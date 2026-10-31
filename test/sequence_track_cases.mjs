// SPDX-License-Identifier: MIT
//
// The genome track on the "Your sequence" tab, end to end: align a sequence, build its custom
// atlas, aggregate per-site supply from it, and find where it differs from NL4-3.
//
//   node test/sequence_track_cases.mjs --version v0.1
//
// The point of these cases is that the track on that tab is driven by the SUBMITTED sequence,
// not the published atlas, so the aggregation has to read the custom atlas's own rate column
// and the variant ticks have to be real differences rather than alignment artefacts.

import fs from "node:fs";
import { MraEngine } from "../js/mra_engine.js";
import { alignToReference } from "../js/align.js";
import { perSiteSupply, groupByPosition, buildFrameIndex, targetFor } from "../js/explore_core.js";

const vArg = process.argv.indexOf("--version");
const V = vArg > 0 ? process.argv[vArg + 1] : "v0.1";
const J = (f) => JSON.parse(fs.readFileSync(`data/${V}/${f}`, "utf8"));

let failures = 0;
function check(name, cond, detail = "") {
  if (cond) console.log(`  ok    ${name}`);
  else { console.log(`  FAIL  ${name}${detail ? "  -- " + detail : ""}`); failures++; }
}

const reference = J("nl43_reference.json");
const engine = new MraEngine({
  bundle: J("mra_web_model.json"), reference, intervals: J("interval_ratios.json"),
  hxb2: J("hxb2_align.json"),
});

/** What the page does: custom atlas -> per-site supply + variant positions. */
function viewFor(seq) {
  const out = engine.customAtlas({
    name: "t", sequence: seq, alignment: alignToReference(seq, reference.sequence) });
  const perSite = perSiteSupply(groupByPosition(out.nt), "rate");
  const baseAt = new Map();
  for (const r of out.nt) if (!baseAt.has(r.ref_position)) baseAt.set(r.ref_position, r.ref_nucleotide);
  const variants = [...baseAt.entries()]
    .filter(([p, b]) => b !== reference.sequence[p - 1]).map(([p]) => p).sort((a, b) => a - b);
  return { out, perSite, variants };
}

console.log(`\nSequence-tab genome track, atlas ${V}`);
console.log("=".repeat(60));

// --- 1. the reference itself has no differences -------------------------------------
const refOnly = reference.sequence.slice(0, 8213);
const base = viewFor(refOnly);
check("NL4-3 submitted against itself shows zero differences",
      base.variants.length === 0, `${base.variants.length} found`);
check("every modelled position is covered and carries supply",
      base.perSite.size === 8213, `${base.perSite.size} positions`);
check("the custom track reads the atlas's own rate column, not a missing one",
      [...base.perSite.values()].every((v) => v.total > 0 && Number.isFinite(v.total)));

// --- 2. introduced substitutions are found, exactly ---------------------------------
const PLANTED = [1600, 2001, 6804];
const arr = refOnly.split("");
for (const p of PLANTED) arr[p - 1] = arr[p - 1] === "A" ? "G" : "A";
const mut = viewFor(arr.join(""));
check("planted substitutions are detected at exactly their positions",
      JSON.stringify(mut.variants) === JSON.stringify(PLANTED),
      mut.variants.slice(0, 8).join(","));
check("no spurious differences elsewhere", mut.variants.length === PLANTED.length,
      `${mut.variants.length} found`);

// Changing a base changes that site's supply: the rates are position- and context-specific,
// so a custom track must not simply reproduce the published one.
const changed = PLANTED.filter((p) => base.perSite.get(p).total !== mut.perSite.get(p).total);
check("a changed base changes that site's supply", changed.length === PLANTED.length,
      `${changed.length} of ${PLANTED.length} sites moved`);

// --- 3. a partial amplicon covers only its own span ---------------------------------
const amp = refOnly.slice(2000, 3000);            // reference 2001-3000
const part = viewFor(amp);
const sites = [...part.perSite.keys()].sort((a, b) => a - b);
check("a 1 kb amplicon covers only its own span, not the whole genome",
      sites.length > 900 && sites.length <= 1000 &&
      sites[0] >= 2001 && sites[sites.length - 1] <= 3000,
      `${sites.length} sites, ${sites[0]}-${sites[sites.length - 1]}`);
check("an unmutated amplicon shows no differences", part.variants.length === 0,
      part.variants.join(","));

// --- 4. frames and targets inside the covered span ----------------------------------
const aaByPosition = new Map();
for (const r of part.out.aa) {
  const k = r.protein + "|" + r.position;
  let l = aaByPosition.get(k);
  if (!l) aaByPosition.set(k, (l = []));
  l.push(r);
}
const frameAt = buildFrameIndex(reference.proteins, aaByPosition);
const spans = new Map();
for (const r of part.out.aa) {
  const c = spans.get(r.protein);
  if (!c) spans.set(r.protein, { min: r.position, max: r.position });
  else { c.min = Math.min(c.min, r.position); c.max = Math.max(c.max, r.position); }
}
const targets = engine.proteinTargets(spans);
check("the amplicon's frames are all inside the covered span",
      [...frameAt.keys()].every((nt) => nt >= 2001 && nt <= 3000));
check("every frame in a partial sequence still maps to a selectable target",
      [...frameAt.values()].every((fr) => fr.every((f) => targetFor(targets, f.protein, f.geneAa))));
// reference 2001-3000 is inside pol, so the subunit label must be a pol subunit
const someFrame = [...frameAt.values()][0][0];
const t = targetFor(targets, someFrame.protein, someFrame.geneAa);
check("a pol amplicon reports a pol subunit, not the bare gene",
      t && ["PR", "RT", "RNH", "IN"].includes(t.key),
      t ? `${t.key} ${engine.ownPositionFor(t, someFrame.geneAa)}` : "none");

// --- 5. the pie at a changed codon uses the submitted codon -------------------------
// Position 6804 was planted above; its codon in the mutant differs from NL4-3, so the
// one-step amino-acid options must differ too, or the tab would be showing NL4-3's answer.
const env = reference.proteins.Env;
const geneAa = env.codon_nt.findIndex((c) => c.includes(6804)) + env.first_aa;
const keyAt = "Env|" + geneAa;
const codonOf = (o) => (o.aa.find((r) => r.protein === "Env" && r.position === geneAa) || {}).codon;
check("the planted change alters the submitted codon at that site",
      codonOf(base.out) !== codonOf(mut.out),
      `${codonOf(base.out)} -> ${codonOf(mut.out)}`);
const setOf = (o) => new Set(o.aa.filter((r) => r.protein === "Env" && r.position === geneAa &&
                                                r.min_nt_changes === 1).map((r) => r.alt_aa));
check("the one-step options at that codon change with it",
      JSON.stringify([...setOf(base.out)].sort()) !== JSON.stringify([...setOf(mut.out)].sort()),
      `${[...setOf(base.out)].sort().join("")} vs ${[...setOf(mut.out)].sort().join("")}`);

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL SEQUENCE-TRACK CASES PASSED");
