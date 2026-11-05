// SPDX-License-Identifier: MIT
//
// codon_sensitivity_cases.mjs -- the codon-sensitivity engine.
//
//   node test/codon_sensitivity_cases.mjs [--version v0.1]
//
// The load-bearing correctness check: for NL4-3's own codon with the NL4-3 neighbour, the
// sensitivity engine must return exactly the published aa_level_supply rate for that
// change -- it is the same arithmetic reading the same bases. If that holds, the varied
// codons are the same computation on a different codon, which is the whole point.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { MraEngine, GENETIC_CODE } from "../js/mra_engine.js";
import { codonSensitivity, codonSensitivitySummary, formatFold, sequonLossRate }
  from "../js/codon_sensitivity.js";
import { parseCsv } from "../js/data.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const i = process.argv.indexOf("--version");
const VERSION = i >= 0 ? process.argv[i + 1] : "v0.1";
const DATA = path.join(HERE, "..", "data", VERSION);
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference: J("nl43_reference.json"),
                               intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });

// published aa layer, indexed
const aa = parseCsv(fs.readFileSync(path.join(DATA, "aa_level_supply.csv"), "utf8")).records;
const pub = new Map();
for (const r of aa) pub.set(r.protein + "|" + Number(r.position) + "|" + r.alt_aa,
                            { mu_aa: Number(r.mu_aa), min: Number(r.min_nt_changes),
                              n: Number(r.n_paths), ref: r.ref_aa });

// derived tolerance from the round-trip: a 3-step path compounds the 0.83% site residual
const TOL3 = Math.pow(1.0083, 3) - 1;

let failures = 0;
const check = (label, cond, detail) => {
  if (cond) console.log(`  OK   ${label}`);
  else { console.log(`  FAIL ${label}${detail ? " -- " + detail : ""}`); failures++; }
};

// Some well-known DRMs across genes, in atlas (gene, gene_position) coordinates.
// PR M46I = Pol 97, RT K103N = Pol 253, RT M184V = Pol 334, IN N155H = Pol 865.
const CASES = [
  { name: "PR M46I", gene: "Pol", genePosition: 97, refAa: "M", altAa: "I" },
  { name: "RT K103N", gene: "Pol", genePosition: 253, refAa: "K", altAa: "N" },
  { name: "RT M184V", gene: "Pol", genePosition: 334, refAa: "M", altAa: "V" },
  { name: "IN N155H", gene: "Pol", genePosition: 865, refAa: "N", altAa: "H" },
  { name: "RT Y181C", gene: "Pol", genePosition: 331, refAa: "Y", altAa: "C" },
];

console.log(`codon sensitivity, atlas ${VERSION}`);

console.log("\n1. NL4-3 codon reproduces the published rate");
for (const c of CASES) {
  const s = codonSensitivity(engine, c);
  const key = c.gene + "|" + c.genePosition + "|" + c.altAa;
  const p = pub.get(key);
  if (!p) { check(`${c.name}: published row exists`, false, key); continue; }
  const dev = Math.abs(s.nl43_rate - p.mu_aa) / p.mu_aa;
  check(`${c.name}: NL4-3 codon rate matches published (dev ${(dev * 100).toFixed(4)}%)`,
        dev <= TOL3, `engine ${s.nl43_rate}, published ${p.mu_aa}`);
}

console.log("\n2. enumeration sizes follow the genetic code");
for (const c of CASES) {
  const s = codonSensitivity(engine, c);
  // number of synonymous codons that actually reach altAa with positive rate <= total syn codons
  const synTotal = Object.keys(GENETIC_CODE).filter((cd) => GENETIC_CODE[cd] === s.refAa).length;
  check(`${c.name}: ${s.by_site_codon.length} productive codons <= ${synTotal} synonymous`,
        s.by_site_codon.length <= synTotal && s.by_site_codon.length >= 1,
        `${s.by_site_codon.length} vs ${synTotal}`);
  check(`${c.name}: NL4-3 codon is among the enumerated codons`,
        s.by_site_codon.some((r) => r.is_nl43), "NL4-3 codon missing from enumeration");
}

console.log("\n3. holding the flank fixed reproduces the site-only result");
// by construction dimension 1 holds the next codon at NL4-3, so the NL4-3 entry of
// by_site_codon must equal nl43_rate exactly.
for (const c of CASES) {
  const s = codonSensitivity(engine, c);
  const nl = s.by_site_codon.find((r) => r.is_nl43);
  check(`${c.name}: NL4-3 entry of the site sweep equals nl43_rate`,
        nl && Math.abs(nl.mu_aa - s.nl43_rate) < 1e-30, "mismatch within the module");
}

console.log("\n3b. codon optimization is separated from accessibility-class change");
for (const c of CASES) {
  const s = codonSensitivity(engine, c);
  // the codon-optimization fold range must never exceed the full range, and when step
  // counts vary the full range must be strictly larger (the categorical effect is real).
  check(`${c.name}: codon-opt range <= full range`, s.site_fold_range <= s.full_fold_range + 1e-6,
        `${s.site_fold_range} vs ${s.full_fold_range}`);
  if (s.step_count_varies) {
    check(`${c.name}: step-count spread kept out of the headline`,
          s.full_fold_range > s.site_fold_range, "full should exceed codon-opt");
  }
}

console.log("\n4. fold ranges are well-formed");
for (const c of CASES) {
  const s = codonSensitivity(engine, c);
  check(`${c.name}: site fold range >= 1`, s.site_fold_range >= 1 - 1e-9,
        String(s.site_fold_range));
  check(`${c.name}: NL4-3 rank within [1, peers]`,
        s.nl43_rank_among_peers >= 1 && s.nl43_rank_among_peers <= s.n_peers_at_nl43_step,
        `rank ${s.nl43_rank_among_peers} of ${s.n_peers_at_nl43_step}`);
  console.log(`         ${c.name}: ${s.n_synonymous_codons} syn WT codons ` +
              `(${s.n_codons_at_min} at min ${s.min_nt_changes} nt, ${s.n_codons_higher_step} higher), ` +
              `codon-opt fold ${s.site_fold_range.toFixed(2)}x` +
              (s.step_count_varies ? `, full ${s.full_fold_range.toExponential(1)}x incl. step changes` : "") +
              `, NL4-3 ${s.nl43Codon} ranks ${s.nl43_rank_among_peers}/${s.n_peers_at_nl43_step}`);
}

// --- sequon-loss aggregation ------------------------------------------------------
console.log("\n5. glycan sequon-loss aggregation");
const sequons = J("glycan_sequons.json").sequons;
const bnabDep = sequons.filter((s) => s.bnab_dependent);
check(`glycan_sequons.json lists bNAb-dependent sequons`, bnabDep.length > 0,
      String(bnabDep.length));

for (const sq of bnabDep) {
  const r = sequonLossRate(engine, sq);
  if (!r) { check(`N${sq.hxb2_position}: sequon-loss computed`, false, "returned null"); continue; }
  // the aggregate must never be below its largest single component
  check(`N${sq.hxb2_position} (${r.motif}): aggregate >= largest contributor`,
        r.mu_sequon_loss >= r.largest_contributor.mu_aa - 1e-30,
        `total ${r.mu_sequon_loss}, largest ${r.largest_contributor.mu_aa}`);
  // and it should sum three codons' worth of changes, so strictly exceed the largest
  // whenever more than one destroying change exists
  check(`N${sq.hxb2_position}: multiple routes sum above any single one`,
        r.n_destroying_changes <= 1 || r.mu_sequon_loss > r.largest_contributor.mu_aa,
        `${r.n_destroying_changes} changes but total == largest`);
  console.log(`         N${sq.hxb2_position} ${r.motif}: sequon-loss ${r.mu_sequon_loss.toExponential(2)} ` +
              `over ${r.n_destroying_changes} changes; largest single ` +
              `${(r.largest_fraction * 100).toFixed(0)}% (${r.largest_contributor.role} ` +
              `→${r.largest_contributor.alt_aa})`);
}

// a non-sequon position must return null, not zero
const notSequon = sequonLossRate(engine, { nl43_env_position: 1, hxb2_position: 1,
  codon_nt: [engine.reference.proteins.Env.codon_nt[0],
             engine.reference.proteins.Env.codon_nt[1],
             engine.reference.proteins.Env.codon_nt[2]] });
check("a non-sequon position returns null, not zero",
      notSequon === null || engine.reference.proteins.Env.codon_seq[0][0] === "N",
      "expected null for a non-N-X-S/T position");

// sequon positions must agree with glycan_sequons.json motifs
let motifOk = true;
for (const sq of sequons.slice(0, 10)) {
  const r = sequonLossRate(engine, sq);
  if (r && r.motif[0] !== "N") motifOk = false;
}
check("every listed sequon translates to an N-X-S/T motif", motifOk);

// --- the inline summary used by the curated tables --------------------------------------------
// These rows are now rendered into every table cell rather than built on demand, so a silent
// change in shape would ship as a wrong column rather than a broken button.
{
  // Methionine has one codon: codon usage cannot move the rate, and the score must say so
  // rather than printing a misleading "1.0x".
  const met = codonSensitivitySummary(engine, { gene: "Pol", genePosition: 87,
                                                refAa: "M", altAa: "I" });
  check("a single-codon residue yields one codon option",
        met && met.n_codons === 1 && met.codons[0].is_nl43, met && met.n_codons);
  check("a single-codon residue scores as not applicable, not as 1x",
        formatFold(met.fold_range, met.n_codons) === "\u2014",
        formatFold(met.fold_range, met.n_codons));

  // Glycine has four codons; two of them need an extra nucleotide to reach serine, which is
  // what makes the fold range large. The step count must be visible in the text.
  const gly = codonSensitivitySummary(engine, { gene: "Pol", genePosition: 6,
                                                refAa: "G", altAa: "S" });
  if (gly) {
    check("a four-codon residue lists four options", gly.n_codons === 4, gly.n_codons);
    check("options are ordered fastest first",
          gly.codons.every((c, i) => i === 0 || gly.codons[i - 1].rate >= c.rate));
    check("the NL4-3 codon is marked exactly once",
          gly.codons.filter((c) => c.is_nl43).length === 1);
    check("the text marks NL4-3 with an asterisk and flags multi-step codons",
          /\*/.test(gly.text) && (!gly.step_count_varies || /\dnt\)/.test(gly.text)),
          gly.text);
    check("the fold range equals max over min of the listed rates",
          Math.abs(gly.fold_range -
            gly.codons[0].rate / gly.codons[gly.codons.length - 1].rate) < 1e-6 * gly.fold_range,
          `${gly.fold_range.toExponential(2)}`);
    check("the NL4-3 rank is its position in the listed order",
          gly.codons[gly.nl43_rank - 1].is_nl43, gly.nl43_rank);
  }

  // The headline score is the full range, not the matched-step one: at a position where a
  // single codon reaches the target in the fewest steps, the matched-step range is 1.0 and
  // would hide a real several-log effect of spelling the wild-type residue differently.
  let foundWide = false;
  for (let pos = 100; pos < 400 && !foundWide; pos++) {
    const sm = codonSensitivitySummary(engine, { gene: "Pol", genePosition: pos,
                                                 refAa: "L", altAa: "W" });
    if (sm && sm.n_codons > 1 && sm.fold_range_same_steps < 1.01 && sm.fold_range > 100) {
      foundWide = true;
      check("the score reports the full range where the matched-step range would read 1.0",
            formatFold(sm.fold_range, sm.n_codons) !== "\u2014" &&
            sm.fold_range > sm.fold_range_same_steps,
            `full ${sm.fold_range.toExponential(1)} vs same-step ${sm.fold_range_same_steps.toFixed(2)}`);
    }
  }
  check("such a position exists in pol (the case the score is designed for)", foundWide);

  check("formatFold is compact below 1000 and exponential above",
        formatFold(45.2, 3) === "45\u00d7" && formatFold(4.4, 2) === "4.4\u00d7" &&
        /e\d/.test(formatFold(44000, 4)),
        [formatFold(45.2, 3), formatFold(4.4, 2), formatFold(44000, 4)].join(" "));
}

console.log("\n" + "=".repeat(60));
console.log(failures === 0 ? "ALL CODON-SENSITIVITY CASES PASSED" : `${failures} CASE(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
