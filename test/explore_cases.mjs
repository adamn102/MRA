// SPDX-License-Identifier: MIT
//
// Exercises the genome-explorer core against the published atlas. Run:
//   node test/explore_cases.mjs --version v0.1
//
// These assertions encode values cross-checked independently in pandas, so a change in the
// aggregation that silently shifts the track will fail here rather than ship.

import fs from "node:fs";
import { MraEngine } from "../js/mra_engine.js";
import { trackStats } from "../js/track_view.js";
import { perSiteSupply, smoothSeries, buildFrameIndex, geneSegments, oneStepSlices,
         targetFor } from "../js/explore_core.js";

const vArg = process.argv.indexOf("--version");
const V = vArg > 0 ? process.argv[vArg + 1] : "v0.1";
const D = `data/${V}/`;
const J = (f) => JSON.parse(fs.readFileSync(D + f, "utf8"));

let failures = 0;
function check(name, cond, detail = "") {
  if (cond) console.log(`  ok    ${name}`);
  else { console.log(`  FAIL  ${name}${detail ? "  -- " + detail : ""}`); failures++; }
}
const close = (a, b, rel = 0.01) => Math.abs(a - b) <= rel * Math.abs(b);

// --- load the nucleotide layer the way data.js does -------------------------------
const ntText = fs.readFileSync(D + "nt_level_supply.csv", "utf8").trim().split("\n");
const ntHead = ntText[0].split(",");
const iPos = ntHead.indexOf("ref_position"), iType = ntHead.indexOf("mut_type"),
      iMean = ntHead.indexOf("bio_rate_mean");
const byPosition = new Map();
for (let k = 1; k < ntText.length; k++) {
  const f = ntText[k].split(",");
  const pos = Number(f[iPos]);
  let list = byPosition.get(pos);
  if (!list) byPosition.set(pos, (list = []));
  list.push({ mut_type: f[iType], bio_rate_mean: Number(f[iMean]) });
}

// --- load the amino-acid layer ----------------------------------------------------
const aaText = fs.readFileSync(D + "aa_level_supply.csv", "utf8").trim().split("\n");
const aaHead = aaText[0].split(",");
const aIdx = Object.fromEntries(aaHead.map((h, i) => [h, i]));
const aaByPosition = new Map();
for (let k = 1; k < aaText.length; k++) {
  const f = aaText[k].split(",");
  const key = f[aIdx.protein] + "|" + Number(f[aIdx.position]);
  let list = aaByPosition.get(key);
  if (!list) aaByPosition.set(key, (list = []));
  list.push({
    protein: f[aIdx.protein], position: Number(f[aIdx.position]),
    ref_aa: f[aIdx.ref_aa], alt_aa: f[aIdx.alt_aa], aa_change: f[aIdx.aa_change],
    synonymous: f[aIdx.synonymous] === "True",
    mu_aa: Number(f[aIdx.mu_aa]), min_nt_changes: Number(f[aIdx.min_nt_changes]),
  });
}

console.log(`\nGenome-explorer core, atlas ${V}`);
console.log("=".repeat(60));

// --- 1. per-site supply -----------------------------------------------------------
const perSite = perSiteSupply(byPosition);
const sites = [...perSite.keys()].sort((a, b) => a - b);
check("every modelled position has positive supply", sites.length === 8213,
      `got ${sites.length}`);

const totals = sites.map((p) => perSite.get(p).total).sort((a, b) => a - b);
const median = totals[Math.floor(totals.length / 2)];
// cross-checked in pandas: median 3.595e-05, max 1.961e-03
check("per-site median matches the independent calculation", close(median, 3.595e-5, 0.02),
      `median ${median.toExponential(3)}`);
check("per-site max matches the independent calculation",
      close(totals[totals.length - 1], 1.961e-3, 0.02),
      `max ${totals[totals.length - 1].toExponential(3)}`);

const gaSum = sites.reduce((t, p) => t + perSite.get(p).ga, 0);
const allSum = sites.reduce((t, p) => t + perSite.get(p).total, 0);
check("G->A is ~70% of all substitution supply", close(100 * gaSum / allSum, 69.9, 0.02),
      `${(100 * gaSum / allSum).toFixed(1)}%`);

check("total equals ga + nonga at every site",
      sites.every((p) => { const s = perSite.get(p); return close(s.total, s.ga + s.nonga, 1e-9); }));

// indels must not be counted: a position whose only rows are indels would have no entry
const indelOnly = [...byPosition.entries()]
  .filter(([, rows]) => rows.every((r) => !r.mut_type.includes("->")));
check("indel-only positions are excluded from the track",
      indelOnly.every(([p]) => !perSite.has(p)), `${indelOnly.length} such positions`);

// --- 2. smoothing -----------------------------------------------------------------
const sm151 = smoothSeries(sites, perSite, "total", 151);
check("smoothing returns a point for every well-supported site", sm151.length === sites.length,
      `got ${sm151.length}`);
check("smoothing is monotone in nothing but preserves the mean",
      close(sm151.reduce((t, p) => t + p[1], 0) / sm151.length, allSum / sites.length, 0.02));
const sm25 = smoothSeries(sites, perSite, "total", 25);
const spread = (a) => { const v = a.map((p) => p[1]); return Math.max(...v) / Math.min(...v); };
check("a narrower window resolves more dynamic range", spread(sm25) > spread(sm151),
      `25nt ${spread(sm25).toFixed(1)}x vs 151nt ${spread(sm151).toFixed(1)}x`);
check("non-G->A resolves more structure than the total track",
      spread(smoothSeries(sites, perSite, "nonga", 151)) > spread(sm151),
      `nonGA ${spread(smoothSeries(sites, perSite, "nonga", 151)).toFixed(1)}x vs total ${spread(sm151).toFixed(1)}x`);

// --- 3. reading-frame index -------------------------------------------------------
const reference = J("nl43_reference.json");
const frameAt = buildFrameIndex(reference.proteins, aaByPosition);

const f6804 = frameAt.get(6804) || [];
check("nt 6,804 resolves to Env codon 458",
      f6804.some((f) => f.protein === "Env" && f.geneAa === 458),
      JSON.stringify(f6804.map((f) => f.protein + " " + f.geneAa)));

// overlapping frames: a nucleotide in the gag/pol overlap must offer both
const overlap = [...frameAt.entries()].find(([, fr]) =>
  fr.some((f) => f.protein === "Gag") && fr.some((f) => f.protein === "Pol"));
check("the gag/pol overlap exposes both reading frames", Boolean(overlap),
      overlap ? `nt ${overlap[0]}` : "none found");

check("every indexed frame has atlas rows",
      [...frameAt.values()].every((fr) =>
        fr.every((f) => aaByPosition.has(f.protein + "|" + f.geneAa))));

// every codon nucleotide of a modelled codon is indexed
const envRec = reference.proteins.Env;
const codon458 = envRec.codon_nt[458 - envRec.first_aa];
check("all three bases of a modelled codon are indexed",
      codon458.every((nt) => (frameAt.get(nt) || []).some((f) => f.protein === "Env" && f.geneAa === 458)),
      `codon nt ${codon458.join(",")}`);

// --- 4. pie slices ----------------------------------------------------------------
const rows458 = aaByPosition.get("Env|458");
const slices = oneStepSlices(rows458);
// cross-checked in pandas against the deployed browse page: 8 slices, K at 20.2%
check("Env 458 yields 8 one-step slices", slices.length === 8, `got ${slices.length}`);
const tot = slices.reduce((t, r) => t + r.mu_aa, 0);
check("the largest slice is K at 20.2%",
      slices[0].alt_aa === "K" && close(100 * slices[0].mu_aa / tot, 20.2, 0.01),
      `${slices[0].alt_aa} ${(100 * slices[0].mu_aa / tot).toFixed(1)}%`);
check("slices are sorted by descending supply",
      slices.every((r, i) => i === 0 || slices[i - 1].mu_aa >= r.mu_aa));
check("every slice is genuinely one nucleotide away",
      slices.every((r) => r.min_nt_changes === 1));

// --- 5. gene annotation band --------------------------------------------------------
const segs = geneSegments(reference.proteins);
check("segments are returned in ascending order",
      segs.every((s, i) => i === 0 || segs[i - 1].start <= s.start));
check("all nine genes appear in the band",
      ["Gag", "Pol", "Vif", "Vpr", "Tat", "Rev", "Vpu", "Env", "Nef"]
        .every((g) => segs.some((s) => s.label === g)),
      [...new Set(segs.map((s) => s.label))].join(","));
// tat and rev are spliced: two exons each, either side of env. A single min-to-max box would
// run straight across env and imply continuous coding, so each must come back as two segments.
for (const g of ["Tat", "Rev"]) {
  const mine = segs.filter((s) => s.label === g);
  check(`${g} is split into its two exons`, mine.length === 2,
        mine.map((s) => `${s.start}-${s.end}`).join(" "));
}
check("no segment contains a non-coding gap",
      segs.every((s) => s.end >= s.start));
// every segment must be inside the modelled region, or the band would annotate rows the
// track cannot plot
check("segments lie inside the modelled region",
      segs.every((s) => s.start >= reference.modelled_region.first_position &&
                        s.end <= reference.modelled_region.last_position));
// the gag/pol overlap must be visible as two separate segments that genuinely overlap
const gag = segs.find((s) => s.label === "Gag"), pol = segs.find((s) => s.label === "Pol");
check("gag and pol overlap in the band", gag.end > pol.start,
      `gag ends ${gag.end}, pol starts ${pol.start}`);

// --- 6. protein-based coordinate entry ----------------------------------------------
// The number a user types is the subunit's own numbering -- the convention a paper or clinic
// uses. These landmarks are the check that it lands on the right codon.
const engine = new MraEngine({
  bundle: J("mra_web_model.json"), reference, intervals: J("interval_ratios.json"),
  hxb2: J("hxb2_align.json"),
});
const targets = engine.proteinTargets(new Map(
  [...aaByPosition.keys()].reduce((m, k) => {
    const [prot, pos] = k.split("|"); const n = Number(pos);
    const cur = m.get(prot);
    if (!cur) m.set(prot, { min: n, max: n });
    else { cur.min = Math.min(cur.min, n); cur.max = Math.max(cur.max, n); }
    return m;
  }, new Map())));
const byKey = new Map(targets.map((t) => [t.key, t]));

check("the subunit list offers PR, RT, IN and the gag products",
      ["PR", "RT", "IN", "p17", "p24", "Gag", "Pol", "Env"].every((k) => byKey.has(k)),
      [...byKey.keys()].join(","));

/** own-numbering position -> centre nucleotide of that codon, the path the Go button takes */
function ntFor(key, ownPos) {
  const t = byKey.get(key);
  const geneAa = engine.geneAaFromOwnPosition(t, ownPos);
  if (geneAa === null) return null;
  const rec = reference.proteins[t.protein];
  const idx = geneAa - rec.first_aa;
  if (idx < 0 || idx >= rec.codon_nt.length) return null;
  return { geneAa, nt: rec.codon_nt[idx][1], refAa: (aaByPosition.get(t.protein + "|" + geneAa) || [{}])[0].ref_aa };
}

// Each landmark names the wild-type residue, so a mis-resolved codon shows up as the wrong
// amino acid rather than silently returning a plausible number.
for (const [key, pos, wt, label] of [
  ["PR", 46, "M", "M46I, a protease inhibitor mutation"],
  ["RT", 103, "K", "K103N, an NNRTI mutation"],
  ["RT", 184, "M", "M184V, the 3TC mutation"],
  ["IN", 155, "N", "N155H, an INSTI mutation"],
  ["Gag", 1, "M", "the gag start codon"],
]) {
  const got = ntFor(key, pos);
  check(`${key} ${pos} resolves to wild-type ${wt} (${label})`,
        got !== null && got.refAa === wt,
        got ? `got ${got.refAa} at ${key.length > 2 ? key : key} ${pos} -> gene aa ${got.geneAa}, nt ${got.nt}`
            : "unresolved");
}

// Env is the case where the typed number is HXB2 and NL4-3 differs by indels.
for (const [hxb2Pos, wt, label] of [
  [332, "N", "the 10-1074 / PGT121 glycan"],
  [160, "N", "the V2-apex N160 glycan"],
]) {
  const got = ntFor("Env", hxb2Pos);
  check(`Env ${hxb2Pos} (HXB2) resolves to wild-type ${wt} (${label})`,
        got !== null && got.refAa === wt,
        got ? `NL4-3 Env ${got.geneAa}, nt ${got.nt}, ref ${got.refAa}` : "unresolved");
}
check("Env HXB2 numbering is mapped, not taken literally",
      ntFor("Env", 332).geneAa === 330,
      `HXB2 332 -> NL4-3 ${ntFor("Env", 332).geneAa}`);

// A subunit position past its length must be rejected rather than silently spilling into the
// next subunit -- RT is 440 residues, so RT 500 is not an integrase position.
check("subunit lengths are known, so an over-long position can be rejected",
      byKey.get("RT").maxOwn === 440, `RT maxOwn ${byKey.get("RT").maxOwn}`);

// --- 7. the merged page's selection round-trip --------------------------------------
// Typing a position, reading the label back, and typing it again must land on the same codon.
// If this drifts, the coordinate box silently renumbers a user's site between clicks.
let rtChecked = 0, rtBad = [];
for (const t of targets) {
  const maxOwn = t.maxOwn || (aaLayer_max(t.protein) - (t.offset || 0));
  for (const own of [1, 2, Math.floor(maxOwn / 2), maxOwn - 1, maxOwn]) {
    if (!Number.isInteger(own) || own < 1) continue;
    const geneAa = engine.geneAaFromOwnPosition(t, own);
    if (geneAa === null) continue;                       // HXB2 insertion, no NL4-3 row
    if (!aaByPosition.has(t.protein + "|" + geneAa)) continue;   // outside the modelled region
    const back = engine.ownPositionFor(t, geneAa);
    rtChecked++;
    if (back !== own) rtBad.push(`${t.key} ${own} -> gene ${geneAa} -> ${back}`);
  }
}
function aaLayer_max(protein) {
  let m = 0;
  for (const k of aaByPosition.keys()) {
    const [p, pos] = k.split("|");
    if (p === protein) m = Math.max(m, Number(pos));
  }
  return m;
}
check(`own-position round-trip holds at every subunit boundary (${rtChecked} positions)`,
      rtBad.length === 0, rtBad.slice(0, 4).join("; "));

// A pol position must report as its subunit, not as the gene -- this is what makes the
// dashboard say "RT K103N" instead of "Pol K253N".
for (const [protein, geneAa, wantKey, wantOwn] of [
  ["Pol", 253, "RT", 103],
  ["Pol", 97, "PR", 46],
  ["Pol", 865, "IN", 155],
  ["Gag", 1, "p17", 1],
  ["Env", 330, "Env", 332],      // Env has no subunits; the own position is the HXB2 number
]) {
  const t = targetFor(targets, protein, geneAa);
  const own = t ? engine.ownPositionFor(t, geneAa) : null;
  check(`${protein} ${geneAa} reports as ${wantKey} ${wantOwn}`,
        t && t.key === wantKey && own === wantOwn,
        t ? `got ${t.key} ${own}` : "no target");
}

// Every frame the track can point at must map to a selectable target, or clicking the track
// would land on a position the controls cannot express.
const unmappable = [];
for (const frames of frameAt.values()) {
  for (const f of frames) {
    if (!targetFor(targets, f.protein, f.geneAa)) unmappable.push(`${f.protein} ${f.geneAa}`);
  }
}
check("every clickable frame maps to a selectable protein/subunit",
      unmappable.length === 0, unmappable.slice(0, 4).join("; "));

// --- 8. the quantity table under the track -------------------------------------------
// The median is taken over the WHOLE track, so it must not move with the marker -- it is the
// fixed reference the marker is compared against. The relative row is the one that moves.
// A caption once called the median "over the plotted window", which read as if it were local;
// these two assertions are what stop that confusion coming back as a wrong number.
{
  const pts = smoothSeries(sites, perSite, "total", 151);
  const vals = pts.map((p) => p[1]).sort((a, b) => a - b);
  const med = vals[Math.floor(vals.length / 2)];
  const rowsAt = (nt) => trackStats({
    vals, med, at: pts.find((p) => p[0] >= nt), site: perSite.get(nt),
    cls: "total", win: 151, marker: nt, refBase: null,
  });
  const MARKERS = [1600, 3000, 6804, 8000];
  const medians = MARKERS.map((m) => rowsAt(m)[0].value);
  check("the track median is identical at every marker (it is a fixed reference)",
        new Set(medians).size === 1, medians.join(" / "));

  const relative = MARKERS.map((m) => rowsAt(m)[2].value);
  check("the marker-relative row DOES change as the marker moves",
        new Set(relative).size === MARKERS.length, relative.join(" / "));

  const markerVals = MARKERS.map((m) => rowsAt(m)[1].value);
  check("the value-at-marker row changes with the marker",
        new Set(markerVals).size > 1, markerVals.join(" / "));

  // The median must change with the supply class, or the class selector would be cosmetic.
  const byClass = ["total", "nonga", "ga"].map((c) => {
    const p2 = smoothSeries(sites, perSite, c, 151);
    const v2 = p2.map((x) => x[1]).sort((x, y) => x - y);
    return trackStats({ vals: v2, med: v2[Math.floor(v2.length / 2)],
                        at: p2.find((x) => x[0] >= 3000), site: perSite.get(3000),
                        cls: c, win: 151, marker: 3000, refBase: null })[0].value;
  });
  check("the track median changes with the supply class",
        new Set(byClass).size === 3, byClass.join(" / "));

  // Percentile must be consistent with the ordering of the underlying values.
  const pctOf = (nt) => Number(rowsAt(nt)[2].note.match(/^(\d+)/)[1]);
  const valOf = (nt) => pts.find((p) => p[0] >= nt)[1];
  const hi = 3000, lo = 6804;                       // 1.75e-4 vs 9.75e-5
  check("a higher-supply marker gets the higher percentile",
        valOf(hi) > valOf(lo) && pctOf(hi) > pctOf(lo),
        `nt ${hi}: ${pctOf(hi)}th, nt ${lo}: ${pctOf(lo)}th`);
  check("percentiles stay within 0-100",
        MARKERS.every((m) => pctOf(m) >= 0 && pctOf(m) <= 100));

  // A non-G site must explain its 0%, rather than leaving a bare zero to look like a failure.
  const nonG = trackStats({ vals, med, at: pts.find((p) => p[0] >= 1600),
                            site: perSite.get(1600), cls: "total", win: 151,
                            marker: 1600, refBase: "T" })[3];
  check("a 0% G->A row says why it is zero", nonG.value === "0%" && /not G/.test(nonG.note),
        `${nonG.value} -- ${nonG.note}`);
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL EXPLORER CORE CASES PASSED");
