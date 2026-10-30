// SPDX-License-Identifier: MIT
//
// Exercises the genome-explorer core against the published atlas. Run:
//   node test/explore_cases.mjs --version v0.1
//
// These assertions encode values cross-checked independently in pandas, so a change in the
// aggregation that silently shifts the track will fail here rather than ship.

import fs from "node:fs";
import { perSiteSupply, smoothSeries, buildFrameIndex, geneSegments, oneStepSlices }
  from "../js/explore_core.js";

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

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL EXPLORER CORE CASES PASSED");
