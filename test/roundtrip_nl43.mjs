// SPDX-License-Identifier: MIT
//
// roundtrip_nl43.mjs -- the decisive test.
//
//   Run the custom-atlas engine on the NL4-3 reference sequence itself and reproduce
//   the published NL4-3 atlas. Same nucleotide rates, same mu_aa, same path counts,
//   same minimum-substitution classification, for every row where the position maps.
//   Any discrepancy is a bug in the engine, not a modelling choice.
//
// This file does the running; test/compare_roundtrip.py does the comparing. It is a
// permanent regression check: re-run it after any change to the engine, the extractor
// or the assets. v0.1 does not leave the lab until it passes.
//
// It deliberately feeds the reference in as if a user had pasted it, alignment and all,
// so the identity path through the aligner is tested too -- an aligner that silently
// shifted the frame by one would otherwise go unnoticed.
//
//   node test/roundtrip_nl43.mjs [--version v0.1]
//
// Writes test/out/engine_nt.csv, engine_aa.csv, engine_paths.csv and run_summary.json.
//
// ---------------------------------------------------------------------------------
// THE ENGINE API THIS TEST PINS DOWN
//
//   js/align.js
//     alignToReference(querySeq, refSeq, options) -> {
//        queryToRef   Int32Array, one entry per query base, 1-based NL4-3 position
//                     or 0 where the base does not align
//        refToQuery   Int32Array of length refSeq.length + 1, inverse mapping, 0 = gap
//        identity     fraction of aligned columns that match
//        alignedCount number of query bases with a reference position
//        unalignedSpans  [{queryStart, queryEnd}] 1-based inclusive
//        refStart, refEnd  the reference span covered
//     }
//
//   js/mra_engine.js
//     new MraEngine({bundle, reference, intervals, hxb2})      // plain objects, no I/O
//     engine.customAtlas({name, sequence, alignment}) -> {
//        nt[]      one row per (reference position, substitution) that could be evaluated
//        aa[]      one row per (protein, position, alt_aa)
//        paths[]   one row per codon -> mutant-codon path
//        summary   counts, identity, dropped-outside-region, fallback usage
//     }
// ---------------------------------------------------------------------------------

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { alignToReference } from "../js/align.js";
import { MraEngine } from "../js/mra_engine.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, "..");

function arg(name, fallback) {
  const i = process.argv.indexOf(name);
  return i >= 0 && process.argv[i + 1] ? process.argv[i + 1] : fallback;
}

const VERSION = arg("--version", "v0.1");
const DATA = path.join(REPO, "data", VERSION);
const OUT = path.join(HERE, "out");

function loadJson(name) {
  const p = path.join(DATA, name);
  if (!fs.existsSync(p)) {
    console.error(`FAIL: missing asset ${p}. Run scripts/build_web_assets.py first.`);
    process.exit(2);
  }
  return JSON.parse(fs.readFileSync(p, "utf8"));
}

// --- CSV writing, no dependencies -------------------------------------------------
function writeCsv(file, rows, columns) {
  const esc = (v) => {
    if (v === null || v === undefined) return "";
    const s = String(v);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  };
  const lines = [columns.join(",")];
  for (const r of rows) lines.push(columns.map((c) => esc(r[c])).join(","));
  fs.writeFileSync(file, lines.join("\n") + "\n");
}

// ----------------------------------------------------------------------------------
console.log(`NL4-3 round-trip, atlas ${VERSION}`);
console.log("=".repeat(60));

const bundle = loadJson("mra_web_model.json");
const reference = loadJson("nl43_reference.json");
const intervals = loadJson("interval_ratios.json");
const hxb2 = loadJson("hxb2_align.json");

// Rule 4 of the version spine, checked here too: a bundle whose embedded version
// disagrees with its directory means assets were copied between releases.
if (bundle.version !== VERSION) {
  console.error(`FAIL: bundle declares version ${bundle.version} but sits in ${VERSION}/.`);
  process.exit(2);
}

// The submitted sequence IS the reference, pasted in as a user would paste it.
const query = reference.sequence;
console.log(`submitted sequence: ${query.length} bp (the NL4-3 reference itself)`);

const t0 = Date.now();
const alignment = alignToReference(query, reference.sequence, { band: 64 });
const tAlign = Date.now() - t0;
console.log(`alignment: identity ${(alignment.identity * 100).toFixed(3)}%, ` +
            `${alignment.alignedCount} bases mapped, ` +
            `${alignment.unalignedSpans.length} unaligned span(s), ${tAlign} ms`);

// An identity alignment must be exactly that. If this trips, nothing downstream means
// anything, so fail here rather than produce a plausible-looking table.
if (alignment.identity !== 1 || alignment.alignedCount !== query.length) {
  console.error("FAIL: the reference did not align to itself perfectly. " +
                `identity=${alignment.identity}, aligned=${alignment.alignedCount}/${query.length}`);
  process.exit(1);
}
for (let q = 1; q <= query.length; q++) {
  if (alignment.queryToRef[q - 1] !== q) {
    console.error(`FAIL: identity alignment is offset at query base ${q}: ` +
                  `maps to reference ${alignment.queryToRef[q - 1]}`);
    process.exit(1);
  }
}

const engine = new MraEngine({ bundle, reference, intervals, hxb2 });
const t1 = Date.now();
const result = engine.customAtlas({ name: "NL4-3", sequence: query, alignment });
const tEngine = Date.now() - t1;

console.log(`engine: ${result.nt.length} nucleotide rows, ${result.aa.length} amino-acid rows, ` +
            `${result.paths.length} paths, ${tEngine} ms`);
console.log(`dropped outside modelled region: ${result.summary.droppedOutsideRegion} positions`);
console.log(`spline fallback used at: ${result.summary.fallbackPositions} positions`);

fs.mkdirSync(OUT, { recursive: true });
writeCsv(path.join(OUT, "engine_nt.csv"), result.nt,
  ["ref_position", "query_position", "ref_nucleotide", "mut_type", "context_3p",
   "rate", "lo95", "hi95", "spline_class", "position_source", "interval_source"]);
writeCsv(path.join(OUT, "engine_aa.csv"), result.aa,
  ["protein", "position", "aa_change", "ref_aa", "alt_aa", "synonymous",
   "mu_aa", "mu_aa_lo95", "mu_aa_hi95", "min_nt_changes", "n_paths",
   "codon", "query_codon_number", "query_nt_start", "hxb2_label"]);
writeCsv(path.join(OUT, "engine_paths.csv"), result.paths,
  ["protein", "position", "codon", "mutant_codon", "ref_aa", "alt_aa",
   "n_nt_changes", "mutation_type_path", "mu_path", "mu_path_lo95", "mu_path_hi95"]);

fs.writeFileSync(path.join(OUT, "run_summary.json"), JSON.stringify({
  version: VERSION,
  model: bundle.model,
  query_length: query.length,
  align_ms: tAlign,
  engine_ms: tEngine,
  ...result.summary,
}, null, 2) + "\n");

// Structural expectations that do not need the published tables to check. These are
// facts about the reference, so a change here means the engine lost codons or types.
const EXPECT = { codons: 2984, paths: 2984 * 63, aaRows: 62518, ntRows: 8213 * 3 };
let structural = true;
for (const [what, want] of [["paths", EXPECT.paths], ["aa", EXPECT.aaRows], ["nt", EXPECT.ntRows]]) {
  const got = result[what].length;
  if (got !== want) {
    console.error(`FAIL: expected ${want} ${what} rows for NL4-3, got ${got}`);
    structural = false;
  }
}
console.log(structural ? "\nstructural row counts OK" : "\nstructural row counts WRONG");
console.log(`wrote ${OUT}/engine_{nt,aa,paths}.csv`);
console.log("\nNow run:  python test/compare_roundtrip.py --version " + VERSION);
process.exit(structural ? 0 : 1);
