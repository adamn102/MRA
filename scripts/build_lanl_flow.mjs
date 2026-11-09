// SPDX-License-Identifier: MIT
//
// build_lanl_flow.mjs -- precompute the LANL batch so "Browse the atlas" can show a supply
// flow chart at any position without the visitor running a batch first.
//
//     node scripts/build_lanl_flow.mjs [--version v0.1]
//
// WHY PRECOMPUTE. The flow chart answers "across circulating virus, where does residue change
// at this position actually get generated?", and that question needs a real population of
// codons -- it does not exist for a single reference sequence. Producing it in the browser
// means aligning and scoring 182 compendium genomes, which takes ~28 s on this machine. That
// is reasonable to ask of someone who has come to the batch tab to run a batch; it is not
// reasonable on the landing page, where the chart should be there when a position is picked.
// So the batch is run once, here, and only its end results ship.
//
// WHAT SHIPS. Per amino-acid position: how many of the 182 sequences cover it, which codons
// they carry and how often, and for each of those codons the rate to each destination
// residue. Nothing else from the batch -- no per-sample matrices, no alignments, no genome
// tracks. The codon-to-residue mapping is dropped too, being the standard genetic code the
// app already has. Rates are rounded to 3 significant figures, far finer than the chart can
// draw, which is what makes the file compress.
//
// This is the same pipeline the batch tab runs, imported from the same modules, so the chart
// on the landing page and the chart on the batch tab are computed identically by
// construction. test/browse_flow_cases.mjs checks the shipped file against a fresh run.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const REPO = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const argv = process.argv.slice(2);
const version = argv.includes("--version") ? argv[argv.indexOf("--version") + 1] : "v0.1";
const DATA = path.join(REPO, "data", version);
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));

const { MraEngine } = await import(path.join(REPO, "js/mra_engine.js"));
const { Analysis } = await import(path.join(REPO, "js/analysis.js"));
const { parseFasta } = await import(path.join(REPO, "js/align.js"));
const { buildAaIndex, aaPositionStats, AA_ORDER } =
  await import(path.join(REPO, "js/batch_core.js"));

const assets = {
  bundle: J("mra_web_model.json"),
  reference: J("nl43_reference.json"),
  intervals: J("interval_ratios.json"),
  hxb2: J("hxb2_align.json"),
};
const engine = new MraEngine(assets);
void engine;                      // constructed to fail fast on a malformed bundle

const fastaPath = path.join(DATA, "lanl_batch_demo.fasta");
const records = parseFasta(fs.readFileSync(fastaPath, "utf8"));
process.stdout.write(`${records.length} sequences from ${path.basename(fastaPath)}\n`);

const analysis = new Analysis(assets);
await analysis.start();           // no Worker under Node: falls back to the main thread
const t0 = Date.now();
const batch = await analysis.batch(records, () => {});
process.stdout.write(`batch folded in ${((Date.now() - t0) / 1000).toFixed(1)} s\n`);

const aaIndex = buildAaIndex(assets.reference.proteins);
// "nonsyn" is the field the flow chart's own page uses; the flow chart reads the destination
// vectors rather than this summary, but the per-position codon counts come from here.
const stats = aaPositionStats(batch.aaMatrix, aaIndex, { field: "nonsyn" });

const sig3 = (v) => Number(v.toPrecision(3));
const rows = [];
for (const s of stats) {
  const codons = [];
  for (const [codon, count] of s.codon_counts) {
    const d = batch.aaMatrix.destinations.get(s.position_index + "|" + codon);
    if (!d) continue;             // a codon that never folded; should not occur
    const flat = [];
    for (let i = 0; i < AA_ORDER.length; i++) {
      if (d[i] > 0) { flat.push(i, sig3(d[i])); }
    }
    codons.push([codon, count, flat]);
  }
  if (codons.length) rows.push([s.protein, s.position, s.n_covered, codons]);
}

const out = {
  source: path.basename(fastaPath),
  n_sequences: records.length,
  aa_order: AA_ORDER.join(""),
  n_positions: rows.length,
  format: "rows: [protein, position, n_covered, [[codon, count, [aaIndex, rate, ...]], ...]]",
  rows,
};
const dest = path.join(DATA, "lanl_flow.json");
fs.writeFileSync(dest, JSON.stringify(out));
const bytes = fs.statSync(dest).size;
process.stdout.write(`wrote ${path.relative(REPO, dest)}  ` +
  `${rows.length} positions, ${(bytes / 1e6).toFixed(2)} MB\n`);
process.stdout.write("Now re-run scripts/build_web_assets.py to refresh the manifest.\n");
