// SPDX-License-Identifier: MIT
//
// demo_views.mjs -- run the engine on a real, divergent sequence and dump the data
// behind the comparison and genome-track views, so the figures in docs/ are generated
// by the same code the portal runs rather than by a separate reimplementation.
//
//   node test/demo_views.mjs test/fixtures/JRFL_pol.fasta
//
// Writes test/out/demo_views.json.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { alignToReference, parseFasta } from "../js/align.js";
import { MraEngine, GENETIC_CODE } from "../js/mra_engine.js";
import { parseCsv } from "../js/data.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const DATA = path.join(HERE, "..", "data", "v0.1");
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));
const reference = J("nl43_reference.json");
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference,
                               intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });

const fastaPath = process.argv[2] || path.join(HERE, "fixtures", "JRFL_pol.fasta");
const record = parseFasta(fs.readFileSync(fastaPath, "utf8"))[0];
console.log(`${record.name}: ${record.sequence.length} bp`);

const alignment = alignToReference(record.sequence, reference.sequence);
console.log(`aligned to NL4-3 ${alignment.refStart}-${alignment.refEnd}, ` +
            `identity ${(alignment.identity * 100).toFixed(2)}%, ` +
            `${(alignment.coverageOfQuery * 100).toFixed(1)}% of query placed, ` +
            `orientation ${alignment.orientation}`);
const out = engine.customAtlas({ name: record.name, sequence: record.sequence, alignment });
console.log(`${out.summary.codonsAnalysed} codons, ${out.aa.length} aa rows, ` +
            `${out.nt.length} nt rows, ${out.summary.droppedOutsideRegion} positions dropped`);

// published NL4-3 layers, for comparison
const aa = parseCsv(fs.readFileSync(path.join(DATA, "aa_level_supply.csv"), "utf8")).records;
const pubAa = new Map();
for (const r of aa) pubAa.set(r.protein + "|" + Number(r.position) + "|" + r.alt_aa, r);
const nt = parseCsv(fs.readFileSync(path.join(DATA, "nt_level_supply.csv"), "utf8")).records;

const BASES = ["A", "C", "G", "T"];
function oneStep(codon) {
  const s = new Set();
  for (let i = 0; i < 3; i++) for (const b of BASES) {
    if (b !== codon[i]) s.add(GENETIC_CODE[codon.slice(0, i) + b + codon.slice(i + 1)]);
  }
  s.delete(GENETIC_CODE[codon]);
  return s;
}

// --- codon differences -----------------------------------------------------------
const codonDiffs = [];
const seen = new Set();
for (const r of out.aa) {
  const key = r.protein + "|" + r.position;
  if (seen.has(key)) continue;
  seen.add(key);
  const prot = reference.proteins[r.protein];
  const refCodon = prot.codon_seq[r.position - prot.first_aa];
  if (refCodon === r.codon) continue;
  const refAa = GENETIC_CODE[refCodon], yourAa = GENETIC_CODE[r.codon];
  const a = oneStep(refCodon), b = oneStep(r.codon);
  codonDiffs.push({
    hxb2_label: r.hxb2_label, protein: r.protein, position: r.position,
    ref_codon: refCodon, your_codon: r.codon, ref_aa: refAa, your_aa: yourAa,
    synonymous: refAa === yourAa,
    gained: [...b].filter((x) => !a.has(x)).sort(),
    lost: [...a].filter((x) => !b.has(x)).sort(),
  });
}

// --- accessibility ----------------------------------------------------------------
const access = [];
for (const r of out.aa) {
  const pub = pubAa.get(r.protein + "|" + r.position + "|" + r.alt_aa);
  const prot = reference.proteins[r.protein];
  const refCodon = prot.codon_seq[r.position - prot.first_aa];
  if (!pub) continue;
  access.push({
    hxb2_label: r.hxb2_label, aa_change: r.aa_change,
    ref_codon: refCodon, your_codon: r.codon,
    mu_custom: r.mu_aa, mu_nl43: Number(pub.mu_aa),
    log2fc: Math.log2(r.mu_aa / Number(pub.mu_aa)),
    min_nt_custom: r.min_nt_changes, min_nt_nl43: Number(pub.min_nt_changes),
  });
}

// --- genome track, one series per type -------------------------------------------
const track = {};
for (const type of ["G->A", "C->T", "A->G"]) {
  track[type] = {
    nl43: nt.filter((r) => r.mut_type === type && Number(r.bio_rate_median) > 0)
            .map((r) => [Number(r.ref_position), Number(r.bio_rate_median)])
            .sort((a, b) => a[0] - b[0]),
    custom: out.nt.filter((r) => r.mut_type === type && r.rate > 0)
                  .map((r) => [r.ref_position, r.rate]).sort((a, b) => a[0] - b[0]),
  };
}

fs.mkdirSync(path.join(HERE, "out"), { recursive: true });
fs.writeFileSync(path.join(HERE, "out", "demo_views.json"), JSON.stringify({
  name: record.name, length: record.sequence.length,
  alignment: { identity: alignment.identity, coverage: alignment.coverageOfQuery,
               refStart: alignment.refStart, refEnd: alignment.refEnd,
               orientation: alignment.orientation,
               unalignedSpans: alignment.unalignedSpans },
  summary: out.summary, codonDiffs, access, track,
}));
const syn = codonDiffs.filter((d) => d.synonymous);
console.log(`codon differences: ${codonDiffs.length} (${syn.length} synonymous, ` +
            `${syn.filter((d) => d.gained.length || d.lost.length).length} of those reshape ` +
            `one-step reachability)`);
console.log(`accessibility rows: ${access.length}; ` +
            `${access.filter((r) => Math.abs(r.log2fc) > 1).length} shift more than two-fold; ` +
            `${access.filter((r) => r.min_nt_custom !== r.min_nt_nl43).length} change the ` +
            `required number of substitutions`);
console.log("wrote test/out/demo_views.json");
