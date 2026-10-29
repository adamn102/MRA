// SPDX-License-Identifier: MIT
// align_timing.mjs -- measure alignment cost against query length, so MAX_BATCH_SEQUENCES
// is set from data rather than guessed. Writes test/out/align_timing.json.
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { alignToReference, buildSeedIndex } from "../js/align.js";
import { MraEngine } from "../js/mra_engine.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const DATA = path.join(HERE, "..", "data", "v0.1");
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));
const reference = J("nl43_reference.json");
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference,
                               intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });
const REF = reference.sequence;
const seedIndex = buildSeedIndex(REF);

function mutate(seq, rate, seed) {
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
  const bases = ["A", "C", "G", "T"];
  return [...seq].map((b) => (rnd() < rate ? bases[Math.floor(rnd() * 4)] : b)).join("");
}

const lengths = [300, 600, 1000, 1500, 2500, 4000, 6000, 8618];
const REPS = 5;
const rows = [];
for (const L of lengths) {
  const start = Math.max(1, Math.floor((REF.length - L) / 2));
  for (let r = 0; r < REPS; r++) {
    const q = mutate(REF.slice(start - 1, start - 1 + L), 0.05, 1000 + r);
    const t0 = performance.now();
    const aln = alignToReference(q, REF, { seedIndex });
    const t1 = performance.now();
    const out = engine.customAtlas({ name: "t", sequence: q, alignment: aln });
    const t2 = performance.now();
    rows.push({ length: L, rep: r, align_ms: t1 - t0, engine_ms: t2 - t1,
                nt_rows: out.nt.length, aa_rows: out.aa.length, codons: out.summary.codonsAnalysed });
  }
  const mine = rows.filter((x) => x.length === L);
  const med = (k) => mine.map((x) => x[k]).sort((a, b) => a - b)[Math.floor(REPS / 2)];
  console.log(`${String(L).padStart(5)} bp: align ${med("align_ms").toFixed(1)} ms, ` +
              `engine ${med("engine_ms").toFixed(1)} ms, ${mine[0].codons} codons`);
}
fs.mkdirSync(path.join(HERE, "out"), { recursive: true });
fs.writeFileSync(path.join(HERE, "out", "align_timing.json"), JSON.stringify(rows, null, 1) + "\n");
console.log(`\nwrote test/out/align_timing.json (${rows.length} measurements)`);
