// SPDX-License-Identifier: MIT
//
// align.worker.js -- runs alignment AND model evaluation off the main thread.
//
// Alignment is the only expensive step: a near-full-length query against the 8.6 kb
// reference takes a few hundred milliseconds, which would visibly freeze the page.
// Evaluating the model afterwards is milliseconds, so it happens here too rather than
// bouncing the alignment back to the main thread and the result forward again.
//
// Results are KEPT IN THE WORKER and summaries are posted as each sequence finishes.
// The page asks for one sequence's full tables at a time. A batch of 100 near-full-length
// genomes would otherwise mean millions of rows crossing the thread boundary at once.
//
// PRIVACY: this worker never makes a network request. It receives the assets the page
// already fetched from this origin and the sequence the user pasted, and it posts
// results back. Nothing leaves the machine.
//
// Messages in:
//   {type: "init",   assets: {bundle, reference, intervals, hxb2}}
//   {type: "run",    sequences: [{name, sequence}]}
//   {type: "detail", index: n}
// Messages out:
//   {type: "ready"}
//   {type: "progress", index, total, name, phase}
//   {type: "summary",  index, summary, alignment}
//   {type: "done",     summaries}
//   {type: "detail",   index, nt, aa, paths}
//   {type: "error",    message, index}

import { alignToReference, buildSeedIndex, MAX_BATCH_SEQUENCES } from "./align.js";
import { MraEngine } from "./mra_engine.js";

let engine = null;
let reference = null;
let seedIndex = null;     // built once per init, reused for every sequence in a batch
const results = [];

self.onmessage = (event) => {
  const msg = event.data;
  try {
    if (msg.type === "init") {
      engine = new MraEngine(msg.assets);
      reference = msg.assets.reference;
      seedIndex = buildSeedIndex(reference.sequence);
      self.postMessage({ type: "ready" });
      return;
    }
    if (msg.type === "detail") {
      const r = results[msg.index];
      if (!r) throw new Error(`no result at index ${msg.index}`);
      self.postMessage({ type: "detail", index: msg.index, nt: r.nt, aa: r.aa, paths: r.paths });
      return;
    }
    if (msg.type === "run") {
      if (!engine) throw new Error("worker received 'run' before 'init'");
      results.length = 0;
      const seqs = msg.sequences.slice(0, MAX_BATCH_SEQUENCES);
      const summaries = [];
      for (let i = 0; i < seqs.length; i++) {
        const { name, sequence } = seqs[i];
        self.postMessage({ type: "progress", index: i, total: seqs.length, name, phase: "aligning" });
        const alignment = alignToReference(sequence, reference.sequence, { seedIndex });
        if (alignment.unalignable) {
          const summary = { name, unalignable: true, reason: alignment.reason };
          results.push({ nt: [], aa: [], paths: [], summary });
          summaries.push(summary);
          self.postMessage({ type: "summary", index: i, summary, alignment: null });
          continue;
        }
        self.postMessage({ type: "progress", index: i, total: seqs.length, name, phase: "evaluating" });
        const out = engine.customAtlas({ name, sequence, alignment });
        // The page needs the alignment geometry for the genome track, but not the two
        // typed arrays, which are large and of no use outside the engine.
        const geometry = {
          identity: alignment.identity,
          coverageOfQuery: alignment.coverageOfQuery,
          alignedCount: alignment.alignedCount,
          unalignedSpans: alignment.unalignedSpans,
          refStart: alignment.refStart,
          refEnd: alignment.refEnd,
          orientation: alignment.orientation,
          band: alignment.band,
        };
        results.push(out);
        summaries.push(out.summary);
        self.postMessage({ type: "summary", index: i, summary: out.summary, alignment: geometry });
      }
      self.postMessage({ type: "done", summaries, truncated: msg.sequences.length > seqs.length });
      return;
    }
    throw new Error(`unknown message type ${msg.type}`);
  } catch (err) {
    self.postMessage({ type: "error", message: err.message, index: msg.index });
  }
};
