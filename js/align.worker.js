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
//   {type: "batch", sequences: [...]}        -- fold-and-discard batch, see below
// Messages out:
//   {type: "ready"}
//   {type: "progress", index, total, name, phase}
//   {type: "summary",  index, summary, alignment}
//   {type: "done",     summaries}
//   {type: "detail",   index, nt, aa, paths}
//   {type: "batch",   matrix buffers, names, summaries, geometry}
//   {type: "error",    message, index}

import { alignToReference, buildSeedIndex, MAX_BATCH_SEQUENCES, MAX_DETAIL_SEQUENCES }
  from "./align.js";
import { makeBatchMatrix, foldSequenceIntoMatrix, codonByPositionFrom,
         buildAaIndex, makeAaMatrix, foldAaIntoMatrix } from "./batch_core.js";
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
    // BATCH. The difference from "run" is that nothing is retained: each sequence's atlas is
    // folded into the per-site matrix and then dropped, so worker memory stays O(sites x
    // samples) rather than O(rows x samples). Retaining 100 full-genome atlases the way "run"
    // does would be ~286 MB here; the matrix is ~8 MB. That is why this is a separate path and
    // not a flag on the existing one.
    if (msg.type === "batch") {
      if (!engine) throw new Error("worker received 'batch' before 'init'");
      const seqs = msg.sequences.slice(0, MAX_BATCH_SEQUENCES);
      const region = reference.modelled_region;
      const nSites = region.last_position - region.first_position + 1;
      const matrix = makeBatchMatrix(seqs.length, nSites, region.first_position);
      // The codon-level matrix is what the tracks, the divergence view and the scatter read:
      // mutability is a property of a codon, not of a nucleotide. The nucleotide matrix is
      // kept alongside because the per-sequence rate summary is a per-nucleotide-class
      // quantity. Both together are ~22 MB at the cap, so neither needs to be dropped.
      const aaIndex = buildAaIndex(reference.proteins);
      const aaMatrix = makeAaMatrix(seqs.length, aaIndex.nPositions);
      const names = [], summaries = [], geometry = [];
      for (let i = 0; i < seqs.length; i++) {
        const { name, sequence } = seqs[i];
        names.push(name);
        self.postMessage({ type: "progress", index: i, total: seqs.length, name, phase: "aligning" });
        const alignment = alignToReference(sequence, reference.sequence, { seedIndex });
        if (alignment.unalignable) {
          summaries.push({ name, unalignable: true, reason: alignment.reason });
          geometry.push(null);
          continue;                       // row stays all-NaN: excluded, not zero
        }
        self.postMessage({ type: "progress", index: i, total: seqs.length, name, phase: "evaluating" });
        const out = engine.customAtlas({ name, sequence, alignment });
        foldSequenceIntoMatrix(matrix, i, out.nt,
                               codonByPositionFrom(out.aa, reference.proteins));
        foldAaIntoMatrix(aaMatrix, i, out.paths, aaIndex);
        summaries.push(out.summary);
        geometry.push({ identity: alignment.identity, coverageOfQuery: alignment.coverageOfQuery,
                        refStart: alignment.refStart, refEnd: alignment.refEnd,
                        orientation: alignment.orientation });
        // `out` goes out of scope here and is collectable before the next sequence.
      }
      // destinations is a Map of Float32Array keyed by "positionIndex|codon"; structured clone
      // handles it, and it is small (one entry per distinct codon actually observed).
      self.postMessage({
        type: "batch", names, summaries, geometry,
        nSamples: matrix.nSamples, nSites: matrix.nSites, firstSite: matrix.firstSite,
        total: matrix.total, ga: matrix.ga, covered: matrix.covered, codonCode: matrix.codonCode,
        typeSum: matrix.typeSum, typeCount: matrix.typeCount, sitesCovered: matrix.sitesCovered,
        aa: {
          nSamples: aaMatrix.nSamples, nPositions: aaMatrix.nPositions,
          nonsyn: aaMatrix.nonsyn, syn: aaMatrix.syn, stop: aaMatrix.stop,
          ga: aaMatrix.ga,
          covered: aaMatrix.covered, codonCode: aaMatrix.codonCode, aaCode: aaMatrix.aaCode,
          destinations: aaMatrix.destinations, sitesCovered: aaMatrix.sitesCovered,
        },
        truncated: msg.sequences.length > seqs.length,
      }, [matrix.total.buffer, matrix.ga.buffer, matrix.covered.buffer, matrix.codonCode.buffer,
          matrix.typeSum.buffer, matrix.typeCount.buffer, matrix.sitesCovered.buffer,
          aaMatrix.nonsyn.buffer, aaMatrix.syn.buffer, aaMatrix.stop.buffer,
          aaMatrix.ga.buffer,
          aaMatrix.covered.buffer, aaMatrix.codonCode.buffer, aaMatrix.aaCode.buffer,
          aaMatrix.sitesCovered.buffer]);
      return;
    }

    if (msg.type === "run") {
      if (!engine) throw new Error("worker received 'run' before 'init'");
      results.length = 0;
      const seqs = msg.sequences.slice(0, MAX_DETAIL_SEQUENCES);
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
