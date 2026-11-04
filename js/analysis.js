// SPDX-License-Identifier: MIT
//
// analysis.js -- drive the worker that aligns and evaluates submitted sequences.
//
// Both the sequence view and the comparison view need the same thing: hand over some
// FASTA records, watch progress, get summaries back, then pull one sequence's full
// tables at a time. That logic lives here so neither page reimplements it.
//
// FALLBACK: module workers are widely supported but not universally. If constructing
// one throws, the same code runs on the main thread instead and the caller is told, so
// the tool degrades to "briefly unresponsive" rather than "broken".

import { alignToReference, buildSeedIndex, MAX_BATCH_SEQUENCES, MAX_DETAIL_SEQUENCES }
  from "./align.js";
import { makeBatchMatrix, foldSequenceIntoMatrix, codonByPositionFrom,
         buildAaIndex, makeAaMatrix, foldAaIntoMatrix } from "./batch_core.js";
import { MraEngine } from "./mra_engine.js";

export { MAX_BATCH_SEQUENCES, MAX_DETAIL_SEQUENCES };

export class Analysis {
  constructor(assets) {
    this.assets = assets;
    this.worker = null;
    this.usedWorker = false;
    this.results = [];      // only populated in the main-thread fallback
  }

  async start() {
    try {
      this.worker = new Worker("js/align.worker.js", { type: "module" });
      await new Promise((resolve, reject) => {
        const onMessage = (e) => {
          if (e.data.type === "ready") { this.worker.removeEventListener("message", onMessage); resolve(); }
          else if (e.data.type === "error") { reject(new Error(e.data.message)); }
        };
        this.worker.addEventListener("message", onMessage);
        this.worker.addEventListener("error", (e) => reject(new Error(e.message || "worker failed")));
        this.worker.postMessage({ type: "init", assets: this.assets });
      });
      this.usedWorker = true;
    } catch (err) {
      // No module-worker support, or the worker could not be constructed.
      this.worker = null;
      this.usedWorker = false;
      this.engine = new MraEngine(this.assets);
      this.seedIndex = buildSeedIndex(this.assets.reference.sequence);
    }
    return this.usedWorker;
  }

  /** @param {function} onProgress called with {index, total, name, phase} */
  run(records, onProgress) {
    const capped = records.slice(0, MAX_DETAIL_SEQUENCES);
    const truncated = records.length > capped.length;
    if (this.usedWorker) {
      return new Promise((resolve, reject) => {
        const summaries = [], geometry = [];
        const onMessage = (e) => {
          const m = e.data;
          if (m.type === "progress") onProgress && onProgress(m);
          else if (m.type === "summary") { summaries[m.index] = m.summary; geometry[m.index] = m.alignment; }
          else if (m.type === "done") {
            this.worker.removeEventListener("message", onMessage);
            resolve({ summaries, geometry, truncated: m.truncated || truncated });
          } else if (m.type === "error") {
            this.worker.removeEventListener("message", onMessage);
            reject(new Error(m.message));
          }
        };
        this.worker.addEventListener("message", onMessage);
        this.worker.postMessage({ type: "run", sequences: capped });
      });
    }
    // Main-thread fallback. Yield to the event loop between sequences so the progress
    // bar can actually paint.
    const self_ = this;
    return (async () => {
      const summaries = [], geometry = [];
      self_.results = [];
      for (let i = 0; i < capped.length; i++) {
        const { name, sequence } = capped[i];
        onProgress && onProgress({ index: i, total: capped.length, name, phase: "aligning" });
        await new Promise((r) => setTimeout(r, 0));
        const alignment = alignToReference(sequence, self_.assets.reference.sequence,
                                           { seedIndex: self_.seedIndex });
        if (alignment.unalignable) {
          const summary = { name, unalignable: true, reason: alignment.reason };
          self_.results.push({ nt: [], aa: [], paths: [], summary });
          summaries.push(summary); geometry.push(null);
          continue;
        }
        onProgress && onProgress({ index: i, total: capped.length, name, phase: "evaluating" });
        await new Promise((r) => setTimeout(r, 0));
        const out = self_.engine.customAtlas({ name, sequence, alignment });
        self_.results.push(out);
        summaries.push(out.summary);
        geometry.push({
          identity: alignment.identity, coverageOfQuery: alignment.coverageOfQuery,
          alignedCount: alignment.alignedCount, unalignedSpans: alignment.unalignedSpans,
          refStart: alignment.refStart, refEnd: alignment.refEnd,
          orientation: alignment.orientation, band: alignment.band,
        });
      }
      return { summaries, geometry, truncated };
    })();
  }

  /**
   * Fold a whole batch into one per-site x per-sample matrix, retaining no atlases.
   *
   * Use this instead of run()+detail() for batch work: run() keeps every atlas so detail()
   * can serve it later, which for 100 full genomes is ~286 MB inside the worker. This path
   * discards each atlas as soon as it has been folded, so the peak is the matrix (~8 MB).
   * The per-sequence tables are therefore NOT available afterwards -- that is the trade.
   *
   * @returns {Promise<{matrix, names, summaries, geometry, truncated}>}
   */
  batch(records, onProgress) {
    const capped = records.slice(0, MAX_BATCH_SEQUENCES);
    const truncated = records.length > capped.length;
    if (this.usedWorker) {
      return new Promise((resolve, reject) => {
        const onMessage = (e) => {
          const m = e.data;
          if (m.type === "progress") { onProgress && onProgress(m); return; }
          if (m.type === "batch") {
            this.worker.removeEventListener("message", onMessage);
            // Rebuild both matrices from the transferred buffers. typeSum/typeCount/
            // sitesCovered must be carried across or the per-sequence rate summary silently
            // loses its denominators -- they were missing here until the codon layer was
            // added, and only the main-thread fallback was exercising that code.
            resolve({
              matrix: { total: m.total, ga: m.ga, covered: m.covered, codonCode: m.codonCode,
                        typeSum: m.typeSum, typeCount: m.typeCount,
                        sitesCovered: m.sitesCovered,
                        nSamples: m.nSamples, nSites: m.nSites, firstSite: m.firstSite },
              aaMatrix: m.aa,
              names: m.names, summaries: m.summaries, geometry: m.geometry,
              truncated: m.truncated || truncated,
            });
          } else if (m.type === "error") {
            this.worker.removeEventListener("message", onMessage);
            reject(new Error(m.message));
          }
        };
        this.worker.addEventListener("message", onMessage);
        this.worker.postMessage({ type: "batch", sequences: capped });
      });
    }
    // Main-thread fallback, yielding between sequences so the progress bar paints.
    const self_ = this;
    return (async () => {
      const ref = self_.assets.reference;
      const region = ref.modelled_region;
      const nSites = region.last_position - region.first_position + 1;
      const matrix = makeBatchMatrix(capped.length, nSites, region.first_position);
      const aaIndex = buildAaIndex(ref.proteins);
      const aaMatrix = makeAaMatrix(capped.length, aaIndex.nPositions);
      const names = [], summaries = [], geometry = [];
      for (let i = 0; i < capped.length; i++) {
        const { name, sequence } = capped[i];
        names.push(name);
        onProgress && onProgress({ index: i, total: capped.length, name, phase: "aligning" });
        await new Promise((r) => setTimeout(r, 0));
        const alignment = alignToReference(sequence, ref.sequence, { seedIndex: self_.seedIndex });
        if (alignment.unalignable) {
          summaries.push({ name, unalignable: true, reason: alignment.reason });
          geometry.push(null);
          continue;
        }
        onProgress && onProgress({ index: i, total: capped.length, name, phase: "evaluating" });
        await new Promise((r) => setTimeout(r, 0));
        const out = self_.engine.customAtlas({ name, sequence, alignment });
        foldSequenceIntoMatrix(matrix, i, out.nt, codonByPositionFrom(out.aa, ref.proteins));
        foldAaIntoMatrix(aaMatrix, i, out.aa, aaIndex);
        summaries.push(out.summary);
        geometry.push({ identity: alignment.identity, coverageOfQuery: alignment.coverageOfQuery,
                        refStart: alignment.refStart, refEnd: alignment.refEnd,
                        orientation: alignment.orientation });
      }
      return { matrix, aaMatrix, names, summaries, geometry, truncated };
    })();
  }

  /** Full tables for one sequence. */
  detail(index) {
    if (!this.usedWorker) {
      const r = this.results[index];
      return Promise.resolve({ nt: r.nt, aa: r.aa, paths: r.paths });
    }
    return new Promise((resolve, reject) => {
      const onMessage = (e) => {
        if (e.data.type === "detail" && e.data.index === index) {
          this.worker.removeEventListener("message", onMessage);
          resolve({ nt: e.data.nt, aa: e.data.aa, paths: e.data.paths });
        } else if (e.data.type === "error") {
          this.worker.removeEventListener("message", onMessage);
          reject(new Error(e.data.message));
        }
      };
      this.worker.addEventListener("message", onMessage);
      this.worker.postMessage({ type: "detail", index });
    });
  }

  /** Release the worker and, with it, every copy of the submitted sequences. */
  dispose() {
    if (this.worker) { this.worker.terminate(); this.worker = null; }
    this.results = [];
  }
}

/** A progress bar bound to an Analysis run. */
export function progressRow(container) {
  container.innerHTML = "";
  const row = document.createElement("div");
  row.className = "progress-row";
  const bar = document.createElement("progress");
  bar.max = 1; bar.value = 0;
  const text = document.createElement("span");
  row.appendChild(bar); row.appendChild(text);
  container.appendChild(row);
  return {
    update({ index, total, name, phase }) {
      bar.max = total; bar.value = index;
      text.textContent = `${phase} ${index + 1} of ${total}: ${name}`;
    },
    finish(message) { bar.value = bar.max; text.textContent = message; },
  };
}
