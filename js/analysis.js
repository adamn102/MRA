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

import { alignToReference, buildSeedIndex, MAX_BATCH_SEQUENCES } from "./align.js";
import { MraEngine } from "./mra_engine.js";

export { MAX_BATCH_SEQUENCES };

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
    const capped = records.slice(0, MAX_BATCH_SEQUENCES);
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
