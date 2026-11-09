// SPDX-License-Identifier: MIT
//
// lanl_flow.js -- the precomputed LANL batch, reshaped into what drawSupplyFlow expects.
//
// The supply-flow chart was built for the batch tab, where the visitor has just submitted an
// alignment and the batch is in memory. "Browse the atlas" has no batch, and running one on
// the landing page would cost ~28 s of alignment before anything could be drawn. So the batch
// over the 182 LANL group M compendium genomes is run once at build time
// (scripts/build_lanl_flow.mjs) and only its end results ship.
//
// This module turns that file back into the two objects the chart reads -- a per-position
// `stat` and a `matrix` carrying destination rates -- so drawSupplyFlow is the same function
// drawing the same chart from the same numbers on both tabs. Nothing about the chart is
// reimplemented here; if the two ever disagreed it would be a bug in this adapter, which is
// why test/browse_flow_cases.mjs checks a position end to end against a live batch.
//
// Fetched lazily, on the first position a visitor selects: it is ~370 KB gzipped, and a
// visitor who only reads the genome track should not pay for it.

import { assetUrl } from "./version.js";
import { GENETIC_CODE } from "./mra_engine.js";

let _cache = null;        // version -> Promise of the parsed index

/**
 * Load and index the precomputed batch. Returns null if the asset is not shipped in this
 * release, so an older version of the atlas degrades to no chart rather than a broken page.
 */
export function loadLanlFlow(version) {
  if (_cache && _cache.version === version) return _cache.promise;
  const promise = (async () => {
    let resp;
    try {
      resp = await fetch(assetUrl(version, "lanl_flow.json"));
    } catch {
      return null;                                   // offline, or blocked
    }
    if (!resp.ok) return null;                       // 404: not shipped in this release
    const raw = await resp.json();
    const byKey = new Map();
    for (const [protein, position, nCovered, codons] of raw.rows) {
      byKey.set(protein + "|" + position, { protein, position, nCovered, codons });
    }
    return { meta: { source: raw.source, nSequences: raw.n_sequences,
                     aaOrder: raw.aa_order, nPositions: raw.n_positions },
             byKey };
  })();
  _cache = { version, promise };
  return promise;
}

/**
 * The `stat` and `matrix` pair for one protein position, or null where the compendium does
 * not cover it.
 *
 * `position_index` is set to 0 and the matrix holds only this position: the chart keys its
 * destination lookups on `position_index + "|" + codon`, so a single-position matrix is a
 * complete and much smaller input than the batch tab's full one. The chart cannot tell the
 * difference, which is the point.
 */
export function lanlFlowAt(index, protein, position) {
  if (!index) return null;
  const rec = index.byKey.get(protein + "|" + position);
  if (!rec || !rec.codons.length) return null;

  const aaOrder = index.meta.aaOrder.split("");
  const destinations = new Map();
  const codonCounts = [];
  const codonAa = new Map();
  for (const [codon, count, flat] of rec.codons) {
    const vec = new Float64Array(aaOrder.length);
    for (let i = 0; i < flat.length; i += 2) vec[flat[i]] = flat[i + 1];
    destinations.set("0|" + codon, vec);
    codonCounts.push([codon, count]);
    // The codon-to-residue map is the standard genetic code rather than something shipped:
    // storing it would be 18,601 redundant copies of a table the engine already exports.
    codonAa.set(codon, GENETIC_CODE[codon] || "?");
  }
  codonCounts.sort((a, b) => b[1] - a[1]);

  const stat = {
    key: protein + "|" + position, position_index: 0,
    protein, position,
    n_covered: rec.nCovered,
    n_codons: codonCounts.length,
    codon_counts: codonCounts,
    codon_aa: codonAa,
    majority_codon: codonCounts[0][0],
    majority_codon_count: codonCounts[0][1],
    majority_codon_freq: codonCounts[0][1] / rec.nCovered,
  };
  return { stat, matrix: { destinations } };
}
