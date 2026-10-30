// SPDX-License-Identifier: MIT
//
// explore_core.js -- the DOM-free logic behind the genome-explorer page.
//
// Kept separate from explore.js so it can be exercised in node without a browser. Every
// function here is pure: same inputs, same outputs, no module state.

/**
 * Per-site substitution supply, summed over the substitution classes present at that site.
 * Insertions and deletions are deliberately excluded -- they are a different process with a
 * different denominator, and folding them into a "substitution supply" figure misreports both.
 *
 * @param {Map<number, Array>} byPosition  nt position -> atlas rows
 * @returns {Map<number, {total:number, ga:number, nonga:number}>}
 */
export function perSiteSupply(byPosition) {
  const out = new Map();
  for (const [pos, rows] of byPosition) {
    let total = 0, ga = 0;
    for (const r of rows) {
      if (!r.mut_type.includes("->")) continue;
      const v = r.bio_rate_mean;
      if (!(v > 0)) continue;
      total += v;
      if (r.mut_type === "G->A") ga += v;
    }
    if (total > 0) out.set(pos, { total, ga, nonga: total - ga });
  }
  return out;
}

/**
 * Centred running mean of a per-site field over `win` sites.
 *
 * Arithmetic rather than geometric: the window value is the expected per-site supply across
 * the window, which is an arithmetic mean of rates. A geometric mean would answer a different
 * question and would be pulled toward the smallest rate in the window.
 *
 * Windows with fewer than max(3, win/3) contributing sites are dropped rather than reported
 * on thin support, which is what makes the two ends of the track taper off instead of spiking.
 */
export function smoothSeries(sites, perSite, field, win) {
  const half = Math.floor(win / 2);
  const need = Math.max(3, Math.floor(win / 3));
  const out = [];
  for (let i = 0; i < sites.length; i++) {
    let sum = 0, n = 0;
    for (let j = i - half; j <= i + half; j++) {
      if (j < 0 || j >= sites.length) continue;
      sum += perSite.get(sites[j])[field]; n++;
    }
    if (n >= need) out.push([sites[i], sum / n]);
  }
  return out;
}

/**
 * Reverse index from nucleotide to every reading frame that codes at it. HIV genes overlap
 * (gag/pol, and tat/rev/vpu/env), so one nucleotide can belong to several frames; the caller
 * decides which to display. Codons absent from the amino-acid layer are skipped, so the index
 * only ever points at positions the atlas can actually answer for.
 */
export function buildFrameIndex(proteins, aaByPosition) {
  const frameAt = new Map();
  for (const [protein, rec] of Object.entries(proteins)) {
    const codons = rec.codon_nt;
    if (!codons) continue;
    for (let i = 0; i < codons.length; i++) {
      const geneAa = rec.first_aa + i;
      if (!aaByPosition.has(protein + "|" + geneAa)) continue;
      for (const nt of codons[i]) {
        let list = frameAt.get(nt);
        if (!list) frameAt.set(nt, (list = []));
        list.push({ protein, geneAa, codonIndex: i });
      }
    }
  }
  return frameAt;
}

/** Gene spans wide enough to be worth shading on the track. */
export function geneSpans(proteins, minWidth = 300) {
  const spans = [];
  for (const [protein, rec] of Object.entries(proteins)) {
    const c = rec.codon_nt;
    if (!c || !c.length) continue;
    const start = Math.min(...c[0]), end = Math.max(...c[c.length - 1]);
    if (end - start > minWidth) spans.push({ start, end, label: protein });
  }
  return spans.sort((a, b) => a.start - b.start);
}

/** One-step amino-acid changes at a codon, largest supply first -- the pie's slices. */
export function oneStepSlices(rows) {
  return rows.filter((r) => r.min_nt_changes === 1).sort((a, b) => b.mu_aa - a.mu_aa);
}
