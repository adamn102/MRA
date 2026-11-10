// SPDX-License-Identifier: MIT
//
// explore_core.js -- the DOM-free logic behind the genome track on the browse page.
//
// Named for the "Explore the genome" tab it was written for. That tab is gone; Browse the
// atlas (index.html) is the landing page, and this module backs its genome track. Kept separate from the page
// module so it can be exercised in node without a browser: every function here is pure --
// same inputs, same outputs, no module state. Exercised by test/explore_cases.mjs.

/**
 * Per-site substitution supply, summed over the substitution classes present at that site.
 * Insertions and deletions are deliberately excluded -- they are a different process with a
 * different denominator, and folding them into a "substitution supply" figure misreports both.
 *
 * @param {Map<number, Array>} byPosition  nt position -> atlas rows
 * @param {string} rateField  the per-row rate column: "bio_rate_mean" in the published
 *        tables, "rate" in a custom atlas built from a submitted sequence
 * @returns {Map<number, {total:number, ga:number, nonga:number}>}
 */
export function perSiteSupply(byPosition, rateField = "bio_rate_mean") {
  const out = new Map();
  for (const [pos, rows] of byPosition) {
    let total = 0, ga = 0;
    for (const r of rows) {
      if (!r.mut_type.includes("->")) continue;
      const v = r[rateField];
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

/**
 * Contiguous coding segments per protein, for the gene annotation band.
 *
 * Taken as runs of consecutive nucleotides rather than first-to-last, because tat and rev are
 * spliced: their first and last codons sit either side of a multi-kilobase intron, so a single
 * box from min to max would draw them straight across env and imply continuous coding. A run
 * per exon is what is actually true. Segments shorter than `minWidth` are dropped as unlabelable.
 */
export function geneSegments(proteins, minWidth = 40) {
  const out = [];
  for (const [protein, rec] of Object.entries(proteins)) {
    const codons = rec.codon_nt;
    if (!codons || !codons.length) continue;
    const nts = [...new Set(codons.flat())].sort((a, b) => a - b);
    let start = nts[0], prev = nts[0];
    for (let i = 1; i <= nts.length; i++) {
      const cur = nts[i];
      if (cur !== prev + 1) {                     // run ended (or array exhausted)
        if (prev - start + 1 >= minWidth) out.push({ start, end: prev, label: protein });
        start = cur;
      }
      prev = cur;
    }
  }
  return out.sort((a, b) => a.start - b.start || a.end - b.end);
}

/** One-step amino-acid changes at a codon, largest supply first -- the pie's slices. */
export function oneStepSlices(rows) {
  return rows.filter((r) => r.min_nt_changes === 1).sort((a, b) => b.mu_aa - a.mu_aa);
}

/**
 * The most specific selectable target for a gene position: the subunit containing it if there
 * is one (so a pol position reports as RT 103 rather than Pol 253), otherwise the gene itself.
 */
export function targetFor(targets, protein, geneAa) {
  let fallback = null;
  for (const t of targets) {
    if (t.protein !== protein) continue;
    if (!t.maxOwn) { fallback = fallback || t; continue; }
    if (geneAa > t.offset && geneAa <= t.offset + t.maxOwn) return t;
  }
  return fallback;
}

/** Group flat atlas rows into the position map perSiteSupply expects. */
export function groupByPosition(rows, key = "ref_position") {
  const out = new Map();
  for (const r of rows) {
    let list = out.get(r[key]);
    if (!list) out.set(r[key], (list = []));
    list.push(r);
  }
  return out;
}

/** Expected de novo copies of a change in N replication events, one cycle. */
export function expectedCopies(mu, n) { return mu * n; }

// ---------------------------------------------------------------------------------
// One-step reachability from the genetic code alone. The codon-difference and accessibility
// tables that used this no longer appear on Your Sequence. The genetic code is passed in
// rather than imported, so this module stays standalone.

const _BASES = ["A", "C", "G", "T"];

/** Amino acids reachable from a codon in exactly one nucleotide substitution (self excluded). */
export function oneStepAminoAcids(codon, code) {
  const out = new Set();
  for (let i = 0; i < 3; i++) {
    for (const b of _BASES) {
      if (b === codon[i]) continue;
      out.add(code[codon.slice(0, i) + b + codon.slice(i + 1)]);
    }
  }
  out.delete(code[codon]);
  return out;
}
