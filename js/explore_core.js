// SPDX-License-Identifier: MIT
//
// explore_core.js -- the DOM-free logic behind the genome track on the browse page.
//
// Named for the "Explore the genome" tab it was written for; that tab has since been merged
// into browse.html, and this module now backs the track there. Kept separate from the page
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

/**
 * De novo emergence of one amino-acid change, as a function of how much the virus replicates.
 *
 * THE MODEL. One replication cycle produces N new genomes from wild-type templates. Each is an
 * independent trial that yields this specific amino-acid change with probability mu, so the
 * number of new mutant genomes is Binomial(N, mu) and
 *
 *     P(at least one generated) = 1 - (1 - mu)^N        E[copies] = N * mu
 *
 * Against log10(N) the probability is a sigmoid: flat while N << 1/mu, rising steeply around
 * N* = 1/mu, saturating at 1. N* is where ONE copy is expected, which is a 63% chance of
 * seeing at least one -- not a guarantee.
 *
 * WHAT N IS. A count of replication events per cycle, not a fraction and not a viral load:
 * N = (population size) x (fraction replicating this cycle). It is dimensionless only in the
 * sense that it counts genomes; it is not bounded by 1.
 *
 * WHAT THIS IS NOT. Generation, not establishment. A variant that is generated may still be
 * lost to drift, purged by selection, or never sampled; this function speaks only to
 * mutational supply. It is also per cycle, so it understates cumulative supply over an
 * infection. It applies to single-nucleotide changes, where one replication event suffices;
 * a multi-step change is not a single-cycle event and must not be read off this curve.
 *
 * @param {object} r  an atlas row with mu_aa and, optionally, its 95% bounds
 */
export function emergenceStats(r) {
  const mu = r.mu_aa;
  if (!(mu > 0)) return { mu: 0, nStar: null, nLo: null, nHi: null, pAtNStar: null };
  return {
    mu,
    nStar: 1 / mu,
    nLo: r.mu_aa_hi95 > 0 ? 1 / r.mu_aa_hi95 : null,   // upper rate -> lower population
    nHi: r.mu_aa_lo95 > 0 ? 1 / r.mu_aa_lo95 : null,
    pAtNStar: 1 - Math.exp(-1),
  };
}

/** P(at least one copy generated) in N replication events. Exact binomial, not Poisson. */
export function emergenceProbability(mu, n) {
  if (!(mu > 0) || !(n > 0)) return 0;
  return 1 - Math.pow(1 - mu, n);
}

/** Expected de novo copies of a change in N replication events, one cycle. */
export function expectedCopies(mu, n) { return mu * n; }
