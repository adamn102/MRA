// SPDX-License-Identifier: MIT
//
// batch_core.js -- aggregation across a batch of submitted sequences.
//
// WHY A MATRIX AND NOT A LIST OF ATLASES. One full-genome atlas is ~62,500 amino-acid rows;
// a hundred of them is ~286 MB, which no browser tab should be asked to hold. The per-site
// supply for the same hundred samples is a Float32 matrix of 8,213 x 100 = 3.13 MB. So the
// runner streams: it asks for one sequence's tables, folds the per-site totals in here, and
// drops the tables before asking for the next. Memory is O(sites x samples), not O(rows x
// samples). The figures above are measured, not assumed (handoff/batch_timing.json).
//
// COVERED IS NOT ZERO. A sample that does not span a site has no opinion about it. Such cells
// stay NaN and are skipped by every statistic -- they must never enter a median as a zero,
// which would drag the batch curve toward the floor wherever amplicons do not overlap.
//
// Every function here is pure and DOM-free so test/batch_cases.mjs can exercise it in node.

/**
 * An empty batch matrix.
 * @param {number} nSamples  rows
 * @param {number} nSites    columns; site s is stored at column (s - firstSite)
 * @param {number} firstSite  reference position of column 0
 */
// The twelve base substitution classes, in a fixed order so a sample's per-type row is always
// laid out identically. Insertions and deletions are deliberately absent: their rate has a
// different denominator and they are excluded from per-site supply throughout the portal.
export const MUT_TYPES = [
  "A->C", "A->G", "A->T", "C->A", "C->G", "C->T",
  "G->A", "G->C", "G->T", "T->A", "T->C", "T->G",
];
const TYPE_INDEX = new Map(MUT_TYPES.map((t, i) => [t, i]));

export function makeBatchMatrix(nSamples, nSites, firstSite) {
  const total = new Float32Array(nSamples * nSites);
  const ga = new Float32Array(nSamples * nSites);
  total.fill(NaN); ga.fill(NaN);
  const codonCode = new Uint8Array(nSamples * nSites);
  codonCode.fill(NO_CODON);
  return { total, ga, nSamples, nSites, firstSite,
           covered: new Uint8Array(nSamples * nSites), codonCode,
           // Per-sample per-type accumulators for the sequence summary table. Summing in
           // Float64 because these add tens of thousands of small rates and a Float32 running
           // total loses the low bits; the count is per type because only the three
           // substitutions away from a site's reference base apply at that site.
           typeSum: new Float64Array(nSamples * MUT_TYPES.length),
           typeCount: new Int32Array(nSamples * MUT_TYPES.length),
           sitesCovered: new Int32Array(nSamples) };
}

// Codons are stored as one byte, not as 800k JS strings: 2 bits per base, so ACG -> 0b000110.
// 255 means "this site is not inside a codon this atlas models". Keeping the whole matrix in
// typed arrays is what lets it be transferred from the worker without a structured-clone copy.
export const NO_CODON = 255;
const B2I = { A: 0, C: 1, G: 2, T: 3 };
const I2B = ["A", "C", "G", "T"];

export function encodeCodon(codon) {
  if (!codon || codon.length !== 3) return NO_CODON;
  const a = B2I[codon[0]], b = B2I[codon[1]], c = B2I[codon[2]];
  if (a === undefined || b === undefined || c === undefined) return NO_CODON;
  return (a << 4) | (b << 2) | c;
}

export function decodeCodon(code) {
  if (code === NO_CODON) return null;
  return I2B[(code >> 4) & 3] + I2B[(code >> 2) & 3] + I2B[code & 3];
}

/**
 * Map every modelled nucleotide position to the codon the sample carries there, from that
 * sample's amino-acid rows. Built per sequence and thrown away with the atlas.
 */
export function codonByPositionFrom(aaRows, proteins) {
  const m = new Map();
  const seen = new Set();
  for (const r of aaRows) {
    const key = r.protein + "|" + r.position;
    if (seen.has(key)) continue;
    seen.add(key);
    const rec = proteins[r.protein];
    if (!rec || !rec.codon_nt) continue;
    const nts = rec.codon_nt[r.position - rec.first_aa];
    if (!nts) continue;
    for (const nt of nts) m.set(nt, r.codon);
  }
  return m;
}

/**
 * Fold one sequence's nucleotide rows into the matrix at row `sampleIndex`.
 * Substitution classes only: indels are excluded from a per-site supply total because their
 * denominator is not the same per-site quantity (the same rule perSiteSupply uses).
 * Returns the number of sites this sample contributed, for the coverage bookkeeping.
 */
export function foldSequenceIntoMatrix(matrix, sampleIndex, ntRows, codonByPosition = null) {
  const { total, ga, nSites, firstSite, covered, codonCode, typeSum, typeCount } = matrix;
  const nTypes = MUT_TYPES.length;
  const acc = new Map();
  for (const r of ntRows) {
    if (!r.mut_type || !r.mut_type.includes("->")) continue;     // skip ins/del
    const v = r.rate != null ? r.rate : r.bio_rate_mean;
    if (v == null || !Number.isFinite(v)) continue;
    const ti = TYPE_INDEX.get(r.mut_type);
    if (ti !== undefined) {
      typeSum[sampleIndex * nTypes + ti] += v;
      typeCount[sampleIndex * nTypes + ti] += 1;
    }
    let cell = acc.get(r.ref_position);
    if (!cell) acc.set(r.ref_position, (cell = { t: 0, g: 0 }));
    cell.t += v;
    if (r.mut_type === "G->A") cell.g += v;
  }
  let n = 0;
  for (const [pos, cell] of acc) {
    const col = pos - firstSite;
    if (col < 0 || col >= nSites) continue;                      // outside the modelled region
    const k = sampleIndex * nSites + col;
    total[k] = cell.t; ga[k] = cell.g; covered[k] = 1;
    if (codonByPosition) codonCode[k] = encodeCodon(codonByPosition.get(pos));
    n++;
  }
  matrix.sitesCovered[sampleIndex] = n;
  return n;
}

/** Median of an already-sorted numeric array. */
function medianSorted(a) {
  if (!a.length) return NaN;
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
}

/** Quantile of an already-sorted array, by linear interpolation. */
function quantileSorted(a, q) {
  if (!a.length) return NaN;
  const i = (a.length - 1) * q, lo = Math.floor(i), hi = Math.ceil(i);
  return lo === hi ? a[lo] : a[lo] + (a[hi] - a[lo]) * (i - lo);
}

/**
 * Per-site statistics across the samples that cover each site.
 *
 * `divergence` is log2(max/min) over the covering samples -- the fold range between the
 * highest- and lowest-supply sample at that exact site. It is 0 when every covering sample
 * agrees, which is the expected result for a batch of identical sequences, and it is reported
 * only where at least `minCovered` samples contribute, because a fold range over two samples
 * is not a batch statement.
 *
 * @returns {Array} one record per site that at least one sample covers
 */
export function batchSiteStats(matrix, { field = "total", minCovered = 2 } = {}) {
  const { nSamples, nSites, firstSite, covered } = matrix;
  const src = matrix[field];
  const out = [];
  for (let col = 0; col < nSites; col++) {
    const vals = [];
    for (let s = 0; s < nSamples; s++) {
      const k = s * nSites + col;
      if (covered[k] && Number.isFinite(src[k])) vals.push(src[k]);
    }
    if (!vals.length) continue;
    vals.sort((a, b) => a - b);
    const min = vals[0], max = vals[vals.length - 1];
    // A fold range needs a non-zero floor; a site whose minimum is 0 is reported as null
    // divergence rather than Infinity, and the zero is visible in the min column.
    const divergence = (vals.length >= minCovered && min > 0) ? Math.log2(max / min) : null;
    out.push({
      ref_position: firstSite + col,
      n_covered: vals.length,
      median: medianSorted(vals),
      min, max,
      q25: quantileSorted(vals, 0.25),
      q75: quantileSorted(vals, 0.75),
      divergence,
    });
  }
  return out;
}

/**
 * The sites where the samples disagree most, by fold range. Sites covered by fewer than
 * `minCovered` samples are excluded: the point of this list is batch-level heterogeneity,
 * and a two-sample fold range is noise, not a batch signal.
 */
export function topDivergentSites(siteStats, { limit = 25, minCovered = 2 } = {}) {
  // Strictly greater than zero, not merely non-null: a batch whose samples all carry the same
  // codon has divergence 0 at every site, and returning those as "most divergent" would make
  // the page's jump-to-divergent control land on an arbitrary position. An empty result is the
  // honest answer, and callers treat it as "no divergence to show".
  return siteStats
    .filter((r) => r.divergence !== null && r.divergence > 0 && r.n_covered >= minCovered)
    .sort((a, b) => b.divergence - a.divergence)
    .slice(0, limit);
}

/**
 * Every sample's value at one site, for the cross-sample comparison plot. Samples that do not
 * cover the site are returned with value null and covered false, so the plot can show them as
 * absent rather than silently omitting them (an omission reads as "no difference").
 */
export function siteAcrossSamples(matrix, refPosition, names, { field = "total" } = {}) {
  const { nSites, firstSite, covered, codonCode } = matrix;
  const src = matrix[field];
  const col = refPosition - firstSite;
  const rows = [];
  if (col < 0 || col >= nSites) return rows;
  for (let s = 0; s < matrix.nSamples; s++) {
    const k = s * nSites + col;
    const isCov = !!covered[k] && Number.isFinite(src[k]);
    rows.push({ sample: names[s] != null ? names[s] : `sample ${s + 1}`, sample_index: s,
                value: isCov ? src[k] : null, covered: isCov,
                codon: decodeCodon(codonCode[k]) });
  }
  return rows;
}

/**
 * Per-sample rollup against the batch: how much of the modelled region the sample covers and
 * how its supply sits relative to the batch median at the sites it shares. `ratioToBatch` is
 * the median of the per-site ratio (sample / batch median), so it is robust to a few extreme
 * sites and is 1 for a sample that tracks the batch exactly.
 */
export function perSampleSummary(matrix, siteStats, names, { field = "total" } = {}) {
  const { nSites, firstSite, covered, nSamples } = matrix;
  const src = matrix[field];
  const medianAt = new Map(siteStats.map((r) => [r.ref_position, r.median]));
  const out = [];
  for (let s = 0; s < nSamples; s++) {
    const ratios = [];
    let nCov = 0;
    for (let col = 0; col < nSites; col++) {
      const k = s * nSites + col;
      if (!covered[k] || !Number.isFinite(src[k])) continue;
      nCov++;
      const med = medianAt.get(firstSite + col);
      if (med > 0) ratios.push(src[k] / med);
    }
    ratios.sort((a, b) => a - b);
    out.push({ sample: names[s] != null ? names[s] : `sample ${s + 1}`, sample_index: s,
               sites_covered: nCov,
               fraction_of_region: nSites ? nCov / nSites : 0,
               ratio_to_batch: ratios.length ? medianSorted(ratios) : null });
  }
  return out;
}

/**
 * Per-sequence rate summary: the bulk rate, each substitution class, and G->A separately.
 *
 * DEFINITIONS, because "the rate for a sequence" is ambiguous and the three quantities below
 * are not interchangeable:
 *
 *  - `bulk_rate` is the mean TOTAL substitution supply per covered site: for each site the
 *    three substitutions away from its reference base are summed, and those per-site totals are
 *    averaged over the sites this sequence covers. It is the expected number of substitutions
 *    per site per replication, so multiplying by the sequence length gives the expected
 *    substitutions per genome per cycle.
 *
 *  - a per-type rate (`A->C` and so on) is the mean rate of that one class over the sites where
 *    it can occur -- i.e. sites whose reference base is the type's source base. It is NOT
 *    averaged over all sites: averaging a G->A rate over non-G sites would divide by a
 *    denominator that includes sites where the substitution is impossible, and the twelve
 *    numbers would no longer be comparable with each other.
 *
 *  - `ga_rate` is the G->A entry restated for prominence, and `ga_share` is the fraction of the
 *    sequence's total supply that is G->A. APOBEC-driven G->A dominates HIV mutational supply,
 *    so this share is the single most useful summary of whether a sequence's profile is
 *    typical, and it is the quantity to check before reading any per-type comparison.
 *
 * Samples that failed to align contribute no sites and come back with nulls, not zeros.
 */
export function perSampleRates(matrix, names) {
  const nTypes = MUT_TYPES.length;
  const { typeSum, typeCount, sitesCovered, nSamples, nSites, covered, total } = matrix;
  const out = [];
  for (let s = 0; s < nSamples; s++) {
    const nCov = sitesCovered[s];
    const row = { sample: names[s] != null ? names[s] : `sample ${s + 1}`, sample_index: s,
                  sites_covered: nCov };
    if (!nCov) {
      row.bulk_rate = null; row.ga_rate = null; row.ga_share = null;
      row.transition_rate = null; row.transversion_rate = null;
      for (const t of MUT_TYPES) row[t] = null;
      out.push(row);
      continue;
    }
    // Bulk rate: mean of the per-site totals over the sites this sample covers.
    let sumTotal = 0;
    for (let col = 0; col < nSites; col++) {
      const k = s * nSites + col;
      if (covered[k] && Number.isFinite(total[k])) sumTotal += total[k];
    }
    row.bulk_rate = sumTotal / nCov;

    let gaSum = 0, tsSum = 0, tsN = 0, tvSum = 0, tvN = 0;
    const TRANSITIONS = new Set(["A->G", "G->A", "C->T", "T->C"]);
    for (let i = 0; i < nTypes; i++) {
      const c = typeCount[s * nTypes + i];
      const mean = c ? typeSum[s * nTypes + i] / c : null;
      row[MUT_TYPES[i]] = mean;
      if (MUT_TYPES[i] === "G->A") gaSum = typeSum[s * nTypes + i];
      if (c) {
        if (TRANSITIONS.has(MUT_TYPES[i])) { tsSum += typeSum[s * nTypes + i]; tsN += c; }
        else { tvSum += typeSum[s * nTypes + i]; tvN += c; }
      }
    }
    row.ga_rate = row["G->A"];
    row.ga_share = sumTotal > 0 ? gaSum / sumTotal : null;
    row.transition_rate = tsN ? tsSum / tsN : null;
    row.transversion_rate = tvN ? tvSum / tvN : null;
    out.push(row);
  }
  return out;
}

/**
 * Group the samples at one site by the codon they carry.
 *
 * WHY CODONS AND NOT SAMPLES. Every sample carrying the same codon in the same 3' context gets
 * exactly the same rate -- the model is identical for all of them. A plot with one row per
 * sample therefore repeats one value dozens of times and, with real sample names, becomes a
 * wall of text whose rows are mostly duplicates. The codon is the unit that actually varies,
 * so it is the unit plotted; the sample names belong in a table where long ones can be read.
 *
 * Groups come back ordered by rate, highest first, so the fastest codon is the top row.
 *
 * @param {Array} rows      siteAcrossSamples() output
 * @param {string|null} refCodon  the reference's codon here, flagged so the user can see
 *                                whether a group matches NL4-3 or departs from it
 */
export function codonGroupsAtSite(rows, refCodon = null) {
  const covered = rows.filter((r) => r.covered);
  const byCodon = new Map();
  for (const r of covered) {
    const key = r.codon || "\u2014";
    let g = byCodon.get(key);
    if (!g) byCodon.set(key, (g = { codon: key, samples: [], rates: [] }));
    g.samples.push(r.sample);
    g.rates.push(r.value);
  }
  const out = [];
  for (const g of byCodon.values()) {
    const rates = g.rates.slice().sort((a, b) => a - b);
    const mid = rates.length % 2 ? rates[rates.length >> 1]
              : 0.5 * (rates[(rates.length >> 1) - 1] + rates[rates.length >> 1]);
    // Rates within a codon group should be identical; they can differ only if the 3' context
    // differs, which happens at a codon's last base when the next base differs too. Carrying
    // the spread lets the caller say so rather than silently showing one of them.
    out.push({ codon: g.codon, n: g.samples.length, samples: g.samples,
               rate: mid, rate_min: rates[0], rate_max: rates[rates.length - 1],
               uniform: rates[rates.length - 1] - rates[0] <= 1e-12 * Math.max(1, rates[0]),
               is_reference: refCodon != null && g.codon === refCodon });
  }
  out.sort((a, b) => b.rate - a.rate);
  return out;
}

// =====================================================================================
// AMINO-ACID LEVEL. Everything above is per nucleotide; everything below is per codon.
//
// WHY THE BATCH WORKS AT CODON LEVEL. The fitted model is identical for every submitted
// sequence -- a batch is not 200 different mutation rates. What differs between samples is
// the codon each carries and its local context, so what a batch measures is how the codon
// usage present in a set of real sequences distributes mutational accessibility over
// sequence space. That question is only meaningful per codon: a nucleotide's rate says
// nothing about which amino acid it reaches, and "mutability" of a site in evolutionary
// terms means the supply of amino-acid change, not of nucleotide change.
//
// ONE STEP, NOT ALL PATHS. Mutability here sums only min_nt_changes === 1 paths. A
// two-nucleotide path is not a per-generation event; summing it into a per-replication rate
// would mix a rate with an accessibility score. Multi-step changes remain in the published
// per-sequence tables, where they are labelled as such.

export const AA_ORDER = "ACDEFGHIKLMNPQRSTVWY*".split("");
const AA_INDEX = new Map(AA_ORDER.map((a, i) => [a, i]));

/**
 * Index the codon positions of the reference once, so every sample folds into the same rows.
 * Key is "protein|position"; the order is genomic, which is also the plotting order.
 */
export function buildAaIndex(proteins) {
  const keys = [], meta = [];
  for (const [protein, rec] of Object.entries(proteins)) {
    if (!rec.codon_nt) continue;
    for (let i = 0; i < rec.codon_nt.length; i++) {
      const nts = rec.codon_nt[i];
      if (!nts) continue;
      keys.push(protein + "|" + (rec.first_aa + i));
      meta.push({ protein, position: rec.first_aa + i, ntStart: nts[0], nts });
    }
  }
  const order = meta.map((m, i) => i).sort((a, b) => meta[a].ntStart - meta[b].ntStart);
  const index = new Map();
  const sorted = order.map((i) => meta[i]);
  order.forEach((i, rank) => index.set(keys[i], rank));
  return { index, meta: sorted, nPositions: sorted.length };
}

/**
 * An empty codon-level matrix.
 *
 * `nonsyn` is the quantity the tracks and the scatter use: the summed rate of every one-step
 * change that alters the residue, excluding stops. Synonymous and stop supply are carried
 * separately rather than folded in, because they answer different questions -- synonymous
 * supply changes the codon without changing the protein, and stop supply is lethal rather
 * than exploratory.
 */
export function makeAaMatrix(nSamples, nPositions) {
  const n = nSamples * nPositions;
  const nonsyn = new Float32Array(n), syn = new Float32Array(n), stop = new Float32Array(n);
  const ga = new Float32Array(n);
  nonsyn.fill(NaN); syn.fill(NaN); stop.fill(NaN); ga.fill(NaN);
  return { nonsyn, syn, stop, ga, nSamples, nPositions,
           covered: new Uint8Array(n),
           codonCode: new Uint8Array(n).fill(NO_CODON),
           aaCode: new Uint8Array(n).fill(255),
           destinations: new Map(),
           sitesCovered: new Int32Array(nSamples) };
}

/** Fold one sample's amino-acid rows into the codon matrix. */
/**
 * Fold one sample's one-step codon paths into the codon matrix.
 *
 * Takes `paths` (from customAtlas) rather than the amino-acid rows, for two reasons. It makes
 * the quantity exact: an amino-acid row's rate is summed over ALL paths reaching that residue,
 * including two- and three-nucleotide ones, so filtering rows by min_nt_changes === 1 would
 * still fold a little multi-step rate into a per-generation number. And a path carries its own
 * mutation type, which is what lets G->A be separated here -- an amino-acid row cannot be
 * classified that way, because several paths with different types can reach the same residue.
 */
export function foldAaIntoMatrix(matrix, sampleIndex, paths, aaIndex) {
  const { nonsyn, syn, stop, ga, nPositions, covered, codonCode, aaCode, destinations } = matrix;
  const acc = new Map();
  for (const r of paths) {
    if (r.n_nt_changes !== 1) continue;              // per-generation events only
    const col = aaIndex.index.get(r.protein + "|" + r.position);
    if (col === undefined) continue;
    let a = acc.get(col);
    if (!a) acc.set(col, (a = { codon: r.codon, ref: r.ref_aa, ns: 0, sy: 0, st: 0, ga: 0,
                                dest: new Float32Array(AA_ORDER.length) }));
    const di = AA_INDEX.get(r.alt_aa);
    if (di !== undefined) a.dest[di] += r.mu_path;
    if (r.alt_aa === "*") a.st += r.mu_path;
    else if (r.alt_aa === r.ref_aa) a.sy += r.mu_path;
    else {
      a.ns += r.mu_path;
      // The G->A class is tracked within the non-synonymous supply, because that is the
      // quantity the tracks plot; separating it answers whether a region's apparent
      // mutability is APOBEC-driven or not.
      if (r.mutation_type_path === "G->A") a.ga += r.mu_path;
    }
  }
  let n = 0;
  for (const [col, a] of acc) {
    const k = sampleIndex * nPositions + col;
    nonsyn[k] = a.ns; syn[k] = a.sy; stop[k] = a.st; ga[k] = a.ga;
    covered[k] = 1;
    codonCode[k] = encodeCodon(a.codon);
    const ai = AA_INDEX.get(a.ref);
    aaCode[k] = ai === undefined ? 255 : ai;
    const dkey = col + "|" + a.codon;
    if (!destinations.has(dkey)) destinations.set(dkey, a.dest);
    n++;
  }
  matrix.sitesCovered[sampleIndex] = n;
  return n;
}

function medianOfSorted(a) {
  if (!a.length) return NaN;
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
}

/**
 * Per-codon statistics across the batch.
 *
 * `entropy` is Shannon entropy of the residue distribution in bits, the conventional measure
 * of site variability in a sequence alignment. `diversity` is 1 - the majority residue's
 * frequency, which is easier to read and is what the scatter plots. Both describe what the
 * submitted sequences ARE; `mutability` describes how readily they change. Keeping the two
 * apart is the whole point of the view.
 */
export function aaPositionStats(matrix, aaIndex, { field = "nonsyn" } = {}) {
  const { nSamples, nPositions, covered, codonCode, aaCode } = matrix;
  // "nonga" is not stored: it is the non-synonymous supply with the G->A part taken out, and
  // deriving it avoids a fourth array that would always equal this subtraction.
  const derived = field === "nonga";
  const src = derived ? matrix.nonsyn : matrix[field];
  const sub = derived ? matrix.ga : null;
  const out = [];
  for (let col = 0; col < nPositions; col++) {
    const vals = [];
    const aaCounts = new Map(), codonCounts = new Map(), codonAa = new Map();
    for (let s = 0; s < nSamples; s++) {
      const k = s * nPositions + col;
      if (!covered[k] || !Number.isFinite(src[k])) continue;
      vals.push(sub ? Math.max(0, src[k] - sub[k]) : src[k]);
      const aa = aaCode[k] === 255 ? "?" : AA_ORDER[aaCode[k]];
      aaCounts.set(aa, (aaCounts.get(aa) || 0) + 1);
      const cd = decodeCodon(codonCode[k]) || "???";
      codonCounts.set(cd, (codonCounts.get(cd) || 0) + 1);
      // The residue this codon encodes, taken from the sample that carried it: avoids
      // importing a genetic code into a module that is deliberately standalone.
      if (!codonAa.has(cd)) codonAa.set(cd, aa);
    }
    if (!vals.length) continue;
    vals.sort((a, b) => a - b);
    const min = vals[0], max = vals[vals.length - 1];
    const n = vals.length;
    const aaSorted = [...aaCounts.entries()].sort((a, b) => b[1] - a[1]);
    const codonSorted = [...codonCounts.entries()].sort((a, b) => b[1] - a[1]);
    let entropy = 0;
    for (const [, c] of aaSorted) { const p = c / n; entropy -= p * Math.log2(p); }
    const m = aaIndex.meta[col];
    out.push({
      key: m.protein + "|" + m.position, position_index: col,
      protein: m.protein, position: m.position, nt_start: m.ntStart,
      n_covered: n,
      median: medianOfSorted(vals), min, max,
      divergence: min > 0 && n >= 2 ? Math.log2(max / min) : null,
      majority_aa: aaSorted[0][0], majority_aa_count: aaSorted[0][1],
      majority_aa_freq: aaSorted[0][1] / n,
      n_aa: aaSorted.length, aa_counts: aaSorted,
      diversity: 1 - aaSorted[0][1] / n,
      entropy,
      majority_codon: codonSorted[0][0], majority_codon_count: codonSorted[0][1],
      majority_codon_freq: codonSorted[0][1] / n,
      n_codons: codonSorted.length, codon_counts: codonSorted, codon_aa: codonAa,
    });
  }
  return out;
}

/**
 * Batch-frequency-weighted supply to each destination amino acid at one position.
 *
 * THIS IS THE QUANTITY A SET OF CIRCULATING SEQUENCES ADDS. A rate computed from NL4-3 is one
 * codon's draw. Weighting each observed codon's rate by how often that codon actually occurs
 * in the submitted set gives the supply as it is distributed in the sampled population --
 * the LANL-weighted rate, computed from the user's own alignment rather than assumed.
 *
 * The weights are the codon frequencies among the samples COVERING this position, so a
 * partially covered position is weighted within its own denominator and is not diluted by
 * samples that cannot speak to it.
 *
 * @returns {Array} per destination residue: weighted mean rate, the range across codons, and
 *                  the fraction of samples whose codon reaches it at all.
 */
export function weightedDestinationSupply(matrix, stat) {
  const col = stat.position_index;
  const n = stat.n_covered;
  const perAa = AA_ORDER.map(() => ({ weighted: 0, min: Infinity, max: -Infinity, nReach: 0 }));
  let used = 0;
  for (const [codon, count] of stat.codon_counts) {
    const d = matrix.destinations.get(col + "|" + codon);
    if (!d) continue;                       // codon never folded (should not happen)
    used += count;
    const w = count / n;
    for (let i = 0; i < AA_ORDER.length; i++) {
      const v = d[i];
      perAa[i].weighted += w * v;
      if (v > 0) {
        perAa[i].nReach += count;
        if (v < perAa[i].min) perAa[i].min = v;
        if (v > perAa[i].max) perAa[i].max = v;
      }
    }
  }
  const out = [];
  for (let i = 0; i < AA_ORDER.length; i++) {
    const a = perAa[i];
    if (a.weighted <= 0 && a.nReach === 0) continue;
    out.push({
      alt_aa: AA_ORDER[i],
      is_stop: AA_ORDER[i] === "*",
      is_synonymous: AA_ORDER[i] === stat.majority_aa,
      weighted_rate: a.weighted,
      rate_min: a.min === Infinity ? null : a.min,
      rate_max: a.max === -Infinity ? null : a.max,
      fold_range: (a.min > 0 && a.max > 0) ? a.max / a.min : null,
      fraction_reaching: used ? a.nReach / used : 0,
    });
  }
  out.sort((x, y) => y.weighted_rate - x.weighted_rate);
  return out;
}

/**
 * How much the mutability at a position depends on which synonymous codon the virus carries.
 *
 * Positions where every sample shares a codon return null: there is no sensitivity to report,
 * and returning 0 would put them in the same bucket as positions genuinely measured to be
 * insensitive. `fraction_on_fastest` is the share of samples sitting on the highest-supply
 * codon -- the quantity that says whether the population is positioned for or against a
 * change being generated.
 */
export function codonSensitivity(matrix, stat) {
  const col = stat.position_index;
  const rows = [];
  for (const [codon, count] of stat.codon_counts) {
    const d = matrix.destinations.get(col + "|" + codon);
    if (!d) continue;
    const codonAa = stat.codon_aa.get(codon);
    let ns = 0;
    for (let i = 0; i < AA_ORDER.length; i++) {
      const a = AA_ORDER[i];
      if (a === "*") continue;
      // Non-synonymous relative to THIS codon's own residue, not the batch majority: a codon
      // encoding a different residue has a different synonymous set.
      if (a === codonAa) continue;
      ns += d[i];
    }
    rows.push({ codon, count, frequency: count / stat.n_covered, aa: codonAa,
                nonsyn_rate: ns });
  }
  if (rows.length < 2) {
    return { n_codons: rows.length, fold_range: null, codons: rows,
             fraction_on_fastest: rows.length ? rows[0].frequency : null, uniform: true };
  }
  rows.sort((a, b) => b.nonsyn_rate - a.nonsyn_rate);
  const hi = rows[0].nonsyn_rate, lo = rows[rows.length - 1].nonsyn_rate;
  return {
    n_codons: rows.length,
    fold_range: lo > 0 ? hi / lo : null,
    fraction_on_fastest: rows[0].frequency,
    codons: rows,
    uniform: false,
  };
}

/**
 * Observed residue diversity against predicted mutability, one point per covered position.
 *
 * WHAT THIS DOES AND DOES NOT LICENSE. Mutability is a prediction about how often a residue
 * change is generated; diversity is an observation of how often one is seen. Where the two
 * disagree, selection is the usual explanation -- a position that mutates readily yet never
 * varies is one where variants are removed. But submitted sequences are related by descent,
 * so diversity also carries phylogeny, sampling and (for consensus sequences) the collapse of
 * within-host variation. The scatter is a description of the submitted set, not an estimate
 * of a selection coefficient.
 */
export function diversityVsMutability(stats, { minCovered = 2 } = {}) {
  return stats
    .filter((s) => s.n_covered >= minCovered && s.median > 0)
    .map((s) => ({
      key: s.key, protein: s.protein, position: s.position,
      mutability: s.median,
      diversity: s.diversity,
      entropy: s.entropy,
      n_covered: s.n_covered,
      majority_aa: s.majority_aa,
      n_aa: s.n_aa,
      n_codons: s.n_codons,
    }));
}

/** Spearman rank correlation, for describing the scatter without assuming linearity. */
export function spearman(xs, ys) {
  const n = xs.length;
  if (n < 3) return null;
  const rank = (v) => {
    const idx = v.map((x, i) => [x, i]).sort((a, b) => a[0] - b[0]);
    const r = new Array(n);
    let i = 0;
    while (i < n) {
      let j = i;
      while (j + 1 < n && idx[j + 1][0] === idx[i][0]) j++;
      const avg = (i + j) / 2 + 1;                      // average rank for ties
      for (let k = i; k <= j; k++) r[idx[k][1]] = avg;
      i = j + 1;
    }
    return r;
  };
  const rx = rank(xs), ry = rank(ys);
  const mx = rx.reduce((a, b) => a + b, 0) / n, my = ry.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx = 0, dy = 0;
  for (let i = 0; i < n; i++) {
    const a = rx[i] - mx, b = ry[i] - my;
    num += a * b; dx += a * a; dy += b * b;
  }
  return (dx > 0 && dy > 0) ? num / Math.sqrt(dx * dy) : null;
}

/**
 * Every sample at one codon position, by name.
 *
 * The figures group by codon because that is what varies; this is where the sample names
 * live, so long identifiers are readable and filterable rather than crowding a plot. A sample
 * that does not cover the position is returned with covered=false rather than omitted, so the
 * table shows the full denominator.
 */
export function samplesAtAaPosition(matrix, stat, names, { field = "nonsyn" } = {}) {
  const { nPositions, covered, codonCode, aaCode } = matrix;
  const derived = field === "nonga";
  const src = derived ? matrix.nonsyn : matrix[field];
  const sub = derived ? matrix.ga : null;
  const col = stat.position_index;
  const out = [];
  for (let s = 0; s < names.length; s++) {
    const k = s * nPositions + col;
    const isCov = !!covered[k] && Number.isFinite(src[k]);
    const v = isCov ? (sub ? Math.max(0, src[k] - sub[k]) : src[k]) : null;
    out.push({
      sample: names[s],
      covered: isCov,
      codon: isCov ? decodeCodon(codonCode[k]) : null,
      aa: isCov && aaCode[k] !== 255 ? AA_ORDER[aaCode[k]] : null,
      value: v,
      ratio: v !== null && stat.median > 0 ? v / stat.median : null,
    });
  }
  return out;
}


/**
 * Summary statistics for one position's mutational output.
 *
 * WHY AN EFFECTIVE NUMBER RATHER THAN RAW ENTROPY. `effective_destinations` is exp(H) over the
 * batch-weighted supply to each non-synonymous destination residue. Entropy in bits is hard to
 * read; its exponential is in natural units -- "this site effectively generates about five
 * distinct residue changes" -- and it combines richness with evenness in one number, so a site
 * reaching twelve residues but sending 98% of its supply to one scores near 1, not near 12.
 * It is the same construction as Wright's effective number of codons, which readers in this
 * field already know.
 *
 * WHY SYNONYMOUS IS EXCLUDED FROM THAT NUMBER BUT REPORTED BESIDE IT. The question exp(H)
 * answers is how concentrated the PROTEIN-changing output is. Folding synonymous routes in
 * would make the statistic mostly track the synonymous fraction, which varies over an order of
 * magnitude between sites (3.8% to 53% across the LANL compendium) and is worth its own
 * column. Stop codons are likewise kept out of the entropy -- a stop is not one residue
 * outcome among twenty -- and reported separately.
 *
 * Returns null where the position has no residue-changing supply at all.
 */
export function siteSupplyStats(matrix, stat) {
  const byDest = new Map();
  let syn = 0, stop = 0, nonsyn = 0;
  for (const [codon, count] of stat.codon_counts) {
    const d = matrix.destinations.get(stat.position_index + "|" + codon);
    if (!d) continue;
    const own = stat.codon_aa.get(codon);
    const freq = count / stat.n_covered;
    for (let i = 0; i < AA_ORDER.length; i++) {
      const v = freq * d[i];
      if (!(v > 0)) continue;
      const alt = AA_ORDER[i];
      if (alt === own) syn += v;
      else if (alt === "*") stop += v;
      else { nonsyn += v; byDest.set(alt, (byDest.get(alt) || 0) + v); }
    }
  }
  const total = syn + stop + nonsyn;
  if (!(total > 0) || !(nonsyn > 0)) return null;

  const ps = [...byDest.values()].map((v) => v / nonsyn);
  let H = 0;
  for (const prob of ps) if (prob > 0) H -= prob * Math.log(prob);
  const k = ps.length;
  let top = null, topShare = 0;
  for (const [a, v] of byDest) if (v / nonsyn > topShare) { top = a; topShare = v / nonsyn; }

  return {
    total_supply: total,
    nonsyn_supply: nonsyn,
    syn_fraction: syn / total,
    stop_fraction: stop / total,
    n_destinations: k,
    effective_destinations: Math.exp(H),
    // Evenness is exp(H)/k: 1 when every reachable residue gets the same supply, towards 0
    // when one dominates. Null at k = 1, where evenness is undefined rather than perfect.
    evenness: k > 1 ? Math.exp(H) / k : null,
    top_destination: top,
    top_share: topShare,
  };
}
