// SPDX-License-Identifier: MIT
//
// align.js -- place a submitted HIV-1 sequence on the NL4-3 reference.
//
// WHY NOT PLAIN NEEDLEMAN-WUNSCH
//   An 8.6 kb query against an 8.6 kb reference is ~74 million cells. That stalls the
//   main thread for seconds and allocates hundreds of megabytes. Alignment is the only
//   real performance constraint in this portal -- evaluating the model afterwards is
//   about 25,000 lookups and some codon arithmetic, i.e. milliseconds.
//
// WHAT THIS DOES INSTEAD
//   1. Seed. Index the reference's 12-mers, look up the query's 12-mers, and vote on
//      the diagonal (reference position - query position). HIV-1 sequences differ from
//      NL4-3 by a few percent plus indels, so the vote is sharp.
//   2. Extend. Run a banded affine-gap (Gotoh) alignment around the winning diagonal.
//      The band is widened and the alignment re-run if the optimal path touches a band
//      edge, so a sequence with a large indel is not silently clipped.
//
// ALIGNMENT MODE: global in the query, free at both ends of the reference.
//   A clinical pol fragment is the common case. It must align without being charged for
//   the ~7 kb of reference it does not cover, so reference prefix and suffix gaps are
//   free. The query, by contrast, is placed in full: every query base either maps to a
//   reference position or is reported in an unaligned span. Nothing is quietly dropped.
//
// Both orientations are tried and the better-scoring one is returned, because users do
// paste reverse-complemented sequences.

// Scoring. Nucleotide-level, tuned for a few-percent-divergent viral sequence: a
// mismatch must be cheaper than a pair of gaps, and opening a gap must be much more
// expensive than extending one so that real indels come out as single events.
export const SCORES = { match: 2, mismatch: -3, gapOpen: -6, gapExtend: -1 };

export const SEED_K = 12;              // k-mer length for diagonal voting
export const DEFAULT_BAND = 64;        // half-width, in diagonals
export const MAX_BAND = 1024;          // give up rather than fall back to full NW
// Batch cap, set from measurement rather than guessed (test/align_timing.mjs,
// figures/alignment_timing.png): a near-full-length 8.6 kb genome costs ~27 ms to align
// and ~106 ms to evaluate, so 100 of them is about 14 s in a worker with progress
// reporting. A typical 1.5 kb clinical pol amplicon is ~22 ms, so 100 of those is ~2 s.
// Note that after the band-edge fix it is MODEL EVALUATION, not alignment, that
// dominates: a full genome enumerates 187,992 codon paths.
export const MAX_BATCH_SEQUENCES = 100;

const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A", N: "N", "-": "-" };

export function reverseComplement(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]] || "N";
  return out;
}

/** Parse FASTA. A bare sequence with no header is accepted as one unnamed record. */
export function parseFasta(text) {
  const records = [];
  let name = null, chunks = [];
  const flush = () => {
    if (chunks.length) {
      records.push({ name: name || `sequence_${records.length + 1}`,
                     sequence: chunks.join("").toUpperCase().replace(/[^A-Z]/g, "") });
    }
    chunks = [];
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) { flush(); name = line.slice(1).trim(); }
    else if (line.trim()) chunks.push(line.trim());
  }
  flush();
  return records;
}

/** Index of reference k-mers -> list of 1-based start positions.
 *  Exported so a caller aligning a batch builds it once rather than per sequence. */
export function buildSeedIndex(reference, k = SEED_K) {
  const index = new Map();
  for (let i = 0; i + k <= reference.length; i++) {
    const kmer = reference.slice(i, i + k);
    let list = index.get(kmer);
    if (!list) index.set(kmer, (list = []));
    if (list.length < 8) list.push(i + 1);   // cap: a repeated k-mer should not dominate
  }
  return index;
}

/** Vote on the dominant diagonal (reference position - query position). */
function voteDiagonal(query, reference, index, k) {
  const votes = new Map();
  let total = 0;
  for (let j = 0; j + k <= query.length; j += 3) {      // every 3rd k-mer is plenty
    const hits = index.get(query.slice(j, j + k));
    if (!hits) continue;
    for (const refPos of hits) {
      const d = refPos - (j + 1);
      votes.set(d, (votes.get(d) || 0) + 1);
      total++;
    }
  }
  if (!total) return { diagonal: null, support: 0, spread: 0 };
  // Coalesce votes into a window so that an indel-split query still gives one peak.
  let best = null, bestCount = -1;
  for (const d of votes.keys()) {
    let count = 0;
    for (let dd = d - 32; dd <= d + 32; dd++) count += votes.get(dd) || 0;
    if (count > bestCount) { bestCount = count; best = d; }
  }
  const ds = [...votes.keys()];
  return { diagonal: best, support: bestCount, spread: Math.max(...ds) - Math.min(...ds) };
}

const NEG = -1e9;

/**
 * Banded affine-gap alignment, global in the query and free at both reference ends.
 * Returns the traceback as query->reference and reference->query maps.
 */
function bandedAlign(query, reference, diagonal, band) {
  const m = query.length, n = reference.length;
  const width = 2 * band + 1;

  // Column index within the band for reference position i at query position j.
  const loFor = (j) => Math.max(1, j + diagonal - band);
  const hiFor = (j) => Math.min(n, j + diagonal + band);

  // Three score rows (previous, current) for match/insert-in-query/delete-in-query.
  const size = width + 2;
  let Mprev = new Float64Array(size).fill(NEG);
  let Xprev = new Float64Array(size).fill(NEG);   // gap in the reference (query insertion)
  let Yprev = new Float64Array(size).fill(NEG);   // gap in the query (reference deletion)
  let Mcur = new Float64Array(size), Xcur = new Float64Array(size), Ycur = new Float64Array(size);

  // Traceback: 2 bits of state per cell. 0 = from M, 1 = from X, 2 = from Y, 3 = start.
  const tbM = new Uint8Array((m + 1) * width);
  const tbX = new Uint8Array((m + 1) * width);
  const tbY = new Uint8Array((m + 1) * width);

  // Row j = 0: a free reference prefix, so every reference position starts at score 0.
  {
    const lo = loFor(0), hi = hiFor(0);
    for (let i = lo; i <= hi; i++) Mprev[i - lo + 1] = 0;
    Mprev.loBase = lo;
  }
  let loPrev = loFor(0);

  let touchedEdge = false;

  for (let j = 1; j <= m; j++) {
    const lo = loFor(j), hi = hiFor(j);
    Mcur.fill(NEG); Xcur.fill(NEG); Ycur.fill(NEG);
    const qb = query[j - 1];
    for (let i = lo; i <= hi; i++) {
      const c = i - lo + 1;                 // column in the current row
      const cPrevDiag = i - 1 - loPrev + 1;  // (i-1, j-1) in the previous row
      const cPrevSame = i - loPrev + 1;      // (i,   j-1) in the previous row

      // M: query base j paired with reference base i
      const sub = qb === reference[i - 1] ? SCORES.match : SCORES.mismatch;
      let bestM = NEG, fromM = 3;
      if (cPrevDiag >= 1 && cPrevDiag < size) {
        const cands = [Mprev[cPrevDiag], Xprev[cPrevDiag], Yprev[cPrevDiag]];
        for (let s = 0; s < 3; s++) if (cands[s] > bestM) { bestM = cands[s]; fromM = s; }
      }
      if (j === 1) {                         // free reference prefix feeds M directly
        if (0 > bestM) { bestM = 0; fromM = 3; }
      }
      Mcur[c] = bestM <= NEG / 2 ? NEG : bestM + sub;
      tbM[j * width + (c - 1)] = fromM;

      // X: gap in the reference, i.e. a base present in the query only
      let bestX = NEG, fromX = 3;
      if (cPrevSame >= 1 && cPrevSame < size) {
        const open = Mprev[cPrevSame] + SCORES.gapOpen + SCORES.gapExtend;
        const ext = Xprev[cPrevSame] + SCORES.gapExtend;
        if (open >= ext) { bestX = open; fromX = 0; } else { bestX = ext; fromX = 1; }
      }
      Xcur[c] = bestX;
      tbX[j * width + (c - 1)] = fromX;

      // Y: gap in the query, i.e. a reference position the query skips
      let bestY = NEG, fromY = 3;
      if (c - 1 >= 1) {
        const open = Mcur[c - 1] + SCORES.gapOpen + SCORES.gapExtend;
        const ext = Ycur[c - 1] + SCORES.gapExtend;
        if (open >= ext) { bestY = open; fromY = 0; } else { bestY = ext; fromY = 2; }
      }
      Ycur[c] = bestY;
      tbY[j * width + (c - 1)] = fromY;
    }
    [Mprev, Mcur] = [Mcur, Mprev];
    [Xprev, Xcur] = [Xcur, Xprev];
    [Yprev, Ycur] = [Ycur, Yprev];
    loPrev = lo;
  }

  // Free reference suffix: the best cell anywhere along the last query row wins.
  const lo = loFor(m), hi = hiFor(m);
  let bestScore = NEG, bestI = -1, bestState = 0;
  for (let i = lo; i <= hi; i++) {
    const c = i - lo + 1;
    const states = [Mprev[c], Xprev[c], Yprev[c]];
    for (let s = 0; s < 3; s++) {
      if (states[s] > bestScore) { bestScore = states[s]; bestI = i; bestState = s; }
    }
  }
  // Only a band edge that is NOT also a reference edge means the band was too narrow.
  // A query covering the whole reference legitimately ends at reference position 1 or n,
  // and treating that as a clipped path sent the aligner through needless re-runs at
  // ever wider bands -- 347 ms instead of 16 ms for a full-length genome.
  if ((bestI === lo && lo > 1) || (bestI === hi && hi < n)) touchedEdge = true;

  // ---- traceback --------------------------------------------------------------
  const queryToRef = new Int32Array(m);
  const refToQuery = new Int32Array(n + 1);
  let matches = 0, alignedColumns = 0;
  let i = bestI, j = m, state = bestState;

  while (j >= 1 && i >= 1) {
    const loj = loFor(j);
    const c = i - loj + 1;
    if (c < 1 || c > 2 * band + 1) { touchedEdge = true; break; }
    const idx = j * width + (c - 1);
    if (state === 0) {
      queryToRef[j - 1] = i;
      refToQuery[i] = j;
      alignedColumns++;
      if (query[j - 1] === reference[i - 1]) matches++;
      const from = tbM[idx];
      i -= 1; j -= 1; state = from === 3 ? 0 : from;
      if (from === 3) break;
    } else if (state === 1) {
      queryToRef[j - 1] = 0;                 // query base with no reference position
      const from = tbX[idx];
      j -= 1; state = from === 3 ? 0 : from;
      if (from === 3) break;
    } else {
      const from = tbY[idx];
      i -= 1; state = from === 3 ? 0 : from;
      if (from === 3) break;
    }
    const jj = Math.max(j, 1);
    if ((i === loFor(jj) && loFor(jj) > 1) || (i === hiFor(jj) && hiFor(jj) < n)) touchedEdge = true;
  }
  // Any query bases left before the traceback terminated are unaligned by definition.
  for (let jj = 1; jj <= j; jj++) queryToRef[jj - 1] = 0;

  return { queryToRef, refToQuery, score: bestScore, matches, alignedColumns, touchedEdge };
}

function summarise(query, reference, aln, extra) {
  const spans = [];
  let start = null;
  for (let j = 1; j <= query.length; j++) {
    const mapped = aln.queryToRef[j - 1] !== 0;
    if (!mapped && start === null) start = j;
    if (mapped && start !== null) { spans.push({ queryStart: start, queryEnd: j - 1 }); start = null; }
  }
  if (start !== null) spans.push({ queryStart: start, queryEnd: query.length });

  let refStart = 0, refEnd = 0;
  for (let p = 1; p <= reference.length; p++) {
    if (aln.refToQuery[p]) { if (!refStart) refStart = p; refEnd = p; }
  }
  return {
    queryToRef: aln.queryToRef,
    refToQuery: aln.refToQuery,
    score: aln.score,
    matches: aln.matches,
    alignedColumns: aln.alignedColumns,
    alignedCount: aln.alignedColumns,
    identity: aln.alignedColumns ? aln.matches / aln.alignedColumns : 0,
    coverageOfQuery: aln.alignedColumns / query.length,
    unalignedSpans: spans,
    refStart, refEnd,
    ...extra,
  };
}

/**
 * Align a submitted sequence to the reference.
 *
 * @returns {object} see the API block at the top of test/roundtrip_nl43.mjs, plus
 *   `orientation` ("forward" | "reverse_complement"), `band`, `seedSupport`, and
 *   `unalignable` with a `reason` when no seed support exists at all.
 */
export function alignToReference(querySeq, referenceSeq, options = {}) {
  const reference = referenceSeq.toUpperCase();
  const band0 = options.band || DEFAULT_BAND;
  const index = options.seedIndex || buildSeedIndex(reference, SEED_K);

  const candidates = [{ orientation: "forward", seq: querySeq.toUpperCase() }];
  if (options.tryReverseComplement !== false) {
    candidates.push({ orientation: "reverse_complement",
                      seq: reverseComplement(querySeq.toUpperCase()) });
  }

  // Vote first for both orientations, then align only the plausible one(s). A banded
  // alignment is the expensive part, and a sequence pasted the right way round has
  // overwhelming forward seed support, so the reverse pass is usually skippable.
  const voted = candidates.map((c) => ({ ...c, vote: voteDiagonal(c.seq, reference, index, SEED_K) }))
                          .filter((c) => c.vote.diagonal !== null);
  const topSupport = Math.max(0, ...voted.map((c) => c.vote.support));
  const worth = voted.filter((c) => c.vote.support * 5 >= topSupport);

  let best = null;
  for (const cand of worth) {
    const vote = cand.vote;
    let band = band0, aln = null;
    // Widen and re-run if the optimal path touched a band edge: a large indel needs a
    // wider band, and silently clipping it would shift the frame downstream.
    for (;;) {
      aln = bandedAlign(cand.seq, reference, vote.diagonal, band);
      if (!aln.touchedEdge || band >= MAX_BAND) break;
      band = Math.min(band * 4, MAX_BAND);
    }
    const summary = summarise(cand.seq, reference, aln,
      { orientation: cand.orientation, band, seedSupport: vote.support,
        seedDiagonal: vote.diagonal, bandTouchedEdge: aln.touchedEdge,
        query: cand.seq });
    if (!best || summary.score > best.score) best = summary;
  }

  if (!best) {
    return {
      unalignable: true,
      reason: `No ${SEED_K}-mer of the submitted sequence matches the NL4-3 reference in ` +
              "either orientation. The sequence is either not HIV-1, too short, or too " +
              "divergent to place on this reference.",
      queryToRef: new Int32Array(querySeq.length),
      refToQuery: new Int32Array(reference.length + 1),
      identity: 0, alignedCount: 0, unalignedSpans: [{ queryStart: 1, queryEnd: querySeq.length }],
      refStart: 0, refEnd: 0,
    };
  }
  return best;
}
