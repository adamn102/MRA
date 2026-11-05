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
// ALIGNMENT MODE: banded LOCAL alignment (Smith-Waterman) with affine gaps.
//   The reference is gag-to-nef (8,618 bp). A submitted sequence may be SHORTER (a
//   clinical pol amplicon) or LONGER (a full LTR-to-LTR genome, ~9.7 kb, whose 5' LTR +
//   leader and 3' LTR lie outside the reference). Local alignment handles both: it finds
//   the best colinear match between the query and the reference and lets everything
//   outside it -- unmatched reference for an amplicon, LTR overhangs for a full genome --
//   fall away for free, reported as unaligned spans. Nothing in the aligned core is
//   quietly dropped; the scope restriction to positions 1-8,213 happens afterwards.
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
// Two caps, because the two batch paths have different binding constraints.
//
// Both are set from measurement rather than guessed (test/align_timing.mjs,
// figures/alignment_timing.png): a near-full-length 8.6 kb genome costs ~27 ms to align
// and ~106 ms to evaluate; a typical 1.5 kb clinical pol amplicon is ~22 ms. Note that
// after the band-edge fix it is MODEL EVALUATION, not alignment, that dominates: a full
// genome enumerates 187,992 codon paths.
//
// MAX_DETAIL_SEQUENCES governs the "Your sequence" path, which RETAINS every atlas so the
// per-sequence tables can be served on demand. That retention is the limit, not the time:
// a full-genome atlas is ~62,500 amino-acid rows, so 100 of them is ~286 MB held in the
// worker and 250 would be ~715 MB. Time at 100 is about 14 s.
export const MAX_DETAIL_SEQUENCES = 100;

// MAX_BATCH_SEQUENCES governs the batch-analysis path, which folds each atlas into the
// per-site matrix and DISCARDS it, so memory is O(sites x samples) and does not grow with
// the tables. Measured at 250: the matrix is 19.6 MB, a batch of 250 pol amplicons takes
// 9.1 s and 250 full-length genomes 35.9 s, with the per-site statistics adding ~30 ms.
// Raising this further is bounded by patience, not memory.
export const MAX_BATCH_SEQUENCES = 250;

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
      // U -> T: an RNA sequence is a legitimate way to write a viral genome, and the model is
      // defined on DNA bases. Mapping here rather than rejecting means a pasted RNA sequence
      // behaves exactly like its DNA spelling; leaving it unmapped silently dropped every U
      // position and returned a plausible-looking atlas at ~78% identity.
      records.push({ name: name || `sequence_${records.length + 1}`,
                     sequence: chunks.join("").toUpperCase()
                                     .replace(/[^A-Z]/g, "").replace(/U/g, "T") });
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
 * Banded local alignment (Smith-Waterman) with affine gaps, seeded on `diagonal`.
 *
 * Returns the single best-scoring colinear alignment between the query and the reference.
 * Both the query and the reference may overhang the alignment at either end for free, so
 * a query longer than the reference (a full genome with LTRs) works exactly as a shorter
 * amplicon does. Returns the traceback as query->reference and reference->query maps.
 */
function bandedAlign(query, reference, diagonal, band) {
  const m = query.length, n = reference.length;
  const width = 2 * band + 1;

  // Reference column range in the band at query position j.
  const loFor = (j) => Math.max(1, j + diagonal - band);
  const hiFor = (j) => Math.min(n, j + diagonal + band);

  const size = width + 2;
  let Mprev = new Float64Array(size).fill(NEG);
  let Xprev = new Float64Array(size).fill(NEG);   // gap in the reference (query insertion)
  let Yprev = new Float64Array(size).fill(NEG);   // gap in the query (reference deletion)
  let Mcur = new Float64Array(size), Xcur = new Float64Array(size), Ycur = new Float64Array(size);

  // Traceback state per cell: 0 = from M, 1 = from X, 2 = from Y, 3 = fresh start (score 0).
  const tbM = new Uint8Array((m + 1) * width);
  const tbX = new Uint8Array((m + 1) * width);
  const tbY = new Uint8Array((m + 1) * width);

  let loPrev = 1;
  let bestScore = 0, bestI = 0, bestJ = 0;   // 0,0,0 = empty alignment (SW floor)
  let edgeAtBest = false;

  for (let j = 1; j <= m; j++) {
    const lo = loFor(j), hi = hiFor(j);
    if (hi < lo) {
      // The band lies entirely outside the reference: this query base is in an overhang
      // (an LTR, say). No cells to compute; the previous row becomes empty so the next
      // in-band row starts fresh.
      Mprev.fill(NEG); Xprev.fill(NEG); Yprev.fill(NEG);
      loPrev = Math.max(1, j + diagonal - band);
      continue;
    }
    Mcur.fill(NEG); Xcur.fill(NEG); Ycur.fill(NEG);
    const qb = query[j - 1];
    for (let i = lo; i <= hi; i++) {
      const c = i - lo + 1;                 // column in the current row
      const cPrevDiag = i - 1 - loPrev + 1;  // (i-1, j-1) in the previous row
      const cPrevSame = i - loPrev + 1;      // (i,   j-1) in the previous row
      const sub = qb === reference[i - 1] ? SCORES.match : SCORES.mismatch;

      // M: extend a previous M/X/Y diagonally, OR start fresh at 0. The fresh-start option
      // is what frees BOTH leading overhangs (query and reference) at no cost, and the
      // SW floor below frees the trailing ones.
      let bestM = 0, fromM = 3;
      if (cPrevDiag >= 1 && cPrevDiag < size) {
        const cands = [Mprev[cPrevDiag], Xprev[cPrevDiag], Yprev[cPrevDiag]];
        for (let st = 0; st < 3; st++) if (cands[st] > bestM) { bestM = cands[st]; fromM = st; }
      }
      let mval = bestM + sub;
      if (mval < 0) { mval = 0; fromM = 3; }   // Smith-Waterman floor
      Mcur[c] = mval;
      tbM[j * width + (c - 1)] = fromM;

      // X: gap in the reference (a base present in the query only)
      let bestX = NEG, fromX = 3;
      if (cPrevSame >= 1 && cPrevSame < size) {
        const open = Mprev[cPrevSame] + SCORES.gapOpen + SCORES.gapExtend;
        const ext = Xprev[cPrevSame] + SCORES.gapExtend;
        if (open >= ext) { bestX = open; fromX = 0; } else { bestX = ext; fromX = 1; }
      }
      Xcur[c] = bestX;
      tbX[j * width + (c - 1)] = fromX;

      // Y: gap in the query (a reference position the query skips)
      let bestY = NEG, fromY = 3;
      if (c - 1 >= 1) {
        const open = Mcur[c - 1] + SCORES.gapOpen + SCORES.gapExtend;
        const ext = Ycur[c - 1] + SCORES.gapExtend;
        if (open >= ext) { bestY = open; fromY = 0; } else { bestY = ext; fromY = 2; }
      }
      Ycur[c] = bestY;
      tbY[j * width + (c - 1)] = fromY;

      if (Mcur[c] > bestScore) {
        bestScore = Mcur[c]; bestI = i; bestJ = j;
        edgeAtBest = (i === lo && lo > 1) || (i === hi && hi < n);
      }
    }
    [Mprev, Mcur] = [Mcur, Mprev];
    [Xprev, Xcur] = [Xcur, Xprev];
    [Yprev, Ycur] = [Ycur, Yprev];
    loPrev = lo;
  }

  let touchedEdge = edgeAtBest;

  // ---- traceback from the global best cell until a fresh start -----------------
  const queryToRef = new Int32Array(m);
  const refToQuery = new Int32Array(n + 1);
  let matches = 0, alignedColumns = 0;
  let i = bestI, j = bestJ, state = 0;

  while (j >= 1 && i >= 1 && bestScore > 0) {
    const loj = loFor(j);
    const c = i - loj + 1;
    if (c < 1 || c > width) { touchedEdge = true; break; }
    const idx = j * width + (c - 1);
    if (state === 0) {
      queryToRef[j - 1] = i;
      refToQuery[i] = j;
      alignedColumns++;
      if (query[j - 1] === reference[i - 1]) matches++;
      const from = tbM[idx];
      if (from === 3) break;                 // fresh start: alignment begins here
      i -= 1; j -= 1; state = from;
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
    if (j >= 1) {
      const jj = j;
      if ((i === loFor(jj) && loFor(jj) > 1) || (i === hiFor(jj) && hiFor(jj) < n)) touchedEdge = true;
    }
  }

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
