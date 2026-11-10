// SPDX-License-Identifier: MIT
//
// codon_sensitivity.js -- how a DRM's mutation rate depends on codon optimization.
//
// The question this answers: two viruses carrying the same amino acid at a resistance
// position, but encoding it with different synonymous codons, generate the resistant
// variant at different rates. For each DRM we report that rate as a function of
//
//   1. the wild-type codon at the site itself, varied over every synonymous codon for
//      the wild-type amino acid, and
//   2. the 3' flanking context, varied over only those values reachable by SYNONYMOUS
//      recoding of the NEXT codon -- so the neighbouring amino acid is held fixed and
//      every value corresponds to a protein-identical virus that could actually exist.
//
// Everything is computed from the shipped model bundle via MraEngine, at a fixed aligned
// NL4-3 position, so the position spline and the amino-acid identities are exactly those
// of the published atlas. Only the codon and its 3' neighbour vary.
//
// A NOTE ON SEPARABILITY
//   For a substitution whose changed nucleotide is the codon's THIRD base, that base is
//   also the 3' context of... nothing within the codon, but the NEXT codon's first base
//   is this site's 3' context. The two dimensions we vary are therefore genuinely
//   distinct: dimension 1 changes bases inside the DRM codon, dimension 2 changes the
//   first base of the next codon. They only interact when a path's changed site sits at
//   the codon's 3' edge and reads the next codon's first base as ITS context -- that
//   coupling is detected and the affected combinations are labelled rather than counted
//   as if independent.

import { GENETIC_CODE } from "./mra_engine.js";

const BASES = ["A", "C", "G", "T"];

function synonymousCodons(aa) {
  const out = [];
  for (const c1 of BASES) for (const c2 of BASES) for (const c3 of BASES) {
    const codon = c1 + c2 + c3;
    if (GENETIC_CODE[codon] === aa) out.push(codon);
  }
  return out;
}

// Geometric summary of a set of rates: the fold range (max/min) and the values.
function foldRange(values) {
  const v = values.filter((x) => x > 0);
  if (v.length < 2) return { fold: 1, min: v[0] || 0, max: v[0] || 0 };
  return { fold: Math.max(...v) / Math.min(...v), min: Math.min(...v), max: Math.max(...v) };
}

/**
 * Codon sensitivity of one amino-acid change at one atlas position.
 *
 * @param {MraEngine} engine
 * @param {object} spec  {gene, genePosition, refAa, altAa}
 * @returns {object} per-codon rates and the fold-range decomposition, or {applicable:false}
 */
export function codonSensitivity(engine, { gene, genePosition, refAa, altAa }) {
  const ref = engine.reference;
  const protRec = ref.proteins[gene];
  if (!protRec) return { applicable: false, reason: `unknown gene ${gene}` };
  const idx = genePosition - protRec.first_aa;
  if (idx < 0 || idx >= protRec.codon_nt.length) {
    return { applicable: false, reason: `position ${genePosition} outside ${gene}` };
  }

  const nl43Codon = protRec.codon_seq[idx];
  const codonNt = protRec.codon_nt[idx];            // three NL4-3 reference positions
  if (GENETIC_CODE[nl43Codon] !== refAa) {
    // NL4-3 does not carry the expected wild-type amino acid here. Sensitivity is still
    // computed over synonymous codons for NL4-3's actual residue, but flagged, because
    // "the wild-type codon" the user has in mind is not present.
    refAa = GENETIC_CODE[nl43Codon];
  }

  // The next codon supplies the 3' context. It may lie outside the modelled region or be
  // unmapped (end of protein); if so the flank cannot be varied and we report site-only.
  const nextIdx = idx + 1;
  const hasNext = nextIdx < protRec.codon_nt.length;
  const nextCodonNt = hasNext ? protRec.codon_nt[nextIdx] : null;
  const nextCodon = hasNext ? protRec.codon_seq[nextIdx] : null;
  const nextAa = nextCodon ? GENETIC_CODE[nextCodon] : null;

  const wtCodons = synonymousCodons(refAa);

  // Reconstruct the local nucleotide window as NL4-3 has it, then substitute candidate
  // codons into it. Context for a site is read from this mutable local sequence, exactly
  // as the engine reads it from a submitted sequence.
  //
  // Window: [DRM codon 3 nt][next codon 3 nt]. The 3' base / 2-mer of any site in the DRM
  // codon is read from the bases immediately 3' of it within this window.
  function localWindow(drmCodon, nextCodonSeq) {
    return (drmCodon + (nextCodonSeq || "")).split("");
  }

  // Rate of refAa->altAa for a given DRM codon and a given next codon, summed over paths,
  // reading context from the local window. Mirrors MraEngine.customAtlas path arithmetic.
  function rateForCodons(drmCodon, nextCodonSeq) {
    const win = localWindow(drmCodon, nextCodonSeq);
    let muAa = 0, muLo = 0, muHi = 0, nPaths = 0, minSteps = Infinity;
    let usedFlank = false, flankCoupled = false;
    for (const c1 of BASES) for (const c2 of BASES) for (const c3 of BASES) {
      const mutant = c1 + c2 + c3;
      if (mutant === drmCodon) continue;
      if (GENETIC_CODE[mutant] !== altAa) continue;
      let muPath = 1, loPath = 1, hiPath = 1, nChanges = 0, ok = true;
      for (let slot = 0; slot < 3; slot++) {
        if (drmCodon[slot] === mutant[slot]) continue;
        nChanges++;
        const base = win[slot];
        const alt = mutant[slot];
        const refPos = codonNt[slot];
        // Context: 1 base 3' for non-G->A, 2-mer for G->A, read from the local window.
        const need = base === "G" && alt === "A" ? 2 : 1;
        let ctx = "";
        for (let k = 1; k <= need; k++) {
          const w = slot + k;
          if (w < win.length) ctx += win[w];
        }
        if (slot === 2 && need >= 1) { usedFlank = true; }     // reads the next codon
        if (slot === 2 && base === "G" && alt === "A") { flankCoupled = true; }
        const r = engine.rateFor(refPos, base, alt, ctx);
        if (!r || r.lo95 === null) { ok = false; break; }
        muPath *= r.rate; loPath *= r.lo95; hiPath *= r.hi95;
      }
      if (!ok) continue;
      muAa += muPath; muLo += loPath; muHi += hiPath; nPaths++;
      if (nChanges < minSteps) minSteps = nChanges;
    }
    return { mu_aa: muAa, lo95: muLo, hi95: muHi, n_paths: nPaths,
             min_nt_changes: minSteps === Infinity ? null : minSteps,
             uses_flank: usedFlank, flank_coupled: flankCoupled };
  }

  // --- dimension 1: the DRM codon, with the next codon held at NL4-3 ----------------
  const bySiteCodon = wtCodons.map((codon) => {
    const r = rateForCodons(codon, nextCodon);
    return { codon, ...r, is_nl43: codon === nl43Codon };
  }).filter((r) => r.mu_aa > 0);

  // --- dimension 2: the 3' flank, DRM codon held at NL4-3 ---------------------------
  // Vary the next codon over its synonymous options only, so the neighbour AA is fixed.
  let byFlank = [];
  if (hasNext && nextAa) {
    const flankCodons = synonymousCodons(nextAa);
    byFlank = flankCodons.map((fc) => {
      const r = rateForCodons(nl43Codon, fc);
      return { next_codon: fc, flank_1base: fc[0], flank_2mer: fc.slice(0, 2),
               ...r, is_nl43: fc === nextCodon };
    }).filter((r) => r.mu_aa > 0 && r.uses_flank);
  }

  // CODON OPTIMIZATION vs CATEGORICAL ACCESSIBILITY
  //   Synonymous codons can differ in how many nucleotide changes they need to reach the
  //   target amino acid: some reach it in 2 steps, others only in 3. A 2-step rate and a
  //   3-step rate differ by orders of magnitude, but that is a CATEGORICAL difference in
  //   accessibility, not codon optimization. So the headline fold range is computed only among codons at the MINIMAL step
  //   count (the genuine codon-optimization effect at equal accessibility), and the
  //   step-count spread is reported separately.
  const steps = bySiteCodon.map((r) => r.min_nt_changes).filter((v) => v != null);
  const minSteps = steps.length ? Math.min(...steps) : null;
  const atMin = bySiteCodon.filter((r) => r.min_nt_changes === minSteps);
  const higherStep = bySiteCodon.filter((r) => r.min_nt_changes > minSteps);

  const siteRange = foldRange(atMin.map((r) => r.mu_aa));         // codon optimization
  const fullRange = foldRange(bySiteCodon.map((r) => r.mu_aa));   // incl. step-count spread
  const flankRange = foldRange(byFlank.map((r) => r.mu_aa));
  const nl43Rate = bySiteCodon.find((r) => r.is_nl43);

  // Rank of NL4-3's codon among the synonymous options AT THE SAME STEP COUNT as NL4-3.
  const nl43Steps = nl43Rate ? nl43Rate.min_nt_changes : null;
  const peers = bySiteCodon.filter((r) => r.min_nt_changes === nl43Steps)
                           .sort((a, b) => a.mu_aa - b.mu_aa);
  const nl43Rank = nl43Rate ? peers.findIndex((r) => r.is_nl43) + 1 : null;

  return {
    applicable: true,
    gene, genePosition, refAa, altAa,
    nl43Codon, nl43_matches_expected: GENETIC_CODE[nl43Codon] === refAa,
    nl43_rate: nl43Rate ? nl43Rate.mu_aa : null,
    n_synonymous_codons: bySiteCodon.length,
    min_nt_changes: minSteps,
    n_codons_at_min: atMin.length,
    n_codons_higher_step: higherStep.length,
    // Codon-optimization fold range: among synonymous codons that reach the target at the
    // minimal step count. This is the number that answers "how much does codon choice
    // matter", holding accessibility class fixed.
    site_fold_range: siteRange.fold,
    site_min: siteRange.min, site_max: siteRange.max,
    // Full range including codons that need MORE steps. Dominated by the step-count
    // difference, not by codon optimization; reported but not the headline.
    full_fold_range: fullRange.fold,
    step_count_varies: higherStep.length > 0,
    nl43_step_count: nl43Steps,
    nl43_rank_among_peers: nl43Rank,
    n_peers_at_nl43_step: peers.length,
    flank_varied: byFlank.length > 0,
    flank_fold_range: flankRange.fold,
    flank_coupled: byFlank.some((r) => r.flank_coupled) || bySiteCodon.some((r) => r.flank_coupled),
    next_codon_nl43: nextCodon, next_aa: nextAa,
    by_site_codon: bySiteCodon.sort((a, b) => b.mu_aa - a.mu_aa),
    by_flank: byFlank.sort((a, b) => b.mu_aa - a.mu_aa),
  };
}


// ---------------------------------------------------------------------------------
// GLYCAN SEQUON-LOSS AGGREGATION
// ---------------------------------------------------------------------------------
// Much bNAb escape works by destroying an N-linked glycosylation sequon, N-X-S/T (X != P).
// The glycan is lost if ANY of three things happens: the asparagine changes, the X
// position becomes proline, or the S/T position becomes anything that is neither S nor T.
// So escape supply is the SUMMED mutability of every such change across the three codons,
// not the rate of one particular substitution. This is the quantity that predicts escape,
// and it is typically larger than any single annotated change because several routes reach
// the same end.
//
// DEFINITION OF "DESTROYS THE SEQUON", stated explicitly so it can be argued with:
//   - N codon -> any amino acid other than N          (removes the acceptor asparagine)
//   - X codon -> P                                     (N-P-S/T is not glycosylated)
//   - S/T codon -> any amino acid other than S or T    (removes the hydroxyl)
//   N-X-C is treated as NON-glycosylating, i.e. a change S/T -> C still counts as loss.
//   This is the common convention; sequon scanners that accept N-X-C would differ here,
//   and that choice is documented in docs/DRM_METHODS.md.

const SEQUON_RULES = {
  N: (aa) => aa !== "N",                 // asparagine lost
  X: (aa) => aa === "P",                 // proline blocks glycosylation
  ST: (aa) => aa !== "S" && aa !== "T",  // hydroxyl lost (N-X-C counts as loss)
};

/**
 * Summed rate of losing an N-X-S/T glycan.
 *
 * @param {MraEngine} engine
 * @param {object} sequon  one entry from glycan_sequons.json:
 *                 {nl43_env_position, hxb2_position, motif, codon_nt:[[..],[..],[..]]}
 * @returns {object|null} null if the position is not a sequon in NL4-3
 */
export function sequonLossRate(engine, sequon) {
  const ref = engine.reference;
  const protRec = ref.proteins.Env;
  if (!protRec || !sequon || !sequon.codon_nt || sequon.codon_nt.length !== 3) return null;

  // Confirm NL4-3 actually carries the sequon here, by translating the three codons.
  const idx = sequon.nl43_env_position - protRec.first_aa;
  const codons = [0, 1, 2].map((k) => protRec.codon_seq[idx + k]);
  const aas = codons.map((c) => GENETIC_CODE[c]);
  if (!(aas[0] === "N" && aas[1] !== "P" && (aas[2] === "S" || aas[2] === "T"))) return null;

  const roleOf = ["N", "X", "ST"];
  const perCodon = [];
  let total = 0, totalLo = 0, totalHi = 0, nChanges = 0;
  let largest = null;

  for (let slot = 0; slot < 3; slot++) {
    const role = roleOf[slot];
    const codon = codons[slot];
    const codonNt = sequon.codon_nt[slot];
    const refAa = aas[slot];
    const destroys = SEQUON_RULES[role];

    // Local window for context: this codon followed by the next codon in Env (or, for the
    // S/T codon, the codon after the sequon). Read from NL4-3 as the glycan's own context.
    const nextSeq = protRec.codon_seq[idx + slot + 1] || "";
    const win = (codon + nextSeq).split("");

    // enumerate every amino-acid change at this codon that destroys the sequon, sum mu_aa
    const targetAas = new Set();
    for (const c1 of BASES) for (const c2 of BASES) for (const c3 of BASES) {
      const aa = GENETIC_CODE[c1 + c2 + c3];
      if (aa !== refAa && destroys(aa)) targetAas.add(aa);
    }

    let codonTotal = 0, codonLo = 0, codonHi = 0, codonPaths = 0;
    const byTarget = [];
    for (const altAa of targetAas) {
      let muAa = 0, muLo = 0, muHi = 0, paths = 0, minSteps = Infinity;
      for (const c1 of BASES) for (const c2 of BASES) for (const c3 of BASES) {
        const mutant = c1 + c2 + c3;
        if (mutant === codon || GENETIC_CODE[mutant] !== altAa) continue;
        let muPath = 1, loPath = 1, hiPath = 1, nc = 0, ok = true;
        for (let k = 0; k < 3; k++) {
          if (codon[k] === mutant[k]) continue;
          nc++;
          const base = win[k], alt = mutant[k];
          const need = base === "G" && alt === "A" ? 2 : 1;
          let ctx = "";
          for (let j = 1; j <= need; j++) if (k + j < win.length) ctx += win[k + j];
          const r = engine.rateFor(codonNt[k], base, alt, ctx);
          if (!r || r.lo95 === null) { ok = false; break; }
          muPath *= r.rate; loPath *= r.lo95; hiPath *= r.hi95;
        }
        if (!ok) continue;
        muAa += muPath; muLo += loPath; muHi += hiPath; paths++;
        if (nc < minSteps) minSteps = nc;
      }
      if (muAa > 0) {
        byTarget.push({ alt_aa: altAa, mu_aa: muAa, min_nt_changes: minSteps, n_paths: paths });
        codonTotal += muAa; codonLo += muLo; codonHi += muHi; codonPaths += paths;
        const change = { role, position_role: role, codon, alt_aa: altAa, mu_aa: muAa,
                         min_nt_changes: minSteps };
        if (!largest || muAa > largest.mu_aa) largest = change;
      }
    }
    perCodon.push({ role, codon, ref_aa: refAa, nt: codonNt,
                    mu_loss: codonTotal, lo95: codonLo, hi95: codonHi,
                    n_paths: codonPaths, by_target: byTarget.sort((a, b) => b.mu_aa - a.mu_aa) });
    total += codonTotal; totalLo += codonLo; totalHi += codonHi; nChanges += codonPaths;
  }

  return {
    applicable: true,
    hxb2_position: sequon.hxb2_position,
    nl43_env_position: sequon.nl43_env_position,
    motif: aas.join(""),
    mu_sequon_loss: total, lo95: totalLo, hi95: totalHi,
    n_destroying_changes: nChanges,
    largest_contributor: largest,
    largest_fraction: largest ? largest.mu_aa / total : null,
    per_codon: perCodon,
  };
}


/**
 * Display-ready codon sensitivity for one amino-acid change: every synonymous codon the
 * REFERENCE residue could be written as, the rate this change would have from each, and one
 * number for how much that choice matters.
 *
 * WHICH FOLD RANGE IS THE HEADLINE. `codonSensitivity` reports two. `site_fold_range` holds
 * accessibility fixed and varies only the codon, which isolates the codon-optimisation effect
 * but reads 1.0x whenever a single codon reaches the target in the fewest steps -- at RT L210W
 * that hides the fact that five of the six leucine codons put the tryptophan four to five logs
 * further away. `full_fold_range` spans every synonymous codon regardless of step count, which
 * is the answer to the question actually being asked here: if this virus spelled the wild-type
 * residue differently, how much would the supply of this resistance change move? We report the
 * full range as the score and keep the per-codon step counts visible beside it, so a reader can
 * see when the spread is driven by an extra nucleotide rather than by context alone.
 *
 * @returns {object|null} null when the change has no codon decomposition at this position
 */
export function codonSensitivitySummary(engine, spec) {
  const s = codonSensitivity(engine, spec);
  if (!s || !s.applicable) return null;
  const codons = s.by_site_codon.map((c) => ({
    codon: c.codon, rate: c.mu_aa, steps: c.min_nt_changes, is_nl43: !!c.is_nl43,
  }));
  const rank = codons.findIndex((c) => c.is_nl43) + 1;   // by_site_codon is sorted fast -> slow
  return {
    codons,
    n_codons: codons.length,
    nl43_codon: s.nl43Codon,
    nl43_rate: s.nl43_rate,
    nl43_rank: rank || null,
    fold_range: s.full_fold_range,
    fold_range_same_steps: s.site_fold_range,
    step_count_varies: s.step_count_varies,
    // One cell for the table and one field for the CSV, built from the same list so the two
    // can never disagree. "*" marks the codon NL4-3 actually carries.
    text: codons.map((c) => `${c.codon}${c.is_nl43 ? "*" : ""} ${formatSmall(c.rate)}` +
                            (c.steps > 1 ? ` (${c.steps}nt)` : "")).join(" \u00b7 "),
  };
}

/** Compact rate for in-cell text: two significant figures, exponent only when needed. */
function formatSmall(v) {
  if (!(v > 0)) return "0";
  return v.toExponential(1).replace("e-", "e\u2212");
}

/** "45x", "4.4e4x", or an em dash where a single codon leaves nothing to compare. */
export function formatFold(fold, nCodons) {
  if (nCodons <= 1 || !(fold > 1.0000001)) return "\u2014";
  if (fold < 1000) return fold.toFixed(fold < 10 ? 1 : 0) + "\u00d7";
  return fold.toExponential(1).replace("e+", "e") + "\u00d7";
}
