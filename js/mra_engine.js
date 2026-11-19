// SPDX-License-Identifier: MIT
//
// mra_engine.js -- the custom-atlas engine.
//
// Plain ES module, no DOM, no I/O, no dependencies. Everything comes in as parsed
// objects so the same code runs in the browser and under node in the regression test.
//
// WHAT THIS DOES AND DOES NOT CLAIM
//   This is not a measurement of mutation rates in a user's virus. It is the MRA model,
//   fitted on NL4-3, applied to a different sequence. Three components contribute, and
//   every row records which ones it rests on:
//
//     codon    -- the submitted sequence's own codon. Two viruses with the same amino
//                 acid at a position but different synonymous codons reach different
//                 sets of amino acids at different rates. This is a real, large effect
//                 determined entirely by the submitted sequence, and it is the main
//                 reason the custom atlas is worth computing.
//     context  -- the 3' base (1 nt for substitutions other than G->A, the 2-mer for
//                 G->A) read from the submitted sequence, not from NL4-3.
//     position -- the fitted position spline, which transfers only where the input
//                 aligns to NL4-3. Where it does not, the genome mean of that spline
//                 class is substituted and the row is flagged.
//
// THE MODEL (confirm against model_metadata.json / postfit_model_checklist.txt)
//   non-G->A substitutions, 11 types:
//       log rate = beta_sub_type_bio[type]
//                + beta_sub_context_delta[type|3' base]        (zero for pooled cells)
//                + spline_shared_sub(aligned position)
//   G->A, modelled apart because APOBEC3-class editing dominates it:
//       log rate = intercept_bio_ga[pooled 3' 2-mer] + spline_GtoA_<class>(position)
//       the 2-mer is read at the G site; CC, CG, TC and TT pool to the level OTHER;
//       the raw 2-mers AA and GA each have a dedicated spline, the rest share SHARED.
//   indels: intercept + type + a dedicated spline, with NO nucleotide context, so they
//       carry no sequence-specific information and are excluded from custom output.
//
// AMINO-ACID MUTABILITY
//   For each codon, all 63 alternative codons are enumerated. A path's rate is the
//   PRODUCT of the per-site rates at the changed positions, and mu_aa is the SUM of
//   path rates over every path reaching that amino acid -- all paths, 1-, 2- and
//   3-step, not only the minimal ones. Both conventions are taken from
//   comprehensive_mutation_supply_analysis.py and aa_level_supply_ci_methodology.txt,
//   which are the authority, and both are verified by test/roundtrip_nl43.mjs.
//   Intervals are the same arithmetic on the marginal bounds: they approximate joint
//   uncertainty and are useful for ranking and sensitivity, not as exact simultaneous
//   credible intervals.

export const GENETIC_CODE = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", CTT: "L", CTC: "L", CTA: "L", CTG: "L",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M", GTT: "V", GTC: "V", GTA: "V", GTG: "V",
  TCT: "S", TCC: "S", TCA: "S", TCG: "S", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  ACT: "T", ACC: "T", ACA: "T", ACG: "T", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", CAT: "H", CAC: "H", CAA: "Q", CAG: "Q",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K", GAT: "D", GAC: "D", GAA: "E", GAG: "E",
  TGT: "C", TGC: "C", TGA: "*", TGG: "W", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
  AGT: "S", AGC: "S", AGA: "R", AGG: "R", GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

const BASES = ["A", "C", "G", "T"];

// Sub-protein abbreviations users actually type and recognise. Keys are the `protein`
// column of the workbook `alignment` sheet.
export const SUBPROTEIN_ABBREV = {
  Protease: "PR", RT: "RT", RNase: "RNH", Integrase: "IN",
  p17: "p17", p24: "p24", p2: "p2", p7: "p7", p1: "p1", p6: "p6",
};

// The gene name used in the atlas tables for each `gene` value on the alignment sheet.
const GENE_TO_PROTEIN = { gag: "Gag", pol: "Pol", env: "Env", nef: "Nef",
                          Vif: "Vif", Vpr: "Vpr", Vpu: "Vpu", Tat: "Tat", Rev: "Rev" };

function isACGT(ch) {
  return ch === "A" || ch === "C" || ch === "G" || ch === "T";
}

export function mutationType(fromBase, toBase) {
  return fromBase + "->" + toBase;
}

export class MraEngine {
  /**
   * @param {object} assets
   * @param {object} assets.bundle     data/<version>/mra_web_model.json
   * @param {object} assets.reference  data/<version>/nl43_reference.json
   * @param {object} assets.intervals  data/<version>/interval_ratios.json
   * @param {object} assets.hxb2       data/<version>/hxb2_align.json
   */
  constructor({ bundle, reference, intervals, hxb2 }) {
    if (!bundle || !reference) throw new Error("MraEngine needs at least bundle and reference");
    this.bundle = bundle;
    this.reference = reference;
    this.intervals = intervals || null;
    this.hxb2 = hxb2 || null;

    this.firstPosition = bundle.modelled_region.first_position;
    this.lastPosition = bundle.modelled_region.last_position;

    // Pooled context cells are multiplied by zero in the fit, so their sampled delta is
    // prior-driven rather than zero. Applying it would be quietly wrong by up to 0.55%.
    this.pooledCells = new Set(bundle.sub_context_pooled_cells || []);
    this.gaPooled2mers = new Set(bundle.ga_pooled_2mers || []);
    this.splineByPosition = bundle.spline_by_position;
    this.splineFirst = bundle.spline_by_position.first_position;
    this.fallback = bundle.unaligned_fallback;

    this._buildIntervalIndex();
    this._buildCoordinateIndex();
  }

  // -- interval lookup -----------------------------------------------------------
  // Published per-site intervals where they exist, a (type, context) mean where they do
  // not. See the docstring of build_interval_ratios in scripts/build_web_assets.py.
  _buildIntervalIndex() {
    this.intervalByPosition = new Map();
    if (!this.intervals) return;
    for (const [mutType, rec] of Object.entries(this.intervals.by_position || {})) {
      const { positions, lo_ratio, hi_ratio } = rec;
      for (let i = 0; i < positions.length; i++) {
        this.intervalByPosition.set(mutType + ":" + positions[i], [lo_ratio[i], hi_ratio[i]]);
      }
    }
  }

  // -- coordinates ---------------------------------------------------------------
  _buildCoordinateIndex() {
    // NL4-3 amino-acid position -> HXB2 amino-acid position, for Env and Vpu only.
    // Every other protein maps identity. The mapping is position-by-position because
    // of indels, so a constant offset would be wrong.
    this.nl43ToHxb2 = {};
    if (this.hxb2 && this.hxb2.hxb2_to_nl43) {
      for (const [protein, table] of Object.entries(this.hxb2.hxb2_to_nl43)) {
        const inverse = new Map();
        for (const [hxb2Pos, nl43Pos] of Object.entries(table)) {
          if (nl43Pos !== null) inverse.set(Number(nl43Pos), Number(hxb2Pos));
        }
        this.nl43ToHxb2[protein] = inverse;
      }
    }
    // Sub-protein spans within a gene, in amino-acid units.
    this.subProteins = {};
    const offsets = (this.hxb2 && this.hxb2.protein_offsets) || [];
    for (const rec of offsets) {
      const protein = GENE_TO_PROTEIN[rec.gene];
      if (!protein || rec.offset === null) continue;      // Vpu's offset is N/A by design
      if (protein === rec.protein) continue;              // the gene product itself
      (this.subProteins[protein] = this.subProteins[protein] || []).push({
        name: rec.protein,
        abbrev: SUBPROTEIN_ABBREV[rec.protein] || rec.protein,
        firstGeneAa: rec.offset + 1,
        lastGeneAa: rec.offset + rec.length,
        offset: rec.offset,
      });
    }
    for (const list of Object.values(this.subProteins)) list.sort((a, b) => a.firstGeneAa - b.firstGeneAa);
  }

  /** HXB2 amino-acid position for a protein position in NL4-3 numbering. */
  hxb2Position(protein, nl43AaPosition) {
    const table = this.nl43ToHxb2[protein];
    if (!table) return nl43AaPosition;                 // identity for everything but Env, Vpu
    const hit = table.get(nl43AaPosition);
    return hit === undefined ? null : hit;             // null = an NL4-3 insertion vs HXB2
  }
  /**
   * The selectable protein/subunit list. Genes run 5' to 3' — Gag, Pol, Vif, Vpr, Tat,
   * Rev, Vpu, Env, Nef — and each gene is followed by its subunits in that same direction.
   * `offset` converts the
   * subunit's own numbering to gene numbering (RT 103 -> Pol 253).
   *
   * @param {Map<string,{min:number,max:number}>} proteinSpans  from loadAminoAcidLayer
   */
  proteinTargets(proteinSpans) {
    const startNt = (protein) => {
      const rec = this.reference.proteins[protein];
      const codon = rec && rec.codon_nt && rec.codon_nt[0];
      return codon && codon.length ? Math.min(...codon) : -1;
    };
    const targets = [];
    const proteins = [...proteinSpans.keys()].sort(
      (a, b) => startNt(a) - startNt(b) || a.localeCompare(b));
    for (const protein of proteins) {
      const span = proteinSpans.get(protein);
      targets.push({ key: protein, label: `${protein}  (1\u2013${span.max})`, protein, offset: 0 });
      const subs = (this.subProteins[protein] || []).slice()
        .sort((a, b) => a.firstGeneAa - b.firstGeneAa);
      for (const sp of subs) {
        targets.push({ key: sp.abbrev, protein,
                       label: `\u2003${sp.abbrev} \u2014 ${sp.name} (1\u2013${sp.lastGeneAa - sp.offset})`,
                       offset: sp.offset, maxOwn: sp.lastGeneAa - sp.offset });
      }
    }
    return targets;
  }

  /**
   * A subunit position in its own numbering -> the gene's amino-acid position in NL4-3
   * numbering, which is what the atlas is keyed on.
   *
   * For gag and pol the subunit's own numbering is already the clinical convention and maps
   * to HXB2 by identity, so adding the offset is the whole job. Env and vpu are different:
   * the number a reader quotes is HXB2, and NL4-3 differs from it by indels (HXB2 Env 460 is
   * NL4-3 Env 458), so the figure has to go through the alignment. Returns null when the
   * HXB2 position has no NL4-3 equivalent.
   */
  geneAaFromOwnPosition(target, ownPosition) {
    const asGeneAa = ownPosition + (target.offset || 0);
    if (this.nl43ToHxb2[target.protein]) {
      return this.nl43PositionFromHxb2(target.protein, asGeneAa);
    }
    return asGeneAa;
  }

  /**
   * The inverse of geneAaFromOwnPosition: a gene amino-acid position (NL4-3, the atlas's own
   * key) back to the number a reader would quote for it. Round-tripping the two is what keeps
   * the coordinate box showing back what was typed into it.
   */
  ownPositionFor(target, geneAa) {
    if (this.nl43ToHxb2[target.protein]) {
      const h = this.hxb2Position(target.protein, geneAa);
      if (h !== null) return h - (target.offset || 0);
    }
    return geneAa - (target.offset || 0);
  }

  /**
   * Reverse lookup: given an HXB2 position, find the NL4-3 position.
   * Returns null if the HXB2 position doesn't exist in this protein.
   */
  nl43PositionFromHxb2(protein, hxb2AaPosition) {
    const table = this.nl43ToHxb2[protein];
    if (!table) return hxb2AaPosition;                 // identity mapping for proteins without HXB2
    
    // Reverse search: iterate through the forward mapping
    for (const [nl43Pos, hxb2Pos] of table.entries()) {
      if (hxb2Pos === hxb2AaPosition) {
        return nl43Pos;
      }
    }
    return null;                                       // HXB2 position not found
  }


  /** Sub-protein and its own numbering, e.g. Pol 253 -> {abbrev: "RT", position: 103}. */
  subProteinFor(protein, aaPosition) {
    const list = this.subProteins[protein];
    if (!list) return null;
    for (const sp of list) {
      if (aaPosition >= sp.firstGeneAa && aaPosition <= sp.lastGeneAa) {
        return { name: sp.name, abbrev: sp.abbrev, position: aaPosition - sp.offset };
      }
    }
    return null;
  }

  /** The clinical-convention label a user will recognise: "RT 103", "PR 46", "Env 160". */
  hxb2Label(protein, aaPosition) {
    const sp = this.subProteinFor(protein, aaPosition);
    if (sp) return sp.abbrev + " " + sp.position;       // Gag/Pol subunits: HXB2 = NL4-3
    const h = this.hxb2Position(protein, aaPosition);
    if (h === null) return protein + " -";              // no HXB2 equivalent (NL4-3 insertion)
    return protein + " " + h;
  }

  /** The mutation written in HXB2 convention, e.g. "RT K103N", "Env N160K", "PR M46I".
   *  This is the clinical/publication form: the residue change on the HXB2-numbered
   *  position, within the subunit for Gag/Pol. `refAa` is NL4-3's residue at the site. */
  hxb2Mutation(protein, aaPosition, refAa, altAa) {
    const sp = this.subProteinFor(protein, aaPosition);
    if (sp) return sp.abbrev + " " + refAa + sp.position + altAa;
    const h = this.hxb2Position(protein, aaPosition);
    if (h === null) return protein + " " + refAa + "ins" + altAa;   // NL4-3 insertion vs HXB2
    return protein + " " + refAa + h + altAa;
  }

  /** Inverse of hxb2Label's sub-protein half: "RT", 103 -> {protein: "Pol", position: 253}. */
  resolveSubProtein(abbrevOrName, position) {
    const want = String(abbrevOrName).toUpperCase();
    for (const [protein, list] of Object.entries(this.subProteins)) {
      for (const sp of list) {
        if (sp.abbrev.toUpperCase() === want || sp.name.toUpperCase() === want) {
          return { protein, position: position + sp.offset, subprotein: sp.abbrev };
        }
      }
    }
    return null;
  }

  // -- the model -----------------------------------------------------------------

  /** The position spline term for a class, or its genome mean where unaligned. */
  splineTerm(splineClass, refPosition) {
    if (refPosition === null || refPosition === 0) {
      // NOT zero: the sum-to-zero constraint applies to the spline coefficients, not to
      // the curve, so the curve's mean over positions is not zero.
      return { value: this.fallback[splineClass], source: "fallback" };
    }
    const idx = refPosition - this.splineFirst;
    const curve = this.splineByPosition[splineClass];
    if (idx < 0 || idx >= curve.length) {
      return { value: this.fallback[splineClass], source: "fallback" };
    }
    return { value: curve[idx], source: "aligned" };
  }

  /**
   * Rate for one substitution at one site.
   *
   * @param {number|null} refPosition  aligned NL4-3 position, or null if unaligned
   * @param {string} refBase           the submitted sequence's base at the site
   * @param {string} altBase
   * @param {string} context           3' base for non-G->A, 3' 2-mer for G->A, "" if
   *                                   it could not be read from the submitted sequence
   * @returns {object|null}  null when the model cannot be evaluated at all
   */
  rateFor(refPosition, refBase, altBase, context) {
    if (!isACGT(refBase) || !isACGT(altBase) || refBase === altBase) return null;
    const mutType = mutationType(refBase, altBase);
    let logRate, splineClass, contextUsed, contextSource;

    if (mutType === "G->A") {
      // The 2-mer is required: a G->A rate without it is not defined by this model.
      if (!context || context.length !== 2 || !isACGT(context[0]) || !isACGT(context[1])) {
        return null;
      }
      const raw2mer = context;
      const level = this.gaPooled2mers.has(raw2mer) ? "OTHER" : raw2mer;
      const intercept = this.bundle.intercept_bio_ga[level];
      if (intercept === undefined) return null;
      // Spline class follows the RAW 2-mer, the intercept follows the POOLED level.
      const cls = raw2mer === "AA" ? "AA" : raw2mer === "GA" ? "GA" : "SHARED";
      splineClass = "GtoA_" + cls;
      const sp = this.splineTerm(splineClass, refPosition);
      logRate = intercept + sp.value;
      contextUsed = level;
      contextSource = "sequence";
      var positionSource = sp.source;
    } else {
      const beta = this.bundle.beta_sub_type_bio[mutType];
      if (beta === undefined) return null;
      splineClass = "shared_sub";
      const sp = this.splineTerm(splineClass, refPosition);
      let delta = 0;
      if (context && context.length >= 1 && isACGT(context[0])) {
        const key = mutType + "|" + context[0];
        // Pooled cells were fitted type-only; their delta must not be applied.
        delta = this.pooledCells.has(key) ? 0 : (this.bundle.beta_sub_context_delta[key] || 0);
        contextUsed = context[0];
        contextSource = this.pooledCells.has(key) ? "pooled" : "sequence";
      } else {
        // 3' base unreadable (end of the submitted sequence, or an ambiguous base).
        // Fall back to the type-only rate and say so.
        contextUsed = "";
        contextSource = "unavailable";
      }
      logRate = beta + delta + sp.value;
      var positionSource = sp.source;
    }

    const rate = Math.exp(logRate);
    const iv = this._interval(mutType, refPosition, contextUsed, rate);
    return {
      mutType, rate, lo95: iv.lo95, hi95: iv.hi95,
      splineClass, context: contextUsed, contextSource,
      positionSource, intervalSource: iv.source,
    };
  }

  _interval(mutType, refPosition, context, rate) {
    if (!this.intervals) return { lo95: null, hi95: null, source: "unavailable" };
    // 1. the published interval for this exact site and type
    if (refPosition) {
      const hit = this.intervalByPosition.get(mutType + ":" + refPosition);
      if (hit) {
        return { lo95: rate * Math.exp(hit[0]), hi95: rate * Math.exp(hit[1]), source: "published" };
      }
    }
    // 2. the mean log-ratio for this (type, context)
    const tc = this.intervals.by_type_context[mutType + "|" + context];
    if (tc) {
      return { lo95: rate * Math.exp(tc.lo_ratio), hi95: rate * Math.exp(tc.hi_ratio),
               source: "type-context" };
    }
    // 3. the mean log-ratio for this type
    const t = this.intervals.by_type[mutType];
    if (t) {
      return { lo95: rate * Math.exp(t.lo_ratio), hi95: rate * Math.exp(t.hi_ratio), source: "type" };
    }
    return { lo95: null, hi95: null, source: "unavailable" };
  }

  /**
   * Insertion or deletion rate at one aligned position.
   *
   * log rate = intercept_bio + beta_mut_bio[kind] + that kind's position spline.
   * There is no nucleotide context, so the rate is the NL4-3 rate at the coordinate.
   * It matches the published median to the rounding of the shipped coefficients.
   */
  indelRate(refPosition, kind) {
    if (kind !== "insertion" && kind !== "deletion") return null;
    const beta = this.bundle.beta_mut_bio && this.bundle.beta_mut_bio[kind];
    if (beta === undefined || typeof this.bundle.intercept_bio !== "number") return null;
    const sp = this.splineTerm(kind, refPosition);
    return {
      rate: Math.exp(this.bundle.intercept_bio + beta + sp.value),
      positionSource: sp.source,
    };
  }

  // -- reading context out of a submitted sequence ---------------------------------
  // Read from the SUBMITTED sequence, not from NL4-3. This is the step where a user's
  // sequence genuinely changes the answer.
  _contextAt(sequence, queryPosition, length) {
    const start = queryPosition;                  // 1-based site -> 0-based next base
    if (start + length > sequence.length) return "";
    const ctx = sequence.slice(start, start + length);
    for (const ch of ctx) if (!isACGT(ch)) return "";
    return ctx;
  }

  // --------------------------------------------------------------------------------
  // THE CUSTOM ATLAS
  // --------------------------------------------------------------------------------
  /**
   * @param {object} job
   * @param {string} job.name       sequence label from the FASTA header
   * @param {string} job.sequence   the submitted nucleotide sequence, upper case
   * @param {object} job.alignment  output of alignToReference()
   */
  customAtlas({ name, sequence, alignment }) {
    const seq = sequence.toUpperCase();
    const refToQuery = alignment.refToQuery;
    const summary = {
      name,
      queryLength: seq.length,
      identity: alignment.identity,
      alignedBases: alignment.alignedCount,
      unalignedSpans: alignment.unalignedSpans.length,
      droppedOutsideRegion: 0,
      fallbackPositions: 0,
      gaContextUnavailable: 0,
      contextUnavailable: 0,
      codonsAnalysed: 0,
      codonsSkippedUnaligned: 0,
      codonsSkippedOutsideRegion: 0,
      codonsSkippedAmbiguous: 0,
      proteins: {},
    };

    // Positions of the submitted sequence that aligned OUTSIDE the modelled region.
    // The brief's scope exception: align against all 8,618 bp so the frame is anchored,
    // then drop everything outside 1-8,213 and say how much was dropped.
    for (let p = 1; p <= this.reference.length; p++) {
      if (p >= this.firstPosition && p <= this.lastPosition) continue;
      if (refToQuery[p]) summary.droppedOutsideRegion++;
    }

    // ---- nucleotide layer --------------------------------------------------------
    const nt = [];
    const fallbackSeen = new Set();
    for (let refPos = this.firstPosition; refPos <= this.lastPosition; refPos++) {
      const q = refToQuery[refPos];
      if (!q) continue;                                   // reference position not covered
      const base = seq[q - 1];
      if (!isACGT(base)) continue;                        // ambiguous base in the submission
      const ctx1 = this._contextAt(seq, q, 1);
      const ctx2 = this._contextAt(seq, q, 2);
      for (const alt of BASES) {
        if (alt === base) continue;
        const ctx = base === "G" && alt === "A" ? ctx2 : ctx1;
        const r = this.rateFor(refPos, base, alt, ctx);
        if (!r) {
          if (base === "G" && alt === "A") summary.gaContextUnavailable++;
          continue;
        }
        if (r.contextSource === "unavailable") summary.contextUnavailable++;
        if (r.positionSource === "fallback") fallbackSeen.add(refPos);
        nt.push({
          ref_position: refPos,
          query_position: q,
          ref_nucleotide: base,
          mut_type: r.mutType,
          context_3p: r.context,
          rate: r.rate,
          lo95: r.lo95,
          hi95: r.hi95,
          spline_class: r.splineClass,
          position_source: r.positionSource,
          interval_source: r.intervalSource,
          context_source: r.contextSource,
        });
      }
    }
    summary.fallbackPositions = fallbackSeen.size;

    // A per-site rate cache so a codon path never recomputes a site it shares with
    // another path. Key is `${refPos}:${mutType}`.
    const siteCache = new Map();
    for (const row of nt) siteCache.set(row.ref_position + ":" + row.mut_type, row);

    // ---- codon paths and amino-acid mutability -----------------------------------
    const paths = [];
    const aaGroups = new Map();      // protein|position|alt_aa -> accumulator

    for (const [protein, rec] of Object.entries(this.reference.proteins)) {
      const pstat = summary.proteins[protein] = {
        codonsInReference: rec.n_codons, codonsAnalysed: 0,
        firstAnalysed: null, lastAnalysed: null,
      };
      let codonOrdinal = 0;
      for (let i = 0; i < rec.codon_nt.length; i++) {
        const aaPosition = rec.first_aa + i;
        const refNt = rec.codon_nt[i];                    // the three NL4-3 positions

        // Every one of the three positions must be inside the modelled region and
        // covered by the submitted sequence, or the codon is not analysable.
        if (refNt.some((p) => p < this.firstPosition || p > this.lastPosition)) {
          summary.codonsSkippedOutsideRegion++; continue;
        }
        const qNt = refNt.map((p) => refToQuery[p]);
        if (qNt.some((q) => !q)) { summary.codonsSkippedUnaligned++; continue; }

        const codon = qNt.map((q) => seq[q - 1]).join("");
        if (![...codon].every(isACGT) || !(codon in GENETIC_CODE)) {
          summary.codonsSkippedAmbiguous++; continue;
        }
        const refAa = GENETIC_CODE[codon];
        codonOrdinal++;
        pstat.codonsAnalysed++;
        if (pstat.firstAnalysed === null) pstat.firstAnalysed = aaPosition;
        pstat.lastAnalysed = aaPosition;
        summary.codonsAnalysed++;

        const hxb2Label = this.hxb2Label(protein, aaPosition);

        // All 63 alternative codons: 1-, 2- and 3-step alike.
        for (const c1 of BASES) for (const c2 of BASES) for (const c3 of BASES) {
          const mutantCodon = c1 + c2 + c3;
          if (mutantCodon === codon) continue;
          const altAa = GENETIC_CODE[mutantCodon];

          // Path rate = PRODUCT of the per-site rates at the changed slots, in codon
          // slot order. Bounds are the same product over the marginal bounds.
          let muPath = 1, muLo = 1, muHi = 1, nChanges = 0, ok = true;
          const typeParts = [];
          for (let slot = 0; slot < 3; slot++) {
            if (codon[slot] === mutantCodon[slot]) continue;
            nChanges++;
            const mt = mutationType(codon[slot], mutantCodon[slot]);
            const site = siteCache.get(refNt[slot] + ":" + mt);
            if (!site || site.lo95 === null) { ok = false; break; }
            typeParts.push(mt);
            muPath *= site.rate; muLo *= site.lo95; muHi *= site.hi95;
          }
          if (!ok) continue;        // a site in this path could not be evaluated

          paths.push({
            protein, position: aaPosition, codon, mutant_codon: mutantCodon,
            ref_aa: refAa, alt_aa: altAa, n_nt_changes: nChanges,
            mutation_type_path: typeParts.join(";"),
            mu_path: muPath, mu_path_lo95: muLo, mu_path_hi95: muHi,
          });

          const key = protein + "|" + aaPosition + "|" + altAa;
          let acc = aaGroups.get(key);
          if (!acc) {
            acc = {
              protein, position: aaPosition, ref_aa: refAa, alt_aa: altAa,
              codon, query_codon_number: codonOrdinal, query_nt_start: qNt[0],
              hxb2_label: hxb2Label,
              mu_aa: 0, mu_aa_lo95: 0, mu_aa_hi95: 0,
              min_nt_changes: Infinity, n_paths: 0,
            };
            aaGroups.set(key, acc);
          }
          // Summed over ALL paths, not only the minimal ones -- verified against the
          // published table, where minimal-only reproduces n_paths for 23,718 of 62,518
          // rows and all-paths reproduces all of them.
          acc.mu_aa += muPath;
          acc.mu_aa_lo95 += muLo;
          acc.mu_aa_hi95 += muHi;
          acc.n_paths += 1;
          if (nChanges < acc.min_nt_changes) acc.min_nt_changes = nChanges;
        }
      }
    }

    const aa = [];
    for (const acc of aaGroups.values()) {
      aa.push({
        ...acc,
        synonymous: acc.alt_aa === acc.ref_aa,
        aa_change: acc.ref_aa + acc.position + acc.alt_aa,        // NL4-3 gene numbering (traceability)
        hxb2_change: this.hxb2Mutation(acc.protein, acc.position, acc.ref_aa, acc.alt_aa),
        accessible_class: acc.min_nt_changes,
        components: "codon+context" + (acc.position ? "+position" : ""),
      });
    }
    // Stable order: protein, then position, then target amino acid.
    aa.sort((x, y) => x.protein.localeCompare(y.protein) || x.position - y.position ||
                      x.alt_aa.localeCompare(y.alt_aa));
    paths.sort((x, y) => x.protein.localeCompare(y.protein) || x.position - y.position ||
                         x.mutant_codon.localeCompare(y.mutant_codon));

    summary.ntRows = nt.length;
    summary.aaRows = aa.length;
    summary.pathRows = paths.length;
    return { nt, aa, paths, summary };
  }
}
