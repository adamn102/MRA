<!-- SPDX-License-Identifier: MIT -->
# Caveats

The limits of what this portal can tell you. Read this before quoting a number from a
custom atlas in a manuscript.

## 1. A custom atlas is a model applied to a sequence, not a measurement of that virus

This is the single most important caveat and it is stated in the interface, in every
export header and here.

ERR-Seq measured mutation rates on NL4-3. The portal applies the model fitted on those
measurements to a sequence you submit. It does not re-measure anything. Three components
contribute, and every output row records which ones it rests on:

| Component | Determined by | Transfers to another virus? |
|---|---|---|
| **Codon identity** | entirely by your sequence | Yes. This is a real, large effect and the main reason to compute a custom atlas. |
| **Local 3′ context** | entirely by your sequence | Yes, to the extent the fitted context effects are general. |
| **Position** | the NL4-3 fit | Only where your sequence aligns to NL4-3, and only on the assumption that positional effects are a property of the genome region rather than of NL4-3 specifically. |

The demonstration figure `figures/compare_demo.png` quantifies the first of these on a real
sequence: for JR-FL *pol* (GenBank U63632, 96.78% identical to NL4-3 over the region it
covers), 51.8% of amino-acid changes at positions whose codon differs shift by more than
two-fold, against 0.06% where the codon is identical. **A fold change in a custom atlas is
a statement about codon and context, not about a different underlying mutation rate.**

## 2. NL4-3-centred estimates may not transfer to other backgrounds

Positional effects are the weakest part of the transfer. The position splines absorb
everything about a genome region that the type and context terms do not explain —
secondary structure, replication-complex behaviour, local sequence features, and any
NL4-3-specific idiosyncrasy. There is no way, from this dataset alone, to say how much of
each. A position-resolved rate for a divergent isolate inherits that uncertainty and the
portal does not quantify it.

## 3. Where a sequence cannot be aligned, the positional term is a genome average

The substituted value is the **genome mean of that spline class**, not zero: the
sum-to-zero constraint applies to the spline coefficients, not to the curve, and the curve
means range from −0.35 to +0.79 on the log scale. Every such position is flagged
`position_source = "fallback"` and counted in the run summary. A result resting mostly on
fallback positions is a result about mutation type and context only.

## 4. Intervals are approximate in two distinct ways

**For novel site-substitution pairs.** Where your sequence carries a base NL4-3 does not
have at that position, no published interval exists for that pair and the portal applies
the mean log interval-ratio for that (type, 3′ base). Within a type these ratios have SD
0.04–0.10, so the bound is right to roughly ±10%. Such rows are flagged
`interval_source = "type-context"`.

**For amino-acid mutability, always.** `mu_path` bounds are the product of the marginal
bounds at the changed sites; `mu_aa` bounds are their sum. These approximate joint
uncertainty. They are appropriate for ranking and sensitivity, and they are **not** exact
simultaneous credible intervals. This convention is inherited from the published atlas
(`aa_level_supply_ci_methodology.txt`), not introduced by the portal.

## 5. Multi-step paths are not one-generation rates

A 2- or 3-substitution `mu_path` is the product of per-generation site rates. Treat it as
a measure of **mutational accessibility over multiple replication events**, not as the
probability of a multi-nucleotide event in one generation. `mu_aa` sums 1-, 2- and 3-step
paths together, so for a change whose minimum is 1 substitution the multi-step
contribution is negligible, but for a change whose minimum is 3 the entire value is of this
kind. `min_nt_changes` is reported on every row precisely so this is visible.

## 6. G→A is dominated by APOBEC3-class editing

G→A is modelled separately for this reason: a 2-mer 3′ context, 13 intercept levels, and
three position splines where every other substitution shares one. Its biology, context
dependence and variance structure differ from the other eleven types. Any genome-wide
summary that pools G→A with the rest is dominated by it. The genome track lets you select
a single mutation type for this reason, and it is worth looking at your result both with
and without G→A.

## 7. Indels are excluded from custom output, deliberately

Their predictor is intercept + type + position spline, with no nucleotide context. For a
submitted sequence the only input is the aligned coordinate, so a "custom" indel rate is
just the NL4-3 indel rate at that coordinate and says nothing about the submitted virus.
They remain in the browse and download layer, where they are a real published result. They
play no part in amino-acid mutability, which is substitution-based.

## 8. Coverage, technical background and the modelled region

The model covers NL4-3 positions 1–8,213. *nef* is partial: 72 of roughly 206 codons.
A query beyond that returns nothing because there is no coverage, and the interface says
so rather than showing an empty table. Rates near the edges of the covered region, and at
positions with low wild-type coverage, carry wider intervals — `wt_count` and
`wt_coverage` are shipped in `nt_level_supply.csv` so you can see which those are.

## 9. What the round-trip test does and does not establish

It establishes that the browser engine reproduces the published atlas exactly: same
classifications, rates within a tolerance derived from the published posterior's own Monte
Carlo residual. It establishes nothing about whether the model is right, and nothing about
accuracy on a non-NL4-3 virus. See [validation_report.md](validation_report.md).

## 10. Rounding

Rates in the shipped CSVs are rounded to 6 significant figures. That is ~1e-6 relative,
four orders of magnitude below the validated tolerance, and irrelevant for any use of these
numbers — but it means a byte-exact comparison against the production tables will differ.


## 11. Drug-resistance and antibody-escape page

- **Supply is not frequency.** The page reports mutational supply (`mu_aa`), the rate a
  variant is generated. Clinical frequency, shown alongside, is supply filtered through
  selection, treatment history, sampling and founder effects. Do not read a high supply as
  a prediction of high prevalence, or the reverse.
- **bNAb annotations are provisional.** They are a literature draft, not a database export,
  marked provisional on every row with citations unfilled. See `docs/BNAB_PROVENANCE.md`.
  The rates are computed by the validated engine; the *annotations* are what await sign-off.
- **NL4-3 is not a neutral Env reference.** It is lab-adapted with an atypical V1/V2/V3, and
  at several bNAb and one integrase DRM position its residue differs from consensus. Those
  rows are flagged; the rate is from NL4-3's residue, not the canonical wild type.
- **Indel-coded DRMs have no rate.** In-frame codon insertions/deletions (Stanford `~`/`#`)
  are real resistance but not substitutions, so no `mu_aa` is given and the per-nucleotide
  indel rate — a frameshift, a different event — is not substituted for it.
- **Codon optimization vs accessibility class.** The headline codon-sensitivity fold range
  is among synonymous codons that reach the target at the same number of nucleotide steps.
  Codons needing more steps are a categorical accessibility difference and are reported
  separately; do not add the two.
- **Sequon-loss aggregation depends on a definition.** The summed glycan-loss rate counts
  every change destroying N-X-S/T, treating N-X-C as non-glycosylating. A different sequon
  convention would change the aggregate.
- **Multi-step DRM rates are accessibility scores.** A DRM needing 2 or 3 nucleotide changes
  has a `mu_aa` that is a product of per-generation site rates — mutational accessibility
  over several replication events, not a one-generation rate. The default view shows only
  one-nucleotide-accessible DRMs for this reason.
- **LANL weighting is observational.** When added, circulating codon usage reflects
  selection and sampling, not only mutation; a frequency-weighted rate describes the
  population, it is not a second rate estimate.
