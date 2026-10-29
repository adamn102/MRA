<!-- SPDX-License-Identifier: MIT -->
# Resistance, escape and codon-sensitivity — methods

Methods for the resistance-and-escape page, written to be lifted into a manuscript
Methods section. Everything here is computed from the shipped atlas assets; nothing is
re-fitted in the browser.

## Drug-resistance mutation catalog

The catalog of drug-resistance mutations (DRMs) was rebuilt from the Stanford HIV Drug
Resistance Database surveillance exports (`SDRD/{PI, NNRTI, RTI, INI, CAI}` files),
retaining DRM-flagged rows only. We used the per-class exports rather than a pre-merged
table because the merged table omitted the NNRTI class entirely. The DRM flags partition
without overlap between classes that share a protein: the NRTI and NNRTI flagged sets over
reverse transcriptase positions 1–560 are disjoint, so concatenating the class files
double-counts no mutation. The result is 244 unique DRMs — 59 protease, 56 NRTI, 57 NNRTI,
54 integrase and 18 capsid — reproducing, where they overlap, the protease DRM set of the
manuscript's existing figures exactly (42 of 42 treatment-naive single-nucleotide changes).

Each DRM was mapped from mature-protein numbering to the atlas's gene coordinates using the
protein offsets from the atlas alignment table (protease = pol +51, RT +150, integrase
+710, capsid = gag +132) and joined to the published amino-acid mutability layer
(`aa_level_supply.csv`). Every substitution DRM joined; no DRM was dropped.

**Mutation rate.** The reported rate is `mu_aa`, the amino-acid mutability from the atlas:
the per-replication rate summed over every nucleotide path from the reference codon to any
codon encoding the resistant amino acid. It is mutational *supply*, not clinical frequency.

**Clinical frequency.** For context we also report the clinical sequence frequency in each
treatment stratum, `(mutants + 0.5) / sequences`, identical to the convention in the
manuscript figures, in both the selected treatment stratum and treatment-naive. Supply and
frequency are shown side by side so selection — the gap between them — is visible.

**Two edge cases are handled explicitly rather than hidden.**

1. *Indel-coded DRMs.* Stanford codes in-frame codon deletions (`~`) and the T69 insertion
   complex (`#`); there are five such RT DRMs. They are kept in the catalog because an NRTI
   list without them would be wrong, but they are given **no** `mu_aa`: they are not
   amino-acid substitutions. The model does estimate a per-nucleotide insertion and
   deletion rate, but that is a single-nucleotide frameshift — a different event from an
   in-frame codon indel — so it is not quoted as the DRM's rate.

2. *NL4-3 ≠ subtype-B consensus.* At integrase 151, NL4-3 carries isoleucine, not the
   consensus valine. The V151A/I/L DRMs computed from NL4-3 are therefore changes from
   isoleucine, not from the clinician's valine. These rows are flagged; the reported rate
   is the rate from NL4-3's residue.

## Codon-sensitivity analysis

For each DRM we quantify how the rate depends on codon usage, in two dimensions.

**At the site.** The wild-type codon at the resistance position is varied over every
synonymous codon for the wild-type amino acid, holding the neighbouring amino acid fixed,
and `mu_aa` is recomputed for each. Because `mu_aa` is a sum over nucleotide paths, two
synonymous codons can differ both in the *number* of nucleotide changes needed to reach the
resistant codon (a categorical difference in accessibility) and in the *rate* at a fixed
number of changes (codon optimization proper). We separate these: the headline fold range
is computed only among synonymous codons that reach the target at the minimal step count,
and codons requiring more steps are reported separately. Conflating them produces
meaningless 10¹⁰-fold ranges driven entirely by step-count differences.

**Around the position.** The 3′ flanking context is varied over only those values reachable
by synonymous recoding of the *next* codon, so the neighbouring amino acid is held fixed and
every context corresponds to a protein-identical virus that could exist. For a substitution
whose changed nucleotide sits at the codon's 3′ edge, the site and flank dimensions couple —
the changed base reads the next codon as its context — and those combinations are flagged
rather than counted as independent.

We also report where NL4-3's own codon sits in the resulting distribution (its rank among
synonymous codons at the same step count), because NL4-3 is the reference the atlas was
measured on, not a neutral choice. Across the 239 substitution DRMs, 52 vary at least
two-fold at equal accessibility; the extreme is integrase G140R, which spans 68-fold across
the four glycine codons, with NL4-3 using the lowest-supply codon.

## Glycan sequon-loss aggregation

Broadly neutralizing antibody escape often works by destroying an N-linked glycosylation
sequon, N-X-S/T (X ≠ P), which can be lost through a change at any of three codons. For
glycan-dependent epitopes we therefore report, in addition to any specific annotated change,
the **summed** mutability of every change that destroys the sequon:

- the asparagine codon → any amino acid other than N;
- the X codon → proline;
- the S/T codon → any amino acid other than serine or threonine.

N-X-C is treated as non-glycosylating, so an S/T → C change counts as loss; a sequon scanner
that accepts N-X-C would differ here. The aggregate is the quantity that predicts escape,
since any sequon-destroying change removes the glycan, and it runs 1.6–4.2× above the single
largest contributing change across the bNAb-dependent glycans of NL4-3 Env. For the V2-apex
(N160) and V3-supersite (N332) glycans, most of the supply is the S/T position mutating to
asparagine — a single high-rate transition — rather than the commonly cited change at the
asparagine itself.

## LANL codon-frequency weighting (optional)

The codon-sensitivity analysis reports the *achievable* range of rates across synonymous
codons. How much of that range circulating virus actually occupies depends on observed codon
usage, which the atlas does not contain. When a LANL-derived codon-frequency table
(`lanl_codon_frequencies.csv`) is present, the view adds a circulating-sequence-weighted
mean rate and the fraction of sequences on an at-or-above-median-supply codon; when it is
absent, those fields show an em dash rather than a silent zero. The weighting multiplies the
per-codon `mu_aa` by the observed codon fraction; it re-fits nothing. See
`scripts/build_lanl_codon_frequencies.py` for the input schema. Observed codon usage is
itself shaped by selection and sampling, so a frequency-weighted rate is a descriptive
summary of circulating virus, not a second mutation-rate estimate.

## Coordinates and provenance

Resistance positions are reported in standard protein numbering (PR, RT, IN, CA); bNAb
positions in HXB2 Env numbering, resolved to NL4-3 through the position-by-position Env
lookup. Every exported file carries the atlas version, the model name and the retrieval
date, and bNAb exports are marked provisional.
