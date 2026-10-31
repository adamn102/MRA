<!-- SPDX-License-Identifier: MIT -->
# Data dictionary — atlas v0.1

Rates throughout are the **control-corrected biological mutation rate per site per
replication**. Counts from wild-type libraries were modelled jointly with a control
library capturing sequencing and PCR error; the released quantity is the biological
component. Rates in the shipped CSVs are rounded to **6 significant figures** (~1e-6
relative, four orders of magnitude finer than the validated tolerance).

Atlas data are licensed **CC-BY-4.0**. Portal code is MIT.

## A warning about `bio_rate_median` versus `bio_rate_mean`

`bio_rate_median` is the posterior median and `bio_rate_mean` the posterior mean. They
are **not interchangeable**. The posterior mean of a rate is not the exponential of the
posterior mean of the linear predictor; the median is, because `exp` is monotone. The
model bundle reproduces the *median*, and the amino-acid layer is built from the
*median*. Substituting the mean changes path rates by up to 22%. Use
`bio_rate_median` unless you specifically want the mean, and say which you used.

## `nt_level_supply.csv` — 41,065 rows

One row per reference position (1–8,213) and mutation type: the 3 substitutions NL4-3
can make at each position, plus insertion and deletion.

| Column | Meaning |
|---|---|
| `ref_position` | NL4-3 reference position, 1-based. Position 1 is the A of the gag start codon. Internal coordinate — report HXB2 protein numbering to readers. |
| `ref_nucleotide` | Reference base at that position. |
| `mut_type` | One of the 12 substitutions (`A->C` … `T->G`), `insertion` or `deletion`. |
| `context_3p` | 3′ context used by the model: the **single 3′ base** for substitutions other than G→A; the **3′ 2-mer at the G site** for G→A, with the sparse 2-mers CC, CG, TC and TT shown as the pooled level `OTHER`; empty for indels, which have no nucleotide context. |
| `bio_rate_median` | Posterior median biological rate. **The canonical column.** |
| `bio_rate_lo95`, `bio_rate_hi95` | 95% credible interval. |
| `bio_rate_mean` | Posterior mean. See the warning above. |
| `wt_count`, `wt_coverage` | Observed mutation count and callable bases in the wild-type libraries. |
| `control_count`, `control_coverage` | The same for the control library. Empty for the three WT-only types (`C->A`, `G->T`, `deletion`). |
| `control_corrected` | Whether this type was fitted jointly with the control library (`True`) or WT-only (`False`). |

## `aa_level_supply.csv` — 62,518 rows

One row per (protein, amino-acid position, target amino acid).

| Column | Meaning |
|---|---|
| `protein` | `Gag`, `Pol`, `Vif`, `Vpr`, `Vpu`, `Tat`, `Rev`, `Env`, `Nef`. Subunits (PR, RT, RNase, IN, p17…) are positions within `Gag`/`Pol`; see `hxb2_align.json` → `protein_offsets`. |
| `position` | Amino-acid position in that protein, NL4-3 numbering. |
| `ref_aa`, `alt_aa` | Reference and target amino acid; `*` is a stop codon. |
| `aa_change` | `ref_aa` + position + `alt_aa`, e.g. `M184V`. |
| `synonymous` | `alt_aa` equals `ref_aa`. |
| `mu_aa` | **Summed** rate over every path from the reference codon to any codon encoding `alt_aa` — including 2- and 3-substitution routes, not only the minimal ones. |
| `mu_aa_lo95`, `mu_aa_hi95` | The same sum over the marginal lower and upper bounds. |
| `min_nt_changes` | Fewest nucleotide substitutions on any path to this amino acid (1, 2 or 3). |
| `n_paths` | Number of paths summed. Counts all alternative codons encoding `alt_aa`, at every step count. |

## `aa_paths/<Protein>.csv` — 187,992 rows in total

Every codon paired with **all 63** alternative codons, for each of 2,984 codons.

| Column | Meaning |
|---|---|
| `position` | Amino-acid position in that protein. |
| `codon`, `mutant_codon` | Reference codon and the alternative. |
| `ref_aa`, `alt_aa`, `aa_change` | As above. |
| `codon_position_1/2/3` | The three NL4-3 **nucleotide** positions of this codon. Not necessarily consecutive: one *rev* codon spans a 2,325 nt splice junction. |
| `n_nt_changes` | 1, 2 or 3. |
| `mutation_type_path` | Semicolon-separated substitutions in codon-slot order, e.g. `T->A;G->A`. |
| `mu_path` | **Product** of `bio_rate_median` at the changed sites. |
| `mu_path_lo95`, `mu_path_hi95` | The same product over the marginal bounds. |
| `synonymous` | As above. |

### What the intervals do and do not mean

`mu_path` bounds are the product of marginal bounds; `mu_aa` bounds are their sum. These
**approximate** joint uncertainty. They are appropriate for ranking changes and for
sensitivity analysis, and they are not exact simultaneous credible intervals. A 2- or
3-substitution `mu_path` is the product of per-generation site rates, so it is a measure
of mutational accessibility over multiple replication events, not a literal
one-generation rate.

## `mra_web_model.json`

| Key | Meaning |
|---|---|
| `version`, `release_date`, `doi` | Stamped inside the file so a download is self-identifying. The portal refuses a bundle whose `version` disagrees with its directory. |
| `model` | Name of the fitted PyMC model. |
| `modelled_region` | `first_position` 1, `last_position` 8,213. |
| `position_scaling` | Mean 4107.0 and SD 2370.917730893402 used to standardise position. |
| `spline_basis` | Cubic B-spline, `df` 30, 31 basis functions. |
| `beta_sub_type_bio` | Log-rate intercept per non-G→A substitution type (11). |
| `beta_sub_context_delta` | Context delta keyed `TYPE|3′ base` (44). |
| `sub_context_pooled_cells` | **Cells whose delta must NOT be applied.** Fitted type-only because they had fewer than 20 WT counts, so the sampled delta is prior-driven, not zero (up to 0.0055 on the log scale, 0.55% in rate). |
| `intercept_bio_ga` | G→A intercept per pooled 3′ 2-mer (13 levels, including `OTHER`). |
| `ga_pooled_2mers` | `CC`, `CG`, `TC`, `TT` → the `OTHER` level. |
| `ga_spline_class_rule` | Raw 2-mer `AA` and `GA` each have a dedicated position spline; all others share `SHARED`. The spline class follows the **raw** 2-mer, the intercept the **pooled** level. |
| `intercept_bio`, `beta_mut_bio` | Indel intercept and type effects. |
| `spline_theta` | 31 coefficients per spline, for the six splines. |
| `spline_by_position` | Each spline pre-evaluated at positions 1–8,213, so the browser does no spline arithmetic. |
| `unaligned_fallback` | **Genome mean of each spline curve**, substituted where a submitted sequence cannot be aligned. Not zero: the sum-to-zero constraint applies to the coefficients, not to the curve, and the curve means range from −0.35 to +0.79. |
| `recovery_report` | Diagnostics from recovering the splines by subtraction. |

### How to evaluate the model

```
substitution other than G->A:
    log rate = beta_sub_type_bio[type]
             + beta_sub_context_delta[type|3' base]     (OMIT if the cell is pooled)
             + spline_by_position.shared_sub[position - 1]

G->A:
    log rate = intercept_bio_ga[pooled 3' 2-mer]
             + spline_by_position.GtoA_{AA|GA|SHARED}[position - 1]

insertion / deletion:
    log rate = intercept_bio + beta_mut_bio[type]
             + spline_by_position.{insertion|deletion}[position - 1]
```

## `interval_ratios.json`

Multiplicative 95% interval ratios to the median, as natural logs: `lo95 = rate *
exp(lo_ratio)`, `hi95 = rate * exp(hi_ratio)`. `by_position` is the published per-site
interval and is exact. `by_type_context` is a mean over sites, used when a submitted
sequence carries a base NL4-3 does not have at that position so no published entry
exists; within a type these log-ratios have SD 0.04–0.10, so the bound is right to about
±10%. Rows using it are flagged `interval_source = "type-context"`. Substitutions only.

## `nl43_reference.json`

Reference sequence (8,618 bp), per-position `annotation_v4` and `bio_class` as a code
dictionary plus index arrays, and per-protein codon frames (`codon_nt`, `codon_seq`)
taken from the published path table rather than re-derived — which is why overlapping
genes and the spliced *tat*/*rev* reading frames are correct without special cases.

## `hxb2_align.json`

`hxb2_to_nl43` maps HXB2 amino-acid position → NL4-3 amino-acid position for **Env and
Vpu only**; `null` marks an HXB2 position with no NL4-3 equivalent. Every other protein
maps identity. Do not assume a constant offset for Env or Vpu — the mapping is
position-by-position because of indels. `protein_offsets` gives subunit spans within a
gene in amino-acid units (Protease = pol 52–150, RT 151–590, RNase 591–710, integrase
711–999).

Whole-genome HXB2 **nucleotide** coordinates are not derivable from these tables and are
not reported. The reference starts at the gag start codon (HXB2 790) with indels
downstream, so a constant offset would be wrong. That would require a separate
full-genome alignment.

## Custom-atlas exports

Files produced by the sequence tool use the same column names as the published tables so
the two are interchangeable downstream, plus:

| Column | Meaning |
|---|---|
| `query_position`, `query_nt_start`, `query_codon_number` | Your own 1-based coordinates in the sequence you submitted. |
| `hxb2_label` | HXB2 protein numbering, e.g. `RT 103`. |
| `position_source` | `aligned` where the fitted position spline applied; `fallback` where the genome mean was substituted because the site could not be aligned. |
| `context_source` | `sequence` where the 3′ context was read from your sequence; `pooled` where the model fitted that cell type-only; `unavailable` where the 3′ base ran off the end of your sequence or was ambiguous. |
| `interval_source` | `published`, `type-context` or `type` — see `interval_ratios.json` above. |
| `components` | Which of codon identity, local context and position the row rests on. |

Indels are **absent** from custom output by design. Their predictor has no nucleotide
context, so a "custom" indel rate is just the NL4-3 indel rate at that coordinate and
says nothing about the submitted virus.


## `drm_catalog.csv` — 1,053 rows, 244 unique DRMs

Drug-resistance mutations rebuilt from the Stanford SDRD class exports and joined to the
amino-acid layer. One row per (mutation, treatment stratum).

| Column | Meaning |
|---|---|
| `drug_class` | `PI`, `NRTI`, `NNRTI`, `INI`, `CAI`. |
| `drug_class_label`, `target` | Human-readable class and its protein target. |
| `rx` | Treatment stratum: `naive`, a class label (`PI`, `NRTI`, `NNRTI`, `INI`, `CAI`) or an individual drug/regimen (`ATV`, `LPV`, `EFV`, `DOR`, `RPV`, `DTG`, `CAB`, `ABC+XTC`, `AZT+XTC`, `TFV+XTC`). |
| `mechanism` | `substitution`, `deletion` or `insertion`. Indel rows have no `mu_aa`. |
| `protein`, `protein_position`, `ref_aa`, `alt_aa` | Mature-protein coordinates (PR, RT, IN, CA). |
| `mutation`, `aa_change` | e.g. `RT M184V`; `aa_change` is in gene (Pol/Gag) numbering. |
| `gene`, `gene_position` | Atlas coordinates, the join key to `aa_level_supply.csv`. |
| `nl43_codon`, `nl43_nt_start` | NL4-3 reference codon and its first nucleotide position. |
| `atlas_ref_aa`, `nl43_matches_consensus` | The amino acid NL4-3 carries, and whether it equals the Stanford consensus. `False` means the rate is from NL4-3's residue, not the clinical wild type (e.g. IN 151). |
| `mu_aa`, `mu_aa_lo95`, `mu_aa_hi95` | Mutation rate and interval. Empty for indel rows. |
| `min_nt_changes`, `n_paths` | As in the amino-acid layer. |
| `all_seqs`, `all_muts`, `freq` | Stanford sequence counts and `freq = (all_muts + 0.5)/all_seqs` in this stratum. |
| `min_prev`, `max_prev`, `ratio_prev` | Stanford cross-subtype prevalence summaries. |
| `<subtype>_seqs`, `<subtype>_muts` | Per-subtype counts (A, B, C, D, F, G, CRF01_AE, CRF02_AG). |
| `nt_indel_rate_min/max`, `rate_note` | For indel rows only: the per-nucleotide indel rate range across the codon, flagged as context, not a substitution rate. |

## `drm_taxonomy.json`

The class→stratum selector tree, with per-class mutation counts, position ranges and the
source SDRD file. `freq_pseudocount` and the subtype list are recorded here.

## `bnab_catalog.csv` — 24 rows

Curated bNAb escape mutations. A row is included only where a cited study reports a measured
effect on neutralisation; structural contacts and glycan positions without a measured effect
are excluded, as are substitutions a source names but reports as not conferring resistance
alone. See `docs/BNAB_PROVENANCE.md` for the inclusion standard and coverage gaps.

6 antibodies (PGDM1400, VRC07-523LS, PGT121, 10-1074, 3BNC117, N6) across 4 epitope classes,
from 5 papers. 20 of the 24 rows carry a rate.

| Column | Meaning |
|---|---|
| `antibody`, `epitope_class` | e.g. `PGDM1400`, `V2-apex`. |
| `hxb2_position`, `wt_aa`, `escape_aa` | The published escape mutation, in HXB2 Env numbering. |
| `nl43_env_position`, `nl43_wt_aa` | Resolved to NL4-3 numbering; blank where the row does not apply. |
| `applies_to_nl43` | `False` where NL4-3 does not carry the published wild-type residue. Such rows keep their annotation and citation but carry **no rate**, and `note` begins `NOT APPLICABLE:` with the reason. |
| `evidence_type` | `in_vivo_selection`, `clinical_rebound` or `literature_reported`. |
| `evidence_strength` | `primary` (reduces neutralisation on its own), `secondary` (no effect alone; escapes only with a partner mutation), `clinical` (emerged in a treated participant and attributed there). |
| `effect` | One line on what the cited study observed. |
| `citation`, `citation_doi` | Author-year and DOI. Rendered as a link on the page. |
| `mu_aa`, `mu_aa_lo95`, `mu_aa_hi95` | Per-replication probability of the substitution from the NL4-3 codon, with its 95% interval. Blank where `applies_to_nl43` is `False`. |
| `min_nt_changes` | Nucleotide changes required from the NL4-3 codon. |
| `n_star` | `1/mu_aa` — replication events for one expected de novo copy. |
| `note` | Caveats, including the `NOT APPLICABLE:` reason where relevant. |

Rows are **one mutation each**. No genotype-level rates are published, even where a source
reports two changes are needed together: combining them would assume sequential accumulation
within a single lineage, which recombination can bypass.

`test/bnab_cases.mjs` checks the schema, that every row cites a resolvable DOI, that the
applicability rule holds, and that every shipped rate equals a live engine recomputation.

## `glycan_sequons.json`

Every N-X-S/T sequon (X≠P) in modelled NL4-3 Env, with the three codon nucleotide
positions of the N, X and S/T codons and a `bnab_dependent` flag. Derived from the
reference sequence.

## `lanl_codon_frequencies.csv` — optional, not shipped in v0.1

When present: `gene`, `protein`, `protein_position`, `codon`, `lanl_fraction`,
`n_sequences`, `alignment_version`. Adds circulating-codon-usage weighting to the
codon-sensitivity view. Absent in v0.1; the weighted columns show an em dash. See
`scripts/build_lanl_codon_frequencies.py`.
