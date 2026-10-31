<!-- SPDX-License-Identifier: MIT -->
# bNAb escape catalog — provenance

## Status

Curated from primary literature. The earlier draft of this file described a provisional,
uncited set written to unblock engineering; that set has been **replaced in full**. Every row
in `data/v0.1/bnab_catalog.csv` now cites a study that measured an effect on neutralisation,
and no row carries a provisional flag.

## Inclusion standard

A substitution is listed only where a cited paper reports a measured effect on neutralisation:
escape selected in vivo, a rebound virus in a treated participant, or an explicit resistance
assignment. Structural contact residues and glycan positions with no measured effect are
deliberately excluded, as are substitutions a source names but reports as *not* conferring
resistance on their own.

`evidence_strength` records which of three tiers a row belongs to:

| Tier | Meaning |
|---|---|
| `primary` | the substitution was shown to reduce neutralisation on its own |
| `secondary` | no effect alone; confers escape only alongside a partner mutation |
| `clinical` | emerged in a treated participant and was attributed to that antibody there |

## Coverage

24 rows, 6 antibodies, 4 epitope classes:

| Antibody | Class | Rows | With a rate |
|---|---|---|---|
| PGDM1400 | V2-apex | 8 | 5 |
| VRC07-523LS | CD4bs | 5 | 5 |
| PGT121 | V3-glycan | 3 | 2 |
| 10-1074 | V3-glycan | 4 | 4 |
| 3BNC117 | CD4bs | 2 | 2 |
| N6 | CD4bs | 2 | 2 |

Sources: Galvez et al. 2026 (10.1016/j.immuni.2026.01.031); Julg et al. 2024
(10.1038/s41591-024-03247-5); Julg et al. 2022 (10.1038/s41591-022-01815-1); Romero et al. 2026
(10.1073/pnas.2524021123); Walker-Sperling et al. 2025 (10.1038/s41467-025-66321-7).

**Not yet covered.** VRC01, 10E8, CAP256-VRC26 and PGT151/VRC34 are clinically advanced but are
not in this catalog: their per-substitution neutralisation data sits in figures and
supplementary tables that were not retrievable as text. They are absent rather than
weakly-evidenced, which is the intended failure mode. Stabell et al. 2026
(10.1038/s41564-026-02347-x) reports a large parallel-selection screen for 3BNC117 and 10-1074
across 15 isolates whose per-mutation results are likewise figure-bound; adding them would
substantially extend the V3-glycan and CD4bs rows.

## NL4-3 applicability

Escape mutations are published in HXB2 numbering against varied viral backgrounds — Galvez et
al. used REJO.c and JR-CSF, not NL4-3. The atlas can only rate a change that starts from the
NL4-3 codon, so each row was checked against the NL4-3 residue at the mapped position.

Where NL4-3 does not carry the published wild type, the row is kept with its citation and
`applies_to_nl43 = False`, and **no rate is shown**. Relabelling such a row to NL4-3's own
residue was considered and rejected: D279K is a different mutation from N279K and need not be
an escape mutation at all. Four rows are in this state, and two of them are informative in
their own right — at HXB2 190 and HXB2 325 NL4-3 already carries the escape residue, so it is
natively "escaped" at those sites.

## One mutation per row

The catalog does not combine mutations into genotype-level rates, even where a source reports
that two are needed together (VRC07-523LS D279A requires N276D or N460D for complete escape).
A combined rate would assume sequential accumulation within one lineage; recombination can
assemble the same genotype without that path, and recombination is not modelled here.

## What the rate column means

`mu_aa` is the per-replication probability of that amino-acid substitution arising from the
NL4-3 codon, with `n_star = 1/mu_aa` the number of replication events giving one expected de
novo copy. It is mutational supply only. It says nothing about whether the escape variant is
fit, is selected, or persists — Galvez et al. show path choice depends on fitness cost and
resistance benefit as well as supply, and the dominant VRC07-523LS path (D279A) has the
*lowest* one-step rate in this catalog.

## Regenerating and checking

`node test/bnab_cases.mjs --version v0.1` asserts the schema, the citation and evidence rules,
the applicability rule, and that every shipped rate equals a live engine recomputation. It also
pins the N6 A281D/A281K contrast, which is the catalog's clearest supply result.
