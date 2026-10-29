<!-- SPDX-License-Identifier: MIT -->
# Broadly neutralizing antibody escape catalog — provenance and status

## Read this first

**These annotations are provisional. They are my reading of the primary bNAb literature,
not a database export, and the `citation` column is deliberately empty.** Every row in
`data/v0.1/bnab_catalog.csv` carries `provisional = TRUE`, and the portal marks every
bNAb result as provisional in the interface. Do not present any bNAb number as
authoritative, and do not clear the provisional flag, until the steps below are done.

Unlike the DRM catalog — which is rebuilt from the Stanford SDRD exports already in this
repository and reproduces the published Figure 6 set exactly — there is **no bNAb source
data in the repository**. This file was written to unblock the engineering: the escape
positions, the glycan dependencies and the sequon-loss computation can all be built and
tested against it now, and the annotations swapped for curated ones later without touching
any code.

## What the catalog contains

23 escape annotations across six epitope classes and 13 antibodies:

| Epitope class | Antibodies | Principal mechanism |
|---|---|---|
| CD4-binding site | VRC01, 3BNC117, N6 | loop D / CD4bs contacts, N197 glycan |
| V3-glycan (N332 supersite) | 10-1074, PGT121, PGT128 | N332 and N301 glycans, V3 base |
| V2-apex | PG9, PGT145, CAP256-VRC26 | N160 glycan, C-strand (166–171) |
| MPER | 10E8, 4E10 | W672, F673, W680 |
| gp120–gp41 interface | PGT151 | N611 glycan |
| fusion peptide | VRC34.01 | FP N-terminus (512), N88 glycan |

Columns: `antibody`, `epitope_class`, `antibody_class`, `mechanism` (contact /
glycan_sequon / conformational), `hxb2_position`, `wt_aa`, `escape_aa` (empty for a
glycan sequon, meaning *any* change that destroys the N-X-S/T motif), the resolved
`nl43_env_position` and `nl43_wt_aa`, `nl43_matches_expected`, `evidence_type`,
`citation` (empty), `provisional`, and a `note`.

## NL4-3 is not a neutral reference for Env

NL4-3 is a lab-adapted subtype-B strain with an atypical V1/V2/V3. At four annotated
positions its residue is **not** the canonical epitope wild type, so a rate computed from
NL4-3 at those positions is the rate *from NL4-3's residue*, which is not the escape event
the literature describes. These rows are flagged `nl43_matches_expected = FALSE`:

| Antibody | HXB2 position | Expected WT | NL4-3 has |
|---|---|---|---|
| VRC01, N6 | 279 (loop D) | N | D |
| 10-1074 | 325 (V3 base) | D | N |
| PG9 | 169 (V2 C-strand) | K | V |

The glycan-dependent annotations are on firmer ground: NL4-3 **does** carry the key
sequons — N160 (V2 apex), N276 (CD4bs), N301 and N332 (V3 supersite), N611 (interface),
N88 (fusion peptide) — all confirmed present by direct translation of the reference and
listed in `glycan_sequons.json`.

## Before the provisional flag can be cleared

1. Replace each annotation with one traceable to a specific paper or to CATNAP/LANL, and
   fill the `citation` column. The escape residues here are the well-known ones, but a
   resource must cite them.
2. Decide the scope deliberately: which antibodies, and for each whether to list specific
   escape substitutions, glycan dependencies, or both. The current set is illustrative,
   not a considered inclusion list.
3. Resolve the four `nl43_matches_expected = FALSE` rows — either note that the escape is
   defined relative to consensus rather than NL4-3, or drop them.
4. Confirm the glycan-sequon mechanism is the right abstraction for each glycan-dependent
   antibody (some also require specific underlying residues, not only the glycan).

## How the rates are computed

A `contact` annotation is treated exactly like a DRM: the amino-acid mutability of the
specific change, from the atlas. A `glycan_sequon` annotation additionally gets a
**sequon-loss rate** — the summed mutability of every change at the N, X and S/T codons
that destroys the N-X-S/T motif — because any sequon-destroying change removes the glycan,
so escape supply is spread over three codons. Both numbers are shown and clearly
distinguished. The definition of which substitutions destroy a sequon is in
`js/codon_sensitivity.js` and documented in `docs/DRM_METHODS.md`; N-X-C is treated as
non-glycosylating, i.e. as loss.
