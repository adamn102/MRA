#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_lanl_codon_frequencies.py -- STUB. Codon usage per protein position from a LANL
alignment, so the codon-sensitivity view can weight by how common each context actually is.

STATUS: not implemented. This file documents exactly what input it needs and what it will
produce, so the step is ready to run once the alignment is in hand -- it is deliberately
not run on invented data.

WHY THIS IS A SEPARATE, OPTIONAL ASSET
    The codon-sensitivity view already answers "how much could codon choice change the
    rate" by enumerating every synonymous codon. What it cannot answer without external
    data is "how common is each of those codons in circulating virus" -- i.e. is NL4-3's
    (often low-supply) codon typical, or does most of the population already sit on a
    high-supply codon. That requires observed codon frequencies, which the MRA does not
    contain. LANL supplies them.

    When data/<version>/lanl_codon_frequencies.csv is present, the portal adds three
    columns to the codon-sensitivity view: the LANL-weighted mean rate, the fraction of
    circulating sequences on a high-supply codon, and a frequency-weighted fold range.
    When it is ABSENT, those columns show an em dash -- never a silent zero.

REQUIRED INPUT
    A codon-level alignment of a LANL-curated sequence set to the NL4-3 reading frames
    used by this atlas, one of:
      (a) a LANL "codon alignment" FASTA per protein (gag, pol, env, ...), already in
          frame and gap-handled, or
      (b) a per-sequence codon table with columns
          (sequence_id, gene, protein_position, codon).
    Plus the alignment's provenance: the LANL download date, the sequence set
    (e.g. "HIV-1 M group, filtered"), and the number of sequences.

OUTPUT -- data/<version>/lanl_codon_frequencies.csv
    gene, protein, protein_position, codon, lanl_fraction, n_sequences, alignment_version
  where lanl_fraction is the fraction of sequences carrying that codon at that position
  among those with a resolvable codon there, and the fractions at a position sum to 1.

CORRECTNESS REQUIREMENTS when this is implemented
    - Map LANL coordinates to this atlas's gene positions explicitly; do not assume LANL
      and NL4-3 numbering agree (they do not in env, because of indels -- the same reason
      the HXB2 env mapping is position-by-position).
    - Drop, and count, positions where the LANL codon cannot be placed in frame.
    - Record alignment_version in every row so a weighted rate is traceable to a specific
      LANL release, exactly as every other asset is versioned.
    - Weighting multiplies the per-codon mu_aa from codonSensitivity() by lanl_fraction;
      it does NOT re-fit anything. State in the UI that the weighting reflects observed
      codon usage, which is itself shaped by selection and sampling, not only mutation.

CAVEATS to carry into the UI (see docs/CAVEATS.md)
    LANL codon usage depends on the representativeness and alignment quality of the
    sequence set, is dominated by subtype B in most LANL downloads, and conflates
    mutational and selective shaping of codon usage. A frequency-weighted rate is a
    descriptive summary of circulating virus, not a second mutation-rate estimate.
"""
import sys

if __name__ == "__main__":
    sys.exit(
        "build_lanl_codon_frequencies.py is a stub. It needs a LANL codon alignment "
        "mapped to this atlas's reading frames; see the module docstring for the exact "
        "input and output schema. The portal runs correctly without it -- the "
        "LANL-weighted columns simply show an em dash until the file exists."
    )
