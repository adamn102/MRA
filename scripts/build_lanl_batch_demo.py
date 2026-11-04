#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_lanl_batch_demo.py -- the worked example for the batch-analysis tab, built from the
LANL HIV-1 compendium genome alignment.

WHY A SUBSET, AND WHICH ONE
    The compendium is a phylogenetic reference set: it deliberately spans HIV-1 groups M, N,
    O and P and includes SIVcpz / SIVgor outgroups, because its purpose is to show the shape
    of the tree. That is the wrong default for this tab. Measured against the NL4-3 reference
    this atlas is built on, the groups separate cleanly:

        group M subtypes and CRFs   85-93% identity   (n = 182)
        groups N, O, P              74-78% identity   (n = 13)
        SIVcpz, SIVgor              74-77% identity   (n =  4)

    A batch mixes those two bands together and reports the result as divergence in mutational
    supply. For group M that reads correctly: the codons differ because circulating virus is
    polymorphic. For groups N/O/P and the SIV outgroups it does not -- the codons differ
    because the lineages separated before HIV-1 group M existed, and two of them are not
    HIV-1 at all. The page would be showing between-lineage distance under a label that says
    codon usage.

    So the shipped demo is group M only. The full compendium stays in data/source/ and the
    --all flag below will emit it, for anyone who wants to see what cross-lineage divergence
    looks like -- but it is not what loads when a visitor clicks the button.

INPUT   data/source/HIV1_COM_2025_genome_DNA.fasta.gz
        The compendium as downloaded: aligned (11,934 columns), gapped, one record per
        reference strain, named subtype.country.year.name.accession.

OUTPUT  data/<version>/lanl_batch_demo.fasta
        Gaps stripped, because the portal aligns every submitted sequence to NL4-3
        independently and an alignment column index means nothing to it. Record names are
        kept verbatim, so every sequence remains traceable to its GenBank accession.

PROVENANCE
    Los Alamos National Laboratory HIV Sequence Database, 2025 compendium, genome DNA
    alignment. https://www.hiv.lanl.gov/ -- the underlying sequences are GenBank records and
    the accession is the last field of each name.
"""
import argparse, gzip, pathlib, sys

# Groups that are not HIV-1 group M. Matched against the first dot-separated field of a
# record name, which is the subtype/CRF label in LANL's naming scheme.
NON_GROUP_M = {"N", "O", "P", "CPZ", "GOR"}


def read_fasta(path):
    """(name, sequence) pairs from a plain or gzipped FASTA."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, buf, out = None, [], []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def group_of(name):
    return name.split(".")[0].upper()


def build_lanl_demo(source, outdir, keep_all=False, wrap=60):
    recs = read_fasta(source)
    if not recs:
        sys.exit(f"no records read from {source}")
    lengths = {len(s) for _, s in recs}
    if len(lengths) != 1:
        print(f"  note: input is not a single alignment ({len(lengths)} distinct lengths)")

    kept, dropped = [], []
    for name, seq in recs:
        (dropped if (not keep_all and group_of(name) in NON_GROUP_M) else kept).append(
            (name, seq.replace("-", "").replace(".", "").upper()))

    if not kept:
        sys.exit("every record was filtered out; refusing to write an empty demo")
    short = [n for n, s in kept if len(s) < 5000]
    if short:
        print(f"  warning: {len(short)} kept record(s) shorter than 5 kb, e.g. {short[0]}")

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "lanl_batch_demo.fasta"
    with open(out, "w") as fh:
        for name, seq in kept:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")
    print(f"  lanl_batch_demo.fasta: {len(kept)} sequences, "
          f"{out.stat().st_size / 1e6:.2f} MB"
          + (f" ({len(dropped)} non-group-M records excluded)" if dropped else " (all records)"))
    return out, len(kept), len(dropped)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source",
                    default="data/source/HIV1_COM_2025_genome_DNA.fasta.gz")
    ap.add_argument("--version", default="v0.1")
    ap.add_argument("--all", action="store_true",
                    help="keep groups N/O/P and the SIV outgroups too (not the shipped default)")
    a = ap.parse_args()
    root = pathlib.Path(__file__).resolve().parent.parent
    build_lanl_demo(root / a.source, root / "data" / a.version, keep_all=a.all)
