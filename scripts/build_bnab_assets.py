#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_bnab_assets.py -- the bNAb escape catalog and the Env glycan sequons.

Imported by build_web_assets.py so a from-scratch rebuild regenerates them.

TWO KINDS OF CONTENT, DIFFERENT SOURCES
    glycan_sequons.json is DERIVED: every N-X-S/T motif in the modelled region of NL4-3
    Env, computed from the reference sequence. It is reproducible and authoritative.

    bnab_catalog.csv is CURATED SOURCE DATA, copied through from
    data/source/bnab_catalog_curated.csv. There is no upstream database to export from:
    the rows are a curation of the primary literature, each citing a study that measured
    an effect on neutralisation, so the committed CSV *is* the source.

    This script previously authored the catalog from a table held here. It no longer does.
    Generating it would silently overwrite the curation on any rebuild -- edit the source
    CSV instead, and see docs/BNAB_PROVENANCE.md for the inclusion standard.
"""
import argparse, csv, json, pathlib, sys
import pandas as pd

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# antibody, epitope_class, antibody_class, mechanism, hxb2_position, wt_aa, escape_aa, note
# mechanism: contact | glycan_sequon | conformational. escape_aa "" for a glycan sequon.
# The ANNOTATIONS table that used to live here authored the retired provisional
# catalog. It was removed when the catalog became curated source data; see the
# module docstring.

def build_bnab(outdir, reference, hxb2, verbose=True):
    env = reference["proteins"]["Env"]
    env_aa = "".join(GENETIC_CODE[c] for c in env["codon_seq"])
    h2nl = {int(h): (int(n) if n is not None else None)
            for h, n in hxb2["hxb2_to_nl43"]["Env"].items()}

    # --- glycan sequons (derived) ---
    sequons = []
    for i in range(len(env_aa) - 2):
        if env_aa[i] == "N" and env_aa[i+1] != "P" and env_aa[i+2] in ("S", "T"):
            nlp = env["first_aa"] + i
            hpos = next((h for h, n in h2nl.items() if n == nlp), None)
            sequons.append({"nl43_env_position": nlp, "hxb2_position": hpos,
                            "motif": env_aa[i:i+3],
                            "codon_nt": [env["codon_nt"][i+k] for k in range(3)]})

    # --- bNAb catalog (curated source, copied through and validated) ---
    source = pathlib.Path(__file__).resolve().parent.parent / "data/source/bnab_catalog_curated.csv"
    if not source.exists():
        sys.exit(f"ERROR: curated bNAb catalog missing at {source}. It is source data, not "
                 f"generated; restore it from version control rather than rebuilding it.")
    with open(source, newline="") as fh:
        rows = list(csv.DictReader(fh))
    REQUIRED = {"antibody", "epitope_class", "hxb2_position", "wt_aa", "escape_aa",
                "applies_to_nl43", "evidence_strength", "citation", "citation_doi"}
    missing = REQUIRED - set(rows[0].keys() if rows else [])
    if missing:
        sys.exit(f"ERROR: curated bNAb catalog is missing columns: {sorted(missing)}")
    uncited = [r for r in rows if not r.get("citation_doi", "").startswith("10.")]
    if uncited:
        sys.exit(f"ERROR: {len(uncited)} bNAb row(s) lack a DOI; every row must cite a source.")

    # A row may only carry a rate if NL4-3 holds the published wild type. Re-derive that here
    # rather than trusting the flag, so a hand edit to the source cannot ship a wrong rate.
    flagged = []
    for r in rows:
        hpos = int(r["hxb2_position"])
        nl = h2nl.get(hpos)
        nl_aa = env_aa[nl-1] if nl else None
        applies = (nl_aa == r["wt_aa"])
        if applies != (str(r["applies_to_nl43"]) == "True"):
            sys.exit(f"ERROR: {r['antibody']} {r['wt_aa']}{hpos}{r['escape_aa']} claims "
                     f"applies_to_nl43={r['applies_to_nl43']} but NL4-3 carries {nl_aa}.")
        if not applies:
            flagged.append((r["antibody"], hpos, r["wt_aa"], nl_aa))
            if r.get("mu_aa"):
                sys.exit(f"ERROR: {r['antibody']} {r['wt_aa']}{hpos} does not apply to NL4-3 "
                         f"but carries a rate. Such rows must have no rate.")

    bnab_sequon_hxb2 = {int(r["hxb2_position"]) for r in rows}
    for s_ in sequons:
        s_["bnab_dependent"] = s_["hxb2_position"] in bnab_sequon_hxb2

    cols = list(rows[0].keys())
    with open(outdir / "bnab_catalog.csv", "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader()
        for r in rows: w.writerow(r)
    (outdir / "glycan_sequons.json").write_text(json.dumps({
        "note": ("N-X-S/T sequons (X != P) in modelled NL4-3 Env, from the reference sequence. "
                 "codon_nt gives the three NL4-3 reference nucleotide positions of the N, X and "
                 "S/T codons. bnab_dependent marks sequons the provisional bNAb catalog relies on."),
        "n_sequons": len(sequons),
        "env_modelled_range_nl43": [env["first_aa"], env["last_aa"]],
        "sequons": sequons}, indent=1))
    if verbose:
        print(f"    bNAb catalog: {len(rows)} curated annotations, "
              f"{len(set(r['antibody'] for r in rows))} antibodies, "
              f"{len(set(r['epitope_class'] for r in rows))} epitope classes")
        print(f"    glycan sequons: {len(sequons)} ({sum(s['bnab_dependent'] for s in sequons)} bNAb-dependent)")
        if flagged:
            print("    NL4-3 != expected WT at: " +
                  ", ".join(f"{ab} {h}({wt}->{got})" for ab, h, wt, got in flagged))
    return rows, sequons


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", default="v0.1")
    args = ap.parse_args()
    repo = pathlib.Path(__file__).resolve().parent.parent
    outdir = repo / "data" / args.version
    reference = json.loads((outdir / "nl43_reference.json").read_text())
    hxb2 = json.loads((outdir / "hxb2_align.json").read_text())
    build_bnab(outdir, reference, hxb2)
    print("\nNOTE: run build_web_assets.py to refresh the manifest. bNAb rows are PROVISIONAL.")


if __name__ == "__main__":
    main()
