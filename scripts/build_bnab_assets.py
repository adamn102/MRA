#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_bnab_assets.py -- the provisional bNAb escape catalog and the Env glycan sequons.

Imported by build_web_assets.py so a from-scratch rebuild regenerates them.

TWO KINDS OF CONTENT, DIFFERENT SOURCES
    glycan_sequons.json is DERIVED: every N-X-S/T motif in the modelled region of NL4-3
    Env, computed from the reference sequence. It is reproducible and authoritative.

    bnab_catalog.csv is HAND-AUTHORED: the escape annotations in the ANNOTATIONS table
    below are a reading of the primary bNAb literature, not a database export. This script
    IS their source of truth -- edit the table here and rebuild. Every row is written with
    provisional=TRUE and an empty citation. See docs/BNAB_PROVENANCE.md before clearing the
    provisional flag.
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
ANNOTATIONS = [
    ("VRC01","CD4bs","3BNC117/VRC01-class","contact",368,"D","R","CD4bs contact Asp368; also impairs CD4 binding."),
    ("VRC01","CD4bs","3BNC117/VRC01-class","contact",279,"N","K","Loop D change reducing CD4bs bNAb binding."),
    ("VRC01","CD4bs","3BNC117/VRC01-class","glycan_sequon",197,"N","","N197 glycan modulates CD4bs bNAb access."),
    ("3BNC117","CD4bs","3BNC117/VRC01-class","contact",458,"G","S","Loop V5/CD4bs periphery; 3BNC117 escape."),
    ("3BNC117","CD4bs","3BNC117/VRC01-class","contact",471,"G","E","CD4bs periphery reducing 3BNC117 neutralisation."),
    ("N6","CD4bs","VRC01-class (broad)","contact",279,"N","D","Loop D; modulates even the broad N6."),
    ("10-1074","V3-glycan","PGT121/10-1074","glycan_sequon",332,"N","","N332 glycan required for 10-1074."),
    ("10-1074","V3-glycan","PGT121/10-1074","contact",325,"D","N","GDIR/324-327 V3 base."),
    ("PGT121","V3-glycan","PGT121/10-1074","glycan_sequon",332,"N","","N332 glycan dependence."),
    ("PGT128","V3-glycan","PGT128","glycan_sequon",301,"N","","N301 glycan contributes to PGT128."),
    ("PGT128","V3-glycan","PGT128","glycan_sequon",332,"N","","N332 glycan dependence."),
    ("PG9","V2-apex","PG9/PG16","glycan_sequon",160,"N","","N160 glycan central to V2 apex; N160K classic escape."),
    ("PG9","V2-apex","PG9/PG16","contact",169,"K","E","C-strand (166-171) cationic contact."),
    ("PGT145","V2-apex","PGT145","glycan_sequon",160,"N","","N160 glycan dependence."),
    ("CAP256-VRC26","V2-apex","CAP256-VRC26","contact",166,"R","S","C-strand contact."),
    ("CAP256-VRC26","V2-apex","CAP256-VRC26","glycan_sequon",160,"N","","N160 glycan dependence."),
    ("10E8","MPER","10E8","contact",672,"W","A","Core MPER contact Trp672."),
    ("10E8","MPER","10E8","contact",680,"W","S","MPER helix 680."),
    ("4E10","MPER","4E10","contact",672,"W","A","Trp672 MPER contact."),
    ("4E10","MPER","4E10","contact",673,"F","A","Phe673 MPER contact."),
    ("PGT151","interface","PGT151","glycan_sequon",611,"N","","N611 glycan at the gp120-gp41 interface."),
    ("VRC34","fusion-peptide","VRC34.01","contact",512,"A","E","Fusion-peptide N-terminus (512)."),
    ("VRC34","fusion-peptide","VRC34.01","glycan_sequon",88,"N","","N88 glycan flanks the fusion peptide."),
]


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

    # --- bNAb catalog (hand-authored) ---
    rows, flagged = [], []
    bnab_sequon_hxb2 = {h for *_, m, h, _, _, _ in ANNOTATIONS if m == "glycan_sequon"}
    for ab, ec, cls, mech, hpos, wt, esc, note in ANNOTATIONS:
        nl = h2nl.get(hpos)
        nl_aa = env_aa[nl-1] if nl else None
        matches = nl_aa == wt
        if not matches:
            flagged.append((ab, hpos, wt, nl_aa))
        rows.append({"antibody": ab, "epitope_class": ec, "antibody_class": cls,
                     "mechanism": mech, "hxb2_position": hpos, "wt_aa": wt, "escape_aa": esc,
                     "nl43_env_position": nl, "nl43_wt_aa": nl_aa,
                     "nl43_matches_expected": matches,
                     "evidence_type": "literature_provisional", "citation": "",
                     "provisional": True, "note": note})
    for s in sequons:
        s["bnab_dependent"] = s["hxb2_position"] in bnab_sequon_hxb2

    cols = ["antibody","epitope_class","antibody_class","mechanism","hxb2_position","wt_aa",
            "escape_aa","nl43_env_position","nl43_wt_aa","nl43_matches_expected",
            "evidence_type","citation","provisional","note"]
    with open(outdir / "bnab_catalog.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow(r)
    (outdir / "glycan_sequons.json").write_text(json.dumps({
        "note": ("N-X-S/T sequons (X != P) in modelled NL4-3 Env, from the reference sequence. "
                 "codon_nt gives the three NL4-3 reference nucleotide positions of the N, X and "
                 "S/T codons. bnab_dependent marks sequons the provisional bNAb catalog relies on."),
        "n_sequons": len(sequons),
        "env_modelled_range_nl43": [env["first_aa"], env["last_aa"]],
        "sequons": sequons}, indent=1))
    if verbose:
        print(f"    bNAb catalog: {len(rows)} provisional annotations, "
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
