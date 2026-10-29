#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_drm_assets.py -- the drug-resistance mutation catalog shipped with the portal.

Imported by build_web_assets.py so that data/<version>/ has a single owner and the
manifest covers everything. Runnable standalone for development.

WHY THIS IS REBUILT FROM THE SDRD CLASS FILES RATHER THAN FROM variant_prevalence.xlsx
    The merged workbook contains no NNRTI mutations at all. Its RT DRM set is entirely
    NRTI/thymidine-analogue -- M41L, K65R, M184V, T215Y, K219Q and so on -- and K103N,
    Y181C, E138K, G190A and V106M are absent. The raw SDRD exports do have them. A
    resource page listing NRTI but not NNRTI resistance would read as an error to any
    clinician, so the catalog is rebuilt from the per-class exports.

    The DRM flags partition cleanly by class: RTI and NNRTI both cover RT positions
    1-560, and the overlap between their DRM-flagged sets is exactly zero. So
    concatenating the five class files double-counts nothing.

    NNRTI vs NNRTI_II: their DRM content is identical (285 rows, 57 unique changes, same
    Rx strata). NNRTI_II differs only by adding background non-DRM positions. This script
    uses NNRTI; switching makes no difference to the catalog.

COORDINATES
    Stanford numbers each mature protein from 1. The portal's atlas is indexed by gene
    (Pol, Gag) amino-acid position, so each class carries the offset from the workbook
    alignment sheet: protease +51, RT +150, integrase +710, capsid(p24) +132. Every DRM
    position in the exports falls inside its subunit, so no DRM lands in the RNase H
    range where Stanford's RT numbering (1-560) runs past RT proper (1-440).
"""
import argparse, glob, json, pathlib, sys
import pandas as pd
import numpy as np

# Stanford class -> (mature protein label, atlas gene, subunit name on the alignment sheet)
CLASS_SPEC = {
    "PI":    {"pattern": "PI.*.csv",    "protein": "PR",  "gene": "Pol", "subunit": "Protease",
              "label": "Protease inhibitors", "target": "protease"},
    "NRTI":  {"pattern": "RTI.*.csv",   "protein": "RT",  "gene": "Pol", "subunit": "RT",
              "label": "Nucleoside RT inhibitors", "target": "reverse transcriptase"},
    "NNRTI": {"pattern": "NNRTI.M*.csv", "protein": "RT", "gene": "Pol", "subunit": "RT",
              "label": "Non-nucleoside RT inhibitors", "target": "reverse transcriptase"},
    "INI":   {"pattern": "INI.*.csv",   "protein": "IN",  "gene": "Pol", "subunit": "Integrase",
              "label": "Integrase strand-transfer inhibitors", "target": "integrase"},
    "CAI":   {"pattern": "CAI.*.csv",   "protein": "CA",  "gene": "Gag", "subunit": "p24",
              "label": "Capsid inhibitors", "target": "capsid"},
}

SUBTYPES = ["A", "B", "C", "D", "F", "G", "CRF01_AE", "CRF02_AG"]

# Stanford codes an in-frame codon deletion as "~" and an insertion as "#". These are real
# and clinically important NRTI mechanisms -- the T69 insertion complex and the 67/69/70
# deletions -- but they are NOT amino-acid substitutions, so they have no mu_aa.
#
# They are kept in the catalog, because an NRTI list without them would be wrong, and they
# are given no rate, because the model does not estimate in-frame codon indel rates. What
# it does estimate is a per-nucleotide-position insertion and deletion rate with NO
# nucleotide context. A single-nucleotide indel frameshifts and is a different event from
# the in-frame codon indel Stanford is describing, so quoting one as the other would be a
# category error. The per-position published rate is carried alongside as context, clearly
# named, and these rows are excluded from the codon-sensitivity analysis -- an indel
# predictor with no context has no codon dimension to vary.
INDEL_CODES = {"~": "deletion", "#": "insertion"}

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}

# Figures 5 and 6 use freq = (all_muts + 0.5) / all_seqs. Kept identical so the portal
# and the manuscript cannot disagree; the pseudocount is named in the data dictionary.
FREQ_PSEUDOCOUNT = 0.5


def build_drm(sdrd_dir, outdir, aa_level, nt_level, reference, hxb2, verbose=True):
    """Write drm_catalog.csv and drm_taxonomy.json. Returns the catalog DataFrame."""
    sdrd_dir = pathlib.Path(sdrd_dir)
    offsets = {rec["protein"]: rec["offset"] for rec in hxb2["protein_offsets"]
               if rec["offset"] is not None}

    frames = []
    for cls, spec in CLASS_SPEC.items():
        hits = sorted(glob.glob(str(sdrd_dir / spec["pattern"])))
        if not hits:
            sys.exit(f"ERROR: no SDRD export matching {spec['pattern']} in {sdrd_dir}")
        path = pathlib.Path(hits[0])
        raw = pd.read_csv(path)
        drm = raw[raw.DRM == 1].copy()
        if drm.empty:
            sys.exit(f"ERROR: {path.name} has no DRM-flagged rows")
        offset = offsets[spec["subunit"]]
        drm["drug_class"] = cls
        drm["drug_class_label"] = spec["label"]
        drm["target"] = spec["target"]
        drm["protein"] = spec["protein"]
        drm["gene"] = spec["gene"]
        drm["gene_position"] = drm.Pos + offset
        drm["source_file"] = path.name
        frames.append(drm)
        if verbose:
            uniq = (drm.Cons + drm.Pos.astype(str) + drm.Mut).nunique()
            print(f"    {cls:6s} {path.name.split('.')[0]:9s} {len(drm):4d} rows, "
                  f"{uniq:3d} unique DRMs, positions {drm.Pos.min()}-{drm.Pos.max()}, "
                  f"Rx {sorted(drm.Rx.unique())}")

    cat = pd.concat(frames, ignore_index=True)

    # The DRM flags must not overlap between classes sharing a protein, or concatenating
    # would double-count. Checked rather than assumed.
    rt = cat[cat.protein == "RT"]
    overlap = (set(zip(rt[rt.drug_class == "NRTI"].Pos, rt[rt.drug_class == "NRTI"].Cons,
                       rt[rt.drug_class == "NRTI"].Mut)) &
               set(zip(rt[rt.drug_class == "NNRTI"].Pos, rt[rt.drug_class == "NNRTI"].Cons,
                       rt[rt.drug_class == "NNRTI"].Mut)))
    assert not overlap, f"NRTI and NNRTI DRM sets overlap at {sorted(overlap)[:5]}"

    cat = cat.rename(columns={"Pos": "protein_position", "Cons": "ref_aa", "Mut": "alt_aa",
                              "Rx": "rx"})
    cat["mutation"] = cat.protein + " " + cat.ref_aa + cat.protein_position.astype(str) + cat.alt_aa
    cat["aa_change"] = cat.ref_aa + cat.gene_position.astype(str) + cat.alt_aa
    cat["freq"] = (cat.all_muts + FREQ_PSEUDOCOUNT) / cat.all_seqs
    cat["mechanism"] = [INDEL_CODES.get(a, "substitution") for a in cat.alt_aa]

    # --- NL4-3 codon coordinates, needed by both the join and the indel handling ----
    codon_of, nt_of = {}, {}
    for gene, rec in reference["proteins"].items():
        for i, (codon, nts) in enumerate(zip(rec["codon_seq"], rec["codon_nt"])):
            codon_of[(gene, rec["first_aa"] + i)] = codon
            nt_of[(gene, rec["first_aa"] + i)] = nts

    # --- join to the published amino-acid layer -----------------------------------
    aa = aa_level.set_index(["protein", "position", "alt_aa"])
    keys = list(zip(cat.gene, cat.gene_position, cat.alt_aa))
    joined = aa.reindex(keys)
    for col in ["mu_aa", "mu_aa_lo95", "mu_aa_hi95", "min_nt_changes", "n_paths"]:
        cat[col] = joined[col].values
    cat["atlas_ref_aa"] = joined["ref_aa"].values

    # Every SUBSTITUTION must join. An unjoined substitution is a coordinate bug and must
    # never be dropped silently.
    subs = cat[cat.mechanism == "substitution"]
    missing = subs[subs.mu_aa.isna()]
    if len(missing):
        print("\n  ERROR: %d substitution rows did not join to the atlas:" % len(missing))
        print(missing[["drug_class", "mutation", "gene", "gene_position"]]
              .drop_duplicates().to_string(index=False))
        sys.exit("ERROR: every DRM substitution must join to the amino-acid layer.")

    # Indel-coded rows: no mu_aa by design. Carry the published per-position single-
    # nucleotide indel rate across the codon's three positions as context only.
    cat["rate_note"] = ""
    indels = cat[cat.mechanism != "substitution"]
    if len(indels):
        nt = nt_level[nt_level.mut_type.isin(["insertion", "deletion"])]
        nt_idx = nt.set_index(["ref_position", "mut_type"]).bio_rate_median
        lo_vals, hi_vals = [], []
        for row in cat.itertuples():
            if row.mechanism == "substitution":
                lo_vals.append(np.nan); hi_vals.append(np.nan); continue
            nts = nt_of.get((row.gene, row.gene_position))
            rates = [nt_idx.get((p, row.mechanism), np.nan) for p in (nts or [])]
            rates = [r for r in rates if r == r]
            lo_vals.append(min(rates) if rates else np.nan)
            hi_vals.append(max(rates) if rates else np.nan)
        cat["nt_indel_rate_min"] = lo_vals
        cat["nt_indel_rate_max"] = hi_vals
        cat.loc[cat.mechanism != "substitution", "rate_note"] = (
            "In-frame codon indel. The model estimates per-nucleotide indel rates with no "
            "sequence context, not in-frame codon indel rates, so no mu_aa is reported and "
            "the per-position range shown is context only, not comparable to mu_aa.")
        if verbose:
            uniq = indels.mutation.nunique()
            print(f"\n  {uniq} indel-coded DRM(s) kept without a rate: "
                  f"{', '.join(sorted(indels.mutation.unique()))}")
            print("      These are real NRTI mechanisms. They have no amino-acid mutability")
            print("      because they are not substitutions, and the per-nucleotide indel")
            print("      rate is a different event (a frameshift), so it is not quoted as one.")

    # --- whether NL4-3 carries the Stanford consensus residue -----------------------
    # This matters scientifically: where NL4-3 differs from subtype-B consensus at a DRM
    # position, the "wild-type" codon a clinician has in mind is not the codon the atlas
    # used, and the reported rate is the rate FROM NL4-3's residue. Flagged, not hidden.
    cat["nl43_codon"] = [codon_of.get((g, p), "") for g, p in zip(cat.gene, cat.gene_position)]
    cat["nl43_nt_start"] = [nt_of.get((g, p), [None])[0] for g, p in zip(cat.gene, cat.gene_position)]
    # Indel rows have no atlas join, so translate their NL4-3 codon directly. Otherwise
    # they would show up as spurious consensus mismatches and dilute a flag that has to
    # mean exactly one thing.
    cat["atlas_ref_aa"] = [a if isinstance(a, str) else GENETIC_CODE.get(c, "")
                           for a, c in zip(cat.atlas_ref_aa, cat.nl43_codon)]
    cat["nl43_matches_consensus"] = cat.atlas_ref_aa == cat.ref_aa

    mismatched = cat[~cat.nl43_matches_consensus][["mutation", "ref_aa", "atlas_ref_aa"]].drop_duplicates()
    if verbose and len(mismatched):
        print(f"\n  NL4-3 differs from the Stanford consensus residue at "
              f"{len(mismatched)} DRM change(s):")
        for r in mismatched.itertuples():
            print(f"      {r.mutation}: Stanford consensus {r.ref_aa}, NL4-3 {r.atlas_ref_aa}")
        print("      These rows are flagged nl43_matches_consensus=False. The reported rate is")
        print("      the rate FROM NL4-3's residue, which is not the clinician's wild type.")

    cols = (["drug_class", "drug_class_label", "target", "rx", "mechanism",
             "protein", "protein_position",
             "ref_aa", "alt_aa", "mutation", "gene", "gene_position", "aa_change",
             "nl43_codon", "nl43_nt_start", "atlas_ref_aa", "nl43_matches_consensus",
             "mu_aa", "mu_aa_lo95", "mu_aa_hi95", "min_nt_changes", "n_paths",
             "all_seqs", "all_muts", "freq", "min_prev", "max_prev", "ratio_prev",
             "nt_indel_rate_min", "nt_indel_rate_max", "rate_note"]
            + [f"{s}_{k}" for s in SUBTYPES for k in ("seqs", "muts")]
            + ["source_file"])
    cols = [c for c in cols if c in cat.columns]
    cat = cat[cols].sort_values(["drug_class", "protein_position", "alt_aa", "rx"])
    cat.to_csv(outdir / "drm_catalog.csv", index=False)

    # --- the selector tree the page is driven by -----------------------------------
    taxonomy = {"note": ("Selector tree for the DRM view. 'classes' are Stanford drug classes; "
                         "'strata' are the treatment strata present in that class's export, "
                         "which mix class-level labels (PI, NRTI, NNRTI, INI, CAI) with "
                         "individual drugs or regimens (ATV, LPV, EFV, DOR, RPV, DTG, CAB, "
                         "ABC+XTC, AZT+XTC, TFV+XTC). 'naive' is treatment-naive and is the "
                         "comparison baseline, not a treatment."),
                "freq_pseudocount": FREQ_PSEUDOCOUNT,
                "subtypes": SUBTYPES,
                "classes": []}
    for cls, spec in CLASS_SPEC.items():
        sub = cat[cat.drug_class == cls]
        strata = sorted(sub.rx.unique())
        taxonomy["classes"].append({
            "class": cls,
            "label": spec["label"],
            "target": spec["target"],
            "protein": spec["protein"],
            "gene": spec["gene"],
            "n_mutations": int(sub.mutation.nunique()),
            "position_range": [int(sub.protein_position.min()), int(sub.protein_position.max())],
            "strata": [{"rx": s, "kind": "naive" if s == "naive"
                        else ("class" if s == cls or s in ("NRTI", "NNRTI") else "drug"),
                        "n_rows": int((sub.rx == s).sum())} for s in strata],
            "source_file": sub.source_file.iloc[0],
        })
    (outdir / "drm_taxonomy.json").write_text(json.dumps(taxonomy, indent=2) + "\n")

    if verbose:
        print(f"\n  catalog: {len(cat)} rows, {cat.mutation.nunique()} unique DRMs across "
              f"{cat.drug_class.nunique()} classes")
        print("  one-nucleotide-accessible: %d of %d unique substitution DRMs"
              % (cat[cat.min_nt_changes == 1].mutation.nunique(),
                 cat[cat.mechanism == "substitution"].mutation.nunique()))
    return cat


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sdrd", required=True, help="directory holding the SDRD class exports")
    ap.add_argument("--version", default="v0.1")
    args = ap.parse_args()
    repo = pathlib.Path(__file__).resolve().parent.parent
    outdir = repo / "data" / args.version
    aa = pd.read_csv(outdir / "aa_level_supply.csv")
    nt = pd.read_csv(outdir / "nt_level_supply.csv")
    reference = json.loads((outdir / "nl43_reference.json").read_text())
    hxb2 = json.loads((outdir / "hxb2_align.json").read_text())
    build_drm(args.sdrd, outdir, aa, nt, reference, hxb2)
    print("\nNOTE: run scripts/build_web_assets.py to refresh the manifest and releases.json.")


if __name__ == "__main__":
    main()
