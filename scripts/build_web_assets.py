#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
build_web_assets.py -- turn the released MRA tables into the portal's static assets.

The portal must be rebuildable from the source tables, so this script is the only
thing that writes data/<version>/. It never modifies an existing release directory.

Inputs (--src, the production output directory):
    posterior_summary.csv                       )  consumed by
    mutation_supply_site_mutation_table.csv      )  extract_mra_web_model.py
    FINAL_mutation_supply_dataset.csv                 published nucleotide layer
    comprehensive_analysis_mra_product/
        aa_level_supply_table.csv                     published amino-acid layer
        aa_mutation_path_table.csv                    codon -> codon paths

Inputs (--ref-root, the repository root of the paper):
    data/pipeline/00-consensus_reference/NL43.fa
    data/pipeline/00-bed_files/NL43_annotations_by_position.csv
    mra/ExcelWorkbook/Extended Data - Modeled AA Mutation Supply Dataset V4.0.xlsm
        sheet `alignment`   protein offsets within genes, amino-acid units
        sheet `hxb2_align`  HXB2 <-> NL4-3 amino-acid mapping, Env and Vpu only

Outputs (data/<version>/):
    mra_web_model.json      fitted model, version stamped inside
    nt_level_supply.csv     published per-site rates, positions 1-8213, incl. indels
    drm_catalog.csv         drug-resistance mutations, rebuilt from the SDRD class exports
    drm_taxonomy.json       the class/drug selector tree the DRM view is driven by
    interval_ratios.json    95% interval ratios, so a custom rate can carry an interval
    aa_level_supply.csv     published mu_aa per (protein, position, alt_aa)
    aa_paths/<Protein>.csv  codon -> codon path shards, loaded lazily
    nl43_reference.json     sequence, per-position protein assignment, codon frames
    hxb2_align.json         HXB2 lookup (Env, Vpu) + protein offsets
    citation.json           machine-readable citation block
    manifest.json           SHA-256 + raw/gzipped size per asset

and rewrites data/releases.json from the release directories that exist.

WHY SOME THINGS ARE READ RATHER THAN DERIVED
    Codon frames come from aa_mutation_path_table.codon_position_{1,2,3}, not from
    re-deriving reading frames. HIV-1 genes overlap, and tat and rev are spliced: one
    rev codon spans a 2,325 nt splice junction. The published table already encodes
    that correctly, so it is the authority.
"""
import argparse, csv, gzip, hashlib, json, pathlib, shutil, subprocess, sys, datetime
import pandas as pd
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from build_drm_assets import build_drm
from build_bnab_assets import build_bnab
from build_lanl_batch_demo import build_lanl_demo

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent

# Significant figures kept for rates in the shipped CSVs. Six sig figs is ~1e-6
# relative, four orders of magnitude finer than the 1% round-trip tolerance, and it
# roughly halves the payload versus full float64 repr. Documented in the data dictionary.
SIGFIGS = 6

MODELLED_FIRST, MODELLED_LAST = 1, 8213


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def sig(x, n=SIGFIGS):
    """Format a float to n significant figures, or empty string for NaN."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return f"{x:.{n}g}"


def sha256_and_sizes(path):
    raw = path.read_bytes()
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "gzip_bytes": len(gzip.compress(raw, 9)),
    }


def require(path, what):
    p = pathlib.Path(path)
    if not p.exists():
        sys.exit(f"ERROR: missing {what}: {p}")
    return p


def read_fasta_single(path):
    seq = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith(">"):
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


# ---------------------------------------------------------------------------
# asset builders
# ---------------------------------------------------------------------------

def build_model_bundle(src, outdir, version, doi, release_date):
    """Run extract_mra_web_model.py, which is the single authority on model extraction.

    Its internal assertions (shared spline must agree across the 11 non-G->A types at
    each position; every spline projection residual under tolerance) are a hard stop.
    """
    script = require(HERE / "extract_mra_web_model.py", "extract_mra_web_model.py")
    out = outdir / "mra_web_model.json"
    cmd = [sys.executable, str(script), "--version", version,
           "--release-date", release_date, "--src", str(src), "--out", str(out)]
    if doi:
        cmd += ["--doi", doi]
    print("  $", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout); print(res.stderr, file=sys.stderr)
        sys.exit("ERROR: model extraction failed -- treat an assertion failure as a hard stop.")
    print("  " + res.stdout.strip().replace("\n", "\n  "))

    bundle = json.loads(out.read_text())
    assert bundle["version"] == version, "extractor did not stamp the requested version"
    assert bundle["modelled_region"]["first_position"] == MODELLED_FIRST
    assert bundle["modelled_region"]["last_position"] == MODELLED_LAST
    for key, curve in bundle["spline_by_position"].items():
        if isinstance(curve, list):
            assert len(curve) == MODELLED_LAST, f"{key} spline has {len(curve)} positions"
    return bundle


def build_nt_layer(src, outdir):
    """Published nucleotide layer: one row per (position, mutation type).

    bio_rate_median is the column the amino-acid layer is built from and the column the
    model bundle reproduces, so it leads. bio_rate_mean is carried for reference because
    it is what the paper's summary statistics use, with a note in the data dictionary
    that the two are not interchangeable (mean of a rate != exp(mean of log rate)).
    """
    df = pd.read_csv(src / "FINAL_mutation_supply_dataset.csv")
    before = len(df)
    df = df[(df.ref_position >= MODELLED_FIRST) & (df.ref_position <= MODELLED_LAST)]
    print(f"  nucleotide rows: {before} -> {len(df)} after restricting to {MODELLED_FIRST}-{MODELLED_LAST}")

    keep = ["ref_position", "ref_nucleotide", "mut_type", "context_3p",
            "bio_rate_median", "bio_rate_lo95", "bio_rate_hi95", "bio_rate_mean",
            "wt_count", "wt_coverage", "control_count", "control_coverage",
            "control_corrected"]
    out = df[keep].copy()
    out["context_3p"] = out.context_3p.fillna("")
    for c in ["bio_rate_median", "bio_rate_lo95", "bio_rate_hi95", "bio_rate_mean"]:
        out[c] = out[c].map(sig)
    path = outdir / "nt_level_supply.csv"
    out.to_csv(path, index=False)
    return df


def build_aa_layer(src, outdir):
    """Published amino-acid layer, and the codon-path table sharded by protein.

    log10_* columns are dropped: they are floored at -10 in the source, which is a
    plotting artifact rather than a published quantity. The portal takes log10 of mu_aa
    at display time instead.
    """
    cap = src / "comprehensive_analysis_mra_product"
    aa = pd.read_csv(cap / "aa_level_supply_table.csv")
    keep = ["protein", "position", "aa_change", "ref_aa", "alt_aa", "synonymous",
            "mu_aa", "mu_aa_lo95", "mu_aa_hi95", "min_nt_changes", "n_paths"]
    out = aa[keep].copy()
    for c in ["mu_aa", "mu_aa_lo95", "mu_aa_hi95"]:
        out[c] = out[c].map(sig)
    out.to_csv(outdir / "aa_level_supply.csv", index=False)
    print(f"  amino-acid rows: {len(out)}")

    paths = pd.read_csv(cap / "aa_mutation_path_table.csv")
    shard_dir = outdir / "aa_paths"
    shard_dir.mkdir(exist_ok=True)
    pkeep = ["position", "codon", "mutant_codon", "ref_aa", "alt_aa", "aa_change",
             "codon_position_1", "codon_position_2", "codon_position_3",
             "n_nt_changes", "mutation_type_path",
             "mu_path", "mu_path_lo95", "mu_path_hi95", "synonymous"]
    for protein, grp in paths.groupby("protein"):
        g = grp[pkeep].copy()
        for c in ["mu_path", "mu_path_lo95", "mu_path_hi95"]:
            g[c] = g[c].map(sig)
        g.to_csv(shard_dir / f"{protein}.csv", index=False)
    print(f"  path shards: {sorted(p.name for p in shard_dir.glob('*.csv'))}")
    return aa, paths


def build_interval_ratios(nt, outdir):
    """interval_ratios.json -- how to put a 95% interval on a CUSTOM rate.

    THE PROBLEM
        mra_web_model.json carries posterior *means* of the coefficients and nothing
        about posterior width, so the model bundle alone cannot produce an interval.
        But the brief requires that no point estimate is ever shown alone, and the
        sequence tool has to work on the bundle plus small assets only.

    THE SOLUTION
        Intervals are carried as multiplicative ratios to the median:
            lo_ratio = log(bio_rate_lo95 / bio_rate_median)
            hi_ratio = log(bio_rate_hi95 / bio_rate_median)
        and applied to the engine's point estimate as median * exp(ratio).

        Two sources, in priority order:

        1. `by_position` -- the published ratio for that exact (position, mutation
           type). Available for the three substitutions NL4-3 actually carries at each
           position. Using it means a submitted sequence that matches NL4-3 at a site
           gets the published interval back to within the bundle's own ~0.8%
           reproduction tolerance. This is what makes the NL4-3 round-trip exact.

        2. `by_type_context` -- the mean log-ratio for that (mutation type, 3' base),
           used when a submitted sequence has a base NL4-3 does not have at that
           position, so the published table has no entry. Within a type the log-ratios
           are tight (SD 0.04-0.10, i.e. the bound moves by about +/-10%), so this is a
           usable approximation -- but it IS an approximation and the engine flags every
           row that used it as interval_source="type-context".

    Indels are excluded: their log-ratios have SD 0.15-0.47, far too variable to
    summarise this way, and indels carry no sequence-specific information so they are
    not part of custom output anyway.
    """
    sub = nt[~nt.mut_type.isin(["insertion", "deletion"])].copy()
    sub["lo_ratio"] = np.log(sub.bio_rate_lo95 / sub.bio_rate_median)
    sub["hi_ratio"] = np.log(sub.bio_rate_hi95 / sub.bio_rate_median)

    by_tc = {}
    for (mt, ctx), g in sub.groupby(["mut_type", sub.context_3p.fillna("")]):
        by_tc[f"{mt}|{ctx}"] = {
            "lo_ratio": round(float(g.lo_ratio.mean()), 5),
            "hi_ratio": round(float(g.hi_ratio.mean()), 5),
            "lo_sd": round(float(g.lo_ratio.std()), 5),
            "hi_sd": round(float(g.hi_ratio.std()), 5),
            "n": int(len(g)),
        }
    # Type-level means too, for the rare case of a 3' context not seen for that type
    # (the pooled cells) -- better a slightly wrong width than no interval at all.
    by_type = {}
    for mt, g in sub.groupby("mut_type"):
        by_type[mt] = {"lo_ratio": round(float(g.lo_ratio.mean()), 5),
                       "hi_ratio": round(float(g.hi_ratio.mean()), 5),
                       "n": int(len(g))}

    by_pos = {}
    for mt, g in sub.groupby("mut_type"):
        g = g.sort_values("ref_position")
        by_pos[mt] = {
            "positions": [int(p) for p in g.ref_position],
            "lo_ratio": [round(float(v), 4) for v in g.lo_ratio],
            "hi_ratio": [round(float(v), 4) for v in g.hi_ratio],
        }

    doc = {
        "note": ("Multiplicative 95% interval ratios to the median rate, as natural logs. "
                 "Apply as lo95 = rate * exp(lo_ratio), hi95 = rate * exp(hi_ratio). "
                 "by_position is the published per-site interval and is exact; "
                 "by_type_context is a mean over sites and is an approximation that the "
                 "engine flags. Substitutions only -- indel interval widths are too "
                 "variable to summarise and indels are not part of custom output."),
        "by_position": by_pos,
        "by_type_context": by_tc,
        "by_type": by_type,
    }
    (outdir / "interval_ratios.json").write_text(json.dumps(doc, separators=(",", ":")))
    n_pos = sum(len(v["positions"]) for v in by_pos.values())
    worst = max(max(v["lo_sd"], v["hi_sd"]) for v in by_tc.values())
    print(f"  interval ratios: {n_pos} published site entries, {len(by_tc)} type x context cells; "
          f"worst within-cell SD of the log ratio {worst:.3f}")
    return doc


def build_reference(ref_root, paths, outdir):
    """nl43_reference.json -- sequence, per-position protein assignment, codon frames.

    Codon frames are taken from the published path table so that gene overlap and the
    tat/rev splice junctions are correct by construction rather than re-derived.
    """
    seq = read_fasta_single(require(
        ref_root / "data/pipeline/00-consensus_reference/NL43.fa", "NL4-3 reference"))
    ann = pd.read_csv(require(
        ref_root / "data/pipeline/00-bed_files/NL43_annotations_by_position.csv",
        "per-position annotations"))
    assert len(seq) == 8618, f"reference is {len(seq)} bp, expected 8618"
    assert len(ann) == 8618, f"annotation table has {len(ann)} rows, expected 8618"
    mismatches = int((ann.base.str.upper().values != np.array(list(seq))).sum())
    assert mismatches == 0, f"{mismatches} positions disagree between NL43.fa and the annotation table"

    # Categorical columns stored as a code dictionary plus an index array: an 8,618-entry
    # array of small integers instead of 8,618 repeated strings.
    ann_levels = sorted(ann.annotation_v4.unique())
    bio_levels = sorted(ann.bio_class.unique())
    ann_code = {v: i for i, v in enumerate(ann_levels)}
    bio_code = {v: i for i, v in enumerate(bio_levels)}

    codons = (paths.drop_duplicates(["protein", "position"])
                   .sort_values(["protein", "position"]))
    codon_records = {}
    for protein, grp in codons.groupby("protein"):
        codon_records[protein] = {
            "first_aa": int(grp.position.min()),
            "last_aa": int(grp.position.max()),
            "n_codons": int(len(grp)),
            # [[nt1, nt2, nt3], ...] in AA order; NL4-3 reference positions
            "codon_nt": [[int(a), int(b), int(c)] for a, b, c in
                         zip(grp.codon_position_1, grp.codon_position_2, grp.codon_position_3)],
            "codon_seq": list(grp.codon),
        }
        # consistency: the codon string must match the reference sequence at those positions
        for (a, b, c), cod in zip(codon_records[protein]["codon_nt"], codon_records[protein]["codon_seq"]):
            got = seq[a - 1] + seq[b - 1] + seq[c - 1]
            assert got == cod, f"{protein} codon at {a},{b},{c}: reference says {got}, table says {cod}"

    doc = {
        "reference_name": "NL4-3 (ERR-Seq consensus reference)",
        "length": len(seq),
        "position_1_is": "the A of the gag start codon",
        "modelled_region": {"first_position": MODELLED_FIRST, "last_position": MODELLED_LAST},
        "note": ("Align submitted sequences against all 8618 bp, then restrict output to the "
                 "modelled region. Trimming the reference before alignment would lose the "
                 "anchoring that keeps the reading frame correct upstream."),
        "sequence": seq,
        "annotation_levels": ann_levels,
        "bio_class_levels": bio_levels,
        "annotation_by_position": [ann_code[v] for v in ann.annotation_v4],
        "bio_class_by_position": [bio_code[v] for v in ann.bio_class],
        "proteins": codon_records,
    }
    (outdir / "nl43_reference.json").write_text(json.dumps(doc, separators=(",", ":")))
    print(f"  reference: {len(seq)} bp, {sum(v['n_codons'] for v in codon_records.values())} codons "
          f"across {len(codon_records)} proteins")
    return doc


def build_hxb2(ref_root, outdir):
    """hxb2_align.json -- HXB2 numbering lookup and protein offsets within genes.

    Only Env and Vpu need a lookup; every other protein maps identity between NL4-3 and
    HXB2 numbering. The mapping for Env and Vpu is position-by-position because of
    indels, so a constant offset would be wrong.
    """
    import openpyxl
    wb_path = require(
        ref_root / "mra/ExcelWorkbook/Extended Data - Modeled AA Mutation Supply Dataset V4.0.xlsm",
        "the atlas workbook")
    wb = openpyxl.load_workbook(wb_path, read_only=True, data_only=True)

    rows = list(wb["alignment"].iter_rows(values_only=True))
    header = [str(h) for h in rows[0]]
    offsets = []
    for r in rows[1:]:
        rec = dict(zip(header, r))
        if rec.get("protein") is None:
            continue
        off = rec["offset"]
        offsets.append({
            "protein": str(rec["protein"]),
            "gene": str(rec["gene"]),
            # offset is in AMINO-ACID units: gene AA position = subprotein AA position + offset
            "offset": None if off in (None, "N/A") else int(off),
            "length": int(rec["length"]),
            "notes": None if rec.get("Notes") is None else str(rec["Notes"]),
        })

    rows = list(wb["hxb2_align"].iter_rows(values_only=True))
    header = [str(h) for h in rows[0]]
    mapping = {}
    for r in rows[1:]:
        rec = dict(zip(header, r))
        if rec.get("protein") is None:
            continue
        prot = str(rec["protein"])
        nl = rec["nl43_pos"]
        mapping.setdefault(prot, {})[str(int(rec["hxb2_pos"]))] = (
            None if nl in (None, "-", "") else int(nl))
    wb.close()

    doc = {
        "note": ("Amino-acid level only. Whole-genome HXB2 NUCLEOTIDE coordinates are not "
                 "derivable from these tables: the atlas reference starts at the gag start "
                 "codon (HXB2 790) and there are indels downstream, so a constant offset "
                 "would be wrong. Report HXB2 at the protein level."),
        "lookup_required_for": sorted(mapping),
        "identity_elsewhere": True,
        "hxb2_to_nl43": mapping,
        "protein_offsets": offsets,
    }
    (outdir / "hxb2_align.json").write_text(json.dumps(doc, separators=(",", ":")))
    print("  hxb2 lookup: " + ", ".join(f"{k} ({len(v)} rows)" for k, v in sorted(mapping.items()))
          + f"; {len(offsets)} protein offsets")
    return doc


def build_citation(outdir, version, doi, release_date, model):
    doc = {
        "atlas": "HIV-1 Mutation Rate Atlas (MRA)",
        "version": version,
        "release_date": release_date,
        "doi": doi,
        "model": model,
        "paper": "[Paper citation to be completed on publication]",
        "code_licence": "MIT",
        "data_licence": "CC-BY-4.0",
        "how_to_cite": ("Cite the atlas and state the version used. Atlas data are CC-BY-4.0: "
                        "attribution is required, not optional."),
    }
    (outdir / "citation.json").write_text(json.dumps(doc, indent=2) + "\n")
    return doc


# ---------------------------------------------------------------------------
# release manifest
# ---------------------------------------------------------------------------

def write_manifest(outdir, version, doi, release_date, model):
    assets = {}
    for p in sorted(outdir.rglob("*")):
        if p.is_dir() or p.name == "manifest.json":
            continue
        assets[str(p.relative_to(outdir))] = sha256_and_sizes(p)
    total_gz = sum(a["gzip_bytes"] for a in assets.values())
    bundle_gz = assets["mra_web_model.json"]["gzip_bytes"]
    manifest = {
        "version": version, "release_date": release_date, "doi": doi, "model": model,
        "generated_by": "scripts/build_web_assets.py",
        "sigfigs_for_rates": SIGFIGS,
        "total_gzip_bytes": total_gz,
        "assets": assets,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    # Payload budget from the brief: bundle ~170 KB gzipped, browse tables 2.6-7.5 MB
    # gzipped. The sequence tool must be usable on the bundle alone, so the bundle is
    # checked separately and hard.
    print(f"\n  PAYLOAD: model bundle {bundle_gz/1024:.0f} KB gzipped; "
          f"all assets {total_gz/1024/1024:.2f} MB gzipped")
    if bundle_gz > 400 * 1024:
        sys.exit(f"ERROR: model bundle is {bundle_gz/1024:.0f} KB gzipped, over the 400 KB ceiling. "
                 "The sequence tool must be usable on the bundle alone.")
    if total_gz > 12 * 1024 * 1024:
        sys.exit(f"ERROR: total assets {total_gz/1024/1024:.1f} MB gzipped, over the 12 MB ceiling.")
    return manifest


def rewrite_releases(repo, version, doi, release_date, model, description, manifest):
    """Regenerate data/releases.json from the release directories that exist.

    Earlier entries are preserved verbatim. Old versions are never deleted: the cost of
    keeping a static directory is nil, the cost of removing one is a broken citation.
    """
    rel_path = repo / "data" / "releases.json"
    doc = json.loads(rel_path.read_text()) if rel_path.exists() else {"releases": []}
    existing = {r["version"]: r for r in doc.get("releases", [])}
    existing[version] = {
        "version": version,
        "release_date": release_date,
        "doi": doi,
        "paper": "[Paper citation to be completed on publication]",
        "model": model,
        "description": description,
        "modelled_region": {"first_position": MODELLED_FIRST, "last_position": MODELLED_LAST},
        "assets": {k: {"sha256": v["sha256"], "gzip_bytes": v["gzip_bytes"]}
                   for k, v in manifest["assets"].items()},
    }

    def sortkey(v):
        return [int(x) for x in v.lstrip("v").split(".")]

    ordered = sorted(existing.values(), key=lambda r: sortkey(r["version"]))
    doc = {
        "_comment": [
            "Stable path, never versioned. Lists every atlas release that has ever existed.",
            "Old versions are never deleted: they are static files, and the cost of removing",
            "one is a broken citation. Regenerated by scripts/build_web_assets.py.",
        ],
        "latest": ordered[-1]["version"],
        "releases": ordered,
    }
    rel_path.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"  releases.json: {[r['version'] for r in ordered]}, latest={doc['latest']}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", required=True, help="release string, e.g. v0.1")
    ap.add_argument("--doi", default=None)
    ap.add_argument("--release-date", default=datetime.date.today().isoformat())
    ap.add_argument("--description", default="Atlas release.")
    ap.add_argument("--src", required=True, help="production output directory")
    ap.add_argument("--ref-root", required=True, help="repository root holding data/pipeline and mra/")
    ap.add_argument("--sdrd", default=None,
                    help="directory of SDRD drug-class exports; omit to skip the DRM catalog")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing release directory (never do this after publication)")
    args = ap.parse_args()

    version = args.version if args.version.startswith("v") else "v" + args.version
    src = require(args.src, "source directory")
    ref_root = require(args.ref_root, "reference repository root")
    outdir = REPO / "data" / version

    if (outdir / "manifest.json").exists() and not args.force:
        sys.exit(f"ERROR: {outdir} already exists and carries a manifest. A release directory is "
                 "never modified in place; cut a new version instead. Use --force only for a "
                 "release that has not left the lab.")
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Building atlas {version} into {outdir}\n")
    print("[1/7] model bundle")
    bundle = build_model_bundle(src, outdir, version, args.doi, args.release_date)
    model = bundle["model"]
    print("\n[2/7] nucleotide layer")
    nt = build_nt_layer(src, outdir)
    print("\n[3/7] interval ratios for custom sequences")
    build_interval_ratios(nt, outdir)
    print("\n[4/7] amino-acid layer")
    aa, paths = build_aa_layer(src, outdir)
    print("\n[5/7] reference sequence and codon frames")
    build_reference(ref_root, paths, outdir)
    print("\n[6/7] HXB2 lookup and protein offsets")
    hxb2 = build_hxb2(ref_root, outdir)
    if args.sdrd:
        print("\n[6b/7] drug-resistance mutation catalog")
        reference = json.loads((outdir / "nl43_reference.json").read_text())
        nt_shipped = pd.read_csv(outdir / "nt_level_supply.csv")
        aa_shipped = pd.read_csv(outdir / "aa_level_supply.csv")
        build_drm(args.sdrd, outdir, aa_shipped, nt_shipped, reference, hxb2)
    else:
        print("\n[6b/7] drug-resistance mutation catalog SKIPPED (--sdrd not given)")
    print("\n[6c/7] provisional bNAb escape catalog and glycan sequons")
    ref_for_bnab = json.loads((outdir / "nl43_reference.json").read_text())
    build_bnab(outdir, ref_for_bnab, hxb2)
    print("\n[6d/7] batch-analysis worked example (LANL compendium, group M)")
    lanl_src = REPO / "data/source/HIV1_COM_2025_genome_DNA.fasta.gz"
    if lanl_src.exists():
        build_lanl_demo(lanl_src, outdir)
    else:
        print(f"  SKIPPED: {lanl_src.relative_to(REPO)} not present")
    # Must follow the worked example and the model bundle: it runs the batch pipeline over
    # that FASTA using the shipped bundle. Node, because it imports the app's own modules --
    # reimplementing the fold in Python would be a second copy of the thing most likely to
    # drift from what the browser actually computes.
    print("\n[6e/7] precomputed LANL batch for the landing-page supply flow")
    if (outdir / "lanl_batch_demo.fasta").exists():
        subprocess.run(["node", str(REPO / "scripts/build_lanl_flow.mjs"),
                        "--version", version], cwd=REPO, check=True)
    else:
        print("  SKIPPED: lanl_batch_demo.fasta not built")
    build_citation(outdir, version, args.doi, args.release_date, model)
    print("\n[7/7] manifest and release list")
    manifest = write_manifest(outdir, version, args.doi, args.release_date, model)
    rewrite_releases(REPO, version, args.doi, args.release_date, model, args.description, manifest)
    print("\nDone.")


if __name__ == "__main__":
    main()
