#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
compare_roundtrip.py -- judge the output of test/roundtrip_nl43.mjs.

    node test/roundtrip_nl43.mjs --version v0.1
    python test/compare_roundtrip.py --version v0.1

Compares the engine's output for the NL4-3 reference against the published atlas that
ships in data/<version>/, which makes this test self-contained for anyone who clones the
repository. Point --src at the production directory to compare against the original
tables instead; the shipped tables are the same numbers rounded to 6 significant
figures, which contributes at most ~1e-6 relative and is four orders of magnitude below
the tolerance.

WHAT MUST MATCH EXACTLY (not to a tolerance -- these are classifications, not estimates)
    reference nucleotide at every position
    3' context assignment, including the G->A 2-mer and its OTHER pooling
    min_nt_changes for every amino-acid change
    n_paths for every amino-acid change
    n_nt_changes and the mutation-type path string for every codon path

WHAT MUST MATCH TO A TOLERANCE
    nucleotide bio_rate and its 95% bounds      -- expected <=0.83%, bug if worse than 1%
    mu_path and its bounds                      -- tolerance established by this test
    mu_aa and its bounds                        -- tolerance established by this test

Exits non-zero on any failure. Writes test/out/roundtrip_report.json.
"""
import argparse, json, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--version", default="v0.1")
ap.add_argument("--tol", type=float, default=0.01,
                help="relative tolerance for the NUCLEOTIDE layer (default 1%%). Path and "
                     "amino-acid layers are judged at (1+tol)^k - 1 for a k-substitution "
                     "path, because a path rate is a product of k site rates.")
ap.add_argument("--src", default=None, help="compare against this production directory instead")
args = ap.parse_args()

OUT = HERE / "out"
DATA = REPO / "data" / args.version
TOL = args.tol

if not (OUT / "engine_nt.csv").exists():
    sys.exit(f"ERROR: {OUT}/engine_nt.csv not found. Run `node test/roundtrip_nl43.mjs` first.")

# --------------------------------------------------------------------------------
# load published side
# --------------------------------------------------------------------------------
if args.src:
    src = pathlib.Path(args.src)
    pub_nt = pd.read_csv(src / "FINAL_mutation_supply_dataset.csv")
    cap = src / "comprehensive_analysis_mra_product"
    pub_aa = pd.read_csv(cap / "aa_level_supply_table.csv")
    pub_paths = pd.read_csv(cap / "aa_mutation_path_table.csv")
    source_label = f"production tables at {src}"
else:
    pub_nt = pd.read_csv(DATA / "nt_level_supply.csv")
    pub_aa = pd.read_csv(DATA / "aa_level_supply.csv")
    pub_paths = pd.concat(
        [pd.read_csv(p).assign(protein=p.stem) for p in sorted((DATA / "aa_paths").glob("*.csv"))],
        ignore_index=True)
    source_label = f"shipped assets in {DATA}"

eng_nt = pd.read_csv(OUT / "engine_nt.csv")
eng_aa = pd.read_csv(OUT / "engine_aa.csv")
eng_paths = pd.read_csv(OUT / "engine_paths.csv")

print(f"NL4-3 round-trip comparison, atlas {args.version}")
print(f"published side: {source_label}")
print(f"tolerance: {TOL*100:g}% relative\n")

failures = []
report = {"version": args.version, "tolerance": TOL, "published_source": source_label, "layers": {}}


def reldev(a, b):
    """Relative deviation |a-b|/|b|, ignoring rows where the published value is 0."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = b != 0
    out = np.full(a.shape, np.nan)
    out[ok] = np.abs(a[ok] - b[ok]) / np.abs(b[ok])
    return out


def check_exact(name, left, right, layer):
    """A classification that must agree row for row."""
    mismatch = int((left.values != right.values).sum())
    status = "OK" if mismatch == 0 else "FAIL"
    print(f"  {status:4s} exact match {name}: {len(left)-mismatch}/{len(left)}")
    layer[name] = {"n": int(len(left)), "mismatches": mismatch}
    if mismatch:
        failures.append(f"{name}: {mismatch} rows disagree")
        bad = left.index[(left.values != right.values)][:5]
        for i in bad:
            print(f"         e.g. row {i}: engine={left.loc[i]!r} published={right.loc[i]!r}")


def step_tolerance(k):
    """Tolerance for a k-substitution path.

    A path rate is the PRODUCT of k per-site rates, so the nucleotide layer's own
    reproduction residual compounds multiplicatively: (1 + d)^k - 1. That residual is
    Monte Carlo noise in the published posterior summaries, not model error and not an
    engine error -- test/out/roundtrip_report.json records the measurement that
    established it. Holding a 3-step path to the 1-step tolerance would be demanding
    that the engine be more accurate than its inputs.
    """
    return (1.0 + TOL) ** k - 1.0


def check_numeric(name, eng_col, pub_col, layer, groups=None, group_tol=None):
    """Compare two columns.

    groups    : {label: boolean mask} to report separately
    group_tol : {label: tolerance} -- when given, pass/fail is judged per group at that
                group's own tolerance instead of the flat TOL.
    """
    d = reldev(eng_col, pub_col)
    finite = d[np.isfinite(d)]
    worst = float(np.nanmax(finite)) if finite.size else float("nan")
    p99 = float(np.nanpercentile(finite, 99)) if finite.size else float("nan")
    entry = {"n": int(finite.size), "max_rel_dev": worst, "p99_rel_dev": p99}

    if group_tol is None:
        frac_in = float((finite <= TOL).mean()) if finite.size else float("nan")
        status = "OK" if worst <= TOL else "FAIL"
        print(f"  {status:4s} {name}: max {worst*100:.4f}%, p99 {p99*100:.4f}%, "
              f"{frac_in*100:.3f}% of {finite.size} rows within {TOL*100:g}%")
        entry["tolerance"] = TOL
        entry["fraction_within_tol"] = frac_in
        if worst > TOL:
            failures.append(f"{name}: max relative deviation {worst*100:.3f}% exceeds {TOL*100:g}%")
    else:
        print(f"       {name}: max {worst*100:.4f}%, p99 {p99*100:.4f}% overall")
        entry["by_group"] = {}
        for g, idx in groups.items():
            dg = d[idx][np.isfinite(d[idx])]
            if not dg.size:
                continue
            tol_g = group_tol[g]
            mx = float(dg.max())
            status = "OK" if mx <= tol_g else "FAIL"
            frac = float((dg <= tol_g).mean())
            print(f"  {status:4s}   {g}-step (n={dg.size}): max {mx*100:.4f}%, "
                  f"tolerance {tol_g*100:.3f}%, {frac*100:.3f}% within")
            entry["by_group"][str(g)] = {"n": int(dg.size), "max_rel_dev": mx,
                                         "p99_rel_dev": float(np.percentile(dg, 99)),
                                         "tolerance": tol_g, "fraction_within_tol": frac}
            if mx > tol_g:
                failures.append(f"{name} at {g} substitution(s): max relative deviation "
                                f"{mx*100:.3f}% exceeds {tol_g*100:.3f}%")
        layer[name] = entry
        return d

    if groups is not None:
        entry["by_group"] = {}
        for g, idx in groups.items():
            dg = d[idx][np.isfinite(d[idx])]
            if dg.size:
                entry["by_group"][str(g)] = {"n": int(dg.size), "max_rel_dev": float(dg.max())}
    layer[name] = entry
    return d



# --------------------------------------------------------------------------------
# 1. nucleotide layer
# --------------------------------------------------------------------------------
print("[1/3] nucleotide layer")
layer = report["layers"]["nucleotide"] = {}
pub_sub = pub_nt[~pub_nt.mut_type.isin(["insertion", "deletion"])].copy()
j = eng_nt.merge(pub_sub, on=["ref_position", "mut_type"], how="outer",
                 suffixes=("_eng", "_pub"), indicator=True)
only_eng = int((j._merge == "left_only").sum())
only_pub = int((j._merge == "right_only").sum())
print(f"  rows: engine {len(eng_nt)}, published {len(pub_sub)}, "
      f"engine-only {only_eng}, published-only {only_pub}")
layer["row_counts"] = {"engine": int(len(eng_nt)), "published": int(len(pub_sub)),
                       "engine_only": only_eng, "published_only": only_pub}
if only_eng or only_pub:
    failures.append(f"nucleotide layer: {only_eng} engine-only and {only_pub} published-only rows")
    print(j[j._merge != "both"].head(5).to_string())

m = j[j._merge == "both"].reset_index(drop=True)
check_exact("ref_nucleotide", m.ref_nucleotide_eng, m.ref_nucleotide_pub, layer)
check_exact("context_3p", m.context_3p_eng.fillna(""), m.context_3p_pub.fillna(""), layer)
check_numeric("bio_rate", m.rate, m.bio_rate_median, layer)
check_numeric("bio_rate_lo95", m.lo95, m.bio_rate_lo95, layer)
check_numeric("bio_rate_hi95", m.hi95, m.bio_rate_hi95, layer)
# Every row of an identity round-trip must have used the published interval, not the
# type-context approximation: if it did not, the engine is not finding the lookup.
approx = int((m.interval_source != "published").sum())
print(f"  {'OK' if approx == 0 else 'FAIL'}  interval source: "
      f"{len(m)-approx}/{len(m)} rows used the published per-site interval")
layer["interval_source_approximated"] = approx
if approx:
    failures.append(f"nucleotide layer: {approx} rows fell back to an approximate interval")

# --------------------------------------------------------------------------------
# 2. codon path layer
# --------------------------------------------------------------------------------
print("\n[2/3] codon path layer")
layer = report["layers"]["paths"] = {}
j = eng_paths.merge(pub_paths, on=["protein", "position", "mutant_codon"], how="outer",
                    suffixes=("_eng", "_pub"), indicator=True)
only_eng = int((j._merge == "left_only").sum()); only_pub = int((j._merge == "right_only").sum())
print(f"  rows: engine {len(eng_paths)}, published {len(pub_paths)}, "
      f"engine-only {only_eng}, published-only {only_pub}")
layer["row_counts"] = {"engine": int(len(eng_paths)), "published": int(len(pub_paths)),
                       "engine_only": only_eng, "published_only": only_pub}
if only_eng or only_pub:
    failures.append(f"path layer: {only_eng} engine-only and {only_pub} published-only rows")
m = j[j._merge == "both"].reset_index(drop=True)
check_exact("codon", m.codon_eng, m.codon_pub, layer)
check_exact("n_nt_changes", m.n_nt_changes_eng, m.n_nt_changes_pub, layer)
check_exact("mutation_type_path", m.mutation_type_path_eng, m.mutation_type_path_pub, layer)
groups = {k: (m.n_nt_changes_pub == k).values for k in sorted(m.n_nt_changes_pub.unique())}
gtol = {k: step_tolerance(k) for k in groups}
check_numeric("mu_path", m.mu_path_eng, m.mu_path_pub, layer, groups, gtol)
check_numeric("mu_path_lo95", m.mu_path_lo95_eng, m.mu_path_lo95_pub, layer, groups, gtol)
check_numeric("mu_path_hi95", m.mu_path_hi95_eng, m.mu_path_hi95_pub, layer, groups, gtol)

# --------------------------------------------------------------------------------
# 3. amino-acid layer
# --------------------------------------------------------------------------------
print("\n[3/3] amino-acid layer")
layer = report["layers"]["amino_acid"] = {}
j = eng_aa.merge(pub_aa, on=["protein", "position", "alt_aa"], how="outer",
                 suffixes=("_eng", "_pub"), indicator=True)
only_eng = int((j._merge == "left_only").sum()); only_pub = int((j._merge == "right_only").sum())
print(f"  rows: engine {len(eng_aa)}, published {len(pub_aa)}, "
      f"engine-only {only_eng}, published-only {only_pub}")
layer["row_counts"] = {"engine": int(len(eng_aa)), "published": int(len(pub_aa)),
                       "engine_only": only_eng, "published_only": only_pub}
if only_eng or only_pub:
    failures.append(f"amino-acid layer: {only_eng} engine-only and {only_pub} published-only rows")
    print(j[j._merge != "both"].head(5)[["protein", "position", "alt_aa", "_merge"]].to_string())
m = j[j._merge == "both"].reset_index(drop=True)
check_exact("ref_aa", m.ref_aa_eng, m.ref_aa_pub, layer)
check_exact("aa_change", m.aa_change_eng, m.aa_change_pub, layer)
check_exact("synonymous", m.synonymous_eng.astype(bool), m.synonymous_pub.astype(bool), layer)
check_exact("min_nt_changes", m.min_nt_changes_eng, m.min_nt_changes_pub, layer)
check_exact("n_paths", m.n_paths_eng, m.n_paths_pub, layer)
# mu_aa sums paths of 1, 2 and 3 steps, so a row is judged at the tolerance of the
# longest path that can contribute to it -- 3 steps for every row, in practice, since
# every amino-acid change has 3-step routes. Grouping by min_nt_changes is for
# reporting; the tolerance applied is the 3-step bound throughout.
groups = {k: (m.min_nt_changes_pub == k).values for k in sorted(m.min_nt_changes_pub.unique())}
gtol = {k: step_tolerance(3) for k in groups}
check_numeric("mu_aa", m.mu_aa_eng, m.mu_aa_pub, layer, groups, gtol)
check_numeric("mu_aa_lo95", m.mu_aa_lo95_eng, m.mu_aa_lo95_pub, layer, groups, gtol)
check_numeric("mu_aa_hi95", m.mu_aa_hi95_eng, m.mu_aa_hi95_pub, layer, groups, gtol)

# --------------------------------------------------------------------------------
report["passed"] = not failures
report["failures"] = failures
(OUT / "roundtrip_report.json").write_text(json.dumps(report, indent=2) + "\n")

print("\n" + "=" * 70)
if failures:
    print(f"ROUND-TRIP FAILED -- {len(failures)} problem(s):")
    for f in failures:
        print("  - " + f)
    print("\nAny discrepancy is a bug in the engine, not a modelling choice.")
    print(f"Full detail: {OUT}/roundtrip_report.json")
    sys.exit(1)
print("ROUND-TRIP PASSED")
print(f"Full detail: {OUT}/roundtrip_report.json")
