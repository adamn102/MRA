#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
extract_mra_web_model.py -- turn the fitted MRA posterior into a small static asset.

WHY THIS EXISTS
    The MRA was fitted in PyMC. The sampled posterior (idata.nc) is 9.7 GB and cannot be
    shipped to a browser. But the model's biological-rate predictor is small: a handful of
    type and context coefficients plus six position P-splines. This script recovers all of
    it from two released files and writes a single JSON of ~125 numbers.

THE MODEL (see model_metadata.json)
    non-G->A substitutions:
        log(rate) = beta_sub_type_bio[type]
                  + powered * beta_sub_context_delta[type, 3' base]
                  + spline_shared(position)
    G->A:
        log(rate) = intercept_bio_ga[pooled 3' 2-mer] + spline_ga[class](position)
        pooled: CC, CG, TC, TT -> "OTHER";  class: "AA" | "GA" | "SHARED"
    indels:
        log(rate) = intercept_bio + beta_mut_bio[type] + spline_{insertion,deletion}(position)

HOW THE SPLINES ARE RECOVERED WITHOUT idata.nc
    posterior_summary.csv carries every coefficient EXCEPT the spline shapes. But the fitted
    per-site rates in mutation_supply_site_mutation_table.csv are a deterministic function of
    those coefficients, so the spline is recoverable by subtraction:

        spline(position) = log(bio_rate_median) - beta_type - beta_context

    Use bio_rate_MEDIAN, not mean. The posterior mean of a rate is not the exponential of the
    posterior mean of the linear predictor (Jensen), but the median is, because exp is monotone.
    Using the mean introduces a bias of a few percent.

    The recovered values are then projected back onto the same patsy cubic B-spline basis the
    model used, giving 31 coefficients per spline that can be evaluated at ANY position --
    including positions where NL4-3 has no observation of that mutation class. This matters:
    6,275 of 8,213 positions are not G in NL4-3, but a submitted sequence may have a G there.

SELF-CHECK
    All 11 non-G->A substitution types share one spline, so at a position carrying several
    types the recovered value must agree across them. The script asserts this. Disagreement
    means the coefficient bookkeeping is wrong.
"""
import argparse, datetime, json, numpy as np, pandas as pd, patsy

# Version metadata is stamped INTO the bundle, not only into its directory path, so a downloaded
# file is self-identifying. The portal cross-checks the embedded version against the directory it
# loaded from and refuses to start if they disagree -- that mismatch means assets were copied
# between release directories, which would silently serve the wrong model.
_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--version", required=True,
                 help="atlas release, e.g. v1.0 -- must match the /data/<version>/ directory")
_ap.add_argument("--doi", default=None, help="DOI of the paper this release corresponds to")
_ap.add_argument("--release-date", default=datetime.date.today().isoformat())
_ap.add_argument("--src", default=".", help="directory holding the two input CSVs")
_ap.add_argument("--out", default="mra_web_model.json")
ARGS = _ap.parse_args()

SRC = ARGS.src
POS_MEAN, POS_STD = 4107.0, 2370.917730893402   # from model_metadata.json
SPLINE_DF, DEGREE = 30, 3
GA_POOLED = {"CC", "CG", "TC", "TT"}            # -> "OTHER"

ps = pd.read_csv(f"{SRC}/posterior_summary.csv", index_col=0)["mean"]
sm = pd.read_csv(f"{SRC}/mutation_supply_site_mutation_table.csv")

def indexed(prefix):
    out = {}
    for name, val in ps.items():
        if name.startswith(prefix + "["):
            key = name.split("[", 1)[1].rstrip("]")
            out[tuple(k.strip() for k in key.split(","))] = float(val)
    return out

beta_type    = {k[0]: v for k, v in indexed("beta_sub_type_bio").items()}
beta_context = indexed("beta_sub_context_delta")
ga_intercept = {k[0]: v for k, v in indexed("intercept_bio_ga").items()}
beta_mut_bio = {k[0]: v for k, v in indexed("beta_mut_bio").items()}
intercept_bio = float(ps["intercept_bio"])

# Basis built on the model's own observation grid so the knots match the fit exactly.
design = patsy.dmatrix(
    f"bs(x, df={SPLINE_DF}, degree={DEGREE}, include_intercept=False)",
    {"x": (sm.ref_position.values - POS_MEAN) / POS_STD},
    return_type="dataframe",
).design_info

def basis(positions):
    x = (np.asarray(positions, float) - POS_MEAN) / POS_STD
    return np.asarray(patsy.build_design_matrices([design], {"x": x})[0])

def project(positions, values, label, tol=0.02):
    B = basis(positions)
    theta, *_ = np.linalg.lstsq(B, values, rcond=None)
    resid = float(np.abs(values - B @ theta).max())
    assert resid < tol, f"{label}: spline projection residual {resid:.3g} exceeds {tol}"
    return theta, resid

splines, report = {}, {}

# --- shared spline for the 11 non-G->A substitutions ---------------------------
sub = sm[sm.mut_type.isin(beta_type)].copy()
sub["resid"] = (
    np.log(sub.bio_rate_median)
    - sub.mut_type.map(beta_type)
    - sub.sub_context_powered.astype(float)
    * [beta_context.get((t, str(c)), 0.0) for t, c in zip(sub.mut_type, sub.context_3p)]
)
by_pos = sub.groupby("ref_position").resid.agg(["mean", "std", "size"])
multi = by_pos[by_pos["size"] > 1]
assert multi["std"].max() < 0.05, (
    f"shared spline disagrees across mutation types (max SD {multi['std'].max():.3g}); "
    "coefficient bookkeeping is wrong"
)
splines["shared_sub"], r = project(by_pos.index.values, by_pos["mean"].values, "shared_sub")
report["shared_sub"] = dict(n_positions=int(len(by_pos)),
                            max_within_position_sd=float(multi["std"].max()), residual=r)

# --- three G->A splines --------------------------------------------------------
ga = sm[sm.mut_type == "G->A"].copy()
ga["resid"] = np.log(ga.bio_rate_median) - ga.context_3p.astype(str).map(ga_intercept)
for cls, grp in ga.groupby(ga.ga_spline_context.astype(str)):
    p = grp.groupby("ref_position").resid.mean()
    splines[f"GtoA_{cls}"], r = project(p.index.values, p.values, f"GtoA_{cls}")
    report[f"GtoA_{cls}"] = dict(n_positions=int(len(p)), residual=r)

# --- indel splines -------------------------------------------------------------
for mt, name in (("insertion", "insertion"), ("deletion", "deletion")):
    grp = sm[sm.mut_type == mt]
    if grp.empty:
        continue
    resid = np.log(grp.bio_rate_median) - intercept_bio - beta_mut_bio[mt]
    p = pd.Series(resid.values, index=grp.ref_position.values).groupby(level=0).mean()
    splines[name], r = project(p.index.values, p.values, name)
    report[name] = dict(n_positions=int(len(p)), residual=r)

lo, hi = int(sm.ref_position.min()), int(sm.ref_position.max())
grid = np.arange(lo, hi + 1)
Bg = basis(grid)

bundle = {
    "version": ARGS.version,
    "release_date": ARGS.release_date,
    "doi": ARGS.doi,
    "model": "mra_joint_control_wt_nb_ga2_sub1_ga3_sum0_bio_splines",
    "units": "biological mutation rate per site per replication (control-corrected)",
    "modelled_region": {"first_position": lo, "last_position": hi,
                        "reference": "NL43.fa, 8618 bp, position 1 = gag start codon"},
    "position_scaling": {"mean": POS_MEAN, "std": POS_STD},
    "spline_basis": {"type": "cubic B-spline (patsy bs)", "df": SPLINE_DF,
                     "degree": DEGREE, "n_basis": int(Bg.shape[1]),
                     "knots_note": "built on the model's observation grid; "
                                   "evaluate via the shipped per-position table, not by rebuilding"},
    "beta_sub_type_bio": beta_type,
    "beta_sub_context_delta": {f"{t}|{c}": v for (t, c), v in beta_context.items()},
    # Cells with fewer than 20 WT counts were pooled to type-only: the model multiplies
    # their context delta by zero. The sampled delta is therefore prior-driven and NOT
    # zero -- up to 0.0055 on the log scale, i.e. 0.55% in rate, half the round-trip
    # tolerance -- so a consumer that applied it would be quietly wrong. Derived from
    # sub_context_powered in the site table rather than read from model_metadata.json,
    # so this script still depends on only its two inputs.
    "sub_context_pooled_cells": sorted(
        f"{t}|{c}" for t, c in set(zip(
            sm.loc[(sm.sub_context_powered == 0) & sm.mut_type.isin(beta_type), "mut_type"],
            sm.loc[(sm.sub_context_powered == 0) & sm.mut_type.isin(beta_type), "context_3p"].astype(str)))
    ),
    # G->A 3' 2-mers are read from the submitted sequence and then pooled: a raw 2-mer in
    # ga_pooled_2mers maps to the intercept level "OTHER".
    "ga_context_is_2mer": True,
    "intercept_bio_ga": ga_intercept,
    "intercept_bio": intercept_bio,
    "beta_mut_bio": beta_mut_bio,
    "ga_pooled_2mers": sorted(GA_POOLED),
    "ga_spline_class_rule": {"AA": "AA", "GA": "GA", "_default": "SHARED"},
    "spline_theta": {k: v.tolist() for k, v in splines.items()},
    # Evaluated per position so the browser needs no spline arithmetic at all.
    "spline_by_position": {"first_position": lo,
                           **{k: np.round(Bg @ v, 6).tolist() for k, v in splines.items()}},
    # Correct fallback where a submitted sequence cannot be aligned: the genome-mean of that
    # class's spline. NOT zero -- sum-to-zero constrains the coefficients, not the curve.
    "unaligned_fallback": {k: float((Bg @ v).mean()) for k, v in splines.items()},
    "recovery_report": report,
}
with open(ARGS.out, "w") as fh:
    json.dump(bundle, fh)
print(f"wrote {ARGS.out}  version={ARGS.version}  release_date={ARGS.release_date}")
print(json.dumps(report, indent=2))
