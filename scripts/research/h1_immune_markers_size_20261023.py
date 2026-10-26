#!/usr/bin/env python3
"""Exploratory marker-by-size extension of H1; current accepted source is read only."""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import scipy.linalg as la
import statsmodels.api as sm

from h1_immune_types_size_20261023 import input_data, fit_one, build_design, write_json

MARKERS = ["CD45", "CD3e", "CD8a", "CD4", "CD20", "CD68", "CD163", "HLADR", "CD35", "LGALS3", "CD44", "VIM"]
BINS = ["<50 um", "50-<100 um", ">=100 um"]
PROTOCOL = {
    "date": "2026-10-23",
    "design": "Exploratory extension, frozen before marker-specific outcome fits; prior broad H1 results and size inventories were already known",
    "markers": MARKERS,
    "scope": "11 immune-associated common markers plus shared-state VIM; all 31 common markers profiled descriptively within curated immune cells; nuclear source Ki67 audited separately",
    "population": "Marker-positive cells WITHIN the six mutually exclusive curated immune phenotypes; marker positivity across all cells is descriptive, not immune identity",
    "size": "Equivalent circle core diameter sqrt(4*core_area_um2/pi), <50,50-<100,>=100 um; immune-independent geometry",
    "threshold": "Current generic scaled fluorescence >=0.10; sensitivity >=0.05 and >=0.20; arithmetic means, not log expression; these thresholds do not redefine phenotype",
    "exposure": "Marker-positive immune core / (marker-positive immune core + peri), per 0.10; undefined if no marker-positive nearby immune cells",
    "outcome": "Beta cell core count; log core area in mm2 offset; no total-cell-fraction outcome",
    "model": "Poisson quasi mean, donor fixed intercepts, equal total object weight per donor; class, log core area, log peri/core area, log1p marker-positive nearby immune, log1p marker-negative nearby immune, log1p independent broad immune density20-100um",
    "secondary": "Same marker models additionally condition on overall immune allocation; unidentified fits stay unestimated; this is a conditional association, not a marker-specific causal effect",
    "inference": "Donor-cluster CR1 sandwich; t reference informative donors-1; informative donor has within-donor exposure variation",
    "gates": {"informative_donors": 6, "core_positive_objects": 30, "eligible_objects": 100},
    "zero_beta_donors": "Exclude/list donors with no Beta outcome count in isolated model, retain individual zero-Beta objects",
    "coverage": "Require every curated immune cell in both compartments measured for the marker for an object to enter that marker model; never fill unknown marker counts with zero",
    "multiplicity": {
        "primary": "12 markers x3 geometric sizes =36 all-cohort threshold0.10 slots",
        "subgroups": "12x3x3 diagnostic groups =108 exploratory threshold0.10 slots",
        "threshold_sensitivity": "12x3x2 thresholds =72 exploratory all-cohort slots",
        "conditional_broad": "12x3 =36 exploratory all-cohort threshold0.10 slots",
        "missing": "Untestable slots internally use p=1 for BH but displayed p/q stay missing"
    },
    "interpretation": "Markers overlap and typing markers are partly defined by identity gates; parallel marker tests are correlated and do not identify independent mechanisms. HLADR/CD35/LGALS3/CD44/VIM are shared states, not lineage/function labels. No traffic, killing, activation, exhaustion, Treg or M1/M2 inference."
}


def marker_frame(base, table, marker, threshold):
    part = table[(table.marker == marker) & np.isclose(table.threshold, threshold)].copy()
    required = ["n_all", "n_all_measured", "n_all_positive", "n_immune", "n_immune_measured", "n_immune_positive"]
    if part.duplicated(["islet_id", "compartment"]).any():
        raise ValueError("Duplicate marker compartment rows")
    wide = part.pivot(index="islet_id", columns="compartment", values=required)
    wide.columns = [f"{v}_{c}" for v, c in wide.columns]
    d = base.merge(wide, left_on="islet_id", right_index=True, how="left", validate="one_to_one")
    coverage_cols = [f"{v}_{c}" for v in ["n_all", "n_all_measured", "n_immune_measured"] for c in ["core", "peri"]]
    if d[coverage_cols].isna().any().any():
        raise ValueError(f"Missing compartment row or measured denominator for {marker}/{threshold}")
    for compartment in ["core", "peri"]:
        if not d[f"n_immune_{compartment}_x" if f"n_immune_{compartment}_x" in d else f"n_immune_{compartment}"].eq(d[f"n_immune_{compartment}_y" if f"n_immune_{compartment}_y" in d else f"n_immune_{compartment}"]).all():
            raise ValueError("Immune population mismatch")
    # Rename possible merge collisions explicitly back to the audited broad counts.
    for comp in ["core", "peri"]:
        if f"n_immune_{comp}_x" in d:
            d = d.rename(columns={f"n_immune_{comp}_x": f"n_immune_{comp}"}).drop(columns=f"n_immune_{comp}_y")
    d["marker_complete"] = (d.n_immune_measured_core == d.n_immune_core) & (d.n_immune_measured_peri == d.n_immune_peri)
    d["type_core"] = d.n_immune_positive_core
    d["type_peri"] = d.n_immune_positive_peri
    d["type_nearby"] = d.type_core + d.type_peri
    d["other_nearby"] = d.n_immune_nearby - d.type_nearby
    if (d[["type_core", "type_peri", "other_nearby"]] < 0).any().any():
        raise ValueError("Invalid marker counts")
    d["type_allocation"] = d.type_core.div(d.type_nearby.where(d.type_nearby > 0))
    return d


def fit_conditional(d, infer=True):
    """Model with broad immune location held fixed, to expose redundancy."""
    x = build_design(d)
    x["broad_allocation_10pp"] = d.n_immune_core / d.n_immune_nearby / .1
    a = x.to_numpy(float)
    rank = np.linalg.matrix_rank(a)
    if rank != a.shape[1]:
        raise ValueError(f"Marker and broad allocation not separately identifiable: rank {rank}/{a.shape[1]}")
    condition = float(np.linalg.cond(a))
    if condition > 1e8:
        raise ValueError(f"Ill-conditioned design: {condition:g}")
    sizes = d.case_id.map(d.groupby("case_id").size()).to_numpy(float)
    w = 1 / sizes
    w /= w.mean()
    y = d.n_beta_core.to_numpy(float)
    fit = sm.GLM(y, a, family=sm.families.Poisson(), offset=np.log(d.core_area_um2.to_numpy(float) / 1e6), freq_weights=w).fit(maxiter=100, tol=1e-9)
    if not fit.converged or not np.isfinite(fit.params).all():
        raise ValueError("Conditional fit did not converge")
    j = x.columns.get_loc("allocation_10pp")
    b = float(fit.params[j])
    detail = {"design_columns": list(x), "condition_number": condition, "coefficients": fit.params.tolist(), "donors": sorted(d.case_id.unique()), "converged": True}
    if not infer:
        return {"log_rr": b, "rr": float(np.exp(b))}, detail
    group = d.case_id.to_numpy()
    donors = sorted(d.case_id.unique())
    g, n, p = len(donors), len(d), a.shape[1]
    mu = np.asarray(fit.fittedvalues)
    scores = np.stack([(a * (w * (y - mu))[:, None])[group == donor].sum(axis=0) for donor in donors])
    bread = la.pinvh((a.T * (w * mu)) @ a)
    cov = g / (g - 1) * (n - 1) / (n - p) * bread @ (scores.T @ scores) @ bread
    se = float(np.sqrt(max(cov[j, j], 0)))
    df_t = int(d.groupby("case_id").type_allocation.nunique().gt(1).sum()) - 1
    if se <= 0 or not np.isfinite(se):
        raise ValueError("Conditional uncertainty not identifiable")
    critical = stats.t.ppf(.975, df_t)
    return {"log_rr": b, "rr": float(np.exp(b)), "se_cluster": se, "lower": float(np.exp(b-critical*se)), "upper": float(np.exp(b+critical*se)), "p_cluster_t": float(2*stats.t.sf(abs(b/se), df_t)), "t_df": df_t}, detail


def run(base, table, out):
    rows, details, lodo = [], {}, []
    for threshold in [.10, .05, .20]:
        cohorts = ["All", "ND", "Aab+", "T1D"] if threshold == .10 else ["All"]
        for marker in MARKERS:
            d_all = marker_frame(base, table, marker, threshold)
            for size_bin in BINS:
                for cohort in cohorts:
                    raw = d_all[d_all.equivalent_diameter == size_bin]
                    if cohort != "All":
                        raw = raw[raw.donor_status == cohort]
                    before = raw[raw.base_eligible & raw.marker_complete & (raw.type_nearby > 0)].copy()
                    z = before.groupby("case_id").n_beta_core.sum()
                    zero = sorted(z[z == 0].index)
                    d = before[~before.case_id.isin(zero)].copy()
                    informative = int(d.groupby("case_id").type_allocation.nunique().gt(1).sum())
                    positive = int((d.type_core > 0).sum())
                    fail = []
                    if informative < 6: fail.append("fewer_than_6_informative_donors")
                    if positive < 30: fail.append("fewer_than_30_core_positive_objects")
                    if len(d) < 100: fail.append("fewer_than_100_objects")
                    models = ["parallel", "conditional_broad"] if threshold == .10 and cohort == "All" else ["parallel"]
                    for model in models:
                        family = "conditional_broad" if model == "conditional_broad" else "threshold_sensitivity" if threshold != .10 else "primary" if cohort == "All" else "subgroups"
                        key = f"{model}|{threshold}|{marker}|{size_bin}|{cohort}"
                        row = {"marker": marker, "threshold": threshold, "size_bin": size_bin, "cohort": cohort, "model": model, "family": family,
                            "source_objects": len(raw), "base_eligible_objects": int(raw.base_eligible.sum()), "incomplete_marker_objects": int((raw.base_eligible & ~raw.marker_complete).sum()),
                            "exposed_before_zero_beta_removal": len(before), "objects": len(d), "donors": d.case_id.nunique(), "informative_donors": informative,
                            "core_positive_objects": positive, "marker_positive_core_cells": int(d.type_core.sum()), "marker_positive_peri_memberships": int(d.type_peri.sum()),
                            "all_zero_beta_donors": ";".join(zero), "status": "DESCRIPTIVE_ONLY", "reason": ";".join(fail),
                            "log_rr": np.nan, "rr": np.nan, "lower": np.nan, "upper": np.nan, "p_cluster_t": np.nan, "q_bh": np.nan}
                        if not fail:
                            try:
                                result, detail = fit_conditional(d) if model == "conditional_broad" else fit_one(d)
                                row.update(result)
                                row.update(status="FITTED_EXPLORATORY", reason="")
                                details[key] = detail
                                if family == "primary":
                                    for donor in sorted(d.case_id.unique()):
                                        lr = {"marker": marker, "size_bin": size_bin, "omitted_donor": donor, "rr": np.nan, "status": "FAILED"}
                                        try:
                                            point, _ = fit_one(d[d.case_id != donor], infer=False)
                                            lr.update(rr=point["rr"], status="FITTED_POINT_ONLY")
                                        except Exception as exc:
                                            lr["reason"] = str(exc)
                                        lodo.append(lr)
                            except Exception as exc:
                                row["reason"] = "model_failure: " + str(exc)
                                details[key] = {"error": str(exc)}
                        rows.append(row)
    result = pd.DataFrame(rows)
    for family, indices in result.groupby("family").groups.items():
        p = result.loc[indices, "p_cluster_t"].to_numpy(float)
        working = np.where(np.isfinite(p), p, 1.)
        order = np.argsort(working, kind="stable")
        qsort = np.minimum.accumulate((working[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
        q = np.empty_like(qsort)
        q[order] = np.minimum(qsort, 1.)
        q[~np.isfinite(p)] = np.nan
        result.loc[indices, "q_bh"] = q
        result.loc[indices, "bh_family_slots"] = len(p)
    result.to_csv(out/"estimates.csv", index=False)
    pd.DataFrame(lodo).to_csv(out/"leave_one_donor_out.csv", index=False)
    write_json(out/"fit_details.json", details)
    summary = {"slots": len(result), "statuses": result.status.value_counts().to_dict(), "families": result.groupby("family").agg(slots=("status", "size"), fitted=("p_cluster_t", lambda s: int(s.notna().sum())), minimum_q=("q_bh", "min")).reset_index().to_dict("records")}
    write_json(out/"summary.json", summary)
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, default=Path("output/h1_20261022/h1_object_features.csv"))
    p.add_argument("--marker-counts", type=Path, default=Path("output/h1_immune_size_20261023/object_marker_counts.csv"))
    p.add_argument("--output", type=Path, default=Path("output/h1_immune_size_20261023/marker_models"))
    p.add_argument("--freeze-only", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    protocol = dict(PROTOCOL, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), features_sha256=hashlib.sha256(a.features.read_bytes()).hexdigest())
    if a.marker_counts.exists():
        protocol["marker_counts_sha256"] = hashlib.sha256(a.marker_counts.read_bytes()).hexdigest()
    write_json(a.output/"protocol.json", protocol)
    if a.freeze_only:
        print("Marker model protocol frozen; no fits run")
        return
    base = input_data(a.features)
    table = pd.read_csv(a.marker_counts, dtype={"islet_id": str, "case_id": str})
    if table.compartment.eq("islet").any():
        table["compartment"] = table.compartment.replace({"islet": "core"})
    run(base, table, a.output)


if __name__ == "__main__":
    main()
