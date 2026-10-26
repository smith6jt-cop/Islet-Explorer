#!/usr/bin/env python3
"""H1 immune-phenotype and size stratification, read-only source data.

Uses the audited 2026-10-22 feature table. No old analysis, app, or source file is
modified. Sparse or unsupported slots remain descriptive with no fabricated p.
"""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg as la
from scipy import stats
import statsmodels.api as sm

STATUSES = ["ND", "Aab+", "T1D"]
TYPES = {
    "Broad immune union": None,
    "CD8a T cell": "CD8a_Tcell",
    "CD4 T cell": "CD4_Tcell",
    "Other T cell": "T_cell",
    "B cell": "B_cell",
    "Macrophage": "Macrophage",
    "Unresolved immune": "Immune",
}
SIZE_BINS = {
    "equivalent_diameter": ["<50 um", "50-<100 um", ">=100 um"],
    "endocrine_cell_count": ["<10 endocrine cells", "10-49 endocrine cells", ">=50 endocrine cells"],
}
PROTOCOL = {
    "date": "2026-10-23",
    "design": "Exploratory stratification of H1; fixed specifications before fitting",
    "observations": "Current paired Islet and EndocrineCluster objects; no Singles or donor 78",
    "size_primary": "Equivalent diameter sqrt(4 * core_area_um2 / pi): <50, 50-<100, >=100 um",
    "size_sensitivity": "Typed endocrine cells in core: <10, 10-49, >=50; not total cells",
    "size_sensitivity_caveat": "Endocrine count includes Beta outcome and hormone-dependent typing; conditional selection is not immune-independent proof",
    "types": list(TYPES),
    "type_identity": "Six mutually exclusive current phenotypes; sum exactly equals broad immune union",
    "exposure": "Type core/(type core + type peri), in 0.10 units; undefined when type nearby pool is zero",
    "outcome": "n_beta_core, with log(core_area_um2/1e6) offset",
    "eligibility": "Positive finite core/peri areas, field20_100_area_um2 >=1000, finite nonnegative independent broad field density, positive type pool",
    "cohorts": ["All", "ND", "Aab+", "T1D"],
    "nuisance": ["donor fixed intercepts", "Islet/EndocrineCluster class", "log core area", "log peri/core area", "log1p type nearby count", "log1p OTHER immune nearby count", "log1p independent broad 20-100 um field density per mm2"],
    "weighting": "Equal total object weight for each included donor within each isolated fit",
    "zero_beta_donors": "Exclude and list donors whose Beta count sum is zero in the type/size/cohort analysis population",
    "inference_gate": {"informative_donors": 6, "core_positive_objects": 30, "eligible_objects_after_zero_beta_donor_removal": 100},
    "informative_donor": "At least two distinct type-allocation values among eligible nonzero-Beta-outcome-donor objects",
    "mean_model": "Poisson quasi maximum likelihood, maxiter=100, tol=1e-9; not a Poisson variance assumption",
    "uncertainty": "Donor-cluster sandwich, finite G and N/P correction; t reference with informative_donors-1 degrees of freedom",
    "failure": "DESCRIPTIVE_ONLY for gate failures, no allocation variation, unidentified/nonconverged fit; no p or CI",
    "multiplicity": {
        "primary_subtypes": "All-cohort geometric-size 6 types x 3 bins =18 slots",
        "primary_broad": "All-cohort geometric-size broad union x3 bins =3 slots, separate family",
        "subgroup_subtypes": "Geometric-size subgroups 3 groups x6 types x3 bins =54 exploratory slots",
        "subgroup_broad": "Geometric-size subgroups 3 groups x3 bins =9 exploratory slots",
        "sensitivity_subtypes": "Endocrine-count all-cohort 18 exploratory slots, separate family",
        "sensitivity_broad": "Endocrine-count all-cohort broad 3 exploratory slots, separate family",
        "sensitivity_subgroup_subtypes": "Endocrine-count subgroup 54 exploratory slots, separate family",
        "sensitivity_subgroup_broad": "Endocrine-count subgroup broad 9 exploratory slots, separate family",
        "missing_slots": "Untestable slots internally contribute p=1 to BH; their displayed p and q stay missing",
    },
    "descriptive_denominators": "Full donor/type/size tables use all base-eligible objects, including type-unexposed objects; model eligibility separately requires type pool >0",
    "density_units": "Cells/mm2 for type densities; Beta descriptive means also supplied as cells/1000 um2",
    "membership": "Peri and overlapping field counts are object-associated memberships, not unique tissue-cell totals",
    "causality": "Current occupancy and Beta density are cross-sectional; effects are not trafficking or killing estimates",
    "leave_one_donor_out": "Point-estimate stability only for fitted all-cohort geometric-size slots; no additional significance tests",
}


def write_json(path: Path, value) -> None:
    def convert(x):
        if isinstance(x, (np.integer, np.floating, np.bool_)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value, indent=2, default=convert) + "\n")


def input_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"case_id": str, "islet_id": str})
    columns = {"case_id", "islet_id", "donor_status", "object_class", "core_area_um2", "peri_area_um2", "n_beta_core", "n_endocrine_core", "n_immune_core", "n_immune_peri", "n_immune_nearby", "field20_100_area_um2", "immune_density_field20_100_mm2"}
    for suffix in list(TYPES.values())[1:]:
        columns.update({f"n_core_{suffix}", f"n_peri_{suffix}"})
    if columns - set(df):
        raise ValueError(f"Missing columns: {sorted(columns-set(df))}")
    if df.duplicated(["case_id", "islet_id"]).any():
        raise ValueError("Duplicate object identities")
    if pd.to_numeric(df.case_id, errors="coerce").eq(78).any():
        raise ValueError("Excluded donor 78 in input")
    if set(df.object_class) != {"Islet", "EndocrineCluster"}:
        raise ValueError(f"Unexpected source classes: {df.object_class.unique()}")
    if df.case_id.nunique() != 20:
        raise ValueError("This frozen analysis requires the audited 20-donor input")
    for region in ["core", "peri"]:
        total = df[[f"n_{region}_{suffix}" for suffix in list(TYPES.values())[1:]]].sum(axis=1)
        if not total.eq(df[f"n_immune_{region}"]).all():
            raise ValueError(f"Exclusive {region} phenotypes do not sum to broad immune union")
    if not (df.n_immune_core + df.n_immune_peri).eq(df.n_immune_nearby).all():
        raise ValueError("Immune pool identity failed")
    df["equivalent_diameter_um"] = np.sqrt(4*df.core_area_um2/np.pi)
    df["equivalent_diameter"] = pd.cut(df.equivalent_diameter_um, [0, 50, 100, np.inf], right=False, labels=SIZE_BINS["equivalent_diameter"]).astype(str)
    df["endocrine_cell_count"] = pd.cut(df.n_endocrine_core, [-np.inf, 10, 50, np.inf], right=False, labels=SIZE_BINS["endocrine_cell_count"]).astype(str)
    finite = np.isfinite(df[["core_area_um2", "peri_area_um2", "n_beta_core", "field20_100_area_um2", "immune_density_field20_100_mm2"]]).all(axis=1)
    df["base_eligible"] = finite & (df.core_area_um2 > 0) & (df.peri_area_um2 > 0) & (df.field20_100_area_um2 >= 1000) & (df.immune_density_field20_100_mm2 >= 0)
    df["beta_density_mm2"] = 1e6*df.n_beta_core/df.core_area_um2
    df["beta_density_1000um2"] = 1000*df.n_beta_core/df.core_area_um2
    return df


def type_data(df: pd.DataFrame, suffix: str | None) -> pd.DataFrame:
    d = df.copy()
    d["type_core"] = d.n_immune_core if suffix is None else d[f"n_core_{suffix}"]
    d["type_peri"] = d.n_immune_peri if suffix is None else d[f"n_peri_{suffix}"]
    d["type_nearby"] = d.type_core + d.type_peri
    d["other_nearby"] = d.n_immune_nearby - d.type_nearby
    if (d[["type_core", "type_peri", "other_nearby"]] < 0).any().any():
        raise ValueError("Negative phenotype counts or OTHER immune pool")
    d["type_allocation"] = d.type_core.div(d.type_nearby.where(d.type_nearby > 0))
    return d


def descriptive_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    donor_rows, mean_rows = [], []
    donor_info = df[["case_id", "donor_status"]].drop_duplicates().sort_values("case_id")
    for size_kind, bins in SIZE_BINS.items():
        for immune_type, suffix in TYPES.items():
            typed = type_data(df, suffix)
            for size_bin in bins:
                sized = typed[typed[size_kind] == size_bin]
                for info in donor_info.itertuples(index=False):
                    raw = sized[sized.case_id == info.case_id]
                    base = raw[raw.base_eligible]
                    exposed = base[base.type_nearby > 0]
                    core_pos = base[base.type_core > 0]
                    n = len(base)
                    donor_rows.append({"case_id": info.case_id, "donor_status": info.donor_status, "size_definition": size_kind, "size_bin": size_bin, "immune_type": immune_type,
                        "source_objects": len(raw), "base_eligible_objects": n, "type_exposed_objects": len(exposed), "core_positive_objects": len(core_pos),
                        "islets": int((base.object_class == "Islet").sum()), "clusters": int((base.object_class == "EndocrineCluster").sum()),
                        "type_core_cells": int(base.type_core.sum()), "type_peri_memberships": int(base.type_peri.sum()), "other_nearby_memberships": int(base.other_nearby.sum()),
                        "core_area_um2": float(base.core_area_um2.sum()), "peri_area_um2": float(base.peri_area_um2.sum()),
                        "type_core_density_mm2": float(1e6*base.type_core.sum()/base.core_area_um2.sum()) if n else np.nan,
                        "type_peri_density_mm2": float(1e6*base.type_peri.sum()/base.peri_area_um2.sum()) if n else np.nan,
                        "type_core_cells_per_positive_object": float(core_pos.type_core.mean()) if len(core_pos) else np.nan,
                        "fraction_all_objects_core_positive": len(core_pos)/n if n else np.nan,
                        "fraction_exposed_objects_core_positive": len(core_pos)/len(exposed) if len(exposed) else np.nan,
                        "allocation_unique_exposed": int(exposed.type_allocation.nunique()), "beta_cells_base": int(base.n_beta_core.sum()),
                        "mean_beta_density_base_cells_per_1000um2": float(base.beta_density_1000um2.mean()) if n else np.nan,
                        "beta_cells_type_exposed": int(exposed.n_beta_core.sum()), "beta_zero_objects_type_exposed": int((exposed.n_beta_core == 0).sum())})
                    for occupancy, part in [("All type-exposed", exposed), ("Core immune absent", exposed[exposed.type_core == 0]), ("Core immune present", exposed[exposed.type_core > 0])]:
                        mean_rows.append({"case_id": info.case_id, "donor_status": info.donor_status, "size_definition": size_kind, "size_bin": size_bin, "immune_type": immune_type, "occupancy": occupancy,
                            "objects": len(part), "mean_beta_density_cells_per_1000um2": float(part.beta_density_1000um2.mean()) if len(part) else np.nan,
                            "median_beta_density_cells_per_1000um2": float(part.beta_density_1000um2.median()) if len(part) else np.nan,
                            "mean_type_core_density_mm2": float((1e6*part.type_core/part.core_area_um2).mean()) if len(part) else np.nan,
                            "mean_type_peri_density_mm2": float((1e6*part.type_peri/part.peri_area_um2).mean()) if len(part) else np.nan})
    return pd.DataFrame(donor_rows), pd.DataFrame(mean_rows)


def build_design(d: pd.DataFrame) -> pd.DataFrame:
    x = pd.DataFrame({"intercept": np.ones(len(d))}, index=d.index)
    for donor in sorted(d.case_id.unique())[1:]:
        x["donor_" + donor] = (d.case_id == donor).astype(float)
    x["class_cluster"] = (d.object_class == "EndocrineCluster").astype(float)
    for name, values in {
        "log_core_area": np.log(d.core_area_um2),
        "log_peri_core_area_ratio": np.log(d.peri_area_um2/d.core_area_um2),
        "log1p_type_nearby": np.log1p(d.type_nearby),
        "log1p_other_immune_nearby": np.log1p(d.other_nearby),
        "log1p_broad_independent_field": np.log1p(d.immune_density_field20_100_mm2),
    }.items():
        sd = values.std(ddof=0)
        if sd > 1e-12:
            x[name] = (values-values.mean())/sd
    x["allocation_10pp"] = d.type_allocation/0.1
    x = x.loc[:, [col == "intercept" or x[col].std(ddof=0) > 1e-12 for col in x]]
    if "allocation_10pp" not in x:
        raise ValueError("Allocation has no variation")
    rank = np.linalg.matrix_rank(x.to_numpy(float))
    if rank != x.shape[1]:
        raise ValueError(f"Unidentified design rank {rank}/{x.shape[1]}")
    return x.astype(float)


def fit_one(d: pd.DataFrame, infer=True) -> tuple[dict, dict]:
    x = build_design(d)
    a = x.to_numpy(float)
    y = d.n_beta_core.to_numpy(float)
    sizes = d.case_id.map(d.groupby("case_id").size()).to_numpy(float)
    w = 1/sizes
    w /= w.mean()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = sm.GLM(y, a, family=sm.families.Poisson(), offset=np.log(d.core_area_um2.to_numpy(float)/1e6), freq_weights=w).fit(maxiter=100, tol=1e-9)
    if not fit.converged or not np.isfinite(fit.params).all():
        raise ValueError("Poisson quasi mean fit did not converge")
    j = x.columns.get_loc("allocation_10pp")
    coef = float(fit.params[j])
    detail = {"design_columns": x.columns.tolist(), "coefficients": fit.params.tolist(), "converged": bool(fit.converged), "warnings": sorted({str(c.message) for c in caught}), "donors": sorted(d.case_id.unique())}
    if not infer:
        return {"log_rr": coef, "rr": float(np.exp(coef))}, detail
    group = d.case_id.to_numpy()
    donors = sorted(d.case_id.unique())
    g, n, p = len(donors), len(d), a.shape[1]
    if n <= p or g < 2:
        raise ValueError("Insufficient observations for finite-cluster covariance")
    mu = np.asarray(fit.fittedvalues)
    cluster_scores = np.stack([(a*(w*(y-mu))[:, None])[group == donor].sum(axis=0) for donor in donors])
    bread = la.pinvh((a.T*(w*mu))@a)
    covariance = (g/(g-1)*(n-1)/(n-p))*bread@(cluster_scores.T@cluster_scores)@bread
    se = float(np.sqrt(max(covariance[j, j], 0)))
    informative = int(d.groupby("case_id").type_allocation.nunique().gt(1).sum())
    df_t = informative-1
    if not np.isfinite(se) or se <= 0 or df_t < 1:
        raise ValueError("Cluster uncertainty is unavailable")
    critical = stats.t.ppf(0.975, df_t)
    result = {"log_rr": coef, "rr": float(np.exp(coef)), "se_cluster": se,
              "lower": float(np.exp(coef-critical*se)), "upper": float(np.exp(coef+critical*se)),
              "p_cluster_t": float(2*stats.t.sf(abs(coef/se), df_t)), "t_df": df_t,
              "pearson_dispersion": float(np.sum((y-mu)**2/np.maximum(mu, 1e-12))/(n-p))}
    detail.update({"covariance_cluster": covariance.tolist(), "t_reference_degrees_of_freedom": df_t,
                   "donor_total_weights": {donor: float(w[group == donor].sum()) for donor in donors}})
    return result, detail


def family_name(size_kind, cohort, immune_type):
    prefix = "primary" if size_kind == "equivalent_diameter" else "sensitivity"
    subgroup = "_subgroup" if cohort != "All" else ""
    kind = "broad" if immune_type == "Broad immune union" else "subtypes"
    return prefix+subgroup+"_"+kind


def estimates(df: pd.DataFrame, do_fit: bool) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    rows, details, leaveout = [], {}, []
    for size_kind, bins in SIZE_BINS.items():
        for size_bin in bins:
            for immune_type, suffix in TYPES.items():
                typed = type_data(df[df[size_kind] == size_bin], suffix)
                for cohort in ["All"]+STATUSES:
                    allsub = typed if cohort == "All" else typed[typed.donor_status == cohort]
                    before = allsub[allsub.base_eligible & (allsub.type_nearby > 0)].copy()
                    zero = before.groupby("case_id").n_beta_core.sum()
                    zero_donors = sorted(zero[zero == 0].index.tolist())
                    d = before[~before.case_id.isin(zero_donors)].copy()
                    informative = int(d.groupby("case_id").type_allocation.nunique().gt(1).sum())
                    positive = int((d.type_core > 0).sum())
                    failures = []
                    if informative < 6:
                        failures.append("fewer_than_6_informative_donors")
                    if positive < 30:
                        failures.append("fewer_than_30_core_positive_objects")
                    if len(d) < 100:
                        failures.append("fewer_than_100_eligible_objects")
                    key = f"{size_kind}|{size_bin}|{immune_type}|{cohort}"
                    row = {"size_definition": size_kind, "size_bin": size_bin, "immune_type": immune_type, "cohort": cohort,
                           "family": family_name(size_kind, cohort, immune_type), "source_objects": len(allsub),
                           "base_eligible_objects": int(allsub.base_eligible.sum()), "type_exposed_before_zero_beta_removal": len(before),
                           "objects": len(d), "donors_before_zero_beta_removal": before.case_id.nunique(), "donors": d.case_id.nunique(),
                           "informative_donors": informative, "core_positive_objects": positive,
                           "type_core_cells": int(d.type_core.sum()), "type_peri_memberships": int(d.type_peri.sum()),
                           "all_zero_beta_donors": ";".join(zero_donors), "status": "DESCRIPTIVE_ONLY", "reason": ";".join(failures),
                           "log_rr": np.nan, "rr": np.nan, "se_cluster": np.nan, "lower": np.nan, "upper": np.nan,
                           "p_cluster_t": np.nan, "q_bh": np.nan, "t_df": np.nan, "pearson_dispersion": np.nan}
                    if not failures and not do_fit:
                        row["status"] = "READY_TO_FIT"
                    if not failures and do_fit:
                        try:
                            result, detail = fit_one(d)
                            row.update(result)
                            row["status"] = "FITTED_EXPLORATORY"
                            row["reason"] = ""
                            details[key] = detail
                            if cohort == "All" and size_kind == "equivalent_diameter":
                                for donor in sorted(d.case_id.unique()):
                                    sub = d[d.case_id != donor]
                                    lr = {"size_definition": size_kind, "size_bin": size_bin, "immune_type": immune_type, "omitted_case_id": donor, "remaining_donors": sub.case_id.nunique(), "log_rr": np.nan, "rr": np.nan, "status": "FAILED", "reason": ""}
                                    try:
                                        point, _ = fit_one(sub, infer=False)
                                        lr.update(point)
                                        lr["status"] = "FITTED_POINT_ONLY"
                                    except Exception as exc:
                                        lr["reason"] = str(exc)
                                    leaveout.append(lr)
                        except Exception as exc:
                            row["reason"] = "model_failure: "+str(exc)
                            details[key] = {"error": str(exc)}
                    rows.append(row)
    result = pd.DataFrame(rows)
    for family, indices in result.groupby("family").groups.items():
        p = result.loc[indices, "p_cluster_t"].to_numpy(float)
        working = np.where(np.isfinite(p), p, 1.0)
        order = np.argsort(working, kind="stable")
        ranked = working[order]
        qsort = np.minimum.accumulate((ranked*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
        q = np.empty_like(qsort)
        q[order] = np.minimum(qsort, 1.0)
        q[~np.isfinite(p)] = np.nan
        result.loc[indices, "q_bh"] = q
        result.loc[indices, "bh_family_slots"] = len(p)
    return result, details, pd.DataFrame(leaveout)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=Path("output/h1_20261022/h1_object_features.csv"))
    parser.add_argument("--output", type=Path, default=Path("output/h1_immune_size_20261023/types"))
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    df = input_data(args.features)
    args.output.mkdir(parents=True, exist_ok=True)
    protocol = dict(PROTOCOL)
    protocol.update({"input": str(args.features.resolve()), "input_sha256": hashlib.sha256(args.features.read_bytes()).hexdigest(),
                     "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     "software": {"numpy": np.__version__, "pandas": pd.__version__, "statsmodels": sm.__version__}})
    # Write frozen protocol before any outcome-model fit.
    write_json(args.output/"protocol.json", protocol)
    donor_counts, donor_means = descriptive_tables(df)
    donor_counts.to_csv(args.output/"donor_type_size_counts.csv", index=False)
    donor_means.to_csv(args.output/"donor_occupancy_means.csv", index=False)
    group_rows, balanced_rows = [], []
    for cohort in ["All"]+STATUSES:
        count_sub = donor_counts if cohort == "All" else donor_counts[donor_counts.donor_status == cohort]
        for (kind, size, immune_type), part in count_sub.groupby(["size_definition", "size_bin", "immune_type"], sort=False):
            row = {"cohort": cohort, "size_definition": kind, "size_bin": size, "immune_type": immune_type, "donors_in_cohort": len(part), "donors_with_exposed_objects": int((part.type_exposed_objects > 0).sum()), "donors_with_core_events": int((part.core_positive_objects > 0).sum())}
            for col in ["source_objects", "base_eligible_objects", "type_exposed_objects", "core_positive_objects", "type_core_cells", "type_peri_memberships", "islets", "clusters"]:
                row[col] = int(part[col].sum())
            group_rows.append(row)
        mean_sub = donor_means if cohort == "All" else donor_means[donor_means.donor_status == cohort]
        for (kind, size, immune_type, occupancy), part in mean_sub.groupby(["size_definition", "size_bin", "immune_type", "occupancy"], sort=False):
            balanced_rows.append({"cohort": cohort, "size_definition": kind, "size_bin": size, "immune_type": immune_type, "occupancy": occupancy,
                                  "donors_with_objects": int((part.objects > 0).sum()), "objects": int(part.objects.sum()),
                                  "equal_donor_mean_beta_density_cells_per_1000um2": float(part.mean_beta_density_cells_per_1000um2.mean()),
                                  "equal_donor_mean_type_core_density_mm2": float(part.mean_type_core_density_mm2.mean()),
                                  "equal_donor_mean_type_peri_density_mm2": float(part.mean_type_peri_density_mm2.mean())})
    pd.DataFrame(group_rows).to_csv(args.output/"cohort_counts.csv", index=False)
    pd.DataFrame(balanced_rows).to_csv(args.output/"donor_balanced_means.csv", index=False)
    result, details, lodo = estimates(df, do_fit=not args.preflight_only)
    result.to_csv(args.output/("support_preflight.csv" if args.preflight_only else "estimates.csv"), index=False)
    if not args.preflight_only:
        write_json(args.output/"fit_details.json", details)
        lodo.to_csv(args.output/"leave_one_donor_out.csv", index=False)
    summary = {"input_objects": len(df), "donors": df.case_id.nunique(), "base_eligible_objects": int(df.base_eligible.sum()),
               "model_slots": len(result), "statuses": result.status.value_counts().to_dict(),
               "families": result.groupby("family").agg(slots=("status", "size"), fitted=("p_cluster_t", lambda s: int(s.notna().sum()))).reset_index().to_dict("records")}
    write_json(args.output/("preflight_summary.json" if args.preflight_only else "summary.json"), summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
