#!/usr/bin/env python3
"""Donor-balanced matched descriptive check for H1, frozen 2026-10-22.

Run with the existing Islet Explorer Python environment. This reads derived H1
features, writes only its own output directory, and never changes source/app data.
Matching uses geometry/exposure, never the Beta outcome. It estimates a binary
core-occupancy contrast, not the continuous-fraction estimand of the main model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


PROTOCOL = {
    "date": "2026-10-22",
    "purpose": "Related matched descriptive sensitivity, not independent validation",
    "exposure": "n_immune_core >= 1 versus n_immune_core == 0",
    "exact_matching": ["case_id", "object_class"],
    "eligibility": [
        "Islet or EndocrineCluster (app Cluster); no Singles",
        "positive finite core and peri areas; nonnegative finite counts",
        "positive unique nearby immune pool",
        "finite nonnegative 20-100 micrometer field immune density",
        "field20_100_area_um2 >= 1000",
    ],
    "calipers": {
        "absolute log2 core area difference": 0.5,
        "donor-z log1p nearby immune density per mm2 difference": 0.5,
        "donor-z log1p 20-100 micrometer field immune density per mm2 difference": 0.5,
    },
    "nearby_immune_density": "1e6 * n_immune_nearby / (core_area_um2 + peri_area_um2)",
    "matching_algorithm": "Maximum-cardinality minimum-squared-caliper-distance 1:1 without replacement",
    "primary_outcome": "Paired difference in Beta cells per 1000 square micrometers",
    "sensitivity_outcome": "Paired difference log1p(Beta density per mm2), units fixed",
    "donor_effect": "Mean of within-donor matched pair differences",
    "minimum_pairs_for_balanced_effect": 5,
    "sparse_donors": "1-4 pairs are shown descriptively but excluded from donor-balanced pooled effect",
    "balance_flag": "Absolute SMD > 0.25 is flagged; donors are not selected based on outcomes",
    "uncertainty": "Donor bootstrap only if at least 6 donors have >=5 pairs; no pair-based inference",
    "bootstrap_seed": 20261022,
    "bootstrap_replicates": 10000,
    "limitations": [
        "An observational comparison, not a causal treatment effect",
        "Matching and main model use the same data and related estimands",
        "The density matching estimand differs from the main model count adjustment",
        "Shared immune bands do not create independent pairs; donor is the uncertainty unit",
        "Calipers and thresholds frozen before inspecting matched outcome differences",
    ],
}

MATCH_COLUMNS = ["log2_area", "z_log1p_nearby_density_mm2", "z_log1p_field_density_mm2"]
BALANCE_COLUMNS = ["log2_area", "log1p_nearby_density_mm2", "log1p_field_density_mm2"]
OBJECT_CLASSES = ["EndocrineCluster", "Islet"]


def donor_z(series: pd.Series) -> pd.Series:
    sd = float(series.std(ddof=0))
    return (series - series.mean()) / sd if sd > 0 else series * 0.0


def smd(treated: pd.Series, control: pd.Series, scale: float) -> float:
    if not len(treated) or not len(control):
        return float("nan")
    diff = float(treated.mean() - control.mean())
    if not np.isfinite(scale) or scale == 0:
        return 0.0 if diff == 0 else float("nan")
    return diff / scale


def match_one(group: pd.DataFrame) -> tuple[list[tuple[int, int, float]], int]:
    t = group[group.core_immune_present].sort_values("islet_id")
    c = group[~group.core_immune_present].sort_values("islet_id")
    if len(t) == 0 or len(c) == 0:
        return [], 0
    delta = np.abs(t[MATCH_COLUMNS].to_numpy()[:, None, :] - c[MATCH_COLUMNS].to_numpy()[None, :, :])
    feasible = np.all(delta <= 0.5 + 1e-12, axis=2)
    distance = np.sum((delta / 0.5) ** 2, axis=2)
    # A dummy penalty larger than every feasible assignment's total cost gives
    # maximum feasible cardinality first, minimum distance second.
    dummy = 4.0 * (min(len(t), len(c)) + 1)
    costs = np.full((len(t), len(c) + len(t)), dummy)
    costs[:, :len(c)] = np.where(feasible, distance, dummy * 10)
    rows, cols = linear_sum_assignment(costs)
    result = [(int(t.index[i]), int(c.index[j]), float(distance[i, j]))
              for i, j in zip(rows, cols) if j < len(c) and feasible[i, j]]
    return result, int(feasible.sum())


def summarize(effects: pd.DataFrame, seed: int) -> dict:
    effective = effects[effects.n_pairs >= 5]
    value = effective.beta_density_difference_mean.to_numpy(dtype=float)
    out = {
        "donors_with_any_pairs": int((effects.n_pairs > 0).sum()),
        "donors_with_at_least_5_pairs": int(len(effective)),
        "total_pairs": int(effects.n_pairs.sum()),
        "pairs_in_effective_donors": int(effective.n_pairs.sum()),
        "donor_balanced_beta_density_difference_cells_per_1000um2": float(value.mean()) if len(value) else None,
        "median_donor_beta_density_difference_cells_per_1000um2": float(np.median(value)) if len(value) else None,
        "donors_negative": int((value < 0).sum()),
        "donors_positive": int((value > 0).sum()),
        "donor_balanced_log1p_density_mm2_difference": float(effective.log1p_beta_density_mm2_difference_mean.mean()) if len(value) else None,
        "bootstrap_95pct_interval": None,
        "uncertainty_status": "Insufficient effective donors for bootstrap" if len(value) < 6 else "Donor percentile bootstrap; observational uncertainty only",
    }
    if len(value) >= 6:
        rng = np.random.default_rng(seed)
        resampled = value[rng.integers(0, len(value), size=(10000, len(value)))].mean(axis=1)
        out["bootstrap_95pct_interval"] = [float(x) for x in np.quantile(resampled, [0.025, 0.975])]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=Path("output/h1_20261022/h1_object_features.csv"))
    parser.add_argument("--output", type=Path, default=Path("output/h1_20261022/independent"))
    args = parser.parse_args()
    if not args.features.is_file():
        raise FileNotFoundError(f"Awaiting feature input: {args.features}")
    df = pd.read_csv(args.features, dtype={"case_id": str, "islet_id": str})
    required = {"case_id", "islet_id", "donor_status", "object_class", "n_beta_core", "core_area_um2", "peri_area_um2", "n_immune_core", "n_immune_peri", "n_immune_nearby", "immune_density_field20_100_mm2", "field20_100_area_um2"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required features: {sorted(missing)}")
    if df.duplicated(["case_id", "islet_id"]).any():
        raise ValueError("Duplicate case_id/islet_id keys")
    unknown_classes = set(df.object_class.dropna().unique()) - set(OBJECT_CLASSES)
    if unknown_classes:
        raise ValueError(f"Unexpected object classes; resolve explicitly: {sorted(unknown_classes)}")
    if pd.to_numeric(df.case_id, errors="coerce").eq(78).any():
        raise ValueError("Excluded donor 78 is present")
    number_cols = sorted(required - {"case_id", "islet_id", "donor_status", "object_class"})
    for col in number_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    valid = np.isfinite(df[number_cols]).all(axis=1)
    valid &= df.object_class.isin(OBJECT_CLASSES)
    valid &= (df.core_area_um2 > 0) & (df.peri_area_um2 > 0)
    valid &= df[["n_beta_core", "n_immune_core", "n_immune_peri", "n_immune_nearby", "immune_density_field20_100_mm2"]].ge(0).all(axis=1)
    valid &= (df.n_immune_nearby > 0) & (df.field20_100_area_um2 >= 1000)
    if 'core_complete' in df:
        columns = ['core_complete', 'neighborhood_complete_20', 'neighborhood_complete_100']
        if df[columns].isna().any().any() or not df[columns].isin([True, False]).all().all():
            raise ValueError('Invalid acquisition eligibility metadata')
        valid &= df[columns].all(axis=1)
    df["eligible"] = valid
    df["core_immune_present"] = df.n_immune_core >= 1
    df["log2_area"] = np.log2(df.core_area_um2.where(df.core_area_um2 > 0))
    df["nearby_density_mm2"] = 1e6 * df.n_immune_nearby / (df.core_area_um2 + df.peri_area_um2)
    df["log1p_nearby_density_mm2"] = np.log1p(df.nearby_density_mm2)
    df["log1p_field_density_mm2"] = np.log1p(df.immune_density_field20_100_mm2)
    df["beta_density_1000um2"] = 1000 * df.n_beta_core / df.core_area_um2
    df["log1p_beta_density_mm2"] = np.log1p(1e6 * df.n_beta_core / df.core_area_um2)
    eligible = df[valid].copy()
    for col in ["log1p_nearby_density_mm2", "log1p_field_density_mm2"]:
        eligible["z_" + col] = eligible.groupby("case_id")[col].transform(donor_z)
    matches, preflight, balances, effects = [], [], [], []
    for (donor, status), source in df.groupby(["case_id", "donor_status"], sort=True):
        local = eligible[eligible.case_id == donor]
        donor_pairs = []
        for object_class in OBJECT_CLASSES:
            part = local[local.object_class == object_class]
            selected, feasible = match_one(part)
            donor_pairs.extend(selected)
            preflight_row = {"case_id": donor, "donor_status": status, "object_class": object_class,
                             "n_source_objects": int((source.object_class == object_class).sum()),
                             "n_eligible": len(part), "n_core_immune_present": int(part.core_immune_present.sum()),
                             "n_core_immune_absent": int((~part.core_immune_present).sum()),
                             "n_feasible_edges": feasible, "n_matched_pairs": len(selected)}
            for col in MATCH_COLUMNS:
                pt, pc = part[part.core_immune_present][col], part[~part.core_immune_present][col]
                preflight_row[col + "_present_min"] = pt.min()
                preflight_row[col + "_present_max"] = pt.max()
                preflight_row[col + "_absent_min"] = pc.min()
                preflight_row[col + "_absent_max"] = pc.max()
                preflight_row[col + "_shared_min"] = max(pt.min(), pc.min()) if len(pt) and len(pc) else np.nan
                preflight_row[col + "_shared_max"] = min(pt.max(), pc.max()) if len(pt) and len(pc) else np.nan
            preflight.append(preflight_row)
        treated_indices = [p[0] for p in donor_pairs]
        control_indices = [p[1] for p in donor_pairs]
        t_all = local[local.core_immune_present]
        c_all = local[~local.core_immune_present]
        for col in BALANCE_COLUMNS:
            scale = float(np.sqrt((t_all[col].var(ddof=1) + c_all[col].var(ddof=1)) / 2))
            balances.append({"case_id": donor, "donor_status": status, "covariate": col,
                             "before_smd": smd(t_all[col], c_all[col], scale),
                             "after_smd": smd(local.loc[treated_indices, col], local.loc[control_indices, col], scale)})
        differences, log_differences = [], []
        for k, (ti, ci, cost) in enumerate(donor_pairs, start=1):
            t, c = local.loc[ti], local.loc[ci]
            difference = float(t.beta_density_1000um2 - c.beta_density_1000um2)
            log_difference = float(t.log1p_beta_density_mm2 - c.log1p_beta_density_mm2)
            differences.append(difference)
            log_differences.append(log_difference)
            row = {"case_id": donor, "donor_status": status, "object_class": t.object_class,
                   "pair_id": f"{donor}_{k}", "treated_islet_id": t.islet_id, "control_islet_id": c.islet_id,
                   "match_distance_squared": cost, "beta_density_difference_cells_per_1000um2": difference,
                   "log1p_beta_density_mm2_difference": log_difference}
            for col in ["n_beta_core", "core_area_um2", "peri_area_um2", "n_immune_core", "n_immune_peri", "n_immune_nearby", "nearby_density_mm2", "immune_density_field20_100_mm2", "field20_100_area_um2", "beta_density_1000um2"] + MATCH_COLUMNS:
                row["treated_" + col] = t[col]
                row["control_" + col] = c[col]
            matches.append(row)
        count = len(differences)
        effects.append({"case_id": donor, "donor_status": status, "n_eligible": len(local),
                        "n_core_immune_present": len(t_all), "n_core_immune_absent": len(c_all),
                        "n_pairs": count, "effective_donor": count >= 5,
                        "support_status": "effective" if count >= 5 else "descriptive_only" if count else "no_matches",
                        "beta_density_difference_mean": float(np.mean(differences)) if count else np.nan,
                        "beta_density_difference_median": float(np.median(differences)) if count else np.nan,
                        "log1p_beta_density_mm2_difference_mean": float(np.mean(log_differences)) if count else np.nan,
                        "treated_beta_density_mean": float(local.loc[treated_indices, "beta_density_1000um2"].mean()) if count else np.nan,
                        "control_beta_density_mean": float(local.loc[control_indices, "beta_density_1000um2"].mean()) if count else np.nan})
    args.output.mkdir(parents=True, exist_ok=True)
    matched_table = pd.DataFrame(matches)
    matched_table.to_csv(args.output / "matched_pairs.csv", index=False)
    if len(matched_table):
        class_effects = matched_table.groupby(["case_id", "donor_status", "object_class"]).agg(
            n_pairs=("pair_id", "count"),
            beta_density_difference_mean=("beta_density_difference_cells_per_1000um2", "mean"),
            beta_density_difference_median=("beta_density_difference_cells_per_1000um2", "median"),
        ).reset_index()
        class_effects.to_csv(args.output / "donor_class_effects.csv", index=False)
    pd.DataFrame(preflight).to_csv(args.output / "overlap_preflight.csv", index=False)
    balance = pd.DataFrame(balances)
    balance.to_csv(args.output / "covariate_balance.csv", index=False)
    effect = pd.DataFrame(effects)
    effect["max_absolute_after_smd"] = effect.case_id.map(balance.assign(abs_smd=balance.after_smd.abs()).groupby("case_id").abs_smd.max())
    effect["balance_flag_gt025"] = effect.max_absolute_after_smd > 0.25
    effect.to_csv(args.output / "donor_effects.csv", index=False)
    pooled = summarize(effect, PROTOCOL["bootstrap_seed"])
    pooled["input_objects"] = len(df)
    pooled["eligible_objects"] = len(eligible)
    pooled["group_descriptive_summaries"] = {str(status): summarize(g, PROTOCOL["bootstrap_seed"]) for status, g in effect.groupby("donor_status")}
    effective = effect[effect.effective_donor]
    leaveout = [{"omitted_case_id": row.case_id,
                 "remaining_donors": len(effective) - 1,
                 "donor_balanced_beta_density_difference_cells_per_1000um2": float(effective.loc[effective.case_id != row.case_id, "beta_density_difference_mean"].mean())}
                for row in effective.itertuples()]
    pd.DataFrame(leaveout).to_csv(args.output / "leave_one_donor_out.csv", index=False)
    protocol = dict(PROTOCOL)
    protocol["feature_path"] = str(args.features.resolve())
    protocol["feature_sha256"] = hashlib.sha256(args.features.read_bytes()).hexdigest()
    protocol["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    (args.output / "summary.json").write_text(json.dumps(pooled, indent=2, allow_nan=False) + "\n")
    lines = ["# H1 matched descriptive sensitivity", "", "This is a related binary-occupancy comparison using the same observations as the main H1 analysis, not independent validation.", "",
             f"Eligible objects: {len(eligible)} / {len(df)}. Matched pairs: {pooled['total_pairs']}. Donors with any pairs: {pooled['donors_with_any_pairs']}. Donors with at least five pairs: {pooled['donors_with_at_least_5_pairs']}.", "",
             "Exact matches are within donor and object class. Frozen calipers are 0.5 in log2 core area and 0.5 donor SD in log1p nearby/field immune density per mm2. The primary matched endpoint is treated-minus-control Beta cell density per 1000 square micrometers, with equal weight per eligible donor. Donors with 1-4 pairs remain descriptive.", "",
             "No pair-level tests were performed. Shared bands, nonrandom immune occupancy, and unmeasured historical injury limit interpretation.", "",
             "See protocol.json for the frozen estimand and exclusions, overlap_preflight.csv for all donor/classes, covariate_balance.csv for matching quality, donor_effects.csv for the effective replication, and matched_pairs.csv for every selected pair.", "",
             "```json", json.dumps(pooled, indent=2), "```", ""]
    (args.output / "README.md").write_text("\n".join(lines))
    print(json.dumps(pooled, indent=2))


if __name__ == "__main__":
    main()
