#!/usr/bin/env python3
"""Check original-donor source records while allowing joint derived refits.

The gate aligns rows by identity, never by order. Original-donor islet additions,
removals and duplicate identities fail. New donors and new fields are permitted;
original identities, geometry, measurements, missingness, phenotype composition
and the unharmonized scVI latent must remain equal. Only explicitly named
trajectory, UMAP, Leiden, donor-batch Harmony and supporting fit fields may differ;
their differences are reported separately and do not fail source-record parity.
Absolute numeric tolerance defaults to zero (rtol=0). New keys are permitted.
No mapping from regenerated schema-4 annotations to frozen islets is invented.

``audit(baseline_root, candidate_root, atol=0.0)`` is read-only. The CLI alone
writes a JSON report, exclusively to an explicit directory outside app data and
both input roots. A missing canonical analogue is a failed/incomparable result.
The checked artifacts are islet-level H5ADs, embedded tables, Barlow features
and an optional workbook. This gate does not check the original single-cell
H5ADs; source-cell identity, measurements and labels require separate parity
evidence before any 21-donor integration.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

DONOR_COLUMNS = ("imageid", "case_id", "Case ID")
ARTIFACTS = ("adata_ins_root.h5ad", "islet_explorer.h5ad", "barlow_islet_features.csv")

# Explicit allowlists from the frozen 15-donor H5ADs. These are results of a
# cohort fit, not measurements or cell-type assignments. Unknown baseline fields
# continue to require exact parity, including X_scVI_mean and phenotype arrays.
DERIVED_OBS_COLUMNS = frozenset({
    "dpt_pseudotime", "dpt_pseudotime_combined", "dpt_pseudotime_small",
    "dpt_pseudotime_large", "DC1", "DC2", "DC1_combined", "DC2_combined",
    "leiden_0.3", "leiden_0.5", "leiden_0.8", "leiden_1.0",
    "leiden_umap_1", "leiden_umap_2",
})
DERIVED_OBSM_KEYS = frozenset({
    "X_diffmap", "X_pca", "X_umap", "X_scVI_harmony_pt",
    "X_scVI_combined_harmony_pt",
})
DERIVED_OBSP_KEYS = frozenset({
    "connectivities", "distances", "neighbors_viz_connectivities",
    "neighbors_viz_distances",
})
DERIVED_VARM_KEYS = frozenset({"PCs"})
DERIVED_UNS_KEYS = frozenset({
    "diffmap_evals", "donor_status_sizes", "iroot", "neighbors",
    "neighbors_viz", "paga", "pca", "split_pseudotime", "umap",
})


def donor_id(value: Any) -> str:
    """Normalize donor keys only for matching; raw column values are still checked."""
    if pd.isna(value):
        raise ValueError("Missing donor identity")
    text = str(value)
    try:
        number = float(text)
        if math.isfinite(number) and number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def donor_column(frame: pd.DataFrame) -> str:
    return next((column for column in DONOR_COLUMNS if column in frame), "")


def table_keys(frame: pd.DataFrame) -> list[str]:
    donor = donor_column(frame)
    if not donor:
        raise ValueError("No donor identity column; original subset is incomparable")
    if "region" in frame:
        return [donor, "region"] + [c for c in ("region_type", "type", "class", "marker") if c in frame]
    for identity in ("combined_islet_id", "islet_id", "islet_key", "name"):
        if identity in frame:
            return [donor, identity]
    if not frame[donor].map(donor_id).duplicated().any():
        return [donor]
    raise ValueError("No stable islet/region identity; original subset is incomparable")


def _check(report: dict, name: str, passed: bool, **details) -> None:
    report["checks"].append({"name": name, "passed": bool(passed), **details})


def _equal_arrays(left, right, atol: float) -> tuple[int, float | None]:
    """Value comparison with explicit missingness; no relative tolerance."""
    left, right = np.asarray(left), np.asarray(right)
    if left.shape != right.shape:
        return max(left.size, right.size, 1), None
    left, right = left.reshape(-1), right.reshape(-1)
    missing_left, missing_right = pd.isna(left), pd.isna(right)
    unequal = missing_left != missing_right
    valid = ~missing_left & ~missing_right
    maximum = None
    if np.issubdtype(left.dtype, np.number) and np.issubdtype(right.dtype, np.number):
        equal = (left[valid] == right[valid]) if atol == 0 else np.isclose(
            left[valid], right[valid], rtol=0, atol=atol, equal_nan=False)
        unequal[valid] = ~equal
        finite = valid & np.isfinite(left) & np.isfinite(right)
        if np.any(finite):
            maximum = float(np.max(np.abs(left[finite].astype(float) - right[finite].astype(float))))
    else:
        unequal[valid] = left[valid] != right[valid]
    return int(np.count_nonzero(unequal)), maximum


def compare_values(report: dict, name: str, left, right, atol: float) -> None:
    if getattr(left, "shape", None) != getattr(right, "shape", None):
        _check(report, name, False, reason="shape_mismatch", baseline_shape=list(np.shape(left)),
               candidate_shape=list(np.shape(right)))
        return
    if sparse.issparse(left) or sparse.issparse(right):
        # Bound memory for pairwise graphs; sparse zero and dense zero are equal.
        mismatches, maximum = 0, 0.0
        for start in range(0, left.shape[0], 128):
            a, b = left[start:start + 128], right[start:start + 128]
            a = a.toarray() if sparse.issparse(a) else a
            b = b.toarray() if sparse.issparse(b) else b
            count, error = _equal_arrays(a, b, atol)
            mismatches += count
            maximum = max(maximum, error or 0.0)
    else:
        mismatches, maximum = _equal_arrays(left, right, atol)
    _check(report, name, mismatches == 0, mismatched_values=mismatches,
           max_absolute_difference=maximum)


def _derived_missing(report: dict, name: str) -> None:
    report["derived_differences"].append({"name": name, "change": "missing_in_candidate"})


def compare_derived_values(report: dict, name: str, left, right, atol: float) -> None:
    """Record an allowed refit difference without turning it into a parity pass."""
    comparison = {"checks": []}
    compare_values(comparison, name, left, right, atol)
    result = comparison["checks"][0]
    if not result["passed"]:
        report["derived_differences"].append({
            "name": name, "change": "changed",
            **{key: value for key, value in result.items() if key not in {"name", "passed"}},
        })


def compare_derived_metadata(report: dict, name: str, left, right, atol: float) -> None:
    if isinstance(left, dict):
        if not isinstance(right, dict):
            report["derived_differences"].append({"name": name, "change": "type_changed"})
            return
        for key, value in left.items():
            child = name + "/" + str(key)
            if key not in right:
                _derived_missing(report, child)
            else:
                compare_derived_metadata(report, child, value, right[key], atol)
    else:
        compare_derived_values(report, name, left, right, atol)


def compare_frame_values(report: dict, name: str, left: pd.DataFrame,
                         right: pd.DataFrame, atol: float,
                         derived_columns: frozenset[str] = frozenset()) -> None:
    missing = sorted(set(left.columns) - set(right.columns) - derived_columns)
    _check(report, name + "/required_columns", not missing, missing_columns=missing)
    for column in left.columns:
        if column in derived_columns:
            if column not in right:
                _derived_missing(report, name + "/" + str(column))
            else:
                compare_derived_values(report, name + "/" + str(column),
                                       left[column].to_numpy(), right[column].to_numpy(), atol)
        elif column in right:
            compare_values(report, name + "/" + str(column), left[column].to_numpy(),
                           right[column].to_numpy(), atol)


def align_original_rows(report: dict, name: str, left: pd.DataFrame,
                        right: pd.DataFrame, donors: set[str], keys: list[str]):
    """Return row positions in baseline order, or None after recording failures."""
    missing = sorted(set(keys) - set(right.columns))
    if missing:
        _check(report, name + "/identity", False, reason="missing_identity_columns", missing_columns=missing)
        return None
    donor = keys[0]
    left_donors, right_donors = left[donor].map(donor_id), right[donor].map(donor_id)
    left_positions = np.flatnonzero(left_donors.isin(donors))
    right_positions = np.flatnonzero(right_donors.isin(donors))
    a, b = left.iloc[left_positions][keys].copy(), right.iloc[right_positions][keys].copy()
    a[donor], b[donor] = left_donors.iloc[left_positions].to_numpy(), right_donors.iloc[right_positions].to_numpy()
    if a.isna().any().any() or b.isna().any().any():
        _check(report, name + "/identity", False, reason="missing_identity_values")
        return None
    # Non-donor identifiers retain their literal string representation.
    a_ids = pd.MultiIndex.from_frame(a.astype(str))
    b_ids = pd.MultiIndex.from_frame(b.astype(str))
    missing_ids, added_ids = a_ids.difference(b_ids), b_ids.difference(a_ids)
    duplicates_a, duplicates_b = int(a_ids.duplicated().sum()), int(b_ids.duplicated().sum())
    valid = not(len(missing_ids) or len(added_ids) or duplicates_a or duplicates_b)
    _check(report, name + "/identity", valid,
           baseline_rows=len(a), candidate_original_rows=len(b),
           missing_original_ids=len(missing_ids), added_original_donor_ids=len(added_ids),
           baseline_duplicate_ids=duplicates_a, candidate_duplicate_ids=duplicates_b,
           missing_examples=[list(x) for x in missing_ids[:5]], added_examples=[list(x) for x in added_ids[:5]])
    if not valid:
        return None
    return left_positions, right_positions[b_ids.get_indexer(a_ids)]


def compare_table(report: dict, name: str, left: pd.DataFrame, right: pd.DataFrame,
                  donors: set[str], atol: float) -> None:
    keys = table_keys(left)
    aligned = align_original_rows(report, name, left, right, donors, keys)
    if aligned is not None:
        a, b = aligned
        compare_frame_values(report, name, left.iloc[a], right.iloc[b], atol)


def _uns_tables(adata) -> dict[str, pd.DataFrame]:
    tables = {str(k): v for k, v in adata.uns.items() if isinstance(v, pd.DataFrame)}
    for key, columns in adata.uns.items():
        if key.endswith("_columns"):
            prefix = key[:-len("_columns")]
            tables[prefix] = pd.DataFrame({str(c): adata.uns[prefix + "_" + str(c)] for c in columns})
    return tables


def compare_metadata(report: dict, name: str, left, right, atol: float) -> None:
    if isinstance(left, dict):
        if not isinstance(right, dict):
            _check(report, name, False, reason="metadata_mapping_missing")
            return
        for key, value in left.items():
            _check(report, name + "/" + str(key) + "/present", key in right)
            if key in right:
                compare_metadata(report, name + "/" + str(key), value, right[key], atol)
    else:
        compare_values(report, name, left, right, atol)


def compare_h5ad(report: dict, name: str, left, right, donors: set[str], atol: float) -> None:
    donor = donor_column(left.obs)
    identity = "combined_islet_id" if "combined_islet_id" in left.obs else "islet_id"
    if identity not in left.obs or not donor:
        raise ValueError("Canonical H5AD lacks donor/islet identity")
    if identity in right.obs:
        unique = left.obs[identity].is_unique and right.obs[identity].is_unique
        _check(report, name + "/global_islet_identity", unique)
        if not unique:
            return
    aligned = align_original_rows(report, name + "/obs", left.obs, right.obs, donors, [donor, identity])
    if aligned is None:
        return
    a, b = aligned
    compare_frame_values(report, name + "/obs", left.obs.iloc[a], right.obs.iloc[b],
                         atol, DERIVED_OBS_COLUMNS)
    compare_values(report, name + "/obs_names", left.obs_names.to_numpy()[a], right.obs_names.to_numpy()[b], atol)
    vars_ok = left.var_names.is_unique and right.var_names.is_unique and left.var_names.isin(right.var_names).all()
    _check(report, name + "/var_identity", vars_ok,
           missing_variables=left.var_names.difference(right.var_names).tolist(),
           duplicate_variables=int(right.var_names.duplicated().sum()))
    if vars_ok:
        v = right.var_names.get_indexer(left.var_names)
        compare_frame_values(report, name + "/var", left.var, right.var.iloc[v], atol)
        compare_values(report, name + "/X", left.X[a], right.X[b][:, v], atol)
        for key in left.layers:
            _check(report, name + "/layers/" + key + "/present", key in right.layers)
            if key in right.layers:
                compare_values(report, name + "/layers/" + key, left.layers[key][a], right.layers[key][b][:, v], atol)
        for slot in ("varm", "varp"):
            x, y = getattr(left, slot), getattr(right, slot)
            for key in x:
                field = name + "/" + slot + "/" + key
                derived = slot == "varm" and key in DERIVED_VARM_KEYS
                if key not in y:
                    if derived:
                        _derived_missing(report, field)
                    else:
                        _check(report, field + "/present", False)
                else:
                    if not derived:
                        _check(report, field + "/present", True)
                    candidate = y[key][v]
                    if slot == "varp":
                        candidate = candidate[:, v]
                    if derived:
                        compare_derived_values(report, field, x[key], candidate, atol)
                    else:
                        compare_values(report, field, x[key], candidate, atol)
    if left.raw is not None:
        _check(report, name + "/raw/present", right.raw is not None)
        if right.raw is not None:
            # Raw has its own variables but the same observation identities.
            compare_h5ad(report, name + "/raw", left.raw.to_adata(), right.raw.to_adata(), donors, atol)
    for slot in ("obsm", "obsp"):
        x, y = getattr(left, slot), getattr(right, slot)
        for key in x:
            field = name + "/" + slot + "/" + key
            derived = key in (DERIVED_OBSM_KEYS if slot == "obsm" else DERIVED_OBSP_KEYS)
            if key not in y:
                if derived:
                    _derived_missing(report, field)
                else:
                    _check(report, field + "/present", False)
            else:
                if not derived:
                    _check(report, field + "/present", True)
                left_value, right_value = x[key][a], y[key][b]
                if slot == "obsp":
                    left_value, right_value = left_value[:, a], right_value[:, b]
                if derived:
                    compare_derived_values(report, field, left_value, right_value, atol)
                else:
                    compare_values(report, field, left_value, right_value, atol)
    original_tables, candidate_tables = _uns_tables(left), _uns_tables(right)
    for table, values in original_tables.items():
        _check(report, name + "/uns/" + table + "/present", table in candidate_tables)
        if table in candidate_tables:
            compare_table(report, name + "/uns/" + table, values, candidate_tables[table], donors, atol)
    table_keys_used = set(original_tables)
    for table, frame in original_tables.items():
        table_keys_used.add(table + "_columns")
        table_keys_used.update(table + "_" + str(column) for column in frame.columns)
    for key, value in left.uns.items():
        if key in table_keys_used:
            continue
        field = name + "/uns/" + key
        derived = key in DERIVED_UNS_KEYS
        if key not in right.uns:
            if derived:
                _derived_missing(report, field)
            else:
                _check(report, field + "/present", False)
            continue
        if not derived:
            _check(report, field + "/present", True)
        if key == "iroot":
            # iroot is a row index: compare the identified islet after reordering.
            try:
                original_root = np.asarray(left.obs[identity].iloc[int(value)])
                candidate_root = np.asarray(right.obs[identity].iloc[int(right.uns[key])])
                compare_derived_values(report, name + "/uns/iroot_islet_identity",
                                       original_root, candidate_root, atol)
            except (IndexError, TypeError, ValueError) as error:
                report["derived_differences"].append({
                    "name": field, "change": "invalid_root_index", "reason": str(error)})
        elif derived:
            compare_derived_metadata(report, field, value, right.uns[key], atol)
        else:
            compare_metadata(report, field, value, right.uns[key], atol)


def _candidate_census(report: dict, baseline, candidate_root: Path, donors: set[str]) -> None:
    """Describe schema-4 IDs without claiming feature equivalence."""
    candidates = (candidate_root / "trajectory/composition_features.parquet",
                  candidate_root / "trajectory/composition_features.csv")
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        return
    frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    column = donor_column(frame)
    if not column or "islet_id" not in frame:
        raise ValueError("Native candidate feature table lacks donor/islet IDs")
    bd = baseline.obs[donor_column(baseline.obs)].map(donor_id)
    cd = frame[column].map(donor_id)
    original_id = "combined_islet_id" if "combined_islet_id" in baseline.obs else "islet_id"
    per_donor = []
    for donor in sorted(donors):
        old = set(baseline.obs.loc[bd == donor, original_id].astype(str))
        new = set(frame.loc[cd == donor, "islet_id"].astype(str))
        per_donor.append({"donor": donor, "baseline_ids": len(old), "candidate_ids": len(new),
                          "shared_ids": len(old & new), "missing_original_ids": len(old - new),
                          "added_ids": len(new - old)})
    original_rows = frame.loc[cd.isin(donors)]
    duplicate_ids = int(original_rows["islet_id"].duplicated().sum())
    missing_ids = int(original_rows["islet_id"].isna().sum())
    preserved = not(duplicate_ids or missing_ids) and all(
        row["missing_original_ids"] == row["added_ids"] == 0 for row in per_donor)
    _check(report, "native_candidate_feature_identity", preserved,
           reason="Native annotations must retain original identities; no inferred remapping is permitted",
           source=str(path), baseline_original_rows=int(bd.isin(donors).sum()),
           candidate_original_rows=int(cd.isin(donors).sum()),
           candidate_duplicate_ids=duplicate_ids, candidate_missing_ids=missing_ids, per_donor=per_donor)


def audit(baseline_root, candidate_root, atol: float = 0.0) -> dict:
    """Read-only islet-artifact gate; ``passed`` does not certify cell-level parity."""
    baseline_root, candidate_root = Path(baseline_root).resolve(), Path(candidate_root).resolve()
    report = {"schema_version": 2, "baseline_root": str(baseline_root),
              "candidate_root": str(candidate_root), "numeric_atol": atol, "numeric_rtol": 0,
              "scope": "islet_level_canonical_artifacts",
              "source_cell_parity_checked": False,
              "source_cell_parity_requirement": (
                  "Check original single-cell IDs, measurements, missingness and labels "
                  "separately before 21-donor integration"),
              "checks": [], "derived_differences": [], "derived_difference_count": 0,
              "passed": False}
    if not math.isfinite(atol) or atol < 0:
        _check(report, "numeric_tolerance", False, reason="atol must be finite and nonnegative")
        return report
    try:
        baseline = ad.read_h5ad(baseline_root / "adata_ins_root.h5ad")
        column = donor_column(baseline.obs)
        if not column:
            raise ValueError("Baseline H5AD lacks donor identity")
        donors = set(baseline.obs[column].map(donor_id))
        report["original_donors"] = sorted(donors)
        _check(report, "baseline_cohort_nonempty", bool(donors))
        for filename in ARTIFACTS:
            original, candidate = baseline_root / filename, candidate_root / filename
            if not original.is_file():
                _check(report, filename + "/baseline_present", False, reason="Required baseline artifact absent")
                continue
            if not candidate.is_file():
                _check(report, filename + "/candidate_present", False, status="incomparable",
                       reason="No canonical analogue; schema-4 features cannot substitute for frozen values")
                continue
            try:
                if filename.endswith(".h5ad"):
                    left = baseline if filename == "adata_ins_root.h5ad" else ad.read_h5ad(original)
                    compare_h5ad(report, filename, left, ad.read_h5ad(candidate), donors, atol)
                else:
                    compare_table(report, filename, pd.read_csv(original), pd.read_csv(candidate), donors, atol)
            except Exception as error:
                _check(report, filename + "/comparison", False, status="incomparable", reason=str(error))
        original, candidate = baseline_root / "master_results.xlsx", candidate_root / "master_results.xlsx"
        if candidate.is_file():
            try:
                x, y = pd.read_excel(original, sheet_name=None), pd.read_excel(candidate, sheet_name=None)
                for name, frame in x.items():
                    _check(report, "master_results.xlsx/" + name + "/present", name in y)
                    if name in y:
                        compare_table(report, "master_results.xlsx/" + name, frame, y[name], donors, atol)
            except Exception as error:
                _check(report, "master_results.xlsx/comparison", False, reason=str(error))
        else:
            report["excel"] = "Candidate has no workbook; mandatory canonical H5AD uns tables cover original rows"
        try:
            _candidate_census(report, baseline, candidate_root, donors)
        except Exception as error:
            _check(report, "native_candidate_feature_identity", False, reason=str(error))
    except Exception as error:
        _check(report, "baseline_readable", False, reason=str(error))
    report["passed"] = bool(report["checks"]) and all(c["passed"] for c in report["checks"])
    report["derived_difference_count"] = len(report["derived_differences"])
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", required=True, type=Path)
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--report-root", required=True, type=Path)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args(argv)
    destination = args.report_root.resolve()
    forbidden = [args.baseline_root.resolve(), args.candidate_root.resolve(),
                 Path(__file__).resolve().parents[2] / "data"]
    if any(destination.is_relative_to(root) for root in forbidden):
        parser.error("--report-root must be outside app data and both source roots")
    report = audit(args.baseline_root, args.candidate_root, args.atol)
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / "original_subset_parity.json"
    with path.open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"passed": report["passed"], "checks": len(report["checks"]),
                      "failures": sum(not c["passed"] for c in report["checks"]),
                      "derived_difference_count": report["derived_difference_count"],
                      "source_cell_parity_checked": report["source_cell_parity_checked"],
                      "report": str(path)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
