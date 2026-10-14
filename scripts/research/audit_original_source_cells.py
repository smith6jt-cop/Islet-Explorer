#!/usr/bin/env python3
"""Read-only, exact source-cell parity for the frozen 15-donor CODEX H5AD.

Cells are matched by (obs/imageid, obs/Object ID), not by the synthetic
obs/_index or row order. Original cells, marker values, missingness, obs source
fields and phenotype labels must be preserved. A single-cell refit may
change X_scVI, X_umap and the two scVI registry columns; those differences are
reported separately. This checks no islet-level app artifacts or new-donor
typing quality. It never loads a whole H5AD through AnnData: each source field
is scanned in chunks, with at most one original-cohort field aligned in memory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


REPO = Path(__file__).resolve().parents[2]
DEFAULT_BASELINE = REPO / "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
IDENTITY_COLUMNS = ("imageid", "Object ID")
DERIVED_OBS = frozenset({"_scvi_batch", "_scvi_labels"})
DERIVED_COORDINATES = ("X_scVI", "X_umap")
REQUIRED_BASELINE_PATHS = (
    "X", "obs/imageid", "obs/Object ID", "obs/X_centroid", "obs/Y_centroid",
    "obs/phenotype", "layers/log", "raw/X", "uns/gates", "obsm/X_scVI",
    "obsm/X_umap",
)
CHUNK_ROWS = 16384


def _text(value) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("utf-8")
    return str(value)


def _bytes(value) -> bytes:
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value)
    if isinstance(value, str):
        return value.encode("utf-8")
    raise ValueError(f"Expected a cell/marker string, found {type(value).__name__}")


def _check(report: dict, name: str, passed: bool, **details) -> None:
    report["checks"].append({"name": name, "passed": bool(passed), **details})


def _difference(report: dict, name: str, **details) -> None:
    report["derived_differences"].append({"name": name, **details})


def _identities(file: h5py.File, donors: set[int] | None = None):
    obs = file["obs"]
    for column in IDENTITY_COLUMNS:
        if column not in obs:
            raise ValueError(f"Missing identity column obs/{column}")
    donor_dataset = obs["imageid"]
    object_dataset = obs["Object ID"]
    if not isinstance(donor_dataset, h5py.Dataset) or donor_dataset.dtype.kind not in "iu":
        raise ValueError("obs/imageid must be an integer dataset")
    if not isinstance(object_dataset, h5py.Dataset) or object_dataset.ndim != 1:
        raise ValueError("obs/Object ID must be a string dataset")
    if donor_dataset.shape != object_dataset.shape:
        raise ValueError("Identity column lengths differ")
    donor_values = donor_dataset[:].astype(np.int64, copy=False)
    selected = np.arange(len(donor_values), dtype=np.int64) if donors is None else np.flatnonzero(
        np.isin(donor_values, list(donors)))
    selected_donors = donor_values[selected]
    # QuPath UUIDs are variable-length HDF5 strings. Retain only selected
    # original-donor IDs, then pack them into a fixed-width sortable array.
    blocks = []
    for start in range(0, len(donor_values), CHUNK_ROWS):
        end = min(start + CHUNK_ROWS, len(donor_values))
        mask = None if donors is None else np.isin(donor_values[start:end], list(donors))
        if mask is not None and not mask.any():
            continue
        values = object_dataset[start:end]
        blocks.extend(_bytes(value) for value in (values if mask is None else values[mask]))
    if len(blocks) != len(selected) or any(not value for value in blocks):
        raise ValueError("Missing or invalid obs/Object ID")
    width = max(map(len, blocks), default=1)
    objects = np.asarray(blocks, dtype=f"S{width}")
    return selected, selected_donors, objects, len(donor_values)


def _identity_label(key) -> str:
    return f"{int(key['donor'])}:{_text(key['object'])}"


def _alignment(report: dict, baseline: h5py.File, candidate: h5py.File):
    b_pos, b_donor, b_object, b_total = _identities(baseline)
    if not b_total:
        raise ValueError("Baseline has no source cells")
    donors = set(int(value) for value in np.unique(b_donor))
    c_pos, c_donor, c_object, c_total = _identities(candidate, donors)
    width = max(b_object.dtype.itemsize, c_object.dtype.itemsize)
    dtype = np.dtype([("donor", "<i8"), ("object", f"S{width}")])
    b_keys, c_keys = np.empty(len(b_pos), dtype=dtype), np.empty(len(c_pos), dtype=dtype)
    b_keys["donor"], b_keys["object"] = b_donor, b_object
    c_keys["donor"], c_keys["object"] = c_donor, c_object
    b_order, c_order = np.argsort(b_keys), np.argsort(c_keys)
    b_sorted, c_sorted = b_keys[b_order], c_keys[c_order]
    b_duplicates = int(np.count_nonzero(b_sorted[1:] == b_sorted[:-1]))
    c_duplicates = int(np.count_nonzero(c_sorted[1:] == c_sorted[:-1]))
    missing = np.setdiff1d(b_sorted, c_sorted)
    added = np.setdiff1d(c_sorted, b_sorted)
    valid = not (b_duplicates or c_duplicates or len(missing) or len(added))
    _check(report, "original_cell_identity", valid,
           baseline_cells=b_total, candidate_cells=c_total,
           original_donors=sorted(donors), candidate_original_donor_cells=len(c_pos),
           candidate_new_donor_cells=c_total - len(c_pos),
           baseline_duplicate_ids=b_duplicates, candidate_duplicate_ids=c_duplicates,
           missing_original_ids=len(missing), added_original_donor_ids=len(added),
           missing_examples=[_identity_label(key) for key in missing[:5]],
           added_examples=[_identity_label(key) for key in added[:5]])
    if not valid:
        return None
    # Row map lets each candidate field be scanned sequentially. Candidate
    # values for original cells scatter into baseline order; no random HDF5 IO.
    candidate_to_baseline = np.full(c_total, -1, dtype=np.int32)
    candidate_to_baseline[c_pos[c_order]] = b_order.astype(np.int32)
    return b_keys, candidate_to_baseline


def _marker_map(report: dict, name: str, baseline, candidate):
    if "var" not in baseline or "_index" not in baseline["var"]:
        raise ValueError(f"Baseline {name} marker index missing")
    if "var" not in candidate or "_index" not in candidate["var"]:
        _check(report, name + "/marker_identity", False, reason="candidate_marker_index_missing")
        return None
    left = [_bytes(value) for value in baseline["var/_index"][:]]
    right = [_bytes(value) for value in candidate["var/_index"][:]]
    duplicates = len(right) - len(set(right))
    missing = sorted(set(left) - set(right))
    valid = len(left) == len(set(left)) and not duplicates and not missing
    _check(report, name + "/marker_identity", valid,
           baseline_markers=len(left), candidate_markers=len(right),
           candidate_duplicate_markers=duplicates,
           missing_markers=[_text(value) for value in missing[:10]])
    return [right.index(value) for value in left] if valid else None


def _count_numeric(left: np.ndarray, right: np.ndarray):
    if left.dtype.kind not in "biuf" or right.dtype.kind not in "biuf":
        raise ValueError("Expected numeric source data")
    left_missing = np.isnan(left) if left.dtype.kind == "f" else np.zeros(left.shape, dtype=bool)
    right_missing = np.isnan(right) if right.dtype.kind == "f" else np.zeros(right.shape, dtype=bool)
    mismatch = (left_missing != right_missing) | (
        ~left_missing & ~right_missing & (left != right))
    return mismatch, int(np.count_nonzero(left_missing != right_missing))


def _record_comparison(report: dict, name: str, mismatch_count: int,
                       missing_count: int, examples: list[str], derived: bool) -> None:
    if derived:
        if mismatch_count:
            _difference(report, name, change="changed", mismatched_values=mismatch_count,
                        changed_missingness=missing_count, example_original_cells=examples)
    else:
        _check(report, name, mismatch_count == 0, mismatched_values=mismatch_count,
               changed_missingness=missing_count, example_original_cells=examples)


def _compare_rows(report: dict, name: str, left: h5py.Dataset, right: h5py.Dataset,
                  row_map: np.ndarray, keys: np.ndarray, *, columns=None,
                  derived: bool = False, baseline_transform=None,
                  candidate_transform=None, output_dtype=None) -> None:
    if not isinstance(left, h5py.Dataset) or not isinstance(right, h5py.Dataset):
        _check(report, name, False, reason="unsupported_non_dense_source_dataset")
        return
    if left.shape[0] != len(keys) or right.shape[0] != len(row_map):
        _check(report, name, False, reason="row_count_mismatch")
        return
    tail = left.shape[1:]
    if right.ndim != left.ndim or (columns is None and right.shape[1:] != tail):
        if derived:
            _difference(report, name, change="shape_changed",
                        baseline_shape=list(left.shape), candidate_shape=list(right.shape))
        else:
            _check(report, name, False, reason="shape_mismatch",
                   baseline_shape=list(left.shape), candidate_shape=list(right.shape))
        return
    if columns is not None and (left.ndim != 2 or len(columns) != left.shape[1]
                                or max(columns, default=-1) >= right.shape[1]):
        _check(report, name, False, reason="marker_dimension_mismatch")
        return
    dtype = output_dtype or right.dtype
    if np.dtype(dtype).kind not in "biuf":
        _check(report, name, False, reason="unsupported_source_dtype")
        return
    aligned = np.empty(left.shape, dtype=dtype)
    for start in range(0, len(row_map), CHUNK_ROWS):
        end = min(start + CHUNK_ROWS, len(row_map))
        target = row_map[start:end]
        mask = target >= 0
        if not mask.any():
            continue
        block = right[start:end]
        if columns is not None:
            block = block[:, columns]
        if candidate_transform is not None:
            block = candidate_transform(block)
        aligned[target[mask]] = block[mask]
    mismatches = missing_changes = 0
    examples: list[str] = []
    for start in range(0, len(keys), CHUNK_ROWS):
        end = min(start + CHUNK_ROWS, len(keys))
        block = left[start:end]
        if baseline_transform is not None:
            block = baseline_transform(block)
        different, missing = _count_numeric(block, aligned[start:end])
        mismatches += int(np.count_nonzero(different))
        missing_changes += missing
        if len(examples) < 5 and different.any():
            rows = np.flatnonzero(different.reshape(len(block), -1).any(axis=1))
            examples.extend(_identity_label(keys[start + row]) for row in rows[:5 - len(examples)])
    _record_comparison(report, name, mismatches, missing_changes, examples, derived)


def _categorical_map(baseline: h5py.Group, candidate: h5py.Group):
    if not all(key in baseline for key in ("categories", "codes")) or not all(
            key in candidate for key in ("categories", "codes")):
        raise ValueError("Categorical codes/categories missing")
    left = [_bytes(value) for value in baseline["categories"][:]]
    right = [_bytes(value) for value in candidate["categories"][:]]
    if len(left) != len(set(left)) or len(right) != len(set(right)):
        raise ValueError("Duplicate categorical labels")
    lookup = {label: code for code, label in enumerate(left)}
    mapping = np.asarray([lookup.get(label, -2) for label in right], dtype=np.int32)

    def map_codes(codes):
        if np.any(codes < -1) or np.any(codes >= len(mapping)):
            raise ValueError("Invalid categorical code")
        return np.where(codes == -1, -1, mapping[np.maximum(codes, 0)])

    return map_codes


def _compare_obs(report: dict, baseline: h5py.File, candidate: h5py.File,
                 row_map: np.ndarray, keys: np.ndarray) -> None:
    for column, left in baseline["obs"].items():
        if column in {"_index", *IDENTITY_COLUMNS}:
            continue
        name = "obs/" + column
        derived = column in DERIVED_OBS
        if column not in candidate["obs"]:
            if derived:
                _difference(report, name, change="missing_in_candidate")
            else:
                _check(report, name, False, reason="missing_in_candidate")
            continue
        right = candidate["obs"][column]
        try:
            if isinstance(left, h5py.Group) and isinstance(right, h5py.Group):
                mapping = _categorical_map(left, right)
                _compare_rows(report, name, left["codes"], right["codes"], row_map, keys,
                              derived=derived, candidate_transform=mapping, output_dtype=np.int32)
            elif isinstance(left, h5py.Dataset) and isinstance(right, h5py.Dataset):
                _compare_rows(report, name, left, right, row_map, keys, derived=derived)
            else:
                if derived:
                    _difference(report, name, change="encoding_changed")
                else:
                    _check(report, name, False, reason="column_encoding_changed")
        except Exception as error:
            if derived:
                _difference(report, name, change="comparison_error", error=str(error))
            else:
                _check(report, name, False, reason="comparison_error", error=str(error))


def _compare_matrix(report: dict, name: str, baseline: h5py.File, candidate: h5py.File,
                    row_map: np.ndarray, keys: np.ndarray, *, columns=None,
                    derived: bool = False) -> None:
    if name not in baseline:
        return
    if name not in candidate:
        if derived:
            _difference(report, name, change="missing_in_candidate")
        else:
            _check(report, name, False, reason="missing_in_candidate")
        return
    try:
        _compare_rows(report, name, baseline[name], candidate[name], row_map, keys,
                      columns=columns, derived=derived)
    except Exception as error:
        _check(report, name, False, reason="comparison_error", error=str(error))


def _compare_gates(report: dict, baseline: h5py.File, candidate: h5py.File) -> None:
    if "uns/gates" not in baseline:
        return
    if "uns/gates" not in candidate:
        _check(report, "uns/gates", False, reason="missing_in_candidate")
        return
    left, right = baseline["uns/gates"], candidate["uns/gates"]
    if not isinstance(left, h5py.Group) or not isinstance(right, h5py.Group):
        _check(report, "uns/gates", False, reason="gate_table_encoding_changed")
        return
    try:
        marker_name = _text(left.attrs["_index"])
        candidate_marker_name = _text(right.attrs["_index"])
        left_markers = [_bytes(value) for value in left[marker_name][:]]
        right_markers = [_bytes(value) for value in right[candidate_marker_name][:]]
        if len(right_markers) != len(set(right_markers)) or not set(left_markers) <= set(right_markers):
            _check(report, "uns/gates/marker_identity", False, reason="missing_or_duplicate_markers")
            return
        _check(report, "uns/gates/marker_identity", True)
        order = [right_markers.index(value) for value in left_markers]
        for donor in left.attrs["column-order"]:
            name = _text(donor)
            if name not in right:
                _check(report, "uns/gates/" + name, False, reason="missing_original_donor_gate")
                continue
            a, b = left[name][:], right[name][:][order]
            mismatch, missing = _count_numeric(a, b)
            _record_comparison(report, "uns/gates/" + name,
                               int(np.count_nonzero(mismatch)), missing, [], False)
    except Exception as error:
        _check(report, "uns/gates", False, reason="comparison_error", error=str(error))


def audit(baseline_path: Path, candidate_path: Path) -> dict:
    """Compare a future candidate to frozen single-cell source records, read-only."""
    baseline_path, candidate_path = Path(baseline_path).resolve(), Path(candidate_path).resolve()
    report = {"schema_version": 1, "scope": "original_single_cell_source_records",
              "baseline": str(baseline_path), "candidate": str(candidate_path),
              "numeric_atol": 0, "numeric_rtol": 0,
              "identity_columns": list(IDENTITY_COLUMNS),
              "synthetic_index_compared": False,
              "joint_scvi_refit_allowed": True,
              "derived_coordinates": list(DERIVED_COORDINATES),
              "checks": [], "derived_differences": [],
              "derived_difference_count": 0,
              "source_cell_parity_checked": False, "passed": False}
    try:
        with h5py.File(baseline_path, "r") as baseline, h5py.File(candidate_path, "r") as candidate:
            missing_baseline = [name for name in REQUIRED_BASELINE_PATHS if name not in baseline]
            if missing_baseline:
                raise ValueError(f"Baseline lacks frozen source fields: {missing_baseline}")
            aligned = _alignment(report, baseline, candidate)
            if aligned is not None:
                keys, row_map = aligned
                _compare_obs(report, baseline, candidate, row_map, keys)
                var_columns = _marker_map(report, "var", baseline, candidate)
                if var_columns is not None:
                    _compare_matrix(report, "X", baseline, candidate, row_map, keys,
                                    columns=var_columns)
                    for layer in baseline.get("layers", {}):
                        _compare_matrix(report, "layers/" + layer, baseline, candidate,
                                        row_map, keys, columns=var_columns)
                if "raw" in baseline:
                    if "raw" not in candidate:
                        _check(report, "raw", False, reason="missing_in_candidate")
                    else:
                        raw_columns = _marker_map(report, "raw/var", baseline["raw"],
                                                  candidate["raw"])
                        if raw_columns is not None:
                            _compare_matrix(report, "raw/X", baseline, candidate, row_map,
                                            keys, columns=raw_columns)
                for coordinate in DERIVED_COORDINATES:
                    _compare_matrix(report, "obsm/" + coordinate, baseline, candidate,
                                    row_map, keys, derived=True)
                _compare_gates(report, baseline, candidate)
                # A measured mismatch is a completed comparison; a missing or
                # malformed source field leaves coverage incomplete.
                report["source_cell_parity_checked"] = all(
                    item["passed"] or "mismatched_values" in item
                    for item in report["checks"])
    except Exception as error:
        _check(report, "audit_readable", False, reason=str(error))
    report["derived_difference_count"] = len(report["derived_differences"])
    report["passed"] = bool(report["checks"]) and all(item["passed"] for item in report["checks"])
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    report_path = args.report.resolve()
    forbidden = (args.baseline.resolve(), args.candidate.resolve())
    if report_path in forbidden or report_path.is_relative_to(REPO / "data"):
        parser.error("--report must not overwrite either H5AD or be inside app data")
    report = audit(args.baseline, args.candidate)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"passed": report["passed"], "checks": len(report["checks"]),
                      "failures": sum(not item["passed"] for item in report["checks"]),
                      "derived_difference_count": report["derived_difference_count"],
                      "source_cell_parity_checked": report["source_cell_parity_checked"],
                      "report": str(report_path)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
