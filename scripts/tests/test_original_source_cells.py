"""Small H5AD fixtures for exact, row-order-independent source-cell parity."""
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))
from audit_original_source_cells import audit, main  # noqa: E402


CELLS = ((1, "cell-a"), (1, "cell-b"), (2, "cell-c"), (3, "cell-d"))
VALUES = {
    "cell-a": {"INS": 1.0, "GCG": np.nan, "x": 10.0, "phenotype": "Beta cell"},
    "cell-b": {"INS": 2.0, "GCG": 3.0, "x": 20.0, "phenotype": "Alpha cell"},
    "cell-c": {"INS": 4.0, "GCG": 5.0, "x": 30.0, "phenotype": "Beta cell"},
    "cell-d": {"INS": 6.0, "GCG": 7.0, "x": 40.0, "phenotype": "Alpha cell"},
}


def write_fixture(path, cells=CELLS[:3], markers=("INS", "GCG"),
                  phenotype_categories=("Beta cell", "Alpha cell")):
    obs = pd.DataFrame({
        "imageid": [donor for donor, _ in cells],
        "Object ID": [cell for _, cell in cells],
        "X_centroid": [VALUES[cell]["x"] for _, cell in cells],
        "Y_centroid": [VALUES[cell]["x"] + 1 for _, cell in cells],
        "Cell Area": [VALUES[cell]["x"] * 2 for _, cell in cells],
        "nullable_measurement": [np.nan if cell == "cell-a" else 1.0 for _, cell in cells],
        "phenotype": pd.Categorical([VALUES[cell]["phenotype"] for _, cell in cells],
                                    categories=phenotype_categories),
        "Donor Status": pd.Categorical(["ND" if donor == 1 else "T1D"
                                        for donor, _ in cells]),
        "_scvi_batch": [donor - 1 for donor, _ in cells],
        "_scvi_labels": [0 for _ in cells],
    }, index=[str(i + 100) for i in range(len(cells))])
    var = pd.DataFrame(index=list(markers))
    x = np.asarray([[VALUES[cell][marker] for marker in markers]
                    for _, cell in cells], dtype=float)
    value = ad.AnnData(X=x.copy(), obs=obs, var=var)
    value.layers["log"] = x + 0.5
    value.raw = ad.AnnData(X=x * 10, obs=obs.copy(), var=var.copy())
    value.obsm["X_scVI"] = np.asarray([[VALUES[cell]["x"], 1.0]
                                       for _, cell in cells], dtype=np.float32)
    value.obsm["X_umap"] = np.asarray([[VALUES[cell]["x"], 2.0]
                                       for _, cell in cells], dtype=np.float32)
    value.uns["gates"] = pd.DataFrame({"1": [1.0, 2.0], "2": [3.0, 4.0],
                                      "3": [5.0, 6.0]},
                                     index=pd.Index(["INS", "GCG"], name="markers"))
    value.write_h5ad(path)


@pytest.fixture
def paths(tmp_path):
    baseline, candidate = tmp_path / "baseline.h5ad", tmp_path / "candidate.h5ad"
    write_fixture(baseline)
    write_fixture(candidate)
    return baseline, candidate


def mutate(path, callback):
    with h5py.File(path, "r+") as file:
        callback(file)


def failures(report):
    return {item["name"] for item in report["checks"] if not item["passed"]}


def test_reordered_original_cells_added_donor_and_markers_pass(paths):
    baseline, candidate = paths
    write_fixture(candidate, cells=(CELLS[2], CELLS[0], CELLS[3], CELLS[1]),
                  markers=("GCG", "INS"),
                  phenotype_categories=("Alpha cell", "Beta cell"))
    report = audit(baseline, candidate)
    assert report["passed"], report
    assert report["source_cell_parity_checked"]
    assert report["derived_difference_count"] == 0
    identity = next(item for item in report["checks"]
                    if item["name"] == "original_cell_identity")
    assert identity["candidate_new_donor_cells"] == 1
    assert report["synthetic_index_compared"] is False


@pytest.mark.parametrize("field, expected", [
    ("X", "X"), ("layers/log", "layers/log"), ("raw/X", "raw/X"),
])
def test_changed_original_measurement_fails(paths, field, expected):
    baseline, candidate = paths
    mutate(candidate, lambda file: file[field].__setitem__((0, 0), 99.0))
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert expected in failures(report)


def test_missingness_change_fails(paths):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["raw/X"].__setitem__((0, 1), 5.0))
    report = audit(baseline, candidate)
    assert not report["passed"]
    changed = next(item for item in report["checks"] if item["name"] == "raw/X")
    assert changed["changed_missingness"] == 1


@pytest.mark.parametrize("field, value", [
    ("obs/X_centroid", 999.0), ("obs/Cell Area", 999.0),
    ("obs/nullable_measurement", 5.0),
])
def test_changed_geometry_or_source_obs_fails(paths, field, value):
    baseline, candidate = paths
    mutate(candidate, lambda file: file[field].__setitem__(0, value))
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert field in failures(report)


def test_changed_phenotype_label_fails(paths):
    baseline, candidate = paths
    def change(file):
        file["obs/phenotype/codes"][0] = 1
    mutate(candidate, change)
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert "obs/phenotype" in failures(report)


def test_missing_protected_obs_field_fails_as_incomplete(paths):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["obs"].__delitem__("Cell Area"))
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert not report["source_cell_parity_checked"]
    assert "obs/Cell Area" in failures(report)


def test_missing_or_added_original_cell_fails(paths):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["obs/Object ID"].__setitem__(1, "cell-c"))
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert "original_cell_identity" in failures(report)
    assert not report["source_cell_parity_checked"]


def test_joint_refit_coordinates_and_registry_are_reported(paths):
    baseline, candidate = paths
    def change(file):
        file["obsm/X_scVI"][0, 0] = 99
        file["obsm/X_umap"][0, 0] = 88
        file["obs/_scvi_batch"][0] = 9
    mutate(candidate, change)
    report = audit(baseline, candidate)
    assert report["passed"], report
    assert {item["name"] for item in report["derived_differences"]} == {
        "obsm/X_scVI", "obsm/X_umap", "obs/_scvi_batch"}
    assert report["derived_difference_count"] == 3


def test_missing_derived_coordinate_is_reported(paths):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["obsm"].__delitem__("X_umap"))
    report = audit(baseline, candidate)
    assert report["passed"], report
    assert {"name": "obsm/X_umap", "change": "missing_in_candidate"} in report[
        "derived_differences"]


def test_original_numeric_gate_change_fails(paths):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["uns/gates/1"].__setitem__(0, 9.0))
    report = audit(baseline, candidate)
    assert not report["passed"]
    assert "uns/gates/1" in failures(report)


def test_cli_writes_explicit_report_outside_source(paths, tmp_path):
    baseline, candidate = paths
    mutate(candidate, lambda file: file["obsm/X_scVI"].__setitem__((0, 0), 99))
    output = tmp_path / "report" / "cells.json"
    summary = StringIO()
    with redirect_stdout(summary):
        assert main(["--baseline", str(baseline), "--candidate", str(candidate),
                     "--report", str(output)]) == 0
    result = json.loads(summary.getvalue())
    saved = json.loads(output.read_text())
    assert result["derived_difference_count"] == saved["derived_difference_count"] == 1
    assert result["source_cell_parity_checked"] is True
    assert saved["passed"]
