"""The CODEX substitutions in scripts/uniform's islet classifier fit (2026-10-15 canary findings)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "uniform"))
pytest.importorskip("cv2")
pytest.importorskip("shapely")
import islet_cell_test as ict  # noqa: E402
import islet_classifier_fit as fit  # noqa: E402


def test_pgp95_does_not_name_a_lineage():
    assert "PGP9.5" in ict.SHARED_WITH_ENDOCRINE
    assert "PGP9.5" not in ict.IDENTITY
    assert all("PGP9.5" not in markers for markers in ict.LINEAGE_MARKERS.values())
    assert ict.LINEAGE_MARKERS["Neural"] == ("B3TUBB", "GAP43")
    assert ict.LEVEL_MARKER not in ict.IDENTITY


def _tissue(n_acinar: int = 4000, n_beta: int = 40, seed: int = 0) -> pd.DataFrame:
    """Raw counts shaped like canary 6505: E-cadherin two orders of magnitude above insulin."""
    rng = np.random.default_rng(seed)
    acinar = pd.DataFrame({"ECAD": rng.lognormal(np.log(7600), 0.6, n_acinar), "INS": rng.lognormal(np.log(17), 0.3, n_acinar)})
    beta = pd.DataFrame({"ECAD": rng.lognormal(np.log(2700), 0.4, n_beta), "INS": rng.lognormal(np.log(600), 0.2, n_beta)})
    return pd.concat([acinar, beta], ignore_index=True)


def test_common_scale_lets_a_dim_hormone_win():
    values = _tissue()
    beta = np.arange(len(values)) >= 4000
    eligible = np.ones(len(values), dtype=bool)
    assert (ict.top_marker(values)[beta] == "ECAD").all()           # raw counts: insulin never wins
    scaled, top = ict.common_scale(values, eligible)
    assert top["ECAD"] == pytest.approx(np.percentile(values["ECAD"], ict.SCALE_QUANTILE))
    winners = ict.top_marker(scaled)
    assert (winners[beta] == "INS").mean() > 0.9
    assert (winners[~beta] == "ECAD").mean() > 0.99


def test_common_scale_fits_on_eligible_cells_and_keeps_empty_columns_raw():
    values = pd.DataFrame({"A": [1.0, 2.0, 1000.0], "Z": [0.0, 0.0, 0.0]})
    scaled, top = ict.common_scale(values, np.array([True, True, False]), quantile=100.0)
    assert top["A"] == 2.0                                          # the ineligible row does not set the scale
    assert scaled["A"].tolist() == [0.5, 1.0, 500.0]
    assert top["Z"] == 0.0 and scaled["Z"].tolist() == [0.0, 0.0, 0.0]


def test_label_cells_compares_the_level_marker_on_raw_values():
    rng = np.random.default_rng(1)
    n = 3000
    values = pd.DataFrame({"ECAD": rng.lognormal(np.log(7600), 0.5, n), "INS": rng.lognormal(np.log(17), 0.3, n),
                           "CHGA": rng.lognormal(np.log(270), 0.3, n)})
    beta = np.zeros(n, dtype=bool)
    beta[:80] = True
    values.loc[beta, "INS"] = rng.lognormal(np.log(600), 0.2, beta.sum())
    values.loc[beta, "ECAD"] = 2700.0
    # the level marker takes part in the highest-marker test (as CD99 does in Senior): keep it below the hormone
    values.loc[beta, "CHGA"] = rng.lognormal(np.log(900), 0.4, beta.sum())
    x = rng.uniform(0, 5000, n)
    y = rng.uniform(0, 5000, n)
    x[beta] = 2500 + rng.uniform(-40, 40, beta.sum())                # one islet
    y[beta] = 2500 + rng.uniform(-40, 40, beta.sum())
    eligible = np.ones(n, dtype=bool)
    lab = fit.label_cells(values, x, y, eligible, np.zeros(n, dtype=bool), set(), 99.0)
    assert not (lab["top"][beta] == "ECAD").any()                   # raw counts would make every one of them ECAD
    assert lab["hormone_top"][beta].mean() > 0.6
    # the brightest-CHGA beta cells have CHGA highest and come in by the level-marker route instead
    assert set(lab["top"][beta & ~lab["hormone_top"]]) <= {"CHGA"}
    assert lab["endocrine"][beta].mean() >= 0.85
    ref = lab["reference"]
    assert lab["level_marker_reference_level"] == pytest.approx(np.percentile(values["CHGA"].to_numpy()[ref], 99.0))
    assert set(lab["highest_marker_scale"]) == {"ECAD", "INS", "CHGA"}


def test_common_scale_floor_raises_a_low_hormone_top_only():
    values = pd.DataFrame({"INS": [10.0, 20.0, 40.0], "ECAD": [100.0, 200.0, 400.0]})
    eligible = np.ones(3, dtype=bool)
    scaled, top = ict.common_scale(values, eligible, quantile=100.0, floor={"INS": 80.0, "ECAD": 50.0})
    assert top == {"INS": 80.0, "ECAD": 400.0}                      # INS raised to the floor; ECAD keeps its own
    assert scaled["INS"].tolist() == [0.125, 0.25, 0.5]


# ---------------------------------------------------------------------------------------------- the runner

import json  # noqa: E402

import cohort_scale  # noqa: E402
import run_uniform as ru  # noqa: E402


def test_dependents_follow_the_step_graph():
    assert ru.ISLET_CHAIN[0] == "1" and "fit" not in ru.ISLET_CHAIN
    assert "C_planes" not in ru.dependents("0")                     # the planes are image pixels, not cells
    assert {"C_cells", "fit", "1", "J_islets"} <= set(ru.dependents("2"))
    assert set(ru.CELL_STEPS) & set(ru.dependents("fit")) == set()


def _project(tmp_path: Path) -> tuple[ru.Runner, Path]:
    out = tmp_path / "out"
    (out / "qupath" / "data" / "7").mkdir(parents=True)
    (out / "qupath" / "project.qpproj").write_text(json.dumps({"images": [{"imageName": "T1D_0001", "entryID": 7}]}))
    qp = out / "qupath" / "data" / "7" / "data.qpdata"
    qp.write_bytes(b"after the cell steps")
    cohort = tmp_path / "cohort.csv"
    cohort.write_text("case_id,donor_status,cohort,image_name,source_project,source_entry_id\n1,T1D,original,T1D_0001,x,1\n")
    return ru.Runner(out, [{"image_name": "T1D_0001"}], cohort), qp


def test_islet_chain_always_starts_from_the_snapshot(tmp_path):
    runner, qp = _project(tmp_path)
    assert runner.start_islets("T1D_0001")                          # first islet run: snapshot taken
    snap = runner.snapshot_of("T1D_0001")
    assert snap.read_bytes() == b"after the cell steps"
    qp.write_bytes(b"after the islet steps")
    assert runner.start_islets("T1D_0001")                          # a rerun restores it
    assert qp.read_bytes() == b"after the cell steps"


def test_islet_steps_without_a_snapshot_are_refused(tmp_path):
    runner, qp = _project(tmp_path)
    rec = runner.record("T1D_0001", "1b")
    rec.parent.mkdir(parents=True)
    rec.write_text(json.dumps({"status": "pass"}))
    assert not runner.start_islets("T1D_0001")
    assert not runner.snapshot_of("T1D_0001").exists()


def test_running_a_step_retires_what_rests_on_it(tmp_path):
    runner, qp = _project(tmp_path)
    assert runner.start_islets("T1D_0001")
    for step in ("1", "4", "J_islets", "C_planes"):
        runner.record("T1D_0001", step).parent.mkdir(parents=True, exist_ok=True)
        runner.record("T1D_0001", step).write_text(json.dumps({"status": "pass"}))
    runner.invalidate("T1D_0001", "fit")
    assert not runner.record("T1D_0001", "1").exists() and not runner.record("T1D_0001", "J_islets").exists()
    assert runner.record("T1D_0001", "C_planes").exists()           # does not rest on the fit
    assert runner.snapshot_of("T1D_0001").exists()                  # the fit does not change the image
    runner.invalidate("T1D_0001", "2")
    assert not runner.snapshot_of("T1D_0001").exists()              # new measurements: the copy is stale
    assert any((runner.snapshot_of("T1D_0001").parent / "history").iterdir())


def test_cohort_scale_must_hold_every_current_cells_file(tmp_path):
    out = tmp_path
    (out / "cells").mkdir()
    for tag in ("A", "B"):
        (out / "cells" / f"{tag}.parquet").write_bytes(tag.encode())
    scale = {"images": {t: cohort_scale.file_id(cohort_scale.cells_path(out, t)) for t in ("A", "B")}}
    assert cohort_scale.holds_current(scale, out, ["A", "B"]) == []
    assert cohort_scale.holds_current(scale, out, ["A", "B", "C"]) == ["C"]
    (out / "cells" / "B.parquet").write_bytes(b"re-exported")
    assert cohort_scale.holds_current(scale, out, ["A", "B"]) == ["B"]


def test_only_insulin_and_somatostatin_are_floored():
    assert ict.FLOORED_HORMONES == ("INS", "SST")                  # alpha cells are in every image: GCG keeps its own top
    assert set(ict.FLOORED_HORMONES) <= set(ict.HORMONES)
    assert cohort_scale.rule() == {"quantile": ict.SCALE_QUANTILE, "floor_markers": ["INS", "SST"]}


def test_a_fit_on_a_scale_built_under_another_rule_is_stale(tmp_path):
    runner, _ = _project(tmp_path)
    cells = cohort_scale.cells_path(runner.out, "T1D_0001")
    cells.parent.mkdir(parents=True)
    cells.write_bytes(b"cells")
    scale_file = runner.cohort_scale_file()
    scale_file.parent.mkdir(parents=True)
    rec = runner.record("T1D_0001", "fit")
    rec.parent.mkdir(parents=True)

    def write(scale: dict) -> None:
        scale_file.write_text(json.dumps(scale))
        rec.write_text(json.dumps({"status": "pass", "cohort_scale": {"sha256": ru.sha256(scale_file)}}))

    images = {"T1D_0001": cohort_scale.file_id(cells)}
    write({**cohort_scale.rule(), "images": images})
    assert runner.passed("T1D_0001", "fit")
    write({"quantile": ict.SCALE_QUANTILE, "floor_markers": list(ict.HORMONES), "images": images})   # the first rule
    assert not runner.passed("T1D_0001", "fit")
    assert not cohort_scale.follows_rule(json.loads(scale_file.read_text()))
