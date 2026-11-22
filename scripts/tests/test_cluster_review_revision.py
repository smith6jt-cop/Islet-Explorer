"""Regression cases for the separate v23 review replay."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "uniform" / "typing"))
import registry
from review_revision import reconsider


def trial(label, lineage, positives, *, runner="", status="", cluster=False, gal3=True):
    row = {
        "object_id": "cell",
        "start_class": "EndocrineCluster" if cluster else "Tissue",
        "label": label,
        "compartment": "Epithelial" if cluster else "Hematopoietic",
        "lineage": lineage,
        "tags": "canary_top_Hematopoietic" if cluster else "",
        "status": status,
        "compartment_basis": "object" if cluster else "canary",
        "runner_up_compartment": runner,
        "spatial_tie_evidence": "",
        "n_interior_pixels": 100,
    }
    cov = {"object_id": "cell", "n_interior_pixels": 100}
    for m in registry.load().typing_markers():
        strength, interior = positives.get(m, (0, 0))
        row["fold__" + m] = strength
        row["pos__" + m] = int(strength > 0)
        row["canary__" + m] = int(strength >= 1)
        cov["covint__" + m] = interior
    ctx = {"object_id": "cell", "x_um": 1.0, "y_um": 2.0}
    for name, field in [
        ("Blood vessel", "vessel"),
        ("Nerve", "nerve"),
        ("Lymphatic", "lymphatic"),
    ]:
        ctx["in_" + name] = False
        ctx["crossed_frac_" + field] = 0.0
    return reconsider(
        *[pd.DataFrame([r]) for r in (row, cov, ctx)], pd.Series({"cell": gal3})
    ).iloc[0]


def test_cluster_recovers_ranked_t_lineage_despite_competing_cd163():
    markers = {
        "CD45": (4, 1),
        "CD3e": (3, 1),
        "CD8a": (4, 1),
        "CD4": (1, 0.6),
        "CD163": (1, 0.1),
        "INS": (0.1, 0.68),
    }
    result = trial("Beta cell", "Endocrine", markers, cluster=True)
    assert (result.proposal, result.named_by) == ("CD8a Tcell", "CD8a")


def test_ins_positive_cluster_without_cd45_keeps_endocrine_call():
    result = trial(
        "Beta cell",
        "Endocrine",
        {"CD3e": (3, 1), "CD8a": (4, 1), "INS": (3, 0.98)},
        cluster=True,
    )
    assert result.disposition == "keep"


def test_cd45_exception_requires_localized_lineage_signal():
    strong = trial("Hematopoietic CD45-", "Myeloid", {"CD68": (2, 0.8)})
    weak = trial(
        "Hematopoietic CD45-",
        "Myeloid",
        {"CD68": (2, 0.4), "SMA": (3, 1)},
        runner="Mesenchymal",
    )
    assert strong.proposal == "Macrophage_LGALS3-"
    assert weak.disposition == "review_required"


def test_isolated_cd45_does_not_acquire_a_t_subtype():
    assert (
        trial("Hematopoietic CD45+", "", {"CD45": (4, 1)}).disposition
        == "review_required"
    )


def test_review_pattern_requires_shared_only_evidence():
    assert (
        trial(
            "Unclassified", "", {"VIM": (0.5, 0.4)}, status="shared_marker_only"
        ).proposal
        == "Fibroblast"
    )
    assert (
        trial("Unclassified", "", {}, status="no_marker_above_background").disposition
        == "review_required"
    )


def test_missing_gal3_keeps_broad_identity():
    sample = trial("Hematopoietic CD45-", "Myeloid", {"CD68": (2, 0.8)}, gal3=False)
    assert sample.proposal == "Macrophage"
