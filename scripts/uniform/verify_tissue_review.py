"""Reconcile the terminal release with its native readbacks, prior cohort and cell-level evidence."""

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import review_project as review
from resolve_tissue_calls import CONFLICT, REASON, TERMINAL, UNDETECTED


def verify(prior, baseline, current, project, replay=None):
    decisions = pd.read_parquet(current / "assignments.parquet")
    release = json.loads((baseline / "uniform_release.json").read_text())
    names = {
        r["imageName"]: str(r["entryID"])
        for r in json.loads((project / "project.qpproj").read_text())["images"]
    }
    if set(names) != set(release["identity_sha256"]):
        raise ValueError("Native donor inventory differs")
    catalog = {
        c["name"]
        for c in json.loads((project / "classifiers/classes.json").read_text())[
            "pathClasses"
        ]
    }
    if not {CONFLICT, UNDETECTED, "Tissue"} <= catalog or any(
        x.startswith("Tissue cell - ") for x in catalog
    ):
        raise ValueError("Invalid native class catalog")
    native_counts, by_image = Counter(), []
    for tag, entry in names.items():
        case = str(int(tag.split("_")[-1]))

        def read(path):
            return pd.read_parquet(path).set_index("object_id", verify_integrity=True)

        old, new = (
            read(root / "identity" / f"{tag}.parquet") for root in (prior, current)
        )
        pd.testing.assert_index_equal(old.index, new.index)
        if set(new.label) - TERMINAL:
            raise ValueError(f"{tag}: unexpected terminal class")
        target = old.label.eq("Tissue")
        pd.testing.assert_frame_equal(old.loc[~target], new.loc[~target])
        changed_cols = [
            "label",
            "excluded_reason",
            "compartment",
            "compartment_basis",
            "compartment_status",
            "lineage",
            "lineage_basis",
            "lineage_status",
            "subtype_status",
            "status",
            "tags",
            "named_by",
            "confirmed",
            "confirming_markers",
            "established_label",
        ]
        pd.testing.assert_frame_equal(
            old.drop(columns=changed_cols), new.drop(columns=changed_cols)
        )
        audited = decisions.loc[decisions.image.eq(tag)].set_index(
            "object_id", verify_integrity=True
        )
        assert set(old.index[target]) == set(audited.index)
        old_status = old.loc[audited.index, "status"]
        for status, cls in (
            ("structure_marker_only", "Stromal"),
            ("shared_marker_only", CONFLICT),
            ("no_marker_above_background", UNDETECTED),
        ):
            assert audited.loc[old_status.eq(status), "destination"].eq(cls).all()
        display = new.label.map(review.display).mask(
            new.excluded_reason.eq(REASON), CONFLICT
        )
        assert display.loc[audited.index].eq(audited.destination).all()
        assert (new.label.eq(UNDETECTED) == new.n_positive.eq(0)).loc[target].all()
        assert not (
            new.start_class.eq("Tissue")
            & new.label.isin(
                ["Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative"]
            )
        ).any()
        receipt = json.loads((current / "samples" / (tag + ".json")).read_text())
        assert str(receipt["entry"]) == entry and receipt["changed"] == len(audited)
        assert receipt["detections"] == sum(receipt["classes"].values())
        counts = Counter(receipt["classes"])
        qc = counts.pop("Excluded - cell QC", 0)
        assert counts == Counter(display)
        native_counts.update(receipt["classes"])
        cell = read(current / "cell_types" / (case + ".parquet"))
        old_cell = read(baseline / "cell_types" / (case + ".parquet"))
        keep = old_cell.index[~new.label.loc[old_cell.index].eq("Excluded")]
        pd.testing.assert_index_equal(cell.index, keep)
        assert cell.phenotype.eq(new.label.loc[cell.index]).all()
        assert not cell.phenotype.isin(["Excluded", "Tissue", CONFLICT]).any()
        static = cell.columns.difference(
            ["phenotype", "compartment", "lineage", "named_by", "tags", "status"]
        )
        pd.testing.assert_frame_equal(cell[static], old_cell.loc[keep, static])
        if replay is not None:
            pd.testing.assert_frame_equal(
                new, read(replay / "identity" / f"{tag}.parquet")
            )
            pd.testing.assert_frame_equal(
                cell, read(replay / "cell_types" / (case + ".parquet"))
            )
        by_image.append(
            {
                "image": tag,
                "detections": receipt["detections"],
                "qc": qc,
                "reassigned": len(audited),
                "retained": len(cell),
            }
        )
    idx = read(current / "identity/Aab_6450.parquet").loc[
        "87b35e6a-1f95-46bc-b1d8-85afe222ad09"
    ]
    assert idx.label == "CD8a Tcell"
    result = {
        "donors": len(names),
        "native_detections": sum(native_counts.values()),
        "native_class_counts": dict(native_counts),
        "analysed": sum(x["retained"] for x in by_image),
        "reassigned": len(decisions),
        "index_cell": idx.label,
        "full_replay_identical": replay is not None,
        "by_image": by_image,
    }
    (current / "reconciliation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        {
            k: v
            for k, v in result.items()
            if k not in ("native_class_counts", "by_image")
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("prior", "baseline", "current", "project"):
        p.add_argument("--" + k, type=Path, required=True)
    p.add_argument("--replay", type=Path)
    verify(**vars(p.parse_args()))
