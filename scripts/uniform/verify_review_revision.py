"""Reconcile a completed review replay, native QuPath readbacks and returned assembly."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def verify(baseline, assembled, evidence, typing):
    release = json.loads((assembled / "uniform_release.json").read_text())
    original = json.loads((baseline / "uniform_release.json").read_text())
    if release["levels_sha256"] != original["levels_sha256"]:
        raise ValueError("Donor calibration changed in assembly")
    changes = pd.read_csv(evidence / "transitions.csv")
    by_case = []
    counts = {}
    for tag, digest in release["identity_sha256"].items():
        case = str(int(tag.split("_")[-1]))
        with (evidence / "identity" / (tag + ".parquet")).open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != digest:
                raise ValueError(
                    f"{tag}: assembled identity hash differs from reviewed file"
                )
        objects = pd.read_parquet(
            assembled / "cell_types" / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        expected = pd.read_parquet(
            evidence / "cell_types" / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        prior = pd.read_parquet(
            baseline / "cell_types" / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        overlap = [
            "crossed_frac_nerve",
            "crossed_frac_lymphatic",
            "crossed_frac_vessel",
        ]
        geometry = pd.read_parquet(
            typing / "context" / (tag + ".parquet"), columns=["object_id", *overlap]
        ).set_index("object_id")
        pd.testing.assert_frame_equal(
            objects[overlap], geometry.loc[objects.index, overlap], check_dtype=False
        )
        pd.testing.assert_frame_equal(
            objects.drop(columns=[c for c in overlap if c not in expected]),
            expected,
            check_dtype=False,
        )
        static = prior.columns.difference(
            ["phenotype", "compartment", "lineage", "status", "named_by", "tags"]
        )
        pd.testing.assert_frame_equal(objects[static], prior[static], check_dtype=False)
        n = int(objects.phenotype.ne(prior.phenotype).sum())
        audit = int(
            changes.loc[
                changes.image.eq(tag) & changes.disposition.eq("retyped"), "cells"
            ].sum()
        )
        native = json.loads((evidence / "samples" / (tag + ".json")).read_text())
        if n != native["changed"] or n != audit:
            raise ValueError(f"{tag}: label changes disagree among sources")
        counts.update(
            {
                (int(case), label): int(size)
                for label, size in objects.phenotype.value_counts().items()
            }
        )
        by_case.append(
            {
                "image": tag,
                "eligible": len(objects),
                "changed": n,
                "qupath_detections": native["detections"],
            }
        )
    distance = pd.read_csv(
        assembled / "cell_distance.csv", usecols=["case_id", "phenotype"]
    )
    if distance.value_counts().to_dict() != counts:
        raise ValueError("Distance analysis phenotype counts differ from UUID tables")
    if (
        len(by_case) != release["donors"]
        or len(distance) != release["cells_analysed"]
        or len(distance) != original["cells_analysed"]
    ):
        raise ValueError("Cohort census differs from release")
    result = {
        "donors": len(by_case),
        "analysed_cells": len(distance),
        "retyped": sum(r["changed"] for r in by_case),
        "native_detections": sum(r["qupath_detections"] for r in by_case),
        "by_donor": by_case,
    }
    (evidence / "reconciliation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "assembled", "evidence", "typing"):
        p.add_argument("--" + name, type=Path, required=True)
    print(verify(**vars(p.parse_args())))
