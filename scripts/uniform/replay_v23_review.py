"""Build cohort-wide candidate identities and QuPath review deltas from frozen v23 evidence."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import review_project as review

sys.path.insert(0, str(Path(__file__).resolve().parent / "typing"))
from review_revision import RETIRED, apply_decisions, reconsider

SENTINEL = "87b35e6a-1f95-46bc-b1d8-85afe222ad09"


def build(typing, baseline, dest):
    cohort = pd.read_csv(Path(__file__).with_name("cohort.csv"), dtype=str)
    dest.mkdir(parents=True, exist_ok=False)
    for name in ("identity", "cell_types", "qupath", "audit", "samples"):
        (dest / name).mkdir()
    totals, transitions = [], []
    release = json.loads((baseline / "uniform_release.json").read_text())

    def digest(path):
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()

    if digest(typing / "levels.csv") != release["levels_sha256"]:
        raise ValueError("Frozen donor levels differ from the v23 release")
    levels = pd.read_csv(typing / "levels.csv").set_index(["image", "marker"])
    from registry import load

    reg = load()
    for tag in cohort.image_name:
        receipt = json.loads((typing / "coverage" / (tag + ".json")).read_text())
        if (
            digest(typing / "identity" / (tag + ".parquet"))
            != release["identity_sha256"][tag]
            or digest(typing / "coverage" / (tag + ".parquet"))
            != receipt["parquet_sha256"]
            or receipt["inputs"]["sha256"]["typing/levels.csv"]
            != release["levels_sha256"]
        ):
            raise ValueError(
                f"{tag}: replay inputs differ from frozen v23 measurements"
            )
        identity, coverage, context = [
            pd.read_parquet(typing / name / (tag + ".parquet"))
            for name in ("identity", "coverage", "context")
        ]
        column = reg.column("LGALS3", levels.at[(tag, "LGALS3"), "measurement"])
        raw_profile = pd.read_parquet(
            typing / "redsea" / (tag + ".parquet"), columns=["object_id", column]
        ).set_index("object_id")
        audit = reconsider(
            identity, coverage, context, np.isfinite(raw_profile[column])
        )
        out = apply_decisions(identity, audit)
        case = str(int(tag.split("_")[-1]))
        table = pd.read_parquet(
            baseline / "cell_types" / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        old = identity.loc[identity.label.ne("Excluded")].set_index(
            "object_id", verify_integrity=True
        )
        if set(table.index) != set(old.index) or not table.phenotype.equals(
            old.label.reindex(table.index)
        ):
            raise ValueError(f"{tag}: assembled v23 UUID/label mismatch")
        indexed = out.set_index("object_id", verify_integrity=True)
        for col in (
            "compartment",
            "lineage",
            "named_by",
            "tags",
            "status",
            "phenotype",
        ):
            table[col] = indexed.loc[
                table.index, "label" if col == "phenotype" else col
            ]
        audit["image"] = tag
        audit["destination"] = audit.proposal.where(
            audit.proposal.ne(""), audit.old_label
        )
        if tag == "Aab_6450" and audit.loc[SENTINEL, "proposal"] != "CD8a Tcell":
            raise ValueError("Sentinel failed calibrated CD8a T cell rule")
        if (
            table.phenotype.isin(RETIRED).sum()
            != audit.disposition.eq("review_required").sum()
        ):
            raise ValueError(f"{tag}: incomplete residual reconciliation")
        table.reset_index().to_parquet(
            dest / "cell_types" / (case + ".parquet"), index=False
        )
        out.to_parquet(dest / "identity" / (tag + ".parquet"), index=False)
        audit.to_parquet(dest / "audit" / (tag + ".parquet"))
        delta = indexed.loc[
            audit.index[audit.disposition.eq("retyped")],
            review.IDENTITY_META + ["label"],
        ].copy()
        delta["old_label"] = audit.old_label
        delta["display"] = delta.label.map(review.display)
        delta["confirmed"] = ""
        delta.to_csv(dest / "qupath" / (tag + ".tsv"), sep="\t")
        # Fixed selection before pixel review: one sample for each transition or
        # residual reason on every donor, controls and the user-reviewed UUID.
        sample = audit.sample(frac=1, random_state=23).drop_duplicates(
            ["old_label", "destination", "old_status", "rule"]
        )
        controls = (
            old[~old.label.isin(RETIRED)]
            .sample(frac=1, random_state=23)
            .drop_duplicates("label")
        )
        selected = sample.index.union(controls.index)
        if tag == "Aab_6450":
            selected = selected.union([SENTINEL])
        pd.Series(selected).to_csv(
            dest / "samples" / (tag + ".txt"), index=False, header=False
        )
        totals.append(
            table.groupby("phenotype")
            .size()
            .rename("cells")
            .reset_index()
            .assign(image=tag)
        )
        transitions.append(
            audit.groupby(
                ["old_label", "destination", "disposition", "rule"], dropna=False
            )
            .size()
            .rename("cells")
            .reset_index()
            .assign(image=tag)
        )
        print(tag, len(table), audit.disposition.value_counts().to_dict(), flush=True)
    pd.concat(totals).to_csv(dest / "phenotype_counts.csv", index=False)
    transitions = pd.concat(transitions)
    denominator = transitions.groupby(["image", "old_label"]).cells.transform("sum")
    transitions["percent_of_source_in_scope"] = 100 * transitions.cells / denominator
    by_destination = transitions.groupby(
        ["old_label", "destination", "disposition"], as_index=False
    ).cells.sum()
    by_destination["percent_of_source_in_scope"] = (
        100
        * by_destination.cells
        / by_destination.groupby("old_label").cells.transform("sum")
    )
    by_destination.to_csv(dest / "cohort_transitions.csv", index=False)
    transitions.to_csv(dest / "transitions.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("typing", "baseline", "dest"):
        p.add_argument("--" + name, type=Path, required=True)
    build(**vars(p.parse_args()))
