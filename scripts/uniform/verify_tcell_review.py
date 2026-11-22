"""Verify the v27 release against v26, the native readback, and a fresh pipeline replay."""

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import review_project as review
from resolve_tcell_subtypes import BROAD, CD4, CD8, KEYS, direct_call
from resolve_tissue_calls import TERMINAL


def read(root, folder, name):
    return pd.read_parquet(root / folder / (name + ".parquet")).set_index(
        "object_id", verify_integrity=True
    )


def verify(prior, current, project, replay=None):
    audit = pd.read_parquet(current / "subtype_assignments.parquet")
    audited = audit.set_index(KEYS, verify_integrity=True)
    controls = pd.read_parquet(current / "heldout_subtypes.parquet")
    if not audited.index.intersection(
        controls.set_index(KEYS, verify_integrity=True).index
    ).empty:
        raise ValueError("A target was used as a reference")
    entries = {
        e["imageName"]: str(e["entryID"])
        for e in json.loads((project / "project.qpproj").read_text())["images"]
    }
    catalog = {
        c["name"]
        for c in json.loads((project / "classifiers/classes.json").read_text())[
            "pathClasses"
        ]
    }
    if (
        review.display(BROAD) in catalog
        or not {review.display(CD4), review.display(CD8)} <= catalog
    ):
        raise ValueError("Invalid native subtype catalog")
    if set(entries) != {p.stem for p in (current / "identity").glob("*.parquet")}:
        raise ValueError("Donor inventory mismatch")
    native_counts = Counter()
    analysed = 0
    for tag, entry in entries.items():
        case = str(int(tag.split("_")[-1]))
        old = read(prior, "identity", tag)
        new = read(current, "identity", tag)
        pd.testing.assert_index_equal(old.index, new.index)
        target = old.label.eq(BROAD)
        decisions = audit.loc[audit.image.eq(tag)].set_index(
            "object_id", verify_integrity=True
        )
        ids = decisions.index
        if set(old.index[target]) != set(ids) or set(new.label) - (TERMINAL - {BROAD}):
            raise ValueError(f"{tag}: invalid target or terminal inventory")
        pd.testing.assert_frame_equal(
            old.loc[~target], new.loc[~target], check_exact=True
        )
        mutable = ["label", "named_by", "subtype_status", "tags"]
        pd.testing.assert_frame_equal(
            old.drop(columns=mutable), new.drop(columns=mutable), check_exact=True
        )
        pd.testing.assert_series_equal(
            new.loc[ids, "label"], decisions.destination, check_names=False
        )
        gate = direct_call(old.loc[ids])
        inferred = gate.eq("")
        assert new.label.loc[ids[~inferred]].eq(gate[~inferred]).all()
        assert (
            new.subtype_status.loc[ids]
            .eq("inferred_below_both_levels")
            .equals(inferred)
        )
        before, after = (read(root, "cell_types", case) for root in (prior, current))
        pd.testing.assert_index_equal(before.index, after.index)
        pd.testing.assert_frame_equal(
            before.loc[~before.phenotype.eq(BROAD)],
            after.loc[~before.phenotype.eq(BROAD)],
            check_exact=True,
        )
        fields = ["phenotype", "named_by", "tags"]
        pd.testing.assert_frame_equal(
            before.drop(columns=fields), after.drop(columns=fields), check_exact=True
        )
        assert after.phenotype.eq(new.loc[after.index, "label"]).all()
        assert set(after.index[before.phenotype.eq(BROAD)]) == set(ids)
        rec = json.loads((current / "samples" / (tag + ".json")).read_text())
        prev = json.loads((prior / "samples" / (tag + ".json")).read_text())
        assert (
            str(rec["entry"]) == entry
            and rec["changed"] == len(ids)
            and rec["detections"] == prev["detections"]
        )
        expected = Counter(prev["classes"])
        assert expected.pop(review.display(BROAD), 0) == len(ids)
        expected.update(decisions.destination.map(review.display))
        assert (
            dict(expected) == rec["classes"]
            and sum(expected.values()) == rec["detections"]
        )
        native_counts.update(expected)
        analysed += len(after)
        if replay:
            pd.testing.assert_frame_equal(
                new, read(replay, "identity", tag), check_exact=True
            )
            pd.testing.assert_frame_equal(
                after, read(replay, "cell_types", case), check_exact=True
            )
    if replay:
        fresh = pd.read_parquet(replay / "subtype_assignments.parquet").set_index(
            KEYS, verify_integrity=True
        )
        pd.testing.assert_frame_equal(
            audited.sort_index(), fresh.sort_index(), check_exact=True
        )
    idx = read(current, "identity", "Aab_6450").at[
        "87b35e6a-1f95-46bc-b1d8-85afe222ad09", "label"
    ]
    assert idx == CD8
    summary = {
        "donors": len(entries),
        "native_detections": sum(native_counts.values()),
        "analysed": analysed,
        "retyped": len(audit),
        "broad_T_cells": native_counts[review.display(BROAD)],
        "transitions": pd.crosstab(audit.rule, audit.destination).to_dict(),
        "native_class_counts": dict(native_counts),
        "index_cell": idx,
        "inferred_winning_vote_below_0.60": int(audit.subtype_vote.lt(0.60).sum()),
        "controls": len(controls),
        "heldout_reference_agreement": float(controls.profile_correct.mean()),
        "weak_reference_agreement": float(
            controls.loc[controls.weak_reference, "profile_correct"].mean()
        ),
        "other_markers_only_balanced_accuracy": float(
            controls.groupby("label").panel_only_correct.mean().mean()
        ),
        "full_replay_identical": replay is not None,
    }
    (current / "reconciliation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        {
            k: v
            for k, v in summary.items()
            if k not in ("native_class_counts", "transitions")
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("prior", "current", "project"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--replay", type=Path)
    verify(**vars(p.parse_args()))
