"""Resolve the broad v26 T lineage to CD4/CD8a while retaining the measurement basis."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import review_project as review
from sklearn.ensemble import ExtraTreesClassifier

BROAD = "T cell"
CD4 = "CD4 Tcell"
CD8 = "CD8a Tcell"
LABELS = [BROAD, CD4, CD8]
KEYS = ["image", "object_id"]
CLUSTERS = ["stage1_cluster", "immune_cluster"]
LOCAL = ["membrane_contrast", "interior_coverage_contrast"]


def aligned_measurements(path, ids, columns):
    return (
        pd.read_parquet(path, columns=["object_id", *columns])
        .set_index("object_id", verify_integrity=True)
        .loc[ids, columns]
    )


def load_panel(prior, typing, historical):
    history = pd.read_parquet(historical, columns=KEYS + CLUSTERS)
    levels = pd.read_csv(typing / "levels.csv").set_index(
        ["image", "marker"], verify_integrity=True
    )
    frames = []
    for path in sorted((prior / "identity").glob("*.parquet")):
        tag = path.stem
        cells = pd.read_parquet(path, filters=[("label", "in", LABELS)]).set_index(
            "object_id", verify_integrity=True
        )
        cols = [
            "label",
            "named_by",
            "pos__CD3e",
            "pos__CD4",
            "pos__CD8a",
            *cells.filter(regex="^fold__").columns,
        ]
        frame = cells[cols].copy()

        rim = aligned_measurements(
            typing / "redsea" / path.name,
            frame.index,
            ["Membrane__CD8a", "Membrane__CD4"],
        )
        cov = aligned_measurements(
            typing / "coverage" / path.name,
            frame.index,
            ["n_interior_pixels", "covint__CD8a", "covint__CD4"],
        )
        if not cov.n_interior_pixels.eq(cells.n_interior_pixels).all():
            raise ValueError(f"{tag}: interior count mismatch")
        # Express the two membrane readings relative to this image's existing marker levels.
        # The difference is a continuous predictor; neither whole-cell gate is changed.
        shift = (
            levels.at[(tag, "CD8a"), "level_log"] - levels.at[(tag, "CD4"), "level_log"]
        )
        frame["membrane_contrast"] = (
            np.log1p(rim.Membrane__CD8a.clip(lower=0))
            - np.log1p(rim.Membrane__CD4.clip(lower=0))
            - shift
        )
        frame["interior_coverage_contrast"] = (
            cov.covint__CD8a - cov.covint__CD4
        ).where(cov.n_interior_pixels.gt(0))
        frame["image"] = tag
        frames.append(frame.reset_index())
    panel = pd.concat(frames, ignore_index=True).merge(
        history, how="left", on=KEYS, validate="one_to_one"
    )
    if (
        panel[["pos__CD4", "pos__CD8a", "pos__CD3e"]].isna().any().any()
        or not np.isfinite(panel[["fold__CD4", "fold__CD8a"]]).all().all()
    ):
        raise ValueError("Missing direct subtype measurement")
    return panel


def direct_call(panel):
    cd4, cd8 = panel.pos__CD4.eq(1), panel.pos__CD8a.eq(1)
    prefer8 = cd8 & (~cd4 | panel.fold__CD8a.ge(panel.fold__CD4))
    return pd.Series(
        np.select([prefer8, cd4], [CD8, CD4], default=""), index=panel.index
    )


def reference_votes(panel):
    is4 = (
        panel.label.eq(CD4)
        & panel.named_by.eq("CD4")
        & panel.pos__CD4.eq(1)
        & panel.pos__CD8a.eq(0)
    )
    is8 = (
        panel.label.eq(CD8)
        & panel.named_by.eq("CD8a")
        & panel.pos__CD8a.eq(1)
        & panel.pos__CD4.eq(0)
    )
    anchor = panel.pos__CD3e.eq(1) & (is4 | is8)
    target = panel.label.eq(BROAD)
    remaining_markers = [
        c for c in panel.filter(regex="^fold__") if c not in ("fold__CD4", "fold__CD8a")
    ]
    out = []
    for image in sorted(panel.image.unique()):
        train = panel.loc[anchor & panel.image.ne(image)]
        test = panel.loc[(anchor | target) & panel.image.eq(image)]
        if test.empty:
            continue
        if set(train.label) != {CD4, CD8}:
            raise ValueError(f"{image}: missing independent subtype references")
        # Give the two reference classes equal training weight; within each class balance donors.
        sizes = train.groupby(["label", "image"]).label.transform("size")
        donors = train.groupby("label").image.transform("nunique")
        weight = 1 / (sizes * donors)
        evidence = test[
            KEYS
            + [
                "label",
                "named_by",
                "pos__CD4",
                "pos__CD8a",
                "fold__CD4",
                "fold__CD8a",
                *LOCAL,
                *CLUSTERS,
            ]
        ].copy()
        for features, column in [
            (remaining_markers, "panel_only_cd8_vote"),
            (remaining_markers + LOCAL, "profile_cd8_vote"),
        ]:
            model = ExtraTreesClassifier(
                n_estimators=96,
                min_samples_leaf=5,
                max_features=0.65,
                random_state=27,
                n_jobs=4,
            )
            model.fit(
                train[features],
                train.label.eq(CD8).astype(int),
                sample_weight=weight.to_numpy() * len(train) / weight.sum(),
            )
            model.n_jobs = 1
            evidence[column] = model.predict_proba(test[features])[:, 1]
        # Transfer votes only from other donors in the same historical cluster. Unmatched cells use the profile alone.
        counts = train.groupby(CLUSTERS).label.agg(
            n="size", cd8=lambda s: s.eq(CD8).sum()
        )
        evidence = evidence.join(counts, on=CLUSTERS, validate="many_to_one")
        evidence["cluster_n"] = evidence.pop("n").fillna(0).astype(int)
        evidence["cluster_cd8_vote"] = evidence.pop("cd8").div(
            evidence.cluster_n.where(evidence.cluster_n.gt(0))
        )
        w = 0.20 * (evidence.cluster_n / 20).clip(upper=1)
        evidence["cd8_vote"] = (
            1 - w
        ) * evidence.profile_cd8_vote + w * evidence.cluster_cd8_vote.fillna(0)
        evidence["inferred_label"] = np.where(
            evidence.cd8_vote.gt(0.5)
            | (evidence.cd8_vote.eq(0.5) & evidence.fold__CD8a.ge(evidence.fold__CD4)),
            CD8,
            CD4,
        )
        out.append(evidence)
    votes = pd.concat(out, ignore_index=True)
    if not np.isfinite(votes.cd8_vote).all():
        raise ValueError("Nonfinite subtype vote")
    return votes


def assign(votes):
    targets = votes.loc[votes.label.eq(BROAD)].copy()
    direct = direct_call(targets)
    targets["destination"] = direct.where(direct.ne(""), targets.inferred_label)
    targets["rule"] = np.select(
        [direct.eq(""), targets.pos__CD4.eq(1) & targets.pos__CD8a.eq(1)],
        ["inferred_below_both_levels", "higher_of_two_positive"],
        default="single_positive",
    )
    targets["subtype_vote"] = np.maximum(targets.cd8_vote, 1 - targets.cd8_vote).where(
        direct.eq("")
    )
    if set(targets.destination) - {CD4, CD8} or targets.duplicated(KEYS).any():
        raise ValueError("Invalid T-cell resolution")
    return targets


def validation(votes):
    controls = votes.loc[votes.label.ne(BROAD)].copy()
    win = np.where(controls.label.eq(CD8), controls.fold__CD8a, controls.fold__CD4)
    controls["gate_distance"] = win
    weak = controls.groupby(["image", "label"]).gate_distance.transform(
        lambda v: v.le(v.quantile(0.25))
    )
    controls["weak_reference"] = weak
    controls["profile_correct"] = controls.inferred_label.eq(controls.label)
    controls["panel_only_correct"] = controls.panel_only_cd8_vote.ge(0.5).eq(
        controls.label.eq(CD8)
    )
    return controls


def build(prior, typing, historical, dest):
    policy = json.loads((prior / "terminal_classes.json").read_text())
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("identity", "cell_types", "qupath", "samples"):
        (dest / name).mkdir(exist_ok=False)
    panel = load_panel(prior, typing, historical)
    votes = reference_votes(panel)
    targets = assign(votes)
    controls = validation(votes)
    targets.to_parquet(dest / "subtype_assignments.parquet", index=False)
    controls.to_parquet(dest / "heldout_subtypes.parquet", index=False)
    counts = []
    for path in sorted((prior / "identity").glob("*.parquet")):
        tag = path.stem
        case = str(int(tag.split("_")[-1]))
        d = pd.read_parquet(path).set_index("object_id", verify_integrity=True)
        tab = pd.read_parquet(prior / "cell_types" / (case + ".parquet")).set_index(
            "object_id", verify_integrity=True
        )
        rows = targets.loc[targets.image.eq(tag)].set_index(
            "object_id", verify_integrity=True
        )
        ids = rows.index
        if set(ids) != set(d.index[d.label.eq(BROAD)]) or set(ids) != set(
            tab.index[tab.phenotype.eq(BROAD)]
        ):
            raise ValueError(f"{tag}: broad T-cell inventories differ")
        inferred = rows.rule.eq("inferred_below_both_levels")
        named = rows.destination.map({CD4: "CD4", CD8: "CD8a"}).mask(
            inferred, "T_subtype_reference"
        )
        d.loc[ids, "label"] = rows.destination
        d.loc[ids, "named_by"] = named
        d.loc[ids, "subtype_status"] = rows.rule
        d.loc[ids, "tags"] = (d.loc[ids, "tags"] + ";t_subtype_" + rows.rule).str.strip(
            ";"
        )
        for col, source in [
            ("phenotype", "label"),
            ("tags", "tags"),
            ("named_by", "named_by"),
        ]:
            tab.loc[ids, col] = d.loc[ids, source]
        if d.label.eq(BROAD).any() or tab.phenotype.eq(BROAD).any():
            raise ValueError(f"{tag}: broad T call survived")
        d.reset_index().to_parquet(dest / "identity" / path.name, index=False)
        tab.reset_index().to_parquet(
            dest / "cell_types" / (case + ".parquet"), index=False
        )
        delta = d.loc[
            ids, ["label", "tags", "named_by", "subtype_status", "excluded_reason"]
        ].copy()
        delta["old_label"], delta["old_display"] = BROAD, review.display(BROAD)
        delta["display"] = rows.destination.map(review.display)
        delta["t_subtype_cd8_vote"] = rows.cd8_vote.where(inferred).round(6)
        delta["t_subtype_cluster_references"] = rows.cluster_n.where(inferred)
        delta.to_csv(dest / "qupath" / (tag + ".tsv"), sep="\t")
        (dest / "samples" / (tag + ".txt")).write_text("")
        counts.append(
            rows.groupby(["destination", "rule"])
            .size()
            .rename("cells")
            .reset_index()
            .assign(image=tag, old_label=BROAD, disposition="retyped")
        )
    pd.concat(counts).to_csv(dest / "transitions.csv", index=False)
    # Keep the accumulated review decisions, including all exclusions required by the assembler.
    previous = pd.read_parquet(prior / "assignments.parquet").set_index(
        KEYS, verify_integrity=True
    )
    current = targets.set_index(KEYS, verify_integrity=True)[
        ["destination", "rule"]
    ].copy()
    current["destination"] = current.destination.map(review.display)
    current["old_label"] = review.display(BROAD)
    shared = current.index.intersection(previous.index)
    if not previous.loc[shared, "destination"].eq(review.display(BROAD)).all():
        raise ValueError("Prior decision is not a broad T cell")
    current.loc[shared, "old_label"] = previous.loc[shared, "old_label"]
    pd.concat([previous.drop(index=shared), current]).reset_index().to_parquet(
        dest / "assignments.parquet", index=False
    )
    policy["detection_classes"].remove(review.display(BROAD))
    policy["retired_detection_classes"] = [review.display(BROAD)]
    policy["retired_phenotypes"] = [BROAD]
    (dest / "terminal_classes.json").write_text(json.dumps(policy, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("prior", "typing", "historical", "dest"):
        parser.add_argument("--" + name, type=Path, required=True)
    build(**vars(parser.parse_args()))
