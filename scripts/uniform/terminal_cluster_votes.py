"""Donor-held-out v26 transfer. Only approved named reference types can receive votes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from recluster_residuals import (
    SEED,
    STAGES,
    STRUCTURES,
    cluster_votes,
    constraints,
    marker_matrix,
)
from sklearn.ensemble import ExtraTreesClassifier
from threadpoolctl import threadpool_limits


def fit(prior, typing, historical, dest):
    from resolve_tissue_calls import TYPES

    ref = pd.read_parquet(
        prior / "heldout_predictions.parquet",
        columns=["image", "object_id", "label", "fold"],
    )
    history = pd.read_parquet(historical, columns=["image", "object_id", *STAGES])
    tables = []
    for image, controls in ref.groupby("image"):
        (f,) = controls.fold.unique()
        d = pd.read_parquet(prior / "identity" / (image + ".parquet")).set_index(
            "object_id", verify_integrity=True
        )
        target = d.label.eq("Tissue") & d.status.isin(["", "endocrine_marker_only"])
        selected = d.loc[d.index.isin(controls.object_id) | target]
        control = controls.set_index("object_id").label
        if not d.label.loc[control.index].equals(control):
            raise ValueError(f"{image}: v25 reference labels differ")
        p = selected[
            [
                "label",
                "n_positive",
                "start_class",
                *selected.filter(regex="^(pos__|fold__)").columns,
            ]
        ].copy()
        cv = pd.read_parquet(typing / "coverage" / (image + ".parquet")).set_index(
            "object_id", verify_integrity=True
        )
        ctx = pd.read_parquet(typing / "context" / (image + ".parquet")).set_index(
            "object_id", verify_integrity=True
        )
        cols = [
            k
            for name, field, _ in STRUCTURES.values()
            for k in ("in_" + name, "crossed_frac_" + field)
        ]
        p = p.join(cv.loc[p.index].filter(regex="^covint__")).join(
            ctx.loc[p.index, cols]
        )
        ids = (
            history.loc[history.image.eq(image)]
            .set_index("object_id", verify_integrity=True)[STAGES]
            .reindex(p.index)
            .fillna(-1)
        )
        p["parent"] = ids.astype(int).astype(str).agg("/".join, axis=1)
        p["fold"], p["image"] = f, image
        p["control"], p["target"] = p.index.isin(control.index), target.loc[p.index]
        tables.append(p.reset_index())
    data = pd.concat(tables, ignore_index=True)
    train = data.control & data.label.isin(TYPES)
    labels = np.array(sorted(data.loc[train, "label"].unique()))
    _, x = marker_matrix(data)
    y = pd.Categorical(data.label, categories=labels).codes
    allowed = constraints(data, labels)
    parent = data.parent.to_numpy()
    parts = []
    for f in sorted(data.fold.unique()):
        tr, te = (
            np.flatnonzero(train & data.fold.ne(f)),
            np.flatnonzero(data.fold.eq(f)),
        )
        local, _ = cluster_votes(x, y, parent, tr, te, len(labels))
        model = ExtraTreesClassifier(
            n_estimators=96,
            min_samples_leaf=3,
            max_features=0.65,
            n_jobs=6,
            random_state=SEED + f,
        )
        model.fit(x[tr], y[tr])
        model.n_jobs = 1  # sum tree predictions in a fixed order
        forest = np.zeros_like(local)
        forest[:, model.classes_] = model.predict_proba(x[te])
        weight = 0.2 * (local.sum(axis=1) > 0)[:, None]
        scores = np.where(allowed[te], (1 - weight) * forest + weight * local, 0)
        top = np.argsort(-scores, axis=1, kind="stable")[:, :2]
        row = np.arange(len(te))
        result = data.iloc[te][
            ["image", "object_id", "label", "target", "control", "fold"]
        ].copy()
        result["vote"] = scores[row, top[:, 0]]
        result["margin"] = result.vote - scores[row, top[:, 1]]
        result["candidate"] = np.where(result.vote.gt(0), labels[top[:, 0]], "")
        parts.append(result)
        print(f"Finished fold {f}", flush=True)
    votes = pd.concat(parts, ignore_index=True)
    votes.to_parquet(dest / "named_votes.parquet", index=False)
    controls = votes.loc[
        votes.control & votes.vote.ge(0.75) & votes.margin.ge(0.20)
    ].copy()
    controls["correct"] = controls.candidate.eq(controls.label)
    controls.groupby(["fold", "candidate"]).correct.agg(["sum", "size"]).to_csv(
        dest / "named_vote_controls.csv"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("prior", "typing", "historical", "dest"):
        p.add_argument("--" + k, type=Path, required=True)
    with threadpool_limits(limits=2):
        fit(**vars(p.parse_args()))
