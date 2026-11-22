"""Resolve v24 residuals with donor-held-out reference transfer within Barlow marker clusters."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import review_project as review
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import ExtraTreesClassifier
from threadpoolctl import threadpool_limits

RETIRED = {"Hematopoietic CD45+", "Hematopoietic CD45-", "Unclassified"}
STRUCTURES = {
    "Endothelium Blood Vessel": ("Blood vessel", "vessel", "CD31"),
    "Endothelium Lymphatic": ("Lymphatic", "lymphatic", "PDPN"),
    "Neural": ("Nerve", "nerve", "B3TUBB"),
}
HORMONES = {
    "Beta cell": "INS",
    "Delta cell": "SST",
    "Alpha cell": "GCG",
    "CHGA+ hormone-negative": "CHGA",
}
STAGES = ["stage1_cluster", "immune_cluster", "endocrine_cluster"]
SEED = 25


def load_sample(typing, prior, historical):
    hist = pd.read_parquet(historical, columns=["image", "object_id", *STAGES])
    frames = []
    for path in sorted((prior / "identity").glob("*.parquet")):
        tag = path.stem
        d = pd.read_parquet(path).set_index("object_id", verify_integrity=True)
        target = d.label.isin(RETIRED)
        # The retired labels contribute no training examples. Sample reference
        # types within donors so rare types enter the fit and validation.
        ref = (
            d[~target & d.label.ne("Excluded")]
            .sample(frac=1, random_state=SEED)
            .groupby("label")
            .head(500)
        )
        selected = pd.concat([ref, d[target]])
        markers = [c[6:] for c in d if c.startswith("fold__")]
        cols = [
            "label",
            "start_class",
            "n_positive",
            *[p + m for p in ("fold__", "pos__") for m in markers],
        ]
        panel = selected[cols].copy()
        panel["image"] = tag
        cov = pd.read_parquet(
            typing / "coverage" / path.name,
            columns=["object_id", *["covint__" + m for m in markers]],
        ).set_index("object_id", verify_integrity=True)
        ctx_cols = [
            k
            for name, field, _ in STRUCTURES.values()
            for k in ("in_" + name, "crossed_frac_" + field)
        ]
        ctx = pd.read_parquet(
            typing / "context" / path.name, columns=["object_id", *ctx_cols]
        ).set_index("object_id", verify_integrity=True)
        panel = panel.join(cov.loc[panel.index]).join(ctx.loc[panel.index])
        old = hist.loc[hist.image.eq(tag)].set_index("object_id", verify_integrity=True)
        panel["matched"] = panel.index.isin(old.index)
        panel = panel.join(old[STAGES]).fillna({s: -1 for s in STAGES})
        panel["parent"] = panel[STAGES].astype(int).astype(str).agg("/".join, axis=1)
        frames.append(panel.reset_index())
    return pd.concat(frames, ignore_index=True)


def cluster_votes(x, y, parents, train, test, nlabels):
    votes = np.zeros((len(test), nlabels))
    child = np.full(len(test), -1)
    # Subdivide each historical population (and the unmatched population)
    # on the current measurements. Labels are used only after fitting.
    for parent in np.unique(parents[test]):
        tr = train[parents[train] == parent]
        loc = np.flatnonzero(parents[test] == parent)
        if len(tr) < 20:
            continue
        n = min(48, max(1, len(tr) // 70))
        km = MiniBatchKMeans(n_clusters=n, random_state=SEED, n_init=3, batch_size=1024)
        codes = km.fit_predict(x[tr])
        table = np.zeros((n, nlabels))
        np.add.at(table, (codes, y[tr]), 1)
        child[loc] = km.predict(x[test[loc]])
        counts = table[child[loc]]
        votes[loc] = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1)
    return votes, child


def constraints(panel, labels):
    mask = np.ones((len(panel), len(labels)), bool)
    for j, label in enumerate(labels):
        # Cluster similarity can resolve a missing pan marker, but a compatible
        # cell-local phenotype marker must support a specific destination.
        if label == "Acinar cell":
            mask[:, j] = panel[["pos__ECAD", "pos__panCK"]].any(
                axis=1
            ) & panel.pos__CK19.eq(0)
        if label == "Ductal cell":
            mask[:, j] = panel.pos__CK19.eq(1)
        if label == "Fibroblast":
            mask[:, j] = panel[["pos__VIM", "pos__CD34", "pos__PDPN"]].any(axis=1)
        if label == "Stromal":
            mask[:, j] = panel.pos__SMA.eq(1)
        if label == "B cell":
            mask[:, j] = panel.pos__CD20.eq(1)
        if label == "T cell":
            mask[:, j] = panel[["pos__CD3e", "pos__CD8a", "pos__CD4"]].any(axis=1)
        if label in ("CD8a Tcell", "CD4 Tcell"):
            marker, other = (
                ("CD8a", "CD4") if label == "CD8a Tcell" else ("CD4", "CD8a")
            )
            preferred = (
                panel["fold__" + marker].ge(panel["fold__" + other])
                if marker == "CD8a"
                else panel["fold__" + marker].gt(panel["fold__" + other])
            )
            mask[:, j] = panel["pos__" + marker].eq(1) & (
                panel["pos__" + other].eq(0) | preferred
            )
        if label in STRUCTURES:
            name, field, m = STRUCTURES[label]
            mask[:, j] = (
                panel["in_" + name].eq(True)
                & panel["crossed_frac_" + field].gt(0.5)
                & panel["pos__" + m].eq(1)
            )
        if label in HORMONES:
            mask[:, j] = panel.start_class.isin(
                ["Islet", "EndocrineCluster", "EndocrineSingle"]
            ) & panel["pos__" + HORMONES[label]].eq(1)
            if label == "CHGA+ hormone-negative":
                mask[:, j] &= ~panel[["pos__INS", "pos__GCG", "pos__SST"]].any(axis=1)
        if label.startswith("Macrophage_LGALS3"):
            mask[:, j] = panel.pos__LGALS3.eq(int(label.endswith("+"))) & panel[
                ["pos__CD68", "pos__CD163"]
            ].any(axis=1)
        if label != "Tissue":
            mask[:, j] &= panel.n_positive.gt(0)
    return mask


def marker_matrix(panel):
    cols = [c for c in panel if c.startswith(("fold__", "covint__", "pos__"))]
    x = panel[cols].fillna(-1).to_numpy(np.float32)
    # Scale continuous fold measurements for Euclidean clustering; bounded
    # coverage and positivity already share a scale across donors.
    for j, c in enumerate(cols):
        if c.startswith("fold__"):
            x[:, j] = np.arcsinh(x[:, j]) / 2
    return cols, x


def group_residuals(panel, destinations, dest):
    """Fit subclusters to residual cells themselves. Characterize each group
    with independent reference votes and the distribution of issued calls."""
    _, x = marker_matrix(panel)
    target = panel.label.isin(RETIRED).to_numpy()
    parent = panel.parent.to_numpy()
    groups = np.full(len(panel), "", dtype=object)
    for group in np.unique(parent[target]):
        ids = np.flatnonzero(parent == group)
        residuals = ids[target[ids]]
        n = min(48, max(1, len(residuals) // 100))
        km = MiniBatchKMeans(n_clusters=n, n_init=3, batch_size=1024, random_state=SEED)
        km.fit(x[residuals])
        groups[ids] = [f"{group}:{c}" for c in km.predict(x[ids])]
    frame = panel[["object_id", "image", "label"]].assign(residual_cluster=groups)
    refs = frame.loc[~target & (groups != "")]
    members = frame.loc[target].merge(
        destinations[["object_id", "image", "destination", "score"]],
        on=["object_id", "image"],
        validate="one_to_one",
    )
    profile = []
    members["compatible_reference_type"] = ""
    members["compatible_reference_share"] = 0.0
    for tag in members.image.unique():
        counts = (
            refs.loc[refs.image.ne(tag)]
            .groupby(["residual_cluster", "label"])
            .size()
            .unstack(fill_value=0)
        )
        which = members.image.eq(tag)
        subset = members.loc[which]
        raw_counts = counts.reindex(subset.residual_cluster, fill_value=0).to_numpy()
        measurements = (
            panel.loc[panel.image.eq(tag)].set_index("object_id").loc[subset.object_id]
        )
        allowed = constraints(measurements, counts.columns)
        compatible = np.where(allowed, raw_counts, 0)
        best = compatible.max(axis=1)
        members.loc[which, "compatible_reference_type"] = np.where(
            best > 0, counts.columns.to_numpy()[compatible.argmax(axis=1)], ""
        )
        members.loc[which, "compatible_reference_share"] = best / np.maximum(
            1, raw_counts.sum(axis=1)
        )
        profile.append(
            pd.DataFrame(
                {
                    "reference_type": counts.idxmax(axis=1),
                    "reference_share": counts.max(axis=1) / counts.sum(axis=1),
                    "reference_cells": counts.sum(axis=1),
                    "image": tag,
                }
            ).reset_index()
        )
    members = members.merge(
        pd.concat(profile),
        on=["image", "residual_cluster"],
        how="left",
        validate="many_to_one",
    )
    members.to_parquet(dest / "residual_cluster_members.parquet", index=False)
    counts = (
        members.groupby(["residual_cluster", "destination"])
        .size()
        .unstack(fill_value=0)
    )
    overview = members.groupby("residual_cluster").agg(
        cells=("object_id", "size"),
        donors=("image", "nunique"),
        median_cell_score=("score", "median"),
    )
    overview["leading_v25_type"] = counts.idxmax(axis=1)
    overview["leading_type_fraction"] = counts.max(axis=1) / overview.cells
    overview.join(counts.add_prefix("count__")).sort_values(
        "cells", ascending=False
    ).to_csv(dest / "residual_clusters.csv")


def infer(panel):
    residual = panel.label.isin(RETIRED).to_numpy()
    labels = np.array(sorted(panel.loc[~residual, "label"].unique()))
    y = pd.Categorical(panel.label, categories=labels).codes
    folds = [
        g
        for _, g in panel[["image"]]
        .drop_duplicates()
        .groupby(lambda i: panel.at[i, "image"].split("_")[0])
    ]
    donor_fold = {
        tag: i % 5 for group in folds for i, tag in enumerate(sorted(group.image))
    }
    fold = panel.image.map(donor_fold).to_numpy()
    cols, x = marker_matrix(panel)
    parents = panel.parent.to_numpy()
    allowed = constraints(panel, labels)
    result = []
    stress = []
    for f in range(5):
        tr = np.flatnonzero((fold != f) & ~residual)
        te = np.flatnonzero(fold == f)
        local, child = cluster_votes(x, y, parents, tr, te, len(labels))
        forest = ExtraTreesClassifier(
            n_estimators=96,
            min_samples_leaf=3,
            max_features=0.65,
            n_jobs=6,
            random_state=SEED + f,
        )
        forest.fit(x[tr], y[tr])
        model = np.zeros_like(local)
        model[:, forest.classes_] = forest.predict_proba(x[te])
        # A lineage panel may be genuinely absent in the residuals. Measure
        # recovery separately after withholding the pan marker from known immune
        # controls. Omit historical cluster votes in this stress test because
        # their memberships were fitted before masking.
        immune_controls = (~residual[te]) & panel.label.iloc[te].isin(
            [
                "CD8a Tcell",
                "CD4 Tcell",
                "T cell",
                "B cell",
                "Macrophage_LGALS3+",
                "Macrophage_LGALS3-",
            ]
        ).to_numpy()
        idx = te[immune_controls]
        for missing in [("CD45",), ("CD45", "CD3e")]:
            xm = x[idx].copy()
            for marker in missing:
                for prefix in ("fold__", "pos__", "covint__"):
                    xm[:, cols.index(prefix + marker)] = -1 if prefix == "fold__" else 0
            prob = np.zeros((len(idx), len(labels)))
            prob[:, forest.classes_] = forest.predict_proba(xm)
            masked_panel = panel.iloc[idx].copy()
            for marker in missing:
                masked_panel["pos__" + marker] = 0
            masked_scores = np.where(constraints(masked_panel, labels), prob, -np.inf)
            pred = labels[masked_scores.argmax(axis=1)]
            pred[masked_scores.max(axis=1) <= 0] = "Tissue"
            stress.append(
                pd.DataFrame(
                    {
                        "image": panel.image.iloc[idx],
                        "reference": panel.label.iloc[idx],
                        "masked": "+".join(missing),
                        "prediction": pred,
                    }
                )
            )
        # If a historical stratum has no training data, use the full forest
        # score. Preserve zero support when geometry rejects all voted classes.
        weight = 0.2 * (local.sum(axis=1) > 0)[:, None]
        scores = np.where(allowed[te], (1 - weight) * model + weight * local, -np.inf)
        empty = scores.max(axis=1) <= 0
        scores[empty] = -np.inf
        scores[empty, np.flatnonzero(labels == "Tissue")[0]] = 0
        top = np.argsort(-scores, axis=1, kind="stable")[:, :2]

        rows = np.arange(len(te))
        res = panel.iloc[te][
            ["image", "object_id", "label", "parent", "matched"]
        ].copy()
        res["fold"] = f
        res["child"] = child
        res["cluster_label"] = np.where(
            local.sum(axis=1) > 0, labels[local.argmax(axis=1)], ""
        )
        res["cluster_vote"] = local.max(axis=1)
        res["destination"] = labels[top[:, 0]]
        res["runner_up"] = np.where(scores[rows, top[:, 1]] > 0, labels[top[:, 1]], "")
        res["no_allowed_vote"] = empty
        res["score"] = scores[rows, top[:, 0]]
        res["margin"] = res.score - np.maximum(0, scores[rows, top[:, 1]])
        res["no_positive_marker"] = panel.iloc[te].n_positive.eq(0).to_numpy()
        res["support_band"] = pd.cut(
            res.score.round(12),  # parallel tree sums can straddle an exact reporting boundary by one ULP
            [-1, 0.5, 0.75, 0.9, 1.001],
            right=False,
            labels=["<0.50", "0.50–0.75", "0.75–0.90", ">=0.90"],
        )
        result.append(res)
        print("Finished fold", f, flush=True)
    return pd.concat(result, ignore_index=True), pd.concat(stress, ignore_index=True)


def write_transitions(targets, dest):
    for keys, name in [
        (["image", "label", "destination"], "transitions.csv"),
        (["label", "destination"], "cohort_transitions.csv"),
    ]:
        tally = targets.groupby(keys).size().rename("cells").reset_index()
        tally["percent_of_source"] = (
            100 * tally.cells / tally.groupby(keys[:-1]).cells.transform("sum")
        )
        tally.rename(columns={"label": "old_label"}).assign(
            disposition="retyped"
        ).to_csv(dest / name, index=False)
    targets.groupby(
        ["label", "destination", "support_band"], observed=True
    ).size().rename("cells").to_csv(dest / "support.csv")
    targets.groupby(
        ["fold", "parent", "child", "destination"], observed=True
    ).size().rename("cells").to_csv(dest / "cluster_composition.csv")


def export(results, prior, baseline, dest):
    targets = results[results.label.isin(RETIRED)]
    if len(targets) != 134001 or targets.duplicated(["image", "object_id"]).any():
        raise ValueError("v24 residual inventory mismatch")
    refs = results[~results.label.isin(RETIRED)].copy()
    refs["correct"] = refs.destination.eq(refs.label)
    refs["cluster_correct"] = refs.cluster_label.eq(refs.label)
    refs.to_parquet(dest / "heldout_predictions.parquet", index=False)
    refs.groupby(["label", "destination"]).size().rename("cells").to_csv(
        dest / "heldout_confusion.csv"
    )
    refs.groupby("support_band", observed=True).correct.agg(["size", "mean"]).to_csv(
        dest / "heldout_by_score.csv"
    )
    refs.groupby("label").correct.agg(["size", "mean"]).to_csv(
        dest / "heldout_by_type.csv"
    )
    targets.to_parquet(dest / "assignments.parquet", index=False)
    write_transitions(targets, dest)
    for name in ("identity", "cell_types", "qupath", "samples"):
        (dest / name).mkdir()
    report = []
    for tag, rows in targets.groupby("image"):
        ids = rows.object_id
        case = str(int(tag.split("_")[-1]))
        d = pd.read_parquet(prior / "identity" / (tag + ".parquet")).set_index(
            "object_id", verify_integrity=True
        )
        if set(ids) != set(d.index[d.label.isin(RETIRED)]):
            raise ValueError(f"{tag}: source residuals differ")
        byid = rows.set_index("object_id")
        label = byid.destination
        meta = {
            "Tissue": ("", ""),
            "Acinar cell": ("Epithelial", "Exocrine"),
            "Ductal cell": ("Epithelial", "Exocrine"),
            "Fibroblast": ("Mesenchymal", "Fibroblast"),
            "Stromal": ("Mesenchymal", "Mural"),
            "B cell": ("Hematopoietic", "B"),
            **{s: ("Epithelial", "Endocrine") for s in HORMONES},
            **{
                s: ("Neural" if s == "Neural" else "Endothelial", v[0])
                for s, v in STRUCTURES.items()
            },
        }
        for s in label.unique():
            if s.startswith("Macrophage"):
                meta[s] = ("Hematopoietic", "Myeloid")
            if s in ("CD8a Tcell", "CD4 Tcell", "T cell"):
                meta[s] = ("Hematopoietic", "T")
            selected = label.index[label.eq(s)]
            d.loc[selected, ["label", "compartment", "lineage"]] = [s, *meta[s]]
        for k in ["compartment_basis", "lineage_basis", "named_by"]:
            d.loc[ids, k] = "cluster_reference"
        for k in [
            "compartment_status",
            "lineage_status",
            "subtype_status",
            "status",
            "confirming_markers",
            "established_label",
        ]:
            d.loc[ids, k] = ""
        d.loc[ids, "confirmed"] = np.nan
        d.loc[ids, "tags"] = "cluster_reference;score_" + byid.support_band.astype(str)
        d.reset_index().to_parquet(dest / "identity" / (tag + ".parquet"), index=False)
        table = pd.read_parquet(
            baseline / "cell_types" / (case + ".parquet")
        ).set_index("object_id", verify_integrity=True)
        if set(table.index[table.phenotype.isin(RETIRED)]) != set(ids) or not table.loc[
            ids, "phenotype"
        ].equals(byid.label):
            raise ValueError(f"{tag}: baseline UUID/phenotype mismatch")
        for k in ["compartment", "lineage", "status", "tags", "named_by", "phenotype"]:
            table.loc[ids, k] = d.loc[ids, "label" if k == "phenotype" else k]
        if table.phenotype.isin(RETIRED).any():
            raise ValueError("retired label survived")
        table.reset_index().to_parquet(
            dest / "cell_types" / (case + ".parquet"), index=False
        )
        delta = d.loc[ids, review.IDENTITY_META + ["label"]].copy()
        delta["old_label"] = byid.label
        delta["display"] = delta.label.map(review.display)
        delta["confirmed"] = ""
        delta["cluster_score"] = byid.score.round(4)
        delta["cluster_runner_up"] = byid.runner_up
        delta["cluster_id"] = (
            byid.fold.astype(str) + "/" + byid.parent + "/" + byid.child.astype(str)
        )
        delta.to_csv(dest / "qupath" / (tag + ".tsv"), sep="\t")
        rows.sort_values("margin").drop_duplicates(
            ["label", "destination"]
        ).object_id.to_csv(dest / "samples" / (tag + ".txt"), header=False, index=False)
        report.append({"image": tag, "eligible": len(table), "retyped": len(rows)})
    summary = {
        "eligible": sum(r["eligible"] for r in report),
        "retyped": len(targets),
        "matched": int(targets.matched.sum()),
        "heldout": len(refs),
        "heldout_agreement": float(refs.correct.mean()),
        "cluster_only_agreement": float(refs.cluster_correct.mean()),
        "donors": report,
    }
    (dest / "summary.json").write_text(json.dumps(summary, indent=2))


def record_method(dest, historical, baseline):
    import sklearn

    with historical.open("rb") as stream:
        historical_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    release = json.loads((baseline / "uniform_release.json").read_text())
    marker_files = release["identity_sha256"]
    targets = pd.read_parquet(dest / "assignments.parquet")
    refs = pd.read_parquet(dest / "heldout_predictions.parquet")
    named = targets.destination.ne("Tissue")
    method = {
        "historical_clusters": str(historical.resolve()),
        "historical_sha256": historical_hash,
        "reference_identity_sha256": marker_files,
        "donor_level_sha256": release["levels_sha256"],
        "numpy": np.__version__,
        "sklearn": sklearn.__version__,
        "seed": SEED,
        "folds": refs.groupby("fold").image.unique().map(list).to_dict(),
        "reference_cells_per_donor_and_phenotype_max": 500,
        "cluster_vote_weight": 0.2,
        "extra_trees_weight": 0.8,
        "score_definition": "Uncalibrated vote share before biological incompatibilities are removed; compatible votes are not renormalized.",
        "tissue": int((~named).sum()),
        "named": int(named.sum()),
        "no_compatible_vote": int(targets.no_allowed_vote.sum()),
        "no_positive_typing_marker": int(targets.no_positive_marker.sum()),
        "specific_call_below_half_vote": int((named & targets.score.lt(0.5)).sum()),
    }
    (dest / "method.json").write_text(json.dumps(method, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("typing", "prior", "baseline", "historical", "dest"):
        p.add_argument("--" + arg, type=Path, required=True)
    archive = p.add_mutually_exclusive_group()
    archive.add_argument("--replay-v25", action="store_true", help="Reproduce the archived v25 intermediate")
    archive.add_argument("--replay-v26", action="store_true", help="Reproduce v26 before the terminal T-cell subtype pass")
    p.add_argument("--endocrine-audit", type=Path, help="Independent native v25 outline and mask audit, when available")
    a = p.parse_args()
    a.dest.mkdir(parents=True, exist_ok=False)
    dest = a.dest if a.replay_v25 else a.dest / "cluster_intermediate"
    dest.mkdir(exist_ok=True)
    with threadpool_limits(limits=2):
        panel = load_sample(a.typing, a.prior, a.historical)
        predictions, stress = infer(panel)
        stress.groupby(["image", "reference", "masked", "prediction"]).size().rename(
            "cells"
        ).to_csv(dest / "masked_controls.csv")
        export(predictions, a.prior, a.baseline, dest)
        group_residuals(panel, predictions.loc[predictions.label.isin(RETIRED)], dest)
        record_method(dest, a.historical, a.baseline)
        if not a.replay_v25:
            from terminal_cluster_votes import fit
            from resolve_tissue_calls import build, combine_project_deltas
            tissue = a.dest if a.replay_v26 else a.dest / "tissue_intermediate"
            tissue.mkdir(exist_ok=True)
            fit(dest, a.typing, a.historical, tissue)
            build(dest, a.typing, a.baseline, a.endocrine_audit, tissue, source_tables=dest/"cell_types")
            combine_project_deltas(dest, tissue)
            if not a.replay_v26:
                from resolve_tcell_subtypes import build as subtype_t_cells
                subtype_t_cells(tissue, a.typing, a.historical, a.dest)
                combine_project_deltas(tissue, a.dest, transition_name="t_subtype_transitions.csv")
