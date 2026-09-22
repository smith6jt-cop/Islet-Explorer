#!/usr/bin/env python3
"""Build the saved all-donor islet UMAP and Leiden labels used by Spatial.

The embedding is descriptive and uses the same typed, endocrine-normalized
composition representation declared for the primary trajectory.  Disease
labels, clinical outcomes, section area and whole-islet marker averages are
never features.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import umap
from sklearn.neighbors import NearestNeighbors

from common import sha256, write_json


SEED = 20260916


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    source, output = args.input_root, args.output_root

    feature_path = source / "trajectory/composition_features.parquet"
    feature_manifest_path = source / "trajectory/feature_manifest.json"
    features = pd.read_parquet(feature_path)
    feature_manifest = json.loads(feature_manifest_path.read_text())
    if feature_manifest.get("representation") != "typed_composition":
        raise ValueError("Spatial embedding requires typed composition features")
    if feature_manifest.get("formula") != "log1p(validated base-type count / classifier-defined core cell count)":
        raise ValueError("Unexpected endocrine-normalized feature definition")
    columns = feature_manifest.get("primary_columns")
    if not columns or not set(columns).issubset(features.columns):
        raise ValueError("Primary typed-composition columns are missing")

    eligible_flag = features["eligible"].fillna(False).astype(bool)
    if "trajectory_feature_supported" not in features:
        raise ValueError("Classifier eligibility and trajectory feature support must be recorded separately")
    supported_flag = eligible_flag & features["trajectory_feature_supported"].fillna(False).astype(bool)
    observed_support = eligible_flag & np.isfinite(features[columns]).all(axis=1)
    if not np.array_equal(supported_flag.to_numpy(), observed_support.to_numpy()):
        raise ValueError("Recorded typed-composition support differs from the primary matrix")
    eligible = features.loc[supported_flag, ["islet_id", "case_id"]].copy().reset_index(drop=True)
    matrix = features.loc[supported_flag, columns].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    values = matrix.to_numpy(float)
    if not np.isfinite(values).all():
        missing = int((~np.isfinite(values)).sum())
        raise ValueError(f"Eligible typed-composition matrix has {missing} missing values")
    center = np.median(values, axis=0)
    scale = np.subtract(*np.percentile(values, [75, 25], axis=0))
    # Drop truly constant features, while retaining sparse identities such as
    # CD8 T cells whose IQR can be zero despite real between-islet variation.
    keep = np.ptp(values, axis=0) > 1e-12
    scale[~np.isfinite(scale) | (scale < 1e-12)] = 1.0
    if keep.sum() < 2:
        raise ValueError("Too few nonconstant typed-composition features for UMAP")
    values = (values[:, keep] - center[keep]) / scale[keep]

    reducer = umap.UMAP(
        n_neighbors=15, min_dist=0.25, metric="cosine", n_components=2,
        random_state=SEED, transform_seed=SEED, n_jobs=1,
    )
    embedding = reducer.fit_transform(values)

    neighbors = NearestNeighbors(n_neighbors=min(16, len(values)), metric="cosine").fit(values)
    distances, indices = neighbors.kneighbors(values)
    # Symmetrise the directed k-nearest-neighbour result.  Keeping edges only
    # when the query row has the smaller index silently discarded asymmetric
    # neighbour relationships.
    pair_weights: dict[tuple[int, int], float] = {}
    for row in range(len(values)):
        for pos in range(1, indices.shape[1]):
            other = int(indices[row, pos])
            edge = tuple(sorted((row, other)))
            weight = float(np.exp(-distances[row, pos]))
            pair_weights[edge] = max(weight, pair_weights.get(edge, 0.0))
    edges = list(pair_weights)
    weights = [pair_weights[edge] for edge in edges]
    graph = ig.Graph(n=len(values), edges=edges, directed=False)
    graph.es["weight"] = weights

    result = eligible.reset_index(drop=True)
    result["leiden_umap_1"] = embedding[:, 0]
    result["leiden_umap_2"] = embedding[:, 1]
    for resolution in (0.3, 0.5, 0.8, 1.0):
        partition = leidenalg.find_partition(
            graph, leidenalg.RBConfigurationVertexPartition,
            weights="weight", resolution_parameter=resolution, seed=SEED,
        )
        result[f"leiden_{resolution:.1f}"] = np.asarray(partition.membership, dtype=int).astype(str)

    destination = output / "analysis/islet_embedding.parquet"
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(destination, index=False)

    kept_columns = [column for column, use in zip(columns, keep) if use]
    dropped_columns = [column for column, use in zip(columns, keep) if not use]
    cluster_counts = {
        column: int(result[column].nunique())
        for column in result.columns if column.startswith("leiden_") and not column.startswith("leiden_umap_")
    }
    write_json(output / "analysis/islet_embedding_manifest.json", {
        "schema_version": 1,
        "representation": "typed_composition",
        "formula": feature_manifest["formula"],
        "feature_columns": kept_columns,
        "dropped_constant_columns": dropped_columns,
        "transformation": "median center and interquartile-range scale on feature-supported classifier-eligible islets; unit scale when a nonconstant sparse feature has zero IQR",
        "missing_value_policy": "classifier eligibility is retained independently; only rows with complete supported typed-composition features are embedded; no imputation",
        "classifier_eligible_islets": int(eligible_flag.sum()),
        "feature_supported_islets": int(supported_flag.sum()),
        "embedded_islets": int(len(result)),
        # Backward-compatible name retains the classifier-level meaning of
        # eligibility.  ``embedded_islets`` is the supported subset count.
        "eligible_islets": int(eligible_flag.sum()),
        "donors": int(result["case_id"].astype(str).nunique()),
        "outcomes_excluded": True,
        "excluded_feature_families": [
            "disease status", "clinical outcomes", "section area",
            "whole-islet marker averages", "existing pseudotime",
        ],
        "umap": {
            "n_neighbors": 15, "min_dist": 0.25, "metric": "cosine",
            "n_components": 2, "seed": SEED,
        },
        "leiden": {
            "graph": "symmetrized cosine 15-nearest-neighbour graph",
            "resolutions": [0.3, 0.5, 0.8, 1.0], "seed": SEED,
            "cluster_counts": cluster_counts,
        },
        "sources": {
            "trajectory/composition_features.parquet": sha256(feature_path),
            "trajectory/feature_manifest.json": sha256(feature_manifest_path),
        },
        "embedding_sha256": sha256(destination),
    })

    summary_path = output / "summary/islets.parquet"
    summary = pd.read_parquet(summary_path)
    old = [c for c in summary if c.startswith("leiden_")]
    summary = summary.drop(columns=old, errors="ignore").merge(result.drop(columns="case_id"), on="islet_id", how="left", validate="one_to_one")
    summary.to_parquet(summary_path, index=False)
    print(f"embedded {len(result)} of {int(eligible_flag.sum())} classifier-eligible islets with {values.shape[1]} nonconstant typed-composition features")


if __name__ == "__main__":
    main()
