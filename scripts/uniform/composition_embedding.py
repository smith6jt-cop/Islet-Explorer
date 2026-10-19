#!/usr/bin/env python3
"""composition_embedding.py -- the served islet-level UMAP and Leiden clusters: endocrine composition of the objects
with ten or more core cells (owner, 2026-10-19).

The protein embedding (mean marker values, PCA, Harmony on donor; `build_protein_embedding.fit_embedding`) showed
little disease structure. Measured on v5, 2026-10-19, over each object's 30 nearest neighbours from other donors:
same disease exceeded chance by +0.025 for protein PCA + Harmony refitted on the objects with ten or more core
cells, +0.068 without Harmony (where 45 % of neighbours were the same donor, against 6.5 % by chance), and +0.098
for endocrine composition with no correction (17 % same donor). Harmony removes donor shifts, and disease is a donor
property, so it removes much of the disease difference too; cell-type shares are counted within each object and
carry no staining batch, so they need no correction. Objects under ten core cells take only a few composition values
(a two-cell object is all-alpha, all-beta or half and half), so they are left out, as they are from the Trajectory
fit (the app's TRAJ_MIN_CELLS).

Features: each endocrine label's share of the object's core endocrine cells (Beta, Alpha, Delta, CHGA+
hormone-negative, Endocrine), standardised; a share that is zero in every fitted object is left out and recorded
(Endocrine on v5: every endocrine cell has a subtype). One graph -- 30 nearest neighbours, Euclidean -- gives Leiden at
0.3 / 0.5 / 0.8 / 1.0 and the UMAP (min_dist 0.5, spread 2), seed 0 (`build_protein_embedding.cluster_graph`).
Objects outside the fit get no coordinates and no cluster (NaN). The protein embedding is kept for the record:
`X_protein_pca` (Harmony), `X_pca` (uncorrected), `X_umap_protein` and `protein_leiden_<r>`.

    python scripts/uniform/composition_embedding.py --input <root>/islet_explorer.h5ad --report report.json
re-fits a built root's H5AD and prints the diagnostics without writing it; `assemble_uniform.py` calls
`fit_composition_embedding` when it builds a root.
"""
from __future__ import annotations

import argparse
import json
import sys
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_protein_embedding import RESOLUTIONS, cluster_graph  # noqa: E402

ENDOCRINE = ('Beta cell', 'Alpha cell', 'Delta cell', 'CHGA+ hormone-negative', 'Endocrine')
MIN_CORE_CELLS = 10
NEIGHBOURS = 30


def endocrine_shares(adata) -> np.ndarray:
    """Each endocrine label's share of the object's core endocrine cells; NaN where the core has none."""
    names = list(adata.uns['phenotype_names'])
    counts = np.asarray(adata.obsm['phenotype_counts'], dtype=float)[:, [names.index(n) for n in ENDOCRINE]]
    total = counts.sum(axis=1, keepdims=True)
    return np.divide(counts, total, out=np.full_like(counts, np.nan), where=total > 0)


def cross_donor_disease(coords, status, donor, k=NEIGHBOURS) -> dict:
    """Same-disease share among each object's k nearest neighbours from other donors, against chance.

    Disease is a donor property, so neighbours from the same donor share it whatever the embedding shows; only the
    other-donor neighbours test it. Chance for an object is its disease's share of the other donors' objects."""
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(return_distance=False)[:, 1:]
    other = donor[nn] != donor[:, None]
    observed = ((status[nn] == status[:, None]) & other).sum() / other.sum()
    n_status, n_donor = pd.Series(status).value_counts(), pd.Series(donor).value_counts()
    chance = np.array([(n_status[s] - n_donor[d]) / (len(status) - n_donor[d]) for s, d in zip(status, donor)])
    weight = other.sum(axis=1)
    return dict(same_donor_neighbour_fraction=float(1 - other.mean()),
                same_donor_chance=float((pd.Series(donor).value_counts(normalize=True) ** 2).sum()),
                other_donor_same_disease=float(observed),
                other_donor_same_disease_chance=float((chance * weight).sum() / weight.sum()))


def fit_composition_embedding(source, min_core_cells=MIN_CORE_CELLS, seed=0, stability=True):
    """Return (copy of `source` with the composition UMAP and Leiden, diagnostics)."""
    shares = endocrine_shares(source)
    fit = (source.obs.total_cells_core.to_numpy() >= min_core_cells) & np.isfinite(shares).all(axis=1)
    if fit.sum() < NEIGHBOURS + 1:
        raise ValueError(f'only {fit.sum()} objects with {min_core_cells} or more core cells')
    values = shares[fit]
    mean, scale = values.mean(axis=0), values.std(axis=0)
    active = scale > 0  # a label with no cells in any fitted object (Endocrine on v5) carries nothing
    if active.sum() < 2:
        raise ValueError('fewer than two endocrine shares vary across the fitted objects')
    scaled = (values[:, active] - mean[active]) / scale[active]
    graph = cluster_graph(scaled, seed=seed)

    result = source.copy()
    # The protein embedding stays for the record under its own names; its all-object graph does not describe the
    # served layout, so the graph and its parameters are dropped.
    if 'X_umap' in result.obsm:
        result.obsm['X_umap_protein'] = result.obsm['X_umap']
    for r in RESOLUTIONS:
        if f'leiden_{r}' in result.obs:
            result.obs[f'protein_leiden_{r}'] = result.obs.pop(f'leiden_{r}')
        result.uns.pop(f'leiden_{r}', None)
    for key in ('connectivities', 'distances'):
        result.obsp.pop(key, None)
    for key in ('neighbors', 'umap'):
        result.uns.pop(key, None)

    umap = np.full((result.n_obs, 2), np.nan)
    umap[fit] = graph.obsm['X_umap']
    result.obsm['X_umap'] = umap
    result.obs['leiden_umap_1'], result.obs['leiden_umap_2'] = umap[:, 0], umap[:, 1]
    result.obs['in_composition_embedding'] = fit
    for r in RESOLUTIONS:
        labels = pd.Series(pd.NA, index=result.obs_names, dtype=object)
        labels[fit] = graph.obs[f'leiden_{r}'].astype(str).to_numpy()
        levels = sorted(labels[fit].unique(), key=int)
        result.obs[f'leiden_{r}'] = pd.Categorical(labels, categories=levels)
    result.uns['composition_embedding'] = dict(
        method='endocrine_composition_neighbors_Leiden_UMAP',
        features=np.asarray([f'{n} share of core endocrine cells' for n, k in zip(ENDOCRINE, active) if k], dtype=object),
        excluded_constant_features=np.asarray([n for n, k in zip(ENDOCRINE, active) if not k], dtype=object),
        min_core_cells=min_core_cells, objects=int(fit.sum()), objects_left_out=int((~fit).sum()),
        scaling_mean=mean, scaling_std=scale, n_neighbors=NEIGHBOURS, metric='euclidean', seed=seed,
        min_dist=0.5, spread=2., resolutions=np.array(RESOLUTIONS), leiden_flavor='igraph', leiden_iterations=2,
        donor_correction='none', graph_shared_by_clustering_and_umap=True,
        software={p: version(p) for p in ('scanpy', 'anndata', 'numpy', 'scikit-learn', 'igraph', 'umap-learn')})
    if 'protein_embedding' in result.uns:
        result.uns['protein_embedding']['role'] = ('record only: X_protein_pca, X_pca, X_umap_protein and '
                                                   'protein_leiden_<r>; the served UMAP and Leiden are '
                                                   'composition_embedding (2026-10-19)')

    obs = result.obs[fit]
    status, donor = obs.donor_status.astype(str).to_numpy(), obs.imageid.astype(str).to_numpy()
    report = dict(objects=int(fit.sum()), objects_left_out=int((~fit).sum()), min_core_cells=min_core_cells,
                  by_object_class={str(k): int(v) for k, v in obs.object_class.value_counts().items()},
                  by_donor_status={str(k): int(v) for k, v in obs.donor_status.value_counts().items()},
                  donors=int(len(set(donor))), umap=cross_donor_disease(graph.obsm['X_umap'], status, donor),
                  resolution_diagnostics={}, seed_stability=[])
    for r in RESOLUTIONS:
        labels = graph.obs[f'leiden_{r}'].astype(str).to_numpy()
        report['resolution_diagnostics'][f'leiden_{r}'] = dict(
            clusters=int(len(set(labels))),
            donor_normalized_mutual_information=float(normalized_mutual_info_score(donor, labels)),
            disease_normalized_mutual_information=float(normalized_mutual_info_score(status, labels)),
            max_cluster_donor_fraction=float(pd.crosstab(labels, donor, normalize='index').max(axis=1).max()))
    if stability:
        for s in range(1, 5):
            trial = cluster_graph(scaled, seed=s, with_umap=False)
            report['seed_stability'].append(dict(seed=s, **{
                f'leiden_{r}': float(adjusted_rand_score(graph.obs[f'leiden_{r}'], trial.obs[f'leiden_{r}']))
                for r in RESOLUTIONS}))
    return result, report


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', type=Path, required=True, help="a built root's islet_explorer.h5ad")
    ap.add_argument('--report', type=Path, help='write the diagnostics here as well')
    a = ap.parse_args()
    _, rep = fit_composition_embedding(ad.read_h5ad(a.input))
    text = json.dumps(rep, indent=2) + '\n'
    if a.report:
        a.report.write_text(text)
    print(text)
