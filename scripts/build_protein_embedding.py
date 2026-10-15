#!/usr/bin/env python3
"""Fit an islet protein PCA graph; use the SAME graph for Leiden and UMAP.

Input X is the mean gate-rescaled protein matrix, not counts or a latent space.
Source X, obs identities/measurements, layers and raw are never transformed in
place. The output is a new H5AD; diagnostics go to a separate fresh directory.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

RESOLUTIONS = (0.3, 0.5, 0.8, 1.0)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def protein_pca(values, components=15):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or len(values) < 3 or not np.isfinite(values).all():
        raise ValueError('At least three islets with finite protein measurements required')
    mean, scale = values.mean(axis=0), values.std(axis=0)
    active = scale > 0
    if not active.any():
        raise ValueError('All protein features are constant')
    scaled = (values[:, active] - mean[active]) / scale[active]
    rank = int(np.linalg.matrix_rank(scaled))
    if rank == 0:
        raise ValueError('Protein matrix has zero rank')
    pca = PCA(n_components=min(components, rank), svd_solver='full')
    scores = pca.fit_transform(scaled)
    loadings = np.zeros((values.shape[1], scores.shape[1]))
    loadings[active] = pca.components_.T
    return scores, dict(mean=mean, scale=scale, active=active, loadings=loadings,
                        variance_ratio=pca.explained_variance_ratio_, rank=rank)


def cluster_graph(scores, seed=0, with_umap=True):
    import scanpy as sc
    graph = ad.AnnData(X=np.asarray(scores, dtype=np.float32))
    graph.obsm['X_protein_pca'] = np.asarray(scores, dtype=np.float64)
    sc.pp.neighbors(graph, n_neighbors=min(30, graph.n_obs - 1),
                    use_rep='X_protein_pca', metric='euclidean', random_state=seed)
    for resolution in RESOLUTIONS:
        sc.tl.leiden(graph, resolution=resolution, key_added=f'leiden_{resolution}',
                     random_state=seed, flavor='igraph', n_iterations=2, directed=False)
    if with_umap:
        sc.tl.umap(graph, min_dist=0.5, spread=2, random_state=seed)
    return graph


def fit_embedding(source, stability=True, donor_correction='none'):
    if not source.obs_names.is_unique or not source.var_names.is_unique:
        raise ValueError('Duplicate islet or marker identities')
    if 'imageid' not in source.obs or source.obs.imageid.isna().any():
        raise ValueError('Missing donor identities')
    values = source.X.toarray() if sparse.issparse(source.X) else np.asarray(source.X)
    scores, scaling = protein_pca(values)
    uncorrected_scores = scores.copy()
    if donor_correction == 'harmony':
        import harmonypy as hm
        metadata = pd.DataFrame({'donor': source.obs.imageid.astype(str).to_numpy()})
        harmony = hm.run_harmony(scores, metadata, 'donor', theta=2,
                                  max_iter_harmony=20, random_state=0, verbose=False)
        scores = np.asarray(harmony.Z_corr)
        if scores.shape != uncorrected_scores.shape:
            raise ValueError('Unexpected Harmony output shape')
    elif donor_correction != 'none':
        raise ValueError('Unknown donor correction')
    graph = cluster_graph(scores)
    result = source.copy()
    result.obsm['X_protein_pca'] = scores
    result.obsm['X_pca'] = uncorrected_scores
    result.obsm['X_umap'] = graph.obsm['X_umap'].copy()
    result.varm['PCs'] = scaling['loadings']
    for key in ('connectivities', 'distances'):
        result.obsp[key] = graph.obsp[key].copy()
    for key in ('neighbors', 'umap'):
        result.uns[key] = graph.uns[key].copy()
    result.uns['pca'] = dict(variance_ratio=scaling['variance_ratio'],
                             params=dict(zero_center=True, use_highly_variable=False))
    for resolution in RESOLUTIONS:
        key = f'leiden_{resolution}'
        result.obs[key] = pd.Categorical(graph.obs[key].to_numpy())
        result.uns[key] = graph.uns[key].copy()
    result.obs['leiden_umap_1'] = result.obsm['X_umap'][:, 0]
    result.obs['leiden_umap_2'] = result.obsm['X_umap'][:, 1]
    provenance = dict(method='standardized_protein_PCA_neighbors_Leiden_UMAP',
        measurement='mean gate-rescaled log1p protein intensity per islet core',
        markers=source.var_names.to_numpy(dtype=str),
        excluded_constant_markers=source.var_names[~scaling['active']].to_numpy(dtype=str),
        scaling_mean=scaling['mean'], scaling_std=scaling['scale'],
        components=scores.shape[1], numerical_rank=scaling['rank'],
        n_neighbors=min(30, source.n_obs - 1), metric='euclidean', seed=0,
        min_dist=0.5, spread=2., resolutions=np.array(RESOLUTIONS),
        leiden_flavor='igraph', leiden_iterations=2,
        donor_correction=donor_correction, harmony_theta=2. if donor_correction == 'harmony' else 0., scvi_required=False,
        graph_shared_by_clustering_and_umap=True,
        source_matrix_sha256=hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest(),
        software={p: version(p) for p in ('scanpy', 'anndata', 'numpy', 'scikit-learn', 'igraph', 'umap-learn')})
    result.uns['protein_embedding'] = provenance
    donor = source.obs.imageid.astype(str).to_numpy()
    report = dict(islets=source.n_obs, donors=len(set(donor)), components=scores.shape[1], donor_correction=donor_correction,
                  scvi_required=False, source_measurements_unchanged=np.array_equal(values, result.X),
                  resolution_diagnostics={}, seed_stability=[])
    for resolution in RESOLUTIONS:
        key = f'leiden_{resolution}'
        labels = result.obs[key].astype(str)
        fractions = pd.crosstab(labels, donor, normalize='index')
        entry = dict(clusters=int(labels.nunique()),
                     donor_normalized_mutual_information=float(normalized_mutual_info_score(donor, labels)),
                     max_cluster_donor_fraction=float(fractions.max(axis=1).max()))
        if key in source.obs and source.obs[key].notna().all():
            entry['saved_cluster_adjusted_rand'] = float(adjusted_rand_score(source.obs[key].astype(str), labels))
        report['resolution_diagnostics'][key] = entry
    # Historical scVI graph/assignments are comparisons, never fitting inputs.
    if 'connectivities' in source.obsp:
        old = source.obsp['connectivities'].tocsr().astype(bool)
        new = result.obsp['connectivities'].tocsr().astype(bool)
        overlap = old.multiply(new).sum(axis=1).A.ravel()
        union = (old + new).astype(bool).sum(axis=1).A.ravel()
        report['saved_input_graph_mean_jaccard'] = float(np.mean(np.divide(overlap, union,
            out=np.zeros_like(overlap, dtype=float), where=union > 0)))
    nn = NearestNeighbors(n_neighbors=min(30, len(scores) - 1)).fit(scores).kneighbors(return_distance=False)
    report['same_donor_neighbor_fraction'] = float(np.mean(donor[nn] == donor[:, None]))
    if stability:
        for seed in range(1, 5):
            trial = cluster_graph(scores, seed=seed, with_umap=False)
            report['seed_stability'].append(dict(seed=seed, **{
                f'leiden_{r}': float(adjusted_rand_score(result.obs[f'leiden_{r}'], trial.obs[f'leiden_{r}']))
                for r in RESOLUTIONS}))
    return result, report


def build(input_path, output_path, report_root, stability=True, donor_correction='none'):
    input_path, output_path = Path(input_path).resolve(strict=True), Path(output_path).resolve()
    report_root = Path(report_root).resolve()
    if output_path.exists() or report_root.exists():
        raise FileExistsError('Output H5AD and report directory must both be new')
    if report_root == output_path.parent or report_root.is_relative_to(output_path.parent):
        raise ValueError('Keep diagnostics outside the app data root')
    source = ad.read_h5ad(input_path)
    result, report = fit_embedding(source, stability=stability, donor_correction=donor_correction)
    report.update(input=str(input_path), input_sha256=sha256(input_path),
                  runner_sha256=sha256(__file__))
    result.uns['protein_embedding']['input_sha256'] = report['input_sha256']
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output_path, compression='gzip')
    report_root.mkdir(parents=True)
    report['output_sha256'] = sha256(output_path)
    (report_root / 'diagnostics.json').write_text(json.dumps(report, indent=2) + '\n')
    means = pd.DataFrame(result.X, index=result.obs_names, columns=result.var_names)
    means.groupby(result.obs['leiden_0.5'], observed=True).mean().to_csv(report_root / 'cluster_marker_means.csv')
    pd.crosstab(result.obs['leiden_0.5'], result.obs.imageid).to_csv(report_root / 'cluster_donor_counts.csv')
    result.obs[[f'leiden_{r}' for r in RESOLUTIONS] + ['imageid', 'leiden_umap_1', 'leiden_umap_2']].to_csv(
        report_root / 'islet_embedding.csv')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--report-root', required=True, type=Path)
    parser.add_argument('--donor-correction', choices=['none','harmony'], default='harmony',
                        help='Harmony on protein PCA is the owner-selected release method (2026-10-03); '
                             '"none" fits the uncorrected comparison')
    args = parser.parse_args()
    print(json.dumps(build(args.input, args.output, args.report_root, donor_correction=args.donor_correction), indent=2))
