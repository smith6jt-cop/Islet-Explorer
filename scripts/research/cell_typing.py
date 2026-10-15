#!/usr/bin/env python3
"""Disease-blind, donor-balanced hierarchical cell typing from corrected intensities.

The discover stage is deliberately separate from biological annotation. A Leiden
cluster is not a validated cell type. The apply stage requires an evidence-bound
annotation file and preserves unmeasured, unsupported and ambiguous cells.

Commands (all paths explicit; no historical-label or original-cohort fallback)::

    cell_typing.py discover --input-root INPUT --output-root OUTPUT --workers 14
    cell_typing.py apply --input-root INPUT --output-root OUTPUT \
        --annotations OUTPUT/typing/cluster_annotations.json
    cell_typing.py discover --input-root INPUT --output-root OUTPUT \
        --parent-typing-root OUTPUT --lineage Immune --workers 14
    cell_typing.py apply --input-root INPUT --output-root OUTPUT \
        --model-root OUTPUT/typing/refinements/Immune \
        --annotations REVIEWED_REFINEMENT.json --parent-typing-root OUTPUT

INPUT/cohort.csv contains unique case_id values. INPUT/corrected/{case_id}-cells.parquet
contains unique (case_id, object_id) keys and all 31 nonnegative corrected markers.
Optional x_um/y_um or centroid_x/centroid_y columns are carried to review candidates.
No disease labels, archived cell labels, geometry or clinical features enter fits.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import pickle
import re
import tempfile
from typing import Any

import numpy as np
import pandas as pd
from common import MARKERS, sha256, write_json

SEED = 20260919
REQUIRED_REFINEMENT_MARKERS = ['CHGA', 'CD45', 'VIM', 'LGALS3', 'CD44', 'PGP9.5', 'GAP43', 'CD35', 'NaKATPase']
ENDOCRINE_EVIDENCE = {
    'Alpha cell': 'GCG', 'Beta cell': 'INS', 'Delta cell': 'SST',
    'Hormone-unresolved endocrine cell': 'CHGA',
}
KEYS = ['case_id', 'object_id']
FIT_IDS = ['log1p_k15', 'log1p_k30', 'arcsinh_k15', 'arcsinh_k30']
COORDINATES = ['x_um', 'y_um', 'centroid_x', 'centroid_y', 'x_source_px', 'y_source_px',
               'X_centroid', 'Y_centroid']


@dataclass(frozen=True)
class TypingConfig:
    seed: int = SEED
    cells_per_donor: int = 20000
    resolution: float = 1.0
    neighbor_counts: tuple[int, ...] = (15, 30)
    transformations: tuple[str, ...] = ('log1p', 'arcsinh')
    consensus_votes: int = 3
    projection_agreement: float = 0.60
    support_quantile: float = 0.99
    batch_size: int = 50000
    refinement: str | None = None
    supplemental: bool = False
    omitted_markers: tuple[str, ...] = ()
    neighbor_backend: str = 'pynndescent'
    minimum_cell_area_um2: float = 5.0
    technical_exclusion_regions: tuple[dict, ...] = ()

    def validate(self):
        if self.refinement and self.supplemental:
            raise ValueError('Choose lineage refinement or one unsupported-cell discovery, not both')
        if self.supplemental and self.cells_per_donor > 5000:
            raise ValueError('Supplemental unsupported-cell discovery is capped at 5000 cells per donor')
        if self.cells_per_donor < 1 or self.batch_size < 1:
            raise ValueError('Positive sample and batch sizes are required')
        if len(self.neighbor_counts) * len(self.transformations) != 4:
            raise ValueError('The four-fit ensemble is required')
        if min(self.neighbor_counts) < 2 or self.consensus_votes != 3:
            raise ValueError('Require k >= 2 and three-of-four consensus')
        if not 0.5 <= self.projection_agreement <= 1 or not 0 < self.support_quantile < 1:
            raise ValueError('Invalid projection support configuration')
        if not set(self.omitted_markers).issubset(MARKERS):
            raise ValueError('Unregistered omitted marker')
        if self.neighbor_backend not in {'pynndescent', 'exact'}:
            raise ValueError('Unsupported neighbor backend')
        if not np.isfinite(self.minimum_cell_area_um2) or self.minimum_cell_area_um2 < 0:
            raise ValueError('Cell segmentation minimum area must be finite and nonnegative')
        for region in self.technical_exclusion_regions:
            if not region.get('case_id') or not region.get('geometry') or not region.get('evidence'):
                raise ValueError('Technical cell exclusions require case, source-pixel geometry and evidence binding')


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def stable_order(case_ids, object_ids, seed=SEED):
    return np.array([hashlib.sha256(f'{seed}\0{c}\0{o}'.encode()).hexdigest()
                     for c, o in zip(case_ids, object_ids)])


def atomic_parquet(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.' + path.name, suffix='.parquet', delete=False) as f:
        temp = Path(f.name)
    try:
        frame.to_parquet(temp, index=False)
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def atomic_pickle(value, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.' + path.name, delete=False) as f:
        temp = Path(f.name)
        pickle.dump(value, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    temp.replace(path)


def read_cells(path, case_id, markers=MARKERS):
    """Read only canonical measurement fields, ignoring disease/legacy labels."""
    import pyarrow.parquet as pq
    columns = pq.read_schema(path).names
    required = KEYS + list(markers)
    if not set(required).issubset(columns):
        raise ValueError(f'Missing corrected cell fields for {case_id}: {sorted(set(required) - set(columns))}')
    frame = pd.read_parquet(path, columns=required + [c for c in COORDINATES + ['Cell Area', 'Nucleus Area'] if c in columns])
    if frame[KEYS].isna().any().any():
        raise ValueError(f'Missing cell key for {case_id}')
    for col in KEYS:
        frame[col] = frame[col].astype(str)
    if set(frame.case_id) != {str(case_id)} or frame.duplicated(KEYS).any():
        raise ValueError(f'Mixed donor or duplicate cell key for {case_id}')
    values = frame[list(markers)].to_numpy(dtype=float)
    if np.any(values[np.isfinite(values)] < 0):
        raise ValueError(f'Negative corrected intensities require an explicit measurement revision: {case_id}')
    return frame.sort_values(KEYS, kind='stable').reset_index(drop=True)


def cohort_inputs(input_root):
    root = Path(input_root).resolve(strict=True)
    manifest = root / 'cohort.csv'
    cohort = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    if 'case_id' not in cohort or cohort.case_id.eq('').any() or cohort.case_id.duplicated().any():
        raise ValueError('Explicit nonempty unique cohort case_id values are required')
    paths = {cid: root / 'corrected' / f'{cid}-cells.parquet' for cid in sorted(cohort.case_id)}
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise ValueError('Missing corrected measurement inputs: ' + ', '.join(missing))
    # Deliberately exclude disease/clinical columns from fit identity. Only IDs and
    # corrected cell bytes determine inputs; complete manifest hash is provenance.
    inputs = {cid: {'path': str(path), 'sha256': sha256(path)} for cid, path in paths.items()}
    return paths, inputs, sha256(manifest)


def segmentation_quality_mask(frame, minimum_cell_area_um2=5.0):
    """Cell-area QC is an eligibility rule, never a clustering feature."""
    if minimum_cell_area_um2 <= 0:
        return np.ones(len(frame), dtype=bool)
    if 'Cell Area' not in frame:
        raise ValueError('Cell Area in calibrated µm² is required for segmentation quality eligibility')
    area = pd.to_numeric(frame['Cell Area'], errors='coerce').to_numpy(float)
    return np.isfinite(area) & (area >= minimum_cell_area_um2)


def cell_eligibility(frame, minimum_cell_area_um2=5.0, technical_exclusion_regions=()):
    """Canonical fit/projection/positivity eligibility; original rows stay intact.

    Technical polygons use zero-origin full-resolution source-pixel coordinates.
    A centroid on an exclusion boundary is also excluded. Uncertain biological
    islet-support flags are not technical exclusion polygons.
    """
    area_supported = segmentation_quality_mask(frame, minimum_cell_area_um2)
    technical = np.zeros(len(frame), dtype=bool)
    if technical_exclusion_regions:
        import shapely
        from shapely.geometry import shape
        if not {'x_source_px', 'y_source_px'}.issubset(frame):
            raise ValueError('Technical exclusion requires explicit source-pixel cell centroids')
        for region in technical_exclusion_regions:
            selected = frame.case_id.astype(str).eq(str(region['case_id'])).to_numpy()
            if not selected.any():
                continue
            geom = shape(region['geometry'])
            if not geom.is_valid or geom.is_empty or geom.area <= 0:
                raise ValueError('Invalid technical cell exclusion polygon')
            x = frame.loc[selected, 'x_source_px'].to_numpy(float)
            y = frame.loc[selected, 'y_source_px'].to_numpy(float)
            if not np.isfinite(x).all() or not np.isfinite(y).all():
                raise ValueError('Technical exclusion cannot be evaluated for missing source-pixel centroids')
            technical[np.flatnonzero(selected)] |= shapely.intersects_xy(geom, x, y)
    measured = np.isfinite(frame[MARKERS].to_numpy(float)).all(axis=1)
    quality = area_supported & ~technical
    return pd.DataFrame({'marker_measurements_complete': measured,
        'cell_area_supported': area_supported, 'technical_region_excluded': technical,
        'segmentation_quality_excluded': ~quality, 'typing_eligible': quality & measured,
        'segmentation_quality_reason': np.where(technical, 'confirmed_technical_image_region',
            np.where(area_supported, '', 'cell_area_below_minimum_or_missing'))}, index=frame.index)


def attach_quality_status(result, frame, minimum_cell_area_um2=5.0, technical_exclusion_regions=()):
    """Preserve raw rows while preventing small fragments from acquiring identity."""
    result = result.copy()
    eligibility = cell_eligibility(frame, minimum_cell_area_um2, technical_exclusion_regions)
    quality = ~eligibility.segmentation_quality_excluded.to_numpy()
    complete = eligibility.marker_measurements_complete.to_numpy()
    for col in eligibility:
        result[col] = eligibility[col].to_numpy()
    result['cell_area_um2'] = pd.to_numeric(frame['Cell Area'], errors='coerce').to_numpy(float)
    result.loc[~complete, 'classification_status'] = 'unmeasured'
    result.loc[~quality, 'classification_status'] = 'segmentation_quality_excluded'
    unavailable = ~complete | ~quality
    for column in ['base_type', 'lineage', 'phenotype_state', 'state', 'candidate_lineages']:
        result.loc[unavailable, column] = None
    if 'phenotype_states' in result:
        result.loc[unavailable, 'phenotype_states'] = pd.Series([[] for _ in range(int(unavailable.sum()))], index=result.index[unavailable], dtype=object)
    result.loc[unavailable, 'projection_supported'] = False
    result.loc[unavailable, 'is_endocrine'] = pd.NA
    result.loc[unavailable, 'confidence'] = 0.
    return result


def balanced_sample(frames, cap, seed=SEED, minimum_cell_area_um2=5.0, technical_exclusion_regions=()):
    """Equal representation among donors having eligible cells, no clinical inputs."""
    eligible = {}
    for cid, frame in frames.items():
        valid = cell_eligibility(frame, minimum_cell_area_um2, technical_exclusion_regions).typing_eligible.to_numpy()
        eligible[str(cid)] = frame.loc[valid].copy()
    represented = {cid: f for cid, f in eligible.items() if len(f)}
    if not represented:
        raise ValueError('No fully measured cells available for discovery')
    target = min(int(cap), min(map(len, represented.values())))
    selected = []
    report = []
    for cid in sorted(eligible):
        frame = eligible[cid]
        n = target if len(frame) else 0
        if n:
            order = np.argsort(stable_order(frame.case_id, frame.object_id, seed), kind='stable')[:n]
            selected.append(frame.iloc[order])
        report.append({'case_id': cid, 'eligible': len(frame), 'sampled': n,
                       'reason': '' if n else 'no_eligible_cells_for_this_discovery_scope'})
    sample = pd.concat(selected, ignore_index=True).sort_values(KEYS, kind='stable').reset_index(drop=True)
    return sample, pd.DataFrame(report)


def fit_transform(values, markers, method):
    """Fit preprocessing solely on discovery landmarks (or a training fold)."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != len(markers) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Transformation requires a finite nonnegative marker matrix')
    constant = np.ptp(values, axis=0) == 0
    if constant.all():
        raise ValueError('No varying markers in discovery sample')
    retained = np.flatnonzero(~constant)
    x = values[:, retained]
    if method == 'log1p':
        cofactor = np.ones(x.shape[1])
        x = np.log1p(x)
        center = np.median(x, axis=0)
        scale = np.subtract(*np.percentile(x, [75, 25], axis=0))
        fallback = scale == 0
        scale[fallback] = np.std(x[:, fallback], axis=0)
    elif method == 'arcsinh':
        cofactor = np.maximum(1., np.quantile(x, .9, axis=0) / 5.)
        x = np.arcsinh(x / cofactor)
        center = np.mean(x, axis=0)
        scale = np.std(x, axis=0)
        fallback = np.zeros(x.shape[1], dtype=bool)
    else:
        raise ValueError('Unknown transformation')
    if np.any(scale <= 0) or not np.isfinite(scale).all():
        raise ValueError('Nonpositive scaling for a varying marker')
    specification = {'method': method, 'input_markers': list(markers),
                     'markers': [markers[i] for i in retained],
                     'dropped_constant': [m for m, c in zip(markers, constant) if c],
                     'retained_indices': retained.tolist(), 'cofactor': cofactor.tolist(),
                     'center': center.tolist(), 'scale': scale.tolist(),
                     'sd_fallback_markers': [markers[retained[i]] for i in np.flatnonzero(fallback)]}
    return ((x - center) / scale).astype(np.float32), specification


def transform(values, specification):
    values = np.asarray(values, dtype=float)
    if values.shape[1] != len(specification['input_markers']):
        raise ValueError('Projection marker alignment differs from fitted preprocessing')
    x = values[:, specification['retained_indices']]
    if not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Nonfinite or negative projection measurements')
    if specification['method'] == 'log1p':
        x = np.log1p(x)
    else:
        x = np.arcsinh(x / np.asarray(specification['cofactor']))
    return ((x - np.asarray(specification['center'])) / np.asarray(specification['scale'])).astype(np.float32)


def make_neighbors(x, k, seed, workers, backend):
    k = min(k, len(x) - 1)
    if k < 2:
        raise ValueError('At least three discovery cells are required')
    if backend == 'exact':
        from sklearn.neighbors import NearestNeighbors
        index = NearestNeighbors(n_neighbors=k + 1, metric='euclidean', n_jobs=workers).fit(x)
        distances, neighbors = index.kneighbors(x)
    else:
        from pynndescent import NNDescent
        index = NNDescent(x, n_neighbors=k + 1, metric='euclidean', random_state=seed,
                          n_jobs=workers, low_memory=True)
        neighbors, distances = index.neighbor_graph
    # Duplicate expression profiles may not have self as the first neighbor.
    clean_neighbors = np.empty((len(x), k), np.int64)
    clean_distances = np.empty((len(x), k), float)
    for i in range(len(x)):
        valid = (neighbors[i] != i) & (neighbors[i] >= 0) & np.isfinite(distances[i])
        n, d = neighbors[i][valid], distances[i][valid]
        order = np.lexsort((n, d))[:k]
        if len(order) < k:
            raise ValueError('Neighbor graph has unsupported discovery rows')
        clean_neighbors[i], clean_distances[i] = n[order], d[order]
    return index, clean_neighbors, clean_distances


def cluster_graph(neighbors, distances, seed, resolution):
    import igraph as ig
    import leidenalg
    from scipy.sparse import coo_matrix
    n, k = neighbors.shape
    local_scale = np.maximum(distances[:, -1], np.finfo(float).eps)
    row = np.repeat(np.arange(n), k)
    col = neighbors.ravel()
    weights = np.exp(-(distances.ravel() ** 2) / (local_scale[row] * local_scale[col]))
    adjacency = coo_matrix((weights, (row, col)), shape=(n, n)).tocsr()
    adjacency = adjacency.maximum(adjacency.T)
    upper = coo_matrix(adjacency)
    keep = upper.row < upper.col
    edges = np.column_stack([upper.row[keep], upper.col[keep]])
    graph = ig.Graph(n=n, edges=edges.tolist(), directed=False)
    partition = leidenalg.find_partition(graph, leidenalg.RBConfigurationVertexPartition,
                                        weights=upper.data[keep].tolist(),
                                        resolution_parameter=resolution, seed=seed)
    labels = np.asarray(partition.membership, dtype=np.int32)
    # Stable cluster numbering by the first sorted discovery-cell index.
    mapping = {label: i for i, label in enumerate(sorted(set(labels), key=lambda c: np.flatnonzero(labels == c)[0]))}
    return np.array([mapping[v] for v in labels], dtype=np.int32), edges, upper.data[keep]


def query_neighbors(index, x, k, backend):
    if backend == 'exact':
        distances, neighbors = index.kneighbors(x, n_neighbors=k)
    else:
        neighbors, distances = index.query(x, k=k)
    return np.asarray(neighbors), np.asarray(distances)


def project_clusters(index, values, labels, specification, support, k, backend, min_agreement):
    x = transform(values, specification)
    neighbors, distances = query_neighbors(index, x, min(k, len(labels)), backend)
    valid = ((neighbors >= 0) & np.isfinite(distances)).all(axis=1)
    safe_neighbors = np.maximum(neighbors, 0)
    neighbor_labels = np.sort(labels[safe_neighbors], axis=1)
    winning = neighbor_labels[:, 0].copy()
    best_count = np.ones(len(x), dtype=np.int32)
    running = best_count.copy()
    for j in range(1, neighbor_labels.shape[1]):
        running = np.where(neighbor_labels[:, j] == neighbor_labels[:, j - 1], running + 1, 1)
        improve = running > best_count
        winning[improve] = neighbor_labels[improve, j]
        best_count[improve] = running[improve]
    agreement = best_count / neighbor_labels.shape[1]
    distance = np.median(distances, axis=1)
    limits = np.array([support[str(int(c))] for c in winning])
    supported = valid & (agreement >= min_agreement) & (distance <= limits)
    winning[~valid] = -1
    agreement[~valid] = 0.
    distance[~valid] = np.nan
    return pd.DataFrame({'cluster_id': winning, 'neighbor_agreement': agreement,
                         'projection_distance': distance, 'projection_supported': supported})


def fit_ensemble_member(sample, config, method, k, workers, destination, input_identity):
    config = TypingConfig(**config) if isinstance(config, dict) else config
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    fit_id = f'{method}_k{k}'
    identity = canonical_hash({'input': input_identity, 'config': asdict(config), 'fit_id': fit_id})
    checkpoint = destination / 'complete.json'
    if checkpoint.exists():
        done = json.loads(checkpoint.read_text())
        if done['identity'] != identity:
            raise ValueError(f'Input/configuration mismatch for existing checkpoint {destination}; use a fresh model root')
        for name, digest in done['artifacts'].items():
            if not (destination / name).is_file() or sha256(destination / name) != digest:
                raise ValueError(f'Incomplete or changed checkpoint artifact {destination / name}')
        return done
    markers = [m for m in MARKERS if m not in config.omitted_markers]
    x, preprocessing = fit_transform(sample[markers].to_numpy(float), markers, method)
    index, neighbors, distances = make_neighbors(x, k, config.seed, workers, config.neighbor_backend)
    clusters, edges, weights = cluster_graph(neighbors, distances, config.seed, config.resolution)
    # Training neighbor median distances have no self-neighbor and are the same
    # statistic used at projection. Each cluster receives its own 99% limit.
    statistic = np.median(distances, axis=1)
    support = {str(int(c)): float(np.quantile(statistic[clusters == c], config.support_quantile))
               for c in sorted(set(clusters))}
    memberships = sample[KEYS].copy()
    memberships['cluster_id'] = clusters
    memberships['training_neighbor_agreement'] = np.mean(clusters[neighbors] == clusters[:, None], axis=1)
    memberships['training_distance'] = statistic
    atomic_parquet(memberships, destination / 'memberships.parquet')
    profiles = sample[markers].copy()
    profiles.insert(0, 'cluster_id', clusters)
    median = profiles.groupby('cluster_id')[markers].median().reset_index()
    median['n_cells'] = profiles.groupby('cluster_id').size().values
    atomic_parquet(median, destination / 'cluster_profiles.parquet')
    by_donor = memberships.groupby(['cluster_id', 'case_id']).size().rename('n_cells').reset_index()
    atomic_parquet(by_donor, destination / 'donor_composition.parquet')
    np.savez_compressed(destination / 'graph.npz', edges=edges, weights=weights,
                        neighbors=neighbors, distances=distances.astype(np.float32))
    if config.neighbor_backend == 'pynndescent':
        index.prepare()
    model = {'index': index, 'cluster_labels': clusters, 'preprocessing': preprocessing,
             'support': support, 'k': min(k, len(sample) - 1), 'backend': config.neighbor_backend,
             'fit_id': fit_id, 'identity': identity, 'seed': config.seed,
             'projection_agreement': config.projection_agreement}
    atomic_pickle(model, destination / 'model.pkl')
    write_json(destination / 'preprocessing.json', preprocessing)
    write_json(destination / 'support.json', {'statistic': 'median_distance_to_k_training_neighbors',
                                            'quantile': config.support_quantile, 'limits': support,
                                            'minimum_neighbor_agreement': config.projection_agreement})
    # Concrete per-cluster review candidates: centers, random representatives and
    # low local agreement, including coordinates when the input supplies them.
    candidates = []
    for c in sorted(set(clusters)):
        ix = np.flatnonzero(clusters == c)
        center = np.median(x[ix], axis=0)
        center_rank = ix[np.argsort(np.linalg.norm(x[ix] - center, axis=1), kind='stable')[:8]]
        random_rank = ix[np.argsort(stable_order(sample.iloc[ix].case_id, sample.iloc[ix].object_id, config.seed + 1))[:8]]
        uncertain_rank = ix[np.argsort(memberships.iloc[ix].training_neighbor_agreement.to_numpy(), kind='stable')[:8]]
        for reason, ranks in [('cluster_center', center_rank), ('random', random_rank), ('low_neighbor_agreement', uncertain_rank)]:
            f = sample.iloc[ranks][KEYS + [c for c in COORDINATES if c in sample]].copy()
            f['fit_id'] = fit_id
            f['cluster_id'] = int(c)
            f['selection_reason'] = reason
            candidates.append(f)
        # Marker tails expose mixed lineages and broad-marker discordance to the
        # reference audit; this does not assign labels by an arbitrary threshold.
        for marker in REQUIRED_REFINEMENT_MARKERS:
            if marker not in markers:
                continue
            v = sample.iloc[ix][marker].to_numpy(float)
            ranks = ix[np.unique(np.r_[np.argsort(v, kind='stable')[:2], np.argsort(v, kind='stable')[-2:]])]
            f = sample.iloc[ranks][KEYS + [c for c in COORDINATES if c in sample]].copy()
            f['fit_id'], f['cluster_id'], f['selection_reason'] = fit_id, int(c), f'marker_tail:{marker}'
            candidates.append(f)
    atomic_parquet(pd.concat(candidates, ignore_index=True), destination / 'reference_candidates.parquet')
    names = ['memberships.parquet', 'cluster_profiles.parquet', 'donor_composition.parquet', 'graph.npz',
             'model.pkl', 'preprocessing.json', 'support.json', 'reference_candidates.parquet']
    done = {'identity': identity, 'fit_id': fit_id, 'n_cells': len(sample), 'n_clusters': int(len(set(clusters))),
            'markers_used': preprocessing['markers'], 'dropped_constant': preprocessing['dropped_constant'],
            'artifacts': {name: sha256(destination / name) for name in names}}
    write_json(checkpoint, done)
    return done


def restrict_to_parent(frame, parent_root, lineage=None, supplemental=False):
    parent_file = Path(parent_root) / 'typing' / 'cells' / f'{frame.case_id.iloc[0]}.parquet'
    p = pd.read_parquet(parent_file)
    for key in KEYS:
        p[key] = p[key].astype(str)
    if p.duplicated(KEYS).any() or set(map(tuple, p[KEYS].values)) != set(map(tuple, frame[KEYS].values)):
        raise ValueError('Parent typing keys differ from corrected measurement keys')
    if supplemental:
        keep = ~p.projection_supported.fillna(False) & p.classification_status.ne('unmeasured')
    else:
        keep = parent_lineage_candidates(p, lineage)
    return frame.merge(p.loc[keep, KEYS], on=KEYS, how='inner', validate='one_to_one')


def parent_lineage_candidates(parent, lineage):
    """Two supported lineage votes nominate a discovery pool, never final labels."""
    certain = parent.classification_status.isin(['validated', 'lineage_supported']) & parent.projection_supported.fillna(False) & parent.lineage.eq(lineage)
    if 'candidate_lineages' not in parent:
        return certain
    candidates = parent.candidate_lineages.map(lambda s: lineage in json.loads(s) if isinstance(s, str) and s else False)
    return certain | candidates


def discover(input_root, output_root, config=TypingConfig(), workers=1, fit_workers=1,
             model_root=None, parent_typing_root=None, exclude_donors=()):
    config.validate()
    if workers < 1 or fit_workers < 1 or fit_workers > 4 or fit_workers > workers:
        raise ValueError('Invalid shared/per-fit worker budget')
    root = Path(output_root).resolve()
    default_model = ('supplemental' if config.supplemental else 'model') if not config.refinement else 'refinements/' + safe_name(config.refinement)
    model_root = Path(model_root) if model_root else root / 'typing' / default_model
    model_root.mkdir(parents=True, exist_ok=True)
    if bool(config.refinement or config.supplemental) != bool(parent_typing_root):
        raise ValueError('Refinement/supplemental discovery requires exactly one parent typing root')
    if config.supplemental and (Path(parent_typing_root) / 'typing' / 'application.json').is_file():
        parent_application = json.loads((Path(parent_typing_root) / 'typing' / 'application.json').read_text())
        if parent_application.get('supplemental_discovery_applied'):
            raise ValueError('Only one supplemental unsupported-cell discovery is permitted')
    paths, inputs, manifest_hash = cohort_inputs(input_root)
    if set(exclude_donors) - set(paths):
        raise ValueError('Held-out donor is not in the explicit cohort')
    frames, exclusions = {}, []
    for cid, path in paths.items():
        if cid in exclude_donors:
            continue
        frame = read_cells(path, cid)
        n_all, n_measured = len(frame), int(np.isfinite(frame[MARKERS].to_numpy(float)).all(axis=1).sum())
        if config.refinement or config.supplemental:
            frame = restrict_to_parent(frame, parent_typing_root, config.refinement, config.supplemental)
        frames[cid] = frame
        exclusions.append({'case_id': cid, 'raw_objects': n_all, 'fully_measured': n_measured,
                           'unmeasured_or_partial': n_all - n_measured, 'discovery_scope_cells': len(frame),
                           'segmentation_quality_excluded_in_scope': int(cell_eligibility(frame, config.minimum_cell_area_um2, config.technical_exclusion_regions).segmentation_quality_excluded.sum())})
    sample, sampling = balanced_sample(frames, config.cells_per_donor, config.seed, config.minimum_cell_area_um2, config.technical_exclusion_regions)
    del frames
    training_inputs = {cid: value['sha256'] for cid, value in inputs.items() if cid not in exclude_donors}
    identity = canonical_hash({'inputs': training_inputs, 'keys': sample[KEYS].values.tolist(), 'config': asdict(config)})
    spec_path = model_root / 'specification.json'
    parent_hashes = None
    if parent_typing_root:
        parent_hashes = {cid: sha256(Path(parent_typing_root) / 'typing' / 'cells' / f'{cid}.parquet') for cid in training_inputs}
        identity = canonical_hash({'discovery': identity, 'parent_typing_hashes': parent_hashes})
    if spec_path.exists() and json.loads(spec_path.read_text())['identity'] != identity:
        raise ValueError('Existing model root has different inputs/configuration; preserve it and use a fresh root')
    versions = {}
    for package in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'pynndescent', 'igraph', 'leidenalg', 'pyarrow']:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = 'not installed'
    specification = {'schema_version': 1, 'identity': identity, 'config': asdict(config),
                     'input_root': str(Path(input_root).resolve()), 'cohort_manifest_sha256': manifest_hash,
                     'training_inputs': training_inputs, 'excluded_donors': sorted(exclude_donors),
                     'all_cohort_donors': sorted(paths), 'parent_typing_hashes': parent_hashes,
                     'n_landmarks': len(sample), 'software_versions': versions,
                     'review_status': 'requires_biological_cluster_annotation',
                     'features': [m for m in MARKERS if m not in config.omitted_markers],
                     'disease_labels_used': False, 'raw_objects_preserved': True}
    write_json(spec_path, specification)
    atomic_parquet(sample, model_root / 'landmarks.parquet')
    atomic_parquet(sampling, model_root / 'sampling.parquet')
    write_json(model_root / 'measurement_exclusions.json', exclusions)
    tasks = [(sample, config, method, k, max(1, workers // fit_workers), model_root / f'{method}_k{k}', identity)
             for method in config.transformations for k in config.neighbor_counts]
    if fit_workers > 1:
        with ProcessPoolExecutor(max_workers=fit_workers) as executor:
            futures = [executor.submit(fit_ensemble_member, *args) for args in tasks]
            completed = [future.result() for future in futures]
    else:
        completed = [fit_ensemble_member(*args) for args in tasks]
    annotation = {'schema_version': 1, 'model_identity': identity,
                  'note': 'Replace pending records only after joint profiles and image evidence are reviewed; no automatic biological acceptance.',
                  'clusters': []}
    for fit in completed:
        for cluster in range(fit['n_clusters']):
            annotation['clusters'].append({'fit_id': fit['fit_id'], 'cluster_id': cluster,
                'decision': 'pending', 'base_type': '', 'lineage': config.refinement or '', 'phenotype_state': '',
                'endocrine_evidence': [], 'marker_dependencies': [], 'review_note': '', 'evidence': []})
    template = model_root / 'cluster_annotations.template.json'
    if not template.exists():
        write_json(template, annotation)
    write_json(model_root / 'discovery_complete.json', {'identity': identity, 'fits': completed,
               'annotation_template': str(template), 'fit_workers': fit_workers, 'worker_budget': workers})
    return specification


def safe_name(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value))


def read_annotations(path, specification, validate_evidence=True):
    path = Path(path)
    annotations = json.loads(path.read_text())
    if annotations.get('model_identity') != specification['identity']:
        raise ValueError('Biological annotations do not belong to this fitted discovery model')
    mapping = {}
    expected_fits = {f'{method}_k{k}' for method in specification['config']['transformations']
                     for k in specification['config']['neighbor_counts']}
    for row in annotations.get('clusters', []):
        key = (str(row['fit_id']), int(row['cluster_id']))
        if key in mapping or key[0] not in expected_fits:
            raise ValueError('Duplicate or unexpected cluster annotation')
        if row.get('decision') not in {'pending', 'rejected', 'accepted'}:
            raise ValueError('Unknown cluster review decision')
        if row['decision'] == 'accepted':
            if not set(row.get('marker_dependencies', [])).issubset(MARKERS):
                raise ValueError('Unregistered biological annotation marker dependency')
            states = row.get('phenotype_states', [])
            if not isinstance(states, list) or any(not isinstance(v, str) or not v for v in states):
                raise ValueError('Overlapping phenotype_states must be a list of nonempty names')
            for name, dependencies in row.get('state_marker_dependencies', {}).items():
                if name not in states and name != row.get('phenotype_state'):
                    raise ValueError('State dependencies refer to an unassigned state')
                if not set(dependencies).issubset(MARKERS):
                    raise ValueError('Unregistered phenotype-state marker dependency')
            level = row.get('annotation_level', 'base')
            if level not in {'base', 'lineage'}:
                raise ValueError('Biological annotation level must be base or lineage')
            if (level == 'base' and not row.get('base_type')) or not row.get('lineage') or not row.get('review_note') or not row.get('evidence'):
                raise ValueError('Accepted biological labels require identity, lineage, review note and concrete evidence')
            if level == 'lineage' and (row.get('base_type') or row.get('phenotype_state') or row.get('phenotype_states')):
                raise ValueError('Intermediate lineage annotation must not assign a final base identity or state')
            if row['base_type'] in {'Unknown', 'Endocrine', 'Immune', 'T cell', 'Structural', 'Endothelial', 'Acinar', 'ECAD+', 'APCs'}:
                raise ValueError('Retired broad labels cannot become final base identities')
            parent_pool = specification['config'].get('refinement')
            if parent_pool and row['lineage'] != parent_pool and row.get('parent_pool') != parent_pool:
                raise ValueError('Cross-lineage biological refinement requires its explicit parent_pool and reviewed evidence')
            if row['base_type'] in ENDOCRINE_EVIDENCE:
                if ENDOCRINE_EVIDENCE[row['base_type']] not in row.get('endocrine_evidence', []):
                    raise ValueError('Endocrine denominator membership lacks its required endocrine evidence')
                if row['lineage'] != 'Endocrine':
                    raise ValueError('Endocrine base identity requires Endocrine lineage')
            elif row['lineage'] == 'Endocrine' and level != 'lineage':
                raise ValueError('Unregistered endocrine base identity would change the agreed denominator')
            if validate_evidence:
                for record in row['evidence']:
                    file = Path(record['path'])
                    if not file.is_absolute():
                        file = path.parent / file
                    if not file.is_file() or sha256(file) != record.get('sha256'):
                        raise ValueError(f'Annotation evidence missing or changed: {file}')
        mapping[key] = row
    return mapping


def consensus_assignments(keys, projections, mappings, minimum_votes=3):
    """Three supported, accepted votes; unreviewed clusters never become labels."""
    if len(projections) != 4:
        raise ValueError('Four ensemble projections are required')
    result = keys.copy()
    columns = {'base_type': [], 'lineage': [], 'phenotype_state': [], 'phenotype_states': [], 'classification_status': [],
               'confidence': [], 'projection_supported': [], 'is_endocrine': [],
               'consensus_votes': [], 'state_consensus_votes': [], 'candidate_lineages': [], 'lineage_vote_counts': []}
    fit_ids = sorted(projections)
    arrays = {}
    for fit_id, frame in projections.items():
        if len(frame) != len(keys):
            raise ValueError('Projection row alignment changed')
        arrays[fit_id] = (frame.cluster_id.to_numpy(), frame.projection_supported.to_numpy(bool),
                          frame.neighbor_agreement.to_numpy(float))
    for i in range(len(keys)):
        votes, states, agreements, lineage_votes, lineage_agreements = {}, {}, {}, {}, {}
        accepted_available = 0
        supported_available = 0
        for fit_id in fit_ids:
            cluster, supported, agreement = (arr[i] for arr in arrays[fit_id])
            if supported:
                supported_available += 1
            annotation = mappings.get((fit_id, int(cluster)))
            if annotation is None or annotation.get('decision') != 'accepted':
                continue
            accepted_available += 1
            if not supported:
                continue
            lineage = annotation['lineage']
            lineage_votes[lineage] = lineage_votes.get(lineage, 0) + 1
            lineage_agreements.setdefault(lineage, []).append(float(agreement))
            if annotation.get('annotation_level') == 'lineage':
                continue
            name = (annotation['base_type'], lineage)
            votes[name] = votes.get(name, 0) + 1
            agreements.setdefault(name, []).append(float(agreement))
            for state in set(annotation.get('phenotype_states') or [annotation.get('phenotype_state', '')]):
                if state:
                    states[(name, state)] = states.get((name, state), 0) + 1
        winner = max(votes, key=lambda v: (votes[v], v)) if votes else None
        count = votes.get(winner, 0)
        valid = count >= minimum_votes
        if valid:
            name, lineage = winner
            candidates = [(n, state) for (base, state), n in states.items() if base == winner and state and n >= minimum_votes]
            state_votes = max([n for n, _ in candidates], default=0)
            state_list = sorted(f'{name}: {s}' for _, s in candidates)
            state = '; '.join(state_list)
            status = 'validated'
            confidence = count / 4. * float(np.mean(agreements[winner]))
            endocrine = name in ENDOCRINE_EVIDENCE
        elif lineage_votes and max(lineage_votes.values()) >= minimum_votes:
            lineage = max(lineage_votes, key=lambda v: (lineage_votes[v], v))
            name, state, endocrine, state_votes, state_list = None, None, pd.NA, 0, []
            count = lineage_votes[lineage]
            status = 'lineage_supported'
            confidence = count / 4. * float(np.mean(lineage_agreements[lineage]))
        elif lineage_votes and max(lineage_votes.values()) >= 2:
            nominees = sorted(v for v, n in lineage_votes.items() if n >= 2)
            lineage = nominees[0] if len(nominees) == 1 else None
            name, state, endocrine, state_votes, state_list = None, None, pd.NA, 0, []
            count = max(lineage_votes.values())
            status, confidence = 'lineage_candidate', count / 4.
        else:
            name, lineage, state, endocrine = None, None, None, pd.NA
            state_list = []
            state_votes = 0
            status = ('unsupported' if supported_available < minimum_votes else
                      'unreviewed' if accepted_available < minimum_votes else 'disagreement')
            confidence = count / 4.
        for col, value in {'base_type': name, 'lineage': lineage, 'phenotype_state': state, 'phenotype_states': state_list,
            'classification_status': status, 'confidence': confidence,
            'projection_supported': supported_available >= minimum_votes,
            'is_endocrine': endocrine, 'consensus_votes': count, 'state_consensus_votes': state_votes,
            'candidate_lineages': json.dumps(sorted(v for v, n in lineage_votes.items() if n >= 2)),
            'lineage_vote_counts': json.dumps(lineage_votes, sort_keys=True)}.items():
            columns[col].append(value)
    for column, values in columns.items():
        result[column] = pd.array(values, dtype='boolean') if column == 'is_endocrine' else values
    result['state'] = result.phenotype_state
    return result


def load_model(path):
    """Only load a model produced by this pipeline, verified against its checksum."""
    path = Path(path)
    done = json.loads((path / 'complete.json').read_text())
    if sha256(path / 'model.pkl') != done['artifacts']['model.pkl']:
        raise ValueError('Cell typing model checksum mismatch')
    with (path / 'model.pkl').open('rb') as stream:
        return pickle.load(stream)


def project_frame(frame, model_root, spec, batch_size=50000, model_cache=None):
    """Project one donor; a caller-owned cache avoids reloading four large indexes.

    Fitted landmarks retain their exact membership and never need neighbor queries.
    The optional cache is scoped to a single application, never shared across runs.
    """
    projections = {}
    measured = cell_eligibility(frame, spec['config'].get('minimum_cell_area_um2', 0.0), spec['config'].get('technical_exclusion_regions', ())).typing_eligible.to_numpy()
    for method in spec['config']['transformations']:
        for k in spec['config']['neighbor_counts']:
            fit_id = f'{method}_k{k}'
            model_path = Path(model_root) / fit_id
            cache_key = str(model_path.resolve())
            cached = None if model_cache is None else model_cache.get(cache_key)
            if cached is None:
                model = load_model(model_path)
                if model['identity'] != json.loads((model_path / 'complete.json').read_text())['identity']:
                    raise ValueError('Model/checkpoint identity mismatch')
                training = pd.read_parquet(model_path / 'memberships.parquet')
                if model_cache is not None:
                    model_cache[cache_key] = (model, training)
            else:
                model, training = cached
            p = pd.DataFrame({'cluster_id': np.full(len(frame), -1, dtype=int), 'neighbor_agreement': 0.,
                              'projection_distance': np.nan, 'projection_supported': False})
            key_index = pd.MultiIndex.from_frame(frame[KEYS])
            lookup = key_index.get_indexer(pd.MultiIndex.from_frame(training[KEYS]))
            hit = lookup >= 0
            query_rows = measured.copy()
            query_rows[lookup[hit]] = False
            ix = np.flatnonzero(query_rows)
            for start in range(0, len(ix), batch_size):
                rows = ix[start:start + batch_size]
                query = frame.iloc[rows][model['preprocessing']['input_markers']].to_numpy(float)
                projected = project_clusters(model['index'], query, model['cluster_labels'], model['preprocessing'],
                    model['support'], model['k'], model['backend'], model['projection_agreement'])
                p.iloc[rows] = projected.to_numpy()
            # Landmark membership is a fitted identity, not an out-of-sample guess.
            if hit.any():
                p.loc[lookup[hit], 'cluster_id'] = training.loc[hit, 'cluster_id'].to_numpy()
                p.loc[lookup[hit], 'neighbor_agreement'] = 1.
                p.loc[lookup[hit], 'projection_distance'] = 0.
                p.loc[lookup[hit], 'projection_supported'] = True
            p['cluster_id'] = p.cluster_id.astype(int)
            p['projection_supported'] = p.projection_supported.astype(bool)
            projections[fit_id] = p
            del model
    return projections, measured


def attach_positivity(result, positivity_file=None):
    positive_columns = [f'positive__{m}' for m in MARKERS]
    if positivity_file is None:
        for column in positive_columns:
            result[column] = pd.array([pd.NA] * len(result), dtype='boolean')
        return result
    positive = pd.read_parquet(positivity_file)
    if not set(KEYS + positive_columns).issubset(positive):
        raise ValueError('Positivity file must explicitly contain all 31 nullable marker calls')
    for key in KEYS:
        positive[key] = positive[key].astype(str)
    if positive.duplicated(KEYS).any() or set(map(tuple, positive[KEYS].values)) != set(map(tuple, result[KEYS].values)):
        raise ValueError('Positivity keys differ from cell typing keys')
    for col in positive_columns:
        # Refuse numeric intensities accidentally passed as thresholded calls.
        if not positive[col].dropna().isin([True, False]).all():
            raise ValueError(f'Invalid nullable marker-positive call: {col}')
        positive[col] = positive[col].astype('boolean')
    return result.merge(positive[KEYS + positive_columns], on=KEYS, validate='one_to_one', how='left')


def merge_refinement(parent, child, lineage):
    """A reviewed refinement may resolve identity/state but cannot erase a parent."""
    for key in KEYS:
        parent[key] = parent[key].astype(str)
    parent = parent.copy()
    child = child.copy()
    index = pd.MultiIndex.from_frame(parent[KEYS])
    rows = index.get_indexer(pd.MultiIndex.from_frame(child[KEYS]))
    if np.any(rows < 0):
        raise ValueError('Refinement includes cells outside parent typing')
    if not parent_lineage_candidates(parent.iloc[rows], lineage).all():
        raise ValueError('Refinement attempts to move cells outside its lineage')
    selected = child.classification_status.eq('validated') & child.projection_supported
    already = parent.iloc[rows]
    same_identity = (already.base_type.to_numpy() == child.base_type.to_numpy()) & (already.lineage.to_numpy() == child.lineage.to_numpy())
    prior_validated = already.classification_status.eq('validated').to_numpy()
    conflict = selected.to_numpy() & (already.classification_status.eq('ambiguous_lineage_refinement').to_numpy() | (prior_validated & ~same_identity))
    consistent = selected.to_numpy() & prior_validated & same_identity
    # Identical biology reached through two provisional candidate routes is one
    # identity, not a conflict or eight-vote ensemble. Keep independent states
    # already supported within either four-fit route and the stronger confidence.
    for position in np.flatnonzero(consistent):
        old = already.iloc[position]
        new = child.iloc[position]
        old_states = old.get('phenotype_states', [])
        new_states = new.get('phenotype_states', [])
        states = sorted(set(old_states if isinstance(old_states, (list, np.ndarray)) else []) |
                        set(new_states if isinstance(new_states, (list, np.ndarray)) else []))
        child.at[child.index[position], 'phenotype_states'] = states
        child.at[child.index[position], 'phenotype_state'] = '; '.join(states)
        child.at[child.index[position], 'state'] = '; '.join(states)
        for column in ['confidence', 'consensus_votes', 'state_consensus_votes']:
            child.at[child.index[position], column] = max(old[column], new[column])
    selected &= ~conflict
    columns = ['base_type', 'lineage', 'phenotype_state', 'phenotype_states', 'state', 'classification_status', 'projection_supported', 'confidence', 'is_endocrine',
               'consensus_votes', 'state_consensus_votes']
    for column in columns:
        parent.loc[parent.index[rows[selected]], column] = child.loc[selected, column].to_numpy()
    parent['refinement_status'] = parent.get('refinement_status', pd.Series('not_selected', index=parent.index))
    parent.loc[parent.index[rows], 'refinement_status'] = child.classification_status.to_numpy()
    if conflict.any():
        changed = parent.index[rows[conflict]]
        parent.loc[changed, 'classification_status'] = 'ambiguous_lineage_refinement'
        parent.loc[changed, ['base_type', 'lineage', 'phenotype_state', 'state']] = None
        parent.loc[changed, 'phenotype_states'] = pd.Series([[] for _ in changed], index=changed, dtype=object)
        parent.loc[changed, 'is_endocrine'] = pd.NA
        parent.loc[changed, 'projection_supported'] = False
        parent.loc[changed, 'confidence'] = 0.
    parent['is_endocrine'] = parent.is_endocrine.astype('boolean')
    return parent


def merge_supplemental(parent, child):
    """Only resolve previously unsupported cells; never overwrite a reviewed type."""
    parent = parent.copy()
    for key in KEYS:
        parent[key] = parent[key].astype(str)
    lookup = pd.MultiIndex.from_frame(parent[KEYS]).get_indexer(pd.MultiIndex.from_frame(child[KEYS]))
    if np.any(lookup < 0):
        raise ValueError('Supplemental cell lies outside parent typing')
    if parent.iloc[lookup].projection_supported.fillna(False).any():
        raise ValueError('Supplemental discovery may only update unsupported cells')
    valid = child.classification_status.eq('validated') & child.projection_supported
    columns = ['base_type', 'lineage', 'phenotype_state', 'phenotype_states', 'state', 'classification_status',
               'confidence', 'projection_supported', 'is_endocrine', 'consensus_votes', 'state_consensus_votes']
    for column in columns:
        parent.loc[parent.index[lookup[valid]], column] = child.loc[valid, column].to_numpy()
    parent['is_endocrine'] = parent.is_endocrine.astype('boolean')
    parent['supplemental_resolved'] = parent.get('supplemental_resolved', pd.Series(False, index=parent.index))
    parent.loc[parent.index[lookup[valid]], 'supplemental_resolved'] = True
    return parent


def catalog_from_annotations(mapping, existing=None):
    types, states = {}, {}
    for record in (existing or {}).get('base_types', []):
        types[record['name']] = record
    for record in (existing or {}).get('states', []):
        states[(record['name'], record['parent'])] = record
    for row in mapping.values():
        if row.get('decision') != 'accepted' or row.get('annotation_level') == 'lineage':
            continue
        name = row['base_type']
        dependencies = set(types.get(name, {}).get('marker_dependencies', [])) | set(row.get('marker_dependencies', [])) | set(row.get('endocrine_evidence', []))
        types[name] = {'name': name, 'lineage': row['lineage'], 'description': row.get('description', ''),
                       'marker_dependencies': sorted(dependencies),
                       'computational_features': MARKERS}
        for raw_state in row.get('phenotype_states') or ([row['phenotype_state']] if row.get('phenotype_state') else []):
            state = f"{name}: {raw_state}"
            states[(state, name)] = {'name': state, 'parent': name, 'description': row.get('state_description', ''),
                                    'marker_dependencies': row.get('state_marker_dependencies', {}).get(raw_state, [])}
    return {'schema_version': 1, 'base_types': [types[name] for name in sorted(types)],
            'base_type_names': sorted(types), 'states': [states[key] for key in sorted(states)],
            'endocrine_population': 'Endocrine union', 'endocrine_components': list(ENDOCRINE_EVIDENCE),
            'unresolved_is_a_biological_type': False}


def apply(input_root, output_root, annotations, model_root=None, positivity_root=None,
          parent_typing_root=None, batch_size=50000):
    output_root = Path(output_root).resolve()
    model_root = Path(model_root) if model_root else output_root / 'typing' / 'model'
    spec = json.loads((model_root / 'specification.json').read_text())
    paths, inputs, manifest_hash = cohort_inputs(input_root)
    if set(paths) != set(spec['all_cohort_donors']):
        raise ValueError('Projection cohort IDs differ from the explicit discovery cohort')
    for cid, digest in spec['training_inputs'].items():
        if inputs[cid]['sha256'] != digest:
            raise ValueError(f'Corrected measurements changed since cell discovery for donor {cid}')
    mapping = read_annotations(annotations, spec)
    lineage = spec['config'].get('refinement')
    supplemental = spec['config'].get('supplemental', False)
    if bool(lineage or supplemental) != bool(parent_typing_root):
        raise ValueError('Refinement/supplemental application requires its explicit parent typing root')
    existing_catalog = None
    parent_supplemental = False
    parent_hashes = None
    if parent_typing_root:
        parent_hashes = {cid: sha256(Path(parent_typing_root) / 'typing' / 'cells' / f'{cid}.parquet') for cid in paths}
        for cid, digest in (spec.get('parent_typing_hashes') or {}).items():
            if parent_hashes.get(cid) != digest:
                raise ValueError(f'Parent cell typing changed since refinement/supplemental discovery: {cid}')
        parent_application = Path(parent_typing_root) / 'typing' / 'application.json'
        if parent_application.is_file():
            parent_supplemental = json.loads(parent_application.read_text()).get('supplemental_discovery_applied', False)
        existing_catalog = json.loads((Path(parent_typing_root) / 'typing' / 'catalog.json').read_text())
    catalog = catalog_from_annotations(mapping, existing_catalog)
    # Distinguish the features actually fit from biological annotation evidence.
    for row in catalog['base_types']:
        row['computational_features'] = spec['features']
    annotation_hash = sha256(annotations)
    positivity_hashes = None
    if positivity_root is not None:
        positivity_hashes = {cid: sha256(Path(positivity_root) / f'{cid}.parquet') for cid in paths}
    run_identity = canonical_hash({'model': spec['identity'], 'annotations': annotation_hash,
                                  'inputs': {cid: item['sha256'] for cid, item in inputs.items()},
                                  'parent_typing_hashes': parent_hashes,
                                  'positivity_hashes': positivity_hashes})
    summary = []
    model_cache = {}
    for cid, path in paths.items():
        frame = read_cells(path, cid)
        if lineage or supplemental:
            frame = restrict_to_parent(frame, parent_typing_root, lineage, supplemental)
        projections, measured = project_frame(frame, model_root, spec, batch_size, model_cache)
        result = consensus_assignments(frame[KEYS], projections, mapping)
        result = attach_quality_status(result, frame, spec['config'].get('minimum_cell_area_um2', 0.0), spec['config'].get('technical_exclusion_regions', ()))
        result['typing_model_identity'] = spec['identity']
        result['annotation_sha256'] = annotation_hash
        result['typing_run_identity'] = run_identity
        if lineage:
            result['refinement_parent_pool'] = lineage
        if not lineage and not supplemental:
            positive_file = None if positivity_root is None else Path(positivity_root) / f'{cid}.parquet'
            result = attach_positivity(result, positive_file)
        else:
            parent = pd.read_parquet(Path(parent_typing_root) / 'typing' / 'cells' / f'{cid}.parquet')
            result = merge_supplemental(parent, result) if supplemental else merge_refinement(parent, result, lineage)
            if 'refinement_model_identity' not in result:
                result['refinement_model_identity'] = None
            selected_rows = pd.MultiIndex.from_frame(result[KEYS]).get_indexer(pd.MultiIndex.from_frame(frame[KEYS]))
            result.loc[result.index[selected_rows], 'refinement_model_identity'] = spec['identity']
        # Save per-fit diagnostics independently from app's compact typed-cell table.
        diagnostic = frame[KEYS].copy()
        for fit_id, projected in projections.items():
            for col in projected:
                diagnostic[f'{fit_id}__{col}'] = projected[col].to_numpy()
        diagnostic_root = model_root / 'projections'
        atomic_parquet(diagnostic, diagnostic_root / f'{cid}.parquet')
        atomic_parquet(result, output_root / 'typing' / 'cells' / f'{cid}.parquet')
        print(json.dumps({'event': 'typing_application_donor_complete', 'case_id': cid,
            'n_objects': len(result), 'statuses': {str(k): int(v) for k, v in result.classification_status.value_counts().items()}}), flush=True)
        for status, n in result.classification_status.value_counts().items():
            summary.append({'case_id': cid, 'classification_status': status, 'n_cells': int(n),
                            'n_endocrine': int(result.loc[result.classification_status.eq(status), 'is_endocrine'].fillna(False).sum())})
    write_json(output_root / 'typing' / 'catalog.json', catalog)
    write_json(output_root / 'typing' / 'application.json', {
        'schema_version': 1, 'identity': run_identity, 'model_identity': spec['identity'],
        'model_root': str(model_root.resolve()), 'annotation_file': str(Path(annotations).resolve()),
        'annotation_sha256': annotation_hash, 'cohort_manifest_sha256': manifest_hash,
        'parent_typing_hashes': parent_hashes, 'positivity_hashes': positivity_hashes,
        'endocrine_components': list(ENDOCRINE_EVIDENCE), 'summary': summary,
        'status': 'complete', 'supplemental_discovery_applied': supplemental or parent_supplemental, 'validated_cells': sum(v['n_cells'] for v in summary if v['classification_status'] == 'validated')})
    return summary


def configure_runtime(cache_root, workers):
    cache_root = Path(cache_root)
    cache_root.mkdir(parents=True, exist_ok=True)
    for key, suffix in [('NUMBA_CACHE_DIR', 'numba'), ('MPLCONFIGDIR', 'matplotlib')]:
        path = cache_root / suffix
        path.mkdir(exist_ok=True)
        os.environ.setdefault(key, str(path.resolve()))
    # Numba uses explicit requested parallelism; BLAS must not nest worker pools.
    os.environ['NUMBA_NUM_THREADS'] = str(workers)
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS']:
        os.environ[name] = '1'
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ['discover', 'apply']:
        command = sub.add_parser(name)
        command.add_argument('--input-root', required=True, type=Path)
        command.add_argument('--output-root', required=True, type=Path)
        command.add_argument('--model-root', type=Path)
        command.add_argument('--parent-typing-root', type=Path)
        command.add_argument('--workers', type=int, default=1)
        command.add_argument('--cache-root', type=Path)
        command.add_argument('--batch-size', type=int, default=50000)
    discovery = sub.choices['discover']
    discovery.add_argument('--fit-workers', type=int, default=1)
    discovery.add_argument('--cells-per-donor', type=int, default=20000)
    discovery.add_argument('--seed', type=int, default=SEED)
    discovery.add_argument('--resolution', type=float, default=1.)
    discovery.add_argument('--minimum-cell-area-um2', type=float, default=5.0)
    discovery.add_argument('--technical-exclusions', type=Path, help='Evidence-bound source-pixel technical exclusion registry JSON')
    discovery.add_argument('--lineage')
    discovery.add_argument('--supplemental', action='store_true', help='One extra discovery pass on unsupported cells; cap 5000/donor')
    discovery.add_argument('--exclude-donor', action='append', default=[])
    discovery.add_argument('--omit-marker', action='append', default=[])
    discovery.add_argument('--neighbor-backend', choices=['pynndescent', 'exact'], default='pynndescent')
    application = sub.choices['apply']
    application.add_argument('--annotations', required=True, type=Path)
    application.add_argument('--positivity-root', type=Path)
    args = parser.parse_args()
    configure_runtime(args.cache_root or args.output_root / 'typing' / 'cache', args.workers)
    if args.command == 'discover':
        technical = ()
        if args.technical_exclusions:
            registry = json.loads(args.technical_exclusions.read_text())
            technical = tuple(registry['regions'])
            for region in technical:
                for evidence in region['evidence']:
                    if sha256(evidence['path']) != evidence['sha256']:
                        raise ValueError('Technical exclusion evidence changed')
        config = TypingConfig(seed=args.seed, cells_per_donor=args.cells_per_donor,
            resolution=args.resolution, minimum_cell_area_um2=args.minimum_cell_area_um2, technical_exclusion_regions=technical, refinement=args.lineage, supplemental=args.supplemental, omitted_markers=tuple(args.omit_marker),
            neighbor_backend=args.neighbor_backend, batch_size=args.batch_size)
        result = discover(args.input_root, args.output_root, config, args.workers, args.fit_workers,
                          args.model_root, args.parent_typing_root, args.exclude_donor)
    else:
        result = apply(args.input_root, args.output_root, args.annotations, args.model_root,
                       args.positivity_root, args.parent_typing_root, args.batch_size)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
