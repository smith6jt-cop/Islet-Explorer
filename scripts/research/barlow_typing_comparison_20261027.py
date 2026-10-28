#!/usr/bin/env python3
"""Barlow et al. (eLife 2024;13:RP100535) cell typing on the uniform run, and its agreement with the pipeline's labels.

Owner, 2026-10-27: "The barlow paper used a much simpler typing method of just preprocessing and then umap/clustering.
In parallel to the current work, do this and see how well they agree with one another."

The procedure follows the authors' deposited notebook (BioStudies S-BIAD859, supplementarycode.zip,
`2025 code/umap-parametric.ipynb`), step by step:
  1. per-cell lateral-bleed-compensated whole-cell means, z-normalized within each donor (one image per donor here);
  2. downsample: MiniBatchKMeans (200 clusters, random_state 0) on all markers, 10 % of each cluster, then 50,000 cells
     (random_state 2);
  3. UMAP (n_neighbors 25, min_dist 0.02, cosine, n_epochs 1500, random_state 0) on the stage-1 markers of the 50,000
     cells; Leiden (leidenalg ModularityVertexPartition, seed 0) on UMAP's graph as an unweighted directed igraph;
  4. every other cell: UMAP transform, then a 3-nearest-neighbor classifier in the 2-D embedding gives its cluster;
  5. immune refinement: cells of the stage-1 immune clusters plus any other cell with CD45 z > 1.5, all of them,
     UMAP + Leiden on the immune markers;
  6. endocrine refinement: stage-1 endocrine cells not in step 5, 50,000 of them (random_state 2) UMAP (n_epochs 250)
     + Leiden on the endocrine markers, the rest transformed and 3-NN assigned.
Departures (recorded in the README): umap-learn's UMAP replaces ParametricUMAP (tensorflow is not installed); markers
map to this panel by antigen (no substitutes); the authors named clusters by eye, here a written rule reads each
cluster's median z-scores (the notebook's centroid statistic).

    python scripts/research/barlow_typing_comparison_20261027.py run      # clustering and labels (hours)
    python scripts/research/barlow_typing_comparison_20261027.py compare [--root <root>]  # agreement, per root
    python scripts/research/barlow_typing_comparison_20261027.py readme                   # README.md from every root
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / 'data/development/uniform-20donor-20261027'
COHORT = REPO / 'scripts/uniform/cohort.csv'
OUT = REPO / 'output/barlow_typing_comparison_20261027'
WORK = OUT / 'work'

# Barlow's channels (notebook lists) -> this panel, by antigen; None = not on this panel.
STAGE1 = {'HOECHST1/10/20, Draq 5': 'DAPI', 'cyc002 blanks (3)': None, 'NaKATPase': 'NaKATPase',
          'Cytokeratin': 'panCK', 'CD31': 'CD31', 'Podoplanin': 'PDPN', 'alphaSma': 'SMA', 'CD20': 'CD20',
          'Granzyme B': None, 'CD15': None, 'CD3': 'CD3e', 'CD68': 'CD68', 'CD163': 'CD163', 'CD45': 'CD45',
          'PD-L1': None, 'Chromogranin A': 'CHGA', 'Synaptophysin': None, 'MUC-1': None, 'CD56': None, 'CD57': None,
          'CD138': None, 'S100A6': None, 'Galectin-3': 'LGALS3'}
IMMUNE = {'cyc002 blanks (3)': None, 'HOECHST1': 'DAPI', 'BCL-2': None, 'Galectin-3': 'LGALS3', 'S100A6': None,
          'CD15': None, 'CD11c': None, 'CD163': 'CD163', 'CD56': None, 'CD45': 'CD45', 'Granzyme B': None, 'CD69': None,
          'CD68': 'CD68', 'CD138': None, 'CD20': 'CD20', 'HLA-DR': 'HLADR', 'CD206': None, 'CD3': 'CD3e', 'CD57': None,
          'FOXP3': None, 'CD4': 'CD4', 'CD8': 'CD8a', 'MPO': None, 'TCR g/d': None, 'VISTA': None, 'CD16': None}
ENDOCRINE = {'Chromogranin A': 'CHGA', 'CD31': 'CD31', 'CD45': 'CD45', 'Synaptophysin': None, 'PD-L1': None,
             'CD56': None, 'Proinsulin': 'INS', 'Glucagon': 'GCG', 'Somatostatin': 'SST'}
# Channels not acquired on every image cannot be z-normalized per donor; left out of the all-marker downsample.
NOT_ON_EVERY_IMAGE = ('KRT8-18', 'epCAM', 'Ghrelin')
CD45_IMMUNE_Z = 1.5            # the notebook's cell [41]: non-immune cells with CD45 > 1.5 join the immune refinement

# Cluster naming (replaces the authors' reading of heatmaps): the label of the defining marker with the highest
# cluster median z-score, when that median is above 0 (the donor mean); Figure 1B/C and the notebook's names.
STAGE1_DEFINING = {'CHGA': 'Endocrine', 'CD31': 'Vasculature', 'PDPN': 'Lymphatics', 'SMA': 'Smooth Muscle',
                   'LGALS3': 'Stroma', 'panCK': 'Ducts/Epithelium', 'CD45': 'Immune', 'CD3e': 'Immune',
                   'CD20': 'Immune', 'CD68': 'Immune', 'CD163': 'Immune'}
STAGE1_NONE = 'Acinar'         # Figure 1B: the Acinar row has no enriched marker
ENDO_DEFINING = {'INS': 'Beta', 'GCG': 'Alpha', 'SST': 'Delta', 'CD31': 'Vasculature', 'CD45': 'Other immune'}
ENDO_NONE = 'Endocrine (unresolved)'

UMAP_KW = dict(n_components=2, n_neighbors=25, random_state=0, min_dist=0.02, metric='cosine', low_memory=False)


def log(msg):
    print(f'{time.strftime("%H:%M:%S")} {msg}', flush=True)


def panel(mapping):
    return [m for m in dict.fromkeys(v for v in mapping.values() if v)]


# ------------------------------------------------------------------------------------------- data
def load_z() -> pd.DataFrame:
    """Analysis-eligible cells, REDSEA whole-cell means (the compensated cell mean, as CellSeg's), z per image."""
    frames = []
    for r in csv.DictReader(COHORT.open()):
        tag = r['image_name']
        src = pd.read_parquet(RUN / 'typing' / 'redsea' / f'{tag}.parquet')
        src = src[src['qc_analysis_eligible']].reset_index(drop=True)
        cols = [c for c in src.columns if c.startswith('Cell__') and c != 'Cell__n_px'
                and c[6:] not in NOT_ON_EVERY_IMAGE]
        x = src[cols].to_numpy(np.float64)
        z = (x - np.nanmean(x, axis=0)) / np.nanstd(x, axis=0)
        f = pd.DataFrame(z.astype(np.float32), columns=[c[6:] for c in cols])
        f.insert(0, 'object_id', src['object_id'].values)
        f.insert(1, 'image', tag)
        f.insert(2, 'case_id', int(r['case_id']))
        f.insert(3, 'status', r['donor_status'])
        frames.append(f)
        log(f'{tag}: {len(f):,} cells')
    z = pd.concat(frames, ignore_index=True)
    bad = z.drop(columns=['object_id', 'image', 'case_id', 'status']).isna().sum()
    if bad.any():
        raise SystemExit(f'non-finite z-scores: {bad[bad > 0].to_dict()}')
    return z


def leiden_on_graph(graph) -> np.ndarray:
    import igraph as ig
    import leidenalg as la
    src, dst = graph.nonzero()
    g = ig.Graph(directed=True)
    g.add_vertices(graph.shape[0])
    g.add_edges(list(zip(src.tolist(), dst.tolist())))
    return np.asarray(la.find_partition(g, la.ModularityVertexPartition, seed=0).membership)


def fit_cluster(x: np.ndarray, n_epochs: int, name: str):
    import umap
    t = time.time()
    emb = umap.UMAP(n_epochs=n_epochs, verbose=False, **UMAP_KW).fit(x)
    lab = leiden_on_graph(emb.graph_)
    log(f'{name}: UMAP + Leiden on {len(x):,} cells, {lab.max() + 1} clusters, {time.time() - t:.0f}s')
    return emb, lab


def transfer(emb, fit_xy, fit_lab, x_rest, name, chunk=100_000) -> tuple[np.ndarray, np.ndarray]:
    """The notebook's fit_predict_data: UMAP transform, then KNeighborsClassifier(3, euclidean) in 2-D."""
    from sklearn.neighbors import KNeighborsClassifier
    knn = KNeighborsClassifier(n_neighbors=3, metric='euclidean').fit(fit_xy, fit_lab)
    xy = np.empty((len(x_rest), 2), np.float32)
    lab = np.empty(len(x_rest), int)
    t = time.time()
    for i in range(0, len(x_rest), chunk):
        u = emb.transform(x_rest[i:i + chunk])
        xy[i:i + chunk] = u
        lab[i:i + chunk] = knn.predict(u)
        log(f'{name}: transformed {min(i + chunk, len(x_rest)):,} / {len(x_rest):,} ({time.time() - t:.0f}s)')
    return xy, lab


def medians(z: pd.DataFrame, lab: np.ndarray, markers) -> pd.DataFrame:
    m = z[markers].copy()
    m['cluster'] = lab
    out = m.groupby('cluster').median()
    out.insert(0, 'cells', m.groupby('cluster').size())
    return out


def name_by_top(med: pd.DataFrame, defining: dict, none_label: str) -> pd.Series:
    d = med[list(defining)]
    top = d.idxmax(axis=1)
    return pd.Series(np.where(d.max(axis=1) > 0, top.map(defining), none_label), index=med.index)


def name_stage1(med: pd.DataFrame) -> pd.Series:
    """name_by_top, except that Immune needs the cluster's CD45 above 0 (Figure 1B/C: every immune row is CD45-high);
    without it the highest non-immune defining marker above 0 names the cluster, else Acinar."""
    out = name_by_top(med, STAGE1_DEFINING, STAGE1_NONE)
    for c in out.index[(out == 'Immune') & (med['CD45'] <= 0)]:
        rest = med.loc[c, [m for m, v in STAGE1_DEFINING.items() if v != 'Immune']]
        out[c] = STAGE1_DEFINING[rest.idxmax()] if rest.max() > 0 else STAGE1_NONE
    return out


def name_immune(med: pd.DataFrame) -> pd.Series:
    """Figure 1C's immune labels on this panel, from each cluster's median z-scores. Every Figure 1C row is
    CD45-high: a cluster with CD45 at or below 0 is Acinar (the paper gated acinar cells out of the immune cluster)."""
    t_markers, mac = ('CD3e', 'CD4', 'CD8a'), ('CD68', 'CD163', 'HLADR')
    out = {}
    for c, r in med.iterrows():
        cand = r[list(t_markers + ('CD20',) + mac)]
        top = cand.idxmax()
        if r['CD45'] <= 0:
            out[c] = 'Acinar'
        elif cand.max() <= 0:
            out[c] = 'Other immune' if r['CD45'] > 0 else 'Acinar'   # the notebook gated acinar out of immune
        elif top in t_markers:
            if r['CD3e'] <= 0:
                out[c] = 'Immune (unresolved)'
            elif max(r['CD8a'], r['CD4']) > 0:
                out[c] = 'CD8 T cells' if r['CD8a'] > r['CD4'] else 'CD4 T cells'
            else:
                out[c] = 'T cells (unresolved)'
        elif top == 'CD20':
            out[c] = 'B Cells'
        else:
            out[c] = 'Macrophage/DC'
    return pd.Series(out)


# ------------------------------------------------------------------------------------------- run
def run(reuse_stage1: bool = False) -> int:
    from sklearn.cluster import MiniBatchKMeans
    WORK.mkdir(parents=True, exist_ok=True)
    z = load_z()
    all_markers = [c for c in z.columns if c not in ('object_id', 'image', 'case_id', 'status')]
    s1, im, en = panel(STAGE1), panel(IMMUNE), panel(ENDOCRINE)
    log(f'{len(z):,} cells; all markers {len(all_markers)}; stage 1 {s1}; immune {im}; endocrine {en}')

    km = MiniBatchKMeans(n_clusters=200, random_state=0, n_init=3).fit_predict(z[all_markers].to_numpy())
    rng = np.random.RandomState(0)
    samp = np.concatenate([rng.permutation(idx)[:int(round(0.1 * len(idx)))]
                           for idx in pd.Series(np.arange(len(z))).groupby(km).apply(lambda s: s.to_numpy())])
    fit_idx = np.sort(np.random.RandomState(2).choice(samp, 50_000, replace=False))
    log(f'downsample: {len(samp):,} cells (10 % of 200 k-means clusters), 50,000 fitted')

    prev = WORK / 'barlow_labels.parquet'
    if reuse_stage1 and prev.exists():     # stage 1 (UMAP fit, transform, 3-NN) from the earlier run of this script
        old = pd.read_parquet(prev, columns=['object_id', 'stage1_cluster', 'clustered', 'u_x', 'u_y'])
        if not old['object_id'].equals(z['object_id']):
            raise SystemExit('the saved stage 1 is for other cells')
        if not np.array_equal(np.flatnonzero(old['clustered'].to_numpy()), fit_idx):
            raise SystemExit('the saved stage 1 was fitted on other cells')
        c1 = old['stage1_cluster'].to_numpy()
        u1 = old[['u_x', 'u_y']].to_numpy(np.float32)
        lab_fit = c1[fit_idx]
        log(f'stage 1 reused from {prev.name}: {lab_fit.max() + 1} clusters')
    else:
        x1 = z[s1].to_numpy()
        emb1, lab_fit = fit_cluster(x1[fit_idx], 1500, 'stage 1')
        rest = np.setdiff1d(np.arange(len(z)), fit_idx)
        xy_rest, lab_rest = transfer(emb1, emb1.embedding_, lab_fit, x1[rest], 'stage 1')
        c1 = np.empty(len(z), int); c1[fit_idx] = lab_fit; c1[rest] = lab_rest
        u1 = np.empty((len(z), 2), np.float32); u1[fit_idx] = emb1.embedding_; u1[rest] = xy_rest
    med1 = medians(z.iloc[fit_idx], lab_fit, s1)
    name1 = name_stage1(med1)
    med1.insert(1, 'label', name1)
    med1.to_csv(WORK / 'stage1_clusters.csv')
    stage1 = name1.reindex(c1).to_numpy()

    maybe = (stage1 == 'Immune') | (z['CD45'].to_numpy() > CD45_IMMUNE_Z)
    imm_idx = np.flatnonzero(maybe)
    emb2, lab2 = fit_cluster(z[im].to_numpy()[imm_idx], 1500, 'immune')
    med2 = medians(z.iloc[imm_idx], lab2, im)
    name2 = name_immune(med2)
    med2.insert(1, 'label', name2)
    med2.to_csv(WORK / 'immune_clusters.csv')

    endo_idx = np.flatnonzero((stage1 == 'Endocrine') & ~maybe)
    n_fit = min(50_000, len(endo_idx))
    efit = np.sort(np.random.RandomState(2).choice(endo_idx, n_fit, replace=False))
    erest = np.setdiff1d(endo_idx, efit)
    xe = z[en].to_numpy()
    emb3, lab3 = fit_cluster(xe[efit], 250, 'endocrine')
    if len(erest):
        _, lab3r = transfer(emb3, emb3.embedding_, lab3, xe[erest], 'endocrine')
    c3 = pd.Series(np.concatenate([lab3, lab3r if len(erest) else []]).astype(int),
                   index=np.concatenate([efit, erest]))
    med3 = medians(z.iloc[efit], lab3, en)
    name3 = name_by_top(med3, ENDO_DEFINING, ENDO_NONE)
    med3.insert(1, 'label', name3)
    med3.to_csv(WORK / 'endocrine_clusters.csv')

    label = stage1.astype(object).copy()
    label[imm_idx] = name2.reindex(lab2).to_numpy()
    label[c3.index.to_numpy()] = name3.reindex(c3.to_numpy()).to_numpy()
    cells = z[['object_id', 'image', 'case_id', 'status']].copy()
    cells['stage1_cluster'] = c1
    cells['stage1_label'] = stage1
    cells['immune_cluster'] = -1
    cells.loc[imm_idx, 'immune_cluster'] = lab2
    cells['endocrine_cluster'] = -1
    cells.loc[c3.index, 'endocrine_cluster'] = c3.to_numpy()
    cells['barlow_label'] = label
    cells['clustered'] = False
    cells.loc[fit_idx, 'clustered'] = True
    cells['u_x'], cells['u_y'] = u1[:, 0], u1[:, 1]
    imm_xy = np.full((len(z), 2), np.nan, np.float32); imm_xy[imm_idx] = emb2.embedding_
    cells['ui_x'], cells['ui_y'] = imm_xy[:, 0], imm_xy[:, 1]
    cells.to_parquet(WORK / 'barlow_labels.parquet', index=False)
    (WORK / 'run.json').write_text(json.dumps({
        'cells': len(z), 'all_markers': all_markers, 'stage1_markers': s1, 'immune_markers': im,
        'endocrine_markers': en, 'downsample_cells': int(len(samp)), 'stage1_fit_cells': 50_000,
        'immune_cells': int(len(imm_idx)), 'endocrine_cells': int(len(endo_idx)), 'endocrine_fit_cells': int(n_fit),
        'stage1_clusters': int(lab_fit.max() + 1), 'immune_clusters': int(lab2.max() + 1),
        'endocrine_clusters': int(lab3.max() + 1), 'labels': cells['barlow_label'].value_counts().to_dict()},
        indent=1))
    log('labels: ' + json.dumps(cells['barlow_label'].value_counts().to_dict()))
    return 0


# ------------------------------------------------------------------------------------------- compare
# Common granularity. The pipeline types no structural cell (STRUCTURE_CELL_CALLS off): its Exocrine holds every
# non-endocrine, non-immune cell, so Barlow's acinar, duct/epithelium and structural labels map to one group.
PIPE_COMMON = {'Beta cell': 'Beta', 'Alpha cell': 'Alpha', 'Delta cell': 'Delta',
               'CHGA+ hormone-negative': 'Other endocrine', 'Endocrine': 'Other endocrine',
               'CD8a Tcell': 'CD8 T', 'CD4 Tcell': 'CD4 T', 'T cell': 'T (no subtype)', 'B cell': 'B',
               'Macrophage': 'Macrophage', 'Exocrine': 'Non-endocrine non-immune'}
NON = 'Non-endocrine non-immune'


def pipeline_mapping(root: Path) -> dict:
    """The root's labels (uniform_release.json) at the common granularity, from the list's order: the label just
    before Exocrine is the CD45 rule's class ('Immune' to v15, 'CD45+' from v16) -> Other immune; every label after
    Exocrine is a cell no rule claimed, named by its top marker ('<marker>+', 'CD45+ (top marker)') -> the
    non-endocrine non-immune bucket, kept apart as `detail`; a macrophage profile label -> Macrophage."""
    labels = json.loads((root / 'uniform_release.json').read_text())['labels']
    i = labels.index('Exocrine')
    cd45_rule, top = labels[i - 1], labels[i + 1:]
    if cd45_rule not in ('Immune', 'CD45+'):
        raise SystemExit(f'{root.name}: unexpected label before Exocrine: {cd45_rule}')
    def common(lab):
        if not isinstance(lab, str):
            return None
        if 'Macrophage' in lab:          # v15/v16 also hold CD68-CD163-Macrophage_..., which their list leaves out
            return 'Macrophage'
        if lab == cd45_rule:
            return 'Other immune'
        if lab in top:
            return NON
        return PIPE_COMMON.get(lab)

    def detail(lab):
        c = common(lab)
        return lab if c == NON else c
    return {'common': common, 'detail': detail, 'cd45_rule_label': cd45_rule, 'top_marker_labels': top}
BARLOW_COMMON = {'Beta': 'Beta', 'Alpha': 'Alpha', 'Delta': 'Delta', 'Endocrine (unresolved)': 'Other endocrine',
                 'CD8 T cells': 'CD8 T', 'CD4 T cells': 'CD4 T', 'T cells (unresolved)': 'T (no subtype)',
                 'B Cells': 'B', 'Macrophage/DC': 'Macrophage', 'Other immune': 'Other immune',
                 'Immune (unresolved)': 'Other immune', 'Acinar': 'Non-endocrine non-immune',
                 'Ducts/Epithelium': 'Non-endocrine non-immune', 'Vasculature': 'Non-endocrine non-immune',
                 'Lymphatics': 'Non-endocrine non-immune', 'Smooth Muscle': 'Non-endocrine non-immune',
                 'Stroma': 'Non-endocrine non-immune', 'Immune': 'Other immune', 'Endocrine': 'Other endocrine'}
ORDER = ['Beta', 'Alpha', 'Delta', 'Other endocrine', 'CD8 T', 'CD4 T', 'T (no subtype)', 'B', 'Macrophage',
         'Other immune', 'Non-endocrine non-immune']


def agreement(df: pd.DataFrame) -> dict:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    a, b = df['pipeline_common'], df['barlow_common']
    per = {}
    for lab in [l for l in ORDER if (a == l).any() or (b == l).any()]:
        pa, pb = a == lab, b == lab
        per[lab] = {'pipeline_cells': int(pa.sum()), 'barlow_cells': int(pb.sum()),
                    'pipeline_in_matching_barlow': float((pa & pb).sum() / pa.sum()) if pa.any() else None,
                    'barlow_in_matching_pipeline': float((pa & pb).sum() / pb.sum()) if pb.any() else None}
    return {'cells': int(len(df)), 'same_label': float((a == b).mean()),
            'ARI': float(adjusted_rand_score(a, b)), 'NMI': float(normalized_mutual_info_score(a, b)),
            'per_label': per}


PROFILE_MARKERS = ('CD68', 'CD163', 'LGALS3', 'CD4')
_SIGN = re.compile(r'(?<=[A-Za-z0-9.])([+-])(?=[A-Z_ (]|$)')


def sup(label) -> str:
    """A label's marker signs as superscripts for text: CD68-CD163+Macrophage_LGALS3-CD4- -> CD68⁻CD163⁺..."""
    return _SIGN.sub(lambda m: '⁺' if m.group(1) == '+' else '⁻', str(label))
TITLES = {'stage1_cluster': 'stage-1 Leiden cluster', 'immune_cluster': 'immune Leiden cluster',
          'pipeline_common': 'pipeline label', 'barlow_common': 'Barlow-method label'}


def macrophage_profiles(both: pd.DataFrame):
    """Cells both methods call macrophages: the pipeline's CD68/CD163/LGALS3/CD4 profile (from its label) against
    the profile of the cell's Barlow immune cluster (each marker + when the cluster median z is above 0)."""
    import re
    mac = both[(both['pipeline_common'] == 'Macrophage') & (both['barlow_common'] == 'Macrophage')
               & (both['immune_cluster'] >= 0)].copy()
    if not mac['phenotype'].str.contains('_').any():
        return None
    med = pd.read_csv(WORK / 'immune_clusters.csv', index_col=0)
    for m in PROFILE_MARKERS:
        sign = (med[m] > 0).map({True: '+', False: '-'})
        mac[f'barlow_{m}'] = mac['immune_cluster'].map(sign)
        pat = re.escape(m) + r'([+-])'
        mac[f'pipeline_{m}'] = mac['phenotype'].str.extract(r'(?:^|[^A-Za-z0-9])' + pat if m == 'CD4' else pat)[0]
    rows = []
    for m in PROFILE_MARKERS:
        a, b = mac[f'pipeline_{m}'], mac[f'barlow_{m}']
        rows.append({'marker': m, 'cells': len(mac), 'same_sign': float((a == b).mean()),
                     'pipeline_plus': int((a == '+').sum()), 'barlow_plus': int((b == '+').sum()),
                     'both_plus': int(((a == '+') & (b == '+')).sum()),
                     'both_minus': int(((a == '-') & (b == '-')).sum())})
    key = lambda pre: mac[[f'{pre}_{m}' for m in PROFILE_MARKERS]].apply(
        lambda r: ''.join(f'{m}{v}' for m, v in zip(PROFILE_MARKERS, r)), axis=1)
    mac['pipeline_profile'], mac['barlow_profile'] = key('pipeline'), key('barlow')
    conf = pd.crosstab(mac['pipeline_profile'], mac['barlow_profile'])
    summary = {'cells_both_macrophage': int(len(mac)),
               'same_full_profile': float((mac['pipeline_profile'] == mac['barlow_profile']).mean()),
               'by_marker': rows}
    return pd.DataFrame(rows), conf, summary


def compare(root: Path) -> int:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells = pd.read_parquet(WORK / 'barlow_labels.parquet')
    pipe = pd.concat([pd.read_parquet(p, columns=['object_id', 'phenotype']) for p in sorted((root / 'cell_types').glob('*.parquet'))])
    df = cells.merge(pipe, on='object_id', how='left', validate='one_to_one')
    mapping = pipeline_mapping(root)
    df['pipeline_common'] = df['phenotype'].map(mapping['common'])
    df['pipeline_detail'] = df['phenotype'].map(mapping['detail'])
    if df.loc[df['phenotype'].notna(), 'pipeline_common'].isna().any():
        raise SystemExit(f'unmapped pipeline labels {sorted(set(df.loc[df.phenotype.notna() & df.pipeline_common.isna(), "phenotype"]))}')
    df['barlow_common'] = df['barlow_label'].map(BARLOW_COMMON)
    if df['barlow_common'].isna().any():
        raise SystemExit(f'unmapped Barlow labels {sorted(set(df.loc[df.barlow_common.isna(), "barlow_label"]))}')
    left_out = df[df['phenotype'].isna()]
    both = df[df['phenotype'].notna()].copy()
    tag = root.name.rsplit('-', 1)[-1]                      # e.g. v13
    tables = OUT / 'tables' / tag
    tables.mkdir(parents=True, exist_ok=True)

    cm = pd.crosstab(both['pipeline_common'], both['barlow_common']).reindex(index=ORDER, columns=ORDER, fill_value=0)
    order = [l for l in ORDER if cm.loc[l].sum() or cm[l].sum()]      # v14 has no T cell without a subtype
    cm = cm.loc[order, order]
    cm.to_csv(tables / 'confusion_common_counts.csv')
    cm.div(cm.sum(axis=1), axis=0).round(4).to_csv(tables / 'confusion_common_row_proportions.csv')
    cm.div(cm.sum(axis=0), axis=1).round(4).to_csv(tables / 'confusion_common_column_proportions.csv')
    full = pd.crosstab(both['phenotype'], both['barlow_label'])
    full.to_csv(tables / 'confusion_full_labels_counts.csv')
    # the non-endocrine non-immune bucket split into Exocrine and each top-marker label (v15 on)
    top = [l for l in mapping['top_marker_labels'] if (both['pipeline_detail'] == l).any()]
    top = sorted(top, key=lambda l: -int((both['pipeline_detail'] == l).sum()))
    ext_rows = [l for l in order if l != NON] + ['Exocrine'] + top
    ext = pd.crosstab(both['pipeline_detail'], both['barlow_common']).reindex(index=ext_rows, columns=order, fill_value=0)
    ext.to_csv(tables / 'confusion_extended_counts.csv')
    ext.div(ext.sum(axis=1), axis=0).round(4).to_csv(tables / 'confusion_extended_row_proportions.csv')
    top_res = {}
    if top:
        tm = both[both['pipeline_detail'].isin(top)]
        tcm = pd.crosstab(tm['pipeline_detail'], tm['barlow_label']).reindex(index=top, fill_value=0)
        tcm = tcm[tcm.sum().sort_values(ascending=False).index]
        tcm.to_csv(tables / 'top_marker_cells_by_barlow_label_counts.csv')
        tcm.div(tcm.sum(axis=1), axis=0).round(4).to_csv(tables / 'top_marker_cells_by_barlow_label_row_proportions.csv')
        top_res = {l: {'cells': int(tcm.loc[l].sum()),
                       'barlow_common': {k: int(v) for k, v in ext.loc[l].items() if v},
                       'barlow_label': {k: int(v) for k, v in tcm.loc[l].items() if v}} for l in top}
    left_out['barlow_label'].value_counts().rename('cells').to_csv(tables / 'barlow_labels_of_cells_the_pipeline_left_out.csv')

    res = {'pipeline_root': str(root.relative_to(REPO)), 'cells_barlow': int(len(df)),
           'cells_compared': int(len(both)), 'cells_pipeline_left_out': int(len(left_out)),
           'cd45_rule_label': mapping['cd45_rule_label'], 'top_marker_labels': top_res,
           'pipeline_labels': {str(k): int(v) for k, v in both['phenotype'].value_counts().items()},
           'left_out_by_barlow_label': {str(k): int(v) for k, v in left_out['barlow_label'].value_counts().items()},
           'overall': agreement(both),
           'by_status': {s: agreement(g) for s, g in both.groupby('status')},
           'by_image': {t: agreement(g) for t, g in both.groupby('image')}}
    (OUT / f'agreement_{tag}.json').write_text(json.dumps(res, indent=1))
    pd.DataFrame([{'image': t, 'status': both.loc[both.image == t, 'status'].iat[0], 'cells': v['cells'],
                   'same_label': v['same_label'], 'ARI': v['ARI'], 'NMI': v['NMI']}
                  for t, v in res['by_image'].items()]).to_csv(tables / 'agreement_by_image.csv', index=False)
    pd.DataFrame([{'status': s, 'cells': v['cells'], 'same_label': v['same_label'], 'ARI': v['ARI'], 'NMI': v['NMI']}
                  for s, v in res['by_status'].items()]).to_csv(tables / 'agreement_by_status.csv', index=False)
    pd.DataFrame(res['overall']['per_label']).T.to_csv(tables / 'per_label_agreement.csv')
    prof = macrophage_profiles(both)
    if prof is not None:
        prof[0].to_csv(tables / 'macrophage_profile_agreement.csv', index=False)
        prof[1].to_csv(tables / 'macrophage_profile_confusion.csv')
        res['macrophage_profile'] = prof[2]
        (OUT / f'agreement_{tag}.json').write_text(json.dumps(res, indent=1))

    figs = OUT / 'figures' / tag
    figs.mkdir(parents=True, exist_ok=True)
    rp = cm.div(cm.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(11, 9))
    im = ax.imshow(rp.to_numpy(), cmap='viridis', vmin=0, vmax=1)
    ax.set_xticks(range(len(order)), order, rotation=45, ha='right', fontsize=12)
    ax.set_yticks(range(len(order)), order, fontsize=12)
    for i in range(len(order)):
        for j in range(len(order)):
            ax.text(j, i, f'{rp.iat[i, j]:.2f}\n{cm.iat[i, j]:,}', ha='center', va='center', fontsize=7,
                    color='white' if rp.iat[i, j] < 0.6 else 'black')
    ax.set_xlabel('Barlow-method label', fontsize=14)
    ax.set_ylabel('Pipeline label', fontsize=14)
    ax.set_title(f'Proportion of each pipeline label ({tag}) by Barlow-method label', fontsize=15)
    fig.colorbar(im, ax=ax, label='Row proportion')
    fig.tight_layout()
    fig.savefig(figs / 'confusion_heatmap.png', dpi=200)
    plt.close(fig)
    if top:
        erp = ext.div(ext.sum(axis=1), axis=0)
        fig, ax = plt.subplots(figsize=(11, 4 + 0.42 * len(ext_rows)))
        im = ax.imshow(erp.to_numpy(), cmap='viridis', vmin=0, vmax=1, aspect='auto')
        ax.set_xticks(range(len(order)), order, rotation=45, ha='right', fontsize=11)
        ax.set_yticks(range(len(ext_rows)), [sup(r) for r in ext_rows], fontsize=11)
        for i in range(len(ext_rows)):
            for j in range(len(order)):
                ax.text(j, i, f'{erp.iat[i, j]:.2f}\n{ext.iat[i, j]:,}', ha='center', va='center', fontsize=6,
                        color='white' if erp.iat[i, j] < 0.6 else 'black')
        ax.set_xlabel('Barlow-method label', fontsize=13)
        ax.set_ylabel('Pipeline label (non-endocrine non-immune split by its label)', fontsize=12)
        ax.set_title(f'Proportion of each pipeline label ({tag}) by Barlow-method label', fontsize=14)
        fig.colorbar(im, ax=ax, label='Row proportion')
        fig.tight_layout()
        fig.savefig(figs / 'confusion_extended_heatmap.png', dpi=200)
        plt.close(fig)

    sub = both[both['clustered']].copy()
    for col, name in (('stage1_cluster', 'umap_by_stage1_cluster.png'), ('pipeline_common', 'umap_by_pipeline_label.png'),
                      ('barlow_common', 'umap_by_barlow_label.png')):
        fig, ax = plt.subplots(figsize=(10, 9))
        cats = sorted(sub[col].unique(), key=str)
        cmap = plt.get_cmap('tab20')
        order = sub[col].map({c: i for i, c in enumerate(cats)})
        for i, c in enumerate(cats):
            s = sub[order == i]
            ax.scatter(s['u_x'], s['u_y'], s=1, color=cmap(i % 20), label=f'{c} ({len(s):,})', rasterized=True)
        ax.legend(markerscale=8, fontsize=9, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel('UMAP 1', fontsize=14)
        ax.set_ylabel('UMAP 2', fontsize=14)
        ax.set_title(f'Stage-1 UMAP of the 50,000 clustered cells, by {TITLES[col]}', fontsize=14)
        fig.tight_layout()
        fig.savefig(figs / name, dpi=200)
        plt.close(fig)
    imm = both[both['ui_x'].notna()]
    imm = imm.sample(min(len(imm), 100_000), random_state=0)
    for col, name in (('immune_cluster', 'immune_umap_by_cluster.png'), ('pipeline_common', 'immune_umap_by_pipeline_label.png')):
        fig, ax = plt.subplots(figsize=(10, 9))
        cats = sorted(imm[col].unique(), key=str)
        cmap = plt.get_cmap('tab20')
        for i, c in enumerate(cats):
            s = imm[imm[col] == c]
            ax.scatter(s['ui_x'], s['ui_y'], s=1, color=cmap(i % 20), label=f'{c} ({len(s):,})', rasterized=True)
        ax.legend(markerscale=8, fontsize=9, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel('UMAP 1', fontsize=14)
        ax.set_ylabel('UMAP 2', fontsize=14)
        ax.set_title(f'Immune-refinement UMAP (up to 100,000 cells), by {TITLES[col]}', fontsize=14)
        fig.tight_layout()
        fig.savefig(figs / name, dpi=200)
        plt.close(fig)
    print(json.dumps({k: res[k] for k in ('pipeline_root', 'cells_barlow', 'cells_compared', 'cells_pipeline_left_out')}))
    print(json.dumps({k: v for k, v in res['overall'].items() if k != 'per_label'}))
    return 0


# ------------------------------------------------------------------------------------------- readme
def md(df: pd.DataFrame, index=True) -> str:
    df = df.copy()
    if index:
        df.index = [sup(i) for i in df.index]
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(sup)
    df.columns = [sup(c) for c in df.columns]
    return df.to_markdown(index=index)


def cluster_table(name: str) -> str:
    d = pd.read_csv(WORK / name, index_col=0)
    mk = [c for c in d.columns if c not in ('cells', 'label')]
    top3 = d[mk].apply(lambda x: ', '.join(f'{m} {x[m]:.2f}' for m in x.sort_values(ascending=False).index[:3]), axis=1)
    return md(pd.DataFrame({'cells (clustered)': d['cells'], 'label': d['label'], 'top three markers (median z)': top3}))


def readme() -> int:
    """README.md from methods.md, the work/ cluster tables and every agreement_v*.json; the highest version is the
    current root and gets the full tables, the earlier roots a summary, every root with top-marker labels the
    Barlow labels of those cells."""
    runs = {p.stem.split('_')[1]: json.loads(p.read_text()) for p in OUT.glob('agreement_v*.json')}
    tags = sorted(runs, key=lambda t: int(t[1:]))
    cur = tags[-1]
    run = json.loads((WORK / 'run.json').read_text())
    t = OUT / 'tables'
    head = f"""# Barlow-method cell typing compared with the pipeline (2026-10-27)

Owner request: "The barlow paper used a much simpler typing method of just preprocessing and then umap/clustering. In parallel to the current work, do this and see how well they agree with one another."

Barlow et al.'s two-stage UMAP/Leiden typing was applied once to the {run['cells']:,} analysis-eligible cells of run `data/development/uniform-20donor-20261027` and compared cell by cell (joined on `object_id`) with each pipeline root's labels. **The current pipeline root is {cur}** (`data/development/uniform-20donor-20261027-{cur}`); the earlier roots ({', '.join(tags[:-1])}) are compared beside it, on the same Barlow labels. Measurements only; no interpretation.

Script: `scripts/research/barlow_typing_comparison_20261027.py` (`run` once; `compare --root <root>` per root; `readme`). Per-cell Barlow labels: `work/barlow_labels.parquet`.

"""
    methods = (OUT / 'methods.md').read_text()
    clusters = f"""
## Clusters

Stage 1 (50,000 clustered cells; {run['stage1_clusters']} Leiden clusters):

{cluster_table('stage1_clusters.csv')}

Immune refinement ({run['immune_cells']:,} cells, all clustered; {run['immune_clusters']} clusters):

{cluster_table('immune_clusters.csv')}

Endocrine refinement ({run['endocrine_cells']:,} cells, {run['endocrine_fit_cells']:,} clustered, the rest transferred; {run['endocrine_clusters']} clusters):

{cluster_table('endocrine_clusters.csv')}

Barlow-method labels over all cells: {', '.join(f'{k} {v:,}' for k, v in run['labels'].items())}.

## Common granularity

| Common label | Pipeline labels | Barlow-method labels |
|---|---|---|
| Beta / Alpha / Delta | Beta cell / Alpha cell / Delta cell | Beta / Alpha / Delta |
| Other endocrine | CHGA⁺ hormone-negative, Endocrine | Endocrine (unresolved) |
| CD8 T / CD4 T | CD8a Tcell / CD4 Tcell | CD8 T cells / CD4 T cells |
| T (no subtype) | T cell (v12, v13 only) | T cells (unresolved) (none named) |
| B | B cell | B Cells |
| Macrophage | Macrophage (v12, v13); from v14 the macrophage profile labels CD68⁺CD163⁻, CD68⁻CD163⁺, CD68⁺CD163⁺ or (v15 on) CD68⁻CD163⁻ Macrophage, each LGALS3⁺ or LGALS3⁻ and CD4⁺ or CD4⁻ (e.g. CD68⁻CD163⁺Macrophage_LGALS3⁺CD4⁻) | Macrophage/DC |
| Other immune | the CD45 rule's class: Immune (to v15), CD45⁺ (from v16) | Other immune, Immune (unresolved) |
| Non-endocrine non-immune | Exocrine (the pipeline types no structural cell: `STRUCTURE_CELL_CALLS` off) and, from v15, the cells no rule claimed, named by their top marker: B3TUBB⁺, CD31⁺, CD34⁺, CD4⁺, CD8a⁺, GCG⁺, INS⁺, LGALS3⁺, PDPN⁺, SMA⁺, SST⁺, VIM⁺ and CD45⁺ (v15) / CD45⁺ (top marker) (v16 on); these are also shown as their own rows below | Acinar, Ducts/Epithelium, Stroma, Lymphatics, Vasculature, Smooth Muscle |

The CD45 rule's class and the top-marker labels are read from each root's `uniform_release.json` label list (the label before Exocrine, and the labels after it). v15 and v16 also hold a thirteenth macrophage label their list leaves out, CD68⁻CD163⁻Macrophage_LGALS3⁻CD4⁻ (v15 21,654 cells, v16 17,499); it is Macrophage here like the twelve.

## Agreement

| Root | Cells compared | Left out by the pipeline | Same common label | Adjusted Rand index | Normalized mutual information |
|---|---|---|---|---|---|
""" + '\n'.join(f"| {v} | {r['cells_compared']:,} | {r['cells_pipeline_left_out']:,} | {r['overall']['same_label']:.4f} | {r['overall']['ARI']:.4f} | {r['overall']['NMI']:.4f} |" for v, r in ((v, runs[v]) for v in tags)) + "\n"

    def per_label_table(r):
        pl = pd.DataFrame(r['overall']['per_label']).T
        pl.columns = ['pipeline cells', 'Barlow cells', 'pipeline label in matching Barlow label', 'Barlow label in matching pipeline label']
        for c in pl.columns[:2]:
            pl[c] = pl[c].astype(int)
        for c in pl.columns[2:]:
            pl[c] = pl[c].astype(float).round(3)
        return pl

    def top_marker_section(v, r):
        if not r['top_marker_labels']:
            return ''
        rp = pd.read_csv(t / v / 'top_marker_cells_by_barlow_label_row_proportions.csv', index_col=0)
        rp.insert(0, 'cells', [r['top_marker_labels'][l]['cells'] for l in rp.index])
        return f"""
### {v}: the cells named by their top marker, by Barlow-method label

Row proportions over the full Barlow labels (counts in `tables/{v}/top_marker_cells_by_barlow_label_counts.csv`):

{md(rp.round(3))}
"""

    def current_section(v, r):
        cm = pd.read_csv(t / v / 'confusion_common_counts.csv', index_col=0)
        rp = pd.read_csv(t / v / 'confusion_common_row_proportions.csv', index_col=0)
        off = cm.stack()
        off = off[[a != b for a, b in off.index]].sort_values(ascending=False).head(8)
        offt = pd.DataFrame([{'pipeline label': a, 'Barlow-method label': b, 'cells': int(n),
                              'proportion of the pipeline label': round(rp.loc[a, b], 3)} for (a, b), n in off.items()])
        st = pd.read_csv(t / v / 'agreement_by_status.csv').round(3)
        im = pd.read_csv(t / v / 'agreement_by_image.csv').round(3)
        lo = r['left_out_by_barlow_label']
        out = f"""
## {v} (current root)

Cells in the Barlow universe that {v} left out (its Excluded cells, {r['cells_pipeline_left_out']:,}), by Barlow-method label: {', '.join(f'{a} {b:,}' for a, b in lo.items())}.

### Confusion matrix (rows pipeline, columns Barlow method)

Counts:

{md(cm)}

Row proportions (proportion of each pipeline label falling in each Barlow-method label):

{md(rp.round(3))}

Largest off-diagonal cells:

{md(offt, index=False)}

### Per label

{md(per_label_table(r))}

### By donor status

{md(st, index=False)}

### By image

{md(im, index=False)}
"""
        if r['top_marker_labels']:
            ext = pd.read_csv(t / v / 'confusion_extended_row_proportions.csv', index_col=0)
            cnt = pd.read_csv(t / v / 'confusion_extended_counts.csv', index_col=0)
            ext.insert(0, 'cells', cnt.sum(axis=1).astype(int))
            out += f"""
### Pipeline labels with the non-endocrine non-immune bucket split (rows), by Barlow-method common label

Row proportions, with each row's cell count (counts in `tables/{v}/confusion_extended_counts.csv`; heatmap `figures/{v}/confusion_extended_heatmap.png`):

{md(ext.round(3))}
"""
            out += top_marker_section(v, r)
        if 'macrophage_profile' in r:
            mp = pd.read_csv(t / v / 'macrophage_profile_agreement.csv').round(3)
            out += f"""
### Macrophage profile

For the {r['macrophage_profile']['cells_both_macrophage']:,} cells both methods call macrophages, each pipeline profile (CD68⁺/⁻, CD163⁺/⁻, LGALS3⁺/⁻, CD4⁺/⁻ from the {v} label) against the profile of the cell's Barlow immune cluster, where a marker is ⁺ when the cluster's median z is above 0. Barlow's clusters resolve a profile only per cluster, not per cell.

{md(mp, index=False)}

Same full four-marker profile: {r['macrophage_profile']['same_full_profile']:.3f}. Every Barlow macrophage cluster has a CD4 median above 0, so the Barlow side is CD4⁺ for all of them. Profile-by-profile counts: `tables/{v}/macrophage_profile_confusion.csv`.
"""
        return out

    def earlier_section(v, r):
        keys = [k for k in ('CD8 T', 'CD4 T', 'T (no subtype)', 'B', 'Macrophage', 'Other immune') if k in r['overall']['per_label']]
        pl = per_label_table(r).loc[keys]
        st = pd.read_csv(t / v / 'agreement_by_status.csv').round(3)
        labs = r['pipeline_labels']
        imm = {k: n for k, n in labs.items() if k in ('CD8a Tcell', 'CD4 Tcell', 'T cell', 'B cell', r['cd45_rule_label']) or 'Macrophage' in k}
        out = f"""
## {v}

Same common label {r['overall']['same_label']:.4f}, ARI {r['overall']['ARI']:.4f}, NMI {r['overall']['NMI']:.4f} on {r['cells_compared']:,} cells ({r['cells_pipeline_left_out']:,} left out by the pipeline). Immune labels in {v}: {', '.join(f'{sup(k)} {n:,}' for k, n in imm.items())}.

{md(pl)}

{md(st, index=False)}
"""
        return out + top_marker_section(v, r)

    body = head + methods + clusters + current_section(cur, runs[cur])
    for v in reversed(tags[:-1]):
        body += earlier_section(v, runs[v])
    body += f"""
## Figures

Per root under `figures/<root>/`: `confusion_heatmap.png` (row proportions with counts); `confusion_extended_heatmap.png` (from v15, the non-endocrine non-immune bucket split by label); `umap_by_stage1_cluster.png`, `umap_by_barlow_label.png`, `umap_by_pipeline_label.png` (the stage-1 UMAP of the 50,000 clustered cells); `immune_umap_by_cluster.png`, `immune_umap_by_pipeline_label.png` (the immune-refinement UMAP, 100,000 of its cells).

## Files

`agreement_<root>.json` (overall, by status, by image, per label, the top-marker cells, macrophage profile), `tables/<root>/` (confusion counts, row and column proportions, extended confusion, full-label confusion, top-marker cells by Barlow label, per label, by status, by image, left-out cells, macrophage profile), `work/` (cluster tables, per-cell labels, run logs, the first naming pass), `methods.md` (the method section of this file).
"""
    (OUT / 'README.md').write_text(body)
    print(f'README.md: roots {tags}, current {cur}, {len(body):,} characters')
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['run', 'compare', 'readme'])
    ap.add_argument('--reuse-stage1', action='store_true', help='keep the saved stage 1; redo naming and stages 2-3')
    ap.add_argument('--root', type=Path, default=REPO / 'data/development/uniform-20donor-20261027-v14')
    a = ap.parse_args()
    return {'run': lambda: run(a.reuse_stage1), 'compare': lambda: compare(a.root.resolve()), 'readme': readme}[a.step]()


if __name__ == '__main__':
    sys.exit(main())
