#!/usr/bin/env python3
"""assemble_uniform.py -- build a complete app data root for all 20 donors from the uniform run.

Every value is regenerated from the uniform workflow (`$OUT`, see README.md): the accepted endocrine objects
(`$OUT/islets/`), the structure objects (`$OUT/structures/`), the cells and their Islet-Explorer-Senior typing
(`$OUT/typing/`). Nothing is appended to the frozen 15-donor release. Owner decisions applied:

  A5  every endocrine object is an analysis row -- Islet_N, Cluster_N and Single_N -- with `object_class`
      (Islet / EndocrineCluster / EndocrineSingle) so the app can filter them.
  A6  nerve, blood-vessel and lymphatic cells come only from the structure objects (the typing).
  D3  cell types are Senior's classes, shown under the app's labels (APP_LABEL below).
  D4  hormone positivity is the image's hormone level on the scaled hormone (Senior's endocrine admission);
      regional immune positivity and the `<marker>_pos` flags are scaled >= 0.10 (Senior's rule level).
  D5  marker values are the typing's: compartment REDSEA, UniFORM-aligned, cohort floor/top scale to [0, 1].

Regions: an object's core is the cells whose centroid lies in it; its peri band is the 20 um ring around it,
inside Tissue, minus every endocrine object; a cell may sit in several bands. Cells ruled out of analysis (no
nucleus) and Excluded cells (exocrine inside an object, or a Single without a unique hormone) are counted in
`census.csv` and kept out of every composition, neighbourhood and display table.

Donor demographics and autoantibodies, and the viewer image links, are carried from the last served root
(`--donor-source`, v3), which recorded them from nPOD; nothing else is read from it.

    python scripts/uniform/assemble_uniform.py --run data/development/uniform-20donor-20261015 \\
        --donor-source data/development/protein-expanded-20261015-v3 \\
        --output data/development/uniform-20donor-20261016-v1 --report-root tmp/uniform-20donor-20261016-v1
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from shapely.geometry import mapping, shape

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / 'typing'))
sys.path.insert(0, str(REPO / 'scripts'))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402
from excluded_donors import assert_no_excluded  # noqa: E402

PX = CT.PIXEL_UM
PERI_UM = 20.0
MARKERS = ['CD31', 'CD8a', 'CD3e', 'SMA', 'CD4', 'CD34', 'HLADR', 'PDPN', 'panCK', 'ECAD', 'CD163', 'SST', 'ColIV',
           'VIM', 'CD20', 'LGALS3', 'B3TUBB', 'GCG', 'KRT14', 'GAP43', 'CD35', 'CHGA', 'PGP9.5', 'INS', 'CD44',
           'CD45', 'NaKATPase', 'BCatenin', 'CD68', 'BActin', 'CK19']          # the 31-marker panel, app order
LABELS = ['Beta cell', 'Alpha cell', 'Delta cell', 'Endocrine', 'CD8a Tcell', 'CD4 Tcell', 'T cell', 'Macrophage',
          'Immune', 'Blood Vessel', 'Neural', 'Lymphatic', 'Exocrine']
ENDOCRINE = ('Beta cell', 'Alpha cell', 'Delta cell', 'Endocrine')
IMMUNE = ('CD8a Tcell', 'CD4 Tcell', 'T cell', 'Macrophage', 'Immune')
POS_MARKERS = ['INS', 'GCG', 'SST', 'CD3e', 'CD8a', 'CD4', 'CD20', 'CD68', 'CD163', 'HLADR', 'CD31', 'CD34', 'SMA',
               'PDPN', 'ECAD', 'panCK', 'CD45', 'CHGA']
REGION_MARKERS = ['CD3e', 'CD163', 'CD68', 'LGALS3']
MORPHOLOGY = {'Nucleus Area': 'Nucleus: Area µm^2', 'Nucleus Length': 'Nucleus: Length µm',
              'Nucleus Circularity': 'Nucleus: Circularity', 'Nucleus Solidity': 'Nucleus: Solidity',
              'Nucleus: Max diameter µm': 'Nucleus: Max diameter µm',
              'Nucleus: Min diameter µm': 'Nucleus: Min diameter µm', 'Cell Area': 'Cell: Area µm^2',
              'Cell Length': 'Cell: Length µm', 'Cell Circularity': 'Cell: Circularity',
              'Cell Solidity': 'Cell: Solidity', 'Cell: Max diameter µm': 'Cell: Max diameter µm',
              'Cell: Min diameter µm': 'Cell: Min diameter µm'}
STRUCTURE_CLASSES = ('Nerve', 'Capillary', 'Lymphatic')


def safe(label: str) -> str:
    return label.replace(' ', '_').replace('+', 'plus')


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def features(path: Path) -> list:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as fh:
        raw = json.load(fh)
    return raw['features'] if isinstance(raw, dict) else raw


def valid(g):
    return g if g.is_valid else shapely.make_valid(g)


# ------------------------------------------------------------------------------------------- the label
def app_labels(calls: pd.DataFrame, immune: pd.DataFrame) -> np.ndarray:
    """Senior's classes under the app's labels (D3)."""
    cls = calls['cell_class'].to_numpy(dtype=object)
    out = np.full(len(cls), '', dtype=object)
    direct = {'Beta': 'Beta cell', 'Alpha': 'Alpha cell', 'Delta': 'Delta cell',
              'Endocrine (unresolved subtype)': 'Endocrine', 'Endocrine (ambiguous subtype)': 'Endocrine',
              'Blood vessel': 'Blood Vessel', 'Nerve': 'Neural', 'Lymphatic': 'Lymphatic', 'Exocrine': 'Exocrine',
              'Excluded': 'Excluded'}
    for k, v in direct.items():
        out[cls == k] = v
    imm = cls == 'Immune'
    sub = calls[['object_id']].merge(immune, on='object_id', how='left')
    ic = sub['immune_class'].fillna('').to_numpy(dtype=object)
    st = sub['immune_subtype'].fillna('').to_numpy(dtype=object)
    definition = calls['immune_definition'].fillna('').to_numpy(dtype=object)
    use = np.where(ic != '', ic, definition)                 # most specific rule; else the chain's definition
    lab = np.full(len(cls), 'Immune', dtype=object)          # CD45-only, supporting tier, no rule fired
    lab[use == 'Imm-T-cell'] = 'T cell'
    lab[(use == 'Imm-T-cell') & (st == 'T-CD8')] = 'CD8a Tcell'
    lab[(use == 'Imm-T-cell') & (st == 'T-CD4-helper')] = 'CD4 Tcell'
    lab[use == 'Imm-macrophage'] = 'Macrophage'
    out[imm] = lab[imm]
    if (out == '').any():
        raise ValueError(f'unmapped classes: {sorted(set(cls[out == ""]))}')
    return out


# ------------------------------------------------------------------------------------------- one image
def donor_fields(donors: pd.DataFrame, case: str) -> dict:
    r = donors.loc[donors.imageid.astype(int) == int(case)].iloc[0]
    return {'Case ID': f'{int(case):04d}', 'Donor Status': r.donor_status, 'Age': float(r.age), 'Gender': str(r.gender),
            'AAb_GADA': float(r.GADA), 'AAb_IA2A': float(r.IA2A), 'AAb_ZnT8A': float(r.ZnT8A), 'AAb_mIAA': float(r.mIAA)}


def build_image(run: Path, out: Path, tag: str, case: str, donors: pd.DataFrame, frequencies: dict) -> dict:
    t0 = time.time()
    typ = CT.typing_dir(run)
    meta = donor_fields(donors, case)
    cid = str(int(case))
    scale = CL.load_scale(run)
    calls = pd.read_parquet(typ / 'calls' / f'{tag}.parquet')
    immune = pd.read_parquet(typ / 'immune' / f'{tag}.parquet', columns=['object_id', 'immune_class', 'immune_subtype'])
    audit = json.loads((typ / 'calls' / f'{tag}.json').read_text())
    refined = pd.read_parquet(typ / 'refine' / f'{tag}.parquet', columns=['object_id', 'refined_class'])
    aligned = pd.read_parquet(typ / 'aligned' / f'{tag}.parquet')
    raw = pd.read_parquet(run / 'cells' / f'{tag}.parquet',
                          columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'] + list(MORPHOLOGY.values()))
    census = {'image': tag, 'case_id': cid, 'objects_detected': int(len(aligned)),
              'not_analysis_eligible': int((~aligned['qc_analysis_eligible']).sum())}

    cells = calls.merge(refined, on='object_id', how='left', validate='one_to_one')
    cells['phenotype'] = app_labels(cells, immune)
    census['excluded'] = int((cells['phenotype'] == 'Excluded').sum())
    excluded_in = cells.loc[cells['phenotype'] == 'Excluded', 'islet_num'].value_counts()
    census['excluded_reasons'] = {str(k): int(v) for k, v in
                                  cells.loc[cells.phenotype == 'Excluded', 'excluded_reason'].value_counts().items()}
    cells = cells[cells['phenotype'] != 'Excluded'].reset_index(drop=True)
    cells = cells.merge(raw.rename(columns={'Object ID': 'object_id'}), on='object_id', how='left', validate='one_to_one')
    vals = cells[['object_id']].merge(aligned[['object_id'] + MARKERS], on='object_id', how='left')
    X = CL.apply_scale(vals[MARKERS].to_numpy(), MARKERS, scale)       # cells x 31, [0, 1], NaN where unmeasured
    labels = cells['phenotype'].to_numpy(dtype=object)
    xy_um = cells[['Centroid X µm', 'Centroid Y µm']].to_numpy(float)
    pts = shapely.points(xy_um / PX)
    census['analysed'] = int(len(cells))

    # objects, tissue, structures (pixel coordinates)
    objs = []
    for f in features(run / 'islets' / f'{tag}.geojson'):
        p = f.get('properties') or {}
        c = (p.get('classification') or {}).get('name')
        if c in CT.ENDOCRINE_OBJECTS and f.get('geometry'):
            objs.append((p['name'], c, valid(shape(f['geometry']))))
    tissue = shapely.union_all([valid(shape(f['geometry'])) for f in features(run / 'inputs' / 'tissue' / f'{tag}.geojson')
                                if ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue'])
    shapely.prepare(tissue)
    structures = {}
    for cls in STRUCTURE_CLASSES:
        g = [valid(shape(f['geometry'])) for f in features(run / 'structures' / cls / f'{cls}__{tag}.geojson')
             if f.get('geometry')]
        structures[cls] = (np.asarray(g, dtype=object), STRtree(g)) if g else None

    names = [o[0] for o in objs]
    geoms = np.asarray([o[2] for o in objs], dtype=object)
    otree = STRtree(geoms)
    in_obj = calls.set_index('object_id').reindex(cells['object_id'])['islet_num'].to_numpy(dtype=object)
    index = {n: i for i, n in enumerate(names)}
    core_owner = np.array([index.get(k, -1) if isinstance(k, str) else -1 for k in in_obj])
    ptree = STRtree(pts)
    hormones = ('INS', 'GCG', 'SST')
    hormone_level = np.array([audit['levels']['hormone_scaled'][h] for h in hormones])
    hi = [MARKERS.index(h) for h in hormones]

    # signed distance to the nearest object boundary (negative inside one)
    bounds = shapely.boundary(geoms)
    nearest = STRtree(bounds).nearest(pts)
    signed = shapely.distance(pts, bounds[nearest]) * PX
    signed[core_owner >= 0] *= -1

    obs, means, counts, morph, comp, targets, markers_rows, feats, lookup, unpaired = ([] for _ in range(10))
    display_name = np.full(len(cells), '', dtype=object)
    display_region = np.full(len(cells), 'tissue', dtype=object)
    gfeatures = []
    cell_dir = out / 'cells'
    for i, (name, ocls, core) in enumerate(objs):
        local = shapely.union_all(geoms[otree.query(core.buffer(PERI_UM / PX), predicate='intersects')])
        band = core.buffer(PERI_UM / PX).difference(local).intersection(tissue)
        c = np.flatnonzero(core_owner == i)
        p = np.setdiff1d(ptree.query(band, predicate='covers'), np.flatnonzero(core_owner >= 0)) if not band.is_empty \
            else np.array([], dtype=int)
        display_name[c], display_region[c] = name, 'core'
        free = p[(display_region[p] == 'tissue') | ((display_region[p] == 'peri') & (nearest[p] == i))]
        display_name[free], display_region[free] = name, 'peri'
        identity = f'{cid}_{name}'
        gfeatures.append({'type': 'Feature', 'geometry': mapping(core),
                          'properties': {'name': name, 'object_class': ocls, 'classification': {'name': 'Islet'}}})
        gfeatures.append({'type': 'Feature', 'geometry': mapping(core.union(band)),
                          'properties': {'name': name + '_exp20um', 'object_class': ocls,
                                         'classification': {'name': 'IsletExpanded'}}})
        # hormone positivity (D4): each scaled hormone at or above the image's hormone level on that hormone
        # (classify.hormone_scaled_levels: the level times the hormone's own tissue p99, on the cohort scale)
        h = X[c][:, hi]
        hv = np.isfinite(h).all(axis=1)
        pos = h >= hormone_level[None, :]
        npos = pos.sum(axis=1)
        rec = dict(meta, source_file='uniform_20donor', islet_id=name, name=name, object_class=ocls,
                   cells_total=len(c), Multi_Pos=int(np.sum(hv & (npos > 1))), Triple_Neg=int(np.sum(hv & (npos == 0))))
        for j, pre in enumerate(['Ins', 'Glu', 'Stt']):
            rec[pre + '_single'] = int(np.sum(hv & pos[:, j] & (npos == 1)))
            rec[pre + '_any'] = int(np.sum(hv & pos[:, j]))
            rec[pre + '_thr'] = float(hormone_level[j])
        comp.append(rec)
        for rtype, g, ids in [('islet_core', core, c), ('islet_band', band, p),
                              ('islet_union', core.union(band), np.union1d(c, p))]:
            region = name + {'islet_core': '', 'islet_band': '_exp20um', 'islet_union': '_union'}[rtype]
            area_region = g.area * PX ** 2
            for cls in STRUCTURE_CLASSES:
                s = structures[cls]
                area, count = 0.0, 0
                if s is not None and not g.is_empty:
                    near = s[1].query(g, predicate='intersects')
                    if len(near):
                        inter = shapely.area(shapely.intersection(s[0][near], g))
                        area, count = float(inter.sum() * PX ** 2), int((inter > 0).sum())
                targets.append(dict(meta, source_file='uniform_20donor', region=region, type=rtype, **{'class': cls},
                                    area_um2=area, region_um2=area_region,
                                    area_density=area / area_region if area_region else np.nan, count=count,
                                    _islet_name=name, _region_type=rtype))
            for m in REGION_MARKERS:
                v = X[ids, MARKERS.index(m)]
                v = v[np.isfinite(v)]
                n = len(v)
                k = int(np.sum(v >= CT.IMMUNE_REQUIRE_MIN))
                markers_rows.append(dict(meta, source_file='uniform_20donor', region=region, region_type=rtype,
                                         marker=m, n_cells=n, mean=float(v.mean()) if n else np.nan,
                                         SD=float(v.std(ddof=1)) if n > 1 else np.nan, k=np.nan,
                                         threshold=CT.IMMUNE_REQUIRE_MIN, pos_count=k, pos_frac=k / n if n else np.nan,
                                         _islet_name=name, _region_type=rtype))
        cen = core.centroid
        lookup.append(dict(case_id=int(cid), imageid=int(cid), islet_key=name, combined_islet_id=identity,
                           object_class=ocls, centroid_x_um=cen.x * PX, centroid_y_um=cen.y * PX,
                           area_um2=core.area * PX ** 2))
        if not len(c) or not len(p):
            unpaired.append(dict(case_id=cid, islet_id=identity, object_class=ocls, core_cells=len(c), peri_cells=len(p),
                                 reason='no analysed cell in the core' if not len(c) else 'no analysed cell in the band'))
            continue
        cc = np.array([np.sum(labels[c] == l) for l in LABELS])
        pc = np.array([np.sum(labels[p] == l) for l in LABELS])
        o = dict(islet_id=identity, base_islet_id=name, combined_islet_id=identity, match_key=f'{name}_{cid}',
                 object_class=ocls, region='islet', donor_status=meta['Donor Status'], age=meta['Age'],
                 gender=meta['Gender'], imageid=int(cid), case_id_str=cid, total_cells=len(c),
                 GADA=meta['AAb_GADA'], ZnT8A=meta['AAb_ZnT8A'], IA2A=meta['AAb_IA2A'], mIAA=meta['AAb_mIAA'],
                 total_cells_core=len(c), total_cells_peri=len(p),
                 excluded_cells_core=int(excluded_in.get(name, 0)),
                 core_area_um2=core.area * PX ** 2)
        for l, v in zip(LABELS, cc):
            o['count_' + l] = int(v)
            o['prop_' + l] = v / len(c)
        # neighbourhood (the frozen metrics' formulas)
        cl, pl = labels[c], labels[p]
        ci, pi = int(np.isin(cl, IMMUNE).sum()), int(np.isin(pl, IMMUNE).sum())
        o.update(immune_count_core=ci, immune_count_peri=pi, immune_frac_core=ci / len(c), immune_frac_peri=pi / len(p),
                 immune_ratio=(pi / len(p)) / (ci / len(c)) if ci else np.nan,
                 cd8_to_macro_ratio=np.sum(pl == 'CD8a Tcell') / np.sum(pl == 'Macrophage')
                 if np.any(pl == 'Macrophage') else np.nan,
                 tcell_density_peri=100 * np.isin(pl, ['CD8a Tcell', 'CD4 Tcell', 'T cell']).sum() / len(p))
        dist = np.linalg.norm(xy_um[p] - xy_um[c].mean(axis=0), axis=1)
        idist = dist[np.isin(pl, IMMUNE)]
        o['min_dist_immune_mean'] = float(idist.min()) if len(idist) else np.nan
        for l in LABELS:
            k, n = int(np.sum(pl == l)), len(p)
            o['peri_count_' + safe(l)] = k
            o['peri_prop_' + safe(l)] = k / n
            exp = frequencies[l] * n
            o['enrich_z_' + safe(l)] = (k - exp) / np.sqrt(exp) if exp > 0 else np.nan
            o['min_dist_' + safe(l)] = float(dist[pl == l].min()) if k else np.nan
        obs.append(o)
        means.append(np.nanmean(X[c], axis=0))
        counts.append(cc)
        morph.append(cells.iloc[c][list(MORPHOLOGY.values())].apply(pd.to_numeric, errors='coerce').mean().to_numpy(float))
        endo = int(sum(cc[LABELS.index(l)] for l in ENDOCRINE))
        denom = max(endo, 1)
        f = dict(islet_id=identity, case_id=int(cid), donor_status=meta['Donor Status'], object_class=ocls,
                 total_cells_core=len(c), total_cells_peri=len(p), n_endocrine_core=endo)
        for l, cv, pv in zip(LABELS, cc, pc):
            f.update({f'core_{l}': np.log1p(cv / denom), f'peri_{l}': np.log1p(pv / denom),
                      f'n_core_{l}': int(cv), f'n_peri_{l}': int(pv)})
        feats.append(f)
        ids = np.r_[c, p]
        d = pd.DataFrame({'X_centroid': xy_um[ids, 0], 'Y_centroid': xy_um[ids, 1], 'phenotype': labels[ids],
                          'cell_region': np.r_[np.repeat('core', len(c)), np.repeat('peri', len(p))],
                          'Cell Area': cells['Cell: Area µm^2'].to_numpy()[ids],
                          'Nucleus Area': cells['Nucleus: Area µm^2'].to_numpy()[ids]})
        d = pd.concat([d, pd.DataFrame(X[ids], columns=MARKERS)], axis=1)
        d.to_csv(cell_dir / f'{identity}.csv', index=False)

    if not obs:
        raise ValueError(f'{tag}: no object with analysed core and band cells')
    for cls in STRUCTURE_CLASSES:
        for f in features(run / 'structures' / cls / f'{cls}__{tag}.geojson'):
            if f.get('geometry'):
                gfeatures.append({'type': 'Feature', 'geometry': f['geometry'], 'properties': {'classification': {'name': cls}}})
    gfeatures.append({'type': 'Feature', 'geometry': mapping(tissue), 'properties': {'classification': {'name': 'Tissue'}}})
    with gzip.open(out / 'gson' / f'{int(cid):04d}.geojson.gz', 'wt') as fh:
        json.dump({'type': 'FeatureCollection', 'features': gfeatures}, fh)

    tissue_tab = pd.DataFrame({'X_centroid': xy_um[:, 0], 'Y_centroid': xy_um[:, 1], 'phenotype': labels,
                               'cell_region': display_region, 'islet_name': display_name})
    tissue_tab.to_csv(out / 'donors' / f'{cid}.csv', index=False)
    distance = pd.DataFrame({'case_id': int(cid), 'donor_status': meta['Donor Status'], 'cell_region': display_region,
                             'phenotype': labels, 'signed_dist': signed.astype(np.float32)})
    for m in POS_MARKERS:
        distance[m + '_pos'] = (X[:, MARKERS.index(m)] >= CT.IMMUNE_REQUIRE_MIN).astype(np.int8)
    kde = distance.loc[display_region != 'tissue', ['case_id', 'donor_status', 'cell_region', 'phenotype',
                                                    'signed_dist']].copy()
    kde['is_immune'] = kde['phenotype'].isin(IMMUNE)
    kde['phenotype'] = kde['phenotype'].map(safe)
    membership = pd.DataFrame({'object_id': cells['object_id'], 'case_id': cid, 'phenotype': labels,
                               'cell_class': cells['cell_class'], 'refined_class': cells['refined_class'],
                               'immune_definition': cells['immune_definition'], 'cell_region': display_region,
                               'islet_name': display_name, 'x_um': xy_um[:, 0], 'y_um': xy_um[:, 1]})
    membership.to_parquet(out / 'cell_types' / f'{cid}.parquet', index=False)
    census.update(objects=len(objs), objects_paired=len(obs), seconds=round(time.time() - t0))
    a = ad.AnnData(np.vstack(means).astype(np.float32), obs=pd.DataFrame(obs).set_index('islet_id', drop=False),
                   var=pd.DataFrame(index=MARKERS))
    a.obsm['phenotype_counts'] = np.asarray(counts, dtype=np.float32)
    a.obsm['phenotype_proportions'] = a.obsm['phenotype_counts'] / a.obs.total_cells.to_numpy()[:, None]
    a.obsm['morphology'] = np.asarray(morph, dtype=np.float32)
    print(f'{tag}: {len(obs)} of {len(objs)} objects, {len(cells):,} analysed cells, {census["seconds"]}s', flush=True)
    return dict(adata=a, comp=comp, targets=targets, markers=markers_rows, features=feats, lookup=lookup,
                unpaired=unpaired, distance=distance, kde=kde, census=census)


def groovy(adata, kind, rows, columns):
    df = pd.DataFrame(rows).reindex(columns=columns)
    adata.uns[f'groovy_{kind}_columns'] = np.asarray(columns, dtype=object)
    adata.uns[f'groovy_{kind}_n_rows'] = len(df)
    for col in columns:
        s = df[col]
        if s.dtype == object:
            adata.uns[f'groovy_{kind}_{col}'] = s.fillna('').astype(str).to_numpy(dtype=object)
        else:
            adata.uns[f'groovy_{kind}_{col}'] = pd.to_numeric(s).to_numpy(dtype=float)


COMP_COLS = ['Case ID', 'source_file', 'islet_id', 'name', 'object_class', 'cells_total', 'Ins_single', 'Glu_single',
             'Stt_single', 'Multi_Pos', 'Triple_Neg', 'Ins_any', 'Glu_any', 'Stt_any', 'Ins_thr', 'Glu_thr', 'Stt_thr',
             'Donor Status', 'Age', 'Gender', 'AAb_GADA', 'AAb_ZnT8A', 'AAb_IA2A', 'AAb_mIAA']
TARGET_COLS = ['Case ID', 'source_file', 'region', 'type', 'class', 'area_um2', 'region_um2', 'area_density', 'count',
               'Donor Status', 'Age', 'Gender', 'AAb_GADA', 'AAb_ZnT8A', 'AAb_IA2A', 'AAb_mIAA', '_islet_name',
               '_region_type']
MARKER_COLS = ['Case ID', 'source_file', 'region', 'region_type', 'marker', 'n_cells', 'mean', 'SD', 'k', 'threshold',
               'pos_count', 'pos_frac', 'Donor Status', 'Age', 'Gender', 'AAb_GADA', 'AAb_ZnT8A', 'AAb_IA2A',
               'AAb_mIAA', '_islet_name', '_region_type']


def cohort_frequencies(run: Path, rows: list[dict]) -> dict:
    """Tissue-wide label frequencies over every analysed cell of the cohort (enrich_z's expectation)."""
    typ = CT.typing_dir(run)
    total = pd.Series(0, index=LABELS, dtype=float)
    for r in rows:
        t = r['image_name']
        calls = pd.read_parquet(typ / 'calls' / f'{t}.parquet')
        imm = pd.read_parquet(typ / 'immune' / f'{t}.parquet', columns=['object_id', 'immune_class', 'immune_subtype'])
        lab = pd.Series(app_labels(calls, imm))
        total = total.add(lab[lab != 'Excluded'].value_counts().reindex(LABELS, fill_value=0), fill_value=0)
    return (total / total.sum()).to_dict()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--donor-source', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--report-root', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    out, reports = a.output.resolve(), a.report_root.resolve()
    if out.exists() or reports.exists():
        raise FileExistsError('use fresh output and report roots')
    rows = CT.cohort_rows()
    assert_no_excluded([r['case_id'] for r in rows], 'Cohort manifest')
    donors = ad.read_h5ad(a.donor_source / 'islet_explorer.h5ad', backed='r').obs.drop_duplicates('imageid')[
        ['imageid', 'age', 'gender', 'donor_status', 'GADA', 'ZnT8A', 'IA2A', 'mIAA']].copy()
    if set(donors.imageid.astype(int)) != {int(r['case_id']) for r in rows}:
        raise ValueError('donor source does not cover exactly the cohort')
    for d in ('cells', 'donors', 'gson', 'cell_types', 'images'):
        (out / d).mkdir(parents=True)
    reports.mkdir(parents=True)
    freqs = cohort_frequencies(a.run, rows)
    results = {}
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {r['image_name']: pool.submit(build_image, a.run, out, r['image_name'], r['case_id'], donors, freqs)
                for r in rows}
        for t, f in futs.items():
            results[t] = f.result()
    order = [r['image_name'] for r in rows]
    combined = ad.concat([results[t]['adata'] for t in order], join='outer', merge='same')
    combined.obs.index.name = None
    for col in combined.obs:
        if combined.obs[col].dtype == object:
            combined.obs[col] = pd.Categorical(combined.obs[col])
    combined.uns['phenotype_names'] = np.asarray(LABELS, dtype=object)
    combined.uns['morphology_names'] = np.asarray(list(MORPHOLOGY), dtype=object)
    combined.uns['case_to_imageid'] = {r['case_id']: int(r['case_id']) for r in rows}
    combined.uns['reference_frequencies'] = freqs
    sys.path.insert(0, str(REPO / 'scripts'))
    from build_protein_embedding import fit_embedding
    # A cytoplasm-measured marker is undefined for a cell whose nucleus fills it (REDSEA writes NaN, as Senior's
    # does), so an object made only of such cells has no mean for it. X keeps the NaN; the embedding is fitted on a
    # copy in which a missing object mean is the marker's cohort mean (a zero after standardising), and the
    # number of objects and values filled is recorded.
    measured = np.asarray(combined.X, dtype=float)
    missing = ~np.isfinite(measured)
    filled = np.where(missing, np.nanmean(measured, axis=0)[None, :], measured)
    fit_input = combined.copy()
    fit_input.X = filled.astype(np.float32)
    imputation = {'objects_with_a_missing_mean': int(missing.any(axis=1).sum()), 'values_filled': int(missing.sum()),
                  'by_marker': {m: int(n) for m, n in zip(MARKERS, missing.sum(axis=0)) if n}}
    (reports / 'embedding_imputation.json').write_text(json.dumps(imputation, indent=2) + '\n')
    uncorrected, udiag = fit_embedding(fit_input, donor_correction='none')
    uncorrected.write_h5ad(reports / 'uncorrected_protein_pca.h5ad', compression='gzip')
    (reports / 'uncorrected_embedding.json').write_text(json.dumps(udiag, indent=2, default=float) + '\n')
    del uncorrected
    combined, diag = fit_embedding(fit_input, donor_correction='harmony')
    combined.X = measured.astype(np.float32)
    combined.uns['protein_embedding']['missing_means_filled_for_fit'] = imputation['values_filled']
    (reports / 'embedding.json').write_text(json.dumps(diag, indent=2, default=float) + '\n')
    trajectory = combined.copy()
    for kind, cols in [('composition', COMP_COLS), ('targets', TARGET_COLS), ('markers', MARKER_COLS)]:
        groovy(combined, kind, [r for t in order for r in results[t][{'composition': 'comp'}.get(kind, kind)]], cols)
    head = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    combined.uns['data_provenance'] = {'builder': 'scripts/uniform/assemble_uniform.py', 'git_commit': head,
                                      'run': str(a.run.resolve()), 'cell_typing': 'Islet-Explorer-Senior port',
                                      'donor_source': str(a.donor_source.resolve())}
    combined.write_h5ad(out / 'islet_explorer.h5ad', compression='gzip')
    trajectory.write_h5ad(out / 'adata_ins_root.h5ad', compression='gzip')
    pd.DataFrame([r for t in order for r in results[t]['features']]).to_csv(out / 'barlow_islet_features.csv', index=False)
    pd.DataFrame([r for t in order for r in results[t]['lookup']]).to_csv(out / 'islet_spatial_lookup.csv', index=False)
    pd.concat([results[t]['distance'] for t in order]).to_csv(out / 'cell_distance.csv', index=False)
    pd.concat([results[t]['kde'] for t in order]).to_csv(out / 'signed_distance_kde.csv', index=False)
    nbr = [c for c in combined.obs.columns if c.startswith(('peri_', 'immune_', 'enrich_z_', 'min_dist_', 'cd8_',
                                                            'tcell_', 'total_cells'))]
    combined.obs[['islet_id', 'imageid', 'object_class'] + nbr].to_csv(out / 'neighborhood_metrics.csv', index=False)
    pd.DataFrame([r for t in order for r in results[t]['unpaired']]).to_csv(reports / 'unpaired_objects.csv', index=False)
    census = pd.DataFrame([results[t]['census'] for t in order])
    census.to_csv(out / 'census.csv', index=False)
    md = pd.read_csv(a.donor_source / 'donor_metadata.csv', dtype={'case_id': str})
    assert_no_excluded(md.case_id, 'donor metadata')
    md.to_csv(out / 'donor_metadata.csv', index=False)
    images = []
    for img in sorted((a.donor_source / 'images').iterdir()):
        target = img.resolve(strict=True)
        (out / 'images' / img.name).symlink_to(target)
        images.append({'image': img.name, 'source': str(target), 'bytes': target.stat().st_size})
    (out / 'image_sources.json').write_text(json.dumps(images, indent=2) + '\n')
    (out / 'protein_phenotype_rules.csv').write_text((HERE / 'app_phenotype_rules.csv').read_text())
    rel = {'status': 'assembled_pending_owner_review', 'donors': len(rows), 'objects': int(combined.n_obs),
           'by_object_class': {str(k): int(v) for k, v in combined.obs.object_class.value_counts().items()},
           'labels': LABELS, 'donor_correction': 'harmony', 'scvi_required': False, 'git_commit': head,
           'run': str(a.run.resolve()), 'scale_sha256': sha256(CT.typing_dir(a.run) / 'scale.parquet'),
           'shifts_sha256': sha256(CT.typing_dir(a.run) / 'shifts.json'),
           'cells_analysed': int(census.analysed.sum()), 'cells_excluded': int(census.excluded.sum()),
           'cells_not_analysis_eligible': int(census.not_analysis_eligible.sum())}
    (out / 'uniform_release.json').write_text(json.dumps(rel, indent=2) + '\n')
    print(json.dumps(rel, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
