#!/usr/bin/env python3
"""assemble_uniform.py -- build a complete app data root for all 20 donors from the uniform run.

Every value is regenerated from the uniform workflow (`$OUT`, see README.md): the accepted endocrine objects
(`$OUT/islets/`), the structure objects (`$OUT/structures/`), the cells and their Islet-Explorer-Senior typing
(`$OUT/typing/`). Nothing is appended to the frozen 15-donor release. Owner decisions applied:

  A5  Islets and Clusters are analysis rows -- Islet_N, Cluster_N -- with `object_class` (Islet /
      EndocrineCluster) so the app can filter them. Singles (Single_N, one endocrine cell) are counted and drawn
      but are not analysis objects (owner, 2026-10-18: "Singles should not be included in any islet-level
      analysis except for quantifying or showing how many there are"): no row, no core or band, no signed
      distance to them, and their cell is an ordinary tissue cell, which may lie in an islet's or cluster's band.
  A6  nerve, blood-vessel and lymphatic cells come only from the structure objects (the typing).
  D3  cell types are the final typing's labels (typing/identity, identity.py; owner 2026-10-29: "replace the old
      chain with what we have so far"); until v17 Senior's classes from the old chain (classify.py), shown under
      the app's labels.
  D4  positivity, hormones included, is each marker above its image's own level on its one measurement
      (typing/levels.csv; identity.py's pos__ columns). KRT14, CD44, NaKATPase and BActin, which the typing does not
      read, take display-only levels by E-cadherin's method (typing/display_levels.csv, display_levels.py; owner
      2026-10-29 "just use the ecad methods"). ColIV becomes structure objects (owner 2026-10-29), so it has no cell
      positivity flag; its value scale's floor is the image's empty-slide null.
  D5  marker values are the typing's: despeckled REDSEA on each marker's measurement (whole cell, or membrane where
      the level selection picked it), on the image's own scale (x - level) / (p99.9 - level) clipped to [0, 1]
      (identity.py's hormone typing scale); nothing pooled (owner 2026-10-28). Until v17: UniFORM-aligned values on
      a cohort floor/top scale.

Regions: an object's core is the cells whose centroid lies in it; its peri band is the 20 um ring around it,
inside Tissue, minus every endocrine object; a cell may sit in several bands. Cells ruled out of analysis (no
nucleus) and Excluded cells are counted in `census.csv` and kept out of every composition, neighbourhood and display
table. Label groups (agent choice 2026-10-29, recorded in uniform_release.json): endocrine = the Endocrine lineage's
labels plus unnamed cells of the objects that are their endocrine cells (EndocrineCluster, EndocrineSingle); immune =
the labels the Hematopoietic compartment names (T, B, macrophage profiles, Hematopoietic CD45+), not Hematopoietic
CD45-, conflicts or doublets.

The UMAP and Leiden clusters are the endocrine composition of the objects with ten or more core cells
(composition_embedding.py, owner 2026-10-19); the Harmony protein embedding stays in the H5AD for the record.

Donor demographics and autoantibodies, and the viewer image links, are carried from the last served root
(`--donor-source`, v3), which recorded them from nPOD; nothing else is read from it.

    python scripts/uniform/assemble_uniform.py --run data/development/uniform-20donor-20261015 \\
        --donor-source data/development/protein-expanded-20261015-v3 \\
        --output data/development/uniform-20donor-20261020-v9 --report-root tmp/uniform-20donor-20261020-v9
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
import identity as ID  # noqa: E402
import levels as LV  # noqa: E402
import registry as RG  # noqa: E402
from excluded_donors import assert_no_excluded  # noqa: E402

PX = CT.PIXEL_UM
PERI_UM = 20.0
MARKERS = ['CD31', 'CD8a', 'CD3e', 'SMA', 'CD4', 'CD34', 'HLADR', 'PDPN', 'panCK', 'ECAD', 'CD163', 'SST', 'ColIV',
           'VIM', 'CD20', 'LGALS3', 'B3TUBB', 'GCG', 'KRT14', 'GAP43', 'CD35', 'CHGA', 'PGP9.5', 'INS', 'CD44',
           'CD45', 'NaKATPase', 'BCatenin', 'CD68', 'BActin', 'CK19']          # the 31-marker panel, app order
# The labels are the final typing's (typing/identity), read from the run (label_groups); until v17 a fixed list of the
# old chain's classes.
ENDOCRINE_STARTS = ('EndocrineCluster', 'EndocrineSingle')     # an unnamed cell of an object that is its endocrine cells
IDENTITY_COLUMNS = ['object_id', 'label', 'excluded_reason', 'compartment', 'lineage', 'status', 'named_by', 'tags']
# The `<marker>_pos` flags of cell_distance.csv cover every panel marker (owner, 2026-10-20: the distance profile
# "needs to have all markers that are common to all donors"; all 31 are measured in every image, and build_image
# stops if one is not). Until v7 they covered 18.
# ColIV is left out since 2026-10-29: it becomes structure objects (owner), not a cell call.
STRUCTURE_MARKERS = ('ColIV',)
POS_MARKERS = [m for m in MARKERS if m not in STRUCTURE_MARKERS]
REGION_MARKERS = ['CD3e', 'CD163', 'CD68', 'LGALS3']
MORPHOLOGY = {'Nucleus Area': 'Nucleus: Area µm^2', 'Nucleus Length': 'Nucleus: Length µm',
              'Nucleus Circularity': 'Nucleus: Circularity', 'Nucleus Solidity': 'Nucleus: Solidity',
              'Nucleus: Max diameter µm': 'Nucleus: Max diameter µm',
              'Nucleus: Min diameter µm': 'Nucleus: Min diameter µm', 'Cell Area': 'Cell: Area µm^2',
              'Cell Length': 'Cell: Length µm', 'Cell Circularity': 'Cell: Circularity',
              'Cell Solidity': 'Cell: Solidity', 'Cell: Max diameter µm': 'Cell: Max diameter µm',
              'Cell: Min diameter µm': 'Cell: Min diameter µm'}
STRUCTURE_CLASSES = ('Nerve', 'Capillary', 'Lymphatic')
ANALYSIS_OBJECTS = ('Islet', 'EndocrineCluster')
SINGLE = 'EndocrineSingle'


def safe(label: str) -> str:
    return label.replace(' ', '_').replace('+', 'plus').replace('/', '_')


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
def label_groups(run: Path, rows: list[dict]) -> dict:
    """The cohort's labels (labels-file order first, then the rest sorted), their tissue-wide frequencies over every
    analysed cell (enrich_z's expectation) and the endocrine, immune and macrophage groups."""
    reg = RG.load()
    lab = ID.labels_table(reg)
    total, comp = pd.Series(dtype=float), {}
    for r in rows:
        d = pd.read_parquet(run / 'typing' / 'identity' / f"{r['image_name']}.parquet",
                            columns=['label', 'compartment', 'lineage', 'named_by', 'status'])
        d = d[d['label'] != 'Excluded']
        total = total.add(d['label'].value_counts(), fill_value=0)
        named = (d['status'] == '') & (d['named_by'] != '') & ~d['label'].str.startswith(('Conflict ', 'Doublet'))
        for (lbl, c, lin), _ in d[named].groupby(['label', 'compartment', 'lineage']).size().items():
            comp.setdefault(lbl, set()).add((c, lin))
    first = [x for x in dict.fromkeys(lab['label']) if x in total.index]
    labels = first + sorted(x for x in total.index if x not in first)
    pan_neg = set(lab.loc[lab['kind'] == 'pan_negative', 'label'])
    endocrine = [x for x in labels if any(lin == 'Endocrine' for _, lin in comp.get(x, ())) or x in ENDOCRINE_STARTS]
    immune = [x for x in labels if x not in pan_neg and comp.get(x) and all(c == 'Hematopoietic' for c, _ in comp[x])]
    macrophages = [x for x in immune if any(lin == 'Myeloid' for _, lin in comp[x])]
    return {'labels': labels, 'frequencies': (total / total.sum()).reindex(labels).to_dict(),
            'endocrine': endocrine, 'immune': immune, 'macrophages': macrophages}


# ------------------------------------------------------------------------------------------- one image
def apply_donor_metadata(donors: pd.DataFrame, md: pd.DataFrame) -> pd.DataFrame:
    """Every donor's age from the nPOD donor table; stop if its status, sex or autoantibodies disagree."""
    status = {'No Diabetes': 'ND', 'Autoantibody Positive': 'Aab+', 'Type 1 Diabetes': 'T1D'}
    md = md.set_index(md.case_id.astype(int))
    donors = donors.copy()
    ids = donors.imageid.astype(int)
    if set(ids) != set(md.index):
        raise ValueError('donor metadata does not cover exactly the cohort')
    for i, cid in zip(donors.index, ids):
        m, d = md.loc[cid], donors.loc[i]
        mine = (str(d.donor_status), str(d.gender)[0].upper(),
                *(bool(float(d[k])) for k in ('GADA', 'IA2A', 'ZnT8A', 'mIAA')))
        theirs = (status[m.donor_type], str(m.sex)[0].upper(),
                  *(str(m['AAb_' + k]).lower() == 'true' for k in ('GADA', 'IA2A', 'ZnT8A', 'mIAA')))
        if mine != theirs:
            raise ValueError(f'donor {cid}: donor source {mine} disagrees with donor metadata {theirs}')
    donors['age'] = [float(md.loc[cid, 'age_years_npod']) for cid in ids]
    return donors


def donor_fields(donors: pd.DataFrame, case: str) -> dict:
    r = donors.loc[donors.imageid.astype(int) == int(case)].iloc[0]
    return {'Case ID': f'{int(case):04d}', 'Donor Status': r.donor_status, 'Age': float(r.age), 'Gender': str(r.gender),
            'AAb_GADA': float(r.GADA), 'AAb_IA2A': float(r.IA2A), 'AAb_ZnT8A': float(r.ZnT8A), 'AAb_mIAA': float(r.mIAA)}


def build_image(run: Path, out: Path, tag: str, case: str, donors: pd.DataFrame, groups: dict) -> dict:
    t0 = time.time()
    typ = CT.typing_dir(run)
    meta = donor_fields(donors, case)
    cid = str(int(case))
    LABELS, ENDOCRINE, IMMUNE, MACROPHAGES = groups['labels'], groups['endocrine'], groups['immune'], groups['macrophages']
    frequencies = groups['frequencies']
    reg = RG.load()
    ident = pd.read_parquet(typ / 'identity' / f'{tag}.parquet')
    context = pd.read_parquet(typ / 'context' / f'{tag}.parquet', columns=['object_id', 'qc_analysis_eligible', 'islet_num'])
    lv = pd.read_csv(typ / 'levels.csv')
    lv = lv[lv['image'] == tag].set_index('marker')
    shown = pd.read_csv(typ / 'display_levels.csv')
    shown = shown[shown['image'] == tag].set_index('marker')
    nulls = LV.empty_nulls(run)
    raw = pd.read_parquet(run / 'cells' / f'{tag}.parquet',
                          columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'] + list(MORPHOLOGY.values()))
    census = {'image': tag, 'case_id': cid, 'objects_detected': int(len(context)),
              'not_analysis_eligible': int((~context['qc_analysis_eligible']).sum())}

    cells = ident[IDENTITY_COLUMNS].merge(context[['object_id', 'islet_num']], on='object_id', how='left',
                                          validate='one_to_one')
    cells['phenotype'] = cells['label']
    census['excluded'] = int((cells['phenotype'] == 'Excluded').sum())
    # Singles are counted from the typed cells before Excluded cells are dropped (one cell per Single).
    single_labels = cells.loc[cells['islet_num'].astype(str).str.startswith('Single_'), 'phenotype'].value_counts()
    excluded_in = cells.loc[cells['phenotype'] == 'Excluded', 'islet_num'].value_counts()
    census['excluded_reasons'] = {str(k): int(v) for k, v in
                                  cells.loc[cells.phenotype == 'Excluded', 'excluded_reason'].value_counts().items()}
    cells = cells[cells['phenotype'] != 'Excluded'].reset_index(drop=True)
    cells = cells.merge(raw.rename(columns={'Object ID': 'object_id'}), on='object_id', how='left', validate='one_to_one')
    # D4/D5: each marker on its measurement, above / scaled to the image's own level (display-only level for the
    # markers the typing does not read); the scale's top is the 99.9th percentile over the image's analysed cells
    meas = {m: (lv.loc[m, 'measurement'] if m in lv.index else shown.loc[m, 'measurement'] if m in shown.index
                else 'whole_cell') for m in MARKERS}
    rs = pd.read_parquet(typ / 'redsea' / f'{tag}.parquet', columns=['object_id', 'qc_analysis_eligible']
                         + sorted({reg.column(m, meas[m]) for m in MARKERS}))
    rs = rs[rs['qc_analysis_eligible']].reset_index(drop=True)
    level_log = {m: float(lv.loc[m, 'level_log']) if m in lv.index else float(shown.loc[m, 'level_log'])
                 if m in shown.index else float(np.log1p(nulls[(tag, m)])) for m in MARKERS}
    missing = [m for m in MARKERS if not np.isfinite(level_log[m])]
    if missing:
        raise ValueError(f'{tag}: no level for {missing}')
    levels = pd.DataFrame({'level_raw': {m: float(np.expm1(v)) for m, v in level_log.items()},
                           'display_only': {m: m not in lv.index for m in MARKERS}})
    analysed = rs[[reg.column(m, meas[m]) for m in MARKERS]].clip(lower=0).to_numpy(float)
    Lr = levels.loc[MARKERS, 'level_raw'].to_numpy(float)
    top = np.nanpercentile(analysed, reg.const('HORMONE_SCALE_TOP_QUANTILE'), axis=0)
    # fewer than 0.1 % of the cells above the level (CD20 in an image with few B cells): the top is the highest value;
    # no cell above the level: every value is 0 on the scale
    top = np.where(top > Lr, top, np.nanmax(analysed, axis=0))
    vals = cells[['object_id']].merge(rs, on='object_id', how='left', validate='one_to_one')
    V = np.column_stack([vals[reg.column(m, meas[m])].clip(lower=0).to_numpy(float) for m in MARKERS])
    with np.errstate(invalid='ignore', divide='ignore'):
        X = np.clip((V - Lr) / np.where(top > Lr, top - Lr, np.inf), 0, 1)   # cells x 31, [0, 1], NaN where unmeasured
    E = np.where(np.isfinite(V), np.log1p(V) - np.log1p(Lr), np.nan)        # > 0 positive (above the image's level)
    pos = ident.set_index('object_id').reindex(cells['object_id'])
    for m in MARKERS:                                                       # the typing's own calls where it has them
        if f'pos__{m}' in pos.columns:
            agree = (E[:, MARKERS.index(m)] > 0) == pos[f'pos__{m}'].to_numpy(bool)
            if not agree[np.isfinite(E[:, MARKERS.index(m)])].all():
                raise ValueError(f'{tag}: {m} positivity differs from the typing (identity pos__{m})')
    labels = cells['phenotype'].to_numpy(dtype=object)
    xy_um = cells[['Centroid X µm', 'Centroid Y µm']].to_numpy(float)
    pts = shapely.points(xy_um / PX)
    census['analysed'] = int(len(cells))

    # objects, tissue, structures (pixel coordinates)
    objs, singles = [], []
    for f in features(run / 'islets' / f'{tag}.geojson'):
        p = f.get('properties') or {}
        c = (p.get('classification') or {}).get('name')
        if c in ANALYSIS_OBJECTS and f.get('geometry'):
            objs.append((p['name'], c, valid(shape(f['geometry']))))
        elif c == SINGLE and f.get('geometry'):
            singles.append((p['name'], valid(shape(f['geometry']))))
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
    in_obj = cells['islet_num'].to_numpy(dtype=object)
    index = {n: i for i, n in enumerate(names)}
    core_owner = np.array([index.get(k, -1) if isinstance(k, str) else -1 for k in in_obj])
    ptree = STRtree(pts)
    hormones = ('INS', 'GCG', 'SST')
    hormone_level = levels.loc[list(hormones), 'level_raw'].to_numpy(float)
    hi = [MARKERS.index(h) for h in hormones]

    # signed distance to the nearest object boundary (negative inside one)
    bounds = shapely.boundary(geoms)
    nearest = STRtree(bounds).nearest(pts)
    signed = shapely.distance(pts, bounds[nearest]) * PX
    signed[core_owner >= 0] *= -1
    # The same to the nearest Islet only (owner, 2026-10-20): Clusters are small and scattered through the tissue,
    # so 90 % of the cells beyond 20 um of any object are measured to a Cluster, which blurs the profile around
    # Islets. Negative inside an Islet; a cell inside a Cluster is outside every Islet, so positive.
    islet_ix = np.flatnonzero(np.array([o[1] == 'Islet' for o in objs], dtype=bool))
    if len(islet_ix):
        ibounds = bounds[islet_ix]
        signed_islet = shapely.distance(pts, ibounds[STRtree(ibounds).nearest(pts)]) * PX
        signed_islet[np.isin(core_owner, islet_ix)] *= -1
    else:
        signed_islet = np.full(len(pts), np.nan)

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
        # hormone positivity (D4): each hormone above the image's level of that hormone (the typing's call);
        # *_thr is that level on the measured value
        h = E[c][:, hi]
        hv = np.isfinite(h).all(axis=1)
        pos = h > 0
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
                e = E[ids, MARKERS.index(m)]
                v, e = v[np.isfinite(v)], e[np.isfinite(e)]
                n = len(v)
                k = int(np.sum(e > 0))
                markers_rows.append(dict(meta, source_file='uniform_20donor', region=region, region_type=rtype,
                                         marker=m, n_cells=n, mean=float(v.mean()) if n else np.nan,
                                         SD=float(v.std(ddof=1)) if n > 1 else np.nan, k=np.nan,
                                         threshold=float(levels.loc[m, 'level_raw']) if m in levels.index else np.nan,
                                         pos_count=k, pos_frac=k / n if n else np.nan,
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
                 cd8_to_macro_ratio=np.sum(pl == 'CD8a Tcell') / np.isin(pl, MACROPHAGES).sum()
                 if np.isin(pl, MACROPHAGES).any() else np.nan,
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
    for name, g in singles:
        gfeatures.append({'type': 'Feature', 'geometry': mapping(g),
                          'properties': {'name': name, 'object_class': SINGLE, 'classification': {'name': SINGLE}}})
    gfeatures.append({'type': 'Feature', 'geometry': mapping(tissue), 'properties': {'classification': {'name': 'Tissue'}}})
    with gzip.open(out / 'gson' / f'{int(cid):04d}.geojson.gz', 'wt') as fh:
        json.dump({'type': 'FeatureCollection', 'features': gfeatures}, fh)

    tissue_tab = pd.DataFrame({'X_centroid': xy_um[:, 0], 'Y_centroid': xy_um[:, 1], 'phenotype': labels,
                               'cell_region': display_region, 'islet_name': display_name})
    tissue_tab.to_csv(out / 'donors' / f'{cid}.csv', index=False)
    distance = pd.DataFrame({'case_id': int(cid), 'donor_status': meta['Donor Status'], 'cell_region': display_region,
                             'phenotype': labels, 'signed_dist': signed.astype(np.float32),
                             'signed_dist_islet': signed_islet.astype(np.float32)})
    unmeasured = [m for m in POS_MARKERS if np.isnan(V[:, MARKERS.index(m)]).all()]
    if unmeasured:
        raise ValueError(f'{tag}: {unmeasured} not measured, so not common to all donors; drop them from POS_MARKERS')
    for m in POS_MARKERS:
        distance[m + '_pos'] = (E[:, MARKERS.index(m)] > 0).astype(np.int8)
    kde = distance.loc[display_region != 'tissue', ['case_id', 'donor_status', 'cell_region', 'phenotype',
                                                    'signed_dist']].copy()
    kde['is_immune'] = kde['phenotype'].isin(IMMUNE)
    kde['phenotype'] = kde['phenotype'].map(safe)
    membership = pd.DataFrame({'object_id': cells['object_id'], 'case_id': cid, 'phenotype': labels,
                               **{c: cells[c].to_numpy() for c in IDENTITY_COLUMNS[3:]}, 'cell_region': display_region,
                               'islet_name': display_name, 'x_um': xy_um[:, 0], 'y_um': xy_um[:, 1]})
    membership.to_parquet(out / 'cell_types' / f'{cid}.parquet', index=False)
    in_single = cells['islet_num'].astype(str).str.startswith('Single_').to_numpy()
    census.update(objects=len(objs), objects_paired=len(obs), tissue_area_mm2=tissue.area * PX ** 2 / 1e6,
                  islets=sum(o[1] == 'Islet' for o in objs), clusters=sum(o[1] == 'EndocrineCluster' for o in objs),
                  singles=len(singles), singles_analysed=int(in_single.sum()),
                  singles_excluded=int(single_labels.get('Excluded', 0)),
                  singles_in_peri_bands=int((in_single & (display_region == 'peri')).sum()),
                  **{'singles_' + safe(l): int(single_labels.get(l, 0)) for l in ENDOCRINE},
                  seconds=round(time.time() - t0))
    a = ad.AnnData(np.vstack(means).astype(np.float32), obs=pd.DataFrame(obs).set_index('islet_id', drop=False),
                   var=pd.DataFrame(index=MARKERS))
    a.obsm['phenotype_counts'] = np.asarray(counts, dtype=np.float32)
    a.obsm['phenotype_proportions'] = a.obsm['phenotype_counts'] / a.obs.total_cells.to_numpy()[:, None]
    a.obsm['morphology'] = np.asarray(morph, dtype=np.float32)
    census['display_only_levels'] = {m: float(levels.loc[m, 'level_raw']) for m in MARKERS if levels.loc[m, 'display_only']}
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--donor-source', type=Path, required=True)
    # From v10 (owner, 2026-10-22): the donor table built from the nPOD sheet (scripts/build_donor_metadata.py
    # --uniform-cohort). It is the served donor_metadata.csv and gives every donor's age; status, sex and
    # autoantibodies must agree with the donor source's.
    ap.add_argument('--donor-metadata', type=Path)
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
    md_path = a.donor_metadata or a.donor_source / 'donor_metadata.csv'
    md = pd.read_csv(md_path, dtype={'case_id': str})
    assert_no_excluded(md.case_id, 'donor metadata')
    if a.donor_metadata:
        donors = apply_donor_metadata(donors, md)
    for d in ('cells', 'donors', 'gson', 'cell_types', 'images'):
        (out / d).mkdir(parents=True)
    reports.mkdir(parents=True)
    groups = label_groups(a.run, rows)
    LABELS, ENDOCRINE, freqs = groups['labels'], groups['endocrine'], groups['frequencies']
    results = {}
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {r['image_name']: pool.submit(build_image, a.run, out, r['image_name'], r['case_id'], donors, groups)
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
    combined.uns['endocrine_labels'] = np.asarray(ENDOCRINE, dtype=object)      # composition_embedding.endocrine_shares
    combined.uns['morphology_names'] = np.asarray(list(MORPHOLOGY), dtype=object)
    combined.uns['case_to_imageid'] = {r['case_id']: int(r['case_id']) for r in rows}
    combined.uns['reference_frequencies'] = freqs
    sys.path.insert(0, str(REPO / 'scripts'))
    from build_protein_embedding import fit_embedding
    from composition_embedding import fit_composition_embedding, MIN_CORE_CELLS
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
    # The served UMAP and Leiden: endocrine composition of the objects with ten or more core cells (owner,
    # 2026-10-19; composition_embedding.py). The protein embedding above stays in the H5AD for the record.
    combined, cdiag = fit_composition_embedding(combined)
    (reports / 'composition_embedding.json').write_text(json.dumps(cdiag, indent=2, default=float) + '\n')
    trajectory = combined.copy()
    for kind, cols in [('composition', COMP_COLS), ('targets', TARGET_COLS), ('markers', MARKER_COLS)]:
        groovy(combined, kind, [r for t in order for r in results[t][{'composition': 'comp'}.get(kind, kind)]], cols)
    head = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    combined.uns['data_provenance'] = {'builder': 'scripts/uniform/assemble_uniform.py', 'git_commit': head,
                                      'run': str(a.run.resolve()),
                                      'cell_typing': 'final typing run (typing/identity, identity.py)',
                                      'donor_source': str(a.donor_source.resolve()),
                                      'donor_metadata': str(md_path.resolve())}
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
    # Per-donor counts of every endocrine object class: the only place Singles enter the app.
    by_donor = census[['case_id', 'tissue_area_mm2', 'islets', 'clusters', 'singles', 'singles_excluded']
                      + ['singles_' + safe(l) for l in ENDOCRINE]].copy()
    by_donor.insert(1, 'donor_status', [donor_fields(donors, r['case_id'])['Donor Status'] for r in rows])
    for k in ('islets', 'clusters', 'singles'):
        by_donor[k + '_per_mm2'] = by_donor[k] / by_donor['tissue_area_mm2']
    by_donor.to_csv(out / 'endocrine_objects_by_donor.csv', index=False)
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
           'singles_counted_not_analysed': int(census.singles.sum()),
           'labels': LABELS, 'label_groups': {k: groups[k] for k in ('endocrine', 'immune', 'macrophages')},
           'label_groups_basis': 'agent choice 2026-10-29 (module docstring); owner review pending',
           'cell_typing': 'typing/identity (identity.py)', 'positivity': 'per-image level, typing/levels.csv',
           'measurements': {m: str(v) for m, v in pd.read_csv(CT.typing_dir(a.run) / 'levels.csv')
                            .groupby('marker')['measurement'].agg(lambda x: '/'.join(sorted(set(x)))).items()},
           'display_only_level_markers': sorted(census['display_only_levels'].iloc[0]),
           'display_levels': 'typing/display_levels.csv (E-cadherin method, owner 2026-10-29)',
           'structure_markers_without_cell_flag': list(STRUCTURE_MARKERS),
           'marker_values': 'despeckled REDSEA, per-image (x - level) / (p99.9 - level) clipped to [0, 1]',
           'embedding': 'endocrine_composition', 'embedding_min_core_cells': MIN_CORE_CELLS,
           'embedding_objects': cdiag['objects'], 'donor_correction': 'none',
           'protein_embedding_donor_correction': 'harmony', 'scvi_required': False, 'git_commit': head,
           'run': str(a.run.resolve()), 'levels_sha256': sha256(CT.typing_dir(a.run) / 'levels.csv'),
           'identity_sha256': {r['image_name']: sha256(CT.typing_dir(a.run) / 'identity' / f"{r['image_name']}.parquet")
                               for r in rows},
           'cells_analysed': int(census.analysed.sum()), 'cells_excluded': int(census.excluded.sum()),
           'cells_not_analysis_eligible': int(census.not_analysis_eligible.sum()),
           'donor_metadata': str(md_path.resolve()), 'donor_metadata_sha256': sha256(md_path)}
    (out / 'uniform_release.json').write_text(json.dumps(rel, indent=2) + '\n')
    print(json.dumps(rel, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
