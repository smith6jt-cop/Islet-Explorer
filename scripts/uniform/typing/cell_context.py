#!/usr/bin/env python3
"""cell_context.py -- where each cell is: islet object, structure object or ordinary tissue (Senior's context).

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/typing/context.py` (memberships, object_label),
`islet_membership.py` (membership: centroid in an endocrine object, the SMALLEST containing outline wins) and
`build_crossed_fraction.py` (crossed_frac_<structure> = area(cell outline INTERSECT that structure's objects) /
area(cell), asked of every cell).

Structures (owner, A6, 2026-10-16: "drop cell classes that indicate neural, blood vessel or lymphatic and use the
objects instead"): the Nerve, Capillary and Lymphatic annotations every image already has -- the original
project's pixel classifiers (Nerve on B3TUBB and INS, Capillary on CD31, Lymphatics on PDPN), run over each
image's Tissue in the original project (the 15) or by `scripts/cohort_expansion/prepare_qupath.py` (the five
additions). They are read from the split files `all-donor-20260920/structural/<case>/<Class>_inherited.geojson.gz`
and copied into `$OUT/structures/<Class>/<Class>__<image>.geojson` (Senior's layout), each with its source,
sha256 and the 2026-09-21 audit verdict, in `$OUT/structures/provenance.json`. The five additions' objects are
instead drawn by `structure_thresholds.py make` (thresholds fitted on the original 15) and are not overwritten
here (Lymphatic there is the original classifier's own objects, structure_thresholds.KEEP_CLASSIFIER). Senior's
"Blood vessel" is Capillary here.

Endocrine objects are the accepted step-1b objects in `$OUT/islets/<image>.geojson` (Islet_N, Cluster_N,
Single_N). Membership is checked against QuPath's own parenting (`$OUT/context/<image>.parquet` Parent) on every
image, as Senior checks its containment against an independent source.

    python scripts/uniform/typing/cell_context.py --out data/development/uniform-20donor-20261015 --images all
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from shapely import STRtree, points
from shapely.geometry import shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402

STRUCTURAL = CT.REPO / 'data/development/all-donor-20260920/structural'
CROSSED_COLUMN = {'Blood vessel': 'crossed_frac_vessel', 'Nerve': 'crossed_frac_nerve',
                  'Lymphatic': 'crossed_frac_lymphatic'}


def sha256(path: Path) -> str:
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


def import_structures(out: Path, rows: list[dict]) -> dict:
    """Copy each image's three structure classes into the run, with provenance. Idempotent."""
    dest_root = out / 'structures'
    made = dest_root / 'thresholded.json'
    thresholded = json.loads(made.read_text())['images'] if made.exists() else {}
    refreshed_path = dest_root / 'redrawn.json'
    refreshed = json.loads(refreshed_path.read_text()).get('images', {}) if refreshed_path.exists() else {}
    record = {}
    for r in rows:
        case, tag = r['case_id'], r['image_name']
        audit = json.loads((STRUCTURAL / f'{case}.json').read_text())['masks']
        record[tag] = {}
        for name in CT.STRUCTURES:
            cls = CT.OBJECT_CLASS[name]
            src = STRUCTURAL / case / f'{cls}_inherited.geojson.gz'
            dest = dest_root / cls / f'{cls}__{tag}.geojson'
            if tag in refreshed and cls in refreshed[tag].get('classes', {}):
                entry = refreshed[tag]['classes'][cls]
                if not dest.exists() or sha256(dest) != entry['sha256']:
                    raise ValueError(f'{tag} {cls}: redrawn structure file missing or stale')
                tissue = out / 'inputs' / 'tissue' / f'{tag}.geojson'
                if sha256(tissue) != refreshed[tag]['tissue_sha256']:
                    raise ValueError(f'{tag} {cls}: structures were drawn against different Tissue')
                record[tag][cls] = {'source': str(refreshed_path), **entry}
                continue
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(src, 'rb') as fi, open(dest, 'wb') as fo:
                    shutil.copyfileobj(fi, fo)
            if tag in thresholded:      # written by structure_thresholds.py make (its record says how)
                record[tag][cls] = {'source': str(made), 'sha256': sha256(dest), **thresholded[tag][cls]}
                continue
            m = audit.get(cls, {})
            record[tag][cls] = {'source': str(src), 'source_sha256': sha256(src), 'origin': m.get('source'),
                                'sha256': sha256(dest), 'features': m.get('features'),
                                'audit_2026_09_21': (m.get('review') or {}).get('decision', 'not reviewed')}
    (dest_root / 'provenance.json').write_text(json.dumps(record, indent=1))
    return record


def contain_smallest(geoms: list, keys: list, pts) -> tuple[np.ndarray, np.ndarray, int]:
    area = np.array([g.area for g in geoms])
    number = np.full(len(pts), None, dtype=object)
    best = np.full(len(pts), np.inf)
    hit = np.zeros(len(pts), dtype=bool)
    if not geoms:
        return hit, number, 0
    p_idx, g_idx = STRtree(geoms).query(pts, predicate='within')
    multi = 0
    for p, g in zip(p_idx, g_idx):
        if hit[p]:
            multi += 1
        if area[g] < best[p]:
            best[p], number[p], hit[p] = area[g], keys[g], True
    return hit, number, multi


def process(out: Path, tag: str) -> dict:
    typ = CT.typing_dir(out)
    cells = pd.read_parquet(out / 'cells' / f'{tag}.parquet',
                            columns=['Object ID', 'Centroid X µm', 'Centroid Y µm', 'Nucleus: Max diameter µm',
                                     'Nucleus: Min diameter µm'])
    qc = pd.read_parquet(typ / 'redsea' / f'{tag}.parquet', columns=['object_id', 'qc_analysis_eligible'])
    cells = qc.merge(cells, left_on='object_id', right_on='Object ID', how='left', validate='one_to_one')
    x = cells['Centroid X µm'].to_numpy(float) / CT.PIXEL_UM
    y = cells['Centroid Y µm'].to_numpy(float) / CT.PIXEL_UM
    pts = points(np.column_stack([x, y]))

    # endocrine objects
    objs, keys, meta = [], [], {}
    for f in features(out / 'islets' / f'{tag}.geojson'):
        props = f.get('properties') or {}
        cls = (props.get('classification') or {}).get('name')
        if cls not in CT.ENDOCRINE_OBJECTS or not f.get('geometry'):
            continue
        name = props.get('name')
        objs.append(valid(shape(f['geometry'])))
        keys.append(name)
        meta[name] = (cls, f.get('id'))
    in_islet, islet_num, multi = contain_smallest(objs, keys, pts)
    frame = pd.DataFrame({'object_id': cells['object_id'], 'qc_analysis_eligible': cells['qc_analysis_eligible'],
                          'x_um': cells['Centroid X µm'], 'y_um': cells['Centroid Y µm'],
                          'in_islet': in_islet, 'islet_num': islet_num})
    frame['endocrine_object_class'] = [meta[k][0] if k else None for k in islet_num]
    frame['endocrine_object_id'] = [meta[k][1] if k else None for k in islet_num]
    with np.errstate(invalid='ignore', divide='ignore'):
        frame['nucleus_elongation'] = (cells['Nucleus: Max diameter µm'] / cells['Nucleus: Min diameter µm']).to_numpy()

    # structures: centroid membership per class, and the share of each cell's area each class covers
    cell_geo = {}
    for f in features(out / 'cell_geojson' / f'{tag}.geojson'):
        cell_geo[f['id']] = f['geometry']
    polys = np.asarray([valid(shape(cell_geo[o])) for o in frame['object_id']], dtype=object)
    cell_area = shapely.area(polys)
    rec = {'image': tag, 'cells': len(frame), 'endocrine_objects': len(objs), 'cells_in_two_objects': multi}
    for name in CT.STRUCTURES:
        cls = CT.OBJECT_CLASS[name]
        geoms = [valid(shape(f['geometry'])) for f in features(out / 'structures' / cls / f'{cls}__{tag}.geojson')
                 if f.get('geometry')]
        inside = np.zeros(len(frame), bool)
        frac = np.zeros(len(frame))
        if geoms:
            tree = STRtree(geoms)
            p_idx, _ = tree.query(pts, predicate='within')
            inside[np.unique(p_idx)] = True
            c_idx, g_idx = tree.query(polys, predicate='intersects')
            inter = shapely.area(shapely.intersection(polys[c_idx], np.asarray(geoms, dtype=object)[g_idx]))
            frac = np.bincount(c_idx, weights=inter, minlength=len(frame)) / np.where(cell_area > 0, cell_area, np.nan)
            frac = np.clip(frac, 0, 1)
        frame[f'in_{name}'] = inside
        frame[CROSSED_COLUMN[name]] = frac
        rec[name] = {'objects': len(geoms), 'centroid_inside': int(inside.sum()),
                     'majority_inside': int((frac > CT.MAJORITY_AREA).sum())}
    label = np.full(len(frame), 'tissue', dtype=object)
    for name in reversed(CT.STRUCTURES):
        label[frame[f'in_{name}'].to_numpy()] = name
    frame['object_context'] = label
    frame['context'] = np.where(frame['in_islet'], 'islet', label)

    # cross-check: containment against QuPath's own parenting of the accepted objects
    parent = pd.read_parquet(out / 'context' / f'{tag}.parquet', columns=['Object ID', 'Parent'])
    parent = frame[['object_id']].merge(parent, left_on='object_id', right_on='Object ID', how='left')['Parent']
    q_in = parent.astype(str).str.match(r'^(Islet|Cluster|Single)_').to_numpy()
    rec['membership_vs_qupath_parent'] = {'agree': int((q_in == frame['in_islet']).sum()),
                                          'centroid_in_parent_not': int((frame['in_islet'] & ~q_in).sum()),
                                          'parent_in_centroid_not': int((~frame['in_islet'] & q_in).sum())}
    rec['in_islet'] = int(frame['in_islet'].sum())
    dest = typ / 'context' / f'{tag}.parquet'
    dest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(dest, index=False)
    dest.with_suffix('.json').write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec), flush=True)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    rows = CT.cohort_rows()
    import_structures(a.out, rows)
    tags = [r['image_name'] for r in rows] if a.images == ['all'] else a.images
    with ProcessPoolExecutor(a.workers) as pool:
        list(pool.map(process, [a.out] * len(tags), tags))
    return 0


if __name__ == '__main__':
    sys.exit(main())
