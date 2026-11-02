#!/usr/bin/env python3
"""Trace classifier signal through morphology and cell snapping without changing data.

For cells retaining endocrine-object labels, replay the preserved pixel classifier
on its exported QuPath planes; verify the saved mask, measure signal support on
the classifier grid, and compare exact script-1 polygon overlap with final cell
membership. Pixel support is diagnostic, not a replacement cell-positivity rule.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
from skimage.draw import polygon

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import islet_classifier_fit as fit

LABELS = ('Islet', 'EndocrineCluster', 'EndocrineSingle')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def features(path):
    raw = json.loads(path.read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def sampling_mask(geometry, grid_shape, center_offset=.5):
    """Sample at pixel centers; also report integer-grid sensitivity for small cells."""
    components = geometry['coordinates'] if geometry['type'] == 'MultiPolygon' else [geometry['coordinates']]
    rings = [np.asarray(p[0]) / 2 - center_offset for p in components]
    vertices = np.concatenate(rings)
    x0, y0 = np.maximum(np.floor(vertices.min(axis=0)).astype(int), 0)
    x1, y1 = np.minimum(np.ceil(vertices.max(axis=0)).astype(int) + 1, grid_shape[::-1])
    if x1 <= x0 or y1 <= y0:
        raise ValueError('Cell lies outside classifier plane')
    mask = np.zeros((y1 - y0, x1 - x0), bool)
    for component in components:
        local = np.zeros_like(mask)
        for i, ring in enumerate(component):
            ring = np.asarray(ring) / 2 - center_offset - (x0, y0)
            rr, cc = polygon(ring[:, 1], ring[:, 0], shape=mask.shape)
            local[rr, cc] = i == 0
        mask |= local
    return (slice(y0, y1), slice(x0, x1)), mask


def one(tag, run, base, out):
    identity_path = run / 'typing/identity' / f'{tag}.parquet'
    d = pd.read_parquet(identity_path)
    d = d[d.label.isin(LABELS)][['object_id', 'label', 'status']].copy()
    context = pd.read_parquet(run / 'typing/context' / f'{tag}.parquet',
                              columns=['object_id', 'x_um', 'y_um', 'endocrine_object_id'])
    d = d.merge(context, on='object_id', how='left', validate='one_to_one')
    params_path = base / 'fit' / f'{tag}.json'
    params = json.loads(params_path.read_text())
    if params['resolution_um'] != .7548:
        raise ValueError(f'{tag}: unexpected classifier grid')
    plane_path = base / 'planes' / f'{tag}.planes.tif'
    planes = fit.read_planes(plane_path)
    smoothed = {m: fit.smooth(a, params['sigma_um'] / params['resolution_um']) for m, a in planes.items()}
    del planes
    thresholds = params['thresholds']
    above = {m: a * np.float32(1 / thresholds[m]) > 1 for m, a in smoothed.items() if thresholds.get(m) is not None}
    direct_hormone = np.logical_or.reduce([a for m, a in above.items() if m in fit.HORMONES])
    direct_any = np.logical_or.reduce(list(above.values()))
    seeded = direct_any & fit.dilate(direct_hormone.astype(np.uint8), params['seed_px']).astype(bool)
    classifier = fit.classifier_mask(smoothed, thresholds, seed_px=params['seed_px'], close_px=params['close_px'],
                                     open_px=params['open_px'], margin_px=params['margin_px'])
    del smoothed
    saved = np.load(base / 'fit' / f'{tag}.mask.npy')
    if not np.array_equal(saved, classifier):
        raise ValueError(f'{tag}: classifier replay differs from saved mask at {np.count_nonzero(saved != classifier)} pixels')
    del saved
    stages = {**{'marker__' + m: a for m, a in above.items()}, 'direct_hormone': direct_hormone,
              'direct_any': direct_any, 'seeded': seeded, 'classifier': classifier}
    script1_path = base / 'script1_outlines' / f'{tag}.geojson'
    objects = [shapely.make_valid(shape(f['geometry'])) for f in features(script1_path)
               if f.get('properties', {}).get('classification', {}).get('name') == 'Islet']
    if not objects:
        raise ValueError(f'{tag}: no original script1 objects')
    tree = shapely.STRtree(objects)
    wanted = set(d.object_id)
    cell_path = run / 'cell_geojson' / f'{tag}.geojson'
    cell_geoms = {f['id']: f['geometry'] for f in features(cell_path) if f['id'] in wanted}
    if set(cell_geoms) != wanted:
        raise ValueError(f'{tag}: missing cell geometry')
    rows = []
    for row in d.itertuples(index=False):
        geo = cell_geoms[row.object_id]
        window, mask = sampling_mask(geo, classifier.shape)
        record = row._asdict()
        record['classifier_grid_cell_pixels'] = int(mask.sum())
        cx = int(np.clip(row.x_um / .7548, 0, classifier.shape[1] - 1))
        cy = int(np.clip(row.y_um / .7548, 0, classifier.shape[0] - 1))
        for name, stage in stages.items():
            record[name + '_pixels'] = int(stage[window][mask].sum())
            record[name + '_fraction'] = float(stage[window][mask].mean()) if mask.any() else None
            record[name + '_at_centroid'] = bool(stage[cy, cx])
        integer_window, integer_mask = sampling_mask(geo, classifier.shape, center_offset=0)
        record['integer_grid_cell_pixels'] = int(integer_mask.sum())
        for name in ('direct_any', 'direct_hormone', 'seeded', 'classifier'):
            record[name + '_integer_grid_pixels'] = int(stages[name][integer_window][integer_mask].sum())
        cell = shapely.make_valid(shape(geo))
        if cell.area <= 0:
            raise ValueError(f'{tag}/{row.object_id}: zero cell polygon area')
        hits = tree.query(cell, predicate='intersects')
        overlaps = [cell.intersection(objects[i]) for i in hits]
        record['script1_max_cell_overlap_fraction'] = max((g.area for g in overlaps), default=0) / cell.area
        record['script1_union_cell_overlap_fraction'] = shapely.union_all(overlaps).area / cell.area if overlaps else 0
        record['script1_at_centroid'] = any(objects[i].covers(shapely.Point(row.x_um / .3774, row.y_um / .3774)) for i in hits)
        rows.append(record)
    result = pd.DataFrame(rows)
    result['image'] = tag
    result.to_parquet(out / f'{tag}.parquet', index=False)
    hashes = {str(p): sha(p) for p in [identity_path, params_path, plane_path, script1_path, cell_path,
                                     base / 'fit' / f'{tag}.mask.npy', plane_path.with_suffix('.json'),
                                     run / 'typing/context' / f'{tag}.parquet']}
    return {'image': tag, 'cells': len(result), 'classifier_replay': 'exact', 'source_sha256': hashes}


def summarize(frame):
    return {'cells': len(frame),
            'no_direct_any_marker_pixels': int(frame.direct_any_pixels.eq(0).sum()),
            'no_direct_any_marker_pixels_both_grid_offsets': int((frame.direct_any_pixels.eq(0) & frame.direct_any_integer_grid_pixels.eq(0)).sum()),
            'empty_cell_sampling_mask': int(frame.classifier_grid_cell_pixels.eq(0).sum()),
            'no_direct_hormone_pixels': int(frame.direct_hormone_pixels.eq(0).sum()),
            'no_seeded_pixels': int(frame.seeded_pixels.eq(0).sum()),
            'no_classifier_pixels': int(frame.classifier_pixels.eq(0).sum()),
            'classifier_only_without_direct_marker': int((frame.direct_any_pixels.eq(0) & frame.classifier_pixels.gt(0)).sum()),
            'less_than_half_cell_in_any_single_script1_object': int(frame.script1_max_cell_overlap_fraction.lt(.5).sum()),
            'no_script1_overlap': int(frame.script1_union_cell_overlap_fraction.eq(0).sum())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--classifier-run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Output exists; choose a fresh directory')
    tags = [p.stem for p in sorted((args.run / 'typing/identity').glob('*.parquet'))]
    if not tags:
        parser.error('No identity files found')
    args.out.mkdir(parents=True)
    records = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(one, tag, args.run, args.classifier_run, args.out) for tag in tags]
        for tag, future in zip(tags, futures):
            records.append(future.result())
            print(tag, 'complete', flush=True)
    combined = pd.concat([pd.read_parquet(args.out / f'{tag}.parquet') for tag in tags], ignore_index=True)
    report = {'status': 'diagnostic_not_acceptance', 'images': records, 'totals': summarize(combined),
              'by_label': {label: summarize(group) for label, group in combined.groupby('label')},
              'limits': ['Pixel support sampled on original classifier grid, not current corrected whole-cell positivity.',
                         'Script1 overlap is exact polygon area; less than one half is consistent with fallback picking, not a replay of the full greedy assignment.']}
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['totals'], indent=2))


if __name__ == '__main__':
    main()
