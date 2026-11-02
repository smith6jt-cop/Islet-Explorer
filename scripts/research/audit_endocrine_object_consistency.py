#!/usr/bin/env python3
"""Audit endocrine objects against their constituent cells' current positivity.

This reads existing identity/context records and final object GeoJSON only. It
does not classify pixels, change thresholds, or relabel cells. An object with no
analysis-eligible cells is reported separately from one whose eligible cells
are all negative for INS/GCG/SST/CHGA. Raw-cell means above pixel-classifier
thresholds, when requested, are diagnostic comparisons, not classifier replay.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
CLASSES = ('Islet', 'EndocrineCluster', 'EndocrineSingle')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def objects(path: Path) -> pd.DataFrame:
    raw = json.loads(path.read_text())
    features = raw['features'] if isinstance(raw, dict) else raw
    rows = []
    for feature in features:
        props = feature.get('properties') or {}
        label = (props.get('classification') or {}).get('name')
        if label not in CLASSES or not feature.get('geometry'):
            continue
        rows.append({'endocrine_object_id': feature.get('id'),
                     'endocrine_object_class': label, 'object_name': props.get('name')})
    frame = pd.DataFrame(rows)
    if frame.empty or frame.endocrine_object_id.isna().any() or frame.endocrine_object_id.duplicated().any():
        raise ValueError(f'{path}: empty objects or missing/duplicate object IDs')
    return frame


def summarize(frame: pd.DataFrame) -> dict:
    populated = frame.analysis_cells.gt(0)
    negative = populated & frame.any_endocrine_positive_cells.eq(0)
    return {'objects': len(frame), 'objects_with_analysis_cells': int(populated.sum()),
            'objects_without_analysis_cells': int((~populated).sum()),
            'objects_with_any_endocrine_positive_cell': int(frame.any_endocrine_positive_cells.gt(0).sum()),
            'objects_with_all_analysis_cells_negative': int(negative.sum()),
            'context_cells': int(frame.context_cells.sum()), 'analysis_cells': int(frame.analysis_cells.sum()),
            'any_endocrine_positive_cells': int(frame.any_endocrine_positive_cells.sum()),
            'all_endocrine_negative_cells': int((frame.analysis_cells - frame.any_endocrine_positive_cells).sum()),
            'analysis_cells_in_fully_negative_objects': int(frame.loc[negative, 'analysis_cells'].sum()),
            'raw_means_above_fit_pixel_threshold_cells':
                int(frame.raw_mean_above_fit_pixel_threshold_cells.sum())
                if 'raw_mean_above_fit_pixel_threshold_cells' in frame else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--fit-root', type=Path, help='optional original fit directory containing <image>.json')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a fresh directory')
    paths = sorted((args.run / 'typing/identity').glob('*.parquet'))
    if not paths:
        parser.error('no identity files found')
    levels_path = args.run / 'typing/levels.csv'
    sources = {str(levels_path.resolve()): sha256(levels_path)}
    frames, checks = [], []
    for path in paths:
        tag = path.stem
        context_path = args.run / 'typing/context' / path.name
        object_path = args.run / 'islets' / f'{tag}.geojson'
        for source in (path, context_path, object_path):
            sources[str(source.resolve())] = sha256(source)
        obj = objects(object_path)
        ident = pd.read_parquet(path, columns=['object_id', 'start_class', 'label'] + ['pos__' + m for m in MARKERS])
        ctx = pd.read_parquet(context_path, columns=['object_id', 'qc_analysis_eligible',
                              'endocrine_object_id', 'endocrine_object_class'])
        if ident.object_id.duplicated().any() or ctx.object_id.duplicated().any():
            raise ValueError(f'{tag}: duplicate cell IDs')
        eligible = ctx.loc[ctx.qc_analysis_eligible, 'object_id']
        if set(eligible) != set(ident.object_id):
            raise ValueError(f'{tag}: identity and QC-eligible context cell universes differ')
        d = ident.merge(ctx, on='object_id', how='left', validate='one_to_one')
        expected_start = d.endocrine_object_class.fillna('Tissue')
        if not d.start_class.eq(expected_start).all():
            raise ValueError(f'{tag}: recorded start class disagrees with object context')
        in_object = d.endocrine_object_id.notna()
        if not d.loc[in_object, 'endocrine_object_id'].isin(obj.endocrine_object_id).all():
            raise ValueError(f'{tag}: context references an absent endocrine object')
        d = d.loc[in_object].copy()
        class_map = obj.set_index('endocrine_object_id').endocrine_object_class
        if not d.endocrine_object_class.eq(d.endocrine_object_id.map(class_map)).all():
            raise ValueError(f'{tag}: context object class disagrees with final geometry')
        if not d[['pos__' + m for m in MARKERS]].isin([0, 1, False, True]).all().all():
            raise ValueError(f'{tag}: invalid or missing endocrine positivity flags')
        d['any_endocrine_positive'] = d[['pos__' + m for m in MARKERS]].any(axis=1)
        count_columns = {m + '_positive_cells': ('pos__' + m, 'sum') for m in MARKERS}
        count_columns.update({'analysis_cells': ('object_id', 'size'),
                              'any_endocrine_positive_cells': ('any_endocrine_positive', 'sum')})
        if args.fit_root:
            fit_path = args.fit_root / f'{tag}.json'
            raw_path = args.run / 'cells' / path.name
            for source in (fit_path, raw_path):
                sources[str(source.resolve())] = sha256(source)
            fit = json.loads(fit_path.read_text())
            raw = pd.read_parquet(raw_path, columns=['Object ID'] + [f'Cell: {m}: Mean' for m in MARKERS])
            raw = raw.rename(columns={'Object ID': 'object_id'})
            d = d.merge(raw, on='object_id', how='left', validate='one_to_one', indicator=True)
            if not d._merge.eq('both').all():
                raise ValueError(f'{tag}: raw means missing for object cells')
            flags = []
            for marker in MARKERS:
                threshold = fit['thresholds'].get(marker)
                flags.append(d[f'Cell: {marker}: Mean'].to_numpy(float) > threshold
                             if threshold is not None else np.zeros(len(d), bool))
            d['raw_mean_above_fit_pixel_threshold'] = np.column_stack(flags).any(axis=1)
            count_columns['raw_mean_above_fit_pixel_threshold_cells'] = ('raw_mean_above_fit_pixel_threshold', 'sum')
        counts = d.groupby('endocrine_object_id').agg(**count_columns)
        ctx_counts = ctx.groupby('endocrine_object_id').size().rename('context_cells')
        result = obj.merge(counts, on='endocrine_object_id', how='left', validate='one_to_one')
        result = result.merge(ctx_counts, on='endocrine_object_id', how='left', validate='one_to_one')
        for column in list(count_columns) + ['context_cells']:
            result[column] = result[column].fillna(0).astype(int)
        result['all_analysis_cells_negative'] = result.analysis_cells.gt(0) & result.any_endocrine_positive_cells.eq(0)
        result['image'] = tag
        frames.append(result)
        checks.append({'image': tag, 'identity_cells': len(ident), 'object_cells': len(d),
                       'exact_context_join': True, 'object_class_matches_geometry': True})
    all_objects = pd.concat(frames, ignore_index=True)
    report = {'run': str(args.run.resolve()), 'description': __doc__, 'totals': summarize(all_objects),
              'by_class': {label: summarize(group) for label, group in all_objects.groupby('endocrine_object_class')},
              'by_image_and_class': [dict(image=tag, object_class=label, **summarize(group))
                                     for (tag, label), group in all_objects.groupby(['image', 'endocrine_object_class'])],
              'checks': checks, 'source_sha256': sources,
              'limitations': ['A negative corrected-cell mean is not proof that the classifier found no signal in an object.',
                              'Context membership uses centroids and the smallest containing object; it does not preserve classifier pixel support.',
                              'Objects with zero eligible cells are separate from fully marker-negative populated objects.']}
    args.out.mkdir(parents=True)
    all_objects.to_csv(args.out / 'objects.csv', index=False)
    all_objects.loc[all_objects.all_analysis_cells_negative].to_csv(args.out / 'fully_negative_objects.csv', index=False)
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['by_class'], indent=2))


if __name__ == '__main__':
    main()
