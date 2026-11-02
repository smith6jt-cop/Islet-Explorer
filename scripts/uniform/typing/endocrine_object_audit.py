#!/usr/bin/env python3
"""Check that every final endocrine object has a marker-positive member cell.

Read the selected levels, REDSEA measurements, context membership and final
object inventory directly, before identity typing. Positivity is the same
strict log1p comparison used by identity.py. At least one analysis-eligible
member must be positive for INS, GCG, SST or CHGA; other member cells may be
negative. Empty objects and EndocrineSingle objects without exactly one
analysis-eligible member are separate membership failures.

This is a consistency gate, not a replay of the pixel classifier. Failure does
not prove absence of raw image signal and does not change levels or labels.
Reports must be written to a fresh directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
CLASSES = ('Islet', 'EndocrineCluster', 'EndocrineSingle')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def object_inventory(path: Path) -> pd.DataFrame:
    raw = json.loads(path.read_text())
    features = raw['features'] if isinstance(raw, dict) else raw
    rows = []
    for feature in features:
        props = feature.get('properties') or {}
        label = (props.get('classification') or {}).get('name')
        if label not in CLASSES:
            continue
        if not feature.get('geometry'):
            raise ValueError(f'{path}: endocrine object has no geometry')
        rows.append({'endocrine_object_id': feature.get('id'),
                     'endocrine_object_class': label, 'object_name': props.get('name')})
    frame = pd.DataFrame(rows, columns=['endocrine_object_id', 'endocrine_object_class', 'object_name'])
    if (frame.endocrine_object_id.isna().any() or frame.endocrine_object_id.eq('').any()
            or frame.endocrine_object_id.duplicated().any()):
        raise ValueError(f'{path}: missing or duplicate endocrine object IDs')
    return frame


def summarize(frame: pd.DataFrame) -> dict:
    negative = frame.all_analysis_cells_negative
    return {'objects': len(frame),
            'objects_with_analysis_cells': int(frame.analysis_cells.gt(0).sum()),
            'objects_without_context_cells': int(frame.context_cells.eq(0).sum()),
            'objects_without_analysis_cells': int(frame.no_analysis_cells.sum()),
            'single_membership_anomalies': int(frame.single_membership_anomaly.sum()),
            'membership_anomalies': int(frame.membership_anomaly.sum()),
            'objects_with_any_endocrine_positive_cell': int(frame.any_endocrine_positive_cells.gt(0).sum()),
            'objects_with_all_analysis_cells_negative': int(negative.sum()),
            'context_cells': int(frame.context_cells.sum()),
            'analysis_cells': int(frame.analysis_cells.sum()),
            'any_endocrine_positive_cells': int(frame.any_endocrine_positive_cells.sum()),
            'all_endocrine_negative_cells': int((frame.analysis_cells - frame.any_endocrine_positive_cells).sum()),
            'analysis_cells_in_fully_negative_objects': int(frame.loc[negative, 'analysis_cells'].sum())}


def audit(run: Path, out: Path, *, levels_path: Path | None = None,
          images: list[str] | None = None) -> dict:
    """Write a read-only audit and return its report, including ``passes``.

    ``levels_path`` lets levels.py audit its new candidate table, including a
    --report trial, without accidentally reading the run's previous levels.
    """
    run, out = run.resolve(), out.resolve()
    if out.exists():
        raise FileExistsError(f'{out}: output already exists; choose a fresh directory')
    levels_path = (levels_path or run / 'typing/levels.csv').resolve()
    tags = images if images is not None else sorted(p.stem for p in (run / 'islets').glob('*.geojson'))
    if not tags or len(tags) != len(set(tags)):
        raise ValueError('images must be nonempty and unique')
    levels = pd.read_csv(levels_path)
    levels = levels.loc[levels.image.isin(tags) & levels.marker.isin(MARKERS)]
    if levels.duplicated(['image', 'marker']).any():
        raise ValueError(f'{levels_path}: duplicate endocrine levels')
    sources = {str(levels_path): sha256(levels_path)}
    frames, checks = [], []
    for tag in tags:
        selected = levels.loc[levels.image.eq(tag)].set_index('marker')
        if set(selected.index) != set(MARKERS):
            raise ValueError(f'{tag}: all four endocrine marker levels are required')
        if not selected.measurement.eq('whole_cell').all():
            raise ValueError(f'{tag}: endocrine positivity must use whole-cell measurements')
        if not np.isfinite(selected.level_log.to_numpy(float)).all():
            raise ValueError(f'{tag}: nonfinite endocrine positivity level')
        rs_path = run / 'typing/redsea' / f'{tag}.parquet'
        ctx_path = run / 'typing/context' / f'{tag}.parquet'
        obj_path = run / 'islets' / f'{tag}.geojson'
        for source in (rs_path, ctx_path, obj_path):
            sources[str(source)] = sha256(source)
        obj = object_inventory(obj_path)
        rs = pd.read_parquet(rs_path, columns=['object_id', 'qc_analysis_eligible'] + ['Cell__' + m for m in MARKERS])
        ctx = pd.read_parquet(ctx_path, columns=['object_id', 'qc_analysis_eligible',
                                               'endocrine_object_id', 'endocrine_object_class'])
        for name, frame in (('REDSEA', rs), ('context', ctx)):
            if frame.object_id.isna().any() or frame.object_id.duplicated().any():
                raise ValueError(f'{tag}: missing or duplicate {name} cell IDs')
            if not frame.qc_analysis_eligible.isin([False, True]).all():
                raise ValueError(f'{tag}: invalid {name} analysis-eligibility flags')
        d = rs.merge(ctx, on='object_id', how='outer', validate='one_to_one',
                     suffixes=('_redsea', '_context'), indicator=True)
        if not d._merge.eq('both').all():
            raise ValueError(f'{tag}: REDSEA and context cell universes differ')
        if not d.qc_analysis_eligible_redsea.eq(d.qc_analysis_eligible_context).all():
            raise ValueError(f'{tag}: REDSEA and context analysis eligibility disagree')
        in_object = d.endocrine_object_id.notna()
        if not in_object.eq(d.endocrine_object_class.notna()).all():
            raise ValueError(f'{tag}: incomplete context object membership')
        class_map = obj.set_index('endocrine_object_id').endocrine_object_class
        if not d.loc[in_object, 'endocrine_object_class'].eq(
                d.loc[in_object, 'endocrine_object_id'].map(class_map)).all():
            raise ValueError(f'{tag}: context references an absent object or mismatched object class')
        d = d.loc[in_object & d.qc_analysis_eligible_redsea.astype(bool)].copy()
        for marker in MARKERS:
            x = d['Cell__' + marker].to_numpy(float)
            d['pos__' + marker] = np.log1p(np.clip(np.nan_to_num(x, nan=0.0), 0, None)) > selected.loc[marker, 'level_log']
        d['any_endocrine_positive'] = d[['pos__' + m for m in MARKERS]].any(axis=1)
        count_columns = {m + '_positive_cells': ('pos__' + m, 'sum') for m in MARKERS}
        count_columns.update({'analysis_cells': ('object_id', 'size'),
                              'any_endocrine_positive_cells': ('any_endocrine_positive', 'sum')})
        counts = d.groupby('endocrine_object_id').agg(**count_columns)
        context_counts = ctx.groupby('endocrine_object_id').size().rename('context_cells')
        result = obj.merge(counts, on='endocrine_object_id', how='left', validate='one_to_one')
        result = result.merge(context_counts, on='endocrine_object_id', how='left', validate='one_to_one')
        for column in list(count_columns) + ['context_cells']:
            result[column] = result[column].fillna(0).astype(int)
        result['all_analysis_cells_negative'] = result.analysis_cells.gt(0) & result.any_endocrine_positive_cells.eq(0)
        result['no_analysis_cells'] = result.analysis_cells.eq(0)
        result['single_membership_anomaly'] = result.endocrine_object_class.eq('EndocrineSingle') & result.analysis_cells.ne(1)
        result['membership_anomaly'] = result.no_analysis_cells | result.single_membership_anomaly
        result['image'] = tag
        frames.append(result)
        checks.append({'image': tag, 'redsea_cells': len(rs), 'context_cells': len(ctx),
                       'analysis_cells': int(rs.qc_analysis_eligible.sum()), 'analysis_cells_in_objects': len(d),
                       'exact_context_join': True, 'matching_analysis_eligibility': True,
                       'object_class_matches_geometry': True})
    all_objects = pd.concat(frames, ignore_index=True)
    totals = summarize(all_objects)
    report = {'run': str(run), 'levels': str(levels_path), 'images': tags,
              'passes': not bool(all_objects.all_analysis_cells_negative.any() or all_objects.membership_anomaly.any()),
              'rule': 'Every final endocrine object has at least one analysis-eligible cell positive for INS, GCG, SST or CHGA; '
                      'EndocrineSingle has exactly one analysis-eligible member.',
              'totals': totals,
              'by_class': {label: summarize(group) for label, group in all_objects.groupby('endocrine_object_class')},
              'by_image_and_class': [dict(image=tag, object_class=label, **summarize(group))
                                     for (tag, label), group in all_objects.groupby(['image', 'endocrine_object_class'])],
              'checks': checks, 'source_sha256': sources,
              'limitations': ['This checks object/cell consistency; it does not replay the raw-image pixel classifier.',
                              'Other cells within a marker-positive object need not be endocrine-marker positive.',
                              'Context membership is centroid-based and selects the smallest containing object.']}
    out.mkdir(parents=True)
    all_objects.to_csv(out / 'objects.csv', index=False)
    all_objects.loc[all_objects.all_analysis_cells_negative].to_csv(out / 'fully_negative_objects.csv', index=False)
    all_objects.loc[all_objects.membership_anomaly].to_csv(out / 'membership_anomalies.csv', index=False)
    (out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--levels', type=Path, help='candidate levels table; default: RUN/typing/levels.csv')
    parser.add_argument('--images', nargs='+', help='default: every image in the final islets directory')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('output already exists; choose a fresh directory')
    report = audit(args.run, args.out, levels_path=args.levels, images=args.images)
    print(json.dumps({'passes': report['passes'], 'totals': report['totals'], 'by_class': report['by_class']}, indent=2))
    return 0 if report['passes'] else 4


if __name__ == '__main__':
    sys.exit(main())
