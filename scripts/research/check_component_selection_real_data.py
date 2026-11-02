#!/usr/bin/env python3
"""Read-only real-data regression for explicit endocrine component selection.

This replays archived hormone evidence and checks existing coverage receipts.
The one-image boundary control composes actual archived fields; it is explicitly
not a current-policy classifier integration run or a scientific accuracy test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/uniform/typing'))
import area_coverage as AC
import component_selection as CS


def check_coverage(old: Path, other: Path, tag: str) -> dict:
    metas = []
    for root in (old, other):
        path = root / 'coverage' / f'{tag}.parquet'
        sidecar = path.with_suffix('.json')
        meta = json.loads(sidecar.read_text())
        parquet = pq.ParquetFile(path)
        assert meta['parquet_sha256'] == CS.sha(path), path
        assert meta['schema_version'] == AC.SCHEMA_VERSION and meta['interior_px'] == AC.INTERIOR_PX
        assert meta['image'] == tag and meta['rows'] == parquet.metadata.num_rows
        assert set(meta['markers']) == {c.removeprefix('covint__') for c in parquet.schema.names
                                        if c.startswith('covint__')}
        for name in ('levels.csv', 'levels_selection.csv'):
            assert meta['inputs']['sha256'][f'typing/{name}'] == CS.sha(root / name), (path, name)
        metas.append(meta)
    current = AC.input_fingerprints(other.parent, tag)
    for key in ('coverage_code', 'pixel_processing_code', f'cell_geojson/{tag}.geojson'):
        assert metas[0]['inputs']['sha256'][key] == metas[1]['inputs']['sha256'][key] == current['sha256'][key], (tag, key)
    assert metas[0]['inputs']['channels'] == metas[1]['inputs']['channels'] == current['channels']
    assert metas[0]['inputs']['image'] == metas[1]['inputs']['image']
    assert metas[0]['inputs']['image']['size'] == current['image']['size']
    if metas[0]['inputs']['image']['path'] == current['image']['path']:
        assert metas[0]['inputs']['image'] == current['image']
    return {'rows': metas[0]['rows'], 'sources': {
        str(root / 'coverage' / f'{tag}.parquet'): metas[i]['parquet_sha256'] for i, root in enumerate((old, other))}}


def check_image(old: Path, other: Path, tag: str, negative_control: bool) -> dict:
    source_path = old / 'identity' / f'{tag}.parquet'
    source = pd.read_parquet(source_path)
    assert source.object_id.is_unique and source.object_id.notna().all()
    redsea_path = other / 'redsea' / f'{tag}.parquet'
    rs = pd.read_parquet(redsea_path, columns=['object_id', 'qc_analysis_eligible'] +
                         [f'Cell__{marker}' for marker in CS.MARKERS])
    rs = rs[rs.qc_analysis_eligible].set_index('object_id')
    assert set(rs.index) == set(source.object_id)
    rs = rs.loc[source.object_id]
    levels = pd.read_csv(old / 'levels.csv')
    levels = levels[levels.image.eq(tag)].set_index('marker').loc[list(CS.MARKERS)]
    canary = pd.read_csv(old / 'canary_levels.csv')
    canary = canary[canary.image.eq(tag)].set_index('marker').loc[list(CS.MARKERS)]
    assert levels.measurement.eq('whole_cell').all()
    values = np.log1p(np.clip(np.nan_to_num(rs[[f'Cell__{marker}' for marker in CS.MARKERS]].to_numpy(float),
                                          nan=0.0), 0, None))
    level = levels.level_log.to_numpy(float)
    p = values > level
    c = values > np.maximum(level, canary.canary_log.to_numpy(float))
    f = (values - level).astype(np.float32)
    for j, marker in enumerate(CS.MARKERS):
        assert np.array_equal(p[:, j], source[f'pos__{marker}'].to_numpy(bool)), (tag, marker, 'positivity')
        assert np.array_equal(c[:, j], source[f'canary__{marker}'].to_numpy(bool)), (tag, marker, 'canary')
        assert np.array_equal(f[:, j], source[f'fold__{marker}'].to_numpy()), (tag, marker, 'fold')
    raw, rawlevel = np.expm1(values), np.expm1(level)
    top = np.percentile(raw, 99.9, axis=0)
    scores = np.clip((raw - rawlevel) / np.where(top > rawlevel, top - rawlevel, np.nan), 0, None)
    indices = {marker: j for j, marker in enumerate(CS.MARKERS)}
    binding = {'endocrine_source': str(old)}
    selected, report = CS.apply(other.parent, tag, source, scores, indices, binding)
    assert selected.label.equals(source.label)
    evidence = [field for field in source if field.startswith(('pos__', 'canary__', 'fold__')) or field == 'Ki67_pos']
    assert selected[evidence].equals(source[evidence])
    for field in CS.DECISIONS:
        assert selected[f'component_source__{field}'].equals(source[field])
    result = {'image': tag, 'cells': len(source), 'named_endocrine_cells': report['preserved_named_endocrine_cells'],
              'unresolved_endocrine_cells': report['preserved_unresolved_endocrine_cells'],
              'all_four_marker_evidence_exact': True, 'selected_hormone_labels_exact': True,
              'all_input_marker_and_state_fields_preserved': True,
              'source_identity_sha256': CS.sha(source_path), 'source_redsea_sha256': CS.sha(redsea_path),
              'coverage': check_coverage(old, other, tag)}
    if negative_control:
        original = pd.read_parquet(other / 'identity' / f'{tag}.parquet').set_index('object_id').loc[source.object_id].reset_index()
        try:
            CS.apply(other.parent, tag, original, scores, indices, binding)
        except ValueError as error:
            result['actual_mismatched_v20_evidence_rejected'] = str(error)
        else:
            raise AssertionError(f'{tag}: actual mismatched v20 endocrine evidence was accepted')
        # Real rows/fields, deliberately composed to exercise selection policy.
        # This does not claim these natural decisions were rerun with mixed levels.
        composed = original.copy()
        endocrine_fields = [f'{prefix}{m}' for m in CS.MARKERS for prefix in ('pos__', 'canary__', 'fold__')]
        composed.loc[:, endocrine_fields] = source[endocrine_fields].to_numpy()
        changed, details = CS.apply(other.parent, tag, composed, scores, indices, binding)
        named = source.label.isin(CS.LABELS)
        unresolved = source.label.isin(CS.UNRESOLVED)
        boundary = composed.label.isin(CS.LABELS) & ~named & ~unresolved
        assert boundary.any(), 'Real control image must exercise at least one component boundary'
        assert changed.loc[boundary, 'label'].eq('Conflict Endocrine component boundary').all()
        assert changed.loc[boundary, 'status'].eq('component_boundary_conflict').all()
        assert changed.loc[unresolved, 'label'].equals(source.loc[unresolved, 'label'])
        assert changed.loc[named, 'label'].equals(source.loc[named, 'label'])
        assert changed[evidence].equals(composed[evidence])
        unaffected = ~(named | unresolved | boundary)
        assert changed.loc[unaffected, 'label'].equals(composed.loc[unaffected, 'label'])
        result['boundary_field_composition_control'] = {
            'not_natural_classifier_integration': True, 'boundary_conflicts': int(boundary.sum()),
            'unchanged_current_Ki67_and_all_marker_fields': True,
            'other_natural_labels_unchanged': True, 'details': details}
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--endocrine-source', type=Path, required=True)
    parser.add_argument('--nonendocrine-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--negative-control-image', default='ND_0112')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    old, other = CS.typing_root(args.endocrine_source), CS.typing_root(args.nonendocrine_source)
    images = sorted(p.stem for p in (old / 'identity').glob('*.parquet'))
    assert len(images) == 20 and args.negative_control_image in images
    rows = []
    for tag in images:
        row = check_image(old, other, tag, tag == args.negative_control_image)
        rows.append(row)
        print(f'{tag}: {row["cells"]} eligible cells; {row["named_endocrine_cells"]} named endocrine; replay/coverage pass', flush=True)
    result = {'schema': 'uniform.component_selection_real_data_check/1', 'passes': True,
              'not_full_pipeline_integration': True, 'scientific_accuracy_established': False,
              'images': rows, 'cells': sum(r['cells'] for r in rows),
              'named_endocrine_cells': sum(r['named_endocrine_cells'] for r in rows),
              'unresolved_endocrine_cells': sum(r['unresolved_endocrine_cells'] for r in rows),
              'source_hashes': {str(p): CS.sha(p) for p in [Path(__file__), Path(CS.__file__),
                    ROOT / 'scripts/uniform/typing/identity.py', old / 'levels.csv', old / 'canary_levels.csv']}}
    assert result['cells'] == 3199606 and result['named_endocrine_cells'] == 110721
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(args.output, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
