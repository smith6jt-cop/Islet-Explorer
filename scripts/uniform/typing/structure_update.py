"""Bind a structure-only revision to its unchanged calibration and cell inputs.

This is not a component-selection or release-acceptance mechanism. The baseline
must be selected explicitly. New geometry remains live; calibration context,
levels, measurements, and existing named endocrine decisions remain fixed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

import calibration_context as CC
import registry as RG

BINDING = 'typing/structure_update.json'
CALIBRATION = (
    'levels.csv', 'levels_selection.csv', 'levels_differentiation.csv',
    'levels_differentiation_selection.csv', 'canary_levels.csv', 'display_levels.csv',
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def pin(run: Path, baseline: Path, images: list[str], provenance: dict) -> dict:
    """Bind an already prepared fresh run; never overwrite data or an old pin."""
    run, baseline = run.resolve(), baseline.resolve()
    target = run / BINDING
    if run == baseline or target.exists():
        raise ValueError('A structure update requires a fresh run and an absent binding')
    if (run / 'typing' / 'identity').exists() and any((run / 'typing' / 'identity').iterdir()):
        raise ValueError('Refusing to bind a run with existing identities; prepare a fresh output root')
    if not images or len(images) != len(set(images)):
        raise ValueError('An explicit unique image list is required')
    tables = {}
    for name in CALIBRATION:
        source = baseline / 'typing' / name
        digest = sha(source)
        if sha(run / 'typing' / name) != digest:
            raise ValueError(f'{name}: structure-only revision changed calibration')
        tables[name] = {'path': str(source), 'sha256': digest}
    rows = {}
    for tag in images:
        source = baseline / 'typing' / 'identity' / f'{tag}.parquet'
        redsea = baseline / 'typing' / 'redsea' / f'{tag}.parquet'
        digest = sha(redsea)
        if sha(run / 'typing' / 'redsea' / f'{tag}.parquet') != digest:
            raise ValueError(f'{tag}: structure-only revision changed measurements/QC')
        rows[tag] = {'identity': {'path': str(source), 'sha256': sha(source)},
                     'redsea_sha256': digest}
    registry = {str(p.resolve()): sha(p) for p in sorted(RG.REGISTRY_DIR.glob('*.csv'))}
    CC.pin(run, baseline / 'typing' / 'context', images=images, source_provenance=provenance)
    record = {'schema_version': 1, 'baseline': str(baseline), 'scope': 'Lymphatic only',
              'source_provenance': provenance, 'calibration': tables, 'images': rows, 'registry': registry}
    with target.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    return record


def validate_calibration(run: Path) -> dict | None:
    path = run / BINDING
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    if record.get('schema_version') != 1 or record.get('scope') != 'Lymphatic only':
        raise ValueError('Unsupported structure-update binding')
    if set(record.get('calibration', {})) != set(CALIBRATION):
        raise ValueError('Incomplete structure-update calibration binding')
    for name, info in record['calibration'].items():
        if sha(run / 'typing' / name) != info['sha256']:
            raise ValueError(f'{name}: recalibration is forbidden in a structure-only revision')
    registry = {str(p.resolve()): sha(p) for p in sorted(RG.REGISTRY_DIR.glob('*.csv'))}
    if registry != record['registry']:
        raise ValueError('Typing policy changed during a structure-only revision')
    return record


def baseline_identity(run: Path, tag: str, object_ids: pd.Series) -> pd.DataFrame | None:
    record = validate_calibration(run)
    if record is None:
        if (run / 'structure_generation.json').exists():
            raise ValueError('Selected Lymphatic geometry requires an explicit structure-update baseline')
        return None
    CC.context_path(run, tag)  # validates pinned and live immutable cell fields
    info = record['images'][tag]
    if sha(run / 'typing' / 'redsea' / f'{tag}.parquet') != info['redsea_sha256']:
        raise ValueError(f'{tag}: structure update changed measurements/QC')
    source = Path(info['identity']['path'])
    if sha(source) != info['identity']['sha256']:
        raise ValueError(f'{tag}: baseline identity changed after binding')
    baseline = pd.read_parquet(source)
    ids = baseline['object_id']
    if (ids.isna().any() or ids.duplicated().any() or object_ids.isna().any()
            or object_ids.duplicated().any() or set(ids) != set(object_ids)):
        raise ValueError(f'{tag}: baseline and candidate eligible UUIDs differ')
    return baseline.set_index('object_id').loc[object_ids].reset_index()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--baseline', type=Path, required=True)
    ap.add_argument('--scientific-baseline', required=True)
    ap.add_argument('--selection-reason', required=True)
    ap.add_argument('--images', nargs='+', required=True)
    a = ap.parse_args()
    tags = [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())] if a.images == ['all'] else a.images
    pin(a.run, a.baseline, tags, {'scientific_baseline': a.scientific_baseline,
                               'selection_reason': a.selection_reason})
    print(f'Bound {len(tags)} images; calibration remains unchanged')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
