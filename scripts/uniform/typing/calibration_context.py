#!/usr/bin/env python3
"""Pin the complete scientific reference context independently of live structures.

The caller explicitly selects a baseline; this module never chooses v19, v20,
or another generation. A pin copies every context column and records its source
and SHA256. Later structural updates may change live structure columns, but
must preserve the cell universe, coordinates, QC and endocrine membership.

``context_path(out, image)`` validates a bound image before returning its pinned
context. Unbound runs without a structure-generation selection retain their
ordinary typing/context path. Selected structures require an explicit pin.
Identity typing should continue to use live context: this binding fixes the
calibration reference, not structural assignment or cell identities.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

SCHEMA = 'uniform.calibration_context/1'
MANIFEST = 'calibration_context.json'
SNAPSHOT_DIR = 'calibration_context'
REQUIRED_COLUMNS = ('object_id', 'x_um', 'y_um', 'qc_analysis_eligible',
                    'endocrine_object_id', 'endocrine_object_class')
OPTIONAL_IMMUTABLE_COLUMNS = ('in_islet', 'islet_num', 'nucleus_elongation')


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _image_name(image: str) -> str:
    if not isinstance(image, str) or not image or Path(image).name != image or image in ('.', '..'):
        raise ValueError(f'invalid image name: {image!r}')
    return image


def _source_provenance(provenance: Mapping[str, Any]) -> dict:
    if not isinstance(provenance, Mapping):
        raise ValueError('source_provenance must be an explicit mapping')
    result = dict(provenance)
    for field in ('scientific_baseline', 'selection_reason'):
        if not isinstance(result.get(field), str) or not result[field].strip():
            raise ValueError(f'source_provenance requires an explicit nonempty {field}')
    # Reject non-JSON values before writing any binding files.
    json.dumps(result, allow_nan=False)
    return result


def _indexed(frame: pd.DataFrame, source: Path) -> pd.DataFrame:
    if 'object_id' not in frame:
        raise ValueError(f'{source}: missing object_id')
    ids = frame.object_id
    if ids.isna().any() or ids.duplicated().any() or ids.astype(str).str.strip().eq('').any():
        raise ValueError(f'{source}: missing or duplicate cell IDs')
    return frame.set_index('object_id')


def _same_rows(reference: pd.DataFrame, candidate: pd.DataFrame, source: Path,
               columns: list[str]) -> None:
    if len(reference) != len(candidate) or not reference.index.isin(candidate.index).all():
        raise ValueError(f'{source}: cell UUID universe differs from calibration baseline')
    candidate = candidate.reindex(reference.index)
    for column in columns:
        if column not in reference or column not in candidate:
            raise ValueError(f'{source}: missing invariant field {column}')
        equal = reference[column].eq(candidate[column]) | (reference[column].isna() & candidate[column].isna())
        if not equal.all():
            raise ValueError(f'{source}: {int((~equal).sum())} cells differ in invariant field {column}')


def _validate_context(out: Path, image: str, pinned_path: Path) -> dict:
    live_path = out / 'typing/context' / f'{image}.parquet'
    pinned_columns = set(pq.read_schema(pinned_path).names)
    live_columns = set(pq.read_schema(live_path).names)
    for source, columns in ((pinned_path, pinned_columns), (live_path, live_columns)):
        missing = (set(REQUIRED_COLUMNS) | {'context'}) - columns
        if missing:
            raise ValueError(f'{source}: missing required context columns {sorted(missing)}')
    if pinned_columns != live_columns:
        raise ValueError(f'{live_path}: context schema differs from the complete calibration baseline')
    immutable = set(REQUIRED_COLUMNS) - {'object_id'}
    immutable.update(c for c in OPTIONAL_IMMUTABLE_COLUMNS if c in pinned_columns or c in live_columns)
    immutable.update(c for c in pinned_columns | live_columns if c.startswith('qc_'))
    read_columns = ['object_id'] + sorted(immutable)
    reference = _indexed(pd.read_parquet(pinned_path, columns=read_columns), pinned_path)
    live = _indexed(pd.read_parquet(live_path, columns=read_columns), live_path)
    for source, frame in ((pinned_path, reference), (live_path, live)):
        if not np.isfinite(frame[['x_um', 'y_um']].to_numpy(float)).all():
            raise ValueError(f'{source}: missing or nonfinite cell coordinates')
        if not frame.qc_analysis_eligible.isin([False, True]).all():
            raise ValueError(f'{source}: invalid analysis-eligibility flags')
        if not frame.endocrine_object_id.notna().eq(frame.endocrine_object_class.notna()).all():
            raise ValueError(f'{source}: incomplete endocrine membership')
    _same_rows(reference, live, live_path, sorted(immutable))
    checks = {'image': image, 'cells': len(reference), 'invariant_columns': sorted(immutable),
              'context_columns': sorted(pinned_columns), 'complete_context_schema_matches': True,
              'live_context': str(live_path), 'live_context_matches': True}
    redsea_path = out / 'typing/redsea' / f'{image}.parquet'
    if redsea_path.exists():
        redsea = _indexed(pd.read_parquet(redsea_path, columns=['object_id', 'qc_analysis_eligible']), redsea_path)
        _same_rows(reference, redsea, redsea_path, ['qc_analysis_eligible'])
        checks['redsea_uuid_and_qc_match'] = True
    cells_path = out / 'cells' / f'{image}.parquet'
    if cells_path.exists():
        cells = pd.read_parquet(cells_path, columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
        cells = cells.rename(columns={'Object ID': 'object_id', 'Centroid X µm': 'x_um', 'Centroid Y µm': 'y_um'})
        _same_rows(reference, _indexed(cells, cells_path), cells_path, ['x_um', 'y_um'])
        checks['raw_cell_uuid_and_coordinates_match'] = True
    return checks


def _binding(out: Path) -> dict | None:
    manifest_path = out / 'typing' / MANIFEST
    if manifest_path.is_symlink() or (out / 'typing' / SNAPSHOT_DIR).is_symlink():
        raise ValueError(f'{out}: calibration binding and snapshot directory must not be symlinks')
    if not manifest_path.exists():
        if (out / 'typing' / SNAPSHOT_DIR).exists():
            raise ValueError(f'{out}: incomplete calibration-context pin; snapshot directory exists without manifest')
        return None
    binding = json.loads(manifest_path.read_text())
    if binding.get('schema') != SCHEMA or not isinstance(binding.get('images'), dict) or not binding['images']:
        raise ValueError(f'{manifest_path}: invalid calibration-context manifest')
    if binding.get('run') != str(out):
        raise ValueError(f'{manifest_path}: binding belongs to another run; explicitly pin the new run')
    _source_provenance(binding.get('source_provenance'))
    return binding


def _validate_image(out: Path, image: str, binding: dict) -> tuple[Path, dict]:
    image = _image_name(image)
    if image not in binding['images']:
        raise ValueError(f'{out}: image {image} is not included in the calibration-context pin')
    record = binding['images'][image]
    expected_relative = f'{SNAPSHOT_DIR}/{image}.parquet'
    if record.get('snapshot_path') != expected_relative:
        raise ValueError(f'{out}: invalid pinned path for {image}')
    path = out / 'typing' / expected_relative
    if path.is_symlink() or path.resolve().parent != (out / 'typing' / SNAPSHOT_DIR).resolve():
        raise ValueError(f'{path}: pinned context must be a copied file in the snapshot directory')
    actual = _sha256(path)
    if actual != record.get('sha256'):
        raise ValueError(f'{path}: pinned calibration context SHA256 changed')
    checks = _validate_context(out, image, path)
    checks['snapshot_sha256'] = actual
    return path, checks


def context_path(out: Path, image: str) -> Path:
    """Return a validated pinned context, or the live path for an unbound run."""
    out, image = Path(out).resolve(), _image_name(image)
    binding = _binding(out)
    if binding is None:
        if (out / 'structure_generation.json').exists():
            raise ValueError(f'{out}: selected structure generation requires an explicit calibration-context pin')
        return out / 'typing/context' / f'{image}.parquet'
    return _validate_image(out, image, binding)[0]


def validate(out: Path, images: list[str] | None = None, *, require_bound: bool = False) -> dict:
    """Verify snapshot hashes and current invariant fields; never repair a pin."""
    out = Path(out).resolve()
    binding = _binding(out)
    if binding is None:
        if require_bound:
            raise ValueError(f'{out}: no explicit calibration-context binding')
        return {'bound': False, 'run': str(out)}
    selected = list(binding['images']) if images is None else images
    if not selected or len(selected) != len(set(selected)):
        raise ValueError('validation images must be nonempty and unique')
    checks = [_validate_image(out, image, binding)[1] for image in selected]
    return {'bound': True, 'run': str(out), 'manifest_sha256': _sha256(out / 'typing' / MANIFEST),
            'source_provenance': binding['source_provenance'], 'checks': checks}


def pin(out: Path, baseline_context_dir: Path, *, images: list[str], source_provenance: Mapping[str, Any],
        expected_sha256: Mapping[str, str] | None = None) -> dict:
    """Copy an explicitly selected context baseline once; rebinding is forbidden.

    ``expected_sha256``, when supplied, must cover exactly the selected images.
    Every source is hashed before and after copying. A failed publication leaves
    an incomplete pin that fails closed, rather than reverting to live context.
    """
    out, baseline = Path(out).resolve(), Path(baseline_context_dir).resolve()
    provenance = _source_provenance(source_provenance)
    if not images or len(images) != len(set(images)):
        raise ValueError('pin images must be explicit, nonempty and unique')
    images = [_image_name(image) for image in images]
    if expected_sha256 is not None and set(expected_sha256) != set(images):
        raise ValueError('expected_sha256 must cover exactly the explicitly selected images')
    manifest_path = out / 'typing' / MANIFEST
    snapshot_dir = out / 'typing' / SNAPSHOT_DIR
    if manifest_path.exists() or manifest_path.is_symlink() or snapshot_dir.exists() or snapshot_dir.is_symlink():
        raise FileExistsError(f'{out}: calibration context already pinned or incomplete; rebinding is prohibited')
    records = {}
    for image in images:
        source = baseline / f'{image}.parquet'
        digest = _sha256(source)
        if expected_sha256 is not None and digest != expected_sha256[image]:
            raise ValueError(f'{source}: source SHA256 differs from the explicitly supplied baseline hash')
        _validate_context(out, image, source)
        records[image] = {'source_path': str(source), 'sha256': digest,
                          'snapshot_path': f'{SNAPSHOT_DIR}/{image}.parquet'}
    snapshot_dir.mkdir(parents=True, exist_ok=False)
    for image, record in records.items():
        source = Path(record['source_path'])
        copied = snapshot_dir / f'{image}.parquet'
        shutil.copyfile(source, copied)
        if _sha256(copied) != record['sha256'] or _sha256(source) != record['sha256']:
            raise ValueError(f'{source}: calibration context changed while copying')
        copied.chmod(0o444)
    binding = {'schema': SCHEMA, 'run': str(out), 'created_utc': datetime.now(timezone.utc).isoformat(),
               'baseline_context_dir': str(baseline), 'source_provenance': provenance, 'images': records}
    with manifest_path.open('x') as handle:
        json.dump(binding, handle, indent=2, allow_nan=False)
        handle.write('\n')
    manifest_path.chmod(0o444)
    return validate(out, require_bound=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    create = commands.add_parser('pin', help='explicitly copy and pin a scientific baseline once')
    create.add_argument('--run', type=Path, required=True)
    create.add_argument('--baseline-context', type=Path, required=True)
    create.add_argument('--images', nargs='+', required=True)
    create.add_argument('--scientific-baseline', required=True)
    create.add_argument('--selection-reason', required=True)
    create.add_argument('--provenance-json', type=Path, help='additional source provenance as a JSON object')
    create.add_argument('--expected-hashes', type=Path, help='JSON mapping image names to baseline context SHA256')
    check = commands.add_parser('validate', help='verify an existing pin against current run inputs')
    check.add_argument('--run', type=Path, required=True)
    check.add_argument('--images', nargs='+')
    check.add_argument('--require-bound', action='store_true')
    args = parser.parse_args()
    if args.command == 'pin':
        provenance = json.loads(args.provenance_json.read_text()) if args.provenance_json else {}
        provenance.update(scientific_baseline=args.scientific_baseline, selection_reason=args.selection_reason)
        expected = json.loads(args.expected_hashes.read_text()) if args.expected_hashes else None
        result = pin(args.run, args.baseline_context, images=args.images, source_provenance=provenance, expected_sha256=expected)
    else:
        result = validate(args.run, args.images, require_bound=args.require_bound)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
