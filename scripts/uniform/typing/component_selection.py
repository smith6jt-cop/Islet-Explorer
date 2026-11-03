#!/usr/bin/env python3
"""Bind and audit an explicit v19-endocrine component selection and resolution.

This is a separate operation from a Lymphatic-only revision. It never recalibrates,
changes measurements, invents a named identity, or grants release acceptance.
The natural current-policy result is retained alongside each selected decision.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import area_coverage as AC
import calibration_context as CC
import registry as RG

SCHEMA = 'uniform.component_selection/2'
AUDIT_SCHEMA = 'uniform.component_selection_audit/2'
STATUS = 'assignments_resolved'
RESOLUTION_POLICY = 'exact_v19_named_membership_existing_nonendocrine_rules'
MANIFEST = 'typing/component_selection.json'
AUDIT = 'typing/component_selection_audit.json'
MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
LABELS = ('Beta cell', 'Alpha cell', 'Delta cell', 'CHGA+ hormone-negative')
UNRESOLVED = ('Islet', 'EndocrineCluster', 'EndocrineSingle', 'Endocrine')
MIXED = ('levels.csv', 'levels_selection.csv', 'levels_options.csv', 'canary_levels.csv')
UNCHANGED = ('levels_differentiation.csv', 'levels_differentiation_selection.csv',
             'display_levels.csv', 'display_levels_options.csv', 'empty_area_nulls.csv')
DECISIONS = ('start_class', 'label', 'excluded_reason', 'named_by', 'compartment',
             'compartment_basis', 'compartment_status', 'n_compartments_canary',
             'n_compartments_plain', 'compartment_depth', 'runner_up_compartment',
             'runner_up_depth', 'lineage', 'lineage_basis', 'lineage_status',
             'subtype_status', 'status', 'tags', 'confirmed', 'confirming_markers',
             'established_label')
SELECTED_IDENTITY = ('label', 'excluded_reason', 'named_by', 'compartment', 'lineage',
                     'subtype_status', 'status', 'established_label')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()


def file_record(path: Path) -> dict:
    return {'path': str(path.resolve(strict=True)), 'sha256': sha(path)}


def check_record(record: dict, path: Path | None = None) -> None:
    source = Path(record['path'])
    if path is not None and path.resolve(strict=True) != source:
        raise ValueError(f'{path}: source path changed after component binding')
    if sha(source) != record['sha256']:
        raise ValueError(f'{source}: component input changed after binding')


def typing_root(path: Path) -> Path:
    path = path.resolve(strict=True)
    return path if path.name == 'typing' else path / 'typing'


def indexed(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    frame = pd.read_parquet(path, columns=columns)
    ids = frame['object_id']
    if ids.isna().any() or ids.duplicated().any() or ids.astype(str).str.strip().eq('').any():
        raise ValueError(f'{path}: missing or duplicate UUID')
    return frame.set_index('object_id')


def equal(a: pd.DataFrame, b: pd.DataFrame, description: str) -> None:
    if set(a.index) != set(b.index) or set(a.columns) != set(b.columns):
        raise ValueError(f'{description}: keys or fields differ')
    b = b.reindex(index=a.index, columns=a.columns)
    same = a.eq(b) | (a.isna() & b.isna())
    if not same.all().all():
        raise ValueError(f'{description}: values differ in {same.columns[~same.all()].tolist()}')


def mixed_table(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    """Exact source-row selection; this does not fit or round any value."""
    columns = list(new.columns) + [c for c in old if c not in new]
    return pd.concat([new.loc[~new.marker.isin(MARKERS)], old.loc[old.marker.isin(MARKERS)]],
                     ignore_index=True).reindex(columns=columns)


def table_equal(actual: Path, expected: pd.DataFrame) -> None:
    a = pd.read_csv(actual, float_precision='round_trip')
    # Options have several rows per marker; sorting all fields also preserves duplicates.
    if set(a.columns) != set(expected.columns) or len(a) != len(expected):
        raise ValueError(f'{actual}: mixed calibration schema/row count differs')
    columns = sorted(a.columns)
    a = a[columns].sort_values(columns, na_position='first').reset_index(drop=True)
    b = expected[columns].sort_values(columns, na_position='first').reset_index(drop=True)
    equal(a, b, str(actual))


def validate_calibration(run: Path, old: Path, other: Path) -> None:
    for name in MIXED:
        a = pd.read_csv(old / name, float_precision='round_trip')
        b = pd.read_csv(other / name, float_precision='round_trip')
        table_equal(run / 'typing' / name, mixed_table(a, b))
    for name in UNCHANGED:
        if sha(run / 'typing' / name) != sha(other / name):
            raise ValueError(f'{name}: nonendocrine calibration changed')


def validate_image(run: Path, old: Path, other: Path, structures: Path, tag: str) -> dict:
    """Validate exact common cell universe and per-marker coverage derivation."""
    op = old / 'context' / f'{tag}.parquet'
    npth = run / 'typing/context' / f'{tag}.parquet'
    oc, nc = indexed(op), indexed(npth)
    immutable = [c for c in oc if c.startswith('qc_') or c in
                 ('x_um', 'y_um', 'endocrine_object_id', 'endocrine_object_class',
                  'in_islet', 'islet_num', 'nucleus_elongation')]
    equal(oc[immutable], nc[immutable], f'{tag}: v19 immutable context')
    equal(nc, indexed(structures / 'typing/context' / f'{tag}.parquet'), f'{tag}: selected live context')
    equal(nc[immutable], indexed(other / 'context' / f'{tag}.parquet')[immutable],
          f'{tag}: nonendocrine baseline immutable context')
    reference = indexed(CC.context_path(run, tag))
    equal(reference, indexed(other / 'context' / f'{tag}.parquet'),
          f'{tag}: frozen nonendocrine calibration context')
    redsea = run / 'typing/redsea' / f'{tag}.parquet'
    if sha(redsea) != sha(other / 'redsea' / f'{tag}.parquet'):
        raise ValueError(f'{tag}: changed REDSEA inputs')
    rs = indexed(redsea, ['object_id', 'qc_analysis_eligible'])
    equal(rs, nc[['qc_analysis_eligible']], f'{tag}: REDSEA universe/QC')
    oi, ni = indexed(old / 'identity' / f'{tag}.parquet'), indexed(other / 'identity' / f'{tag}.parquet')
    eligible = set(rs.index[rs.qc_analysis_eligible])
    if set(oi.index) != eligible or set(ni.index) != eligible:
        raise ValueError(f'{tag}: identity UUIDs differ from eligible universe')
    ca = indexed(old / 'coverage' / f'{tag}.parquet')
    cb = indexed(other / 'coverage' / f'{tag}.parquet')
    cc = indexed(run / 'typing/coverage' / f'{tag}.parquet')
    metadata = []
    for root, frame in ((old, ca), (other, cb)):
        meta = json.loads((root / 'coverage' / f'{tag}.json').read_text())
        if (meta.get('schema_version') != AC.SCHEMA_VERSION or meta.get('image') != tag
                or meta.get('interior_px') != AC.INTERIOR_PX
                or meta.get('rows') != len(frame)
                or set(meta.get('markers', [])) != {c.removeprefix('covint__') for c in frame if c.startswith('covint__')}
                or meta.get('parquet_sha256') != sha(root / 'coverage' / f'{tag}.parquet')):
            raise ValueError(f'{tag}: source coverage provenance does not match payload')
        for name in ('levels.csv', 'levels_selection.csv'):
            if meta['inputs']['sha256'][f'typing/{name}'] != sha(root / name):
                raise ValueError(f'{tag}: source coverage calibrated against different {name}')
        metadata.append(meta)
    current = AC.input_fingerprints(run, tag)
    for key in ('coverage_code', 'pixel_processing_code', f'cell_geojson/{tag}.geojson'):
        if metadata[0]['inputs']['sha256'][key] != metadata[1]['inputs']['sha256'][key]:
            raise ValueError(f'{tag}: source coverage pixels, geometry or processing differ')
        if current['sha256'][key] != metadata[0]['inputs']['sha256'][key]:
            raise ValueError(f'{tag}: candidate coverage uses different geometry or processing')
    if not metadata[0]['inputs']['channels'] == metadata[1]['inputs']['channels'] == current['channels']:
        raise ValueError(f'{tag}: source coverage channel order differs')
    if metadata[0]['inputs']['image'] != metadata[1]['inputs']['image']:
        raise ValueError(f'{tag}: archived source image fingerprints differ')
    old_image = metadata[0]['inputs']['image']
    if old_image['size'] != current['image']['size']:
        raise ValueError(f'{tag}: relocated image size differs from archived inputs')
    if old_image['path'] == current['image']['path'] and old_image != current['image']:
        raise ValueError(f'{tag}: source image identity changed at the same path')
    selected = [f'{prefix}{m}' for m in MARKERS for prefix in ('cov__', 'covint__')]
    equal(ca[[AC.INTERIOR_COUNT]], cb[[AC.INTERIOR_COUNT]], f'{tag}: interior raster')
    expected = cb.copy()
    expected.loc[:, selected] = ca.reindex(cb.index)[selected]
    equal(cc, expected, f'{tag}: component coverage')
    AC.read_validated(run, tag, pd.Series(list(eligible)),
                      [c.removeprefix('covint__') for c in cc if c.startswith('covint__')])
    return {'cells': len(oi), 'detections': len(rs),
            'named_endocrine_cells': int(oi.label.isin(LABELS).sum()),
            'unresolved_endocrine_cells': int(oi.label.isin(UNRESOLVED).sum())}


def pin(run: Path, endocrine_source: Path, nonendocrine_source: Path, structure_source: Path,
        images: list[str], selection_reason: str) -> dict:
    run = run.resolve(strict=True)
    old, other = typing_root(endocrine_source), typing_root(nonendocrine_source)
    structures = structure_source.resolve(strict=True)
    if not selection_reason.strip() or not images or len(images) != len(set(images)):
        raise ValueError('Explicit selection reason and unique images are required')
    if (run / MANIFEST).exists() or (run / 'typing/structure_update.json').exists():
        raise ValueError('Component selection requires an absent binding and no structure-only scope')
    if (run / 'typing/identity').exists() and any((run / 'typing/identity').iterdir()):
        raise ValueError('Component selection requires a fresh identity output directory')
    if run == old.parent or run == other.parent or run == structures:
        raise ValueError('Component outputs must be separate from both sources')
    validate_calibration(run, old, other)
    files = {}
    for name in MIXED + UNCHANGED:
        files[f'candidate/{name}'] = file_record(run / 'typing' / name)
        files[f'nonendocrine/{name}'] = file_record(other / name)
        if name in MIXED:
            files[f'endocrine/{name}'] = file_record(old / name)
    for p in sorted(RG.REGISTRY_DIR.glob('*.csv')):
        files[f'registry/{p.name}'] = file_record(p)
    for name in ('identity.py', 'component_selection.py', 'area_coverage.py', 'calibration_context.py'):
        files[f'code/{name}'] = file_record(Path(__file__).parent / name)
    # Context is explicitly frozen independently of live DG geometry, even though
    # this operation selects existing thresholds and does not refit them.
    for tag in images:
        CC.context_path(run, tag)
    files['calibration_context'] = file_record(run / 'typing/calibration_context.json')
    derivation = json.loads((run / 'typing/coverage/component_derivation.json').read_text())
    if (derivation.get('schema') != 'uniform.component_coverage_derivation/1'
            or derivation.get('pixels_recomputed') is not False
            or set(derivation.get('images', {})) != set(images)):
        raise ValueError('Missing or inconsistent explicit coverage derivation receipt')
    files['coverage_derivation'] = file_record(run / 'typing/coverage/component_derivation.json')
    for name in ('structure_generation.json', 'structure_metrics_manifest.json'):
        files[f'structure_source/{name}'] = file_record(structures / name)
        files[f'candidate_structure/{name}'] = file_record(run / name)
    rows = {}
    for tag in images:
        if Path(tag).name != tag:
            raise ValueError(f'Invalid image name {tag}')
        stats = validate_image(run, old, other, structures, tag)
        inputs = {}
        inputs[f'structure_source/context/{tag}.parquet'] = file_record(structures / 'typing/context' / f'{tag}.parquet')
        for source, root in (('endocrine', old), ('nonendocrine', other), ('candidate', run / 'typing')):
            names = [f'context/{tag}.parquet', f'context/{tag}.json',
                     f'coverage/{tag}.parquet', f'coverage/{tag}.json']
            if source != 'candidate':
                names.append(f'identity/{tag}.parquet')
            if source != 'endocrine':
                names.append(f'redsea/{tag}.parquet')
            for name in names:
                inputs[f'{source}/{name}'] = file_record(root / name)
        rows[tag] = {**stats, 'inputs': inputs}
    record = {'schema': SCHEMA, 'scope': 'explicit endocrine component selection',
              'status': STATUS, 'release_ready': False, 'resolution_policy': RESOLUTION_POLICY,
              'run': str(run), 'endocrine_source': str(old), 'nonendocrine_source': str(other),
              'structure_source': str(structures),
              'selection_reason': selection_reason, 'endocrine_markers': list(MARKERS),
              'endocrine_labels': list(LABELS), 'files': files, 'images': rows,
              'decision_provenance': 'Archived v19 identity tuple only for selected named endocrine cells; unresolved/boundary cells use '
                  'existing nonendocrine rules or explicit Unclassified. Original decision metadata is component_source__*. Marker flags, '
                  'counts, competition diagnostics and Ki67 remain current. Natural decisions are natural__*. '
                  'Selection basis/status and resolution targets are explicit; current confirmation is unset for selected rows.'}
    with (run / MANIFEST).open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    return record


def validate(run: Path, image: str | None = None) -> dict | None:
    run = run.resolve()
    path = run / MANIFEST
    if not path.exists():
        return None
    if path.is_symlink() or (run / 'typing/structure_update.json').exists():
        raise ValueError('Component binding must be local and cannot coexist with a structure-only binding')
    record = json.loads(path.read_text())
    if (record.get('schema') != SCHEMA or record.get('run') != str(run)
            or record.get('scope') != 'explicit endocrine component selection'
            or record.get('status') != STATUS or record.get('release_ready') is not False
            or record.get('resolution_policy') != RESOLUTION_POLICY
            or record.get('endocrine_markers') != list(MARKERS) or record.get('endocrine_labels') != list(LABELS)):
        raise ValueError('Unsupported or relocated component-selection binding')
    if not record.get('images') or not record.get('files'):
        raise ValueError('Incomplete component binding')
    old, other = Path(record['endocrine_source']), Path(record['nonendocrine_source'])
    structures = Path(record['structure_source'])
    expected_files = {f'candidate/{x}': run / 'typing' / x for x in MIXED + UNCHANGED}
    expected_files.update({f'nonendocrine/{x}': other / x for x in MIXED + UNCHANGED})
    expected_files.update({f'endocrine/{x}': old / x for x in MIXED})
    expected_files.update({f'registry/{p.name}': p for p in sorted(RG.REGISTRY_DIR.glob('*.csv'))})
    expected_files.update({f'code/{x}': Path(__file__).parent / x for x in
                          ('identity.py', 'component_selection.py', 'area_coverage.py', 'calibration_context.py')})
    expected_files['calibration_context'] = run / 'typing/calibration_context.json'
    expected_files['coverage_derivation'] = run / 'typing/coverage/component_derivation.json'
    for name in ('structure_generation.json', 'structure_metrics_manifest.json'):
        expected_files[f'structure_source/{name}'] = structures / name
        expected_files[f'candidate_structure/{name}'] = run / name
    if set(expected_files) != set(record['files']):
        raise ValueError('Component binding has missing or unexpected file records')
    for key, expected_path in expected_files.items():
        check_record(record['files'][key], expected_path)
    if image is not None:
        if image not in record['images']:
            raise ValueError(f'{image}: not bound by component selection')
        CC.context_path(run, image)
        inputs = record['images'][image]['inputs']
        expected = {f'structure_source/context/{image}.parquet': structures / 'typing/context' / f'{image}.parquet'}
        for source, root in (('endocrine', old), ('nonendocrine', other), ('candidate', run / 'typing')):
            names = [f'context/{image}.parquet', f'context/{image}.json',
                     f'coverage/{image}.parquet', f'coverage/{image}.json']
            if source != 'candidate':
                names.append(f'identity/{image}.parquet')
            if source != 'endocrine':
                names.append(f'redsea/{image}.parquet')
            expected.update({f'{source}/{name}': root / name for name in names})
        if set(inputs) != set(expected):
            raise ValueError(f'{image}: incomplete component input inventory')
        for key, source in expected.items():
            check_record(inputs[key], source)
    return record


def source_identity(record: dict, tag: str, object_ids: pd.Series) -> pd.DataFrame:
    source = indexed(Path(record['endocrine_source']) / 'identity' / f'{tag}.parquet')
    if object_ids.isna().any() or object_ids.duplicated().any() or set(source.index) != set(object_ids):
        raise ValueError(f'{tag}: source identity UUID mismatch')
    return source.reindex(object_ids).reset_index()


def resolution_mask(record: dict, tag: str, natural: pd.DataFrame) -> np.ndarray:
    """Only unresolved source/current definitions and unselected endocrine calls."""
    source = source_identity(record, tag, natural.object_id)
    selected = source.label.isin(LABELS).to_numpy()
    return ~selected & (source.label.isin(UNRESOLVED).to_numpy()
                        | natural.label.isin((*LABELS, *UNRESOLVED)).to_numpy())


def apply(run: Path, tag: str, natural: pd.DataFrame, scores: np.ndarray,
          marker_indices: dict[str, int], record: dict, *, resolved: pd.DataFrame,
          resolution_target: np.ndarray) -> tuple[pd.DataFrame, dict]:
    """Preserve exact named source calls; resolve only the declared nonselected set.

    Hormone measurements/flags stay intact. The caller reruns the existing rules
    without object-forced endocrine placement or new endocrine names for this
    set. Unsupported anatomical labels become explicit Unclassified, retaining
    the original object class, source calls, natural calls and reasons.
    """
    source = source_identity(record, tag, natural.object_id)
    if source.label.isna().any():
        raise ValueError(f'{tag}: missing source identity')
    fields = [f'{prefix}{m}' for m in MARKERS for prefix in ('pos__', 'canary__', 'fold__')]
    equal(source[fields], natural[fields], f'{tag}: archived endocrine positivity/canary/fold evidence')
    equal(source[['start_class']], natural[['start_class']], f'{tag}: endocrine start class')
    selected = source.label.isin(LABELS).to_numpy()
    source_unresolved = source.label.isin(UNRESOLVED).to_numpy()
    p = natural[[f'pos__{m}' for m in MARKERS[:3]]].to_numpy(bool)
    values = np.where(p, np.nan_to_num(scores[:, [marker_indices[m] for m in MARKERS[:3]]], nan=0.0), -np.inf)
    ordered = np.sort(values, axis=1)
    any_hormone = np.isfinite(ordered[:, -1])
    unique = any_hormone & (ordered[:, -1] > ordered[:, -2])
    expected = np.full(len(natural), '', dtype=object)
    expected[unique] = np.asarray(LABELS[:3], dtype=object)[values.argmax(axis=1)[unique]]
    chga = ~any_hormone & natural['pos__CHGA'].astype(bool).to_numpy()
    expected[chga] = LABELS[3]
    if np.any(expected[selected] != source.label.to_numpy()[selected]):
        bad = natural.loc[selected & (expected != source.label.to_numpy()), 'object_id'].tolist()
        raise ValueError(f'{tag}: {len(bad)} selected v19 identities fail hormone evidence replay; first {bad[:5]}')
    target = resolution_mask(record, tag, natural)
    if not np.array_equal(target, resolution_target):
        raise ValueError(f'{tag}: resolution target differs from the declared source/natural boundary')
    if list(resolved.columns) != list(natural.columns):
        raise ValueError(f'{tag}: resolution changed the identity schema')
    immutable = [c for c in natural if c not in DECISIONS or c == 'start_class']
    equal(natural[immutable], resolved[immutable], f'{tag}: resolution altered measurements or marker states')
    equal(natural.loc[~target], resolved.loc[~target], f'{tag}: resolution altered an untargeted cell')
    if resolved.loc[target, 'label'].isin(LABELS).any():
        raise ValueError(f'{tag}: nonselected resolution created a new endocrine identity')
    result = resolved.copy()
    for field in DECISIONS:
        result[f'natural__{field}'] = natural[field]
        result[f'component_source__{field}'] = source[field]
    result['component_source_label'] = source.label.to_numpy()
    result['component_source_compartment'] = source.compartment.to_numpy()
    result['component_source_lineage'] = source.lineage.to_numpy()
    result['component_resolution_target'] = target.astype(np.int8)
    result['component_decision_source'] = np.where(target, 'nonselected_existing_rules', 'current_mixed_policy')
    result.loc[selected, list(SELECTED_IDENTITY)] = source.loc[selected, list(SELECTED_IDENTITY)].to_numpy()
    for field in ('compartment_basis', 'lineage_basis', 'compartment_status', 'lineage_status'):
        result.loc[selected, field] = 'v19_selected_component'
    result.loc[selected, 'confirmed'] = np.nan
    result.loc[selected, 'confirming_markers'] = ''
    result.loc[selected, 'component_decision_source'] = 'v19_named_endocrine'
    unsupported = target & result.label.isin((*UNRESOLVED, 'Tissue')).to_numpy()
    result.loc[unsupported, 'label'] = 'Unclassified'
    # A context class or a non-naming compartment is not a final biological call.
    # Its evidence remains in the original natural/source fields and marker flags.
    for field in ('compartment', 'lineage', 'named_by', 'confirming_markers'):
        result.loc[unsupported, field] = ''
    for field in ('compartment_basis', 'lineage_basis', 'compartment_status', 'lineage_status'):
        result.loc[unsupported, field] = 'unclassified'
    result.loc[unsupported, 'confirmed'] = np.nan
    result.loc[unsupported, 'component_decision_source'] = 'nonselected_unclassified'
    result.loc[unsupported, 'tags'] = result.loc[unsupported, 'tags'].map(
        lambda x: ';'.join(filter(None, [str(x), 'selected_endocrine_membership_other_identity_unsupported'])))
    actual = result.label.isin(LABELS).to_numpy()
    if not np.array_equal(actual, selected) or not result.loc[selected, 'label'].equals(source.loc[selected, 'label']):
        raise ValueError(f'{tag}: selected endocrine UUID/label contract failed')
    unresolved_final = int(result.label.isin(UNRESOLVED).sum())
    boundary_conflicts = int(result.label.eq('Conflict Endocrine component boundary').sum())
    if unresolved_final or boundary_conflicts:
        raise ValueError(f'{tag}: unresolved endocrine definitions remain after component resolution')
    report = {'image': tag, 'cells': len(result), 'status': STATUS, 'release_ready': False,
              'resolution_policy': RESOLUTION_POLICY,
              'preserved_named_endocrine_cells': int(selected.sum()),
              'source_unresolved_endocrine_cells': int(source_unresolved.sum()),
              'unresolved_endocrine_cells': unresolved_final, 'boundary_conflicts': boundary_conflicts,
              'resolved_nonselected_cells': int(target.sum()),
              'unclassified_cells': int(result.label.eq('Unclassified').sum()),
              'source_literal_unclassified_cells': int(source.label.eq('Unclassified').sum()),
              'natural_literal_unclassified_cells': int(natural.label.eq('Unclassified').sum()),
              'source_unnamed_anatomical_cells': int(source.label.isin((*UNRESOLVED, 'Tissue')).sum()),
              'natural_unnamed_anatomical_cells': int(natural.label.isin((*UNRESOLVED, 'Tissue')).sum()),
              'unclassified_reasons': [dict(source_label=str(a), natural_label=str(b), reason=str(c), cells=int(n))
                  for (a, b, c), n in pd.DataFrame({'source': source.loc[unsupported, 'label'],
                      'natural': natural.loc[unsupported, 'label'], 'reason': result.loc[unsupported, 'status']})
                      .value_counts().items()],
              'natural_named_endocrine_cells': int(natural.label.isin(LABELS).sum()),
              'natural_vs_selected_label_changes': int(result.label.ne(natural.label).sum()),
              'natural_unselected_endocrine': int((natural.label.isin(LABELS).to_numpy() & ~selected).sum()),
              'resolution_label_transitions': [dict(before=str(a), after=str(b), cells=int(n))
                  for (a, b), n in pd.DataFrame({'before': natural.loc[target, 'label'],
                                               'after': result.loc[target, 'label']}).value_counts().items()],
              'natural_policy': 'Current mixed-calibration policy with selected v19 endocrine protected from structure overwrite; before component resolution',
              'selected_endocrine_hormone_evidence_replayed': True,
              'resolution_marker_state_invariants_passed': True,
              'untargeted_identity_invariants_passed': True}
    return result, report

def write_audit(run: Path, reports: list[dict]) -> dict:
    record = validate(run)
    if record is None or {r['image'] for r in reports} != set(record['images']):
        raise ValueError('Component audit requires every bound image exactly once')
    if len(reports) != len(record['images']):
        raise ValueError('Duplicate image in component audit')
    counts = ('cells', 'preserved_named_endocrine_cells', 'source_unresolved_endocrine_cells',
              'unresolved_endocrine_cells', 'boundary_conflicts', 'resolved_nonselected_cells', 'unclassified_cells',
              'natural_named_endocrine_cells',
              'natural_vs_selected_label_changes', 'natural_unselected_endocrine')
    images = {}
    for report in reports:
        tag = report['image']
        validate(run, tag)
        images[tag] = {**report, 'identity': file_record(run / 'typing/identity' / f'{tag}.parquet')}
    audit = {'schema': AUDIT_SCHEMA, 'status': STATUS,
             'release_ready': False, 'manifest_sha256': sha(run / MANIFEST), 'images': images,
             **{key: sum(r[key] for r in reports) for key in counts},
             'scientific_accuracy_established': False,
             'resolution_policy': RESOLUTION_POLICY,
             'note': 'Exact selected component preservation and mathematical replay; nonselected unresolved calls '
                     'use existing nonendocrine rules or explicit Unclassified. This is not biological validation.'}
    with (run / AUDIT).open('x') as stream:
        json.dump(audit, stream, indent=2)
        stream.write('\n')
    return audit


def validate_audit(run: Path) -> dict | None:
    record = validate(run)
    if record is None:
        return None
    audit = json.loads((run / AUDIT).read_text())
    if (audit.get('schema') != AUDIT_SCHEMA
            or audit.get('manifest_sha256') != sha(run / MANIFEST)
            or audit.get('status') != STATUS or audit.get('release_ready') is not False
            or audit.get('scientific_accuracy_established') is not False
            or audit.get('resolution_policy') != RESOLUTION_POLICY
            or set(audit.get('images', {})) != set(record['images'])):
        raise ValueError('Incomplete or mismatched component-selection audit')
    for tag, row in audit['images'].items():
        if (row.get('image') != tag or row.get('status') != STATUS
                or row.get('release_ready') is not False or row.get('cells') != record['images'][tag]['cells']):
            raise ValueError(f'{tag}: inconsistent per-image diagnostic audit')
        validate(run, tag)
        check_record(row['identity'], run / 'typing/identity' / f'{tag}.parquet')
        if (row.get('resolution_policy') != RESOLUTION_POLICY
                or row.get('resolution_marker_state_invariants_passed') is not True
                or row.get('untargeted_identity_invariants_passed') is not True
                or row.get('unresolved_endocrine_cells') != 0 or row.get('boundary_conflicts') != 0):
            raise ValueError(f'{tag}: incomplete component resolution evidence')
        if row.get('selected_endocrine_hormone_evidence_replayed') is not True:
            raise ValueError(f'{tag}: missing independent hormone-evidence replay')
        if row['preserved_named_endocrine_cells'] != record['images'][tag]['named_endocrine_cells']:
            raise ValueError(f'{tag}: selected endocrine count differs from source binding')
        if row['source_unresolved_endocrine_cells'] != record['images'][tag]['unresolved_endocrine_cells']:
            raise ValueError(f'{tag}: unresolved source census differs from binding')
        if (row['preserved_named_endocrine_cells'] + row['resolved_nonselected_cells'] > row['cells']
                or row['unclassified_cells'] > row['cells'] - row['preserved_named_endocrine_cells']):
            raise ValueError(f'{tag}: component resolution census is inconsistent')
    for key in ('cells', 'preserved_named_endocrine_cells', 'source_unresolved_endocrine_cells',
                'unresolved_endocrine_cells', 'boundary_conflicts', 'resolved_nonselected_cells', 'unclassified_cells',
              'natural_named_endocrine_cells',
                'natural_vs_selected_label_changes', 'natural_unselected_endocrine'):
        if audit.get(key) != sum(row[key] for row in audit['images'].values()):
            raise ValueError(f'Component audit has an inconsistent {key} total')
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--endocrine-source', required=True, type=Path)
    parser.add_argument('--nonendocrine-source', required=True, type=Path)
    parser.add_argument('--structure-source', required=True, type=Path)
    parser.add_argument('--images', nargs='+', required=True)
    parser.add_argument('--selection-reason', required=True)
    args = parser.parse_args()
    images = [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())] if args.images == ['all'] else args.images
    record = pin(args.run, args.endocrine_source, args.nonendocrine_source, args.structure_source,
                 images, args.selection_reason)
    print(f'Bound {len(record["images"])} images for exact named endocrine component resolution; release_ready=false')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
