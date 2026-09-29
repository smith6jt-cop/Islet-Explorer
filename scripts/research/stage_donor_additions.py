#!/usr/bin/env python3
"""Stage six raw donors beside independent copies of the frozen app artifacts.

This is an offline evidence bundle, not a runtime candidate. No phenotype,
normalized-expression, scVI, composition or trajectory fields are synthesized.
Source H5AD hashes must match their original import sidecars before copying.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

import anndata as ad
import numpy as np
from scipy import sparse

from audit_original_parity import ARTIFACTS, donor_column, donor_id
from common import MARKERS, read_cohort, require_fresh_output, sha256, write_json

ADDITIONS = {'78', '6090', '6232', '6520', '6539', '6558'}
BASELINE_FILES = (*ARTIFACTS, 'master_results.xlsx')
EXPECTED_TOTALS = dict(raw_cells=943068, fully_measured_cells=942867,
                       wholly_unmeasured_cells=201, partially_measured_cells=0)
UNAVAILABLE = {
    'legacy_phenotype': 'Original initial-lineage method not recovered; no accepted replacement',
    'legacy_normalized_expression': 'Raw marker means are not legacy normalized expression',
    'scvi': 'No compatible saved scVI representation for these additions',
    'islet_composition': 'Requires validated typing and separately reviewed geometry/membership',
    'trajectory': 'Requires validated compatible features/projection preserving original coordinates',
}


def inspect_raw(path, case_id):
    """Read all raw values in bounded blocks without modifying source objects."""
    cells = ad.read_h5ad(path, backed='r')
    try:
        if list(cells.var_names) != MARKERS:
            raise ValueError(f'{case_id}: ordered raw marker panel mismatch')
        if not cells.n_obs or not cells.obs_names.is_unique or cells.obs_names.isna().any():
            raise ValueError(f'{case_id}: empty or duplicate raw cell identities')
        if 'imageid' not in cells.obs or set(cells.obs.imageid.map(donor_id)) != {case_id}:
            raise ValueError(f'{case_id}: cross-donor raw cell source')
        if 'phenotype' in cells.obs or any('scvi' in key.lower() for key in cells.obsm):
            raise ValueError(f'{case_id}: expected unphenotyped raw source without scVI')
        complete = unmeasured = partial = 0
        for start in range(0, cells.n_obs, 50000):
            block = cells.X[start:start + 50000]
            values = block.toarray() if sparse.issparse(block) else np.asarray(block)
            finite = np.isfinite(values)
            if np.isinf(values).any() or (values[finite] < 0).any():
                raise ValueError(f'{case_id}: infinite or negative raw intensity')
            complete += int(finite.all(axis=1).sum())
            unmeasured += int((~finite.any(axis=1)).sum())
            partial += int((finite.any(axis=1) & ~finite.all(axis=1)).sum())
        return dict(case_id=case_id, raw_cells=cells.n_obs, markers=cells.n_vars,
                    fully_measured_cells=complete, wholly_unmeasured_cells=unmeasured,
                    partially_measured_cells=partial, intensity_slot='X',
                    intensity_scale='raw_cell_mean',
                    centroid_units=str(cells.uns.get('centroid_coordinate_units', 'not_recorded')),
                    source_classification_role='archival QuPath annotation, not legacy phenotype')
    finally:
        cells.file.close()


def copy_verified(source, destination, expected=None):
    source, destination = Path(source).resolve(strict=True), Path(destination)
    before = sha256(source)
    if expected is not None and before != expected:
        raise ValueError(f'Source hash differs from import provenance: {source}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation also protects unexpected links at individual targets.
    with source.open('rb') as src, destination.open('xb') as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    if sha256(destination) != before or sha256(source) != before:
        raise ValueError(f'Source or copy changed while staging: {source}')
    return dict(source=str(source), sha256=before, bytes=destination.stat().st_size)


def stage(cohort_path, input_root, baseline_root, output_root):
    baseline_root = Path(baseline_root).resolve(strict=True)
    baseline_files = [baseline_root / name for name in BASELINE_FILES]
    for path in baseline_files:
        if not path.is_file():
            raise ValueError(f'Missing frozen artifact: {path}')
    original = ad.read_h5ad(baseline_root / 'adata_ins_root.h5ad', backed='r')
    try:
        originals = set(original.obs[donor_column(original.obs)].map(donor_id))
        if len(originals) != 15 or originals & ADDITIONS:
            raise ValueError('Expected the frozen 15-donor baseline')
    finally:
        original.file.close()
    cohort = read_cohort(cohort_path, input_root)
    if set(cohort.case_id) != originals | ADDITIONS:
        raise ValueError('Cohort must be the original 15 donors plus the six explicit additions')
    additions = cohort.loc[cohort.case_id.isin(ADDITIONS)].copy()
    if not additions.cell_slot.eq('X').all():
        raise ValueError('Additions must select raw H5AD X explicitly')
    source_dirs = [Path(row.cells).parent for row in additions.itertuples()]
    out = require_fresh_output(output_root, [*baseline_files, *source_dirs])
    # Inspect all provenance before creating the output directory.
    expected = {}
    for row in additions.itertuples():
        sidecar = Path(row.cells).with_suffix('.json')
        digest = json.loads(sidecar.read_text()).get('output_sha256')
        if not digest or sha256(row.cells) != digest:
            raise ValueError(f'{row.case_id}: raw source differs from its import sidecar')
        expected[row.case_id] = digest
    out.mkdir(parents=True)
    write_json(out / 'expansion_stage.json', dict(status='incomplete', app_ready=False))
    files, census = {}, []
    for name in BASELINE_FILES:
        files[name] = copy_verified(baseline_root / name, out / name)
    files['provenance/cohort-v1.csv'] = copy_verified(cohort_path, out / 'provenance/cohort-v1.csv')
    for row in additions.itertuples():
        cid = row.case_id
        relative = f'additions/raw/{cid}.h5ad'
        files[relative] = copy_verified(row.cells, out / relative, expected[cid])
        record = inspect_raw(out / relative, cid)
        record.update(donor_source=row.donor_source, donor_status=row.donor_status,
                      legacy_fields_available=False, unavailable=UNAVAILABLE.copy())
        census.append(record)
        for source, relative in [
            (Path(row.cells).with_suffix('.json'), f'additions/provenance/{cid}-import.json'),
            (row.geometry, f'additions/source_geometry/{cid}.geojson.gz'),
            (row.metadata, f'additions/provenance/{cid}-metadata.csv'),
        ]:
            files[relative] = copy_verified(source, out / relative)
        print(f'{cid}: {record["raw_cells"]:,} raw cells staged; legacy fields unavailable', flush=True)
    totals = {key: sum(r[key] for r in census) for key in
              ['raw_cells', 'fully_measured_cells', 'wholly_unmeasured_cells', 'partially_measured_cells']}
    if totals != EXPECTED_TOTALS:
        raise ValueError(f'Raw universe differs from the verified handoff: {totals}')
    # Images/regional exports remain explicit external references. They are not
    # copied or certified by this raw-cell checkpoint.
    additions.to_csv(out / 'additions/source_manifest.csv', index=False)
    write_json(out / 'additions/availability.json', census)
    manifest = dict(
        schema_version=1, kind='frozen-original-plus-six-raw-donors', status='staged_raw_only',
        created_at_utc=datetime.now(timezone.utc).isoformat(), runner_sha256=sha256(__file__),
        baseline_root=str(baseline_root), original_donors=sorted(originals),
        additional_donors=sorted(ADDITIONS), totals=totals, files=files,
        app_ready=False, scientific_acceptance=False, original_artifacts_rebuilt=False,
        additional_geometry_status='archival source only; requires separate review',
        external_references_status='images and regional exports not revalidated by this stage',
        next_checkpoint='independent raw verification and zero-tolerance original parity',
    )
    write_json(out / 'expansion_stage.json', manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', required=True, type=Path)
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--baseline-root', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args()
    result = stage(args.cohort, args.input_root, args.baseline_root, args.output_root)
    print(json.dumps({'status': result['status'], 'totals': result['totals']}, indent=2))
