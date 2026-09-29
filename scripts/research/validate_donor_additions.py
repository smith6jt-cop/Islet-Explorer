#!/usr/bin/env python3
"""Read-only verification of the phase-2 raw bundle and original preservation."""
import argparse
import json
from pathlib import Path
import pandas as pd

from audit_original_parity import audit
from common import require_fresh_output, sha256, write_json
from stage_donor_additions import ADDITIONS, BASELINE_FILES, EXPECTED_TOTALS, UNAVAILABLE, inspect_raw


def validate(baseline_root, bundle_root):
    root = Path(bundle_root).resolve(strict=True)
    manifest = json.loads((root / 'expansion_stage.json').read_text())
    if (manifest.get('kind') != 'frozen-original-plus-six-raw-donors'
            or manifest.get('status') != 'staged_raw_only'
            or manifest.get('app_ready') is not False
            or manifest.get('scientific_acceptance') is not False):
        raise ValueError('Expected a complete raw-only staging bundle')
    required = set(BASELINE_FILES) | {f'additions/raw/{cid}.h5ad' for cid in ADDITIONS}
    if not required.issubset(manifest['files']):
        raise ValueError('Bundle omits required original or additional artifacts')
    if set(manifest['additional_donors']) != ADDITIONS:
        raise ValueError('Addition cohort differs')
    if {p.name for p in (root / 'additions/raw').iterdir()} != {f'{cid}.h5ad' for cid in ADDITIONS}:
        raise ValueError('Raw bundle contains missing or unexpected additional files')
    sources = pd.read_csv(root / 'additions/source_manifest.csv', dtype=str)
    if len(sources) != 6 or set(sources.case_id) != ADDITIONS or not sources.cell_slot.eq('X').all():
        raise ValueError('Raw source manifest omits, duplicates or changes additional donors')
    for row in sources.itertuples():
        if Path(row.cells).resolve() != Path(manifest['files'][f'additions/raw/{row.case_id}.h5ad']['source']):
            raise ValueError(f'Raw source manifest disagrees with copied source: {row.case_id}')
    for name, record in manifest['files'].items():
        path = root / name
        if not path.resolve().is_relative_to(root) or path.is_symlink():
            raise ValueError(f'Bundle artifact escapes independent snapshot: {name}')
        if (sha256(path) != record['sha256']
                or sha256(record['source']) != record['sha256']):
            raise ValueError(f'Staged or source artifact changed: {name}')
        if path.samefile(record['source']):
            raise ValueError(f'Bundle must use independent copies: {name}')
    recorded = json.loads((root / 'additions/availability.json').read_text())
    if len(recorded) != 6 or {r['case_id'] for r in recorded} != ADDITIONS:
        raise ValueError('Availability omits or duplicates additional donors')
    census = [inspect_raw(root / f'additions/raw/{cid}.h5ad', cid) for cid in sorted(ADDITIONS)]
    for row in census:
        saved = next(r for r in recorded if r['case_id'] == row['case_id'])
        if any(saved.get(key) != value for key, value in row.items()):
            raise ValueError(f'Raw census differs: {row["case_id"]}')
        if saved.get('legacy_fields_available') is not False or saved.get('unavailable') != UNAVAILABLE:
            raise ValueError('Raw data cannot certify legacy field availability')
        source = sources.set_index('case_id').loc[row['case_id']]
        if any(saved.get(key) != source[key] for key in ['donor_source', 'donor_status']):
            raise ValueError('Donor provenance differs from selected source')
    totals = {key: sum(row[key] for row in census) for key in EXPECTED_TOTALS}
    if totals != EXPECTED_TOTALS or totals != manifest['totals']:
        raise ValueError('Additional raw-cell universe changed')
    parity = audit(baseline_root, root)
    if set(manifest['original_donors']) != set(parity.get('original_donors', [])):
        raise ValueError('Recorded original donor cohort differs from baseline')
    return dict(passed=parity['passed'], raw_additions_verified=True, totals=totals,
                app_ready=False, scientific_acceptance=False, original_parity=parity)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-root', required=True, type=Path)
    parser.add_argument('--bundle-root', required=True, type=Path)
    parser.add_argument('--report-root', required=True, type=Path)
    args = parser.parse_args()
    report = require_fresh_output(args.report_root, [args.baseline_root, args.bundle_root])
    result = validate(args.baseline_root, args.bundle_root)
    report.mkdir(parents=True)
    write_json(report / 'validation.json', result)
    print(json.dumps({key: value for key, value in result.items() if key != 'original_parity'}, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
