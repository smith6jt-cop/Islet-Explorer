#!/usr/bin/env python3
"""Stage a new runtime revision while preserving all prior data and fitted coordinates."""
import argparse
import json
import os
import shutil
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd
from common import require_cohort, sha256, write_json


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    a = p.parse_args()
    source, target = a.input_root.resolve(strict=True), a.output_root.resolve()
    if target.exists():
        raise ValueError('Use a new revision directory; existing candidates are immutable')
    manifest = json.loads((source/'manifest.json').read_text())
    cohort = pd.read_csv(source/'cohort.csv', dtype=str)
    if len(cohort) != 21 or manifest['expected_images'] != 21:
        raise ValueError('Expected the explicit all-donor candidate')
    require_cohort(pd.read_csv(source/'trajectory/features.csv', dtype={'case_id': str}), cohort.case_id)
    skip_roots = {'research', 'preview', 'evidence', 'logs'}
    skip_files = {'manifest.json', 'candidate.json', 'capacity.json', 'candidate_validation.json', 'trajectory/counts.csv'}
    records = []
    target.mkdir(parents=True)
    for path in sorted(source.rglob('*')):
        relative = path.relative_to(source)
        if relative.parts[0] in skip_roots or str(relative) in skip_files:
            continue
        destination = target/relative
        if path.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
        elif path.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            if relative.parts[0] in {'images', 'raw', 'gson', 'regional'}:
                os.link(path, destination)
                records.append(dict(path=str(relative), bytes=path.stat().st_size,
                                    device=path.stat().st_dev, inode=path.stat().st_ino))
            else:
                shutil.copy2(path, destination)
    for name in ('logs', 'evidence'):
        (target/name).mkdir(exist_ok=True)
    candidate = json.loads((source/'candidate.json').read_text())
    candidate.update(schema_version='3', revision_of=str(source), research_results_in_app=False)
    write_json(target/'candidate.json', candidate)
    write_json(target/'revision.json', dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        source_root=str(source), source_manifest_sha256=sha256(source/'manifest.json'),
        cohort_sha256=sha256(target/'cohort.csv'), schema_version='3',
        linked_files=records, instruction='Immutable linked inputs; write revised artifacts as new files, never in place.',
        scientific_acceptance=False, ufrc_accepted=False))
    print(f'Staged {len(records)} immutable files at {target}', flush=True)


if __name__ == '__main__':
    main()
