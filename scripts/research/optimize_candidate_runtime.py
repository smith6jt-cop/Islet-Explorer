#!/usr/bin/env python3
"""Materialize compact tables and exact single-islet downloads before serving."""
import argparse
import gzip
import json
import time
from pathlib import Path

import pandas as pd
from common import sha256, write_json
from finalize_candidate import artifact


def refresh_manifest(root, report):
    manifest=json.loads((root/'manifest.json').read_text())
    manifest['artifacts']['trajectory']=artifact(root, 'trajectory')
    manifest['artifacts']['islet_boundaries']=artifact(root, 'islet_boundaries', False)
    if (root/'provenance').exists():manifest['artifacts']['provenance']=artifact(root, 'provenance')
    manifest['artifacts'].pop('source_geometry', None)
    write_json(root/'manifest.json',manifest)
    capacity=json.loads((root/'capacity.json').read_text())
    capacity['runtime_bytes']=sum(a['bytes'] for a in manifest['artifacts'].values())
    capacity['recommended_two_releases_bytes']=2*capacity['runtime_bytes']
    capacity['candidate_logical_bytes']=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    capacity['runtime_optimization']=report
    write_json(root/'capacity.json',capacity)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--defer-manifest', action='store_true', help='Build during the independent image hash pass; refresh only after that process completes')
    parser.add_argument('--refresh-manifest-only', action='store_true')
    args = parser.parse_args()
    root = args.input_root.resolve(strict=True)
    if not (root/'candidate.json').exists():
        raise ValueError('Explicit candidate identity required')
    if (root/'preview/service.json').exists():raise ValueError('Do not optimize a running candidate')
    if args.defer_manifest and args.refresh_manifest_only:raise ValueError('Conflicting manifest modes')
    if args.refresh_manifest_only:
        report=json.loads((root/'runtime_optimization.json').read_text())
        for name, info in report['tables'].items():
            if sha256(root/'trajectory'/f'{name}.csv')!=info['csv_sha256']:
                raise ValueError('Source table changed after optimization')
        refresh_manifest(root,report)
        return
    start = time.monotonic()
    schema=str(json.loads((root/'candidate.json').read_text()).get('schema_version','3'))
    files = ['composition_features', 'marker_summaries', 'coordinates', 'counts', 'endocrine_counts', 'edges'] if schema=='4' else ['features', 'coordinates', 'counts', 'endocrine_counts', 'edges']
    tables = {}
    for name in files:
        path = root/'trajectory'/f'{name}.csv'
        if not path.exists() and path.with_suffix('.parquet').exists():
            # Schema 4 writers can produce the canonical compressed table directly.
            continue
        table = pd.read_csv(path, dtype={'case_id': str})
        target = path.with_suffix('.parquet')
        if target.exists():
            existing=pd.read_parquet(target)
            pd.testing.assert_frame_equal(table,existing,check_dtype=False)
            continue
        table.to_parquet(target, index=False, compression='zstd')
        tables[name] = dict(rows=len(table), csv_bytes=path.stat().st_size,
                           parquet_bytes=target.stat().st_size, csv_sha256=sha256(path))
    destination = root/'islet_boundaries'
    destination.mkdir(exist_ok=False)
    cohort = pd.read_csv(root/'cohort.csv', dtype=str)
    expected = pd.read_csv(root/'trajectory/coordinates.csv', usecols=['islet_id']).islet_id
    expected = set(expected)
    written = set()
    for cid in cohort.case_id:
        with gzip.open(root/'gson'/f'{int(cid):04d}.geojson.gz', 'rt') as stream:
            collection = json.load(stream)
        for feature in collection['features']:
            props = feature.get('properties', {})
            key = cid+'_'+str(props.get('name', ''))
            if key not in expected or props.get('classification', {}).get('name') != 'Islet':
                continue
            if key in written or '/' in key or '\\' in key:
                raise ValueError('Duplicate or unsafe islet identifier')
            with (destination/(key+'.geojson')).open('w') as stream:
                json.dump(dict(type='FeatureCollection', features=[feature]), stream, separators=(',', ':'), allow_nan=False)
            written.add(key)
    if written != expected:raise ValueError('Missing exact islet outlines: '+str(sorted(expected-written)[:10]))
    report = dict(tables=tables, single_islet_outlines=len(written),
                  outline_bytes=sum(p.stat().st_size for p in destination.iterdir()),
                  seconds=time.monotonic()-start, previous_geometry_not_reused=True)
    write_json(root/'runtime_optimization.json', report)
    # This is a staging-only operation: refresh only the artifacts changed here.
    # The large unchanged images, geometry indexes and cells retain their hashes
    # from this candidate's full independent hash pass.
    if (root/'manifest.json').exists() and not args.defer_manifest:refresh_manifest(root, report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
