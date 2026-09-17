#!/usr/bin/env python3
"""Hash the new runtime and original input provenance; no deployment attestation."""
import argparse
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq
from common import sha256, write_json, require_cohort

def artifact(root,relative,startup=True):
    path=root/relative
    if path.is_dir():
        files=sorted(p for p in path.rglob('*') if p.is_file())
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=4) as pool:
            records=list(pool.map(lambda p:dict(path=str(p.relative_to(path)),bytes=p.stat().st_size,sha256=sha256(p)),files))
        import hashlib
        digest=hashlib.sha256('\n'.join(r['path']+'  '+r['sha256'] for r in records).encode()).hexdigest()
        return dict(relative_path=relative,type='tree',sha256=digest,bytes=sum(r['bytes'] for r in records),files=records,startup_verify=startup)
    result=dict(relative_path=relative,type='file',sha256=sha256(path),bytes=path.stat().st_size,startup_verify=startup)
    if relative.endswith('.parquet'):result['row_count']=pq.read_metadata(path).num_rows
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--schema',choices=['2','3'],default='3');p.add_argument('--add-cell-artifact',action='store_true');a=p.parse_args();root=Path(a.input_root).resolve(strict=True)
    if a.add_cell_artifact:
        # Append a newly completed optional artifact to this candidate's already
        # verified manifest. Other inputs must remain immutable; this is not a
        # general refresh operation and never attests to changed existing files.
        manifest=json.loads((root/'manifest.json').read_text())
        if 'cells' in manifest['artifacts']:raise ValueError('Cell artifact is already declared')
        manifest['artifacts']['cells']=artifact(root,'cells',False)
        write_json(root/'manifest.json',manifest)
        capacity=json.loads((root/'capacity.json').read_text())
        capacity['runtime_bytes']+=manifest['artifacts']['cells']['bytes']
        capacity['recommended_two_releases_bytes']=2*capacity['runtime_bytes']
        capacity['candidate_logical_bytes']=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        write_json(root/'capacity.json',capacity)
        return
    cohort=pd.read_csv(root/'cohort.csv',dtype=str)
    if len(cohort)!=21:raise ValueError('Wrong cohort')
    for path in ['trajectory/features.csv','trajectory/coordinates.csv','donor_metadata.csv']:
        require_cohort(pd.read_csv(root/path,dtype={'case_id':str}),cohort.case_id)
    mapping={'core_area':'summary/core_area.parquet','targets_all':'summary/targets.parquet','markers_all':'summary/markers.parquet',
             'comp':'summary/islets.parquet','combined_pairs':'summary/combined_pairs.parquet','geometry':'geometry',
             'islet_spatial_lookup':'islet_spatial_lookup.csv','phenotype_rules':'phenotype_rules.csv','donor_metadata':'donor_metadata.csv',
             'images':'images','cohort':'cohort.csv','trajectory':'trajectory','cells':'cells'}
    if a.schema=='2':mapping['research']='research'
    else:
        for name in ['counts.csv','endocrine_counts.csv','phenotype_availability.csv','marker_panel.json','count_manifest.json']:
            if not (root/'trajectory'/name).is_file():raise ValueError('Missing schema-3 trajectory artifact: '+name)
        if (root/'research').exists():raise ValueError('Research reports belong outside the schema-3 runtime')
    artifacts={}
    for name,rel in mapping.items():
        print('Hashing',name,flush=True);artifacts[name]=artifact(root,rel,name not in ['images','geometry','cells'])
    sources=pd.read_csv(root/'source_manifest.csv',dtype=str)
    paths=sorted(set(sources.cells)|set(sources.geometry)|set(sources.metadata))
    evidence=[dict(path=str(p),sha256=sha256(p),bytes=Path(p).stat().st_size) for p in paths]
    unavailable={'legacy_spatial_cell_distances':'not rebuilt for the expanded cohort; no frozen fallback'}
    if a.schema=='2':unavailable.update(phenotype_ratios='upstream gate validation incomplete; no certified denominator',hypotheses='per-test status in research/results.csv')
    else:
        counts=json.loads((root/'trajectory/count_manifest.json').read_text())
        if counts['status']!='available':unavailable['phenotype_ratios']='see trajectory/phenotype_availability.csv and endocrine_counts.csv for donor/type-specific missing evidence'
        if (root/'trajectory/quality_review.json').exists():unavailable['islet_identity']='candidate islet identity requires review; existing coordinates are provisional'
    manifest=dict(schema_version=a.schema,created_at_utc=datetime.now(timezone.utc).isoformat(),builder='scripts/research/finalize_candidate.py',
                  scientific_acceptance=False,ufrc_accepted=False,expected_images=21,artifacts=artifacts,canonical_sources=evidence,
                  unavailable_analyses=unavailable,
                  research_results_in_app=False,
                  materialization_mode='hardlink raw inputs/images; newly written derived data; never mutate linked inputs')
    write_json(root/'manifest.json',manifest)
    runtime_bytes=sum(x['bytes'] for x in artifacts.values())
    usage=dict(runtime_bytes=runtime_bytes,images_bytes=artifacts['images']['bytes'],
               candidate_logical_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),
               recommended_two_releases_bytes=2*runtime_bytes,
               cpu=platform.processor(),trajectory=json.loads((root/'trajectory/validation.json').read_text()),
               analysis_seconds=sum(d['seconds'] for d in json.loads((root/'build_report.json').read_text())['donor_reports']),
               note='Logical sizes include hardlinked source files; production copy needs full runtime bytes. No image/container acceptance is implied.')
    write_json(root/'capacity.json',usage)

if __name__=='__main__':main()
