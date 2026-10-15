#!/usr/bin/env python3
"""Manifest-driven schema-4 staging and checkpointed endocrine boundary runs.

Only immutable inputs are linked. Derived summaries, memberships, labels and
coordinates are rebuilt. Production prediction on training donors is recorded
as such and is never substituted for the held-out classifier evaluation.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd
import shapely
from shapely.geometry import shape, mapping

from common import MARKERS, read_cohort, require_fresh_output, sha256, write_json

REPO = Path(__file__).resolve().parents[2]


def link_input(source, destination):
    source, destination = Path(source).resolve(strict=True), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if os.path.samefile(source, destination):
            return
        raise FileExistsError(destination)
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def feature_list(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        data = json.load(stream)
    return data['features'] if isinstance(data, dict) else data


def class_name(feature):
    cls = feature.get('properties', {}).get('classification', {})
    return cls.get('name', '') if isinstance(cls, dict) else str(cls)


def stage(args):
    prior = args.input_root.resolve(strict=True)
    out = require_fresh_output(args.output_root, [prior])
    cohort = read_cohort(prior/'source_manifest.csv', REPO)
    out.mkdir(parents=True)
    for folder in ['provenance', 'logs', 'runtime', 'raw', 'corrected', 'typing',
                   'analysis', 'geometry', 'trajectory', 'images', 'boundaries/sources']:
        (out/folder).mkdir(parents=True, exist_ok=True)
    original = REPO/'data/development/endocrine-pixel-refinement-20260917'
    records = []
    for row in cohort.to_dict('records'):
        cid = row['case_id']
        raw = prior/'raw'/f'{cid}-cells.parquet'
        link_input(raw, out/'raw'/raw.name)
        image = Path(row['image'])
        link_input(image, out/'images'/image.name)
        folder = out/'boundaries/sources'/cid
        folder.mkdir()
        reference = original/'references'/cid/'reference_manifest.json'
        if reference.exists():
            record = json.loads(reference.read_text())
            spec = {key: record[key] for key in ['case_id','source_image',
                    'image_shape_cyx','pixel_size_um','source_channel_names']}
            fs = feature_list(record['reference_annotations'])
            tissue = shapely.union_all([shape(f['geometry']) for f in fs if class_name(f)=='Tissue'])
            excludes = [shape(f['geometry']) for f in fs if class_name(f) in ('Exclude','Ignore')]
            if excludes:
                tissue = tissue.difference(shapely.union_all(excludes))
            if tissue.is_empty or not tissue.is_valid:
                raise ValueError(f'{cid}: invalid reference tissue')
            write_json(folder/'valid_tissue.geojson', mapping(tissue))
            spec.update(source_reference=str(reference), source_reference_sha256=sha256(reference),
                        previously_verified_image_sha256=None)
        else:
            source = original/'target_sources'/cid
            spec = json.loads((source/'source.json').read_text())
            shutil.copy2(source/'valid_tissue.geojson', folder/'valid_tissue.geojson')
        source_image = Path(spec['source_image']).resolve(strict=True)
        stat = source_image.stat()
        if 'source_bytes' in spec and (spec['source_bytes']!=stat.st_size or spec['source_mtime_ns']!=stat.st_mtime_ns):
            raise ValueError(f'{cid}: reconciled image fingerprint changed')
        if spec['pixel_size_um'] != .3774 or not set(['INS','GCG','SST']).issubset(spec['source_channel_names']):
            raise ValueError(f'{cid}: unsupported calibration or endocrine channels')
        spec.update(source_image=str(source_image), source_bytes=stat.st_size,
                    source_mtime_ns=stat.st_mtime_ns, full_image_hash_recomputed=False,
                    valid_tissue_sha256=sha256(folder/'valid_tissue.geojson'),
                    SAM_prompts_included=False, cell_gates_included=False)
        write_json(folder/'source.json', spec)
        row.update(raw_cells=str(out/'raw'/raw.name), corrected_cells=str(out/'corrected'/raw.name),
                   image=str(out/'images'/image.name), previous_geometry=row['geometry'],
                   geometry=str(out/'geometry'/f'{cid}.geojson.gz'),
                   boundary_source=str(folder/'source.json'), annotation_review='agent_audit_pending')
        records.append(row)
    result = pd.DataFrame(records)
    result.to_csv(out/'cohort.csv', index=False)
    result.to_csv(out/'source_manifest.csv', index=False)
    shutil.copy2(prior/'donor_metadata.csv', out/'donor_metadata.csv')
    for name in ['candidate.json','manifest.json','revision.json','network_service.json']:
        if (prior/name).exists():
            shutil.copy2(prior/name, out/'provenance'/('previous_'+name))
    subprocess.run(['git','diff','--binary'], cwd=REPO, stdout=(out/'provenance/pre_execution.patch').open('w'), check=True)
    subprocess.run(['git','status','--short'], cwd=REPO, stdout=(out/'provenance/pre_execution_status.txt').open('w'), check=True)
    write_json(out/'candidate.json', dict(schema_version='4',donors=21,
        expected_images=sorted(Path(r['image']).name for r in records),
        cohort_sha256=sha256(out/'cohort.csv'),scientific_status='reprocessing',
        ufrc_accepted=False,research_results_in_app=False))
    write_json(out/'revision.json', dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        parent=str(prior),status='reprocessing',raw_object_preservation=True,
        geometry_policy='all21_frozen_INS_GCG_SST_pixel_classifier_IO60_adaptive',
        cell_typing_policy='donor_balanced_hierarchical_four_ensemble_all31',
        endocrine_denominator='unique validated core alpha/beta/delta/CHGA-supported hormone-unresolved endocrine',
        trajectory_policy='typed_base_cell_counts_per_core_endocrine_cell',
        cpu_budget=56,memory_budget_gib=384,source_project_modified=False))
    print(out, flush=True)


def boundary_job(row, out, model, python, qupath, threads):
    cid = row['case_id']
    prediction = out/'boundaries/predictions'/cid
    adaptive = out/'boundaries/adaptive'/cid
    checkpoint = out/'boundaries/checkpoints'/f'{cid}.json'
    source = Path(row['boundary_source'])
    dependency = dict(source_sha256=sha256(source),model_sha256=sha256(model/'endocrine_forest.json'),
                      helper_sha256=sha256(REPO/'scripts/cohort_expansion/size_adaptive_islet_geometry.groovy'),
                      reference_sha256=sha256(Path.home()/'IO60panc2nd/scripts/1_panc_objects.groovy'),
                      minimum_votes=160)
    if checkpoint.exists():
        previous=json.loads(checkpoint.read_text())
        if previous.get('dependencies')!=dependency:
            raise ValueError(f'{cid}: checkpoint dependency changed; fresh output root required')
        if previous.get('status')=='complete':
            for record in previous['outputs']:
                if sha256(record['path'])!=record['sha256']:
                    raise ValueError(f'{cid}: checkpoint output changed')
            return previous
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    prefs = out/'runtime'/f'java_{cid}'
    prefs.mkdir(parents=True, exist_ok=True)
    env['JAVA_TOOL_OPTIONS']=f'-Xmx12G -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism={threads} -Djava.util.prefs.userRoot={prefs} -Djava.util.prefs.systemRoot={prefs}'
    started=time.monotonic()
    def run(command,label):
        with (out/'logs'/f'boundary_{cid}_{label}.log').open('a') as log:
            log.write(json.dumps(command)+'\n');log.flush()
            subprocess.run(command,cwd=REPO,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    if not (prediction/'prediction_manifest.json').exists():
        run([python,str(REPO/'scripts/cohort_expansion/refine_endocrine_pixel_classifier.py'),
             'predict','--model-root',str(model),'--specification',str(source),
             '--output-root',str(prediction),'--minimum-votes','160','--deployment',
             '--threads',str(threads)],'predict')
    if not (prediction/'vectorization_manifest.json').exists():
        run([qupath,'script','--args',str(prediction),
             str(REPO/'scripts/cohort_expansion/vectorize_endocrine_pixel_regions.groovy')],'trace')
    if not (adaptive/'prediction_manifest.json').exists():
        run([qupath,'script','--args',str(prediction),'--args',str(adaptive),'--args',
             str(REPO/'scripts/cohort_expansion/size_adaptive_islet_geometry.groovy'),'--args',
             str(Path.home()/'IO60panc2nd/scripts/1_panc_objects.groovy'),
             str(REPO/'scripts/cohort_expansion/build_size_adaptive_islets.groovy')],'adaptive')
    outputs=[adaptive/'review_candidates.geojson',adaptive/'prediction_manifest.json',adaptive/'parent_audit.json']
    record=dict(case_id=cid,status='complete',dependencies=dependency,
                outputs=[dict(path=str(p),sha256=sha256(p)) for p in outputs],
                boundary_audit_status='pending',runtime_seconds=time.monotonic()-started)
    write_json(checkpoint,record)
    return record


def boundaries(args):
    out=args.output_root.resolve(strict=True)
    rows=pd.read_csv(out/'cohort.csv',dtype=str).to_dict('records')
    if args.donor:
        rows=[r for r in rows if r['case_id'] in args.donor]
    if args.jobs*args.threads>56:
        raise ValueError('Exceeds the 56-core processing budget')
    failed=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={pool.submit(boundary_job,row,out,args.model_root.resolve(),args.python,args.qupath,args.threads):row['case_id'] for row in rows}
        for f in as_completed(futures):
            cid=futures[f]
            try:
                record=f.result();print(f"{cid}: boundaries complete ({record['runtime_seconds']:.0f}s)",flush=True)
            except Exception as exc:
                failed.append(dict(case_id=cid,error=str(exc)));print(f'{cid}: FAILED {exc}',flush=True)
    write_json(out/'boundaries/run_status.json',dict(requested=len(rows),completed=len(rows)-len(failed),failures=failed))
    if failed:
        raise SystemExit(1)


def refine_all(args):
    """Reuse frozen pixel predictions after a geometry-only implementation fix."""
    root=args.output_root.resolve(strict=True)
    rows=pd.read_csv(root/'cohort.csv',dtype=str).to_dict('records')
    def one(row):
        cid=row['case_id'];out=Path(row['boundary_candidates']).parent
        source=root/'boundaries/predictions'/cid
        evidence=root/'boundaries/assessability'/cid/'exclusions.geojson'
        expected=dict(prediction_sha256=sha256(source/'prediction_manifest.json'),
            geometry_helper_sha256=sha256(REPO/'scripts/cohort_expansion/size_adaptive_islet_geometry.groovy'),
            technical_exclusions_sha256=sha256(evidence) if evidence.exists() else None)
        check=out/'rebuild_checkpoint.json'
        if check.exists():
            previous=json.loads(check.read_text())
            if previous['dependencies']!=expected:raise ValueError(f'{cid}: geometry rebuild dependencies changed')
            return cid
        prefs=root/'runtime'/f'adaptive_v2_{cid}';prefs.mkdir(parents=True,exist_ok=True)
        env=os.environ.copy();env['JAVA_TOOL_OPTIONS']=f'-Xmx8G -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism=1 -Djava.util.prefs.userRoot={prefs} -Djava.util.prefs.systemRoot={prefs}'
        command=[args.qupath,'script','--args',str(source),'--args',str(out),'--args',
            str(REPO/'scripts/cohort_expansion/size_adaptive_islet_geometry.groovy'),'--args',
            str(Path.home()/'IO60panc2nd/scripts/1_panc_objects.groovy')]
        if evidence.exists():command.extend(['--args',str(evidence)])
        command.append(str(REPO/'scripts/cohort_expansion/build_size_adaptive_islets.groovy'))
        with (root/'logs'/f'adaptive_v2_{cid}.log').open('w') as log:
            subprocess.run(command,cwd=REPO,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        sys.path.insert(0,str(REPO/'scripts/cohort_expansion'))
        from verify_size_adaptive_output import verify
        result=verify(out)
        write_json(check,dict(dependencies=expected,geometry_check=result,status='complete'))
        return cid
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for f in as_completed([pool.submit(one,r) for r in rows]):print('Adaptive rebuilt',f.result(),flush=True)


def export_native(args):
    root=args.output_root.resolve(strict=True);rows=pd.read_csv(root/'cohort.csv',dtype=str).to_dict('records')
    def one(row):
        cid=row['case_id'];out=Path(row['native_core_root'])
        if (out/'manifest.json').exists():
            previous=json.loads((out/'manifest.json').read_text())
            if previous['source_candidates_sha256']!=sha256(row['boundary_candidates']):raise ValueError('Native candidates changed')
            return cid
        prefs=root/'runtime'/f'native_v2_{cid}';prefs.mkdir(parents=True,exist_ok=True)
        env=os.environ.copy();env['JAVA_TOOL_OPTIONS']=f'-Xmx5G -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism=1 -Djava.util.prefs.userRoot={prefs} -Djava.util.prefs.systemRoot={prefs}'
        command=[args.qupath,'script','--args',row['boundary_candidates'],'--args',
            str(root/'native_cells'/cid/'native_centroids.tsv'),'--args',cid,'--args',str(out),
            str(REPO/'scripts/cohort_expansion/export_regenerated_core_geometry.groovy')]
        with (root/'logs'/f'native_v2_{cid}.log').open('w') as log:
            subprocess.run(command,cwd=REPO,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        return cid
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for f in as_completed([pool.submit(one,r) for r in rows]):print('Native rebuilt',f.result(),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('stage');s.add_argument('--input-root',type=Path,required=True);s.set_defaults(func=stage)
    b=sub.add_parser('boundaries');b.add_argument('--model-root',type=Path,required=True)
    b.add_argument('--python',default=sys.executable);b.add_argument('--qupath',default='/home/smith6jt/QuPath/bin/QuPath')
    b.add_argument('--jobs',type=int,default=7);b.add_argument('--threads',type=int,default=4)
    b.add_argument('--donor',action='append');b.set_defaults(func=boundaries)
    c=sub.add_parser('refine');c.set_defaults(func=refine_all)
    d=sub.add_parser('export-native');d.set_defaults(func=export_native)
    for parser in [c,d]:
        parser.add_argument('--jobs',type=int,default=7)
        parser.add_argument('--qupath',default='/home/smith6jt/QuPath/bin/QuPath')
    for parser in [s,b,c,d]:parser.add_argument('--output-root',type=Path,required=True)
    args=p.parse_args();args.func(args)


if __name__=='__main__':main()
