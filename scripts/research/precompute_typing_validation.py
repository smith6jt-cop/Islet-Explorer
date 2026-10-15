#!/usr/bin/env python3
"""Checkpoint label-independent validation discovery while image review proceeds.

This does not produce validated identities or complete reconstruction evidence.
The hierarchical validation runner must subsequently reconstruct refinements and
score held-out predictions. Only the broad, unsupervised fit is precomputed.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
from cell_typing import TypingConfig,discover,canonical_hash
from cell_typing_validation import validation_variants
from common import write_json,sha256
import pandas as pd


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True);p.add_argument('--budget',type=Path,required=True)
    p.add_argument('--only',action='append',default=[]);a=p.parse_args()
    root=a.input_root.resolve(strict=True);out=a.output_root.resolve();out.mkdir(parents=True,exist_ok=True)
    specification=json.loads((root/'typing/model/specification.json').read_text())
    config=TypingConfig(**specification['config']);donors=pd.read_csv(root/'cohort.csv',dtype=str).case_id.tolist()
    variants=validation_variants(donors,config)
    if a.only:variants=[v for v in variants if v[0] in a.only]
    record=dict(status='running',source_model_identity=specification['identity'],
                source_specification_sha256=sha256(root/'typing/model/specification.json'),
                scope='Broad unsupervised discovery only; does not complete hierarchical validation',completed=[],failed=[])
    for name,variant,excluded,_ in variants:
        while True:
            budget=json.loads(a.budget.read_text());workers=int(budget['workers'])
            if budget.get('stop_after_current',False):
                record['status']='stopped_after_checkpoint';write_json(out/'precompute_manifest.json',record);return
            if workers>0:break
            time.sleep(10)
        if workers>40:raise ValueError('Coordinate CPU budget before exceeding40 typing workers')
        started=time.monotonic();print(f'Start {name}; {workers} workers',flush=True)
        try:
            result=discover(root,out/name,variant,workers,min(4,workers),exclude_donors=excluded)
            record['completed'].append(dict(variant=name,model_identity=result['identity'],
                                           config_hash=canonical_hash(asdict(variant)),workers=workers,seconds=time.monotonic()-started))
            print(f'Finished {name}; {time.monotonic()-started:.1f} seconds',flush=True)
        except Exception as error:
            record['failed'].append(dict(variant=name,error=str(error),error_type=type(error).__name__))
            print(f'FAILED {name}: {error}',flush=True)
        write_json(out/'precompute_manifest.json',record)
    record['status']='broad_discovery_complete' if not record['failed'] else 'broad_discovery_incomplete'
    write_json(out/'precompute_manifest.json',record)
if __name__=='__main__':main()
