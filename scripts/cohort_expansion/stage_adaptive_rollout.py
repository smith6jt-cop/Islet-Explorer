#!/usr/bin/env python3
"""Stage six independent native copies and bind the owner's HPAP boundary approval.

The approved review hierarchy is copied exactly, including any saved edits.
Approval applies only to its Endocrine pixel candidate objects, not phenotypes
or automatically generated objects in other donors.
"""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

from prepare_review_project import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--approved-review',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--approval-text',required=True)
    args=p.parse_args()
    reviewed=args.approved_review.resolve(strict=True)
    out=args.output_root.resolve()
    if out.exists():raise ValueError('Use a fresh project; preserve existing review edits')
    source_manifest=json.loads((reviewed/'review_manifest.json').read_text())
    if len(source_manifest['donors'])!=1 or source_manifest['donors'][0]['case_id']!='78':
        raise ValueError('This recorded visual approval is specifically for HPAP-078')
    native=reviewed/'data/42/data.qpdata'
    native_hash=digest(native)
    project_document=json.loads((reviewed/'project.qpproj').read_text())
    reviewed_entry=copy.deepcopy(next(e for e in project_document['images'] if e['entryID']==42))
    project_hash=digest(reviewed/'project.qpproj')
    subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_review_project.py')),
                    '--source-root',str(args.source_root.resolve()),'--output-root',str(out)],check=True)
    staged=json.loads((out/'review_manifest.json').read_text())
    donor=next(r for r in staged['donors'] if r['case_id']=='78')
    donor['processing_source_data']=donor['source_data']
    donor['processing_source_data_sha256']=donor['data_sha256']
    for name in ['data.qpdata','server.json','summary.json']:
        shutil.copy2(reviewed/'data/42'/name,out/'data/42'/name)
    if digest(out/'data/42/data.qpdata')!=native_hash or digest(native)!=native_hash:
        raise ValueError('The reviewed hierarchy changed during staging')
    donor.update(source_data=str(native),data_sha256=native_hash,data_bytes=native.stat().st_size,
                 annotation_review='owner_approved_endocrine_pixel_boundaries')
    project=json.loads((out/'project.qpproj').read_text())
    project['images']=[reviewed_entry if e['entryID']==42 else e for e in project['images']]
    for entry in project['images']:
        entry.setdefault('metadata',{})['Boundary review']='Owner-approved HPAP-078 boundaries; promotion pending' if entry['entryID']==42 else 'Fixed protocol transfer pending; donor review required'
    (out/'project.qpproj').write_text(json.dumps(project,indent=2)+'\n')
    (out/'staging_manifest.json').write_text(json.dumps(staged,indent=2)+'\n')
    # Keep the preparation manifest as evidence of the original six-source copy.
    approval=dict(version=1,recorded_at_utc=datetime.now(timezone.utc).isoformat(),case_id='78',
        decision='owner_visually_approved',scope='Saved Endocrine pixel candidate boundary objects in HPAP-078 only',
        user_statement=args.approval_text,source_review_project=str(reviewed/'project.qpproj'),
        source_review_project_sha256=project_hash,source_review_data=str(native),source_review_data_sha256=native_hash,
        source_review_object_class='Endocrine pixel candidate',
        includes_phenotype_gate_validation=False,includes_other_donor_visual_approval=False,
        original_review_project_modified=False,staged_project=str(out/'project.qpproj'))
    (out/'boundary_approval.json').write_text(json.dumps(approval,indent=2)+'\n')
    print(json.dumps(approval),flush=True)


if __name__=='__main__':main()
