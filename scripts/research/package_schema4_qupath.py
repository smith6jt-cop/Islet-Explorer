#!/usr/bin/env python3
"""Package native cell UUIDs, regenerated geometry and schema4 typing for QuPath."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import subprocess
import pandas as pd
from common import MARKERS,sha256,write_json

REPO=Path(__file__).resolve().parents[2]


def donor(root,out,row,entry_id):
    cid=row['case_id'];folder=out/'data'/str(entry_id);folder.mkdir(parents=True,exist_ok=True)
    if (folder/'data.qpdata').exists():raise FileExistsError('Preserve native project edits')
    raw_native=json.loads((root/'native_cells'/cid/'manifest.json').read_text())
    native_root=Path(row['native_core_root']);proof=json.loads((native_root/'manifest.json').read_text())
    cells=pd.read_parquet(root/'typing/cells'/f'{cid}.parquet');corrected=pd.read_parquet(root/'corrected'/f'{cid}-cells.parquet')
    for f in [cells,corrected]:f['case_id']=f.case_id.astype(str)
    keys=['case_id','object_id']
    if len(cells)!=len(corrected) or len(cells)!=int(raw_native['cells']):raise ValueError('Native/corrected/typing universes differ')
    frame=cells.merge(corrected[keys+MARKERS],on=keys,validate='one_to_one')
    if len(frame)!=len(cells):raise ValueError('Unmatched native typing records')
    frame['native_uuid']=frame.object_id.str.split(':',n=1).str[-1]
    frame['Typing: endocrine member']=frame.is_endocrine.astype('Float64')
    for col in ['base_type','lineage','phenotype_state']:
        if col not in frame:frame[col]=''
    renames={m:'Corrected: '+m for m in MARKERS}
    for m in MARKERS:
        col='positive__'+m
        if col in frame:
            frame[col]=frame[col].astype('Float64');renames[col]='Marker positive: '+m
    frame=frame.rename(columns=renames)
    cols=['native_uuid','classification_status','base_type','lineage','phenotype_state','Typing: endocrine member']+list(renames.values())
    table=out/'inputs'/f'{cid}-cells.tsv';table.parent.mkdir(exist_ok=True)
    frame[cols].to_csv(table,index=False,sep='\t',na_rep='')
    job=dict(case_id=cid,candidate_root=str(root),image=row['image'],source_native=raw_native['source_native'],
             source_pixel_size_um=json.loads(Path(row['boundary_source']).read_text())['pixel_size_um'],
             source_sha256=raw_native['source_sha256'],cells=len(frame),islets=proof['cores'],
             geometry=row['geometry'],geometry_sha256=sha256(row['geometry']),
             core_memberships=str(native_root/'native_core_memberships.csv'),
             cell_table=str(table),cell_table_sha256=sha256(table),output_data=str(folder/'data.qpdata'))
    job_file=out/'inputs'/f'{cid}.json';write_json(job_file,job)
    env=os.environ.copy();prefs=out/'runtime'/cid;prefs.mkdir(parents=True,exist_ok=True)
    env['JAVA_TOOL_OPTIONS']=f'-Xmx12G -XX:ActiveProcessorCount=2 -Djava.awt.headless=true -Djava.util.concurrent.ForkJoinPool.common.parallelism=2 -Djava.util.prefs.userRoot={prefs} -Djava.util.prefs.systemRoot={prefs}'
    log=out/'logs'/f'{cid}.log';log.parent.mkdir(exist_ok=True)
    with log.open('w') as stream:
        subprocess.run(['/home/smith6jt/QuPath/bin/QuPath','script','--args',str(job_file),
                        str(REPO/'scripts/cohort_expansion/build_schema4_review.groovy')],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
    result=json.loads((folder/'verification.json').read_text());print(cid,result['status'],flush=True)
    return dict(entryID=entry_id,randomizedName=f'Schema4_{cid}',imageName=Path(row['image']).name,
                serverBuilder=json.loads((folder/'server.json').read_text()),
                metadata={'Donor':cid,'Disease status':row['donor_status'],'Donor source':row['donor_source'],
                          'Candidate':str(root),'Geometry':'Regenerated endocrine-pixel boundaries; sample image audit with region flags',
                          'Cells':'Native UUID/ROI/raw means retained; corrected means and new identities added'})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    root=a.input_root.resolve(strict=True);out=a.output_root.resolve()
    if out.exists():raise FileExistsError('Use a fresh independent native review project')
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);out.mkdir(parents=True)
    entries=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs=[pool.submit(donor,root,out,row,i+1) for i,row in enumerate(cohort.to_dict('records'))]
        for f in as_completed(jobs):entries.append(f.result())
    now=int(datetime.now(timezone.utc).timestamp()*1000)
    write_json(out/'project.qpproj',dict(version='0.7.0',createTimestamp=now,modifyTimestamp=now,
               uri=(out/'project.qpproj').as_uri(),lastID=len(entries),images=sorted(entries,key=lambda x:x['entryID'])))
    write_json(out/'classifiers/classes.json',dict(pathClasses=[]))
    checks=[json.loads((out/'data'/str(e['entryID'])/'verification.json').read_text()) for e in entries]
    write_json(out/'review_manifest.json',dict(status='all_native_roundtrips_verified',donors=len(entries),
               cells=sum(r['cells'] for r in checks),islets=sum(r['islets'] for r in checks),
               source_candidate=str(root),cohort_sha256=sha256(root/'cohort.csv'),
               project_sha256=sha256(out/'project.qpproj'),source_projects_modified=False,
               user_review_claimed=False,verification=checks))
    print(out/'project.qpproj',flush=True)
if __name__=='__main__':main()
