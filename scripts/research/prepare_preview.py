#!/usr/bin/env python3
"""Stage all selected images and existing regional measurements before analysis."""
import argparse
from pathlib import Path
import shutil
import pandas as pd
import anndata as ad
from common import read_cohort, write_json, sha256
from build_candidate import materialize

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--output-root',required=True);p.add_argument('--cohort',required=True);a=p.parse_args()
    root=Path(a.input_root);out=Path(a.output_root);out.mkdir(parents=True,exist_ok=True)
    cohort=read_cohort(a.cohort,root);meta=[];tables={k:[] for k in ['composition','targets','markers']}
    baselines={path:ad.read_h5ad(path,backed='r') for path in set(cohort.loc[cohort.cell_slot=='raw','cells'])}
    for row in cohort.itertuples(index=False):
        m=pd.read_csv(row.metadata,dtype={'case_id':str}).query('case_id == @row.case_id').copy()
        if len(m)!=1: raise ValueError('Metadata join must be one to one')
        if row.cell_slot=='raw':
            baseline=baselines[row.cells]
            old=baseline.obs.loc[baseline.obs.imageid.astype(str)==row.case_id].iloc[0];m['age']=float(old.Age);m['sex']=old.Gender
        else: m['age']=m.age_years_key
        m['donor_source']=row.donor_source;m['donor_status']=row.donor_status;m['inclusion_status']='selected';meta.append(m)
        img=Path(row.image)
        for src in [img,img.with_suffix('.xml'),Path(str(img)+'.offsets.json')]:
            if src.exists(): materialize(src,out/'images'/src.name)
        materialize(row.geometry,out/'gson'/f'{int(row.case_id):04d}.geojson.gz')
        for kind in tables:
            src=Path(row.regional_root)/f'islet_{kind}_{int(row.case_id):04d}.tsv'
            if not src.is_file(): raise ValueError(f'Missing selected regional export {src}')
            d=pd.read_csv(src,sep='\t');d['Case ID']=row.case_id;d['Donor Status']=row.donor_status;d['Age']=float(m.age.iloc[0]);d['Gender']=m.sex.iloc[0]
            for col in ['AAb_GADA','AAb_IA2A','AAb_mIAA','AAb_ZnT8A']:
                d[col]=m[col].iloc[0] if col in m else None
            d['source_file']=src.name;tables[kind].append(d)
            materialize(src,out/'regional'/src.name)
    for baseline in baselines.values():baseline.file.close()
    pd.concat(meta,ignore_index=True).to_csv(out/'donor_metadata.csv',index=False)
    with pd.ExcelWriter(out/'master_results.xlsx',engine='openpyxl') as w:
        for kind,frames in tables.items(): pd.concat(frames,ignore_index=True).to_excel(w,sheet_name={'composition':'Islet_Composition','targets':'Islet_Targets','markers':'Islet_Markers'}[kind],index=False)
    shutil.copyfile(root/'data/phenotype_rules.csv',out/'phenotype_rules.csv')
    shutil.copyfile(a.cohort,out/'cohort.csv')
    (out/'trajectory').mkdir(exist_ok=True)
    if not (out/'trajectory/validation.json').exists(): write_json(out/'trajectory/validation.json',{'status':'building','interpretation':'Offline cohort measurements and trajectory validation are running.'})
    write_json(out/'candidate.json',dict(schema_version='3',cohort_sha256=sha256(out/'cohort.csv'),donors=21,
          groups={'ND':7,'Aab+':7,'T1D':7},scientific_status='validation_in_progress',ufrc_accepted=False,
          expected_images=sorted(Path(p).name for p in cohort.image)))

if __name__=='__main__':main()
