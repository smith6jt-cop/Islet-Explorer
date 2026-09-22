#!/usr/bin/env python3
"""Materialize empty schema-4 trajectory tables for an unavailable fit.

The runtime requires stable table contracts even when the scientific fit is
blocked upstream. Empty coordinates/edges/layout make the app render the
explicit unavailable state without fabricating an axis or reusing an old fit.
"""
from pathlib import Path
import argparse, json
import pandas as pd

def main():
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,required=True); p.add_argument('--validation',type=Path,required=True)
    a=p.parse_args(); root=a.root; tr=root/'trajectory'; tr.mkdir(parents=True,exist_ok=True)
    status=json.loads(a.validation.read_text())
    status['coordinates_available']=False
    status['interpretation']='No outcome-blind trajectory axis was fit because the primary typed core/peri eligibility gate is unavailable for part of the cohort.'
    (tr/'validation.json').write_text(json.dumps(status,indent=2,sort_keys=True)+'\n')
    features=pd.read_parquet(tr/'composition_features.parquet')
    coord_cols=['islet_id','case_id','donor_status','donor_source','islet_key','area_um2','measured_cells','eligible','exclusion_reason']
    coords=features[[c for c in coord_cols if c in features]].copy()
    coords['component']=-1; coords['supported']=False; coords['support_distance']=float('nan')
    coords['pseudotime_raw']=float('nan'); coords['pseudotime_percentile']=float('nan')
    coords.to_csv(tr/'coordinates.csv',index=False)
    pd.DataFrame(columns=['from','to','weight']).to_csv(tr/'edges.csv',index=False)
    pd.DataFrame({'islet_id':coords['islet_id'],'component':-1,'x':float('nan'),'y':float('nan')}).to_csv(tr/'graph_layout.csv',index=False)

if __name__=='__main__': main()
