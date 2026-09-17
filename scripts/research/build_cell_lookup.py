#!/usr/bin/env python3
"""Raw marker drill-down for all donors; do not invent unvalidated phenotypes."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import MARKERS, write_json
from build_candidate import geometry_data, membership

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);a=p.parse_args();root=Path(a.input_root);out=root/'cells';out.mkdir(exist_ok=True)
    donors=pd.read_csv(root/'cohort.csv',dtype=str);records=[]
    for cid in donors.case_id:
        cells=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet')
        cores,geoms,_=geometry_data(root/'gson'/f'{int(cid):04d}.geojson.gz')
        tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()));union=shapely.union_all([g for _,g in cores if g.is_valid])
        x=cells.X_centroid.to_numpy()/.3774;y=cells.Y_centroid.to_numpy()/.3774;index=shapely.STRtree(shapely.points(x,y))
        fields=['object_id','X_centroid','Y_centroid','marker_measurements_complete','Cell Area','Nucleus Area']+MARKERS
        total=0
        for name,g in cores:
            if not g.is_valid: continue
            rows=[]
            for comp,region in [('core',g.intersection(tissue)),('peri',g.buffer(20/.3774).difference(union).intersection(tissue))]:
                ids=membership(index,region,x,y);frame=cells.iloc[ids][fields].copy()
                frame['cell_region']=comp;frame['phenotype']=pd.NA
                rows.append(frame)
            table=pd.concat(rows,ignore_index=True);total+=len(table)
            table.to_csv(out/f'{cid}_{name}.csv',index=False)
        records.append(dict(case_id=cid,islets=len(cores),cell_membership_rows=total))
        print(cid,'cell lookup complete',flush=True)
    write_json(root/'cell_lookup.json',dict(donors=records,units='raw Cell mean intensity; centroid micrometres',
               phenotype='unavailable pending full parent gate validation; archival labels remain in raw source files',
               regions='tissue-clipped core and 20 um peri; neighboring cores excluded from peri',unmeasured='preserved as NaN'))

if __name__=='__main__':main()
