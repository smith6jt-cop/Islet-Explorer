#!/usr/bin/env python3
"""Raw marker drill-down for all donors; do not invent unvalidated phenotypes."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import MARKERS, write_json, sha256
from build_candidate import geometry_data, membership

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--case',action='append');a=p.parse_args();root=Path(a.input_root);out=root/'cells';out.mkdir(exist_ok=True)
    donors=pd.read_csv(root/'cohort.csv',dtype=str)
    if a.case:
        if not set(a.case).issubset(donors.case_id):raise ValueError('Donor outside explicit cohort')
        donors=donors[donors.case_id.isin(a.case)]
    previous=json.loads((root/'cell_lookup.json').read_text()) if (root/'cell_lookup.json').exists() else {}
    records=[r for r in previous.get('donors',[]) if r['case_id'] not in set(donors.case_id)]
    for cid in donors.case_id:
        cells=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet')
        geometry=root/'gson'/f'{int(cid):04d}.geojson.gz'
        provenance=root/'provenance'/f'{cid}-native-geometry.json'
        native=json.loads(provenance.read_text()) if provenance.exists() else {}
        if native and sha256(geometry)!=native['geometry_sha256']:raise ValueError('Native geometry hash mismatch')
        roi_types=native.get('core_roi_type_by_name',{})
        cores,geoms,_=geometry_data(geometry)
        tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()));union=shapely.union_all([g for _,g in cores if g.is_valid])
        x=cells.X_centroid.to_numpy()/.3774;y=cells.Y_centroid.to_numpy()/.3774;index=shapely.STRtree(shapely.points(x,y))
        fields=['object_id','X_centroid','Y_centroid','marker_measurements_complete','Cell Area','Nucleus Area']+MARKERS
        total=core_total=0
        for name,g in cores:
            if not g.is_valid: continue
            rows=[]
            for comp,region in [('core',g.intersection(tissue)),('peri',g.buffer(20/.3774).difference(union).intersection(tissue))]:
                policy='qupath_geometry' if comp=='core' and roi_types.get(name)=='GeometryROI' else 'qupath_polygon'
                ids=membership(index,region,x,y,boundary_policy=policy);frame=cells.iloc[ids][fields].copy()
                if comp=='core':core_total+=len(ids)
                frame['cell_region']=comp;frame['phenotype']=pd.NA
                rows.append(frame)
            table=pd.concat(rows,ignore_index=True);total+=len(table)
            target=out/f'{cid}_{name}.csv'
            if target.exists():raise ValueError('Refusing to overwrite existing cell artifact: '+str(target))
            table.to_csv(target,index=False)
        records.append(dict(case_id=cid,islets=len(cores),cell_membership_rows=total,core_membership_rows=core_total,geometry_sha256=sha256(geometry)))
        print(cid,'cell lookup complete',flush=True)
    write_json(root/'cell_lookup.json',dict(donors=records,units='raw Cell mean intensity; centroid micrometres',
               phenotype='unavailable pending full parent gate validation; archival labels remain in raw source files',
               regions='tissue-clipped core and 20 um peri; neighboring cores excluded from peri',unmeasured='preserved as NaN'))

if __name__=='__main__':main()
