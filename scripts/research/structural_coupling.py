#!/usr/bin/env python3
"""Area-weighted excess nerve proximity under a compartment uniform-location null.

The null expectation is analytic: P(distance to capillary <=20 um | region)
equals the capillary-buffer area fraction of that region. It fixes the region
and capillary map and conditions on nerve abundance, but does not preserve nerve
autocorrelation. This limitation is explicit in every saved result.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import read_cohort, write_json, offline_report_root
from build_candidate import geometry_data

def clip_parts(region, entry):
    if entry is None: return shapely.GeometryCollection()
    parts,tree=entry;ids=tree.query(region,predicate='intersects')
    return shapely.union_all(shapely.intersection(parts[ids],region)) if len(ids) else shapely.GeometryCollection()

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--output-root',required=True);p.add_argument('--report-root',required=True);p.add_argument('--cohort',required=True);a=p.parse_args()
    cohort=read_cohort(a.cohort,a.input_root);out=Path(a.output_root);reports=offline_report_root(out,a.report_root);rows=[]
    assessed=pd.read_parquet(out/'analysis/regions.parquet').query('geometry_valid and not overlap and annotation_review != "pending"')
    selected=set(assessed.islet_id)
    for row in cohort.itertuples(index=False):
        cores,geoms,_=geometry_data(row.geometry);tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
        allcores=shapely.union_all([g for _,g in cores if g.is_valid]);px=.3774
        for name,g in cores:
            if not g.is_valid or f'{row.case_id}_{name}' not in selected: continue
            for band in [10,20,30]:
                regions={'core':g.intersection(tissue),'peri':g.buffer(band/px).difference(allcores).intersection(tissue)}
                for compartment,region in regions.items():
                    nerve=clip_parts(region,geoms.get('Nerve'))
                    capillary=clip_parts(region.buffer(20/px),geoms.get('Capillary'))
                    near=region.intersection(capillary.buffer(20/px))
                    observed=nerve.intersection(near).area/nerve.area if nerve.area>0 else np.nan
                    expected=near.area/region.area if region.area>0 else np.nan
                    rows.append(dict(islet_id=f'{row.case_id}_{name}',case_id=row.case_id,band_um=band,compartment=compartment,
                                     nerve_proximity_observed=observed,nerve_proximity_null=expected,excess_proximity=observed-expected))
        print(f'{row.case_id}: coupling complete',flush=True)
    pd.DataFrame(rows).to_parquet(out/'analysis/coupling.parquet',index=False)
    write_json(reports/'spatial_null.json',dict(method='analytic uniform-location null conditional on compartment and observed structural area',
               distance_um=20,limitation='Does not preserve nerve autocorrelation; spatially constrained shape-preserving randomization remains a required sensitivity before biological acceptance.'))

if __name__=='__main__':main()
