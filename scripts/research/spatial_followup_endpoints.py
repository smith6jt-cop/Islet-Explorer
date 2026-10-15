#!/usr/bin/env python3
"""Exact-membership, positive-mask endpoints for lymphatic follow-up analyses."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from joblib import Parallel,delayed
from common import sha256,write_json
from composition import canonical_cells
from build_candidate import geometry_data
from structural_coupling import clip_parts

SUBSET_ALIASES={'CD8 T cell':'CD8','CD8 Tcell':'CD8','CD8a Tcell':'CD8',
                'CD4 T cell':'CD4','CD4 Tcell':'CD4','B cell':'B','Macrophage':'macrophage',
                'Macrophage-associated cell':'macrophage'}


def subset_mapping(catalog):
    lookup={}
    for entry in catalog['base_types']:
        subset=entry.get('immune_subset') or SUBSET_ALIASES.get(entry['name'])
        if subset is not None:
            if subset not in {'CD8','CD4','B','macrophage'}:raise ValueError('Unknown immune subset catalog mapping')
            lookup[entry['name']]=subset
    return lookup


def donor_endpoints(root,row,catalog,regions,memberships):
    cid=str(row['case_id']);directory=root/'analysis/spatial_followup_parts'/cid;directory.mkdir(parents=True,exist_ok=True)
    typ=root/'typing/cells'/f'{cid}.parquet';raw=root/'corrected'/f'{cid}-cells.parquet'
    fingerprint=hashlib.sha256((sha256(typ)+sha256(raw)+sha256(row['geometry'])+sha256(root/'typing/catalog.json')+
        sha256(root/'analysis/regions.parquet')+sha256(root/'analysis/memberships.parquet')+sha256(__file__)).encode()).hexdigest()
    checkpoint=directory/'complete.json'
    if checkpoint.exists() and json.loads(checkpoint.read_text())['fingerprint']==fingerprint:
        return pd.read_parquet(directory/'neighborhoods.parquet'),pd.read_parquet(directory/'localization.parquet')
    r=regions.loc[regions.case_id.eq(cid)&regions.band_um.eq(20)&regions.geometry_valid&~regions.overlap].copy()
    if 'primary_assessable' in r:r=r.loc[r.primary_assessable.fillna(False)]
    if 'annotation_status' in r:r=r.loc[r.annotation_status.ne('technical_region_unassessable')]
    required={'lymphatic_validated','pixel_size_um'}
    if not required.issubset(r):raise ValueError('Regional lymphatic validation and physical calibration required')
    r=r.loc[r.lymphatic_validated.fillna(False)]
    calibrations=r.pixel_size_um.dropna().unique()
    if len(calibrations)>1:raise ValueError('Mixed physical calibration within donor')
    neighborhoods=[];localization=[]
    if len(r):
        px=float(calibrations[0]) if len(calibrations) else np.nan
        if not np.isfinite(px) or px<=0:raise ValueError('Invalid physical calibration')
        cores,geoms,_=geometry_data(row['geometry']);tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
        allcores=shapely.union_all([g for _,g in cores if g.is_valid]);allnear=allcores.buffer(20/px).intersection(tissue)
        cells=canonical_cells(pd.read_parquet(typ),pd.read_parquet(raw));lookup=subset_mapping(catalog)
        if not {'x_source_px','y_source_px'}.issubset(cells):raise ValueError('Exact source-pixel centroids are required for spatial follow-up')
        cells['immune_subset']=cells.base_type.map(lookup)
        members=memberships.loc[memberships.case_id.eq(cid)&memberships.compartment.eq('peri')&memberships.band_um.eq(20)].copy()
        members.object_id=members.object_id.astype(str)
        joined=members.merge(cells,on=['case_id','object_id'],validate='many_to_one');groups=joined.groupby('islet_id')
        peri_regions=[]
        for name,g in cores:
            iid=cid+'_'+name
            if not g.is_valid or not ((r.islet_id==iid)&r.compartment.eq('peri')).any():continue
            peri=g.buffer(20/px).difference(allcores).intersection(tissue)
            if peri.area<=0:continue
            peri_regions.append(peri)
            # Distances use actual positive pixels inside or immediately around
            # the assessed compartment, not cleaned lymphatic object envelopes.
            lymph=clip_parts(peri.buffer(20/px),geoms.get('Lymphatic'))
            near=peri.intersection(lymph.buffer(20/px));near_fraction=near.area/peri.area
            d=groups.get_group(iid) if iid in groups.groups else joined.iloc[:0]
            points=shapely.points(d.x_source_px.to_numpy(float),d.y_source_px.to_numpy(float))
            is_near=shapely.covers(near,points)
            for subset in ['CD8','CD4','B','macrophage']:
                selected=d.typed&d.immune_subset.eq(subset);known=subset in set(lookup.values())
                neighborhoods.append(dict(islet_id=iid,case_id=cid,subset=subset,
                    near_count=float(np.sum(selected.to_numpy()&is_near)) if known else np.nan,
                    total_count=float(selected.sum()) if known else np.nan,near_area_fraction=near_fraction,
                    near_area_um2=near.area*px**2,valid_area_um2=peri.area*px**2,
                    unresolved_cells=int((~d.typed&d.measured).sum()),validated=bool(known),distance_um=20))
        peri_union=shapely.union_all(peri_regions);distant=tissue.difference(allnear)
        for compartment,region in [('peri',peri_union),('distant',distant)]:
            localization.append(dict(case_id=cid,compartment=compartment,positive_area_um2=clip_parts(region,geoms.get('Lymphatic')).area*px**2,
                                     valid_area_um2=region.area*px**2,validated=region.area>0,pixel_size_um=px,
                                     definition='Union of assessable non-overlapping peri20 regions; distant valid tissue excludes all islet cores and their peri20 bands'))
    n=pd.DataFrame(neighborhoods,columns=['islet_id','case_id','subset','near_count','total_count','near_area_fraction','near_area_um2','valid_area_um2','unresolved_cells','validated','distance_um'])
    l=pd.DataFrame(localization,columns=['case_id','compartment','positive_area_um2','valid_area_um2','validated','pixel_size_um','definition'])
    n.to_parquet(directory/'neighborhoods.parquet',index=False);l.to_parquet(directory/'localization.parquet',index=False)
    write_json(checkpoint,dict(fingerprint=fingerprint,n_neighborhood_rows=len(n),n_localization_rows=len(l),
        status='materialized' if len(r) else 'no_validated_lymphatic_regions'))
    return n,l


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True,type=Path);p.add_argument('--workers',type=int,default=1);a=p.parse_args()
    root=a.input_root;cohort=pd.read_csv(root/'cohort.csv',dtype=str);catalog=json.loads((root/'typing/catalog.json').read_text())
    regions=pd.read_parquet(root/'analysis/regions.parquet');regions.case_id=regions.case_id.astype(str)
    memberships=pd.read_parquet(root/'analysis/memberships.parquet');memberships.case_id=memberships.case_id.astype(str)
    results=Parallel(n_jobs=a.workers)(delayed(donor_endpoints)(root,row,catalog,regions,memberships) for row in cohort.to_dict('records'))
    pd.concat([r[0] for r in results],ignore_index=True).to_parquet(root/'analysis/lymphatic_neighborhoods.parquet',index=False)
    pd.concat([r[1] for r in results],ignore_index=True).to_parquet(root/'analysis/lymphatic_localization.parquet',index=False)

if __name__=='__main__':main()
