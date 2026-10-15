#!/usr/bin/env python3
"""Area-weighted nerve proximity with a shape-preserving compartment null.

The primary null rotates and translates the observed nerve geometry, accepting
only complete placements inside its compartment. An explicitly requested
legacy analytic null remains available for archived descriptive reproduction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import read_cohort, write_json, offline_report_root, SEED, sha256
from build_candidate import geometry_data

def clip_parts(region, entry):
    if entry is None: return shapely.GeometryCollection()
    parts,tree=entry;ids=tree.query(region,predicate='intersects')
    return shapely.union_all(shapely.intersection(parts[ids],region)) if len(ids) else shapely.GeometryCollection()

def shape_preserving_proximity(region,nerve,capillary,pixel_size_um,distance_um=20,replicates=999,seed=SEED,max_proposals=None):
    """Rigid rotations/translations, accepted only wholly inside the compartment.

    Rejection sampling preserves the entire observed nerve shape and component
    relationships. An insufficient set of feasible placements is reported as
    unsupported, never replaced with an analytic or truncated-shape null.
    """
    from shapely.affinity import rotate,translate
    if not np.isfinite(pixel_size_um) or pixel_size_um<=0:raise ValueError('Explicit positive image calibration required')
    if not region.is_valid or not nerve.is_valid or not capillary.is_valid:raise ValueError('Invalid spatial geometry')
    empty=dict(null_method='shape_preserving_compartment_rigid',null_successes=0,nerve_area_um2=float(nerve.area*pixel_size_um**2),valid_area_um2=float(region.area*pixel_size_um**2),nerve_proximity_observed=np.nan,
               nerve_proximity_null=np.nan,excess_proximity=np.nan,spatial_p_value=np.nan)
    if region.area<=0 or nerve.area<=0:return dict(empty,null_status='no_assessable_nerve_area')
    if nerve.difference(region).area>max(1e-9,nerve.area*1e-9):raise ValueError('Nerve shape is not clipped to its compartment')
    near=capillary.buffer(distance_um/pixel_size_um).intersection(region)
    observed=nerve.intersection(near).area/nerve.area
    rng=np.random.default_rng(seed);minx,miny,maxx,maxy=region.bounds;origin=nerve.centroid
    values=[];proposals=0;limit=max_proposals if max_proposals is not None else max(10000,replicates*100)
    while len(values)<replicates and proposals<limit:
        proposals+=1;rotated=rotate(nerve,float(rng.uniform(0,360)),origin=(origin.x,origin.y))
        shape=translate(rotated,xoff=float(rng.uniform(minx,maxx)-origin.x),yoff=float(rng.uniform(miny,maxy)-origin.y))
        if not region.covers(shape):continue
        values.append(shape.intersection(near).area/nerve.area)
    result=dict(empty,nerve_proximity_observed=float(observed),null_successes=len(values),null_proposals=proposals,
                null_status='supported' if len(values)==replicates else 'insufficient_feasible_shape_preserving_placements')
    if len(values)==replicates:
        null=np.asarray(values);result.update(nerve_proximity_null=float(null.mean()),excess_proximity=float(observed-null.mean()),
            spatial_p_value=float((1+np.sum(np.abs(null-null.mean())>=abs(observed-null.mean())))/(1+len(null))),
            null_ci_low=float(np.quantile(null,.025)),null_ci_high=float(np.quantile(null,.975)))
    return result


def positive_geometry(path):
    from pathlib import Path
    import gzip
    path=Path(path);opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as stream:data=json.load(stream)
    if data.get('type')=='FeatureCollection':parts=[shapely.geometry.shape(f['geometry']) for f in data['features']]
    else:parts=[shapely.geometry.shape(data.get('geometry',data))]
    if any(not g.is_valid for g in parts):raise ValueError(f'Invalid positive-pixel polygons: {path}')
    return shapely.union_all(parts)


def run_shape_preserving(cohort,out,mask_manifest,workers=1,replicates=999):
    from joblib import Parallel,delayed
    masks=pd.read_csv(mask_manifest,dtype={'case_id':str})
    required={'case_id','pixel_size_um','nerve_validated','capillary_validated'}
    if not required.issubset(masks) or masks.case_id.duplicated().any():raise ValueError('Incomplete positive-pixel mask manifest')
    if set(masks.case_id)!=set(cohort.case_id):raise ValueError('Structural mask manifest differs from selected cohort')
    assessments=pd.read_parquet(out/'analysis/regions.parquet');assessments.case_id=assessments.case_id.astype(str)
    tasks=[]
    for row in cohort.itertuples(index=False):
        info=masks.loc[masks.case_id.eq(row.case_id)].iloc[0]
        if str(info.nerve_validated).lower()!='true' or str(info.capillary_validated).lower()!='true':continue
        cores,geoms,_=geometry_data(row.geometry);tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
        allcores=shapely.union_all([g for _,g in cores if g.is_valid]);px=float(info.pixel_size_um)
        nerve=positive_geometry(info.nerve_positive_geometry) if 'nerve_positive_geometry' in masks and pd.notna(info.nerve_positive_geometry) else None
        capillary=positive_geometry(info.capillary_positive_geometry) if 'capillary_positive_geometry' in masks and pd.notna(info.capillary_positive_geometry) else None
        if nerve is None and geoms.get('Nerve') is None:continue
        if capillary is None and geoms.get('Capillary') is None:continue
        for name,g in cores:
            iid=f'{row.case_id}_{name}'
            assessed=assessments.loc[assessments.islet_id.eq(iid)&assessments.geometry_valid]
            if not g.is_valid:continue
            for band in [10,20,30]:
                for comp,region in [('core',g.intersection(tissue)),('peri',g.buffer(band/px).difference(allcores).intersection(tissue))]:
                    if not ((assessed.band_um==band)&(assessed.compartment==comp)).any():continue
                    seed=int.from_bytes(hashlib.sha256(f'{SEED}|{iid}|{band}|{comp}'.encode()).digest()[:8],'little')
                    tasks.append((dict(islet_id=iid,case_id=row.case_id,band_um=band,compartment=comp,overlap=bool(assessed.loc[(assessed.band_um==band)&(assessed.compartment==comp),'overlap'].iloc[0])),region,
                                  nerve.intersection(region) if nerve is not None else clip_parts(region,geoms.get('Nerve')),
                                  capillary.intersection(region.buffer(20/px)) if capillary is not None else clip_parts(region.buffer(20/px),geoms.get('Capillary')),px,seed))
    from pathlib import Path
    checkpoint=out/'analysis/spatial_null_checkpoints';checkpoint.mkdir(parents=True,exist_ok=True)
    codehash=sha256(__file__)
    def one(task):
        record,region,nerve,capillary,px,seed=task
        fingerprint=hashlib.sha256(region.wkb+nerve.wkb+capillary.wkb+f'{px}|{replicates}|{seed}|{codehash}'.encode()).hexdigest()
        path=checkpoint/(fingerprint+'.json')
        if path.exists():result=json.loads(path.read_text())
        else:
            result=shape_preserving_proximity(region,nerve,capillary,px,replicates=replicates,seed=seed)
            write_json(path,{k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in result.items()})
        return dict(record,**result)
    columns=['islet_id','case_id','band_um','compartment','null_method','null_successes','excess_proximity','null_status']
    return pd.DataFrame(Parallel(n_jobs=workers)(delayed(one)(t) for t in tasks)) if tasks else pd.DataFrame(columns=columns)


def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--output-root',required=True);p.add_argument('--report-root',required=True);p.add_argument('--cohort',required=True);p.add_argument('--mask-manifest');p.add_argument('--workers',type=int,default=1);p.add_argument('--replicates',type=int,default=999);p.add_argument('--legacy-analytic',action='store_true');a=p.parse_args()
    cohort=read_cohort(a.cohort,a.input_root);out=Path(a.output_root);reports=offline_report_root(out,a.report_root);rows=[]
    if not a.legacy_analytic:
        if not a.mask_manifest:raise ValueError('Positive-pixel structural mask manifest required for primary shape-preserving spatial null')
        result=run_shape_preserving(cohort,out,a.mask_manifest,a.workers,a.replicates)
        result.to_parquet(out/'analysis/coupling.parquet',index=False)
        write_json(reports/'spatial_null.json',dict(method='shape_preserving_compartment_rigid',replicates_required=999,
                   requested=a.replicates,assessed_regions=len(result),supported_regions=int(result.null_status.eq('supported').sum()),
                   limitation='Rigid shape placement can be unidentifiable in narrow or irregular compartments; unsupported regions retain missing excess proximity.'))
        return
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
