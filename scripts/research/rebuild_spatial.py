#!/usr/bin/env python3
"""Rebuild physical compartments and UUID memberships after native core export.

No inherited cell parents, median marker features or old islet summaries are
used. Missing structural masks are missing measurements, never empty masks.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import gzip
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape,mapping
from common import MARKERS,sha256,write_json
from reprocess_candidate import feature_list,class_name
from build_candidate import membership

STRUCTURES=['Nerve','Capillary','Lymphatic','ColIV','Duct']


def overlap_flags(objects,pixel_size_um):
    """Flag positive area intersections of these exact, assessable regions."""
    tree=shapely.STRtree(objects)
    return np.array([any(j!=i and g.intersection(objects[j]).area*pixel_size_um**2>1e-6
                         for j in tree.query(g,predicate='intersects'))
                     for i,g in enumerate(objects)],dtype=bool)


def duct_distance(core,tissue,index,pixel_size_um):
    """Observed distance is censored when an unseen outside-field duct may be nearer."""
    boundary_distance=core.distance(tissue.boundary)*pixel_size_um
    if index is None:
        return dict(duct_distance_um=np.nan,duct_distance_observed_um=np.nan,
                    duct_distance_lower_bound_um=np.nan,duct_distance_censored=True,
                    duct_distance_reason='mask_unavailable')
    objects,tree=index
    observed=core.distance(objects[tree.nearest(core)])*pixel_size_um if len(objects) else np.inf
    censored=bool(observed>boundary_distance)
    return dict(duct_distance_um=observed if not censored else np.nan,
                duct_distance_observed_um=observed if np.isfinite(observed) else np.nan,
                duct_distance_lower_bound_um=min(observed,boundary_distance),
                duct_distance_censored=censored,
                duct_distance_reason='outside_assessable_field_may_be_nearer' if censored else '')


def reviewed_regions(root,cid):
    registry_path=root/'boundaries/assessment_registry.json'
    if not registry_path.exists():raise ValueError('Explicit boundary assessment registry required')
    registry=json.loads(registry_path.read_text())
    path=Path(registry['flagged_regions_path'])
    if sha256(path)!=registry['flagged_regions_sha256']:raise ValueError('Boundary flags changed')
    records=[f for f in feature_list(path) if str(f['properties']['case_id'])==cid]
    flags={f['properties']['candidate_id']:f['properties'] for f in records}
    technical=[shape(f['geometry']) for f in records if f['properties']['decision']=='technical_region_unassessable']
    return flags,registry_path,technical


def local_union(region,index):
    if index is None:return None
    objects,tree=index;ids=tree.query(region,predicate='intersects')
    return shapely.union_all(shapely.intersection(objects[ids],region)) if len(ids) else shapely.GeometryCollection()


def spatial_donor(root,cid):
    root=Path(root);start=time.monotonic()
    row=pd.read_csv(root/'cohort.csv',dtype=str).set_index('case_id').loc[cid]
    folder=root/'analysis/by_donor'/cid;folder.mkdir(parents=True,exist_ok=True)
    native=Path(row.native_core_root) if 'native_core_root' in row else root/'boundaries/native'/cid
    proof=json.loads((native/'manifest.json').read_text())
    if sha256(native/'cores.geojson.gz')!=proof['cores_sha256']:raise ValueError('Native cores changed')
    if sha256(native/'native_core_memberships.csv')!=proof['memberships_sha256']:raise ValueError('Native membership changed')
    source=json.loads(Path(row.boundary_source).read_text());px=float(source['pixel_size_um'])
    if px<=0 or not np.isfinite(px):raise ValueError('Missing physical calibration')
    features=feature_list(native/'cores.geojson.gz')
    previous=[f for f in feature_list(row.previous_geometry) if class_name(f) in STRUCTURES]
    core=np.asarray([shape(f['geometry']) for f in features],dtype=object)
    if not len(core) or not shapely.is_valid(core).all():raise ValueError('Invalid regenerated core')
    names=[f['properties']['name'] for f in features]
    policies=[{'PolygonROI':'qupath_polygon','GeometryROI':'qupath_geometry'}[f['properties']['native_roi_type']] for f in features]
    adaptive=Path(row.boundary_candidates).parent if 'boundary_candidates' in row else root/'boundaries/adaptive'/cid
    adaptive_manifest=json.loads((adaptive/'prediction_manifest.json').read_text())
    tissue_path=Path(adaptive_manifest.get('assessable_tissue_path',Path(row.boundary_source).parent/'valid_tissue.geojson'))
    tissue=shape(json.loads(tissue_path.read_text()))
    flags,review_path,technical=reviewed_regions(root,cid)
    if technical:tissue=tissue.difference(shapely.union_all(technical))
    if set(flags)-set(names):raise ValueError('Boundary review names do not match regenerated native cores')
    masks={}
    for name in STRUCTURES:
        geoms=[shapely.make_valid(shape(f['geometry'])) for f in previous if class_name(f)==name]
        # An explicit new structural manifest can declare an assessed empty
        # mask. Absence of a class in the inherited export cannot do that.
        masks[name]=(np.asarray(geoms,dtype=object),shapely.STRtree(geoms)) if geoms else None
    registry_path=root/'structural'/f'{cid}.json'
    registry=json.loads(registry_path.read_text()) if registry_path.exists() else {}
    for name,record in registry.get('masks',{}).items():
        if name not in STRUCTURES+['DuctEnvelope']:raise ValueError('Unknown structural class')
        path=Path(record['path'])
        if sha256(path)!=record['sha256']:raise ValueError('Structural mask changed')
        objects=[shape(f['geometry']) for f in feature_list(path) if class_name(f)==name]
        if not objects and not record.get('assessed_empty',False):raise ValueError('Undeclared empty mask')
        if not all(g.is_valid for g in objects):raise ValueError('Invalid structural mask')
        masks[name]=(np.asarray(objects,dtype=object),shapely.STRtree(objects))
        previous=[f for f in previous if class_name(f)!=name]+[f for f in feature_list(path) if class_name(f)==name]
    raw=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet')
    if not {'x_source_px','y_source_px'}.issubset(raw):raise ValueError('Exact native centroids required')
    corrected=pd.read_parquet(root/'corrected'/f'{cid}-cells.parquet',columns=['object_id']+MARKERS)
    if corrected.object_id.duplicated().any() or set(corrected.object_id)!=set(raw.object_id):raise ValueError('Corrected/raw UUID mismatch')
    corrected=corrected.set_index('object_id').loc[raw.object_id]
    measured=np.isfinite(corrected.to_numpy(float)).all(axis=1)
    x=raw.x_source_px.to_numpy(float);y=raw.y_source_px.to_numpy(float)
    points=shapely.points(x,y);index=shapely.STRtree(points);core_index=shapely.STRtree(core)
    clipped_cores=np.asarray([g.intersection(tissue) for g in core],dtype=object)
    core_overlap=overlap_flags(clipped_cores,px)
    peris={};peri_overlaps={}
    for band in [10,20,30]:
        these=[]
        for g in core:
            envelope=g.buffer(band/px,quad_segs=16).intersection(tissue)
            nearby=core_index.query(envelope,predicate='intersects')
            these.append(envelope.difference(shapely.union_all(core[nearby])))
        peris[band]=np.asarray(these,dtype=object)
        peri_overlaps[band]=overlap_flags(peris[band],px)
    core_members=[];members=[];islets=[];regions=[];spatial=[]
    for i,(name,g) in enumerate(zip(names,core)):
        iid=f'{cid}_{name}';clipped=clipped_cores[i]
        cm=membership(index,clipped,x,y,policies[i])
        # Native parity checks the retained raw boundary. Analysis membership
        # additionally clips confirmed technical image regions from this core.
        native_cm=membership(index,g,x,y,policies[i])
        core_members.extend((raw.object_id.iloc[j],name) for j in native_cm)
        flag=flags.get(name,{})
        review_status=flag.get('decision','regenerated_systematic_sample_review')
        features[i]['properties'].update(primary_assessable=not flag.get('exclude_primary',False),
                                          annotation_status=review_status)
        base=dict(case_id=cid,islet_id=iid,islet_key=name,pixel_size_um=px,
            donor_status=row.donor_status,donor_source=row.donor_source,
            area_um2=g.area*px**2,centroid_x_um=g.centroid.x*px,centroid_y_um=g.centroid.y*px,
            geometry_valid=True,core_tissue_coverage=clipped.area/g.area,
            overlap=bool(core_overlap[i] or peri_overlaps[20][i]),core_overlap=bool(core_overlap[i]),
            **{f'peri_overlap_{b}um':bool(peri_overlaps[b][i]) for b in [10,20,30]},
            annotation_status=review_status,primary_assessable=not flag.get('exclude_primary',False),
            geometry_sha256=proof['cores_sha256'])
        boundary_measures={}
        for tolerance in [1,2,3]:
            local=local_union(g.boundary.buffer(tolerance/px),masks['ColIV'])
            coverage=g.boundary.intersection(local.buffer(tolerance/px)).length/g.boundary.length if local is not None and g.boundary.length>0 else np.nan
            boundary_measures[f'coliv_boundary_coverage_{tolerance}um']=coverage
        boundary_measures['coliv_boundary_coverage']=boundary_measures['coliv_boundary_coverage_2um']
        boundary_measures.update(duct_distance(g,tissue,masks.get('DuctEnvelope',masks['Duct']),px))
        totals={}
        for band in [10,20,30]:
            peri=peris[band][i]
            for comp,region in [('core',clipped),('peri',peri)]:
                ids=cm if comp=='core' else membership(index,region,x,y,'qupath_geometry')
                if comp=='peri' or band==20:
                    members.extend(dict(case_id=cid,object_id=raw.object_id.iloc[j],islet_id=iid,
                                        compartment=comp,band_um=band) for j in ids)
                record=dict(base,compartment=comp,band_um=band,valid_area_um2=region.area*px**2,
                            raw_cells=len(ids),measured_cells=int(measured[ids].sum()))
                record['overlap']=bool(core_overlap[i] if comp=='core' else peri_overlaps[band][i])
                for structure in STRUCTURES:
                    obj=local_union(region,masks[structure])
                    available=obj is not None and region.area>0
                    area=obj.area*px**2 if available else np.nan
                    record[structure.lower()+'_area_um2']=area
                    record[structure.lower()+'_area_fraction']=area/record['valid_area_um2'] if available else np.nan
                    certified=bool(registry.get('masks',{}).get(structure,{}).get('validated',False))
                    record[structure.lower()+'_validated']=certified
                    if comp=='peri' or band==20:
                        spatial.append(dict(case_id=cid,islet_id=iid,compartment=comp,band_um=band,
                            structure=structure,positive_area_um2=area,valid_area_um2=record['valid_area_um2'],
                            available=available,validated=certified,source='structural_registry' if structure in registry.get('masks',{}) else 'inherited_positive_pixel_classifier_polygons'))
                record.update(boundary_measures)
                regions.append(record)
                if band==20:totals.update({f'{comp}_valid_area_um2':record['valid_area_um2'],
                    f'total_cells_{comp}':len(ids),f'measured_cells_{comp}':int(measured[ids].sum())})
        islets.append(dict(base,**totals,peri_area_um2=totals['peri_valid_area_um2']))
    expected=pd.read_csv(native/'native_core_memberships.csv',dtype=str)
    actual=set(core_members);reference=set(map(tuple,expected[['object_id','islet_key']].to_numpy()))
    if actual!=reference:
        raise ValueError(f'{cid}: Python/native core membership mismatch ({len(actual-reference)} extra, {len(reference-actual)} absent)')
    output_geometry=root/'geometry'/f'{cid}.geojson.gz'
    tissue_feature=dict(type='Feature',geometry=mapping(tissue),properties=dict(name='Assessable tissue',
                        classification=dict(name='Tissue',color=[180,180,180]),source_sha256=sha256(tissue_path)))
    peri_features=[dict(type='Feature',geometry=mapping(peris[20][i]),properties=dict(objectType='annotation',
                        name=name+'_exp20um',classification=dict(name='IsletExpanded',color=[0,180,180]),
                        compartment='peri',band_um=20,source_islet=name,
                        definition='Assessable 20um outer band minus all cores; technical image regions clipped'))
                   for i,name in enumerate(names) if not peris[20][i].is_empty]
    with gzip.open(output_geometry,'wt') as stream:json.dump(dict(type='FeatureCollection',features=features+peri_features+[tissue_feature]+previous),stream)
    (root/'gson').mkdir(exist_ok=True)
    target=root/'gson'/f'{int(cid):04d}.geojson.gz'
    if target.exists():target.unlink()
    target.symlink_to(output_geometry)
    outputs={'islets':islets,'memberships':members,'regions':regions,'spatial_measurements':spatial}
    for name,data in outputs.items():pd.DataFrame(data).to_parquet(folder/f'{name}.parquet',index=False)
    write_json(folder/'complete.json',dict(case_id=cid,status='spatial_rebuilt_native_parity_passed',
        islets=len(islets),core_memberships=len(reference),memberships=len(members),
        exact_native_core_parity=True,pixel_size_um=px,
        geometry_sha256=sha256(output_geometry),raw_sha256=sha256(root/'raw'/f'{cid}-cells.parquet'),
        structural_registry_sha256=sha256(registry_path) if registry_path.exists() else None,
        boundary_assessment_sha256=sha256(review_path),assessable_tissue_sha256=sha256(tissue_path),
        native_core_manifest_sha256=sha256(native/'manifest.json'),runtime_seconds=time.monotonic()-start))
    return cid,len(islets)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--workers',type=int,default=7);p.add_argument('--donor',action='append')
    args=p.parse_args();root=args.input_root.resolve(strict=True)
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);donors=args.donor or cohort.case_id.tolist()
    if set(donors)-set(cohort.case_id):raise ValueError('Unknown donor')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for f in as_completed([pool.submit(spatial_donor,root,cid) for cid in donors]):print(f.result(),flush=True)
    if set(donors)==set(cohort.case_id):
        for name in ['islets','memberships','regions','spatial_measurements']:
            frame=pd.concat([pd.read_parquet(root/'analysis/by_donor'/cid/f'{name}.parquet') for cid in cohort.case_id],ignore_index=True)
            frame.to_parquet(root/'analysis'/f'{name}.parquet',index=False)
        frame=pd.read_parquet(root/'analysis/islets.parquet')
        frame[['case_id','islet_key','islet_id','centroid_x_um','centroid_y_um','area_um2']].rename(columns={'islet_id':'combined_islet_id'}).to_csv(root/'islet_spatial_lookup.csv',index=False)


if __name__=='__main__':main()
