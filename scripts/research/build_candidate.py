#!/usr/bin/env python3
"""Build a separate 21-donor candidate from the explicitly selected inputs.

Original raw objects, images and exports are hardlinked read-only by convention;
derived tables are new files. No phenotype transfer model or expression-scale
mixing is permitted. Large per-donor computations checkpoint independently.
"""
from __future__ import annotations
import argparse
import gzip
import json
import os
import shutil
import time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
from common import MARKERS, RETAINED, ENDOCRINE, read_cohort, sha256, write_json, require_cohort
from phenotypes import read_rules, audit

def materialize(src, dst):
    src,dst=Path(src),Path(dst); dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():
        if not os.path.samefile(src,dst): raise ValueError(f'Conflicting materialized file: {dst}')
        return
    try: os.link(src,dst)
    except OSError: shutil.copy2(src,dst)

def geometry_data(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as stream:
        doc=json.load(stream)
    fs=doc['features'] if isinstance(doc,dict) else doc
    classes={}; cores=[]
    for f in fs:
        p=f['properties']; cls=p.get('classification',{}).get('name','')
        g=shape(f['geometry'])
        if cls == 'Islet': cores.append((p['name'],g))
        elif cls in ['Tissue','Capillary','Nerve','Lymphatic','Exclude']:
            classes.setdefault(cls,[]).append(g)
    if not cores or 'Tissue' not in classes: raise ValueError('Core/tissue geometry required')
    if len(set(n for n,g in cores)) != len(cores): raise ValueError('Duplicate core identifiers')
    # Source invalid polygons remain archived; repair only structural union inputs,
    # recording their number. Invalid cores are excluded, never silently repaired.
    bad={k:sum(not g.is_valid for g in gs) for k,gs in classes.items()}
    result={}
    for k,gs in classes.items():
        parts=np.array([shapely.make_valid(g) for g in gs],dtype=object)
        result[k]=(parts,shapely.STRtree(parts)) if k in ['Nerve','Capillary','Lymphatic'] else shapely.union_all(parts)
    return cores,result,bad

def positive_area(region, structural):
    if structural is None: return np.nan
    if region.is_empty: return 0.
    parts,tree=structural
    nearby=tree.query(region,predicate='intersects')
    # Union the intersected local pieces so an overlapping mask cannot inflate
    # positive area. The spatial index avoids scanning a whole-slide union.
    return shapely.union_all(shapely.intersection(parts[nearby],region)).area if len(nearby) else 0.

def membership(index, geom, x, y, boundary_policy='qupath_polygon'):
    if geom.is_empty: return np.array([],dtype=int)
    if boundary_policy == 'qupath_geometry':
        # GeometryROI.contains includes both JTS INTERIOR and BOUNDARY.
        return np.sort(index.query(geom,predicate='covers'))
    if boundary_policy != 'qupath_polygon':
        raise ValueError('Unknown native ROI containment policy')
    ids=index.query(geom,predicate='contains')
    # Use the audited QuPath PolygonROI winding convention for exact ties.
    ties=index.query(geom.boundary,predicate='intersects')
    if len(ties) and geom.geom_type == 'Polygon' and not geom.interiors:
        import sys
        sys.path.insert(0,str(Path(__file__).parents[1]/'cohort_expansion'))
        from audit_preparation import polygon_boundary_contains
        ids=np.r_[ids,[i for i in ties if polygon_boundary_contains(list(geom.exterior.coords),x[i],y[i])]].astype(int)
    return np.sort(ids)

def process(row, out, rules, gate_table, gate_provenance, core_roi_types=None):
    cid=row.case_id; start=time.monotonic()
    source=ad.read_h5ad(row.cells,backed='r')
    try:
        if list(source.var_names)!=MARKERS: raise ValueError('Ordered marker panel mismatch')
        sel=np.flatnonzero(source.obs.imageid.astype(str).to_numpy()==cid)
        if not len(sel): raise ValueError(f'No raw cells for {cid}')
        obs=source.obs.iloc[sel].copy()
        values=np.asarray((source.raw.X if row.cell_slot=='raw' else source.X)[sel],dtype=np.float32)
    finally: source.file.close()
    measured=np.isfinite(values).all(axis=1)
    if np.any(np.isfinite(values).any(axis=1) & ~measured) or np.any(values[measured]<0):
        raise ValueError('Partial missingness or negative raw means')
    logs=np.log1p(values)
    # The original log layer and the additions are both raw log1p MFI. No scVI
    # X, disease covariates, hormone ratios, or per-donor quantile normalization.
    labels=obs.get('phenotype',pd.Series('',index=obs.index)).astype(str).to_numpy()
    gates=gate_table.loc[gate_table.imageid.astype(str)==cid].set_index('markers').gate.to_dict()
    gp=gate_provenance.loc[gate_provenance.imageid.astype(str)==cid]
    validated=gp.loc[gp.validated.eq(True),'markers'].tolist()
    # Provisional PDPN evidence stays false in the registry until reviewed.
    reports,masks=audit(labels,logs,MARKERS,rules,gates,validated)
    for r in reports: r['case_id']=cid
    cores,geoms,bad=geometry_data(row.geometry)
    tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
    px=.3774; area_scale=px**2
    x=obs.X_centroid.to_numpy(float)/px; y=obs.Y_centroid.to_numpy(float)/px
    points=shapely.points(x,y); idx=shapely.STRtree(points)
    valid_cores=[g for _,g in cores if g.is_valid and not g.is_empty]
    union=shapely.union_all(valid_cores)
    # Flag overlaps at 30 um so all primary and sensitivity footprints are
    # non-overlapping within the same assessable set.
    footprints=[g.buffer(30/px).intersection(tissue) if g.is_valid else shapely.GeometryCollection() for _,g in cores]
    tree=shapely.STRtree(footprints)
    overlap=[any(j!=i and g.intersection(footprints[j]).area>1e-8 for j in tree.query(g,predicate='intersects'))
             for i,g in enumerate(footprints)]
    rows=[]; counts=[]; memberships=[]; distance=[]
    for i,(name,g) in enumerate(cores):
        valid=g.is_valid and not g.is_empty and g.area>0
        core=g.intersection(tissue) if valid else shapely.GeometryCollection()
        base=dict(islet_id=f'{cid}_{name}',case_id=cid,islet_key=name,donor_status=row.donor_status,
                  donor_source=row.donor_source,area_um2=g.area*area_scale,
                  centroid_x_um=g.centroid.x*px,centroid_y_um=g.centroid.y*px,
                  geometry_valid=valid,overlap=overlap[i],annotation_review=row.annotation_review,
                  core_tissue_coverage=core.area/g.area if valid else 0)
        for band in [10,20,30]:
            peri=g.buffer(band/px).difference(union).intersection(tissue) if valid else shapely.GeometryCollection()
            for comp,region in [('core',core),('peri',peri)]:
                roi_type=core_roi_types[name] if comp=='core' and core_roi_types is not None else 'PolygonROI'
                policy={'PolygonROI':'qupath_polygon','GeometryROI':'qupath_geometry'}[roi_type]
                members=membership(idx,region,x,y,policy); use=members[measured[members]]
                record={**base,'compartment':comp,'band_um':band,'valid_area_um2':region.area*area_scale,
                        'raw_cells':len(members),'measured_cells':len(use),
                        'unmeasured_cells':len(members)-len(use)}
                for marker,v in zip(MARKERS,np.median(logs[use],axis=0) if len(use) else np.full(31,np.nan)):
                    record[f'median_log1p_{marker}']=v
                for cls in ['Nerve','Capillary','Lymphatic']:
                    pos=positive_area(region,geoms.get(cls))*area_scale
                    record[f'{cls.lower()}_area_um2']=pos
                    record[f'{cls.lower()}_area_fraction']=pos/record['valid_area_um2'] if region.area>0 else np.nan
                if band==20:
                    # Preserve archival counts separately from scientifically validated counts.
                    for label in RETAINED:
                        report=next(r for r in reports if r['phenotype']==label)
                        ok=valid and report['validated']
                        counts.append(dict(islet_id=base['islet_id'],case_id=cid,compartment=comp,
                                           phenotype=label,count=int(masks[label][use].sum()) if ok else np.nan,
                                           archival_count=int((labels[use]==label).sum()) if row.cell_slot=='raw' else np.nan,
                                           validated=ok, measured_cells=len(use), valid_area_um2=region.area*area_scale,
                                           reason='' if ok else ('invalid_geometry' if not valid else report['reason']+': '+report['missing_validated_gates'])))
                    union_mask=np.logical_or.reduce([masks[t] for t in ENDOCRINE])
                    endo_ok=all(next(r['validated'] for r in reports if r['phenotype']==t) for t in ENDOCRINE)
                    record['endocrine_count']=int(union_mask[use].sum()) if valid and endo_ok else np.nan
                    if comp=='core':
                        memberships.extend((base['islet_id'],str(obs.index[j])) for j in members)
                rows.append(record)
        if (i+1)%200==0: print(f'{cid}: {i+1}/{len(cores)} islets',flush=True)
    frame=pd.DataFrame(rows)
    frame.to_parquet(out/'analysis'/f'{cid}.parquet',index=False)
    pd.DataFrame(counts).to_parquet(out/'analysis'/f'{cid}-counts.parquet',index=False)
    pd.DataFrame(reports).to_csv(out/'analysis'/f'{cid}-phenotype-audit.csv',index=False)
    pd.DataFrame(memberships,columns=['islet_id','object_id']).to_parquet(out/'analysis'/f'{cid}-memberships.parquet',index=False)
    cells=pd.DataFrame(values,columns=MARKERS)
    cells.insert(0,'object_id',obs.index.astype(str));cells['case_id']=cid
    cells['X_centroid']=obs.X_centroid.to_numpy();cells['Y_centroid']=obs.Y_centroid.to_numpy()
    cells['archival_phenotype']=labels;cells['marker_measurements_complete']=measured
    for col in ['Classification','Parent','Cell Area','Nucleus Area']:
        if col in obs: cells[col]=obs[col].astype(str).to_numpy() if col in ['Classification','Parent'] else obs[col].to_numpy()
    cells.to_parquet(out/'raw'/f'{cid}-cells.parquet',index=False)
    excluded=cells.loc[~measured,['object_id','case_id','X_centroid','Y_centroid']].copy()
    excluded['reason']='all_31_marker_means_missing';excluded.to_csv(out/'analysis'/f'{cid}-unmeasured.csv',index=False)
    write_json(out/'analysis'/f'{cid}-complete.json',dict(case_id=cid,raw_cells=len(obs),measured_cells=int(measured.sum()),
               unmeasured_cells=int((~measured).sum()),islets=len(cores),invalid_structural_polygons=bad,
               seconds=time.monotonic()-start))
    print(f'{cid}: complete ({len(obs)} cells, {len(cores)} islets)',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--cohort',required=True);p.add_argument('--input-root',required=True)
    p.add_argument('--output-root',required=True);p.add_argument('--gate-registry',type=Path);p.add_argument('--resume',action='store_true');args=p.parse_args()
    out=Path(args.output_root).resolve();root=Path(args.input_root).resolve()
    if out==root or out==root/'data': raise ValueError('Output must be a separate candidate')
    cohort=read_cohort(args.cohort,root)
    if out.exists() and not args.resume: raise ValueError('Output exists; use --resume for the same manifest only')
    for sub in ['analysis','raw','images','gson','regional','trajectory','summary','logs']:
        (out/sub).mkdir(parents=True,exist_ok=True)
    identity={'cohort_sha256':sha256(args.cohort),'builder_sha256':sha256(__file__),
              'rules_sha256':sha256(root/'data/phenotype_rules.csv'),
              'phenotype_engine_sha256':sha256(Path(__file__).with_name('phenotypes.py')),
              'gate_registry_sha256':sha256(args.gate_registry) if args.gate_registry else None}
    identity_path=out/'build_identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text())!=identity: raise ValueError('Resume inputs/code differ')
    write_json(identity_path,identity)
    shutil.copyfile(args.cohort,out/'cohort.csv')
    cohort.to_csv(out/'source_manifest.csv',index=False)
    rules=read_rules(root/'data/phenotype_rules.csv')
    # No implicit original-15 gate table. Measurement-only builds explicitly
    # leave assignments missing until an all-donor evidence registry is supplied.
    gates=pd.DataFrame(columns=['imageid','markers','gate','validated'])
    if args.gate_registry:
        registry=pd.read_csv(args.gate_registry,dtype={'case_id':str})
        require_cohort(registry,cohort.case_id)
        if registry.duplicated(['case_id','marker']).any():raise ValueError('Duplicate donor/marker gates')
        gates=registry.rename(columns={'case_id':'imageid','marker':'markers','gate_log1p':'gate'})
        if not gates.validated.isin([True,False]).all():raise ValueError('Gate validation must be explicit Boolean values')
    for row in cohort.itertuples(index=False):
        image=Path(row.image)
        for src in [image,image.with_suffix('.xml'),Path(str(image)+'.offsets.json')]:
            if src.exists(): materialize(src,out/'images'/src.name)
        materialize(row.cells,out/'raw'/Path(row.cells).name)
        materialize(row.geometry,out/'gson'/f'{int(row.case_id):04d}.geojson.gz')
        for src in Path(row.regional_root).glob(f'*_{int(row.case_id):04d}.tsv'):
            materialize(src,out/'regional'/src.name)
        if not (out/'analysis'/f'{row.case_id}-complete.json').exists(): process(row,out,rules,gates,gates)
    frames=[pd.read_parquet(out/'analysis'/f'{cid}.parquet') for cid in cohort.case_id]
    regional=pd.concat(frames,ignore_index=True);require_cohort(regional,cohort.case_id)
    regional.to_parquet(out/'analysis/regions.parquet',index=False)
    core=regional.query("band_um == 20 and compartment == 'core'").copy()
    peri=regional.query("band_um == 20 and compartment == 'peri'").copy()
    feature_cols=[f'median_log1p_{m}' for m in MARKERS]
    wide=core.drop(columns=feature_cols).merge(core[['islet_id']+feature_cols].rename(columns={c:'core__'+c[13:] for c in feature_cols}),on='islet_id',validate='one_to_one')
    wide=wide.merge(peri[['islet_id']+feature_cols].rename(columns={c:'peri__'+c[13:] for c in feature_cols}),on='islet_id',validate='one_to_one')
    wide['eligible']=wide.geometry_valid & (wide.measured_cells>=10) & (wide.core_tissue_coverage>.99)
    wide['exclusion_reason']=np.where(~wide.geometry_valid,'invalid_geometry',np.where(wide.measured_cells<10,'fewer_than_10_measured_core_cells',np.where(wide.core_tissue_coverage<=.99,'core_outside_tissue','')))
    empty=~np.isfinite(wide[[c for c in wide if c.startswith(('core__','peri__'))]]).all(axis=1)
    wide.loc[empty,'eligible']=False;wide.loc[empty,'exclusion_reason']='empty_or_unmeasured_compartment'
    wide.to_csv(out/'trajectory/features.csv',index=False)
    core[['case_id','islet_key','islet_id','centroid_x_um','centroid_y_um','area_um2']].rename(columns={'islet_id':'combined_islet_id'}).to_csv(out/'islet_spatial_lookup.csv',index=False)
    allcounts=pd.concat([pd.read_parquet(out/'analysis'/f'{cid}-counts.parquet') for cid in cohort.case_id],ignore_index=True)
    allcounts.to_csv(out/'trajectory/counts.csv',index=False)
    audits=pd.concat([pd.read_csv(out/'analysis'/f'{cid}-phenotype-audit.csv',dtype={'case_id':str}) for cid in cohort.case_id],ignore_index=True)
    audits.to_csv(out/'analysis/phenotype_audit.csv',index=False)
    audits.to_csv(out/'trajectory/phenotype_availability.csv',index=False)
    endocrine=core[['islet_id','case_id','endocrine_count','measured_cells','geometry_valid']].rename(columns={'measured_cells':'measured_core_cells'}).copy()
    endocrine['validated']=endocrine.geometry_valid & endocrine.endocrine_count.notna()
    missing={cid:','.join(sorted(set(';'.join(group.loc[group.phenotype.isin(ENDOCRINE),'missing_validated_gates'].fillna('')).split(';'))-{''}))
             for cid,group in audits.groupby('case_id')}
    endocrine['reason']=np.where(endocrine.validated,'',np.where(~endocrine.geometry_valid,'invalid_geometry',
                                 'gate_review_required: '+endocrine.case_id.map(missing)))
    endocrine['definition']='unique core alpha/beta/delta cell union'
    endocrine['rule_sha256']=sha256(root/'data/phenotype_rules.csv')
    endocrine.to_csv(out/'trajectory/endocrine_counts.csv',index=False)
    write_json(out/'trajectory/marker_panel.json',dict(version=1,markers=MARKERS,measurement='per-islet median log1p core intensity'))
    write_json(out/'trajectory/count_manifest.json',dict(version=1,cohort_sha256=sha256(out/'cohort.csv'),
        rules_sha256=sha256(root/'data/phenotype_rules.csv'),gate_registry_sha256=identity['gate_registry_sha256'],
        counts_sha256=sha256(out/'trajectory/counts.csv'),endocrine_counts_sha256=sha256(out/'trajectory/endocrine_counts.csv'),
        archival_labels_required=False,retained_types=RETAINED,endocrine_union=ENDOCRINE,
        validated_donor_type_pairs=int(audits.validated.sum()),supported_endocrine_donors=sorted(endocrine.loc[endocrine.validated,'case_id'].unique().tolist()),
        status='available' if endocrine.validated.all() else 'gate_review_required'))
    metadata=[]
    for row in cohort.itertuples(index=False):
        d=pd.read_csv(row.metadata,dtype={'case_id':str});d=d.loc[d.case_id==row.case_id]
        if len(d)!=1: raise ValueError('Missing/ambiguous donor metadata')
        d=d.copy();d['donor_source']=row.donor_source;d['donor_status']=row.donor_status;d['inclusion_status']='selected'
        if row.cell_slot=='raw':
            a=ad.read_h5ad(row.cells,backed='r');m=a.obs.loc[a.obs.imageid.astype(str)==row.case_id].iloc[0]
            d['age']=float(m.Age);d['sex']=m.Gender;a.file.close()
        else: d['age']=d.age_years_key
        metadata.append(d)
    pd.concat(metadata,ignore_index=True).to_csv(out/'donor_metadata.csv',index=False)
    shutil.copyfile(root/'data/phenotype_rules.csv',out/'phenotype_rules.csv')
    totals=[json.loads((out/'analysis'/f'{cid}-complete.json').read_text()) for cid in cohort.case_id]
    write_json(out/'build_report.json',dict(donors=21,groups=cohort.donor_status.value_counts().to_dict(),
               raw_cells=sum(r['raw_cells'] for r in totals),unmeasured_cells=sum(r['unmeasured_cells'] for r in totals),
               islets=len(core),eligible_islets=int(wide.eligible.sum()),quantification_pixel_um=.3774,
               harmonization='same ordered panel and raw Cell mean units; median(log1p(raw MFI)); no biological batch correction',
               phenotype_status='retained parent-chain gate audit; uncertified labels remain unavailable',donor_reports=totals))

if __name__=='__main__': main()
