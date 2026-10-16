#!/usr/bin/env python3
"""Assemble the restored app plus the five REVIEWED protein-processed donors.

All outputs are new. Approval records bind the exact rules and geometry hashes;
no historical candidate, scVI model, original islet regeneration or default-path
builder is invoked. The frozen single-cell source remains a read-only partition.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import json
from pathlib import Path
import shutil
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import shapely
from scipy import sparse
from shapely.geometry import mapping, shape

from build_protein_embedding import fit_embedding, sha256
from prepare_protein_expansion import read_column
from protein_phenotyping import RULES, classify
from stage_original_gate_calls import DONORS
from excluded_donors import assert_no_excluded
from compute_neighborhood_metrics import IMMUNE_TYPES, POS_MARKERS_CODEX
sys.path.insert(0, str(Path(__file__).resolve().parent / 'research'))
from build_candidate import geometry_data, membership, positive_area
from audit_original_parity import audit

BASELINE_FILES = ('adata_ins_root.h5ad','islet_explorer.h5ad','barlow_islet_features.csv',
    'islet_spatial_lookup.csv','donor_metadata.csv','neighborhood_metrics.csv',
    'cell_distance.csv','signed_distance_kde.csv','annotations.tsv','phenotype_rules.csv',
    'phenotype_gates_per_image.csv','phenotype_gates_provenance.csv')
BASELINE_DIRS = ('cells','donors','gson','json')
ENDOCRINE = ('Alpha cell','Beta cell','Delta cell','Endocrine')


def validate_review(prepared_root, geometry_root, acceptance):
    prepared=json.loads((prepared_root/'preparation.json').read_text())
    geometry=json.loads((geometry_root/'manifest.json').read_text())
    review=json.loads(Path(acceptance).read_text())
    for name,ids in [('Owner acceptance',review.get('boundaries',{})),('Prepared inputs',[r['case_id'] for r in prepared['donors']]),
                     ('Geometry snapshot',[r['case_id'] for r in geometry['donors']])]:
        assert_no_excluded(ids,name)
    rules_hash=sha256(RULES)
    if review.get('phenotype_rules_sha256') != rules_hash or review.get('phenotypes_accepted') is not True:
        raise ValueError('The exact proposed phenotype rules require owner acceptance')
    if not review.get('reviewer') or not review.get('reviewed_at'):
        raise ValueError('Record who reviewed these outputs and when')
    if set(review.get('boundaries',{})) != set(DONORS):
        raise ValueError('An explicit boundary decision is required for every addition')
    if set(r['case_id'] for r in prepared['donors']) != set(DONORS):
        raise ValueError('Prepared cohort is incomplete')
    if set(r['case_id'] for r in geometry['donors']) != set(DONORS):
        raise ValueError('Geometry cohort is incomplete')
    if sha256(prepared['original_source']) != prepared['original_source_sha256']:
        raise ValueError('Frozen single-cell source changed')
    for row in prepared['donors']:
        if sha256(prepared_root/f"{row['case_id']}.h5ad") != row['prepared_sha256']:
            raise ValueError('Prepared donor changed')
    for row in geometry['donors']:
        cid=row['case_id']; decision=review['boundaries'][cid]
        if decision.get('accepted') is not True or decision.get('geometry_sha256') != row['geometry_sha256']:
            raise ValueError(f'{cid}: exact boundary snapshot needs acceptance')
        for key,digest_key in [('geometry','geometry_sha256'),('centroids','centroids_sha256'),('native_data','native_sha256')]:
            if sha256(row[key]) != row[digest_key]:raise ValueError(f'{cid}: reviewed {key} changed')
    return prepared,geometry,review


def copy_baseline(baseline,out,skip=()):
    out.mkdir(parents=True,exist_ok=False)
    for name in BASELINE_FILES:
        if name not in skip and (baseline/name).is_file():shutil.copy2(baseline/name,out/name)
    for name in BASELINE_DIRS:
        shutil.copytree(baseline/name,out/name)


def reference_frequencies(frozen):
    with h5py.File(frozen) as f:
        labels=f['obs/phenotype/categories'].asstr()[:]
        counts=np.bincount(f['obs/phenotype/codes'][:],minlength=len(labels))
    return dict(zip(labels,counts/counts.sum()))


def get_metadata(cid,cohort,clinical,cells):
    row=cohort.loc[cohort.case_id.eq(cid)].iloc[0]
    extra=clinical.loc[clinical.case_id.eq(cid)]
    extra=extra.iloc[0].to_dict() if len(extra) else {}
    meta={'Case ID':cid,'Donor Status':row.donor_status,
        'Age':float(cells.obs.Age.iloc[0]),'Gender':str(cells.obs.Gender.iloc[0])}
    for k in ('AAb_GADA','AAb_IA2A','AAb_ZnT8A','AAb_mIAA'):
        value=extra.get(k,np.nan)
        meta[k]=np.nan if pd.isna(value) else float(str(value).lower()=='true') if isinstance(value,str) else float(value)
    return meta


def donor_regions(snapshot,raw,centroids):
    """Use exact native centroid coordinates and ROI containment semantics."""
    px=float(snapshot['pixel_size_um'])
    if not np.isfinite(px) or px<=0:raise ValueError('Invalid calibration')
    ids=raw.obs['Object ID'].astype(str)
    if ids.duplicated().any() or centroids.index.duplicated().any() or set(ids)!=set(centroids.index):
        raise ValueError('Native and exported cell UUID census differs')
    native=centroids.loc[ids]
    x,y=native.x_source_px.to_numpy(),native.y_source_px.to_numpy()
    if not np.allclose(x*px,raw.obs.X_centroid,atol=1e-4,rtol=0) or not np.allclose(y*px,raw.obs.Y_centroid,atol=1e-4,rtol=0):
        raise ValueError('Native centroid calibration/alignment differs from raw source')
    cores,geoms,bad=geometry_data(snapshot['geometry'])
    if any(bad.values()):raise ValueError('Invalid structural/tissue geometry requires explicit repair review')
    with gzip.open(snapshot['geometry'],'rt') as stream:features=json.load(stream)['features']
    core_types={f['properties']['name']:f['properties'].get('native_roi_type','PolygonROI') for f in features if f['properties'].get('classification',{}).get('name')=='Islet'}
    tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
    shapely.prepare(tissue)
    def clip_to_tissue(g):
        return g if tissue.covers(g) else g.intersection(tissue)
    if any(not g.is_valid or g.is_empty or g.area<=0 for _,g in cores):raise ValueError('Invalid reviewed core')
    clipped=[clip_to_tissue(g) for _,g in cores]
    tree=shapely.STRtree(clipped)
    overlaps=np.zeros(len(clipped))
    for i,g in enumerate(clipped):
        for j in tree.query(g,predicate='intersects'):
            if j<=i:continue
            area=g.intersection(clipped[j]).area*px**2
            # Native PolygonROI float vertices can leave sub-pixel slivers.
            # Preserve reviewed vertices; reject >0.0001 um2 or any shared cell.
            if area>1e-4:raise ValueError('Overlapping reviewed cores require resolution')
            overlaps[i]+=area;overlaps[j]+=area
    points=shapely.points(x,y);index=shapely.STRtree(points)
    result=[];core_owner=np.full(len(x),-1,dtype=int)
    for i,((name,original),core) in enumerate(zip(cores,clipped)):
        if core.is_empty:raise ValueError('Reviewed core outside valid tissue')
        footprint=original.buffer(20/px)
        # Only intersecting cores can affect this band. Avoid overlaying the
        # whole slide's geometry hundreds of times, preserving the same regions.
        local_cores=shapely.union_all([clipped[j] for j in tree.query(footprint,predicate='intersects')])
        peri=clip_to_tissue(footprint.difference(local_cores))
        roi_type=core_types[name]
        if roi_type not in ('PolygonROI','GeometryROI'):raise ValueError('Unsupported native ROI')
        policy='qupath_polygon' if roi_type=='PolygonROI' else 'qupath_geometry'
        core_ids=membership(index,core,x,y,policy)
        peri_ids=np.sort(index.query(peri,predicate='covers'))
        peri_ids=np.setdiff1d(peri_ids,core_ids)
        if np.any(core_owner[core_ids]>=0):raise ValueError('Cell belongs to multiple cores')
        core_owner[core_ids]=i
        result.append(dict(name=name,core=core,peri=peri,core_ids=core_ids,peri_ids=peri_ids,
                           native_roundoff_overlap_um2=overlaps[i]))
    return result,geoms,points,core_owner,px


def neighborhood(labels,core_ids,peri_ids,coords,frequencies):
    c,p=labels[core_ids],labels[peri_ids];cn,pn=len(c),len(p)
    out=dict(total_cells_core=cn,total_cells_peri=pn)
    core_immune=int(np.isin(c,IMMUNE_TYPES).sum());peri_immune=int(np.isin(p,IMMUNE_TYPES).sum())
    out.update(immune_count_core=core_immune,immune_count_peri=peri_immune,
        immune_frac_core=core_immune/cn,immune_frac_peri=peri_immune/pn,
        immune_ratio=(peri_immune/pn)/(core_immune/cn) if core_immune else np.nan,
        cd8_to_macro_ratio=np.sum(p=='CD8a Tcell')/np.sum(p=='Macrophage') if np.any(p=='Macrophage') else np.nan,
        tcell_density_peri=100*np.isin(p,['CD8a Tcell','CD4 Tcell','T cell']).sum()/pn)
    distance=np.linalg.norm(coords[peri_ids]-coords[core_ids].mean(axis=0),axis=1)
    immune_dist=distance[np.isin(p,IMMUNE_TYPES)]
    out['min_dist_immune_mean']=float(immune_dist.min()) if len(immune_dist) else np.nan
    for label,freq in frequencies.items():
        key=label.replace(' ','_').replace('+','plus');count=int(np.sum(p==label))
        out['peri_count_'+key]=count;out['peri_prop_'+key]=count/pn
        if label!='Unknown':
            expected=freq*pn
            out['enrich_z_'+key]=(count-expected)/np.sqrt(expected) if expected>0 else np.nan
            out['min_dist_'+key]=float(distance[p==label].min()) if count else np.nan
    return out


def process_donor(cid,cells,raw,snapshot,meta,frequencies,out,morphology_names=()):
    centroids=pd.read_csv(snapshot['centroids'],sep='\t',dtype={'object_id':str}).set_index('object_id')
    regions,geoms,points,core_owner,px=donor_regions(snapshot,raw,centroids)
    positions=raw.obs_names.get_indexer(cells.obs_names)
    if np.any(positions<0):raise ValueError('Prepared cell identity missing from raw')
    raw_to_qc=np.full(raw.n_obs,-1,dtype=int);raw_to_qc[positions]=np.arange(cells.n_obs)
    labels,reasons=classify(cells.X,list(cells.var_names))
    cells.obs['phenotype']=pd.Categorical(labels)
    cells.obs['phenotype_rule_reason']=reasons
    cells.obs['archival_parent']=cells.obs['Parent'].astype(str)
    cells.obs['Parent']='Tissue'
    raw_x=np.asarray(raw.X);raw_markers=list(raw.var_names)
    hormones=['INS','GCG','SST'];hi=[raw_markers.index(m) for m in hormones]
    measured=np.isfinite(raw_x).all(axis=1)
    thresholds=np.mean(raw_x[measured][:,hi],axis=0)+1.5*np.std(raw_x[measured][:,hi],axis=0,ddof=1)
    coords=cells.obs[['X_centroid','Y_centroid']].to_numpy(float)
    names=list(frequencies);counts=[];means=[];morph=[];observations=[];features=[];lookup=[]
    composition=[];targets=[];markers=[];unpaired=[];memberships=[]
    native_qc_points=points[positions]
    boundaries=np.array([r['core'].boundary for r in regions],dtype=object)
    nearest=shapely.STRtree(boundaries).nearest(native_qc_points)
    signed=shapely.distance(native_qc_points,boundaries[nearest])*px
    signed[core_owner[positions]>=0]*=-1
    # Tissue display has one nearest-islet context; per-islet overlays keep each
    # exact region membership, including a cell in two overlapping peri bands.
    display_name=np.full(cells.n_obs,'',dtype=object);display_region=np.full(cells.n_obs,'tissue',dtype=object)
    geometry_features=[]
    for i,r in enumerate(regions):
        name=r['name'];identity=f'{cid}_{name}'
        raw_c,raw_p=r['core_ids'],r['peri_ids']
        c=raw_to_qc[raw_c];c=c[c>=0];p=raw_to_qc[raw_p];p=p[p>=0]
        for ids,region in [(c,'core'),(p,'peri')]:
            choose=(core_owner[positions[ids]]==i) if region=='core' else ((core_owner[positions[ids]]<0)&(nearest[ids]==i))
            display_name[ids[choose]]=name;display_region[ids[choose]]=region
        for comp,g in [('Islet',r['core']),('IsletExpanded',r['core'].union(r['peri']))]:
            geometry_features.append(dict(type='Feature',geometry=mapping(g),properties=dict(name=name if comp=='Islet' else name+'_exp20um',classification=dict(name=comp))))
        positive=raw_x[raw_c][:,hi]>=thresholds
        valid=np.isfinite(raw_x[raw_c][:,hi]).all(axis=1)
        npos=positive.sum(axis=1)
        record=dict(meta,islet_id=int(name.removeprefix('Islet_')),name=name,
            cells_total=len(raw_c),Multi_Pos=int(np.sum(valid&(npos>1))),Triple_Neg=int(np.sum(valid&(npos==0))))
        for j,prefix in enumerate(['Ins','Glu','Stt']):
            record[prefix+'_single']=int(np.sum(valid&positive[:,j]&(npos==1)))
            record[prefix+'_any']=int(np.sum(valid&positive[:,j]));record[prefix+'_thr']=thresholds[j]
        # Unmeasured cells remain in the raw denominator, with explicit counts.
        record['unmeasured_cells']=int((~valid).sum());composition.append(record)
        for compartment,g,ids in [('islet_core',r['core'],raw_c),('islet_band',r['peri'],raw_p),('islet_union',r['core'].union(r['peri']),np.union1d(raw_c,raw_p))]:
            region=name+{'islet_core':'','islet_band':'_exp20um','islet_union':'_union'}[compartment]
            for structure in ['Nerve','Capillary','Lymphatic']:
                area=positive_area(g,geoms.get(structure))*px**2
                source_structure=geoms.get(structure)
                count=np.nan
                if source_structure is not None:
                    parts,tree=source_structure
                    nearby=tree.query(g,predicate='intersects')
                    count=int(np.sum(shapely.area(shapely.intersection(parts[nearby],g))>0))
                targets.append(dict(meta,region=region,type=compartment,**{'class':structure},area_um2=area,
                    region_um2=g.area*px**2,area_density=area/(g.area*px**2) if g.area else np.nan,count=count,
                    _islet_name=name,_region_type=compartment))
            for marker in ['CD3e','CD163','CD68','LGALS3']:
                values=raw_x[ids,raw_markers.index(marker)];values=values[np.isfinite(values)]
                mean=float(values.mean()) if len(values) else np.nan
                sd=float(values.std(ddof=1)) if len(values)>1 else np.nan
                threshold=mean+1.5*sd;pos=int(np.sum(values>=threshold)) if np.isfinite(threshold) else np.nan
                markers.append(dict(meta,region=region,region_type=compartment,marker=marker,n_cells=len(values),
                    mean=mean,SD=sd,k=1.5,threshold=threshold,pos_count=pos,pos_frac=pos/len(values) if len(values) else np.nan,
                    _islet_name=name,_region_type=compartment))
        lookup.append(dict(case_id=int(cid),imageid=int(cid),islet_key=name,combined_islet_id=identity,
            centroid_x_um=r['core'].centroid.x*px,centroid_y_um=r['core'].centroid.y*px,area_um2=r['core'].area*px**2))
        if not len(c) or not len(p):
            unpaired.append(dict(case_id=cid,islet_id=identity,core_cells=len(c),peri_cells=len(p),reason='empty_qc_core_or_peri'))
            continue
        cc=np.array([np.sum(labels[c]==label) for label in names]);pc=np.array([np.sum(labels[p]==label) for label in names])
        obs=dict(imageid=int(cid),islet_id=identity,base_islet_id=name,combined_islet_id=identity,match_key=f'{name}_{cid}',
            case_id_str=cid,age=meta['Age'],gender=meta['Gender'],donor_status=meta['Donor Status'],region='islet',total_cells=len(c),
            GADA=meta['AAb_GADA'],IA2A=meta['AAb_IA2A'],ZnT8A=meta['AAb_ZnT8A'],mIAA=meta['AAb_mIAA'])
        obs.update(neighborhood(labels,c,p,coords,frequencies))
        for j,label in enumerate(names):obs['count_'+label]=cc[j];obs['prop_'+label]=cc[j]/len(c)
        observations.append(obs);means.append(np.mean(cells.X[c],axis=0).astype(np.float32));counts.append(cc)
        # Morphology uses the same ordered fields as the frozen aggregate.
        source_rows=positions[c]
        morph.append(raw.obs.iloc[source_rows].reindex(columns=morphology_names).apply(
            pd.to_numeric,errors='coerce').mean().to_numpy(dtype=np.float32))
        endo=int(sum(cc[names.index(t)] for t in ENDOCRINE));denom=max(endo,1)
        feature=dict(islet_id=identity,case_id=int(cid),donor_status=meta['Donor Status'],total_cells_core=len(c),total_cells_peri=len(p),n_endocrine_core=endo)
        for label,cv,pv in zip(names,cc,pc):
            feature.update({f'core_{label}':np.log1p(cv/denom),f'peri_{label}':np.log1p(pv/denom),f'n_core_{label}':cv,f'n_peri_{label}':pv})
        features.append(feature)
        ids=np.r_[c,p];detail=cells.obs.iloc[ids][['X_centroid','Y_centroid','Cell Area','Nucleus Area']].copy()
        detail['phenotype']=labels[ids];detail['cell_region']=np.r_[np.repeat('core',len(c)),np.repeat('peri',len(p))]
        detail=pd.concat([detail.reset_index(drop=True),pd.DataFrame(cells.X[ids],columns=cells.var_names)],axis=1)
        detail.to_csv(out/'cells'/f'{identity}.csv',index=False)
        memberships.extend(dict(cell_id=str(cells.obs_names[j]),islet_id=identity,region=region) for region,ids in [('core',c),('peri',p)] for j in ids)
    if not observations:raise ValueError(f'{cid}: no paired QC-supported islets')
    for region,suffix in [('core',''),('peri','_exp20um')]:
        use=display_region==region;cells.obs.loc[use,'Parent']=display_name[use]+suffix
    cells.obs['Dist to Closest Islet']=signed
    tissue=cells.obs[['X_centroid','Y_centroid']].copy();tissue['phenotype']=labels;tissue['cell_region']=display_region;tissue['islet_name']=display_name
    tissue.to_csv(out/'donors'/f'{cid}.csv',index=False)
    distance=pd.DataFrame(dict(case_id=int(cid),donor_status=meta['Donor Status'],cell_region=display_region,phenotype=labels,signed_dist=signed.astype(np.float32)))
    for marker in POS_MARKERS_CODEX:distance[marker+'_pos']=(cells.X[:,list(cells.var_names).index(marker)]>.5).astype(int)
    kde=distance.loc[display_region!='tissue',['case_id','donor_status','cell_region','phenotype','signed_dist']].copy()
    kde['is_immune']=kde.phenotype.isin(IMMUNE_TYPES)
    with gzip.open(snapshot['geometry'],'rt') as stream:source_geo=json.load(stream)['features']
    geometry_features.extend(f for f in source_geo if f['properties'].get('classification',{}).get('name')!='Islet')
    with gzip.open(out/'gson'/f'{int(cid):04d}.geojson.gz','wt') as stream:json.dump(dict(type='FeatureCollection',features=geometry_features),stream)
    cells.write_h5ad(out/'additions'/f'{cid}.h5ad',compression='gzip')
    pd.DataFrame(memberships).to_parquet(out/'additions'/f'{cid}-memberships.parquet',index=False)
    a=ad.AnnData(np.vstack(means),obs=pd.DataFrame(observations).set_index('islet_id',drop=False),var=cells.var.copy())
    a.obsm['phenotype_counts']=np.asarray(counts,dtype=np.float32)
    a.obsm['phenotype_proportions']=a.obsm['phenotype_counts']/a.obs.total_cells.to_numpy()[:,None]
    if morphology_names:a.obsm['morphology']=np.asarray(morph,dtype=np.float32)
    return a,dict(composition=composition,targets=targets,markers=markers,features=features,lookup=lookup,
                  unpaired=unpaired,distance=distance,kde=kde)


def append_h5ad(baseline, additions):
    """Append without recomputing any original source array or metadata."""
    additions=additions[:,baseline.var_names].copy()
    original=baseline.obs.copy()
    for col in original:
        if isinstance(original[col].dtype,pd.CategoricalDtype):original[col]=original[col].astype(object)
    obs=pd.concat([original,additions.obs],axis=0,sort=False)
    obs.index.name=baseline.obs.index.name
    for col in obs:
        if obs[col].dtype==object and col!=obs.index.name:obs[col]=pd.Categorical(obs[col])
    result=ad.AnnData(np.vstack([baseline.X,additions.X]),obs=obs,var=baseline.var.copy())
    result.uns=copy.deepcopy(baseline.uns)
    for key,array in baseline.obsm.items():
        extra=additions.obsm.get(key,np.full((additions.n_obs,array.shape[1]),np.nan))
        if extra.shape[1:]!=array.shape[1:]:raise ValueError(f'Incompatible preserved {key}')
        result.obsm[key]=np.concatenate([array,extra])
    for key,array in baseline.layers.items():
        result.layers[key]=np.vstack([array,np.full(additions.shape,np.nan)])
    for key,array in baseline.varm.items():result.varm[key]=array.copy()
    for key,array in baseline.varp.items():result.varp[key]=array.copy()
    # Old graphs are historical derived fields, not padded into a new graph.
    return result


def append_csv(path,rows):
    columns=pd.read_csv(path,nrows=0).columns
    new=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    # Append bytes, avoiding any float round-trip on original records.
    new.reindex(columns=columns).to_csv(path,mode='a',header=False,index=False)


def append_groovy(adata,kind,rows):
    prefix='groovy_'+kind+'_';columns=list(adata.uns[prefix+'columns'])
    new=pd.DataFrame(rows)
    new['source_file']='reviewed_protein_expansion'
    for col in columns:
        old=np.asarray(adata.uns[prefix+col]);extra=new[col] if col in new else pd.Series([np.nan]*len(new))
        if old.dtype.kind in 'SUO':extra=extra.fillna('').astype(str).to_numpy();old=old.astype(str)
        else:extra=pd.to_numeric(extra,errors='raise').to_numpy(dtype=float)
        adata.uns[prefix+col]=np.concatenate([old,extra])
    adata.uns[prefix+'n_rows']=len(adata.uns[prefix+columns[0]])


def assemble(args):
    baseline=args.baseline.resolve(strict=True);prepared_root=args.prepared.resolve(strict=True)
    prepared,geometry,review=validate_review(prepared_root,args.geometry,args.acceptance)
    out=args.output.resolve();reports=args.report_root.resolve()
    if out.exists() or reports.exists():raise FileExistsError('Use fresh output and report roots')
    if out==baseline or reports.is_relative_to(out):raise ValueError('Unsafe output or report root')
    cohort=pd.read_csv(args.cohort,dtype=str)
    assert_no_excluded(cohort.case_id,'Cohort manifest')
    if set(cohort.case_id)!=set(DONORS)|set(ad.read_h5ad(baseline/'islet_explorer.h5ad').obs.imageid.astype(str)):
        raise ValueError(f'Expected original 15 plus exactly the {len(DONORS)} additions')
    clinical=pd.read_csv(args.clinical,dtype={'case_id':str})
    frequencies=reference_frequencies(prepared['original_source'])
    copy_baseline(baseline,out,skip=['islet_explorer.h5ad','adata_ins_root.h5ad'])
    (out/'additions').mkdir();reports.mkdir(parents=True)
    for cid in DONORS:
        shutil.copy2(prepared_root/f'{cid}-census.parquet',out/'additions'/f'{cid}-census.parquet')
    shutil.copy2(args.acceptance,out/'owner_acceptance.json')
    shutil.copy2(RULES,out/'protein_phenotype_rules.csv')
    shutil.copy2(args.geometry/'manifest.json',out/'boundary_sources.json')
    (out/'images').mkdir()
    image_inventory=[]
    for row in cohort.itertuples(index=False):
        source=Path(row.image).resolve(strict=True)
        (out/'images'/Path(row.image).name).symlink_to(source)
        image_inventory.append(dict(case_id=row.case_id,image=Path(row.image).name,
            source=str(source),bytes=source.stat().st_size))
    (out/'image_sources.json').write_text(json.dumps(image_inventory,indent=2)+'\n')
    base=ad.read_h5ad(baseline/'islet_explorer.h5ad');donors=[];tables=[]
    for cid in DONORS:
        cells=ad.read_h5ad(prepared_root/f'{cid}.h5ad')
        row=cohort.loc[cohort.case_id.eq(cid)].iloc[0]
        raw=ad.read_h5ad(row.cells)
        if sha256(row.cells)!=prepared['raw_source_sha256'][cid]:raise ValueError('Raw source changed')
        snapshot=next(r for r in geometry['donors'] if r['case_id']==cid)
        a,t=process_donor(cid,cells,raw,snapshot,get_metadata(cid,cohort,clinical,cells),frequencies,out,
                          list(base.uns['morphology_names']))
        donors.append(a);tables.append(t);print(f'{cid}: {a.n_obs} paired islets assembled',flush=True)
    additions=ad.concat(donors,join='outer',merge='same')
    combined=append_h5ad(base,additions)
    for kind in ['composition','targets','markers']:append_groovy(combined,kind,[r for t in tables for r in t[kind]])
    for key in ['features','lookup']:
        append_csv(out/({'features':'barlow_islet_features.csv','lookup':'islet_spatial_lookup.csv'}[key]),[r for t in tables for r in t[key]])
    for key,file in [('distance','cell_distance.csv'),('kde','signed_distance_kde.csv')]:append_csv(out/file,pd.concat([t[key] for t in tables],ignore_index=True))
    append_csv(out/'donor_metadata.csv',clinical.loc[clinical.case_id.isin(DONORS)])
    nbr_columns=pd.read_csv(baseline/'neighborhood_metrics.csv',nrows=0).columns
    append_csv(out/'neighborhood_metrics.csv',additions.obs.reindex(columns=nbr_columns))
    if args.donor_correction=='harmony':
        uncorrected,uncorrected_diagnostics=fit_embedding(combined,donor_correction='none')
        uncorrected.write_h5ad(reports/'uncorrected_protein_pca.h5ad',compression='gzip')
        (reports/'uncorrected_embedding.json').write_text(json.dumps(uncorrected_diagnostics,indent=2)+'\n')
        del uncorrected
    combined,diagnostics=fit_embedding(combined,donor_correction=args.donor_correction)
    combined.uns['protein_expansion']=dict(rules_sha256=sha256(RULES),acceptance_sha256=sha256(args.acceptance),
        reference_frequencies=frequencies,original_cells_sha256=prepared['original_source_sha256'],
        source_phenotypes='original labels preserved; additions use reviewed reconstruction',
        original_single_cell_partition=prepared['original_source'],scvi_required=False)
    combined.write_h5ad(out/'islet_explorer.h5ad',compression='gzip')
    trajectory=append_h5ad(ad.read_h5ad(baseline/'adata_ins_root.h5ad'),additions)
    for key in ['leiden_0.3','leiden_0.5','leiden_0.8','leiden_1.0','leiden_umap_1','leiden_umap_2']:
        trajectory.obs[key]=combined.obs[key].copy()
    for key in ['X_umap','X_pca','X_protein_pca']:trajectory.obsm[key]=combined.obsm[key].copy()
    for key in ['connectivities','distances']:trajectory.obsp[key]=combined.obsp[key].copy()
    for key in ['neighbors','umap','pca','protein_embedding']:trajectory.uns[key]=copy.deepcopy(combined.uns[key])
    trajectory.varm['PCs']=combined.varm['PCs'].copy()
    trajectory.write_h5ad(out/'adata_ins_root.h5ad',compression='gzip')
    # No stale original-only Excel fallback is copied into an expanded root.
    with pd.ExcelWriter(out/'master_results.xlsx') as writer:
        for kind,sheet in [('composition','Islet_Composition'),('targets','Islet_Targets'),('markers','Islet_Markers'),('markers','LGALS3')]:
            old=pd.read_excel(baseline/'master_results.xlsx',sheet_name=sheet)
            extra=pd.DataFrame([r for t in tables for r in t[kind]])
            if kind=='markers':extra=extra.loc[extra.marker.eq('LGALS3') if sheet=='LGALS3' else ~extra.marker.eq('LGALS3')]
            extra['source_file']='reviewed_protein_expansion'
            pd.concat([old,extra.reindex(columns=old.columns)],ignore_index=True).to_excel(writer,sheet_name=sheet,index=False)
    pd.DataFrame([r for t in tables for r in t['unpaired']]).to_csv(reports/'unpaired_islets.csv',index=False)
    pd.DataFrame([r for t in tables for r in t['composition']]).to_csv(reports/'raw_composition_census.csv',index=False)
    (reports/'embedding.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    parity=audit(baseline,out)
    (reports/'original_subset_parity.json').write_text(json.dumps(parity,indent=2)+'\n')
    if not parity['passed']:raise ValueError('Original source-artifact preservation failed; do not publish')
    if sha256(prepared['original_source'])!=prepared['original_source_sha256']:raise ValueError('Original source changed')
    result=dict(status='assembled_pending_browser_verification',donors=len(cohort),islets=combined.n_obs,
        original_source_unchanged=True,original_single_cell_storage=f'immutable source partition plus {len(DONORS)} additions',
        original_subset_parity=True,scvi_required=False,donor_correction=args.donor_correction,
        zero_endocrine_islets=int(pd.read_csv(out/'barlow_islet_features.csv').n_endocrine_core.eq(0).sum()))
    (out/'protein_release.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','prepared','geometry','acceptance','cohort','clinical','output','report-root']:
        p.add_argument('--'+key,type=Path,required=True)
    # 2026-10-03 owner decision: Harmony is the release method. With harmony the
    # uncorrected fit is also written to the reports directory for comparison.
    p.add_argument('--donor-correction',choices=['none','harmony'],default='harmony')
    assemble(p.parse_args())
