#!/usr/bin/env python3
"""Full-census, zero-inclusive tissue-grid analysis inputs for original H8.

No new segmentation, typing, app materialization or image thresholding. Counts
reuse source object UUIDs and the corrected counting/edge audit. Tissue fields
use existing aligned/scaled marker measurements, never the Beta endpoint.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
from shapely.affinity import scale
import h1_features_20261022 as hf
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'uniform/typing'))
import codex_typing as ct
import classify as cl

SRC=Path('data/development/uniform-20donor-20261015')
APP=Path('data/development/uniform-20donor-20261022-v10')
COR=Path('output/h1_edge_corrected_20261024')
OUT=Path('output/h5_h10_original_plan_20261026/h8')
MARKERS=['CK19','panCK','SMA','VIM','ColIV']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def source_objects(tag,case,eligible,paired):
    p=SRC/'islets'/f'{tag}.geojson'
    geoms={};measurements={}
    for f in json.loads(p.read_text())['features']:
        if f['id'] not in set(eligible.object_id):continue
        geoms[f['id']]=shapely.make_valid(scale(shape(f['geometry']),xfact=hf.PIXEL_UM,yfact=hf.PIXEL_UM,origin=(0,0)))
        measurements[f['id']]=(f['properties'].get('measurements') or {}).get('Whole cells',np.nan)
    audit=pd.read_parquet(COR/'eligibility/cells'/f'{case}.parquet')
    measured=audit.served & audit.edge_analysis_eligible
    audit['measured']=measured
    audit['beta']=measured & audit.phenotype.eq('Beta cell')
    audit['endocrine']=measured & audit.phenotype.isin(hf.ENDOCRINE)
    audit['immune']=measured & audit.phenotype.isin(hf.IMMUNE)
    counts=audit.loc[audit.endocrine_object_id.notna()].groupby('endocrine_object_id').agg(
        source_cells=('object_id','size'),measured_cells=('measured','sum'),n_beta=('beta','sum'),
        n_endocrine=('endocrine','sum'),n_immune=('immune','sum'))
    objects=eligible.merge(counts,left_on='object_id',right_index=True,how='left',validate='one_to_one')
    for col in counts:objects[col]=objects[col].fillna(0).astype(int)
    objects['annotation_whole_cells']=objects.object_id.map(measurements)
    objects['unknown_cells']=np.maximum(objects.source_cells-objects.measured_cells,
                                        objects.annotation_whole_cells.fillna(0)-objects.measured_cells).astype(int)
    objects['islet_id']=str(case)+'_'+objects.object_name
    objects['paired_row_available']=objects.islet_id.isin(paired.islet_id)
    objects['centroid_x_um']=objects.object_id.map(lambda k:geoms[k].centroid.x)
    objects['centroid_y_um']=objects.object_id.map(lambda k:geoms[k].centroid.y)
    cluster=objects.object_class.eq('EndocrineCluster')
    small=objects.n_endocrine.between(2,9)
    beta=objects.n_beta.gt(0)
    objects['small_beta_cluster_observed']=cluster&small&beta
    objects['small_beta_cluster_lower']=cluster&objects.n_endocrine.ge(2)&(objects.n_endocrine+objects.unknown_cells).le(9)&beta
    objects['small_beta_cluster_upper']=cluster&objects.n_endocrine.le(9)&(objects.n_endocrine+objects.unknown_cells).ge(2)&(beta|objects.unknown_cells.gt(0))
    objects['small_beta_cluster_wholecells']=cluster&objects.source_cells.between(2,9)&beta
    objects['small_beta_cluster_area']=cluster&objects.core_area_um2.lt(1000)&beta
    objects['small_beta_islet']=objects.object_class.eq('Islet')&small&beta
    objects['large_beta_islet']=objects.object_class.eq('Islet')&objects.n_endocrine.ge(10)&beta
    objects['beta_single']=objects.object_class.eq('EndocrineSingle')&beta
    # Direct UUID source counts must agree with the paired corrected features.
    check=objects.merge(paired[['islet_id','n_beta_core','n_endocrine_core','n_cells_core']],on='islet_id')
    discrepancies=check.loc[(check.n_beta!=check.n_beta_core)|(check.n_endocrine!=check.n_endocrine_core)|(check.measured_cells!=check.n_cells_core)]
    if len(discrepancies):
        discrepancies.to_csv(OUT/f'count_discrepancies_{case}.csv',index=False)
        raise ValueError(f'{case}: {len(discrepancies)} source counts differ from corrected feature table')
    return objects,geoms,audit


def run():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'by_donor').mkdir(exist_ok=True)
    manifest=ct.manifest(SRC);scale_table=cl.load_scale(SRC)
    edge=pd.read_csv(COR/'eligibility/object_eligibility.csv')
    paired=pd.read_parquet(COR/'h1_object_features.parquet')
    hashes={str(COR/'eligibility/object_eligibility.csv'):sha(COR/'eligibility/object_eligibility.csv'),
            str(COR/'h1_object_features.parquet'):sha(COR/'h1_object_features.parquet'),
            str(SRC/'typing/scale.parquet'):sha(SRC/'typing/scale.parquet')}
    object_frames=[];tile_frames=[];donor_rows=[]
    for tag,entry in sorted(manifest.items()):
        case=int(entry['case_id']);status=entry['donor_status']
        objects,geoms,audit=source_objects(tag,case,edge.loc[edge.case_id.eq(case)].copy(),paired)
        objects['donor_status']=status
        _,_,tissue,_,gp=hf.geometry(APP,case)
        plane=json.loads((SRC/'planes'/f'{tag}.planes.json').read_text())
        W=plane['width']*plane['downsample']*hf.PIXEL_UM;H=plane['height']*plane['downsample']*hf.PIXEL_UM
        acquisition=shapely.box(0,0,W,H)
        edge_objects=shapely.union_all([geoms[k] for k in objects.loc[~objects.core_complete,'object_id']])
        # Original islet construction used a 5 um acquisition-edge inset.
        valid=tissue.intersection(acquisition.buffer(-5)).difference(edge_objects)
        objects['centroid_in_valid_tissue']=shapely.covers(valid,shapely.points(objects.centroid_x_um,objects.centroid_y_um))
        objects['density_eligible']=objects.core_complete & objects.centroid_in_valid_tissue
        cp=COR/'analysis_cells_by_donor'/f'{case}.parquet';ap=SRC/'typing/aligned'/f'{tag}.parquet'
        cells=pd.read_parquet(cp)
        cells=cells.loc[cells.edge_analysis_eligible & shapely.covers(valid,shapely.points(cells.x_um,cells.y_um))].copy()
        aligned=pd.read_parquet(ap,columns=['object_id',*MARKERS])
        aligned[MARKERS]=cl.apply_scale(aligned[MARKERS].to_numpy(),MARKERS,scale_table)
        cells=cells.merge(aligned,on='object_id',how='left',validate='one_to_one',indicator=True)
        if cells._merge.ne('both').any():raise ValueError(f'{case}: UUID missing from aligned source')
        extra=set(audit.loc[audit.endocrine_object_id.isna(),'object_id'])
        cells['extraendocrine']=cells.object_id.isin(extra)
        # Continuous cell-associated signals, not positivity gates or mass.
        cells['field_marker_measured']=cells[MARKERS].notna().all(axis=1)&cells.extraendocrine
        cells['field_marker_missing']=~cells[MARKERS].notna().all(axis=1)&cells.extraendocrine
        cells['epithelial_signal']=cells[['CK19','panCK']].mean(axis=1).fillna(0)*cells.field_marker_measured
        cells['stromal_signal']=cells[['SMA','VIM','ColIV']].mean(axis=1).fillna(0)*cells.field_marker_measured
        cells['immune']=cells.phenotype.isin(hf.IMMUNE)
        donorrec=dict(case_id=case,donor_status=status,geometric_objects=len(objects),
          paired_objects=int(objects.paired_row_available.sum()),unpaired_objects=int((~objects.paired_row_available&objects.object_class.ne('EndocrineSingle')).sum()),
          incomplete_objects=int((~objects.core_complete).sum()),eligible_objects=int(objects.density_eligible.sum()),
          full_tissue_area_mm2=tissue.area/1e6,valid_tissue_area_mm2=valid.area/1e6,
          valid_tissue_cells=len(cells),unmeasured_core_objects=int(objects.measured_cells.eq(0).sum()))
        classes=['small_beta_cluster_observed','small_beta_cluster_lower','small_beta_cluster_upper',
                 'small_beta_cluster_wholecells','small_beta_cluster_area','small_beta_islet','large_beta_islet','beta_single']
        for cls in classes:donorrec[cls]=int(objects.loc[objects.density_eligible,cls].sum())
        donor_rows.append(donorrec)
        for width in (100,200,400):
            nx=int(np.ceil(W/width));ny=int(np.ceil(H/width))
            xx,yy=np.meshgrid(np.arange(nx),np.arange(ny));ix=xx.ravel();iy=yy.ravel()
            boxes=shapely.box(ix*width,iy*width,(ix+1)*width,(iy+1)*width)
            areas=shapely.area(shapely.intersection(boxes,valid))
            fullareas=shapely.area(shapely.intersection(boxes,tissue))
            keep=areas>1e-8
            frame=pd.DataFrame(dict(case_id=case,donor_status=status,grid_um=width,tile_number=np.arange(nx*ny),
                tile_x=ix,tile_y=iy,x_min_um=ix*width,y_min_um=iy*width,
                valid_area_um2=areas,full_tissue_area_um2=fullareas,acquisition_edge_tile=(ix==0)|(iy==0)|((ix+1)*width>W)|((iy+1)*width>H)))
            ci=(np.floor(cells.y_um/width).astype(int)*nx+np.floor(cells.x_um/width).astype(int)).to_numpy()
            for name,weights in [('measured_cells',None),('immune_cells',cells.immune),('extraendocrine_cells',cells.extraendocrine),
                                 ('field_marker_measured_cells',cells.field_marker_measured),('field_marker_missing_cells',cells.field_marker_missing),
                                 ('epithelial_signal',cells.epithelial_signal),('stromal_signal',cells.stromal_signal)]:
                frame[name]=np.bincount(ci,weights=None if weights is None else weights.to_numpy(),minlength=nx*ny)
            oi=(np.floor(objects.centroid_y_um/width).astype(int)*nx+np.floor(objects.centroid_x_um/width).astype(int)).to_numpy()
            good=objects.density_eligible.to_numpy()
            for name in classes:
                frame[name]=np.bincount(oi[good],weights=objects.loc[good,name].to_numpy(int),minlength=nx*ny).astype(int)
                if int(frame[name].sum())!=donorrec[name]:raise ValueError('Tile count does not conserve donor census')
            for name,cond in [('all_geometric_objects',np.ones(len(objects),bool)),('edge_excluded_objects',~objects.core_complete),
                              ('unmeasured_core_objects',objects.measured_cells.eq(0)),('unpaired_objects',~objects.paired_row_available&objects.object_class.ne('EndocrineSingle'))]:
                frame[name]=np.bincount(oi[cond],minlength=nx*ny)
            frame=frame.loc[keep].copy()
            frame['tile_id']=str(case)+'_'+str(width)+'_'+frame.tile_x.astype(str)+'_'+frame.tile_y.astype(str)
            for name in ['immune_cells','epithelial_signal','stromal_signal',*classes]:
                frame[name+'_density_mm2']=frame[name]*1e6/frame.valid_area_um2
            assert np.isclose(frame.valid_area_um2.sum(),valid.area,rtol=1e-8,atol=1e-4)
            tile_frames.append(frame)
        object_frames.append(objects)
        objects.to_csv(OUT/'by_donor'/f'{case}_objects.csv',index=False)
        for p in [cp,ap,gp,SRC/'islets'/f'{tag}.geojson',COR/'eligibility/cells'/f'{case}.parquet']:hashes[str(p)]=sha(p)
        print(tag,'objects',len(objects),'eligible',int(objects.density_eligible.sum()),'smallBetaClusters',donorrec['small_beta_cluster_observed'],'unpaired',donorrec['unpaired_objects'],flush=True)
    allobjects=pd.concat(object_frames,ignore_index=True);tiles=pd.concat(tile_frames,ignore_index=True)
    allobjects.to_parquet(OUT/'full_object_census.parquet',index=False);allobjects.to_csv(OUT/'full_object_census.csv',index=False)
    tiles.to_parquet(OUT/'all_tissue_tiles.parquet',index=False);tiles.to_csv(OUT/'all_tissue_tiles.csv',index=False)
    pd.DataFrame(donor_rows).to_csv(OUT/'donor_census.csv',index=False)
    (OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print('Totals',len(allobjects),'objects',len(tiles),'tiles across3grids',flush=True)

if __name__=='__main__':run()
