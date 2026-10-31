#!/usr/bin/env python3
"""Bind non-INS H6/H8 measurements and extend the existing H7 ROI coverage."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,hashlib
import numpy as np
import pandas as pd
import shapely
import h7_nerve_features_20261026 as h7
import app_hypotheses_shared_20261026 as shared
ROOT=shared.ROOT;OUT=ROOT/'output/app_hypotheses_20261026/h6_h8';OLD=ROOT/'output/h5_h10_original_plan_20261026'
BASE=ROOT/'output/app_hypotheses_20261026/shared/object_measurements.parquet'
SRC=ROOT/'data/development/uniform-20donor-20261015';APP=ROOT/'data/development/uniform-20donor-20261022-v10'
MEAS=['peri_assessed_area_um2','nerve_area_um2','cd31_area_um2','nerve_area_fraction','cd31_area_fraction',
 'gap43_supported_nerve_fraction','gap43_supported_nerve_area_um2','nerve_vessel_overlap',
 'mean_nerve_vessel_distance_capped20_um','proximity_available']+[f'{p}nerve_near_cd31_{r}um' for r in [5,10,20] for p in ['', 'gap43_']]

def extra(task):
 tag,entry,rows=task;case=int(entry['case_id']);path=OUT/'h7_extra_by_donor'/f'{case}.csv'
 if path.exists():return case,len(pd.read_csv(path))
 cached=json.loads((OLD/'h7/by_donor'/f'{case}.json').read_text())
 cores,_,tissue,_,gp=h7.hf.geometry(APP,case)
 endo=h7.geoms(SRC/'islets'/f'{tag}.geojson');etree=shapely.STRtree(endo)
 vessels=h7.geoms(SRC/'structures/Capillary'/f'Capillary__{tag}.geojson');vtree=shapely.STRtree(vessels)
 img=h7.Channels(entry['image_path']);names=h7.ct.channels(SRC,entry['entry_id'])
 b3=img.arr[names.index('B3TUBB')] if img.memmap else img(names.index('B3TUBB'))
 gap=img.arr[names.index('GAP43')] if img.memmap else img(names.index('GAP43'))
 records=[]
 for r in rows.itertuples():
  result,_=h7.extract(r,cores[r.base_islet_id][1],tissue,endo,etree,vessels,vtree,b3,gap,
      cached['b3tubb_threshold_raw'],cached['gap43_threshold_raw']);records.append(result)
 pd.DataFrame(records).to_csv(path,index=False)
 path.with_suffix('.json').write_text(json.dumps({'case_id':case,'objects':len(records),'cached_threshold_source':str(OLD/'h7/by_donor'/f'{case}.json'),'b3tubb_threshold_raw':cached['b3tubb_threshold_raw'],'gap43_threshold_raw':cached['gap43_threshold_raw'],'same_existing_helper':'h7_nerve_features_20261026.extract','geometry_sha256':h7.sha(gp),'INS_channel_used':False},indent=2)+'\n')
 del b3,gap;img.close();return case,len(records)

def main():
 for s in ['data','models','figures','h7_extra_by_donor']:(OUT/s).mkdir(parents=True,exist_ok=True)
 base=shared.read(BASE)
 old=pd.read_csv(OLD/'h7/object_neurovascular_features.csv',usecols=['islet_id',*MEAS])
 missing=base[base.core_complete&base.neighborhood_complete_20&~base.islet_id.isin(old.islet_id)].copy()
 manifest=h7.ct.manifest(SRC)
 tasks=[(tag,e,missing[missing.case_id.eq(str(e['case_id']))]) for tag,e in sorted(manifest.items())]
 tasks=[x for x in tasks if len(x[2])]
 with ProcessPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(extra,t) for t in tasks]
  for f in as_completed(futures):print('H7 existing-method additional ROI',f.result(),flush=True)
 additional=pd.concat([pd.read_csv(OUT/'h7_extra_by_donor'/f'{int(e["case_id"])}.csv') for _,e,_ in tasks],ignore_index=True)
 measurements=pd.concat([old,additional[['islet_id',*MEAS]]],ignore_index=True)
 assert measurements.islet_id.is_unique
 seven=base.merge(measurements,on='islet_id',how='left',validate='one_to_one')
 seven['h7_primary_eligible']=seven.core_complete&seven.neighborhood_complete_20&seven.nerve_near_cd31_20um.notna()
 seven['log1p_beta_density']=np.log1p(seven.beta_density_mm2)
 seven['log_core_area']=np.log(seven.core_area_um2);seven['log_perimeter']=np.log(seven.core_perimeter_um)
 seven['log_nearby_immune_density']=np.log1p(1e6*seven.n_immune_nearby/(seven.core_area_um2+seven.peri_area_um2))
 assert seven.loc[seven.core_complete&seven.neighborhood_complete_20,'peri_assessed_area_um2'].notna().all()
 for field in ['nerve_area_fraction','cd31_area_fraction']:assert seven[field].dropna().between(0,1).all()
 seven.to_parquet(OUT/'data/h7_objects.parquet',index=False);seven.to_csv(OUT/'data/h7_objects.csv',index=False)
 sixraw=pd.read_csv(OLD/'h6/h6_object_features.csv',low_memory=False)
 cols=['islet_id','n_permutations','n_endocrine_nodes','n_beta_nodes','n_alpha_nodes','n_delta_nodes','n_edges','beta_fraction','delta_fraction','other_endocrine_fraction','log_core_area','log_nearby_immune','log_core_immune','graph_eligible','complete_graph']
 cols += [c for c in sixraw if c.startswith(('combined_','alpha_observed','alpha_null','alpha_excess','alpha_spatial','delta_observed','delta_null','delta_excess','delta_spatial'))]
 six=base.merge(sixraw[cols],on='islet_id',how='inner',validate='one_to_one')
 six['combined_excess_percentage_points']=100*six.combined_excess
 six['observed_adjacency_percent']=100*six.combined_observed
 six['permuted_adjacency_percent']=100*six.combined_null_mean
 six['h6_primary_eligible']=six.graph_eligible & six.combined_spatially_informative
 assert len(six)==721 and six.h6_primary_eligible.all()
 six.to_parquet(OUT/'data/h6_objects.parquet',index=False);six.to_csv(OUT/'data/h6_objects.csv',index=False)
 census=shared.read(OLD/'h8/full_object_census.parquet')
 repcols=['islet_id','umap_1','umap_2','leiden_0.5','pseudotime','pt_rank','in_composition_embedding']
 census=census.merge(base[repcols],on='islet_id',how='left',validate='one_to_one')
 census['in_composition_embedding']=census.in_composition_embedding.fillna(False).astype(bool)
 tiles=shared.read(OLD/'h8/all_tissue_tiles.csv');tiles=tiles[tiles.grid_um.eq(200)].copy()
 census['tile_x']=np.floor(census.centroid_x_um/200).astype(int);census['tile_y']=np.floor(census.centroid_y_um/200).astype(int)
 census['tile_id']=census.case_id.astype(str)+'_200_'+census.tile_x.astype(str)+'_'+census.tile_y.astype(str)
 tilefields=['tile_id','small_beta_cluster_observed_density_mm2','immune_cells_density_mm2','epithelial_signal_density_mm2','stromal_signal_density_mm2','valid_area_um2']
 census=census.merge(tiles[tilefields],on='tile_id',how='left',validate='many_to_one')
 census['matched_positive_area_tile']=census.valid_area_um2.notna()
 state_counts=census[census.density_eligible&census.matched_positive_area_tile].groupby('tile_id').agg(
  linked_eligible_objects=('islet_id','size'),embedded_objects=('in_composition_embedding','sum')).reset_index()
 state_counts['nonembedded_objects']=state_counts.linked_eligible_objects-state_counts.embedded_objects
 bystate=census[census.density_eligible&census.in_composition_embedding&census.matched_positive_area_tile].groupby(['tile_id','leiden_0.5']).size().unstack(fill_value=0)
 bystate.columns=['saved_cluster_'+str(int(c))+'_objects' for c in bystate.columns]
 tiles=tiles.merge(state_counts,on='tile_id',how='left').merge(bystate,on='tile_id',how='left')
 for c in ['linked_eligible_objects','embedded_objects','nonembedded_objects',*bystate.columns]:tiles[c]=tiles[c].fillna(0).astype(int)
 census.to_parquet(OUT/'data/h8_census_with_states.parquet',index=False);census.to_csv(OUT/'data/h8_census_with_states.csv',index=False)
 tiles.to_csv(OUT/'data/h8_tiles_with_states.csv',index=False)
 audit={'H6_objects':len(six),'H6_donors':six.case_id.nunique(),'H6_INS_availability_not_used':True,'H7_preexisting_objects':len(old),'H7_new_ROIs':len(additional),'H7_all_measured_objects':len(measurements),'H7_primary_proximity_objects':int(seven.h7_primary_eligible.sum()),'H7_zero_Beta_primary_objects':int((seven.h7_primary_eligible&seven.n_beta_core.eq(0)).sum()),'H7_donors':seven.loc[seven.h7_primary_eligible,'case_id'].nunique(),'H8_full_census_objects':len(census),'H8_tiles200':len(tiles),'H8_nonembedded_eligible_objects':int((census.density_eligible&~census.in_composition_embedding).sum()),'H8_unmatched_positive_area_tile_objects':int((~census.matched_positive_area_tile).sum()),'H7_structure_fractions':'Binary raster unions; no summed annotation area used','quantitative_INS_exported':False}
 (OUT/'input_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2),flush=True)
if __name__=='__main__':main()
