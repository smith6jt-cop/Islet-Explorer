#!/usr/bin/env python3
"""Validate numerical reuse, measurement units and native artifact completeness."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/app_hypotheses_20261026/h6_h8'
OLD=ROOT/'output/h5_h10_original_plan_20261026'
audit={}
six=pd.read_parquet(OUT/'data/h6_objects.parquet')
old6=pd.read_csv(OLD/'h6/h6_object_features.csv')
fields=['combined_excess','combined_observed','combined_null_mean','n_beta_core','core_area_um2']
z=six.merge(old6[['islet_id',*fields]],on='islet_id',suffixes=('_new','_old'),validate='1:1')
audit['H6_source_measurement_max_abs_differences']={c:float(np.max(np.abs(z[c+'_new']-z[c+'_old']))) for c in fields}
assert all(np.allclose(z[c+'_new'],z[c+'_old'],rtol=1e-12,atol=1e-10) for c in fields)
audit['H6_percent_units_verified']=bool(np.allclose(six.observed_adjacency_percent,100*six.combined_observed) and np.allclose(six.permuted_adjacency_percent,100*six.combined_null_mean) and np.allclose(six.combined_excess_percentage_points,100*six.combined_excess))
seven=pd.read_parquet(OUT/'data/h7_objects.parquet');q=seven[seven.h7_primary_eligible]
model=pd.read_csv(OUT/'models/h7_primary_20um_objects.csv',low_memory=False)
audit['H7_cohort_matches_model']=set(q.islet_id)==set(model.islet_id)
audit['H7_by_state']=q.groupby('donor_status').agg(n_objects=('islet_id','size'),n_donors=('case_id','nunique'),zero_Beta=('n_beta_core',lambda x:int((x==0).sum()))).reset_index().to_dict('records')
audit['H7_beta_density_max_absolute_error']=float(np.max(np.abs(q.beta_density_mm2-1e6*q.n_beta_core/q.core_area_um2)))
audit['H7_outcome_max_absolute_error']=float(np.max(np.abs(q.log1p_beta_density-np.log1p(q.beta_density_mm2))))
audit['H7_union_fraction_bounds']={c:[float(q[c].min()),float(q[c].max())] for c in ['nerve_area_fraction','cd31_area_fraction','nerve_near_cd31_20um']}
census=pd.read_parquet(OUT/'data/h8_census_with_states.parquet');tiles=pd.read_csv(OUT/'data/h8_tiles_with_states.csv')
oldtiles=pd.read_csv(OLD/'h8/all_tissue_tiles.csv');oldtiles=oldtiles[oldtiles.grid_um.eq(200)]
audit['H8_saved_tile_fields_unchanged']=True
for c in oldtiles:
 if c=='tile_id':continue
 a=tiles.set_index('tile_id').sort_index()[c];b=oldtiles.set_index('tile_id').sort_index()[c]
 if not a.equals(b):
  if not (pd.api.types.is_numeric_dtype(a) and np.allclose(a,b,equal_nan=True)):
   audit['H8_saved_tile_fields_unchanged']=False
audit['H8_census_objects']=len(census);audit['H8_all_tiles']=len(tiles)
audit['H8_eligible_objects']=int(census.density_eligible.sum())
audit['H8_nonembedded_eligible_objects']=int((census.density_eligible&~census.in_composition_embedding).sum())
audit['H8_linked_tile_count_equals_eligible_census']=int(tiles.linked_eligible_objects.sum())==int(census.density_eligible.sum())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit['H8_stats_tables_copied_unchanged']={p.name:sha(p)==sha(OLD/'h8'/p.name.removeprefix('h8_')) for p in (OUT/'models').glob('h8_*.csv')}
audit['quantitative_INS_columns']={p.name:[c for c in pd.read_csv(p,nrows=0) if any(k in c.lower() for k in ['ins_','insulin','mean_ins','ins.mean'])] for p in (OUT/'data').glob('*.csv')}
assert not any(audit['quantitative_INS_columns'].values())
audit['primary_comparisons_matches_tests']=(OUT/'primary_comparisons.csv').read_bytes()==(OUT/'primary_tests.csv').read_bytes()
manifest=json.loads((OUT/'figure_manifest.json').read_text());hashes={}
for r in manifest['figures']:
 for p in [r['png'],r['pdf'],r['svg'],r['legend'],*r['data_csv']]:
  p=Path(p);assert p.exists() and p.stat().st_size>0,str(p);hashes[str(p.relative_to(ROOT))]=sha(p)
audit['native_figure_sets']=len(manifest['figures'])
audit['all_manifest_artifacts_exist']=True
assert audit['H6_percent_units_verified'] and audit['H7_cohort_matches_model'] and audit['H8_saved_tile_fields_unchanged'] and audit['H8_linked_tile_count_equals_eligible_census']
assert all(audit['H8_stats_tables_copied_unchanged'].values())
(OUT/'numerical_QA.json').write_text(json.dumps(audit,indent=2)+'\n')
(OUT/'artifact_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(json.dumps({k:v for k,v in audit.items() if k!='quantitative_INS_columns'},indent=2))
