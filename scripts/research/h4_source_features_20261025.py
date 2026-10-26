#!/usr/bin/env python3
"""Compute H4 source allocation from existing v10 cell measurements, without retyping."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import io
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/development/uniform-20donor-20261022-v10'
OUT = ROOT / 'output/h1_h4_literature_methods_20261025/h4'
BASE = ROOT / 'output/h1_edge_corrected_20261024/h1_object_features.parquet'
PRIMARY = ROOT / 'output/h1_edge_corrected_20261024/models/primary_analysis_objects.csv'
ENDO = ['Beta cell', 'Alpha cell', 'Delta cell', 'CHGA+ hormone-negative', 'Endocrine']

def extract(rec):
    path = DATA / 'cells' / (rec['islet_id'] + '.csv')
    raw = path.read_bytes()
    d = pd.read_csv(io.BytesIO(raw), float_precision='round_trip')
    assert d.cell_region.eq('core').sum() == rec['n_cells_core']
    assert d.cell_region.eq('peri').sum() == rec['n_cells_peri']
    assert d.loc[d.cell_region.eq('core'), 'phenotype'].eq('Beta cell').sum() == rec['n_beta_core']
    vals = d[['LGALS3','ECAD','panCK','CK19']].to_numpy(float)
    assert not ((vals[np.isfinite(vals)] < 0) | (vals[np.isfinite(vals)] > 1)).any()
    macrophage = d.phenotype.eq('Macrophage')
    epithelial = d.phenotype.eq('Exocrine') & d[['ECAD','panCK']].ge(.10).any(axis=1)
    ck19 = d.phenotype.eq('Exocrine') & d.CK19.ge(.10)
    all_exocrine = d.phenotype.eq('Exocrine')
    endocrine = d.phenotype.isin(ENDO)
    rows = []
    for variant, epi, cut in [('continuous',epithelial,None),('lgals3_005',epithelial,.05),
                              ('lgals3_010',epithelial,.10),('lgals3_020',epithelial,.20),
                              ('ck19_epithelial',ck19,None),('all_exocrine',all_exocrine,None)]:
        row = dict(islet_id=rec['islet_id'], variant=variant)
        for name, mask in [('myeloid',macrophage),('epithelial',epi),('endocrine',endocrine)]:
            v = d.loc[mask,'LGALS3']
            row['n_'+name+'_total'] = int(mask.sum())
            row['n_'+name] = int(v.notna().sum())
            row['n_'+name+'_missing'] = int(v.isna().sum())
            row['n_'+name+'_positive010'] = int(v.ge(.10).sum())
            # Source summaries contain measured cells only; missing source-cell
            # means are omitted, never treated as measured zeros. Their count
            # remains explicit and a complete-pool sensitivity is fitted.
            measured = v.dropna()
            row['signal_'+name] = (float(measured.sum()) if cut is None else
                                   float(measured.where(measured.ge(cut),0).sum()))
            row['mean_cell_area_'+name] = (float(d.loc[mask & d.LGALS3.notna(),'Cell Area'].mean())
                                            if len(measured) else np.nan)
        pool = row['signal_myeloid'] + row['signal_epithelial']
        row['signal_pool'] = pool
        row['myeloid_fraction'] = row['signal_myeloid']/pool if np.isfinite(pool) and pool>0 else np.nan
        row['n_epithelial_ck19_supported'] = int((epi & d.CK19.ge(.10)).sum())
        row['source_cells_overlap'] = int((macrophage & epi).sum())
        row['source_pool_complete'] = row['n_myeloid_missing'] + row['n_epithelial_missing'] == 0
        rows.append(row)
    return rows, {str(path.relative_to(ROOT)):hashlib.sha256(raw).hexdigest()}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    d = pd.read_parquet(BASE)
    d = d.loc[d.core_complete & d.neighborhood_complete_20 & d.neighborhood_complete_100].copy()
    assert d.case_id.nunique()==20 and not d.case_id.astype(str).eq('78').any()
    primary = pd.read_csv(PRIMARY)
    d['h1_primary'] = d.islet_id.isin(primary.islet_id)
    rows, hashes = [], {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,(r,h) in enumerate(pool.map(extract,d.to_dict('records'),chunksize=20)):
            rows.extend(r); hashes.update(h)
            if (i+1)%500==0: print(f'{i+1}/{len(d)} objects',flush=True)
    f = pd.DataFrame(rows).merge(d,on='islet_id',validate='many_to_one')
    f['neighborhood_area_mm2'] = (f.core_area_um2+f.peri_area_um2)/1e6
    f['beta_density_mm2'] = f.n_beta_core/f.core_area_um2*1e6
    for label in ['myeloid','epithelial','endocrine']:
        f[label+'_density_mm2'] = f['n_'+label]/f.neighborhood_area_mm2
        f['signal_'+label+'_density_mm2'] = f['signal_'+label]/f.neighborhood_area_mm2
    f['signal_pool_density_mm2'] = f.signal_pool/f.neighborhood_area_mm2
    f['immune_density_mm2'] = f.n_immune_nearby/f.neighborhood_area_mm2
    f['source_pool_nonzero'] = f.signal_pool.gt(0)&np.isfinite(f.signal_pool)
    assert f.source_cells_overlap.eq(0).all()
    f.to_csv(OUT/'object_source_features.csv',index=False)
    f.to_parquet(OUT/'object_source_features.parquet',index=False)
    audit=f.groupby(['variant','donor_status']).agg(objects=('islet_id','size'),
        assessable=('source_pool_nonzero','sum'),donors=('case_id','nunique'),
        myeloid_cells=('n_myeloid','sum'),epithelial_cells=('n_epithelial','sum'),
        myeloid_missing=('n_myeloid_missing','sum'),epithelial_missing=('n_epithelial_missing','sum'),
        myeloid_signal=('signal_myeloid','sum'),epithelial_signal=('signal_epithelial','sum'))
    audit.to_csv(OUT/'measurement_summary.csv')
    for p in [BASE,PRIMARY]:hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    (OUT/'source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(audit.to_string())

if __name__=='__main__':main()
