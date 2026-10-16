#!/usr/bin/env python3
"""Check a finished protein release against frozen files and raw additions."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from assemble_protein_expansion import BASELINE_DIRS, ENDOCRINE
from build_protein_embedding import sha256
from protein_phenotyping import qc_mask
from stage_original_gate_calls import DONORS
from excluded_donors import assert_no_excluded


def validate(baseline, candidate, raw_root, comparison, report):
    a=ad.read_h5ad(candidate/'islet_explorer.h5ad')
    old=ad.read_h5ad(baseline/'islet_explorer.h5ad')
    assert a.n_obs==old.n_obs+sum(a.obs.imageid.astype(str).isin([cid]).sum() for cid in DONORS)
    expected_donors=old.obs.imageid.nunique()+len(DONORS)
    assert a.obs_names.is_unique and a.obs.imageid.nunique()==expected_donors
    assert_no_excluded(a.obs.imageid.astype(str).unique(),'Release islets')
    assert_no_excluded([p.stem for p in (candidate/'donors').glob('*.csv')],'Release donor tables')
    assert np.isfinite(a.X).all()
    trajectory=ad.read_h5ad(candidate/'adata_ins_root.h5ad')
    for key in ['leiden_0.3','leiden_0.5','leiden_0.8','leiden_1.0','leiden_umap_1','leiden_umap_2']:
        np.testing.assert_array_equal(trajectory.obs[key],a.obs[key])
    np.testing.assert_allclose(a.obsm['phenotype_counts'].sum(axis=1),a.obs.total_cells,rtol=0,atol=0)
    assert json.loads((comparison/'original_subset_parity.json').read_text())['passed']
    preserved=0
    for directory in BASELINE_DIRS:
        for source in (baseline/directory).rglob('*'):
            if source.is_file():
                assert sha256(source)==sha256(candidate/directory/source.relative_to(baseline/directory)),source
                preserved+=1
    for name in ['cell_distance.csv','signed_distance_kde.csv','neighborhood_metrics.csv',
                 'donor_metadata.csv','islet_spatial_lookup.csv','barlow_islet_features.csv']:
        with (baseline/name).open('rb') as left,(candidate/name).open('rb') as right:
            while block:=left.read(8*1024*1024):assert block==right.read(len(block)),name
    features=pd.read_csv(candidate/'barlow_islet_features.csv',index_col='islet_id')
    donors=[]
    for cid in DONORS:
        raw=ad.read_h5ad(raw_root/f'{cid}.h5ad')
        cells=ad.read_h5ad(candidate/'additions'/f'{cid}.h5ad')
        assert set(raw.var_names)==set(cells.raw.var_names)
        raw=raw[:,cells.raw.var_names].copy()
        use=qc_mask(raw.X)
        census=pd.read_parquet(candidate/'additions'/f'{cid}-census.parquet')
        assert census.cell_id.astype(str).tolist()==raw.obs_names.astype(str).tolist()
        np.testing.assert_array_equal(census.analysis_included.to_numpy(dtype=bool),use)
        assert list(cells.obs_names)==list(raw.obs_names[use])
        np.testing.assert_array_equal(cells.raw.X,raw.X[use])
        np.testing.assert_allclose(cells.layers['log'],np.log1p(np.asarray(cells.raw.X,dtype=np.float64)),rtol=0,atol=0)
        assert np.isfinite(cells.X).all() and cells.obs.phenotype.notna().all()
        assert set(cells.obs.phenotype).isdisjoint({'Acinar','ECAD+','Structural'})
        part=features.loc[features.case_id.astype(str).eq(cid)]
        denom=part.n_endocrine_core.clip(lower=1)
        np.testing.assert_array_equal(part.n_endocrine_core,sum(part['n_core_'+label] for label in ENDOCRINE))
        for column in part:
            if column.startswith(('core_','peri_')):
                np.testing.assert_allclose(part[column],np.log1p(part['n_'+column]/denom),rtol=1e-13,atol=1e-13)
        donors.append(dict(case_id=cid,raw_cells=raw.n_obs,typed_cells=cells.n_obs,qc_excluded=int((~use).sum()),paired_islets=len(part),
                           phenotype_counts=cells.obs.phenotype.value_counts().to_dict()))
    u=ad.read_h5ad(comparison/'uncorrected_protein_pca.h5ad')
    np.testing.assert_array_equal(a.X,u.X)
    np.testing.assert_array_equal(a.obsm['X_pca'],u.obsm['X_pca'])
    assert a.uns['protein_embedding']['donor_correction']=='harmony'
    assert u.uns['protein_embedding']['donor_correction']=='none'
    assert not np.array_equal(a.obsm['X_protein_pca'],u.obsm['X_protein_pca'])
    images=json.loads((candidate/'image_sources.json').read_text())
    assert len(images)==expected_donors
    assert_no_excluded([row['case_id'] for row in images],'Release images')
    for row in images:
        path=candidate/'images'/row['image']
        assert path.resolve()==Path(row['source']) and path.stat().st_size==row['bytes']
    result=dict(passed=True,islets=a.n_obs,donors=donors,original_islets=old.n_obs,
        original_files_byte_identical=preserved,original_csv_prefixes_byte_identical=True,
        raw_addition_values_exact=True,composition_formulas_checked=True,
        harmony_and_uncorrected_retained=True,images=len(images))
    report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='donors'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','candidate','raw-root','comparison','report']:p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();validate(args.baseline,args.candidate,args.raw_root,args.comparison,args.report)
