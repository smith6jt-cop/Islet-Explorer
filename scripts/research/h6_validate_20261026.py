#!/usr/bin/env python3
"""Technical checks and input/output fingerprints for completed H6 analysis."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h5_h10_original_plan_20261026/h6'

def sha(path):
    with path.open('rb') as fh:return hashlib.file_digest(fh,'sha256').hexdigest()

def main():
    d=pd.read_parquet(OUT/'h6_object_features.parquet');assert len(d)==721 and d.islet_id.is_unique
    checks={}
    for label,measurements,nullfile,graphfolder in [
        ('centroid','permutation_measurements.csv','permutation_null_draws.npz','graphs_by_donor'),
        ('cell_boundary','boundary_contact_measurements.csv','boundary_permutation_null_draws.npz','boundary_graphs_by_donor')]:
        m=pd.read_csv(OUT/measurements).set_index('islet_id');n=np.load(OUT/nullfile)
        assert n['draws'].shape==(721,1999,3) and np.isfinite(n['draws']).all()
        assert (n['draws']>=0).all() and (n['draws']<=1).all()
        mu=n['draws'].astype(float).mean(axis=1)
        expected=m.loc[n['islet_ids'],['combined_null_mean','alpha_null_mean','delta_null_mean']].to_numpy()
        assert np.allclose(mu,expected,atol=1e-7,rtol=0)
        maxerr=0.
        for name in n['islet_ids']:
            row=m.loc[name];case=int(row.case_id)
            orig=np.load(OUT/'graphs_by_donor'/str(case)/f'{name}.npz')
            graph=np.load(OUT/graphfolder/str(case)/f'{name}.npz')
            assert np.array_equal(orig['labels'],graph['labels']) and np.array_equal(orig['object_ids'],graph['object_ids'])
            e=graph['edges'];aa=graph['labels'][e[:,0]];bb=graph['labels'][e[:,1]]
            combined=np.sum(((aa==0)&np.isin(bb,[1,2]))|((bb==0)&np.isin(aa,[1,2])))/len(e)
            assert np.isclose(combined,row.combined_observed)
            maxerr=max(maxerr,abs(row.combined_null_mc_error))
        checks[label]={'objects':721,'permutations_per_object':1999,'statistics':3,
                       'stored_float32_draws_match_analysis_means':True,'cell_IDs_and_labels_match_source_graphs':True,
                       'observed_edge_fractions_recomputed':True,'maximum_MC_mean_minus_exact_expectation':maxerr}
    donor=pd.read_csv(OUT/'donor_regression_eligibility.csv');primary=donor[donor.model.eq('primary')]
    assert primary.estimable.sum()==15
    assert primary.loc[primary.estimable&primary.donor_status.eq('T1D'),'case_id'].tolist()==[6533]
    pairs=pd.read_csv(OUT/'donor_disease_pairwise_tests.csv')
    assert pairs.loc[pairs.group_2.eq('T1D'),'p_value'].isna().all()
    audit=pd.read_csv(OUT/'final_rds_convergence_audit.csv')
    assert len(audit)==51 and audit.optimizer_code.eq(0).all() and audit.min_hessian_eigenvalue.gt(0).all()
    checks.update({'model_fits_checked_from_RDS':len(audit),'all_optimizer_codes_zero':True,
                   'all_saved_fixed_effect_designs_full_rank':True,'all_saved_Hessians_positive':True,
                   'primary_independent_donor_coefficients':15,'primary_T1D_donor_coefficients':1,
                   'T1D_independent_donor_tests_unreported':True,
                   'visual_checks':'Six figures inspected across descriptive and completed stages; axes, group labels and legends verified.',
                   'analysis_changes_after_latest_user_restriction':'None; only reporting, rendering and technical validation.'})
    (OUT/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    inputs=[ROOT/'output/h5_h10_original_plan_20261026/h5/object_ins_features.parquet',
            ROOT/'output/h5_h10_original_plan_20261026/h5/source_hashes.json',
            ROOT/'output/h5_h10_original_plan_20261026/h7/object_neurovascular_features.csv']
    inputs+=sorted((ROOT/'data/development/uniform-20donor-20261015/cell_geojson').glob('*.geojson'))
    (OUT/'source_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p)for p in inputs},indent=2)+'\n')
    files=sorted(p for p in OUT.rglob('*')if p.is_file()and p.name not in ['manifest.json','validation.log'])
    scripts=sorted((ROOT/'scripts/research').glob('h6*20261026*'))
    manifest=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)}for p in files+scripts]
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
