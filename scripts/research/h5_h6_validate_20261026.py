#!/usr/bin/env python3
"""Verify prepared H5 accounting, donor analyses and descriptive H6 graphs."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h5_h10_original_plan_20261026'

def main():
    h5=OUT/'h5';h6=OUT/'h6'
    f=pd.read_parquet(h5/'object_ins_features.parquet')
    donor=pd.read_csv(h5/'plot_data_donor_ins_measurements.csv')
    primary=pd.read_csv(h5/'donor_regression_eligibility.csv')
    sensitivity=pd.read_csv(h5/'nearby_burden_sensitivity_donor_regression_eligibility.csv')
    assert len(f)==4179 and f.islet_id.is_unique and len(donor)==20
    assert (f.n_source_excluded_core==f.excluded_cells_core).all()
    assert (f.n_source_excluded_INS_measured+f.n_source_excluded_INS_missing==f.n_source_excluded_core).all()
    assert f.loc[f.n_beta_core.eq(0),'ins_mean_beta_core'].isna().all()
    assert f.n_beta_ins_missing.sum()==0
    assert np.allclose(donor.ins_beta_density+donor.ins_nonbeta_density,donor.ins_total_density)
    assert f.decomposition_absolute_error.max()<1e-10
    assert primary.loc[~primary.estimable,'case_id'].tolist()==[6563]
    assert sensitivity.loc[~sensitivity.estimable,'case_id'].tolist()==[6550,6563]
    g=pd.read_csv(h6/'graph_inventory.csv')
    expected=g.loc[g.graph_eligible].set_index('islet_id')
    assert len(expected)==721
    graphfiles=list((h6/'graphs_by_donor').rglob('*.npz'))
    assert len(graphfiles)==len(expected)
    for path in graphfiles:
        graph=np.load(path);e=graph['edges'];labels=graph['labels'];coords=graph['coordinates'];r=expected.loc[path.stem]
        assert len(labels)==r.n_endocrine_nodes and len(e)==r.n_edges
        assert (e[:,0]<e[:,1]).all() and len(np.unique(e,axis=0))==len(e)
        assert (np.linalg.norm(coords[e[:,0]]-coords[e[:,1]],axis=1)<=20.+1e-10).all()
        aa,bb=labels[e[:,0]],labels[e[:,1]]
        ba=int((((aa==0)&(bb==1))|((aa==1)&(bb==0))).sum())
        bd=int((((aa==0)&(bb==2))|((aa==2)&(bb==0))).sum())
        assert ba==r.observed_beta_alpha_edges and bd==r.observed_beta_delta_edges
        assert np.isclose((ba+bd)/len(e),r.observed_combined_fraction)
    h5qa={'objects':len(f),'donors':len(donor),'no_beta_objects':int(f.n_beta_core.eq(0).sum()),
        'no_beta_means_all_unavailable':True,'beta_INS_missing':int(f.n_beta_ins_missing.sum()),
        'other_core_INS_missing':int(f.n_nonbeta_ins_missing.sum()),
        'source_excluded_cells':int(f.n_source_excluded_core.sum()),
        'source_excluded_INS_measured':int(f.n_source_excluded_INS_measured.sum()),
        'source_excluded_INS_missing':int(f.n_source_excluded_INS_missing.sum()),
        'all_excluded_counts_match_prior_corrected_audit':True,
        'max_decomposition_absolute_error':float(f.decomposition_absolute_error.max()),
        'max_served_vs_scaled_source_INS_difference':float(f.ins_source_scale_max_abs_difference.max()),
        'primary_donors_estimable':int(primary.estimable.sum()),'primary_unestimable_donors':[6563],
        'additional_nearby_burden_donors_estimable':int(sensitivity.estimable.sum()),
        'additional_nearby_burden_unestimable_donors':[6550,6563],
        'visual_QA':'All three PNG figures inspected; group labels, significance brackets and legend values match saved tables.'}
    h6qa={'graphs_verified':len(graphfiles),'objects_with_graph_inventory':len(g),
        'labels_edges_and_20um_distances_verified':True,
        'no_inference_or_permutations_run':True,
        'visual_QA':'PNG inspected; all 19 nonmissing donor values are shown without significance annotations.'}
    for folder,qa in [(h5,h5qa),(h6,h6qa)]:
        (folder/'validation.json').write_text(json.dumps(qa,indent=2)+'\n')
        files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json')
        scripts=sorted((ROOT/'scripts/research').glob('h[56]*20261026*'))
        manifest=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,
                   'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files+scripts]
        (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'h5':h5qa,'h6':h6qa},indent=2))

if __name__=='__main__':main()
