#!/usr/bin/env python3
"""Validate numerical output, complete donor refits, and unchanged inputs."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from h1_followups_20261024 import OUT,SOURCE,h1
from h234_existing_data_20261024 import CORRECTED,sha


def run():
    results=pd.read_csv(OUT/'h1_followup_models.csv')
    deletes=pd.read_csv(OUT/'h1_donor_deletion_fits.csv')
    assert len(results)==20 and results.status.eq('ok').all()
    assert len(deletes)==200 and deletes.status.eq('ok').all()
    for r in results.itertuples():
        j=deletes.loc[deletes.group.eq(r.group)&deletes.model.eq(r.model)&deletes.variant.eq(r.variant)]
        assert len(j)==r.donors
        se=np.sqrt((r.donors-1)/r.donors*np.sum((j.log_effect-j.log_effect.mean())**2))
        assert np.isclose(se,r.se_jackknife,rtol=1e-10)
        assert np.isclose(r.lower,np.exp(r.log_effect-stats.t.ppf(.975,r.df)*se))
        assert np.isclose(r.upper,np.exp(r.log_effect+stats.t.ppf(.975,r.df)*se))
    primary=results.loc[results.group.eq('Pooled')&results.model.eq('count')].iloc[0]
    assert np.isclose(primary.effect,.9486867298,atol=1e-9)
    raw=pd.read_parquet(SOURCE);d=h1.eligible(raw)
    assert len(d)==2937 and d.case_id.nunique()==20 and '78' not in set(d.case_id)
    assert d.islet_id.is_unique
    t=d.loc[d.donor_status.eq('T1D')]
    assert len(t)==1088 and t.n_beta_core.gt(0).sum()==80
    common=pd.read_csv(OUT/'h1_common_adjustment_deletions.csv')
    assert len(common)==20 and np.isfinite(common[h1.STATUS].to_numpy()).all()
    spatial=pd.read_csv(OUT/'h3_object_spatial_null.csv')
    assert spatial.islet_id.nunique()==89 and len(spatial)==267
    assert spatial.null_draws.eq(999).all()
    assert np.allclose(spatial.excess,spatial.observed-spatial.expected)
    assert spatial[['observed','expected']].ge(0).all().all()
    assert spatial[['observed','expected']].le(1).all().all()
    h4=pd.read_csv(OUT/'h4_lgals3_measurement_audit.csv')
    assert h4.query('population == "Macrophage"').missing.sum()==0
    assert h4.islet_id.nunique()==2937
    readiness=pd.read_csv(OUT/'h4_donor_readiness.csv')
    assert readiness.unique_nearby_macrophages.sum()==4047
    evidence=pd.read_csv(OUT/'structural_evidence_audit.csv')
    assert evidence.loc[evidence.structure.eq('ColIV'),'decision'].eq('rejected').all()
    proof=pd.read_csv(OUT/'h3_current_mask_geometry_verification.csv')
    assert len(proof)==5 and proof.geometry_sets_identical.all()
    sources={str(SOURCE):sha(SOURCE)}
    for path in [OUT/'h234_source_hashes.json',OUT/'h4_source_hashes.json']:
        sources.update(json.loads(path.read_text()))
    for r in proof.itertuples():sources[r.current_path]=r.current_sha256
    for p,h in sources.items():
        assert sha(p)==h,f'Input changed: {p}'
    previous=json.loads((CORRECTED/'delivery_manifest.json').read_text())
    for p,h in previous['files'].items():
        assert sha(CORRECTED/p)==h,f'Previous corrected output changed: {p}'
    for p,h in previous['scripts'].items():
        assert sha(p)==h,f'Previous corrected analysis script changed: {p}'
    summary=dict(status='passed',h1_main_fits=20,h1_donor_deletion_fits=200,
        common_adjustment_deletion_fits=20,previous_output_files_unchanged=len(previous['files']),
        previous_scripts_unchanged=len(previous['scripts']),current_inputs_unchanged=len(sources),
        h3_objects_with_999_valid_rotations=89,current_mask_geometry_identity_verified=5,
        h4_unique_macrophages=4047,excluded_donor_78_absent=True,
        source_hashes=sources)
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2))
    scripts=[Path('scripts/research')/n for n in ['h1_followups_20261024.py',
        'h1_common_adjustment_sensitivity_20261024.py','h234_existing_data_20261024.py',
        'h4_measurement_audit_20261024.py','h3_followup_figures_20261024.py',
        'validate_h1_followups_20261024.py']]
    scripts.append(Path('scripts/tests/test_h1_followups_20261024.py'))
    manifest=dict(scripts={str(p.resolve()):sha(p) for p in scripts},
        files={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*')
               if p.is_file() and p.name!='delivery_manifest.json'})
    (OUT/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':run()
