"""Numerical checks for follow-up estimands and the conditional spatial null."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from shapely.geometry import Point, box
import shapely

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
import h1_followups_20261024 as f
from h234_existing_data_20261024 import rigid_rotations


def synthetic(n=1200):
    rng=np.random.default_rng(317)
    d=pd.DataFrame(dict(case_id=np.repeat([str(i) for i in range(6)],n//6),
        object_class=np.where(rng.random(n)>.5,'Islet','EndocrineCluster'),
        core_area_um2=np.exp(rng.normal(8,1,n)),peri_area_um2=np.exp(rng.normal(8.4,.6,n)),
        field_density=np.exp(rng.normal(4,.7,n)),n_immune_nearby=rng.integers(2,20,n),
        immune_localization_fraction=rng.uniform(0,1,n),beta_density_mm2=rng.uniform(100,10000,n)))
    d['n_beta_core']=np.maximum(1,np.round(d.beta_density_mm2*d.core_area_um2/1e6))
    return d


def test_original_poisson_point_estimate_is_preserved():
    d=synthetic();d['has_core_immune']=1
    b,_=f.fit(d,'count')
    rows,_,_=f.h1.fit_model(d,'comparison')
    assert abs(b-rows[0]['log_rr'])<1e-7


def test_reciprocal_area_offset_recovers_injected_effect():
    d=synthetic()
    x=np.log2(1+d.beta_density_mm2/1000)
    eta=np.log(d.core_area_um2/d.peri_area_um2)+.35*x-.7
    probability=1/(1+np.exp(-eta))
    # Sample actual counts to avoid a degenerate perfect-prediction dataset.
    d['n_immune_core']=np.random.default_rng(43).binomial(d.n_immune_nearby,probability)
    b,details=f.fit(d,'reciprocal')
    assert abs(b-.35)<.08
    assert 'allocation_10pp' not in details['columns']


def test_positive_density_recovers_injected_mean():
    d=synthetic()
    mean=np.exp(7+.12*d.immune_localization_fraction/.1)
    d['beta_density_mm2']=np.random.default_rng(17).gamma(20,mean/20)
    b,_=f.fit(d,'positive_density')
    assert abs(b-.12)<.01


def test_presence_reports_constant_donors():
    d=synthetic()
    d.loc[d.case_id.eq('0'),'n_beta_core']=0
    # Explicit construction avoids treating discarded complete-separation
    # groups as if their donor intercept were estimated reliably.
    d.loc[d.case_id.eq('1'),'n_beta_core']=1
    for g in ['2','3','4','5']:
        ids=d.index[d.case_id.eq(g)];d.loc[ids[::2],'n_beta_core']=0;d.loc[ids[1::2],'n_beta_core']=1
    cleaned,dropped=f.prepare(d,'presence')
    assert set(dropped)=={'0','1'}
    assert cleaned.case_id.nunique()==4


def test_rigid_null_preserves_pattern_and_support():
    core=box(-2,-1,2,1);tissue=box(-3,-.8,3,.9)
    xy=np.array([[-.4,0],[.5,.2],[.6,-.2]])
    rot,tried,accepted=rigid_rotations(xy,np.zeros(2),core,tissue,np.random.default_rng(5),target=100)
    assert rot.shape==(100,3,2)
    assert np.all(shapely.covers(core,shapely.points(rot)))
    assert np.all(shapely.covers(tissue,shapely.points(rot)))
    for sample in rot:
        assert np.allclose(pdist(sample),pdist(xy))
        assert np.allclose(np.linalg.norm(sample,axis=1),np.linalg.norm(xy,axis=1))


def test_null_unavailable_is_not_observed_zero_effect():
    xy=np.array([[100.,0]])
    rot,_,_=rigid_rotations(xy,np.zeros(2),Point(0,0).buffer(1),box(-2,-2,2,2),
                            np.random.default_rng(5),target=10,maximum=100)
    assert len(rot)==0
