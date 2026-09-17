#!/usr/bin/env python3
"""Complete estimable structural coupling, heterogeneity and prediction tests."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import balanced_accuracy_score
from common import MARKERS, SEED, sha256, write_json, offline_report_root
from run_hypotheses import continuous, donor_bootstrap, bh_ten, followups
from prediction import permutation_test

def dispersion(work):
    cols=MARKERS;d=work.copy()
    scales=d.groupby('case_id')[cols].quantile(.75)-d.groupby('case_id')[cols].quantile(.25)
    scale=scales.median().to_numpy();scale[scale<1e-10]=1
    Y=d[cols].to_numpy()/scale
    X=np.c_[np.ones(len(d)),np.log(d.area_um2)]
    weight=1/d.groupby('case_id').case_id.transform('size').to_numpy()
    beta=np.linalg.lstsq(X*np.sqrt(weight[:,None]),Y*np.sqrt(weight[:,None]),rcond=None)[0]
    residual=Y-X@beta;medians=[]
    for donor,indices in d.reset_index(drop=True).groupby('case_id').groups.items():
        values=residual[list(indices)];center=np.median(values,axis=0)
        value=np.sqrt(np.mean(np.sum((values-center)**2,axis=1)/len(cols)))
        row=d.loc[d.case_id==donor].iloc[0]
        medians.append(dict(case_id=donor,donor_status=row.donor_status,age=row.age,sex=row.sex,
                            response=np.log(max(value,1e-10)),n_islets=len(indices),
                            median_cell_area=float(d.loc[d.case_id==donor].median_cell_area.median()),
                            median_endocrine_cells=float(d.loc[d.case_id==donor].endocrine_measured_cells.median())))
    return pd.DataFrame(medians)

def donor_model(table,technical=False):
    X=pd.DataFrame({'Intercept':1.,'disease':(table.donor_status=='T1D').astype(float),
                    'age':table.age,'sex':table.sex.astype(str).str.lower().str.startswith('m').astype(float)})
    if technical:
        X['cell_area']=table.median_cell_area;X['log_endocrine_count']=np.log(table.median_endocrine_cells)
    if np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError('Rank-deficient donor dispersion model')
    return sm.OLS(table.response,X).fit()

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',required=True);p.add_argument('--bootstrap',type=int,default=200);p.add_argument('--permutations',type=int,default=200);a=p.parse_args()
    root=Path(a.input_root);out=offline_report_root(root,a.report_root);r=pd.read_parquet(root/'analysis/regions.parquet')
    donors=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str})
    results=pd.read_csv(out/'results.csv').set_index('id');rows=r.query('band_um==20 and geometry_valid and not overlap and annotation_review!="pending" and core_tissue_coverage>.99').copy()
    rows=rows.loc[rows.donor_status.isin(['ND','T1D'])]
    def save_result(h,d,values,unit):
        for k,v in values.items():
            if k in results:results.loc[h,k]=v
        results.loc[h,'status']='estimated_pending_scientific_review';results.loc[h,'reason']=''
        results.loc[h,'unit']=unit;results.loc[h,'n_donors']=d.case_id.nunique();results.loc[h,'n_islets']=d.islet_id.nunique()
        d.to_csv(out/f'{h}-endpoint.csv',index=False)
        rep=d.sort_values('islet_id').groupby('donor_status').head(1)
        rep[[c for c in ['islet_id','case_id','islet_key','donor_status'] if c in rep]].to_csv(out/f'{h}-representative-islets.csv',index=False)
        results['q_value']=bh_ten(results.p_value);results.reset_index().to_csv(out/'results.csv',index=False)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7,3))
        lo,hi=values.get('ci_low',np.nan),values.get('ci_high',np.nan);effect=values['effect']
        if np.isfinite(lo) and np.isfinite(hi):ax.plot([lo,hi],[0,0],color='#3366cc',linewidth=3)
        ax.scatter([effect],[0],s=50);ax.axvline(0,color='grey',linestyle='--')
        ax.set(title=h+' — donor-aware effect and bootstrap interval',xlabel=unit,yticks=[])
        fig.tight_layout()
        for suffix in ['png','pdf']:fig.savefig(out/f'{h}.{suffix}',dpi=150)
        plt.close(fig)
        print(h,values.get('effect'),values.get('p_value'),flush=True)
    coupling_path=root/'analysis/coupling.parquet'
    if not coupling_path.exists():raise ValueError('Spatial coupling artifact is required; no partial-profile fallback')
    if coupling_path.exists():
        coupling=pd.read_parquet(coupling_path)
        merged=rows.merge(coupling,on=['islet_id','case_id','band_um','compartment'],validate='one_to_one')
        merged['excess_area']=merged.excess_proximity*merged.nerve_area_um2
        merged['nerve_supported_area']=np.where(merged.excess_proximity.notna(),merged.nerve_area_um2,0.)
        agg=merged.groupby('islet_id').agg(excess_area=('excess_area','sum'),nerve_area=('nerve_supported_area','sum'),
                      capillary_area=('capillary_area_um2','sum'),area=('valid_area_um2','sum')).reset_index()
        agg['response']=agg.excess_area/agg.nerve_area;agg['nerve_density']=agg.nerve_area/agg.area;agg['capillary_density']=agg.capillary_area/agg.area
        d=rows.query('compartment=="core"').merge(agg,on='islet_id',validate='one_to_one').dropna(subset=['response'])
        d=d.merge(donors[['case_id','age','sex']],on='case_id',validate='many_to_one');d['log_area']=np.log(d.area_um2)
        value=continuous(d,['nerve_density','capillary_density'],'disease',a.bootstrap)
        save_result('H5',d,value,'difference in excess proximity fraction; analytic compartment-uniform null')
        results.loc['H5','reason']='Uniform-location null conditions on compartment and abundance; nerve autocorrelation is not preserved. Shape-preserving randomization is a sensitivity still awaiting validation.'
    profiles=pd.read_parquet(root/'analysis/endocrine_profiles.parquet')
    base=rows.query('compartment=="core"')
    d=base.merge(profiles,on=['islet_id','case_id'],validate='one_to_one').query('profile_eligible').merge(donors[['case_id','age','sex']],on='case_id',validate='many_to_one')
    d=d.loc[d.groupby('case_id').case_id.transform('size')>=3].reset_index(drop=True)
    table=dispersion(d);fit=donor_model(table)
    draws=donor_bootstrap(d,lambda sub:donor_model(dispersion(sub)).params['disease'],a.bootstrap)
    if len(draws)>=.8*a.bootstrap:
        ci=np.quantile(draws,[.025,.975]);save_result('H9',d,dict(effect=fit.params['disease'],ci_low=ci[0],ci_high=ci[1],p_value=fit.pvalues['disease'],bootstrap_success=len(draws)),'log donor dispersion across all 31 endocrine markers')
        table.to_csv(out/'H9-donor-dispersion.csv',index=False)
        write_json(out/'H9-technical-sensitivity.json',dict(effect_adjusted_for_cell_area_and_endocrine_count=float(donor_model(table,True).params['disease']),
                   acquisition='Donor and acquisition are inseparable; technical covariate sensitivity does not prove a biological difference.'))
    # H10 operates on donor summaries so no donor can leak across a fold.
    f=pd.read_csv(root/'trajectory/features.csv',dtype={'case_id':str})
    f=f.loc[f.eligible & f.donor_status.isin(['ND','T1D'])]
    f=f.merge(base[['islet_id']].drop_duplicates(),on='islet_id',validate='one_to_one')
    comps=[]
    for cid in donors.case_id:
        c=pd.read_csv(root/'regional'/f'islet_composition_{int(cid):04d}.tsv',sep='\t');c['islet_id']=cid+'_'+c.name
        denominator=c.cells_total-c.Triple_Neg
        for marker,col in [('INS','Ins_any'),('GCG','Glu_any'),('SST','Stt_any')]:c['fraction_'+marker]=c[col]/denominator.where(denominator>0)
        comps.append(c[['islet_id','fraction_INS','fraction_GCG','fraction_SST']])
    f=f.merge(pd.concat(comps),on='islet_id',validate='one_to_one')
    if coupling_path.exists():
        c=coupling.query('band_um==20').pivot(index='islet_id',columns='compartment',values='excess_proximity').add_prefix('coupling_')
        f=f.merge(c,on='islet_id',validate='one_to_one')
    spatial=[c for c in f if c.endswith('_area_fraction') or c.startswith('coupling_')]
    feature_cols=[f'{comp}__{marker}' for comp in ['core','peri'] for marker in MARKERS]
    basal=['area_um2','fraction_INS','fraction_GCG','fraction_SST']
    aggregate=f.groupby(['case_id','donor_status'])[basal+feature_cols+spatial].median().reset_index()
    aggregate=aggregate.merge(donors[['case_id','age','sex']],on='case_id',validate='one_to_one')
    aggregate['male']=aggregate.sex.str.lower().str.startswith('m').astype(float);aggregate['log_area']=np.log(aggregate.area_um2)
    basal=['age','male','log_area','fraction_INS','fraction_GCG','fraction_SST'];full=basal+feature_cols+spatial
    y=(aggregate.donor_status=='T1D').astype(int).to_numpy()
    aggregate.to_csv(out/'H10-donor-features.csv',index=False)
    v=permutation_test(aggregate[basal].to_numpy(),aggregate[full].to_numpy(),y,a.permutations,a.bootstrap)
    save_result('H10',f,v,'nested donor-held-out balanced accuracy improvement')
    prediction_table=aggregate[['case_id','donor_status']].copy();prediction_table['predicted_baseline']=v['predictions_baseline'];prediction_table['predicted_full']=v['predictions_full']
    prediction_table.to_csv(out/'H10-held-out.csv',index=False)
    pd.DataFrame({'permutation_delta':v['null']}).to_csv(out/'H10-permutations.csv',index=False)
    pd.DataFrame(v['coefficients'],columns=full).to_csv(out/'H10-fold-coefficients.csv',index=False)
    write_json(out/'H10-specification.json',dict(baseline_features=basal,full_features=full,unit='donor median of assessable islets',
                   baseline_balanced_accuracy=v['baseline_balanced_accuracy'],full_balanced_accuracy=v['full_balanced_accuracy'],
                   permutations=a.permutations,bootstrap_success=v['bootstrap_success'],
                   preprocessing='Within every inner and outer training fold: median imputation, standard scaling, L2 C selected from .01/.1/1 by inner donor-held-out balanced accuracy; ties prefer smaller C.',
                   missingness=aggregate[full].isna().sum().to_dict()))
    results['q_value']=bh_ten(results.p_value);results.reset_index().to_csv(out/'results.csv',index=False)
    # Retain completed F1. Other branches need actual prerequisite tests/stability.
    follow=followups(results.reset_index());old=pd.read_csv(out/'followups.csv');follow=pd.concat([old.loc[old.id=='F1'],follow.loc[follow.id!='F1']],ignore_index=True)
    coefficient=v['coefficients'];nonbaseline=np.array(full)[len(basal):]
    top=np.argsort(np.abs(coefficient[:,len(basal):]),axis=1)[:,-5:]
    frequency=np.bincount(top.ravel(),minlength=len(nonbaseline))/len(top)
    stability=pd.DataFrame({'feature':nonbaseline,'top_five_fold_fraction':frequency})
    stability.to_csv(out/'H10-feature-stability.csv',index=False)
    stable=bool((frequency>=.8).any() and v['effect']>0 and results.loc['H10','q_value']<.05)
    follow.loc[follow.id=='F4','branch']='INS_independent_stable_signature' if stable else 'prespecified_size_plus_structural_density_signature'
    follow.loc[follow.id=='F4','status']='unavailable: validated insulin-positive tissue-area burden is absent; no cell-count substitution'
    follow.to_csv(out/'followups.csv',index=False)
    write_json(out/'extended-run.json',dict(runner_sha256=sha256(__file__),cohort_sha256=sha256(root/'cohort.csv'),completed_tests=int(results.p_value.notna().sum()),
               scientific_acceptance=False,bootstrap_requested=a.bootstrap,pseudotime_used=False))

if __name__=='__main__':main()
