#!/usr/bin/env python3
"""Approved H6 fixed-geometry/count-preserving label permutation measurements."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h5_h10_original_plan_20261026/h6'
SEED=20261026
PERMUTATIONS=1999

def evaluate_graph(args):
    name,case,labels,edges=args
    stable=int.from_bytes(hashlib.sha256(name.encode()).digest()[:8],'little')
    rng=np.random.default_rng(np.random.SeedSequence([SEED,stable]))
    def fractions(v):
        aa,bb=v[edges[:,0]],v[edges[:,1]]
        ba=np.sum(((aa==0)&(bb==1))|((aa==1)&(bb==0)))
        bd=np.sum(((aa==0)&(bb==2))|((aa==2)&(bb==0)))
        return np.array([ba+bd,ba,bd])/len(edges)
    observed=fractions(labels)
    null=np.empty((PERMUTATIONS,3))
    for i in range(PERMUTATIONS):null[i]=fractions(rng.permutation(labels))
    n=len(labels);nb=(labels==0).sum();na=(labels==1).sum();nd=(labels==2).sum()
    exact=2*nb*np.array([na+nd,na,nd])/(n*(n-1))
    row={'islet_id':name,'case_id':case,'n_permutations':PERMUTATIONS}
    for j,kind in enumerate(['combined','alpha','delta']):
        mu=null[:,j].mean();sd=null[:,j].std(ddof=1)
        row.update({f'{kind}_observed':observed[j],f'{kind}_null_mean':mu,
                    f'{kind}_null_sd':sd,f'{kind}_null_lower':np.quantile(null[:,j],.025),
                    f'{kind}_null_upper':np.quantile(null[:,j],.975),
                    f'{kind}_excess':observed[j]-mu,f'{kind}_null_exact_mean':exact[j],
                    f'{kind}_null_mc_error':mu-exact[j],f'{kind}_spatially_informative':sd>1e-12})
    return row,null.astype(np.float32)

def main():
    spec={'approval':'Original spatial label-permutation measurement explicitly approved by user via root',
          'seed':SEED,'permutations_per_object':PERMUTATIONS,
          'random_stream':'SHA256(object identifier) plus master seed through NumPy SeedSequence',
          'primary_contact':'Observed combined Beta–Alpha/Beta–Delta edge fraction minus permutation mean',
          'primary_graph':'Existing 20 micrometre endocrine-centroid graph',
          'preserved':'Every cell location, graph edge, endocrine cell count and cell-type count within each object',
          'primary_eligibility':'Existing >=10 endocrine cells, >=2 Beta, >=2 Alpha+Delta and >0 graph edges, plus nonzero combined null variance and finite within-Beta INS',
          'primary_adjustment':['Beta fraction','Delta fraction','other-endocrine fraction (Alpha reference)',
                                'log core area','log1p nearby immune count','log1p core immune count'],
          'standardization':'Continuous predictors within donor; outcome globally over eligible objects',
          'models':'Donor random intercept and contact slope, disease main effects, common contact slope and contact-by-disease model. Retain independent random-slope covariance sensitivity if correlated fit singular.',
          'donor_replication_check':'Separate full-rank donor regressions with positive residual df; no nonconstant covariates dropped to recover sparse donors. One-sample signed-rank versus zero and disease Kruskal–Wallis plus three BH-adjusted rank-sum pairs.',
          'biological_test':'Adjusted contact main association and its direction, separately from disease heterogeneity tests',
          'sensitivities':['Separate Alpha and Delta contact excess','Beta density outcome','Contact by nearby immune count interaction','Raw and RedSEA-corrected INS','Independent nerve/vessel density covariates from H7','True cell-boundary contact graphs if available'],
          'composition_size_comparison':'Within-donor matching strata formed by pooled tertiles of Beta fraction, Delta fraction and log core area; compare positive versus nonpositive combined excess where both occur. Descriptive stratum differences retained; no claim of exact matching or immune adjustment.',
          'meaningful_effect_exclusion':'No equivalence bound was prespecified; nonsignificance cannot exclude a meaningful effect.'}
    (OUT/'analysis_specification_before_fit.json').write_text(json.dumps(spec,indent=2)+'\n')
    inventory=pd.read_csv(OUT/'graph_inventory.csv').set_index('islet_id')
    tasks=[]
    for path in sorted((OUT/'graphs_by_donor').rglob('*.npz')):
        g=np.load(path);tasks.append((path.stem,int(path.parent.name),g['labels'],g['edges']))
    rows=[];nulls=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for i,(row,null) in enumerate(pool.map(evaluate_graph,tasks)):
            rows.append(row);nulls.append(null)
            if (i+1)%100==0:print(f'Permuted {i+1}/{len(tasks)} objects',flush=True)
    measurements=pd.DataFrame(rows)
    measurements.to_csv(OUT/'permutation_measurements.csv',index=False)
    np.savez_compressed(OUT/'permutation_null_draws.npz',islet_ids=measurements.islet_id.to_numpy(dtype=str),
                        statistic_names=np.array(['combined','alpha','delta']),draws=np.stack(nulls))
    features=pd.read_parquet(OUT.parent/'h5/object_ins_features.parquet')
    d=measurements.merge(inventory.reset_index(),on=['islet_id','case_id'],validate='one_to_one')
    d=d.merge(features.drop(columns=['donor_status']),on=['islet_id','case_id'],validate='one_to_one')
    d['beta_fraction']=d.n_beta_nodes/d.n_endocrine_nodes
    d['delta_fraction']=d.n_delta_nodes/d.n_endocrine_nodes
    d['other_endocrine_fraction']=(d.n_endocrine_nodes-d.n_beta_nodes-d.n_alpha_nodes-d.n_delta_nodes)/d.n_endocrine_nodes
    d['log_core_area']=np.log(d.core_area_um2)
    d['log_nearby_immune']=np.log1p(d.n_immune_nearby)
    d['log_core_immune']=np.log1p(d.n_immune_core)
    d['h6_primary_eligible']=d.combined_spatially_informative&np.isfinite(d.ins_mean_beta_core)
    d.to_parquet(OUT/'h6_object_features.parquet',index=False)
    d.to_csv(OUT/'h6_object_features.csv',index=False)
    print(d.groupby('donor_status').agg(graphs=('islet_id','size'),primary=('h6_primary_eligible','sum')).to_string(),flush=True)

if __name__=='__main__':main()
