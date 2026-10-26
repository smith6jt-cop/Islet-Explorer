#!/usr/bin/env python3
"""Post-fit diagnostic separating inference changes from common-adjustment assumptions."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from h1_followups_20261024 import OUT, SOURCE, h1
from h1_heterogeneity_20261022 import fit_group_slopes


def slopes(d):
    rows,_,_=fit_group_slopes(d)
    return rows['log_rr_10pp'].to_numpy(float)


def run():
    d=h1.eligible(pd.read_parquet(SOURCE));b=slopes(d)
    parts=[];deletes=[]
    for group in h1.STATUS:
        gs=sorted(d.loc[d.donor_status.eq(group),'case_id'].unique());values=[]
        for g in gs:
            v=slopes(d.loc[d.case_id.ne(g)])
            values.append(v);deletes.append(dict(deleted_donor=g,deleted_group=group,
                **dict(zip(h1.STATUS,v))))
        a=np.array(values);a-=a.mean(axis=0)
        parts.append((len(gs)-1)/len(gs)*(a.T@a))
    cov=sum(parts)
    # Conservative df=5 rather than borrowing all 20 donors for a subgroup.
    df=5;q=stats.t.ppf(.975,df);rows=[]
    for j,group in enumerate(h1.STATUS):
        se=np.sqrt(cov[j,j]);rows.append(dict(group=group,effect=np.exp(b[j]),
            lower=np.exp(b[j]-q*se),upper=np.exp(b[j]+q*se),p=2*stats.t.sf(abs(b[j]/se),df),df=df))
    pd.DataFrame(rows).to_csv(OUT/'h1_common_adjustment_sensitivity.csv',index=False)
    contrasts=[]
    for name,c in [('T1D vs ND',[-1,0,1]),('T1D vs Aab+',[0,-1,1]),('Aab+ vs ND',[-1,1,0])]:
        c=np.array(c);delta=c@b;se=np.sqrt(c@cov@c)
        contrasts.append(dict(contrast=name,ratio_of_RRs=np.exp(delta),
            lower=np.exp(delta-q*se),upper=np.exp(delta+q*se),p=2*stats.t.sf(abs(delta/se),df),df=df))
    c=pd.DataFrame(contrasts);c['p_holm']=multipletests(c.p,method='holm')[1]
    c.to_csv(OUT/'h1_common_adjustment_contrasts.csv',index=False)
    pd.DataFrame(deletes).to_csv(OUT/'h1_common_adjustment_deletions.csv',index=False)
    (OUT/'h1_common_adjustment_covariance.json').write_text(json.dumps(dict(
        stratum_covariances=[a.tolist() for a in parts],covariance=cov.tolist(),
        group_order=h1.STATUS,df=df),indent=2))
    print(pd.DataFrame(rows).to_string(index=False));print(c.to_string(index=False))


if __name__=='__main__':run()
