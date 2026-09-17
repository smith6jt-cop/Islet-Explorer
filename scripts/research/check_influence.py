#!/usr/bin/env python3
"""Donor influence and HPAP exclusion for the additional estimable hypotheses."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import balanced_accuracy_score
from run_hypotheses import design, within_between
from extended_tests import dispersion, donor_model
from prediction import nested_predict
from common import write_json, offline_report_root

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--report-root',required=True);a=p.parse_args();root=Path(a.input_root);out=offline_report_root(root,a.report_root);results=[]
    h5=pd.read_csv(out/'H5-endpoint.csv',dtype={'case_id':str});h9=pd.read_csv(out/'H9-endpoint.csv',dtype={'case_id':str})
    for h,d in [('H5',h5),('H9',h9)]:
        for donor in sorted(d.case_id.unique()):
            sub=d.loc[d.case_id!=donor].copy()
            try:
                if h=='H9':effect=donor_model(dispersion(sub)).params['disease']
                else:
                    X=design(sub)
                    for col in ['nerve_density','capillary_density']:
                        z=within_between(sub,col,col);X[col+'_within']=z[col+'_within'];X[col+'_between']=z[col+'_between']
                    effect=sm.WLS(sub.response,X,weights=1/sub.groupby('case_id').case_id.transform('size')).fit().params['disease']
                results.append(dict(id=h,excluded_donor=donor,exclude_hpap=donor=='78',effect=float(effect),n_donors=sub.case_id.nunique()))
            except (ValueError,np.linalg.LinAlgError) as e:results.append(dict(id=h,excluded_donor=donor,status=str(e)))
    spec=json.loads((out/'H10-specification.json').read_text());full=spec['full_features'];basal=spec['baseline_features']
    f=pd.read_csv(out/'H10-endpoint.csv',dtype={'case_id':str});meta=pd.read_csv(root/'donor_metadata.csv',dtype={'case_id':str})
    columns=['area_um2']+[c for c in full if c not in ['age','male','log_area']]
    aggregate=f.groupby(['case_id','donor_status'])[columns].median().reset_index().merge(meta[['case_id','age','sex']],on='case_id',validate='one_to_one')
    aggregate['male']=aggregate.sex.str.lower().str.startswith('m').astype(float);aggregate['log_area']=np.log(aggregate.area_um2)
    aggregate.to_csv(out/'H10-donor-features.csv',index=False)
    # Complete nested reconstruction on each remaining donor set; preprocessing
    # and C selection are repeated rather than reusing the original prediction.
    for donor in sorted(aggregate.case_id):
        sub=aggregate.loc[aggregate.case_id!=donor];y=(sub.donor_status=='T1D').astype(int).to_numpy()
        x,_,_=nested_predict(sub[basal].to_numpy(),y);z,_,_=nested_predict(sub[full].to_numpy(),y)
        results.append(dict(id='H10',excluded_donor=donor,exclude_hpap=donor=='78',effect=float(balanced_accuracy_score(y,z)-balanced_accuracy_score(y,x)),n_donors=len(sub)))
    pd.DataFrame(results).to_csv(out/'additional-donor-influence.csv',index=False)
    write_json(out/'acquisition_limitations.json',dict(acquisitions_per_donor=1,separable=False,
                  consequences='Acquisition effects cannot be distinguished from donor effects. Cohort-level stability and donor-held-out prediction do not establish platform-independent biological validity.'))

if __name__=='__main__':main()
