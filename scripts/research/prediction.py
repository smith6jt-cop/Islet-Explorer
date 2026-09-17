"""Nested donor-held-out classification and donor-label permutation.

One feature row per donor avoids leakage and donor-size weighting. All imputation,
scaling and regularization selection occur strictly inside the training folds.
"""
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from common import SEED

def estimator(c):
    return make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True),StandardScaler(),
                         LogisticRegression(C=c,penalty='l2',solver='liblinear',class_weight='balanced',random_state=SEED,max_iter=1000))

def nested_predict(X,y,groups=None):
    n=len(y);pred=np.zeros(n,int);coefficients=[];choices=[]
    groups=np.arange(n) if groups is None else np.asarray(groups)
    for outer in np.unique(groups):
        train=np.flatnonzero(groups!=outer);test=np.flatnonzero(groups==outer)
        scores=[]
        for c in [.01,.1,1.]:
            inner_pred=[];truth=[]
            for held in np.unique(groups[train]):
                inner=train[groups[train]!=held];validation=train[groups[train]==held]
                if len(np.unique(y[inner]))<2:raise ValueError('Insufficient donors in inner training fold')
                model=estimator(c).fit(X[inner],y[inner]);inner_pred.extend(model.predict(X[validation]));truth.extend(y[validation])
            scores.append(balanced_accuracy_score(truth,inner_pred))
        c=[.01,.1,1.][int(np.argmax(scores))]
        model=estimator(c).fit(X[train],y[train]);pred[test]=model.predict(X[test])
        coefficients.append(model[-1].coef_[0]);choices.append(c)
    return pred,np.asarray(coefficients),choices

def _permuted(baseline,full,perm):
    x,_,_=nested_predict(baseline,perm);z,_,_=nested_predict(full,perm)
    return balanced_accuracy_score(perm,z)-balanced_accuracy_score(perm,x)

def _bootstrap(baseline,full,y,draw):
    try:
        x,_,_=nested_predict(baseline[draw],y[draw],draw);z,_,_=nested_predict(full[draw],y[draw],draw)
        return balanced_accuracy_score(y[draw],z)-balanced_accuracy_score(y[draw],x)
    except ValueError:return np.nan

def permutation_test(baseline,full,y,permutations=200,bootstrap=200):
    a,_,_=nested_predict(baseline,y);b,coef,choices=nested_predict(full,y)
    ba=balanced_accuracy_score(y,a);bb=balanced_accuracy_score(y,b);delta=bb-ba
    rng=np.random.default_rng(SEED)
    from joblib import Parallel, delayed
    null=Parallel(n_jobs=4)(delayed(_permuted)(baseline,full,rng.permutation(y)) for _ in range(permutations))
    print('H10 full-procedure permutations complete',flush=True)
    draws=Parallel(n_jobs=4)(delayed(_bootstrap)(baseline,full,y,rng.choice(len(y),len(y),replace=True)) for _ in range(bootstrap))
    draws=np.asarray(draws);draws=draws[np.isfinite(draws)]
    ci=np.quantile(draws,[.025,.975]) if len(draws)>=.8*bootstrap else [np.nan,np.nan]
    return dict(effect=delta,p_value=(1+sum(v>=delta-1e-12 for v in null))/(1+permutations),
                baseline_balanced_accuracy=ba,full_balanced_accuracy=bb,null=null,predictions_baseline=a,
                predictions_full=b,coefficients=coef,regularization=choices,
                ci_low=ci[0],ci_high=ci[1],bootstrap_success=len(draws),bootstrap_draws=draws)
