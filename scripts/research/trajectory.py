#!/usr/bin/env python3
"""Outcome-blind landmark diffusion ordering, with fixed Nyström projection.

Ten nontrivial diffusion eigenvectors; multiscale distance weights lambda/(1-
lambda). Root is a hash-selected endpoint, never a biological direction.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, kruskal
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
from common import MARKERS, SEED, sha256, write_json
from phenotypes import dependencies

MODULES = [set('INS GCG SST CHGA'.split()), set('B3TUBB PGP9.5 GAP43'.split()),
           set('CD31 CD34 SMA ColIV'.split()),set('PDPN'.split()),
           set('CD45 CD3e CD8a CD4 CD20 CD35 CD68 CD163 HLADR'.split()),
           set('CK19 panCK ECAD BCatenin KRT14 VIM LGALS3 CD44'.split()),set('BActin NaKATPase'.split())]

def excluded_for_outcome(outcome, rules=None):
    excluded=set('INS GCG SST CHGA'.split())
    if outcome in MARKERS: excluded.add(outcome)
    elif rules is not None:
        excluded.update(set().union(*(set(r) for r in dependencies(outcome,rules))))
    else: raise ValueError('Unknown outcome requires explicit phenotype rules')
    # Expand entire marker modules including every inherited gating dependency.
    for module in MODULES:
        if module & excluded: excluded.update(module)
    return excluded

def features(excluded=()):
    omit=set('INS GCG SST CHGA'.split())|set(excluded)
    return [f'{comp}__{m}' for comp in ['core','peri'] for m in MARKERS if m not in omit]

def landmarks(frame, cap=200, seed=SEED):
    n=min(cap,int(frame.groupby('case_id').size().min()))
    if n<1: raise ValueError('No eligible landmarks in a donor')
    chosen=[]
    for donor,g in frame.groupby('case_id',sort=True):
        g=g.sort_values('islet_id')
        # Stable ranks break area ties by ID, and round-robin draws keep small
        # strata represented without using disease labels or marker outcomes.
        quartile=pd.Series(0,index=g.index) if len(g)==1 else pd.qcut(g.area_um2.rank(method='first'),4,labels=False,duplicates='drop')
        buckets=[]
        for q in sorted(quartile.unique()):
            ids=list(g.loc[quartile==q].index)
            buckets.append(sorted(ids,key=lambda i:hashlib.sha256(f'{seed}|{frame.at[i,"islet_id"]}'.encode()).digest()))
        draw=[]
        while len(draw)<n:
            for bucket in buckets:
                if bucket and len(draw)<n: draw.append(bucket.pop(0))
        chosen.extend(draw)
    return np.array(sorted(chosen),int)

def _spectral(graph,n_components=10):
    # Anisotropic alpha=1 density correction followed by symmetric Markov
    # conjugation. The saved psi are right eigenfunctions of the Markov chain.
    q=np.asarray(graph.sum(axis=1)).ravel()
    kernel=sparse.diags(1/q)@graph@sparse.diags(1/q)
    d=np.asarray(kernel.sum(axis=1)).ravel()
    sym=sparse.diags(1/np.sqrt(d))@kernel@sparse.diags(1/np.sqrt(d))
    count=min(n_components+1,graph.shape[0]-1)
    if count<3: raise ValueError('Diffusion component too small')
    vals,u=eigsh(sym,k=count,which='LA',v0=np.linspace(1,2,len(d)),tol=1e-10)
    order=np.argsort(vals)[::-1]; vals=vals[order][1:];u=u[:,order][:,1:]
    nonzero=np.abs(vals)>1e-10;vals=vals[nonzero];u=u[:,nonzero]
    if len(vals)<2:raise ValueError('Insufficient nonzero diffusion components for stable Nystrom extension')
    psi=u/np.sqrt(d[:,None])
    # Sign canonicalization is numeric only, independent of any clinical label.
    signs=np.sign(psi[np.abs(psi).argmax(axis=0),np.arange(psi.shape[1])]);signs[signs==0]=1
    return vals,psi*signs,q,d

def fit(frame,k=15,cap=200,seed=SEED,excluded=(),feature_columns=None,legacy_intensity=False):
    if legacy_intensity:
        cols=features(excluded)
    else:
        cols=list(feature_columns) if feature_columns is not None else sorted(c for c in frame if c.startswith(('core__count__','peri__count__')))
        if any(not c.startswith(('core__count__','peri__count__')) for c in cols) or len(set(cols))!=len(cols):
            raise ValueError('Feature manifest must contain unique post-typing composition columns only')
        if not cols:
            raise ValueError('Post-typing composition features required; archived intensity fits require explicit legacy_intensity=True')
        if excluded:
            raise ValueError('Outcome exclusion requires rebuilding typing and the endocrine denominator, not column removal')
    missing=set(cols)-set(frame)
    if missing: raise ValueError(f'Missing required trajectory features: {sorted(missing)}')
    if frame.islet_id.duplicated().any(): raise ValueError('Duplicate trajectory IDs')
    f=frame.loc[frame.eligible & np.isfinite(frame[cols]).all(axis=1)].sort_values('islet_id').reset_index(drop=True)
    if f.empty: raise ValueError('No eligible islets')
    li=landmarks(f,cap,seed); X=f[cols].to_numpy(float); L=X[li]
    center=np.median(L,axis=0); scale=np.quantile(L,.75,axis=0)-np.quantile(L,.25,axis=0)
    keep=np.ptp(L,axis=0)>1e-12
    scale[scale<1e-12]=1
    cols_used=np.array(cols)[keep].tolist()
    if keep.sum()<3: raise ValueError('Fewer than three variable features')
    Z=(X[:,keep]-center[keep])/scale[keep]; L=Z[li]
    if len(L)<=k: raise ValueError('Not enough landmarks for fixed k')
    nbr=NearestNeighbors(n_neighbors=k+1,metric='cosine',algorithm='brute').fit(L)
    dist,inds=nbr.kneighbors(L)
    # Remove self by identity even when duplicate feature vectors reorder ties.
    pairs=[[(j,d) for j,d in zip(js,ds) if j!=i][:k] for i,(js,ds) in enumerate(zip(inds,dist))]
    inds=np.array([[j for j,d in row] for row in pairs]);dist=np.array([[d for j,d in row] for row in pairs])
    sigma=np.maximum(dist[:,-1],1e-10)
    weights=np.exp(-dist**2/(sigma[:,None]*sigma[inds]))
    graph=sparse.csr_matrix((weights.ravel(),(np.repeat(np.arange(len(L)),k),inds.ravel())),shape=(len(L),len(L)))
    graph=graph.maximum(graph.T);graph.setdiag(1)
    n_comp,comp=connected_components(graph,directed=False)
    primary=int(np.argmax(np.bincount(comp)));ix=np.flatnonzero(comp==primary)
    vals,psi,q,d=_spectral(graph[ix][:,ix])
    metric=psi*(vals/np.maximum(1-vals,1e-8))
    medoid=int(np.argmin(cdist(metric,metric).sum(axis=1)))
    endpoint_a=int(np.argmax(np.linalg.norm(metric-metric[medoid],axis=1)))
    endpoint_b=int(np.argmax(np.linalg.norm(metric-metric[endpoint_a],axis=1)))
    ids=f.iloc[li].islet_id.to_numpy(dtype=str)
    endpoints=[endpoint_a,endpoint_b]
    root=min(endpoints,key=lambda j:hashlib.sha256(ids[ix[j]].encode()).digest())
    raw=np.linalg.norm(metric-metric[root],axis=1)
    norm=raw.max() or 1.;raw/=norm
    model=dict(columns=cols_used,center=center[keep],scale=scale[keep],landmark_vectors=L,
               landmark_ids=ids,landmark_components=comp,primary=primary,primary_indices=ix,
               eigenvalues=vals,eigenfunctions=psi,q=q,degree=d,sigma=sigma,
               root_metric=metric[root],normalizer=norm,k=k,graph=graph,root_id=ids[ix[root]],
               endpoint_ids=[ids[ix[j]] for j in endpoints],medoid_id=ids[ix[medoid]],
               dropped=[c for c,v in zip(cols,keep) if not v],seed=seed,cap=cap,
               representation='legacy_intensity' if legacy_intensity else 'typed_composition')
    coords=project(frame,model)
    # Nyström is defined off-sample. Preserve the exact fitted landmark rows.
    positions=coords.set_index('islet_id').index.get_indexer(ids[ix])
    coords.loc[positions,'pseudotime_raw']=raw
    coords.loc[positions,'supported']=True
    coords.loc[positions,'component']=primary
    for j in range(psi.shape[1]): coords.loc[positions,f'DC{j+1}']=psi[:,j]*vals[j]
    coords['landmark']=coords.islet_id.isin(ids)
    model['primary_landmark_raw']=raw
    model['primary_landmark_coordinates']=psi*vals
    coords['pseudotime_percentile']=np.nan
    good=np.isfinite(coords.pseudotime_raw)
    model['reference_raw']=np.sort(coords.loc[good,'pseudotime_raw'].to_numpy())
    coords.loc[good,'pseudotime_percentile']=reference_percentile(coords.loc[good,'pseudotime_raw'],model['reference_raw'])
    return coords,model

def reference_percentile(values,reference):
    """Frozen mid-distribution ECDF, including identical treatment of ties."""
    ref=np.asarray(reference,float);x=np.asarray(values,float)
    if not len(ref):return np.full(x.shape,np.nan)
    return (np.searchsorted(ref,x,side='left')+np.searchsorted(ref,x,side='right'))/(2*len(ref))

def project(frame,model):
    cols=model['columns'];out=frame[[c for c in ['islet_id','case_id','donor_status','donor_source','islet_key','area_um2','measured_cells','raw_cells','eligible','trajectory_feature_supported','exclusion_reason'] if c in frame]].reset_index(drop=True).copy()
    out['component']=-1;out['supported']=False;out['support_distance']=np.nan;out['pseudotime_raw']=np.nan
    for j in range(10): out[f'DC{j+1}']=np.nan
    eligible=frame.eligible.to_numpy(bool)&np.isfinite(frame[cols]).all(axis=1).to_numpy()
    pos=np.flatnonzero(eligible)
    if not len(pos): return out
    Z=(frame.loc[eligible,cols].to_numpy()-model['center'])/model['scale']
    ds,ns=NearestNeighbors(n_neighbors=model['k'],metric='cosine',algorithm='brute').fit(model['landmark_vectors']).kneighbors(Z)
    comps=model['landmark_components'][ns];assignment=comps[:,0]
    # Support requires all neighbors in one component and distance within twice
    # the nearest landmark's local k-neighbor radius. Frozen, outcome-blind rule.
    supported=(comps==assignment[:,None]).all(axis=1)&(ds[:,0]<=2*model['sigma'][ns[:,0]])
    out.loc[pos,'component']=assignment;out.loc[pos,'support_distance']=ds[:,0]
    out.loc[pos,'supported']=supported
    good=supported&(assignment==model['primary']);pos=pos[good];ns=ns[good];ds=ds[good]
    back=np.full(len(model['landmark_vectors']),-1,int);back[model['primary_indices']]=np.arange(len(model['primary_indices']))
    local=back[ns]
    sig=np.maximum(ds[:,-1],1e-10)
    w=np.exp(-ds**2/(sig[:,None]*model['sigma'][ns]))/model['q'][local]
    w/=w.sum(axis=1,keepdims=True)
    vals=model['eigenvalues'];psi=np.sum(w[:,:,None]*model['eigenfunctions'][local],axis=1)/vals
    metric=psi*(vals/np.maximum(1-vals,1e-8))
    out.loc[pos,'pseudotime_raw']=np.linalg.norm(metric-model['root_metric'],axis=1)/model['normalizer']
    for j in range(len(vals)): out.loc[pos,f'DC{j+1}']=psi[:,j]*vals[j]
    # Re-projecting an unchanged fitted landmark must preserve its exact
    # eigenfunction coordinate, not introduce a self-neighbor approximation.
    if 'primary_landmark_raw' in model:
        primary_ids=model['landmark_ids'][model['primary_indices']]
        lookup={str(v):i for i,v in enumerate(primary_ids)}
        for row,iid in enumerate(out.islet_id.astype(str)):
            j=lookup.get(iid)
            if j is None or not eligible[row]:continue
            vector=(frame.iloc[row][cols].to_numpy(float)-model['center'])/model['scale']
            if not np.array_equal(vector,model['landmark_vectors'][model['primary_indices'][j]]):continue
            out.loc[row,'pseudotime_raw']=model['primary_landmark_raw'][j]
            out.loc[row,'supported']=True;out.loc[row,'component']=model['primary']
            for dc in range(len(vals)):out.loc[row,f'DC{dc+1}']=model['primary_landmark_coordinates'][j,dc]
    if 'reference_raw' in model:
        ref=model['reference_raw'];out['pseudotime_percentile']=np.nan
        finite=out.pseudotime_raw.notna()
        out.loc[finite,'pseudotime_percentile']=reference_percentile(out.loc[finite,'pseudotime_raw'],ref)
    return out

def graph_layout(frame,coordinates,model):
    """Independent spectral layout per component, for graph orientation only.

    Component panels are packed on a grid; their axes are not diffusion
    coordinates or a comparable pseudotime. Eligible non-landmarks receive a
    fixed neighbor-weighted display position in their assigned component.
    """
    from scipy.sparse.csgraph import laplacian
    comp=model['landmark_components'];components=sorted(np.unique(comp));xy=np.full((len(comp),2),np.nan)
    width=max(1,int(np.ceil(np.sqrt(len(components)))))
    for slot,component in enumerate(components):
        ix=np.flatnonzero(comp==component);g=model['graph'][ix][:,ix]
        if len(ix)<4:
            local=np.c_[np.linspace(-.5,.5,len(ix)),np.zeros(len(ix))]
        else:
            lap=laplacian(g,normed=True)
            values,vectors=eigsh(lap,k=3,which='SM',v0=np.linspace(1,2,len(ix)),tol=1e-8)
            local=vectors[:,np.argsort(values)[1:3]]
            for j in range(2):
                if local[np.abs(local[:,j]).argmax(),j]<0:local[:,j]*=-1
                radius=np.max(np.abs(local[:,j]));local[:,j]/=radius if radius>0 else 1
        xy[ix]=local+np.array([3*(slot%width),3*(slot//width)])
    output=coordinates[['islet_id','case_id','component']].copy();output['graph_x']=np.nan;output['graph_y']=np.nan
    output['layout_kind']='componentwise spectral graph layout';output['graph_landmark']=output.islet_id.isin(model['landmark_ids'])
    valid=frame.eligible.to_numpy(bool)&np.isfinite(frame[model['columns']]).all(axis=1).to_numpy()
    if valid.any():
        vectors=(frame.loc[valid,model['columns']].to_numpy(float)-model['center'])/model['scale']
        distances,neighbors=NearestNeighbors(n_neighbors=model['k'],metric='cosine',algorithm='brute').fit(model['landmark_vectors']).kneighbors(vectors)
        for row,ds,ns in zip(np.flatnonzero(valid),distances,neighbors):
            same=comp[ns]==comp[ns[0]];weights=1/np.maximum(ds[same],1e-10);weights/=weights.sum()
            point=np.sum(xy[ns[same]]*weights[:,None],axis=0)
            output.loc[row,['graph_x','graph_y']]=point
    lookup=output.set_index('islet_id').index.get_indexer(model['landmark_ids'])
    output.loc[lookup,['graph_x','graph_y']]=xy
    return output


def concordance(reference,other):
    m=reference[['islet_id','pseudotime_raw']].merge(other[['islet_id','pseudotime_raw']],on='islet_id',suffixes=('_a','_b')).dropna()
    r=float(spearmanr(m.pseudotime_raw_a,m.pseudotime_raw_b).statistic) if len(m)>3 else np.nan
    return dict(shared_supported=len(m),spearman=r if np.isfinite(r) else None,absolute_spearman=abs(r) if np.isfinite(r) else None)

def _validation_one(frame,reference,family,name,params,training,feature_columns):
    start=time.monotonic()
    try:
        with threadpool_limits(limits=1):
            fitted,model=fit(training,feature_columns=feature_columns,**params)
            projection=project(frame,model)
        if family == 'lodo':
            target=projection.set_index('islet_id').index.get_indexer(fitted.islet_id)
            for column in ['pseudotime_raw','supported','component']:
                projection.loc[target,column]=fitted[column].to_numpy()
        score=concordance(reference,projection)
        score.update(coverage=float(projection.pseudotime_raw.notna().sum()/max(1,frame.eligible.sum())),
                     components=int(len(np.unique(model['landmark_components']))),root=model['root_id'],success=True)
        if family=='lodo':
            held=projection.loc[projection.case_id==name]
            score['held_out_coverage']=float(held.pseudotime_raw.notna().sum()/max(1,held.eligible.sum()))
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as e:
        score=dict(error=str(e),absolute_spearman=None,coverage=0.,success=False)
    return dict(family=family,replicate=name,seconds=time.monotonic()-start,**score)


def validate(frame,reference,bootstrap=200,feature_columns=None,workers=1,checkpoint_root=None):
    from joblib import Parallel,delayed
    donors=sorted(frame.case_id.unique());rng=np.random.default_rng(SEED)
    specifications=[('lodo',str(d),dict(),frame.loc[frame.case_id!=d]) for d in donors]
    specifications += [('sensitivity',name,params,frame) for name,params in [
        ('k10',{'k':10}),('k30',{'k':30}),('cap100',{'cap':100}),
        ('seed20260917',{'seed':20260917}),('seed20260918',{'seed':20260918})]]
    for b in range(bootstrap):
        sample=[]
        for j,d in enumerate(rng.choice(donors,len(donors),replace=True)):
            g=frame.loc[frame.case_id==d].copy();g.case_id=f'{d}@{j}';g.islet_id=g.islet_id+f'@{j}';sample.append(g)
        specifications.append(('bootstrap',str(b),{},pd.concat(sample,ignore_index=True)))
    checkpoint=Path(checkpoint_root) if checkpoint_root else None
    if checkpoint:checkpoint.mkdir(parents=True,exist_ok=True)
    columns=list(feature_columns) if feature_columns is not None else sorted(c for c in frame if c.startswith(('core__count__','peri__count__')))
    payload=pd.util.hash_pandas_object(frame[['islet_id','case_id','area_um2','eligible']+columns],index=False).to_numpy().tobytes()
    fingerprint=hashlib.sha256(payload+json.dumps(columns).encode()+sha256(__file__).encode()).hexdigest()
    def run(spec):
        family,name,params,training=spec;path=checkpoint/f'{family}-{name}.json' if checkpoint else None
        if path and path.exists():
            saved=json.loads(path.read_text())
            if saved.get('input_fingerprint')==fingerprint:return saved['result']
        record=_validation_one(frame,reference,family,name,params,training,feature_columns)
        if path:write_json(path,dict(input_fingerprint=fingerprint,result=record))
        print(f'{family} {name}: {record.get("absolute_spearman")}',flush=True)
        return record
    return Parallel(n_jobs=workers)(delayed(run)(spec) for spec in specifications)


def apply_identity_review(frame, review, geometry_hashes):
    required={'islet_id','case_id','decision','geometry_sha256','review_note'}
    if not required.issubset(review) or review.islet_id.duplicated().any():
        raise ValueError('Islet identity decisions must uniquely identify every source region')
    if set(zip(review.islet_id,review.case_id)) != set(zip(frame.islet_id,frame.case_id)):
        raise ValueError('Islet identity review does not match the feature cohort')
    if not review.decision.isin(['confirmed_islet','non_islet','uncertain','unreviewed']).all():
        raise ValueError('Unknown islet identity decision')
    if not (review.geometry_sha256==review.case_id.map(geometry_hashes)).all():
        raise ValueError('Reviewed geometry differs from the feature source; rebuild geometry-dependent measurements')
    if review.loc[review.decision=='confirmed_islet','review_note'].fillna('').str.strip().eq('').any():
        raise ValueError('Confirmed islets require a recorded review note')
    out=frame.merge(review[['islet_id','decision']],on='islet_id',validate='one_to_one',sort=False)
    out['prescreen_eligible']=out.eligible
    out['eligible']=out.eligible & out.decision.eq('confirmed_islet')
    out.loc[out.decision!='confirmed_islet','exclusion_reason']='islet_identity_'+out.loc[out.decision!='confirmed_islet','decision']
    if set(out.loc[out.eligible,'case_id']) != set(frame.case_id):
        raise ValueError('Every cohort donor needs reviewed eligible islets before a new all-donor fit')
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--output-root',required=True)
    p.add_argument('--islet-review',help='Optional geometry-bound review; no blanket fresh review requirement')
    p.add_argument('--bootstrap',type=int,default=200);p.add_argument('--workers',type=int,default=1)
    args=p.parse_args();source=Path(args.input_root);out=Path(args.output_root);out.mkdir(parents=True,exist_ok=True)
    feature_path=source/'trajectory/composition_features.parquet';manifest_path=source/'trajectory/feature_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    if manifest.get('schema_version')!=4 or manifest.get('representation')!='typed_composition':
        raise ValueError('Runtime v4 requires the explicit post-typing composition feature manifest')
    if manifest.get('features_sha256')!=sha256(feature_path):raise ValueError('Composition features changed after manifest creation')
    frame=pd.read_parquet(feature_path);frame['case_id']=frame.case_id.astype(str);columns=manifest['columns']
    if (out/'coordinates.csv').exists():raise ValueError('Preserve the previous fit; select a new trajectory output directory')
    if args.islet_review:
        review=pd.read_csv(args.islet_review,dtype={'case_id':str})
        sources=pd.read_csv(source/'source_manifest.csv',dtype=str)
        geometry_hashes={r.case_id:sha256(r.geometry) for r in sources.itertuples(index=False)}
        frame=apply_identity_review(frame,review,geometry_hashes)
    primary_columns=manifest.get('primary_columns')
    if not primary_columns or not set(primary_columns).issubset(frame):
        raise ValueError('Primary classifier-core composition features are missing')
    feature_supported=frame.eligible & np.isfinite(frame[primary_columns]).all(axis=1)
    if 'trajectory_feature_supported' in frame and not np.array_equal(feature_supported.to_numpy(),frame.trajectory_feature_supported.fillna(False).to_numpy(bool)):
        raise ValueError('Recorded trajectory feature support differs from the primary feature matrix')
    eligible_counts=frame.loc[feature_supported].groupby('case_id').size()
    # Fit the fixed specification on donors with enough eligible islets for
    # k=15 landmark neighbors. Donors without that support remain explicit
    # unsupported observations; they are never forced onto the axis.
    fit_donors=eligible_counts[eligible_counts >= 2].index.astype(str)
    unsupported=sorted(set(frame.case_id)-set(fit_donors))
    if len(fit_donors) < 2:
        write_json(out/'validation.json',dict(status='unavailable',reason='Insufficient supported donors for k=15',unsupported_donors=unsupported))
        raise ValueError(f'Insufficient supported donors for trajectory fit: {unsupported}')
    frame=frame[frame.case_id.isin(fit_donors)].copy()
    # Eligibility comes from the classifier geometry and core-cell census.
    # Typing support controls only whether an eligible row can be represented
    # in this typed-composition fit; it never changes the eligible flag.
    fit_columns=list(primary_columns)
    if not feature_supported.any(): raise ValueError('No feature-supported eligible islets')
    with threadpool_limits(limits=1):coords,model=fit(frame,feature_columns=fit_columns)
    coords.to_csv(out/'coordinates.csv',index=False);coords.to_parquet(out/'coordinates.parquet',index=False)
    layout=graph_layout(frame,coords,model);layout.to_csv(out/'graph_layout.csv',index=False);layout.to_parquet(out/'graph_layout.parquet',index=False)
    sparse.save_npz(out/'graph.npz',model['graph']);coo=sparse.triu(model['graph'],k=1).tocoo()
    pd.DataFrame({'from':model['landmark_ids'][coo.row],'to':model['landmark_ids'][coo.col],'weight':coo.data}).to_csv(out/'edges.csv',index=False)
    arrays={k:v for k,v in model.items() if isinstance(v,np.ndarray)};np.savez_compressed(out/'model.npz',**arrays)
    spec={k:v for k,v in model.items() if not isinstance(v,(np.ndarray,sparse.spmatrix))}
    spec.update(schema_version=4,method='anisotropic cosine kNN diffusion, multiscale distance, Nystrom extension',
                orientation='arbitrary hash endpoint; not elapsed disease time',runner_sha256=sha256(__file__),
                features_sha256=sha256(feature_path),feature_manifest_sha256=sha256(manifest_path),
                versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','scikit-learn','shapely']},
                preprocessing='post-typing log1p(validated base-type count/classifier-defined core cell count); landmark median/IQR scaling',
                eligibility='no image-edge contact and at least one canonical cell in the endocrine pixel-classifier core; typing support cannot revoke eligibility',
                reference_percentile='frozen empirical mid-distribution: (left insertion + right insertion)/(2*N)',
                graph_layout='Componentwise spectral graph layout with fixed neighbor interpolation, not comparable diffusion coordinates',
                projection_support='all k neighbors in same component; nearest distance <= 2 times local kth-neighbor radius')
    write_json(out/'specification.json',spec)
    pd.DataFrame({'islet_id':model['landmark_ids'],'component':model['landmark_components']}).to_csv(out/'landmarks.csv',index=False)
    write_json(out/'validation.json',dict(status='running',bootstrap_required=200,eligible=int(frame.eligible.sum())))
    records=validate(frame,coords,args.bootstrap,feature_columns=fit_columns,workers=args.workers,checkpoint_root=out/'validation_checkpoints')
    state_path=source/'trajectory/state_composition_features.parquet'
    if manifest.get('state_sensitivity_columns'):
        if not state_path.exists() or sha256(state_path)!=manifest.get('state_features_sha256'):raise ValueError('State sensitivity manifest lacks its matching materialized count features')
        state_frame=pd.read_parquet(state_path)
        state_cols=[c for c in manifest['state_sensitivity_columns'] if c in state_frame and np.isfinite(state_frame.loc[state_frame.eligible,c]).all()]
        if state_cols: records.append(_validation_one(state_frame,coords,'sensitivity','state_resolved',{},state_frame,state_cols))
    pd.DataFrame(records).to_csv(out/'validation_replicates.csv',index=False)
    lodo=[r for r in records if r['family']=='lodo'];boots=[r for r in records if r['family']=='bootstrap']
    median=float(np.median([r['absolute_spearman'] if r.get('absolute_spearman') is not None else 0 for r in lodo]))
    coverage=float(coords.pseudotime_raw.notna().sum()/frame.eligible.sum())
    status='unstable' if median<.8 or coverage<.8 else 'stable_by_prespecified_checks'
    successful=sum(r['success'] for r in boots);lodo_success=sum(r['success'] for r in lodo)
    if successful<200 or lodo_success<len(lodo):status='incomplete_validation'
    groups=[g.pseudotime_raw.dropna().to_numpy() for _,g in coords.groupby('case_id')];groups=[g for g in groups if len(g)]
    association=float(kruskal(*groups).statistic) if len(groups)>1 else None
    groupcols=[c for c in ['case_id','donor_source','donor_status'] if c in coords]
    coords.groupby(groupcols,dropna=False).agg(islets=('islet_id','size'),supported=('pseudotime_raw','count'),median_pseudotime=('pseudotime_raw','median')).reset_index().to_csv(out/'donor_associations.csv',index=False)
    write_json(out/'validation.json',dict(status=status,median_lodo_absolute_spearman=median,supported_primary_fraction=coverage,
               eligible=int(frame.eligible.sum()),feature_supported=int(feature_supported.sum()),supported=int(coords.pseudotime_raw.notna().sum()),
               fit_donors=sorted(fit_donors),unsupported_donors=unsupported,
               disconnected_components=int(len(np.unique(model['landmark_components']))),bootstrap_completed=successful,
               bootstrap_requested=args.bootstrap,bootstrap_failed=len(boots)-successful,lodo_completed=lodo_success,lodo_failed=len(lodo)-lodo_success,
               donor_kruskal_statistic_descriptive=association,
               acquisition_association='one acquisition per donor: acquisition and donor effects cannot be separated',
               interpretation='Ordering with arbitrary orientation. Stability does not establish absence of bias or elapsed disease time.'))

if __name__=='__main__': main()
