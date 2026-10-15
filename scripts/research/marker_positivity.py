#!/usr/bin/env python3
"""Propose, image-review and apply pooled marker positivity calibrations.

A mixture split is a review candidate, never a validated biological gate. The
apply command requires separate image evidence and an explicit reviewed gate
for every available marker. Unresolved markers remain nullable, including in
empty compartments. No disease label or historical phenotype enters fitting.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.mixture import GaussianMixture
from common import MARKERS,sha256,write_json
from cell_typing import cell_eligibility

SEED=20260920
KEYS=['case_id','object_id']


def fit_candidate(marker,values):
    values=np.asarray(values,float)
    finite=values[np.isfinite(values)&(values>=0)]
    positive=finite[finite>0]
    record=dict(marker=marker,measured=len(finite),zero_fraction=float(np.mean(finite==0)),
                status='awaiting_image_review',threshold=None)
    if len(positive)<100 or len(np.unique(positive))<10:
        return dict(record,status='insufficient_distribution',components=None)
    x=np.log1p(positive)[:,None]
    one=GaussianMixture(1,n_init=1,random_state=SEED).fit(x)
    two=GaussianMixture(2,n_init=3,random_state=SEED,max_iter=200).fit(x)
    order=np.argsort(two.means_.ravel());means=two.means_.ravel()[order]
    variance=two.covariances_.ravel()[order];weights=two.weights_[order]
    grid=np.linspace(means.mean(),max(float(x.max()),means[1]+6*np.sqrt(variance[1])),20000)
    logs=np.log(weights)[None,:]-.5*(np.log(2*np.pi*variance)[None,:]+(grid[:,None]-means[None,:])**2/variance[None,:])
    posterior=np.exp(logs[:,1]-logsumexp(logs,axis=1))
    supported=grid[(posterior>=.95)&(grid>=means[1])]
    threshold=float(np.expm1(supported[0])) if len(supported) else None
    separation=float((means[1]-means[0])/np.sqrt(np.mean(variance)))
    record.update(threshold=threshold,positive_value_mixture=True,posterior_candidate=.95,
                  components=dict(log1p_means=means.tolist(),variances=variance.tolist(),weights=weights.tolist()),
                  bic_improvement=float(one.bic(x)-two.bic(x)),separation=separation,
                  quantiles=dict(zip(['q01','q10','q50','q90','q99','q999'],np.quantile(finite,[.01,.1,.5,.9,.99,.999]).tolist())),
                  caveat='Positive-tail Gaussian mixture only proposes source-image reference locations; bimodality and posterior do not validate positivity.')
    return record


def discover(root,out,workers):
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);frames=[];inputs={}
    eligibility_path=root/'typing/eligibility/technical_exclusions.json'
    technical=json.loads(eligibility_path.read_text())['regions']
    for cid in sorted(cohort.case_id):
        path=root/'corrected'/f'{cid}-cells.parquet';inputs[cid]=sha256(path)
        frame=pd.read_parquet(path,columns=KEYS+['Cell Area','x_source_px','y_source_px']+MARKERS)
        valid=cell_eligibility(frame,5,technical).typing_eligible
        frame=frame.loc[valid].copy()
        order=frame.object_id.map(lambda value:hashlib.sha256(f'{SEED}:{value}'.encode()).hexdigest())
        frames.append(frame.loc[order.sort_values().index[:20000]])
    target=min(map(len,frames));sample=pd.concat([f.iloc[:target] for f in frames],ignore_index=True)
    out.mkdir(parents=True,exist_ok=True)
    if (out/'candidate_specification.json').exists():raise FileExistsError('Preserve previous marker calibration')
    sample.to_parquet(out/'landmarks.parquet',index=False)
    candidates=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs=[pool.submit(fit_candidate,marker,sample[marker].to_numpy()) for marker in MARKERS]
        for job in as_completed(jobs):
            record=job.result();candidates.append(record);print(record['marker'],record['threshold'],flush=True)
    candidates.sort(key=lambda row:MARKERS.index(row['marker']))
    write_json(out/'candidate_specification.json',dict(schema_version=1,status='awaiting_image_review',
               seed=SEED,donors=sorted(cohort.case_id),landmarks_per_donor=target,cell_area_min_um2=5,
               eligibility_sha256=sha256(eligibility_path),technical_exclusion_regions=technical,
               input_hashes=inputs,landmark_sha256=sha256(out/'landmarks.parquet'),markers=candidates))
    write_json(out/'review_template.json',dict(candidate_specification_sha256=sha256(out/'candidate_specification.json'),
               reviewer=None,markers={m:dict(validated=False,threshold=None,evidence=[],reason='Image review pending') for m in MARKERS}))


def render(root,out,marker):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import tifffile
    import shapely
    from shapely.geometry import shape
    from audit_regenerated_boundaries import contour
    specification=json.loads((out/'candidate_specification.json').read_text())
    record=next(r for r in specification['markers'] if r['marker']==marker)
    threshold=record['threshold']
    reference_only=threshold is None
    if reference_only:threshold=record['quantiles']['q99']
    frame=pd.read_parquet(out/'landmarks.parquet')
    # Fixed intensity strata and distinct donor preference. These are candidates
    # near the proposed transition, plus clear high/low examples, not labels.
    selected=[];used=set();donor_uses={}
    for label,target in [('low',.25),('below',.75),('threshold',1.),('above',1.5),('high',3.),('extreme',6.)]:
        group=frame.assign(distance=np.abs(np.log1p(frame[marker])-np.log1p(threshold*target)))
        pool=group.nsmallest(500,'distance').copy()
        pool['donor_uses']=pool.case_id.map(donor_uses).fillna(0)
        pool=pool.loc[~pool.object_id.isin(used)].sort_values(['donor_uses','distance','object_id'])
        r=pool.iloc[0].to_dict();r['stratum']=label;selected.append(r)
        used.add(r['object_id']);donor_uses[r['case_id']]=donor_uses.get(r['case_id'],0)+1
    cohort=pd.read_csv(root/'cohort.csv',dtype=str).set_index('case_id')
    inputs=json.loads((root/'preprocessing/inputs.json').read_text())
    fig,axes=plt.subplots(6,3,figsize=(9,17),squeeze=False);evidence=[]
    for j,r in enumerate(selected):
        cid=r['case_id'];source=json.loads(Path(cohort.loc[cid,'boundary_source']).read_text())
        names=[{'PGP9_5':'PGP9.5'}.get(m,m) for m in source['source_channel_names']]
        channels=['DAPI',marker,'NaKATPase']
        indices=[names.index(m) for m in channels]
        # Native acquisition order is authoritative for this source image;
        # converted candidate OME order is not assumed identical.
        try:pixels=tifffile.memmap(source['source_image'],mode='r')
        except ValueError:pixels=tifffile.imread(source['source_image'],key=indices,maxworkers=1);indices=list(range(3))
        px=float(source['pixel_size_um']);half=14/px
        cx,cy=r['x_source_px'],r['y_source_px']
        x0,x1=max(0,int(cx-half)),min(pixels.shape[2],int(cx+half+1))
        y0,y1=max(0,int(cy-half)),min(pixels.shape[1],int(cy+half+1))
        geometry_file=Path(next(x['geometry'] for x in inputs['donors'] if x['case_id']==cid))
        bare=r['object_id'].split(':',1)[-1]
        cell=next(shape(f['geometry']) for f in json.loads(geometry_file.read_text())['features'] if str(f['id'])==bare)
        levels=[]
        for k,(channel,idx) in enumerate(zip(channels,indices)):
            lo,hi=np.quantile(pixels[idx,::16,::16],[.001,.9995]);levels.append([float(lo),float(hi)])
            patch=pixels[idx,y0:y1,x0:x1]
            ax=axes[j,k];ax.imshow(patch,cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
            contour(ax,cell,np.array([x0,y0]),color='cyan');ax.set_xticks([]);ax.set_yticks([])
            ax.set_title(channel,fontsize=9)
            if k==0:ax.set_ylabel(f'{r["stratum"]}; donor {cid}\ncorrected {r[marker]:.1f}; {r["Cell Area"]:.1f} µm²',fontsize=8)
        evidence.append(dict(case_id=cid,object_id=r['object_id'],stratum=r['stratum'],
                             corrected_value=r[marker],channel_ranges=levels,
                             source_specification_sha256=sha256(Path(cohort.loc[cid,'boundary_source']))))
    caption='distribution reference (no mixture gate)' if reference_only else 'unvalidated mixture candidate'
    fig.suptitle(f'{marker}: {caption} {threshold:.2f}\nRaw source pixels, DAPI and exact cell boundary; no disease labels; proposal only',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.96));path=out/'review'/f'{marker}.png';path.parent.mkdir(exist_ok=True)
    fig.savefig(path,dpi=130);plt.close(fig)
    write_json(out/'review'/f'{marker}.json',dict(marker=marker,candidate_threshold=None if reference_only else threshold,
               distribution_reference=threshold,reference_only=reference_only,
               panel_path=str(path.resolve()),panel_sha256=sha256(path),examples=evidence,
               acceptance='pending_direct_image_inspection'))
    print(marker,flush=True)


def reviewed_specification(out,review):
    document=json.loads(review.read_text());candidate=out/'candidate_specification.json'
    if document['candidate_specification_sha256']!=sha256(candidate):raise ValueError('Review is bound to other marker calibration')
    if set(document['markers'])!=set(MARKERS):raise ValueError('Explicit review disposition required for all31 markers')
    if not document.get('reviewer'):raise ValueError('Missing reviewer')
    for marker,record in document['markers'].items():
        if not record.get('reason'):raise ValueError('Missing review explanation')
        if not record['validated']:continue
        if not isinstance(record['threshold'],(int,float)) or record['threshold']<=0:raise ValueError('Invalid marker gate')
        if not record['evidence']:raise ValueError('Marker positivity requires image evidence')
        for evidence in record['evidence']:
            if sha256(evidence['path'])!=evidence['sha256']:raise ValueError('Marker evidence changed')
    return document


def apply(root,out,review):
    specification=reviewed_specification(out,review)
    candidates=json.loads((out/'candidate_specification.json').read_text())
    for cid,expected in candidates['input_hashes'].items():
        path=root/'corrected'/f'{cid}-cells.parquet'
        if sha256(path)!=expected:raise ValueError('Corrected measurements changed since gate review')
        frame=pd.read_parquet(path,columns=KEYS+['Cell Area','x_source_px','y_source_px']+MARKERS)
        result=frame[KEYS].copy();eligible=cell_eligibility(frame,5,candidates['technical_exclusion_regions']).typing_eligible
        for marker,record in specification['markers'].items():
            calls=pd.Series(pd.NA,index=frame.index,dtype='boolean')
            if record['validated']:calls.loc[eligible]=frame.loc[eligible,marker].ge(record['threshold'])
            result['positive__'+marker]=calls
        result.to_parquet(out/f'{cid}.parquet',index=False)
    write_json(out/'manifest.json',dict(status='reviewed_calibration_applied',review_path=str(review.resolve()),
               review_sha256=sha256(review),source_specification_sha256=sha256(out/'candidate_specification.json'),
               validated_markers=[m for m,r in specification['markers'].items() if r['validated']],
               unavailable_markers=[m for m,r in specification['markers'].items() if not r['validated']],
               cell_area_min_um2=5,disease_used=False))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['discover','render','apply'])
    p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--marker',action='append');p.add_argument('--review',type=Path)
    a=p.parse_args();root=a.input_root.resolve(strict=True);out=a.output_root.resolve()
    if a.command=='discover':discover(root,out,a.workers)
    elif a.command=='render':
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for f in as_completed([pool.submit(render,root,out,m) for m in (a.marker or MARKERS)]):f.result()
    else:
        if a.review is None:p.error('--review required for apply')
        apply(root,out,a.review)
if __name__=='__main__':main()
