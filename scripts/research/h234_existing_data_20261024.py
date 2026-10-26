#!/usr/bin/env python3
"""Read-only structural evidence audit and H3 spatial test; no classifiers."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/tmp/h234-mpl')
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
from shapely.affinity import scale
from scipy import stats
import h1_features_20261022 as hf
import h1_analysis_20261022 as h1

OUT=Path('output/h1_followups_h2_h4_20261024')
ROOT=Path('data/development/uniform-20donor-20261022-v10')
OLD=Path('data/development/all-donor-20260920/structural')
SRC=Path('data/development/uniform-20donor-20261015')
CORRECTED=Path('output/h1_edge_corrected_20261024')


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read_geoms(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:j=json.load(f)
    gs=[shape(x['geometry']) for x in j['features']]
    return [shapely.make_valid(scale(g,xfact=hf.PIXEL_UM,yfact=hf.PIXEL_UM,origin=(0,0))) for g in gs]


def evidence_audit(d):
    provenance=json.loads((SRC/'structures/provenance.json').read_text())
    rows=[]; hashes={}
    for case in sorted(d.case_id.unique()):
        p=OLD/f'{case}.json'; registry=json.loads(p.read_text());hashes[str(p)]=sha(p)
        status=d.loc[d.case_id.eq(case),'donor_status'].iloc[0]
        current=next(v for k,v in provenance.items() if int(k.split('_')[-1])==int(case))
        for label,mask in registry['masks'].items():
            review=mask.get('review',{})
            rows.append(dict(case_id=case,donor_status=status,structure=label,
                decision=review.get('decision','unavailable'),validated=mask.get('validated',False),
                current_mask_review=current.get(label,{}).get('audit_2026_09_21','not_bound_to_current_mask'),
                mask_path=mask['path'],mask_sha256=mask['sha256'],rationale=review.get('rationale','')))
    pd.DataFrame(rows).to_csv(OUT/'structural_evidence_audit.csv',index=False)
    return pd.DataFrame(rows), hashes


def rigid_rotations(xy, center, core, tissue, rng, target=999, maximum=100000):
    """Conditional angular null; every accepted configuration stays in support."""
    centered=xy-center
    samples=[]; tried=0; accepted=0
    while accepted<target and tried<maximum:
        n=min(2000,maximum-tried)
        a=rng.uniform(0,2*np.pi,n);c=np.cos(a);s=np.sin(a)
        rot=np.stack([c[:,None]*centered[None,:,0]-s[:,None]*centered[None,:,1],
                      s[:,None]*centered[None,:,0]+c[:,None]*centered[None,:,1]],axis=-1)+center
        pts=shapely.points(rot)
        good=(shapely.covers(core,pts)&shapely.covers(tissue,pts)).all(axis=1)
        samples.extend(rot[good]);accepted+=int(good.sum());tried+=n
    return np.array(samples[:target]),tried,accepted


def h3(d,audit,hashes):
    accepted=audit.loc[audit.structure.eq('Capillary')&audit.validated]
    records=[]; counts=[]; geometry_proof=[]
    for ev in accepted.to_dict('records'):
        case=ev['case_id'];sub=d.loc[d.case_id.eq(case)&d.n_immune_core.gt(0)]
        path=Path(ev['mask_path'])
        if sha(path)!=ev['mask_sha256']:raise ValueError('Accepted mask digest mismatch')
        hashes[str(path)]=sha(path)
        vessels=read_geoms(path)
        if ev['current_mask_review']=='accepted':
            current_path=next((SRC/'structures/Capillary').glob(f'*_{int(case):04d}.geojson'))
            current=read_geoms(current_path)
            identical={g.wkb for g in vessels}=={g.wkb for g in current}
            if not identical:raise ValueError('Current mask differs from reviewed geometry')
            hashes[str(current_path)]=sha(current_path)
            geometry_proof.append(dict(case_id=case,current_path=str(current_path),
                current_sha256=sha(current_path),reviewed_path=str(path),geometry_sets_identical=identical))
        tree=shapely.STRtree(vessels)
        cores,_,tissue,_,gp=hf.geometry(ROOT,int(case));hashes[str(gp)]=sha(gp)
        cp=CORRECTED/'analysis_cells_by_donor'/f'{case}.parquet'
        mp=CORRECTED/'object_memberships_by_donor'/f'{case}.parquet'
        hashes[str(cp)]=sha(cp);hashes[str(mp)]=sha(mp)
        cells=pd.read_parquet(cp);mem=pd.read_parquet(mp)
        immune=mem.loc[mem.cell_region.eq('core')&mem.phenotype.isin(hf.IMMUNE)].merge(
            cells[['object_id','x_um','y_um']],on='object_id',validate='many_to_one')
        for rec in sub.to_dict('records'):
            core=cores[rec['base_islet_id']][1]
            pts=immune.loc[immune.islet_id.eq(rec['islet_id']),['x_um','y_um']].to_numpy()
            if len(pts)!=rec['n_immune_core']:raise ValueError('Core immune count mismatch')
            row=dict(case_id=case,donor_status=ev['donor_status'],islet_id=rec['islet_id'],
                     core_immune=len(pts), current_mask_accepted=ev['current_mask_review']=='accepted')
            if not (shapely.covers(core,shapely.points(pts))&shapely.covers(tissue,shapely.points(pts))).all():
                counts.append(dict(row,status='observed_points_outside_null_support'));continue
            crossings=[]
            for k in tree.query(core,predicate='intersects'):
                vessel=vessels[k]
                if vessel.intersection(core).area>1e-6 and vessel.difference(core).area>1e-6:
                    cross=vessel.intersection(core.boundary)
                    if not cross.is_empty:crossings.append(cross)
            row['crossing_components']=len(crossings)
            if not crossings:
                counts.append(dict(row,status='no_boundary_crossing'));continue
            crossing=shapely.union_all(crossings)
            center=np.array([core.centroid.x,core.centroid.y])
            seed=int(hashlib.sha256(rec['islet_id'].encode()).hexdigest()[:8],16)
            rotations,tried,naccepted=rigid_rotations(pts,center,core,tissue,np.random.default_rng(seed))
            row.update(null_draws=len(rotations),null_candidates=tried,null_acceptance=naccepted/tried)
            if len(rotations)!=999:
                counts.append(dict(row,status='insufficient_null_orientations'));continue
            observed_dist=shapely.distance(shapely.points(pts),crossing)
            null_dist=shapely.distance(shapely.points(rotations),crossing)
            for radius in [5,10,20]:
                observed=float(np.mean(observed_dist<=radius))
                null_per_draw=np.mean(null_dist<=radius,axis=1)
                records.append(dict(row,radius_um=radius,observed=observed,
                    expected=float(null_per_draw.mean()),excess=observed-float(null_per_draw.mean()),
                    null_sd=float(null_per_draw.std(ddof=1)),
                    null_mcse=float(null_per_draw.std(ddof=1)/np.sqrt(len(rotations)))))
            counts.append(dict(row,status='evaluated'))
        print(f'H3 donor {case}: {len(sub)} core-positive objects checked',flush=True)
    result=pd.DataFrame(records);result.to_csv(OUT/'h3_object_spatial_null.csv',index=False)
    pd.DataFrame(counts).to_csv(OUT/'h3_spatial_eligibility.csv',index=False)
    donors=result.groupby(['case_id','donor_status','current_mask_accepted','radius_um']).agg(
        objects=('islet_id','size'),core_immune=('core_immune','sum'),observed=('observed','mean'),
        expected=('expected','mean'),excess=('excess','mean')).reset_index()
    donors.to_csv(OUT/'h3_donor_effects.csv',index=False)
    summary=[]
    for selection in ['current_reviewed_masks','all_retained_reviewed_masks']:
        z=donors.loc[donors.current_mask_accepted] if selection=='current_reviewed_masks' else donors
        for group in ['Pooled',*h1.STATUS]:
            for radius in [5,10,20]:
                g=z.loc[z.radius_um.eq(radius)&(True if group=='Pooled' else z.donor_status.eq(group))]
                n=len(g);mu=g.excess.mean()
                se=g.excess.std(ddof=1)/np.sqrt(n) if n>1 else np.nan
                q=stats.t.ppf(.975,n-1) if n>1 else np.nan
                summary.append(dict(selection=selection,group=group,radius_um=radius,donors=n,
                    objects=g.objects.sum(),core_immune=g.core_immune.sum(),excess=mu,
                    lower=mu-q*se,upper=mu+q*se,
                    p=2*stats.t.sf(abs(mu/se),n-1) if n>1 and se>0 else np.nan))
    pd.DataFrame(summary).to_csv(OUT/'h3_summary.csv',index=False)
    pd.DataFrame(geometry_proof).to_csv(OUT/'h3_current_mask_geometry_verification.csv',index=False)


def run():
    d=h1.eligible(pd.read_parquet(CORRECTED/'h1_object_features.parquet'))
    audit,hashes=evidence_audit(d)
    h3(d,audit,hashes)
    (OUT/'h234_source_hashes.json').write_text(json.dumps(hashes,indent=2))


if __name__=='__main__':run()
