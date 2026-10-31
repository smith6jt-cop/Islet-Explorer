#!/usr/bin/env python3
"""Existing-cell INS decomposition and fixed endocrine-contact graph preparation.

Does not retype cells, rebuild images, or execute spatial label permutations.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import io
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/development/uniform-20donor-20261022-v10'
SOURCE=ROOT/'data/development/uniform-20donor-20261015'
OUT=ROOT/'output/h5_h10_original_plan_20261026'
BASE=ROOT/'output/h1_edge_corrected_20261024/h1_object_features.parquet'
ENDO=['Beta cell','Alpha cell','Delta cell','CHGA+ hormone-negative','Endocrine']
QC=['CHGA','NaKATPase','BActin','ECAD','BCatenin','CD44','HLADR']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def safe_mean(v):
    v=np.asarray(v,float);v=v[np.isfinite(v)]
    return float(v.mean()) if len(v) else np.nan
def safe_median(v):
    v=np.asarray(v,float);v=v[np.isfinite(v)]
    return float(np.median(v)) if len(v) else np.nan

def donor_extract(case,objects,oldhash):
    case=int(case);hashes={}
    apppath=DATA/'cell_types'/f'{case}.parquet'
    src=pd.read_parquet(apppath);hashes[str(apppath.relative_to(ROOT))]=sha(apppath)
    ix=pd.MultiIndex.from_frame(src[['x_um','y_um','phenotype']]);assert ix.is_unique
    lookup=pd.Series(src.object_id.to_numpy(),index=ix)
    paths=list((SOURCE/'typing/aligned').glob(f'*_{case:04d}.parquet'))
    assert len(paths)==1
    alignedpath=paths[0];tag=alignedpath.stem
    aline=pd.read_parquet(alignedpath,columns=['object_id','INS','ECAD','panCK']).set_index('object_id')
    redpath=SOURCE/'typing/redsea'/f'{tag}.parquet'
    red=pd.read_parquet(redpath,columns=['object_id','Cytoplasm__INS']).set_index('object_id')
    rawpath=SOURCE/'cells'/f'{tag}.parquet'
    raw=pd.read_parquet(rawpath,columns=['Object ID','Cytoplasm: INS: Mean']).set_index('Object ID')
    callpath=SOURCE/'typing/calls'/f'{tag}.parquet'
    calls=pd.read_parquet(callpath,columns=['object_id','cell_class','excluded_reason','islet_num'])
    for p in [alignedpath,redpath,rawpath,callpath]:hashes[str(p.relative_to(ROOT))]=sha(p)
    scale=pd.read_csv(SOURCE/'typing/scale.csv').set_index('marker')
    scaled={m:np.clip((aline[m]-scale.loc[m,'floor'])/scale.loc[m,'span'],0,1) for m in ['INS','ECAD','panCK']}
    ids=src.loc[src.phenotype.eq('Exocrine')&~src.cell_region.eq('core'),'object_id']
    epithelial=(scaled['ECAD'].reindex(ids).ge(.5)|scaled['panCK'].reindex(ids).ge(.5)).to_numpy()
    ids=ids.loc[epithelial]
    background={'case_id':case,'image_tag':tag,'background_cells':len(ids),
        'ins_background_scaled_median':safe_median(scaled['INS'].reindex(ids)),
        'ins_background_aligned_median':safe_median(aline.INS.reindex(ids)),
        'ins_background_raw_median':safe_median(raw['Cytoplasm: INS: Mean'].reindex(ids))}
    excluded=calls.loc[calls.cell_class.eq('Excluded')&calls.islet_num.notna()].copy()
    excluded['INS_scaled']=scaled['INS'].reindex(excluded.object_id).to_numpy()
    exgroups={name:d for name,d in excluded.groupby('islet_num')}
    rows=[];nodes=[];graphs=[]
    graphdir=OUT/'h6/graphs_by_donor'/str(case);graphdir.mkdir(parents=True,exist_ok=True)
    for rec in objects:
        path=DATA/'cells'/f"{rec['islet_id']}.csv";rawbytes=path.read_bytes()
        digest=hashlib.sha256(rawbytes).hexdigest()
        assert oldhash[str(path)]==digest
        hashes[str(path.relative_to(ROOT))]=digest
        d=pd.read_csv(io.BytesIO(rawbytes),float_precision='round_trip')
        d=d.loc[d.cell_region.eq('core')].copy()
        assert len(d)==rec['n_cells_core']
        key=pd.MultiIndex.from_frame(d[['X_centroid','Y_centroid','phenotype']].rename(columns={'X_centroid':'x_um','Y_centroid':'y_um'}))
        ids=lookup.reindex(key);assert ids.notna().all() and not ids.duplicated().any()
        d['object_id']=ids.to_numpy()
        av=scaled['INS'].reindex(d.object_id).to_numpy()
        # The served matrix is stored at float32 precision before CSV export.
        assert np.allclose(d.INS.to_numpy(),av,equal_nan=True,rtol=5e-7,atol=1e-7)
        d['INS_raw']=raw['Cytoplasm: INS: Mean'].reindex(d.object_id).to_numpy()
        d['INS_corrected']=red['Cytoplasm__INS'].reindex(d.object_id).to_numpy()
        beta=d.phenotype.eq('Beta cell');nonbeta=~beta;hn=d.phenotype.eq('CHGA+ hormone-negative')
        assert beta.sum()==rec['n_beta_core']
        row={'islet_id':rec['islet_id'],'case_id':case,**background}
        delta=np.abs(d.INS.to_numpy()-av)
        row['ins_source_scale_max_abs_difference']=float(np.nanmax(delta)) if np.isfinite(delta).any() else np.nan
        for name,mask in [('all',np.ones(len(d),bool)),('beta',beta),('nonbeta',nonbeta),('CHGA_hormone_negative',hn)]:
            v=d.loc[mask,'INS']
            row[f'n_{name}_ins_measured']=int(v.notna().sum())
            row[f'n_{name}_ins_missing']=int(v.isna().sum())
            row[f'ins_sum_{name}_core']=float(v.sum())
            row[f'ins_mean_{name}_core']=safe_mean(v)
            row[f'ins_median_{name}_core']=safe_median(v)
        row['n_CHGA_hormone_negative']=int(hn.sum())
        row['ins_mean_beta_raw']=safe_mean(d.loc[beta,'INS_raw'])
        row['ins_mean_beta_corrected']=safe_mean(d.loc[beta,'INS_corrected'])
        row['ins_mean_beta_margin4']=safe_mean(d.loc[beta & d.INS.ge(4*d[['GCG','SST']].max(axis=1)),'INS'])
        row['n_beta_margin4']=int((beta & d.INS.ge(4*d[['GCG','SST']].max(axis=1))).sum())
        row['beta_mean_cell_area_um2']=safe_mean(d.loc[beta,'Cell Area'])
        for m in QC:row[f'beta_mean_{m}']=safe_mean(d.loc[beta,m])
        beta_component=0 if beta.sum()==0 else row['n_beta_ins_measured']*row['ins_mean_beta_core']
        row['decomposition_absolute_error']=abs(row['ins_sum_all_core']-beta_component-row['ins_sum_nonbeta_core'])
        assert row['decomposition_absolute_error']<1e-10
        ex=exgroups.get(rec['base_islet_id'])
        row['n_source_excluded_core']=len(ex) if ex is not None else 0
        row['n_source_excluded_INS_measured']=int(ex.INS_scaled.notna().sum()) if ex is not None else 0
        row['n_source_excluded_INS_missing']=int(ex.INS_scaled.isna().sum()) if ex is not None else 0
        row['ins_source_excluded_sum']=float(ex.INS_scaled.sum()) if ex is not None else 0.
        marg=ex.loc[ex.excluded_reason.eq('hormone_unresolved_margin')] if ex is not None else None
        row['n_source_unresolved_margin']=len(marg) if marg is not None else 0
        row['n_source_unresolved_INS_measured']=int(marg.INS_scaled.notna().sum()) if marg is not None else 0
        row['n_source_unresolved_INS_missing']=int(marg.INS_scaled.isna().sum()) if marg is not None else 0
        row['ins_source_unresolved_margin_sum']=float(marg.INS_scaled.sum()) if marg is not None else 0.
        rows.append(row)
        e=d.loc[d.phenotype.isin(ENDO)].sort_values('object_id').copy().reset_index(drop=True)
        coords=e[['X_centroid','Y_centroid']].to_numpy(float)
        labels=pd.Categorical(e.phenotype,categories=ENDO).codes.astype(np.int8)
        edges=cKDTree(coords).query_pairs(20.,output_type='ndarray').astype(np.int32) if len(e)>1 else np.zeros((0,2),np.int32)
        n_beta=int((labels==0).sum());n_ad=int(np.isin(labels,[1,2]).sum())
        if len(edges):
            aa=labels[edges[:,0]];bb=labels[edges[:,1]]
            ba=int((((aa==0)&(bb==1))|((aa==1)&(bb==0))).sum())
            bd=int((((aa==0)&(bb==2))|((aa==2)&(bb==0))).sum())
        else:ba=bd=0
        eligible=(len(e)>=10 and n_beta>=2 and n_ad>=2 and len(edges)>0)
        graphs.append({'islet_id':rec['islet_id'],'case_id':case,'donor_status':rec['donor_status'],
          'n_endocrine_nodes':len(e),'n_beta_nodes':n_beta,'n_alpha_nodes':int((labels==1).sum()),
          'n_delta_nodes':int((labels==2).sum()),'n_edges':len(edges),'observed_beta_alpha_edges':ba,
          'observed_beta_delta_edges':bd,'observed_combined_fraction':(ba+bd)/len(edges) if len(edges) else np.nan,
          'graph_eligible':eligible,'complete_graph':len(edges)==len(e)*(len(e)-1)//2})
        if eligible:
            np.savez_compressed(graphdir/f"{rec['islet_id']}.npz",coordinates=coords,labels=labels,edges=edges,
                                 object_ids=e.object_id.astype(str).to_numpy(dtype=str))
        e['islet_id']=rec['islet_id'];e['case_id']=case
        nodes.append(e[['islet_id','case_id','object_id','phenotype','X_centroid','Y_centroid','INS','INS_raw','INS_corrected']])
    print(f'{case}: {len(rows)} INS objects, {sum(g["graph_eligible"] for g in graphs)} prepared contact graphs',flush=True)
    return rows,pd.concat(nodes,ignore_index=True),graphs,hashes,background

def run(args):return donor_extract(*args)
def main():
    (OUT/'h5').mkdir(parents=True,exist_ok=True);(OUT/'h6').mkdir(parents=True,exist_ok=True)
    base=pd.read_parquet(BASE)
    base=base.loc[base.core_complete&base.neighborhood_complete_20&base.neighborhood_complete_100].copy()
    audit=json.loads((ROOT/'output/h1_edge_corrected_20261024/feature_audit.json').read_text())
    tasks=[(int(case),d.to_dict('records'),audit['source_sha256']) for case,d in base.groupby('case_id')]
    rows=[];nodes=[];graphs=[];hashes={};bgs=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for rr,nn,gg,hh,bg in pool.map(run,tasks):rows+=rr;nodes.append(nn);graphs+=gg;hashes.update(hh);bgs.append(bg)
    f=pd.DataFrame(rows).merge(base,on=['islet_id','case_id'],validate='one_to_one')
    f['beta_density_mm2']=f.n_beta_core/f.core_area_um2*1e6
    f['CHGA_hormone_negative_density_mm2']=f.n_CHGA_hormone_negative/f.core_area_um2*1e6
    f['h5_primary_eligible']=f.n_beta_core.gt(0)&f.n_immune_nearby.gt(0)&np.isfinite(f.ins_mean_beta_core)
    f.to_parquet(OUT/'h5/object_ins_features.parquet',index=False);f.to_csv(OUT/'h5/object_ins_features.csv',index=False)
    pd.DataFrame(bgs).to_csv(OUT/'h5/donor_background.csv',index=False)
    pd.concat(nodes,ignore_index=True).to_parquet(OUT/'h6/core_endocrine_cells.parquet',index=False)
    pd.DataFrame(graphs).to_csv(OUT/'h6/graph_inventory.csv',index=False)
    hashes[str(BASE.relative_to(ROOT))]=sha(BASE)
    hashes[str((SOURCE/'typing/scale.csv').relative_to(ROOT))]=sha(SOURCE/'typing/scale.csv')
    (OUT/'h5/source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(f.groupby('donor_status').agg(objects=('islet_id','size'),beta_positive=('n_beta_core',lambda s:(s>0).sum()),
        h5_primary=('h5_primary_eligible','sum'),missing_beta_ins=('n_beta_ins_missing','sum')).to_string())

if __name__=='__main__':main()
