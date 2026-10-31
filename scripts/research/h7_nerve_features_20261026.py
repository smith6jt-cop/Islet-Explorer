#!/usr/bin/env python3
"""INS-independent B3TUBB structures and observed CD31 proximity.

Per-image thresholds and small-object filtering follow Richardson2023.
No inherited INS-dependent nerve classifier or spatial randomization is used.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,argparse,json,hashlib
from pathlib import Path
import cv2,numpy as np,pandas as pd,shapely
from shapely.geometry import shape
from shapely.affinity import scale,translate
from scipy.ndimage import distance_transform_edt
from skimage.filters import threshold_otsu
import h1_features_20261022 as hf
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'uniform/typing'))
import codex_typing as ct
from redsea_cells import Channels
SRC=Path('data/development/uniform-20donor-20261015')
APP=Path('data/development/uniform-20donor-20261022-v10')
INPUT=Path('output/h5_h10_original_plan_20261026/h5/object_ins_features.parquet')
OUT=Path('output/h5_h10_original_plan_20261026/h7');PX=hf.PIXEL_UM
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pixel_mask(geom,shape_yx,x0=0,y0=0,step=1):
    mask=np.zeros(shape_yx,np.uint8)
    g=translate(scale(geom,xfact=1/PX,yfact=1/PX,origin=(0,0)),xoff=-x0,yoff=-y0)
    if step!=1:g=scale(g,xfact=1/step,yfact=1/step,origin=(0,0))
    def paint(p):
        if p.is_empty:return
        if p.geom_type=='Polygon':
            cv2.fillPoly(mask,[np.round(np.asarray(p.exterior.coords)*16).astype(np.int32)],1,shift=4)
            for ring in p.interiors:cv2.fillPoly(mask,[np.round(np.asarray(ring.coords)*16).astype(np.int32)],0,shift=4)
        elif hasattr(p,'geoms'):
            for part in p.geoms:paint(part)
        # Zero-area line/point remnants from clipping are not image area.
    paint(g);return mask.astype(bool)
def geoms(path):
    return np.array([shapely.make_valid(scale(shape(f['geometry']),xfact=PX,yfact=PX,origin=(0,0)))
      for f in json.loads(Path(path).read_text())['features']],dtype=object)
def component_filter(binary):
    _,labels,stats,_=cv2.connectedComponentsWithStats(binary.astype(np.uint8),connectivity=8)
    keep=stats[:,cv2.CC_STAT_AREA]*PX**2>=2;keep[0]=False
    return keep[labels]
def extract(row,core,tissue,endo,etree,vessels,vtree,b3,gap,b3t,gapt):
    bounds=core.buffer(42).bounds;H,W=b3.shape
    x0=max(0,int(np.floor(bounds[0]/PX)));y0=max(0,int(np.floor(bounds[1]/PX)))
    x1=min(W,int(np.ceil(bounds[2]/PX))+1);y1=min(H,int(np.ceil(bounds[3]/PX))+1)
    raw=np.asarray(b3[y0:y1,x0:x1]);graw=np.asarray(gap[y0:y1,x0:x1])
    region=core.buffer(20).difference(core).intersection(tissue)
    em=pixel_mask(hf.union_local(etree,endo,core.buffer(42)),raw.shape,x0,y0)
    tm=pixel_mask(tissue.intersection(core.buffer(42)),raw.shape,x0,y0)
    peri=pixel_mask(region,raw.shape,x0,y0)&~em
    nerve=component_filter((raw>b3t)&tm&~em)
    vm=pixel_mask(hf.union_local(vtree,vessels,core.buffer(42)),raw.shape,x0,y0)&tm
    distance=distance_transform_edt(~vm,sampling=PX) if vm.any() else np.full(raw.shape,np.inf)
    nmask=nerve&peri;gnmask=nmask&(graw>gapt);n=int(nmask.sum());area=int(peri.sum())*PX**2
    rec=dict(islet_id=row.islet_id,case_id=int(row.case_id),donor_status=row.donor_status,
        peri_assessed_area_um2=area,nerve_area_um2=n*PX**2,cd31_area_um2=int((vm&peri).sum())*PX**2,
        nerve_area_fraction=n*PX**2/area if area else np.nan,
        cd31_area_fraction=int((vm&peri).sum())*PX**2/area if area else np.nan,
        gap43_supported_nerve_fraction=float(gnmask.sum()/n) if n else np.nan,
        gap43_supported_nerve_area_um2=int(gnmask.sum())*PX**2,
        nerve_vessel_overlap=float(vm[nmask].mean()) if n else np.nan,
        mean_nerve_vessel_distance_capped20_um=float(np.minimum(distance[nmask],20).mean()) if n else np.nan,
        proximity_available=bool(n and area))
    for r in (5,10,20):
        rec[f'nerve_near_cd31_{r}um']=float((distance[nmask]<=r).mean()) if n else np.nan
        rec[f'gap43_nerve_near_cd31_{r}um']=float((distance[gnmask]<=r).mean()) if gnmask.any() else np.nan
    return rec,(raw,nerve,vm,peri)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cases',type=int,nargs='*');a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'by_donor').mkdir(exist_ok=True)
    d=pd.read_parquet(INPUT);d['case_id']=d.case_id.astype(int);manifest=ct.manifest(SRC)
    for tag,entry in sorted(manifest.items()):
        case=int(entry['case_id'])
        if a.cases is not None and case not in a.cases:continue
        dest=OUT/'by_donor'/f'{case}.csv';jp=dest.with_suffix('.json')
        if dest.exists() and jp.exists():print(tag,'complete',flush=True);continue
        sub=d.loc[d.case_id.eq(case)&d.neighborhood_complete_20].copy()
        cores,_,tissue,_,gp=hf.geometry(APP,case)
        ep=SRC/'islets'/f'{tag}.geojson';vp=SRC/'structures/Capillary'/f'Capillary__{tag}.geojson'
        endo=geoms(ep);etree=shapely.STRtree(endo);vessels=geoms(vp);vtree=shapely.STRtree(vessels)
        img=Channels(entry['image_path']);names=ct.channels(SRC,entry['entry_id'])
        b3=img.arr[names.index('B3TUBB')] if img.memmap else img(names.index('B3TUBB'))
        gap=img.arr[names.index('GAP43')] if img.memmap else img(names.index('GAP43'))
        tm=pixel_mask(tissue,b3[::8,::8].shape,step=8)
        b3t=float(threshold_otsu(np.asarray(b3[::8,::8])[tm]));gapt=float(threshold_otsu(np.asarray(gap[::8,::8])[tm]))
        rows=[];mid=sub.sort_values('core_area_um2').iloc[len(sub)//2].islet_id if len(sub) else ''
        for row in sub.itertuples():
            rec,arrays=extract(row,cores[row.base_islet_id][1],tissue,endo,etree,vessels,vtree,b3,gap,b3t,gapt)
            rows.append(rec)
            if row.islet_id==mid and tag in ['ND_6516','Aab_6450','T1D_6533']:
                raw,nerve,vm,peri=arrays
                np.savez_compressed(OUT/f'example_{case}.npz',raw=raw,nerve=nerve,cd31=vm,peri=peri,islet_id=row.islet_id,pixel_um=PX)
        frame=pd.DataFrame(rows);frame.to_csv(dest,index=False);st=Path(entry['image_path']).stat()
        jp.write_text(json.dumps(dict(case_id=case,tag=tag,b3tubb_threshold_raw=b3t,gap43_threshold_raw=gapt,
            source_image=entry['image_path'],source_image_bytes=st.st_size,source_image_mtime_ns=st.st_mtime_ns,
            channel_index_b3tubb=names.index('B3TUBB'),channel_index_gap43=names.index('GAP43'),
            minimum_component_area_um2=2,ins_channel_used=False,null_randomization_used=False,
            geometry_sha256=sha(gp),endocrine_geometry_sha256=sha(ep),cd31_geometry_sha256=sha(vp),
            objects=len(frame),nerve_present_objects=int(frame.proximity_available.sum())),indent=2)+'\n')
        print(tag,'objects',len(frame),'nerve present',int(frame.proximity_available.sum()),'threshold',b3t,flush=True)
        del b3,gap;img.close()
    if a.cases is not None:return
    features=pd.concat([pd.read_csv(p) for p in sorted((OUT/'by_donor').glob('*.csv'))],ignore_index=True)
    if features.case_id.nunique()!=20:raise ValueError('Incomplete cohort')
    result=d.merge(features.drop(columns=['case_id','donor_status']),on='islet_id',how='left',validate='one_to_one')
    result['proximity_exclusion_reason']=np.where(~result.neighborhood_complete_20,'incomplete_20um_acquisition_support',
      np.where(result.nerve_area_um2.fillna(0).le(0),'no_measured_peri_B3TUBB_structure',''))
    result.to_parquet(OUT/'object_neurovascular_features.parquet',index=False)
    result.to_csv(OUT/'object_neurovascular_features.csv',index=False)
    (OUT/'feature_input_provenance.json').write_text(json.dumps(dict(input=str(INPUT),sha256=sha(INPUT),objects=len(result),
      measured_proximity_objects=int(result.nerve_near_cd31_20um.notna().sum()),
      method='Per-image B3TUBB Otsu threshold; remove components below2um2; periTissue outside all endocrine cores; observed area proximity to existingCD31. No INS channel or spatial null.'),indent=2)+'\n')
if __name__=='__main__':main()
