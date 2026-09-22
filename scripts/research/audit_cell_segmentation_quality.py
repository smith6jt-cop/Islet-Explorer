#!/usr/bin/env python3
"""Render source-pixel evidence for cell fragments; no automatic biological labels."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile
from common import sha256,write_json
from audit_regenerated_boundaries import contour


def render(root,out,cid):
    source=json.loads((root/'boundaries/sources'/cid/'source.json').read_text())
    source_image=Path(source['source_image'])
    px=float(source['pixel_size_um'])
    inputs=json.loads((root/'preprocessing/inputs.json').read_text())
    geom_file=Path(next(r['geometry'] for r in inputs['donors'] if r['case_id']==cid))
    columns=['case_id','object_id','x_source_px','y_source_px','Cell Area','Nucleus Area','INS','GCG','SST','CHGA']
    frame=pd.read_parquet(root/'corrected'/f'{cid}-cells.parquet',columns=columns)
    frame['area_um2']=frame['Cell Area']
    frame['nucleus_area_px2']=frame['Nucleus Area']/px**2
    frame['selection_hash']=frame.object_id.map(lambda s:hashlib.sha256(f'20260920:{s}'.encode()).hexdigest())
    frame['hormone_sum']=frame[['INS','GCG','SST']].sum(axis=1,min_count=1)
    selected=[]
    for lo,hi in [(0,1),(1,5),(5,10),(20,80)]:
        group=frame.loc[frame.area_um2.ge(lo)&frame.area_um2.lt(hi)].copy()
        if not len(group):continue
        ranks=[group.sort_values('area_um2').iloc[len(group)//2],group.sort_values('selection_hash').iloc[0],group.sort_values('hormone_sum',ascending=False).iloc[0]]
        used=set()
        for r,reason in zip(ranks,['median_area','deterministic_random','largest_hormone_sum']):
            if r.object_id in used:continue
            used.add(r.object_id);record=r.to_dict();record['selection_reason']=f'{lo}–{hi} µm²: {reason}';selected.append(record)
    table=pd.DataFrame(selected);wanted=set(table.object_id.str.split(':',n=1).str[-1]);geometries={str(f['id']):shape(f['geometry']) for f in json.loads(geom_file.read_text())['features'] if str(f['id']) in wanted}
    if len(geometries)!=len(wanted):raise ValueError('Cell geometry UUID mismatch')
    markers=['DAPI','INS','GCG','SST','CHGA','NaKATPase'];indices=[source['source_channel_names'].index(m) for m in markers]
    try: pixels=tifffile.memmap(source_image,mode='r')
    except ValueError: pixels=tifffile.imread(source_image,key=indices,maxworkers=1);indices=list(range(len(markers)))
    levels=np.array([np.quantile(pixels[i,::16,::16],[.001,.9995]) for i in indices])
    folder=out/cid;folder.mkdir(parents=True,exist_ok=True);table.to_csv(folder/'selected_cells.csv',index=False)
    evidence=[]
    for page,start in enumerate(range(0,len(table),6)):
        rows=table.iloc[start:start+6];fig,axes=plt.subplots(len(rows),4,figsize=(12,2.6*len(rows)),squeeze=False)
        for j,r in enumerate(rows.itertuples()):
            g=geometries[r.object_id.split(':',1)[-1]];cx,cy=r.x_source_px,r.y_source_px;half=max(12/px,(max(g.bounds[2]-g.bounds[0],g.bounds[3]-g.bounds[1])/2)+5/px)
            x0,x1=max(0,int(cx-half)),min(pixels.shape[2],int(cx+half+1));y0,y1=max(0,int(cy-half)),min(pixels.shape[1],int(cy+half+1));patch=np.array([pixels[i,y0:y1,x0:x1] for i in indices]);normalized=np.clip((patch-levels[:,0,None,None])/np.maximum(levels[:,1,None,None]-levels[:,0,None,None],1),0,1)
            for k,display in enumerate([normalized[0],np.moveaxis(normalized[1:4],0,-1),normalized[4],normalized[5]]):
                ax=axes[j,k];ax.imshow(display,cmap='gray',vmin=0,vmax=1,interpolation='nearest');contour(ax,g,np.array([x0,y0]));ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([]);ax.set_title(['DAPI + cell boundary','INS / GCG / SST + boundary','CHGA + boundary','NaKATPase + boundary'][k],fontsize=8)
                if k==0:ax.set_ylabel(f'{r.selection_reason}\nCell {r.area_um2:.3f} µm²; nuclear {r.nucleus_area_px2*px**2:.3f} µm²\n{r.object_id[-12:]}',fontsize=8)
                ax.plot([3,3+5/px],[y1-y0-4]*2,color='white',lw=2);ax.text(3,y1-y0-7,'5 µm',color='white',fontsize=7)
        fig.suptitle(f'Image reference {cid}: cell segmentation quality review; no disease category shown\nFixed image-wide channel ranges. Tiny objects are retained in raw files; these images do not assign cell types.',fontsize=10);fig.tight_layout(rect=(0,0,1,.96));file=folder/f'cell_area_quality_{page+1}.png';fig.savefig(file,dpi=130);plt.close(fig);evidence.append({'path':str(file.resolve()),'sha256':sha256(file),'object_ids':rows.object_id.tolist()})
    record={'case_id':cid,'status':'rendered_awaiting_visual_review','source_image':str(source_image),'source_specification_sha256':sha256(root/'boundaries/sources'/cid/'source.json'),'geometry_sha256':sha256(geom_file),'examined_channels':markers,'fixed_full_image_display_ranges':dict(zip(markers,levels.tolist())),'selected_cells':len(table),'evidence':evidence,'source_pixels_modified':False,'biological_labels_assigned':False,'owner_approval_claimed':False};write_json(folder/'manifest.json',record);print(json.dumps({'case_id':cid,'panels':len(evidence),'cells':len(table)}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);p.add_argument('--case-id',action='append',required=True);a=p.parse_args()
    for cid in a.case_id:render(a.input_root,a.output_root,cid)
if __name__=='__main__':main()
