#!/usr/bin/env python3
"""Efficient, disease-blinded source-image reference panels for typing clusters.

Only requested TIFF tiles are decoded. Full image-wide display ranges come from
the coarsest existing pyramid level and remain fixed across cells. Marker means
and cluster profiles are REDSEA-corrected; images show immutable source pixels.
Rendering never validates a biological annotation.
"""
from __future__ import annotations
import argparse
from collections import OrderedDict
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile
from common import MARKERS,sha256,write_json
from audit_regenerated_boundaries import contour
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'external/Phenocycler_Analysis'))
from phenocycler.spatial import read_ome_region

MODULES={
 'identity':['DAPI','INS','GCG','SST','CHGA','CD45','PGP9.5','VIM','NaKATPase'],
 'immune':['DAPI','CD45','CD3e','CD4','CD8a','CD20','CD68','CD163','HLADR','CD35'],
 'structural':['DAPI','CD31','CD34','PDPN','SMA','ColIV','CK19','panCK','ECAD','BCatenin'],
 'states':['DAPI','B3TUBB','PGP9.5','GAP43','VIM','LGALS3','CD44','KRT14','NaKATPase','BActin']}


def select_candidates(table,n,selection):
    rows=[]
    for cluster,group in table.groupby('cluster_id',sort=True):
        if selection=='center_random':
            # Interleave central and random references; prefer donor diversity.
            a=group.loc[group.selection_reason.eq('cluster_center')]
            b=group.loc[group.selection_reason.eq('random')]
            group=pd.concat([a.iloc[[i]] if i<len(a) else a.iloc[:0] for i in range(min(2,len(a)))] + [b,a],ignore_index=True)
        else:group=group.loc[group.selection_reason.eq(selection)]
        group=group.drop_duplicates(['case_id','object_id'])
        chosen=[];donors=set()
        for r in group.to_dict('records'):
            if r['case_id'] not in donors:
                chosen.append(r);donors.add(r['case_id'])
            if len(chosen)>=n:break
        for r in group.to_dict('records'):
            if len(chosen)>=n:break
            if not any(r['object_id']==x['object_id'] and r['case_id']==x['case_id'] for x in chosen):chosen.append(r)
        rows.extend(chosen)
    return pd.DataFrame(rows)


class ImageReader:
    def __init__(self,image):
        self.image=Path(image);self.tf=tifffile.TiffFile(image)
        sidecar=self.image.with_suffix('').with_suffix('.ome.xml')
        xml=self.tf.ome_metadata
        if not xml:raise ValueError('Named OME metadata required')
        doc=ET.fromstring(xml);self.names=[e.attrib['Name'] for e in doc.iter() if e.tag.endswith('Channel')]
        if sidecar.is_file():
            side_names=[e.attrib['Name'] for e in ET.fromstring(sidecar.read_text()).iter() if e.tag.endswith('Channel')]
            if side_names!=self.names:raise ValueError('Sidecar channel labels disagree with authoritative embedded OME metadata')
        self.names=[{'PGP9_5':'PGP9.5','KRT8_18':'KRT8-18'}.get(m,m) for m in self.names]
        self.ranges={};self.tile_cache=OrderedDict()
        if len(self.names)!=self.tf.series[0].shape[0] or len(set(self.names))!=len(self.names):raise ValueError('OME channel identity mismatch')
    def display_range(self,marker):
        if marker not in self.ranges:
            ch=self.names.index(marker);level=self.tf.series[0].levels[-1]
            a=level.pages[ch].asarray(maxworkers=1)
            lo,hi=np.quantile(a,[.001,.9995]);self.ranges[marker]=[float(lo),float(max(hi,lo+1))]
        return self.ranges[marker]
    def crop(self,marker,bounds):
        x0,y0,x1,y1=bounds;ch=self.names.index(marker);page=self.tf.series[0].pages[ch]
        if self.tf.series[0].keyframe.is_tiled:return read_ome_region(self.tf,0,ch,y0,y1,x0,x1)
        # Converted candidates are tiled. A declared untiled input stays explicit.
        a=tifffile.memmap(self.image,mode='r');return np.asarray(a[ch,y0:y1,x0:x1])
    def close(self):self.tf.close()


def render(input_root,model_root,output_root,fit_ids,cluster_ids,modules,n,selection):
    input_root,model_root,output_root=map(Path,[input_root,model_root,output_root]);output_root.mkdir(parents=True,exist_ok=True)
    spec=json.loads((model_root/'specification.json').read_text())
    requests=[]
    for fit_id in fit_ids:
        path=model_root/fit_id/'reference_candidates.parquet'
        table=pd.read_parquet(path);table['case_id']=table.case_id.astype(str)
        if cluster_ids:table=table.loc[table.cluster_id.isin(cluster_ids)]
        chosen=select_candidates(table,n,selection)
        if len(chosen):requests.append(chosen)
    if not requests:raise ValueError('No requested cluster references')
    refs=pd.concat(requests,ignore_index=True)
    donors=json.loads((input_root/'preprocessing/inputs.json').read_text())['donors'];donors={r['case_id']:r for r in donors}
    objects={};readers={}
    try:
        for cid,group in refs.groupby('case_id'):
            cellfile=input_root/'corrected'/f'{cid}-cells.parquet';f=pd.read_parquet(cellfile,columns=['case_id','object_id','x_source_px','y_source_px','Cell Area']+MARKERS);wanted=set(group.object_id);f=f.loc[f.object_id.isin(wanted)].set_index('object_id')
            geomfile=Path(donors[cid]['geometry']);raw={str(g['id']):g for g in json.loads(geomfile.read_text())['features'] if f'{cid}:{g["id"]}' in wanted}
            readers[cid]=ImageReader(donors[cid]['image']);geometry_digest=sha256(geomfile)
            px=float(json.loads((input_root/'boundaries/sources'/cid/'source.json').read_text())['pixel_size_um'])
            for oid,r in f.iterrows():
                g=shape(raw[oid.split(':',1)[-1]]['geometry']);cx,cy=r.x_source_px,r.y_source_px;half=max(12/px,max(g.bounds[2]-g.bounds[0],g.bounds[3]-g.bounds[1])/2+4/px);height,width=readers[cid].tf.series[0].shape[-2:];bounds=[max(0,int(cx-half)),max(0,int(cy-half)),min(width,int(cx+half+1)),min(height,int(cy+half+1))];objects[(cid,oid)]={'geometry':g,'bounds':bounds,'means':r.to_dict(),'pixel_um':px,'geometry_sha256':geometry_digest}
        panels=[]
        for (fit_id,cluster),rows in refs.groupby(['fit_id','cluster_id'],sort=True):
            for module in modules:
                markers=MODULES[module];fig,axes=plt.subplots(len(rows),len(markers),figsize=(2*len(markers),2.3*len(rows)),squeeze=False);evidence=[]
                for j,r in enumerate(rows.itertuples()):
                    o=objects[(r.case_id,r.object_id)];reader=readers[r.case_id];x0,y0,x1,y1=o['bounds']
                    for k,m in enumerate(markers):
                        patch=reader.crop(m,o['bounds']);lo,hi=reader.display_range(m);ax=axes[j,k];ax.imshow(patch,cmap='gray',vmin=lo,vmax=hi,interpolation='nearest');contour(ax,o['geometry'],np.array([x0,y0]));ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([]);value=o['means'].get(m);ax.set_title(m+(f'\ncorrected MFI {value:.1f}' if value is not None else ''),fontsize=8)
                        if k==0:ax.set_ylabel(f'Image {r.case_id}\n{r.object_id[-12:]}\n{r.selection_reason}\n{o["means"]["Cell Area"]:.1f} µm²',fontsize=8)
                        ax.plot([3,3+5/o['pixel_um']],[y1-y0-4]*2,c='white',lw=2)
                    evidence.append({'case_id':r.case_id,'object_id':r.object_id,'source_pixel_bounds':o['bounds'],'geometry_sha256':o['geometry_sha256'],'selection_reason':r.selection_reason})
                fig.suptitle(f'{fit_id}, cluster {cluster}: {module} reference audit; no disease labels\nCyan is the native cell boundary; fixed image-wide channel ranges; 5 µm scale bars. Images are raw, printed cell MFI is REDSEA corrected.',fontsize=11);fig.tight_layout(rect=(0,0,1,.93));folder=output_root/fit_id;folder.mkdir(exist_ok=True);path=folder/f'cluster_{cluster:03d}_{module}_{selection.replace(":","_")}.png';fig.savefig(path,dpi=110);plt.close(fig);record={'fit_id':fit_id,'cluster_id':int(cluster),'module':module,'markers':markers,'model_identity':spec['identity'],'panel':str(path.resolve()),'panel_sha256':sha256(path),'references':evidence,'status':'rendered_awaiting_biological_review'};write_json(path.with_suffix('.json'),record);panels.append(record);print(json.dumps({'fit_id':fit_id,'cluster_id':int(cluster),'module':module,'panel':str(path)}),flush=True)
        write_json(output_root/'render_manifest.json',{'status':'rendered_awaiting_biological_review','model_identity':spec['identity'],'model_root':str(model_root.resolve()),'input_root':str(input_root.resolve()),'panels':panels,'display_ranges':{cid:r.ranges for cid,r in readers.items()},'biological_labels_assigned':False,'source_images_modified':False})
    finally:
        for reader in readers.values():reader.close()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True);p.add_argument('--model-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);p.add_argument('--fit-id',action='append',required=True);p.add_argument('--cluster-id',type=int,action='append',default=[]);p.add_argument('--module',choices=MODULES,action='append',default=[]);p.add_argument('--examples',type=int,default=4);p.add_argument('--selection',default='center_random');a=p.parse_args()
    if not 1<=a.examples<=12:raise ValueError('Review panels support 1–12 references per cluster')
    render(a.input_root,a.model_root,a.output_root,a.fit_id,a.cluster_id,a.module or ['identity'],a.examples,a.selection)
if __name__=='__main__':main()
