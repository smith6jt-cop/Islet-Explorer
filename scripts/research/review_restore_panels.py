#!/usr/bin/env python3
"""Compact donor-by-marker source-image audit pages for fresh RESTORE calls."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import tifffile
from shapely.geometry import shape
from common import sha256,write_json
from restore_positivity import pair_name,lineage
from render_restore_positivity import select_examples
from audit_regenerated_boundaries import contour

IMMUNE=['CD45','CD3e','CD8a','CD4','CD20','CD68','CD163','HLADR','CD35','CD44','LGALS3']


def donor_panels(root,screen,out,cid,markers):
    folder=screen/'donors'/cid;baseline=dict(lineage.FIRST_PASS_PAIRS+lineage.LEVEL2_PAIRS)
    sample=pd.read_parquet(folder/'sample.parquet');cohort=pd.read_csv(root/'cohort.csv',dtype={'case_id':str}).set_index('case_id')
    acquisition_path=Path(cohort.loc[cid,'boundary_source']);acquisition=json.loads(acquisition_path.read_text())
    acquisition['source_channel_names']=[m.replace('PGP9_5','PGP9.5') for m in acquisition['source_channel_names']]
    definitions={};needed=set();selected_keys=set()
    for marker in markers:
        reference=baseline[marker];name=pair_name(marker,reference);pair_path=folder/(name+'.json')
        record=json.loads(pair_path.read_text());metrics=record['metrics'];threshold=metrics.get('candidate_full_maximum')
        selected=[]
        if metrics['state'] in ['VALID','PAIR_REVIEW_PENDING'] and threshold is not None and threshold>0:
            labels=pd.read_parquet(folder/(name+'_sample_labels.parquet'))
            selected=[r for r in select_examples(sample,labels,marker,threshold) if r['selection'] in ['brightest reference control','just above divisor','high target']]
        definitions[marker]=dict(reference=reference,name=name,pair_path=pair_path,record=record,metrics=metrics,threshold=threshold,selected=selected)
        for r in selected:
            selected_keys.add(r['object_id'].removeprefix(cid+':'))
            for channel in ['DAPI',marker,reference]:needed.add(channel)
    inputs=json.loads((root/'preprocessing/inputs.json').read_text())
    geometry=Path(next(r['geometry'] for r in inputs['donors'] if r['case_id']==cid))
    cells={str(f['id']):shape(f['geometry']) for f in json.loads(geometry.read_text())['features'] if str(f['id']) in selected_keys}
    if set(cells)!=selected_keys:raise ValueError('Missing exact native cell geometry')
    ranges={};patches={};offsets={};half=14/float(acquisition['pixel_size_um'])
    # Decode each required channel once per donor. Memmapped raw sources avoid
    # copies; compressed native TIFFs are read one channel at a time (bounded RAM).
    with tifffile.TiffFile(acquisition['source_image']) as image:
        try:mapped=tifffile.memmap(acquisition['source_image'],mode='r')
        except ValueError:mapped=None
        for channel in sorted(needed):
            index=acquisition['source_channel_names'].index(channel)
            plane=mapped[index] if mapped is not None else image.series[0].pages[index].asarray(maxworkers=1)
            ranges[channel]=np.quantile(plane[::16,::16],[.001,.9995])
            for marker,definition in definitions.items():
                if channel not in ['DAPI',marker,definition['reference']]:continue
                for example in definition['selected']:
                    key=example['object_id'];cx,cy=example['x_source_px'],example['y_source_px']
                    x0,x1=max(0,int(cx-half)),min(plane.shape[1],int(cx+half+1));y0,y1=max(0,int(cy-half)),min(plane.shape[0],int(cy+half+1))
                    patches[(key,channel)]=np.array(plane[y0:y1,x0:x1]);offsets[key]=np.array([x0,y0])
            del plane
    destination=out/'donors'/cid;destination.mkdir(parents=True,exist_ok=True);records=[]
    for marker,definition in definitions.items():
        selected=definition['selected'];threshold=definition['threshold'];reference=definition['reference']
        fig,axes=plt.subplots(1,9,figsize=(13.5,2.15),squeeze=False);axes=axes[0]
        for ax in axes:ax.set_xticks([]);ax.set_yticks([])
        if not selected:
            for ax in axes:ax.set_axis_off()
            axes[4].text(.5,.5,f'Donor {cid} / {marker}: unavailable\n{definition["metrics"]["state"]}',ha='center',va='center')
        for j,example in enumerate(selected):
            for k,channel in enumerate(['DAPI',marker,reference]):
                ax=axes[j*3+k];key=example['object_id'];levels=ranges[channel]
                ax.imshow(patches[(key,channel)],cmap='gray',vmin=levels[0],vmax=levels[1],interpolation='nearest')
                contour(ax,cells[key.removeprefix(cid+':')],offsets[key],color='cyan')
                ax.set_title(channel,fontsize=8)
                if k==1:ax.set_xlabel(f'{example["selection"]}\n{example["redsea__"+marker]:.0f}',fontsize=7)
        threshold_text=f'{threshold:.1f}' if threshold is not None else 'unavailable'
        fig.suptitle(f'Donor {cid} | {marker} ← {reference} | fresh divisor {threshold_text} | current image review pending',fontsize=9)
        fig.tight_layout(rect=(0,0,1,.87),w_pad=.05);path=destination/(marker+'.png');fig.savefig(path,dpi=110);plt.close(fig)
        evidence=dict(case_id=cid,marker=marker,reference=reference,current_review='pending',
            computational_state=definition['metrics']['state'],threshold=threshold,
            evaluation_path=str(definition['pair_path'].resolve()),evaluation_sha256=sha256(definition['pair_path']),
            source_acquisition_sha256=sha256(acquisition_path),native_source_sha256=definition['record']['fingerprint']['native_sha256'],
            corrected_sha256=definition['record']['fingerprint']['corrected_sha256'],panel_path=str(path.resolve()),panel_sha256=sha256(path),
            examples=[{k:r[k] for k in ['object_id','selection','x_source_px','y_source_px','Cell Area','redsea__'+marker,'redsea__'+reference]} for r in selected],
            channel_ranges={c:ranges[c].tolist() for c in ['DAPI',marker,reference] if c in ranges})
        write_json(destination/(marker+'.json'),evidence);records.append(evidence)
    print(cid,len(records),'audit strips',flush=True);return records


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--screen-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--marker',action='append');p.add_argument('--donor',action='append')
    a=p.parse_args();donors=sorted(pd.read_csv(a.input_root/'cohort.csv',dtype=str).case_id,key=int)
    if a.donor:donors=[c for c in donors if c in a.donor]
    markers=a.marker or IMMUNE;a.output_root.mkdir(parents=True,exist_ok=True);records=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(donor_panels,a.input_root,a.screen_root,a.output_root,cid,markers) for cid in donors]
        for future in as_completed(futures):records.extend(future.result())
    pages=[]
    for marker in markers:
        for start in range(0,len(donors),7):
            selected=donors[start:start+7];images=[Image.open(a.output_root/'donors'/cid/(marker+'.png')).convert('RGB') for cid in selected]
            canvas=Image.new('RGB',(max(i.width for i in images),sum(i.height for i in images)),color='white');y=0
            for picture in images:canvas.paste(picture,(0,y));y+=picture.height
            target=a.output_root/f'{marker}_page{start//7+1}.png';canvas.save(target)
            pages.append(dict(marker=marker,donors=selected,path=str(target.resolve()),sha256=sha256(target)))
    write_json(a.output_root/'manifest.json',dict(status='pending_actual_image_inspection',donors=donors,markers=markers,
        panel_count=len(records),pages=pages,runner_sha256=sha256(__file__),
        interpretation='Actual source-cell evidence. Computational separation and these rendered panels do not themselves constitute biological acceptance.'))


if __name__=='__main__':main()
