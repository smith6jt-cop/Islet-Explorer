#!/usr/bin/env python3
"""Render current donor-specific RESTORE controls and exact source-cell images.

These are review evidence, not automatically accepted gates. The biological
reference pair and computational result remain separate from image inspection.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tifffile
from shapely.geometry import shape
from common import sha256,write_json
from restore_positivity import pair_name,lineage
from audit_regenerated_boundaries import contour


def select_examples(sample,labels,target,threshold):
    frame=sample.merge(labels,on=['object_id','x_source_px','y_source_px'],validate='one_to_one')
    column='redsea__'+target
    requests=[('brightest reference control',frame.loc[frame.full_reference_control].sort_values(column,ascending=False)),
              ('target anchor',frame.loc[frame.full_target_population].sort_values(column,ascending=False)),
              ('just below divisor',frame.loc[frame[column]<=threshold].sort_values(column,ascending=False)),
              ('just above divisor',frame.loc[frame[column]>threshold].sort_values(column)),
              ('low target',frame.iloc[(frame[column]-threshold*.1).abs().argsort()]),
              ('high target',frame.iloc[(frame[column]-threshold*3).abs().argsort()])]
    used=set();selected=[]
    for label,pool in requests:
        available=pool.loc[~pool.object_id.isin(used)]
        if available.empty:continue
        record=available.iloc[0].to_dict();record['selection']=label;selected.append(record);used.add(record['object_id'])
    return selected


def render(root,screen,cid,target,reference):
    folder=screen/'donors'/cid;name=pair_name(target,reference);pair_path=folder/(name+'.json')
    record=json.loads(pair_path.read_text());metrics=record['metrics'];threshold=metrics.get('candidate_full_maximum')
    if threshold is None or threshold<=0:raise ValueError('No positive candidate control maximum to review')
    sample=pd.read_parquet(folder/'sample.parquet');labels=pd.read_parquet(folder/(name+'_sample_labels.parquet'))
    selected=select_examples(sample,labels,target,threshold)
    if not selected:raise ValueError('No assessable source-cell examples')
    cohort=pd.read_csv(root/'cohort.csv',dtype={'case_id':str}).set_index('case_id')
    acquisition_path=Path(cohort.loc[cid,'boundary_source']);acquisition=json.loads(acquisition_path.read_text())
    channels=['DAPI',target,reference,'NaKATPase'];names=[m.replace('PGP9_5','PGP9.5') for m in acquisition['source_channel_names']];indices=[names.index(m) for m in channels]
    try:pixels=tifffile.memmap(acquisition['source_image'],mode='r')
    except ValueError:pixels=tifffile.imread(acquisition['source_image'],key=indices,maxworkers=1);indices=list(range(len(channels)))
    if pixels.ndim!=3:raise ValueError('Source image must retain explicit channel,Y,X axes')
    inputs=json.loads((root/'preprocessing/inputs.json').read_text())
    geometry=Path(next(r['geometry'] for r in inputs['donors'] if r['case_id']==cid))
    keys={r['object_id'].removeprefix(cid+':') for r in selected}
    cells={str(f['id']):shape(f['geometry']) for f in json.loads(geometry.read_text())['features'] if str(f['id']) in keys}
    if set(cells)!=keys:raise ValueError('Image evidence lacks exact source UUID cell boundaries')
    ranges=[np.quantile(pixels[index,::16,::16],[.001,.9995]) for index in indices]
    output=screen/'review'/cid;output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(len(selected),len(channels),figsize=(12,2.8*len(selected)),squeeze=False)
    for row,example in enumerate(selected):
        half=14/float(acquisition['pixel_size_um']);cx,cy=example['x_source_px'],example['y_source_px']
        x0,x1=max(0,int(cx-half)),min(pixels.shape[2],int(cx+half+1));y0,y1=max(0,int(cy-half)),min(pixels.shape[1],int(cy+half+1))
        cell=cells[example['object_id'].removeprefix(cid+':')]
        for col,(marker,index,levels) in enumerate(zip(channels,indices,ranges)):
            ax=axes[row,col];ax.imshow(pixels[index,y0:y1,x0:x1],cmap='gray',vmin=levels[0],vmax=levels[1],interpolation='nearest')
            contour(ax,cell,np.array([x0,y0]),color='cyan');ax.set_xticks([]);ax.set_yticks([])
            if row==0:ax.set_title(marker)
            if col==0:ax.set_ylabel(f'{example["selection"]}\n{target}: {example["redsea__"+target]:.1f}\n{reference}: {example["redsea__"+reference]:.1f}',fontsize=8)
    fig.suptitle(f'Donor {cid}: {target} ← {reference}; fresh divisor {threshold:.2f}\nNative source intensities, fixed per-channel display range, exact UUID boundaries; review pending',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.955));path=output/(name+'_cells.png');fig.savefig(path,dpi=140);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5));joined=sample.merge(labels,on=['object_id','x_source_px','y_source_px'],validate='one_to_one')
    colors=np.full(len(joined),'#bbbbbb',object);colors[joined.full_reference_control]='#2b83ba';colors[joined.full_target_population]='#d95f02'
    axes[0].scatter(joined['redsea__'+reference],joined['redsea__'+target],c=colors,s=4,alpha=.5,rasterized=True)
    axes[0].axhline(threshold,color='red',linewidth=1);separator=metrics.get('reference_separator')
    if separator is not None:axes[0].axvline(separator,color='black',linewidth=1)
    axes[0].set(xscale='symlog',yscale='symlog',xlabel=reference+' corrected mean',ylabel=target+' corrected mean',title='Reference controls blue; target anchors orange')
    values=joined['redsea__'+target]
    axes[1].hist(np.log1p(values),bins=100,color='#777777');axes[1].axvline(np.log1p(threshold),color='red')
    axes[1].set(xlabel='log1p corrected '+target,ylabel='sample cells',title=f'Computational state: {metrics["state"]}')
    fig.tight_layout();diagnostic=output/(name+'_diagnostics.png');fig.savefig(diagnostic,dpi=150);plt.close(fig)
    columns=['object_id','selection','x_source_px','y_source_px','Cell Area','redsea__'+target,'redsea__'+reference]
    pd.DataFrame(selected)[columns].to_csv(output/(name+'_examples.csv'),index=False)
    write_json(output/(name+'.json'),dict(case_id=cid,target=target,reference=reference,threshold=threshold,
        state='awaiting_image_review',evaluation_path=str(pair_path.resolve()),evaluation_sha256=sha256(pair_path),
        source_acquisition_path=str(acquisition_path.resolve()),source_acquisition_sha256=sha256(acquisition_path),
        corrected_sha256=record['fingerprint']['corrected_sha256'],native_source_sha256=record['fingerprint']['native_sha256'],
        panel_path=str(path.resolve()),panel_sha256=sha256(path),diagnostic_path=str(diagnostic.resolve()),diagnostic_sha256=sha256(diagnostic),
        channel_display_ranges={m:levels.tolist() for m,levels in zip(channels,ranges)},examples=len(selected)))
    print(cid,target,path,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--screen-root',type=Path,required=True);parser.add_argument('--donor',action='append',required=True)
    parser.add_argument('--marker',action='append',required=True);parser.add_argument('--reference')
    args=parser.parse_args();baseline=dict(lineage.FIRST_PASS_PAIRS+lineage.LEVEL2_PAIRS)
    for cid in args.donor:
        for marker in args.marker:
            render(args.input_root,args.screen_root,cid,marker,args.reference or baseline[marker])


if __name__=='__main__':main()
