#!/usr/bin/env python3
"""Render reference training examples and the fixed HPAP-078 transfer review."""
import argparse
import gzip
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import shapely
from shapely.geometry import shape
import tifffile

from audit_endocrine_boundaries import raw_crop
from prepare_endocrine_candidates import digest,write_json
from review_hpap_boundary_strategy import geometry_lines
from refine_endocrine_pixel_classifier import MARKERS


def save_figure(fig,out,name):
    for ext in ['png','pdf']:fig.savefig(out/f'{name}.{ext}',dpi=150)
    plt.close(fig)


def load_geometries(path,classification=None):
    fs=json.loads(path.read_text())['features']
    return [shape(f['geometry']) for f in fs if classification is None or f['properties'].get('classification',{}).get('name')==classification]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['samples-root','model-root','prediction-root','baseline-root','sam-root','selection','output-root']:
        p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    sampling=json.loads((args.samples_root/'sampling_manifest.json').read_text())
    trained=json.loads((args.model_root/'training_manifest.json').read_text())
    prediction=json.loads((args.prediction_root/'prediction_manifest.json').read_text())
    if prediction['case_id']!='78' or '78' in trained['training_donors']:raise ValueError('HPAP must be outside training')
    if digest(args.model_root/'endocrine_forest.json')!=prediction['model_sha256']:raise ValueError('Model provenance differs')
    source_pages=[]
    for r in sampling['donors']:
        cid=r['case_id'];rawdir=Path(r['raw_samples_root'])
        raw=np.load(rawdir/'endocrine_overview.npy')
        samples=np.load(rawdir/'samples.npz')
        old=np.asarray(Image.open(rawdir/'source_islet_regions.png'))>0
        supported=np.asarray(Image.open(rawdir/'supported_training_regions.png'))>0
        high=np.maximum(np.quantile(raw[old],[.995],axis=0).ravel(),1.)
        low=np.quantile(raw.reshape(-1,3),.05,axis=0)
        rgb=np.clip((raw-low)/np.maximum(high-low,1),0,1)
        indices=samples['linear_indices'][:6000];dominant=samples['dominant_marker']
        selected=[]
        for i,m in enumerate(MARKERS):
            choices=indices[dominant==i]
            if len(choices):
                linear=int(choices[len(choices)//2]);y,x=divmod(linear,int(samples['shape_yx'][1]))
                selected.append((x//2,y//2,f'{m}-dominant training sample'))
        yy,xx=np.nonzero(old & ~supported)
        if len(xx):selected.append((int(xx[len(xx)//2]),int(yy[len(yy)//2]),'Excluded reference interior'))
        fig,axes=plt.subplots(2,len(selected),figsize=(3*len(selected),6),squeeze=False)
        for col,(x,y,title) in enumerate(selected):
            x0=max(0,x-64);x1=min(raw.shape[1],x+65);y0=max(0,y-64);y1=min(raw.shape[0],y+65)
            for row in range(2):
                ax=axes[row,col];ax.imshow(rgb[y0:y1,x0:x1]);ax.set_axis_off()
                if row==0:ax.set_title(title,fontsize=9)
                else:
                    ref=old[y0:y1,x0:x1];sup=supported[y0:y1,x0:x1]
                    if ref.any() and not ref.all():ax.contour(ref,levels=[.5],colors=['yellow'],linewidths=.7)
                    if sup.any() and not sup.all():ax.contour(sup,levels=[.5],colors=['cyan'],linewidths=.6)
            axes[0,col].plot(x-x0,y-y0,'+',color='white',markersize=7)
        fig.suptitle(f'Original donor {cid}: reference-derived training locations\nINS red, GCG green, SST blue; yellow = old SAM islet; cyan = signal-supported interior\nThese are weak reference labels, not manually verified pixel ground truth.',fontsize=10)
        fig.tight_layout(rect=(0,0,1,.9));name=f'training_{cid}';save_figure(fig,out,name);source_pages.append(name)
    validation=pd.read_csv(args.model_root/'leave_one_donor_out.csv',dtype={'case_id':str})
    fig,axes=plt.subplots(1,2,figsize=(13,4))
    x=np.arange(len(validation))
    for col,prefix,title in [(0,'clear_reference_','Clear reference pixels'),(1,'all_reference_','All sampled reference interiors/background')]:
        axes[col].plot(x,validation[prefix+'sensitivity'],'o-',label='Interior / foreground agreement')
        axes[col].plot(x,validation[prefix+'specificity'],'s-',label='Background agreement')
        axes[col].set_xticks(x,validation.case_id,rotation=60);axes[col].set_ylim(0,1.02)
        axes[col].set_title(title);axes[col].legend(fontsize=8)
    fig.suptitle('Each donor excluded from fitting; agreement with source-derived labels, not biological accuracy')
    fig.tight_layout();save_figure(fig,out,'donor_held_out_reference_agreement')
    folder=args.sam_root/'78';spec=json.loads((folder/'prompts.json').read_text())
    pixels=tifffile.memmap(spec['source_image'],mode='r')
    old_by_name={f['properties']['name']:shape(f['geometry']) for f in json.loads((folder/'review_candidates.geojson').read_text())['features']}
    baseline=load_geometries(args.baseline_root/'raw_pixel_regions.geojson','Islet')
    raw_new=load_geometries(args.prediction_root/'raw_pixel_regions.geojson')
    boundaries=load_geometries(args.prediction_root/'review_candidates.geojson')
    variants=[baseline,raw_new,boundaries];trees=[shapely.STRtree(gs) for gs in variants]
    unions=[shapely.union_all(gs) for gs in variants]
    table=pd.read_csv(args.selection);pages=[];records=[]
    for start in range(0,len(table),3):
        chosen=table.iloc[start:start+3];fig,axes=plt.subplots(len(chosen),5,figsize=(16,3.2*len(chosen)),squeeze=False)
        for row,example in enumerate(chosen.itertuples()):
            old=old_by_name[example.candidate_id];raw,_,bounds=raw_crop(pixels,spec,old,pad_um=20)
            x0,y0,x1,y1=bounds;rgb=np.clip((raw-spec['display_min'])/(np.array(spec['display_max'])-spec['display_min']),0,1)
            for col,title in enumerate(['Raw hormones','Withdrawn SAM','Existing Jesus classifier','Original15 forest pixels','Forest + bounded cleanup']):
                ax=axes[row,col];ax.imshow(rgb)
                if col==1:geometry_lines(ax,old,np.array([x0,y0]),'yellow')
                if col>=2:
                    idx=col-2
                    for i in trees[idx].query(shapely.box(*bounds),predicate='intersects'):
                        geometry_lines(ax,variants[idx][i],np.array([x0,y0]),'cyan')
                ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(title,fontsize=10)
                if col==0:ax.set_ylabel(example.candidate_id,fontsize=8)
            records.append(dict(candidate_id=example.candidate_id,old_area_um2=old.area*spec['pixel_size_um']**2,
                existing_classifier_overlap_fraction=old.intersection(unions[0]).area/old.area,
                refined_pixel_overlap_fraction=old.intersection(unions[1]).area/old.area,
                refined_boundary_overlap_fraction=old.intersection(unions[2]).area/old.area,
                boundary_components_touching_old=len(trees[2].query(old,predicate='intersects')),review_decision='unreviewed'))
        fig.suptitle('HPAP-078 excluded from forest training — fixed 30-location transfer review\nCyan regions need review; no accepted geometry, phenotype assignment or downstream analysis.',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.93));name=f'hpap_comparison_{start//3+1:02d}';save_figure(fig,out,name);pages.append(name)
        print(f'HPAP comparison {start+len(chosen)}/{len(table)}',flush=True)
    pd.DataFrame(records).to_csv(out/'hpap_comparison_metrics.csv',index=False)
    with gzip.open(spec['source_geometry'],'rt') as stream:
        old={f['properties'].get('name'):shape(f['geometry']) for f in json.load(stream)['features'] if f['properties'].get('classification',{}).get('name')=='Islet'}
    controls=[]
    for name in ['Islet_130','Islet_206','Islet_1148','Islet_124','Islet_161','Islet_179']:
        controls.append(dict(previous_islet=name,refined_pixel_fraction=old[name].intersection(unions[1]).area/old[name].area,
                             refined_boundary_fraction=old[name].intersection(unions[2]).area/old[name].area))
    pd.DataFrame(controls).to_csv(out/'hpap_previous_negative_controls.csv',index=False)
    # Full field catches islands/halos away from the old candidate locations.
    raw=np.stack([pixels[spec['source_channel_names'].index(m),::8,::8] for m in MARKERS],-1)
    rgb=np.clip((raw-spec['display_min'])/(np.array(spec['display_max'])-spec['display_min']),0,1)
    fig,axes=plt.subplots(1,3,figsize=(16,6));axes[0].imshow(rgb);axes[0].set_title('Full HPAP-078 hormones')
    for ax,name,title in zip(axes[1:],['pixel_classes','boundary_regions'],['Reference-trained foreground','Boundaries after cleanup']):
        ax.imshow(np.asarray(Image.open(args.prediction_root/f'{name}.png'))[::4,::4],cmap='gray',vmin=0,vmax=1);ax.set_title(title)
    for ax in axes:ax.set_axis_off()
    fig.tight_layout();save_figure(fig,out,'hpap_whole_image')
    summary=dict(status='reference_trained_classifier_transfer_review',training_donors=trained['training_donors'],
        training_pixels=trained['training_pixels'],HPAP_in_training=False,
        original_islet_annotations=sum(json.loads(Path(r['source_reference_manifest']).read_text())['islet_annotations'] for r in sampling['donors']),
        held_out_clear_reference_balanced_agreement_median=float(validation.clear_reference_balanced_accuracy.median()),
        held_out_all_reference_sensitivity_median=float(validation.all_reference_sensitivity.median()),
        held_out_all_reference_specificity_median=float(validation.all_reference_specificity.median()),
        fixed_hpap_comparison_locations=len(table),negative_controls=controls,raw_foreground_regions=len(raw_new),
        boundary_regions=len(boundaries),accepted_islets=0,new_boundary_delivery=False,
        validation_interpretation=trained['validation_interpretation'],model_sha256=prediction['model_sha256'])
    write_json(out/'review_manifest.json',summary)
    def page(name):return f'<a href="{name}.pdf"><img loading="lazy" src="{name}.png"></a>'
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Original15 endocrine pixel refinement</title>'
        '<style>body{font:16px system-ui;max-width:1600px;margin:2rem auto;padding:1rem}img{max-width:100%}</style>'
        '<h1>Endocrine pixel classifier trained from original 15 donors</h1>'
        '<p>INS/GCG/SST only; independent per-image tissue normalization. Source islets locate signal-supported training pixels. '
        'All 5,382 source outlines carry SAM provenance, so their edges and dark interiors are not foreground training labels. '
        'HPAP-078 and the other five new donors were excluded from training. No biological accuracy or analytical acceptance is claimed.</p>'
        '<p><a href="review_manifest.json">Review record</a> · <a href="hpap_comparison_metrics.csv">Fixed comparison metrics</a> · '
        '<a href="hpap_previous_negative_controls.csv">Previous negative controls</a></p>'
        '<h2>Donor-held-out reference agreement</h2>'+page('donor_held_out_reference_agreement')+
        '<h2>HPAP-078 full image</h2>'+page('hpap_whole_image')+
        ''.join(f'<h2>HPAP comparison {i+1}</h2>'+page(name) for i,name in enumerate(pages))+
        '<h2>Training reference examples</h2>'+''.join(page(name) for name in source_pages)+'</html>')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
