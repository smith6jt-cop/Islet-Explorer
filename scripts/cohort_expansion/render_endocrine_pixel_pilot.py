#!/usr/bin/env python3
"""Compare an actual QuPath endocrine-only pixel classifier with withdrawn SAM."""
import argparse
import gzip
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


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prediction-root',type=Path,required=True)
    parser.add_argument('--candidate-root',type=Path,required=True)
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args()
    pred=args.prediction_root.resolve(strict=True)
    manifest=json.loads((pred/'prediction_manifest.json').read_text())
    if manifest['case_id']!='78' or manifest['SAM_used']:
        raise ValueError('Expected the independent HPAP pixel-classifier pilot')
    for name,key in [('raw_pixel_regions.geojson','raw_region_sha256'),('pixel_classes.png','pixel_raster_sha256'),('endocrine_classifier.json','classifier_sha256')]:
        if digest(pred/name)!=manifest[key]:raise ValueError('Pixel output changed')
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    folder=args.candidate_root.resolve(strict=True)/'78'
    spec=json.loads((folder/'prompts.json').read_text())
    if digest(folder/'prompts.json')!=manifest['source_specification_sha256']:
        raise ValueError('Source specification differs')
    table=pd.read_csv(args.selection)
    features={f['properties']['name']:f for f in json.loads((folder/'review_candidates.geojson').read_text())['features']}
    native=json.loads((pred/'raw_pixel_regions.geojson').read_text())['features']
    foreground=[shape(f['geometry']) for f in native if f['properties'].get('classification',{}).get('name')=='Islet']
    if not foreground:raise ValueError('No endocrine pixel regions')
    tree=shapely.STRtree(foreground)
    union=shapely.union_all(foreground)
    tissue=shape(json.loads((folder/'valid_tissue.geojson').read_text()))
    if union.difference(tissue).area>1e-4:raise ValueError('Pixel objects extend outside valid tissue')
    pixels=tifffile.memmap(spec['source_image'],mode='r')
    display=json.loads((folder/'raw_channel_review.json').read_text())
    records=[];pages=[]
    for start in range(0,len(table),4):
        sample=table.iloc[start:start+4]
        fig,axes=plt.subplots(len(sample),4,figsize=(13,3*len(sample)),squeeze=False)
        for row,example in enumerate(sample.itertuples()):
            old=shape(features[example.candidate_id]['geometry'])
            raw,_,bounds=raw_crop(pixels,spec,old,pad_um=20)
            x0,y0,x1,y1=bounds
            rgb=np.clip((raw-np.asarray(spec['display_min']))/(np.asarray(spec['display_max'])-spec['display_min']),0,1)
            chga=pixels[spec['source_channel_names'].index('CHGA'),y0:y1,x0:x1]
            relevant=tree.query(shapely.box(*bounds),predicate='intersects')
            for col,title in enumerate(['Raw INS/GCG/SST','Withdrawn SAM boundary','Endocrine-only pixel classifier','CHGA context']):
                ax=axes[row,col]
                if col==3:ax.imshow(chga,cmap='gray',vmin=display['display_min']['CHGA'],vmax=display['display_max']['CHGA'])
                else:ax.imshow(rgb)
                if col==1:geometry_lines(ax,old,np.array([x0,y0]),'yellow')
                if col==2:
                    for i in relevant:geometry_lines(ax,foreground[i],np.array([x0,y0]),'cyan')
                ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(title,fontsize=10)
                if col==0:ax.set_ylabel(f'{example.candidate_id}\n{example.hormone_seed_cells} prior hormone seeds\nPrior proximity flag: {example.low_signal_proximity}',fontsize=8)
            records.append(dict(candidate_id=example.candidate_id,source_pixel_bounds=bounds,
                prior_proximity_flag=bool(example.low_signal_proximity),old_area_um2=example.area_um2,
                pixel_foreground_fraction_within_old_mask=old.intersection(union).area/old.area,
                pixel_components_touching_old_mask=len(tree.query(old,predicate='intersects')),
                review_decision='unreviewed'))
        name=f'comparison_{start//4+1:02d}';pages.append(name)
        fig.suptitle('HPAP-078: existing INS/GCG/SST pixel classifier, unchanged\nCyan = native classifier foreground contours; no SAM or cell gates; no postprocessing or biological acceptance.',fontsize=10)
        fig.tight_layout(rect=(0,0,1,.95))
        for ext in ['png','pdf']:fig.savefig(out/f'{name}.{ext}',dpi=160)
        plt.close(fig)
        print(f'{name}: {len(records)}/{len(table)} examples',flush=True)
    pd.DataFrame(records).to_csv(out/'comparison_metrics.csv',index=False)
    with gzip.open(spec['source_geometry'],'rt') as stream:
        previous={f['properties'].get('name'):shape(f['geometry']) for f in json.load(stream)['features'] if f['properties'].get('classification',{}).get('name')=='Islet'}
    controls=[]
    for name in ['Islet_130','Islet_206','Islet_1148','Islet_124','Islet_161','Islet_179']:
        old=previous[name]
        controls.append(dict(previous_islet=name,pixel_foreground_fraction=old.intersection(union).area/old.area))
    pd.DataFrame(controls).to_csv(out/'old_negative_controls.csv',index=False)
    # Whole-image coverage provides context beyond the candidate-based sample.
    step=8
    raw=np.stack([pixels[spec['source_channel_names'].index(m),::step,::step] for m in spec['rgb_channels']],-1)
    rgb=np.clip((raw-np.asarray(spec['display_min']))/(np.asarray(spec['display_max'])-spec['display_min']),0,1)
    classes=np.asarray(Image.open(pred/'pixel_classes.png'))
    if classes.ndim!=2 or not set(np.unique(classes)).issubset({0,1}):raise ValueError('Expected class-index raster')
    fig,axes=plt.subplots(1,2,figsize=(14,7))
    axes[0].imshow(rgb);axes[0].set_title('Full HPAP-078 hormone image')
    axes[1].imshow(classes,cmap='gray',vmin=0,vmax=1);axes[1].set_title('Full-image classifier foreground (including outside tissue)')
    for ax in axes:ax.set_axis_off()
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(out/f'whole_image.{ext}',dpi=150)
    plt.close(fig)
    areas=np.array([g.area*spec['pixel_size_um']**2 for g in foreground])
    summary=dict(status='pixel_classifier_comparison_not_delivery',case_id='78',sample_count=len(records),
        full_image_classified=True,raw_foreground_regions=len(foreground),
        raw_regions_below_80um2=int((areas<80).sum()),raw_regions_at_least_80um2=int((areas>=80).sum()),
        total_foreground_area_um2=float(areas.sum()),all_native_regions_valid=bool(shapely.is_valid(foreground).all()),
        native_regions_clipped_to_valid_tissue=True,native_coordinate_transform_used=True,
        raw_source_cell_objects_modified=False,new_boundary_delivery=False,accepted_islets=0,
        classifier_sha256=manifest['classifier_sha256'],classifier_training_provenance='Unknown',
        negative_control_comparisons=controls,
        comparison_code_sha256=digest(Path(__file__)))
    write_json(out/'comparison_manifest.json',summary)
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>HPAP endocrine pixel classifier</title>'
        '<style>body{font:16px system-ui;max-width:1400px;margin:2rem auto;padding:1rem}img{max-width:100%}</style>'
        '<h1>HPAP-078 endocrine-only pixel classifier</h1><p>Existing Jesus.json, using only INS/GCG/SST, tested unchanged on the whole image. '
        'Cyan contours are native foreground regions. They are not accepted whole-islet boundaries. The classifier training history is unknown. '
        'No SAM refinement, cell gates, new phenotype assignments or native-project edits were applied.</p>'
        '<p><a href="../endocrine_classifier.json">QuPath classifier JSON</a> · <a href="comparison_metrics.csv">Comparison metrics</a></p>'
        '<h2>Full image</h2><a href="whole_image.pdf"><img src="whole_image.png"></a>'+
        ''.join(f'<h2>Comparison {i+1}</h2><a href="{page}.pdf"><img loading="lazy" src="{page}.png"></a>' for i,page in enumerate(pages))+'</html>')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
