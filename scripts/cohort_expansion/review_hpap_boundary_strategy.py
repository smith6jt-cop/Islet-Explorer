#!/usr/bin/env python3
"""Compare boundary strategies over a deterministic HPAP sample, not one example.

No native objects or full-cohort geometry are written. The sample includes
flagged and unflagged masks across seed-count ranges plus the reported failure
and two large reference structures. No sampled mask is treated as ground truth.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage
import shapely
from shapely.geometry import shape
import tifffile

from audit_endocrine_boundaries import raw_crop
from endocrine_signal_support import hormone_support, signal_distances
from endocrine_region_envelopes import region_envelopes
from prepare_endocrine_candidates import digest,write_json


def select_examples(table, per_stratum=4):
    table = table.copy()
    table['seed_stratum'] = pd.cut(table.hormone_seed_cells,[0,2,9,49,np.inf],labels=['1-2','3-9','10-49','50+']).astype(str)
    table['selection_hash'] = table.candidate_id.map(lambda s:hashlib.sha256(('20260917:'+s).encode()).hexdigest())
    chosen = table.sort_values('selection_hash').groupby(['seed_stratum','low_signal_proximity'],sort=True).head(per_stratum).copy()
    chosen['selection_reason']='deterministic stratified sample'
    controls=table[table.candidate_id.eq('New_78_0242')].copy()
    controls['selection_reason']='owner-reported failure'
    large=table.sort_values(['area_um2','candidate_id'],ascending=[False,True]).head(2).copy()
    large['selection_reason']='large-structure reference, not ground truth'
    return pd.concat([controls,large,chosen]).drop_duplicates('candidate_id').reset_index(drop=True)


def geometry_lines(ax,geom,offset,color):
    for part in shapely.get_parts(geom):
        for ring in [part.exterior,*part.interiors]:
            xy=np.asarray(ring.coords)-offset
            ax.plot(xy[:,0],xy[:,1],color=color,linewidth=.6)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-root',type=Path,required=True)
    parser.add_argument('--support-audit',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args()
    root=args.candidate_root.resolve(strict=True)
    folder=root/'78'
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    spec=json.loads((folder/'prompts.json').read_text())
    verification=json.loads((root/'source-image-verification.json').read_text())
    source=next(d for d in verification['images'] if d['case_id']=='78')
    path=Path(spec['source_image']).resolve(strict=True)
    if str(path)!=source['path'] or path.stat().st_size!=source['bytes'] or path.stat().st_mtime_ns!=source['mtime_ns']:
        raise ValueError('Source fingerprint differs')
    geometry_path=folder/'review_candidates.geojson'
    if digest(geometry_path)!=json.loads((folder/'candidate_summary.json').read_text())['overlay_sha256']:
        raise ValueError('Candidate geometry differs')
    table=pd.read_csv(args.support_audit/'78/mask_support.csv')
    if len(table)!=437 or table.candidate_id.duplicated().any():
        raise ValueError('Expected the complete R5 HPAP candidate audit')
    selected=select_examples(table)
    selected.to_csv(out/'selected_examples.csv',index=False)
    features={f['properties']['name']:f for f in json.loads(geometry_path.read_text())['features']}
    pixels=tifffile.memmap(path,mode='r')
    seeds=pd.read_csv(folder/'seed_cells.csv')[['X_centroid','Y_centroid']].to_numpy()/spec['pixel_size_um']
    thresholds=np.array([spec['operational_thresholds'][m] for m in spec['rgb_channels']])
    display=json.loads((folder/'raw_channel_review.json').read_text())
    records=[];page_names=[]
    for start in range(0,len(selected),4):
        sample=selected.iloc[start:start+4]
        fig,axes=plt.subplots(len(sample),4,figsize=(13,3*len(sample)),squeeze=False)
        for row,example in enumerate(sample.itertuples()):
            geom=shape(features[example.candidate_id]['geometry'])
            raw,old_mask,bounds=raw_crop(pixels,spec,geom,pad_um=20)
            x0,y0,x1,y1=bounds
            nearby=seeds[((seeds>=[x0,y0])&(seeds<[x1,y1])).all(axis=1)]-[x0,y0]
            rgb=np.clip((raw-np.asarray(spec['display_min']))/(np.asarray(spec['display_max'])-spec['display_min']),0,1)
            # The proposed contours have no SAM geometry as an input.
            support=hormone_support(raw,thresholds,spec['pixel_size_um'])
            labels=region_envelopes(support,nearby,spec['pixel_size_um'])
            relevant=np.unique(labels[old_mask]);relevant=relevant[relevant>0]
            primary=np.isin(labels,relevant)
            # Only scope overlap metrics to this old location; every contextual
            # component is still displayed, including disconnected small foci.
            distances=signal_distances(support,spec['pixel_size_um'])
            base=primary.sum()
            variation={}
            for multiplier in [.5,2.]:
                alternative=region_envelopes(hormone_support(raw,thresholds*multiplier,spec['pixel_size_um']),nearby,spec['pixel_size_um'])
                ids=np.unique(alternative[old_mask]);ids=ids[ids>0]
                result=np.isin(alternative,ids)
                variation[f'pixel_cutoff_{multiplier:g}x_area_um2']=float(result.sum()*spec['pixel_size_um']**2)
            chga=pixels[spec['source_channel_names'].index('CHGA'),y0:y1,x0:x1]
            for col,title in enumerate(['Raw INS/GCG/SST','Withdrawn SAM outline','Experimental signal envelopes','CHGA context']):
                ax=axes[row,col]
                if col==3:ax.imshow(chga,cmap='gray',vmin=display['display_min']['CHGA'],vmax=display['display_max']['CHGA'])
                else:ax.imshow(rgb)
                if col==1:geometry_lines(ax,geom,np.array([x0,y0]),'yellow')
                if col==2 and labels.any() and not (labels>0).all():
                    ax.contour(labels>0,levels=[.5],colors=['cyan'],linewidths=.6)
                ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(title,fontsize=10)
                if col==0:ax.set_ylabel(f'{example.candidate_id}\n{example.hormone_seed_cells} hormone seeds\nPrior proximity flag: {example.low_signal_proximity}',fontsize=8)
            records.append(dict(candidate_id=example.candidate_id,selection_reason=example.selection_reason,
                prior_proximity_flag=bool(example.low_signal_proximity),source_pixel_bounds=bounds,
                seed_stratum=example.seed_stratum,raw_sam_area_um2=example.area_um2,
                new_envelope_area_um2=float(base*spec['pixel_size_um']**2),
                components_at_old_location=len(relevant),context_components=int(labels.max()),
                fraction_old_mask_covered=float((primary&old_mask).sum()/old_mask.sum()),
                envelope_within_3um_signal_fraction=float((distances[primary]<=3).mean()) if base else None,
                original_seed_count=example.hormone_seed_cells,**variation,
                review_decision='unreviewed',biological_reference_available=False))
        page=f'comparison_{start//4+1:02d}';page_names.append(page)
        fig.suptitle('HPAP-078: broader boundary-strategy review\nFixed per-image display ranges; experimental contours are not accepted islet boundaries.',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.95))
        for suffix in ['png','pdf']:fig.savefig(out/f'{page}.{suffix}',dpi=160)
        plt.close(fig)
        print(f'{page}: {len(records)}/{len(selected)} examples',flush=True)
    pd.DataFrame(records).to_csv(out/'comparison_metrics.csv',index=False)
    write_json(out/'pilot_manifest.json',dict(status='experimental_comparison_not_delivery',created_at_utc=datetime.now(timezone.utc).isoformat(),
        case_id='78',source_image_fingerprint=source,source_geometry_sha256=digest(geometry_path),
        input_support_audit_sha256=digest(args.support_audit/'78/mask_support.csv'),sample_count=len(selected),
        total_r5_masks=437,owner_reports_multiple_segmentation_errors=True,
        selection='Four per seed-count x previous-proximity-flag stratum by SHA256 identifier order, plus owner failure and two large references',
        pixel_cutoffs=spec['operational_thresholds'],pixel_cutoffs_validated=False,
        smoothing_sigma_um=.6,closing_radius_um=2,seed_proximity_um=4,
        interior_handling='Fill enclosed holes separately within each connected component; external gaps remain outside',
        native_objects_modified=False,new_delivery=False,accepted_islets=0,
        limitations=['The prior proximity screen is not a complete error detector.',
                     'Using cell-mean cutoffs on pixels is provisional; pixel threshold calibration is still required.',
                     'A filled signal envelope does not establish biological islet identity.',
                     'This candidate-based sample cannot measure missed islets elsewhere in the image.',
                     'Coordinates and cropping are taken from previous detections for comparison only; no full-image replacement has been built.'],
        code=[dict(path=str(p.resolve()),sha256=digest(p)) for p in [Path(__file__),Path(__file__).with_name('endocrine_region_envelopes.py')]]))
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>HPAP-078 broader boundary review</title>'
        '<style>body{font:16px system-ui;max-width:1400px;margin:2rem auto;padding:1rem}img{max-width:100%}</style>'
        '<h1>HPAP-078: broader boundary-strategy review</h1><p>R5 is withdrawn. These comparisons sample flagged and unflagged '
        'masks across seed counts. Cyan shows an experimental envelope from raw hormone foreground, without SAM. '
        'No replacement boundaries are accepted or written to QuPath. Pixel thresholds remain provisional.</p>'
        '<p><a href="selected_examples.csv">Sample selection</a> · <a href="comparison_metrics.csv">Metrics and review decisions</a></p>'+
        ''.join(f'<h2>Page {i+1}</h2><a href="{page}.pdf"><img loading="lazy" src="{page}.png"></a>' for i,page in enumerate(page_names))+'</html>')


if __name__=='__main__':main()
