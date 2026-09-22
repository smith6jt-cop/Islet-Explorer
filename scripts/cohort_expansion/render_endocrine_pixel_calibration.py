#!/usr/bin/env python3
"""Render source operating-point tradeoffs and four diagnostic HPAP comparisons."""
import argparse
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

from audit_endocrine_boundaries import raw_crop
from review_hpap_boundary_strategy import geometry_lines


def save(fig,out,name):
    for suffix in ['png','pdf']:fig.savefig(out/f'{name}.{suffix}',dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['calibration-root','prediction-root','sensitivity-root','sam-root','output-root']:
        parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args();out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    table=pd.read_csv(args.calibration_root/'operating_point_summary.csv')
    result=json.loads((args.calibration_root/'calibration_result.json').read_text())
    target_percent=100*result['target_mean_background_reference_FPR']
    outcome=(f'no tested cutoff met the {target_percent:g}% background-reference target'
             if result['selected_minimum_votes'] is None else
             f'source-reference target first met at {result["selected_minimum_votes"]} votes')
    unanimous=table.iloc[-1]
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].plot(table.cutoff,100*table.mean_background_FPR,label='All background references')
    axes[0].plot(table.cutoff,100*table.mean_clear_background_FPR,label='Clear background references')
    axes[0].axhline(target_percent,color='red',linestyle='--',label=f'{target_percent:g}% reference target')
    axes[0].set_ylabel('Mean background-reference positive fraction (%)');axes[0].legend(fontsize=8)
    axes[1].plot(table.cutoff,100*table.median_foreground_recall,label='Median across donors')
    axes[1].plot(table.cutoff,100*table.min_foreground_recall,label='Lowest donor')
    axes[1].set_ylabel('Clear foreground-reference recall (%)');axes[1].legend(fontsize=8)
    for ax in axes:ax.set_xlabel('Minimum forest vote fraction');ax.grid(alpha=.2)
    fig.suptitle(f'Source-donor calibration: {outcome}\nWeak reference labels; these curves are not manual boundary accuracy estimates.',fontsize=11)
    fig.tight_layout();save(fig,out,'source_operating_point_tradeoff')
    spec=json.loads((args.sam_root/'78/prompts.json').read_text())
    pixels=tifffile.memmap(spec['source_image'],mode='r')
    previous={f['properties']['name']:shape(f['geometry']) for f in json.loads((args.sam_root/'78/review_candidates.geojson').read_text())['features']}
    variants=[[shape(f['geometry']) for f in json.loads((root/'review_candidates.geojson').read_text())['features']]
              for root in [args.prediction_root,args.sensitivity_root]]
    trees=[shapely.STRtree(g) for g in variants]
    examples=['New_78_0242','New_78_0215','New_78_0116','New_78_0333']
    fig,axes=plt.subplots(4,3,figsize=(12,12))
    for row,name in enumerate(examples):
        raw,_,(x0,y0,x1,y1)=raw_crop(pixels,spec,previous[name],pad_um=20)
        rgb=np.clip((raw-spec['display_min'])/(np.array(spec['display_max'])-spec['display_min']),0,1)
        for col in range(3):
            ax=axes[row,col];ax.imshow(rgb)
            if col:
                for i in trees[col-1].query(shapely.box(x0,y0,x1,y1),predicate='intersects'):
                    geometry_lines(ax,variants[col-1][i],np.array([x0,y0]),'cyan')
            ax.set_xlim(0,x1-x0);ax.set_ylim(y1-y0,0);ax.set_xticks([]);ax.set_yticks([])
            if row==0:ax.set_title(['Raw INS/GCG/SST','80/160 trees','160/160 trees'][col])
        axes[row,0].set_ylabel(name)
    fig.suptitle('Diagnostic operating-point comparison; neither version is accepted replacement geometry',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.97));save(fig,out,'operating_point_examples')
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Endocrine classifier operating points</title>'
        '<style>body{font:16px system-ui;max-width:1300px;margin:2rem auto;padding:1rem}img{max-width:100%}</style>'
        '<h1>Conservative operating-point sensitivity</h1><p>The 160/160-tree version is a diagnostic sensitivity, '
        f'not an accepted primary cutoff. Calibration outcome: {outcome}. '
        f'The lowest observed background-reference fraction was {100*result["minimum_achieved_FPR"]:.3f}%; '
        f'median clear-foreground recall at unanimity was {100*unanimous.median_foreground_recall:.1f}%. '
        'These figures concern weak reference labels, not independent biological accuracy.</p>'
        '<p>The four image examples illustrate the reported failure, weak staining and two large structures. '
        'They were selected for diagnostic comparison, not statistical evaluation.</p>'+
        ''.join(f'<a href="{name}.pdf"><img src="{name}.png"></a>' for name in ['source_operating_point_tradeoff','operating_point_examples'])+'</html>')


if __name__=='__main__':main()
