#!/usr/bin/env python3
"""Assess source-donor operating points and save an explicitly unselected sensitivity.

No HPAP labels or disease variables enter cutoff assessment. Failure to meet the
reference-background target remains a failure, even if a sensitivity looks better.
"""
import argparse
import json
from pathlib import Path
import time

import cv2
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
from shapely.geometry import shape

from prepare_endocrine_candidates import digest,write_json
from refine_endocrine_pixel_classifier import (clean_regions,probability,rasterize,
                                               refuse_invalid_samples,validate_training_donors)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['samples-root','model-root','prediction-root','output-root','sensitivity-output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();refuse_invalid_samples(args.samples_root)
    source=json.loads((args.samples_root/'sampling_manifest.json').read_text())
    validate_training_donors(source['donors'])
    model_info=json.loads((args.model_root/'training_manifest.json').read_text())
    if digest(args.samples_root/'sampling_manifest.json')!=model_info['sampling_manifest_sha256']:
        raise ValueError('Samples differ from fitted model')
    out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=False)
    target=args.sensitivity_output.resolve();target.mkdir(parents=True,exist_ok=False)
    trees=model_info['config']['forest_trees'];rows=[]
    for record in source['donors']:
        cid=record['case_id'];file=args.samples_root/cid/'samples.npz'
        if digest(file)!=record['sample_sha256']:raise ValueError('Held-out sample changed')
        data=np.load(file);model=cv2.ml.RTrees_load(str(args.model_root/f'held_out_{cid}_forest.json'))
        clear=probability(model,data['x']);all_refs=probability(model,data['eval_x'])
        for count in range(trees//2,trees+1):
            cutoff=count/trees
            rows.append(dict(case_id=cid,minimum_votes=count,cutoff=cutoff,
                clear_foreground_recall=float((clear[data['y']==1]>=cutoff).mean()),
                clear_background_FPR=float((clear[data['y']==0]>=cutoff).mean()),
                all_background_FPR=float((all_refs[data['eval_y']==0]>=cutoff).mean()),
                all_interior_recall=float((all_refs[data['eval_y']==1]>=cutoff).mean())))
        print(f'Assessed held-out donor {cid}',flush=True)
    table=pd.DataFrame(rows);table.to_csv(out/'held_out_operating_points.csv',index=False)
    summary=table.groupby(['minimum_votes','cutoff']).agg(
        mean_background_FPR=('all_background_FPR','mean'),median_foreground_recall=('clear_foreground_recall','median'),
        min_foreground_recall=('clear_foreground_recall','min'),mean_clear_background_FPR=('clear_background_FPR','mean')).reset_index()
    summary.to_csv(out/'operating_point_summary.csv',index=False)
    eligible=summary[summary.mean_background_FPR<=.005]
    calibration=dict(status='source_reference_target_met' if len(eligible) else 'no_cutoff_met_background_reference_target',
        target_mean_background_reference_FPR=.005,selected_minimum_votes=int(eligible.iloc[0].minimum_votes) if len(eligible) else None,
        minimum_achieved_FPR=float(summary.mean_background_FPR.min()),
        minimum_achieved_at_votes=int(summary.loc[summary.mean_background_FPR.idxmin(),'minimum_votes']),
        HPAP_labels_used=False,validation_scope='Weak source-reference labels; no manual boundary accuracy estimate',
        operating_point_table_sha256=digest(out/'held_out_operating_points.csv'))
    write_json(out/'calibration_result.json',calibration)
    # Unanimity is a sensitivity regardless of the calibration outcome. It is not
    # retrospectively renamed the accepted primary operating point.
    started=time.monotonic()
    manifest=json.loads((args.prediction_root/'prediction_manifest.json').read_text())
    if manifest['model_sha256']!=model_info['forest_sha256']:raise ValueError('Prediction/model mismatch')
    p=np.load(args.prediction_root/'foreground_probability.npy',mmap_mode='r')
    # Recover integer votes exactly from float16 storage (half-unit error < .5 vote).
    votes=np.rint(p.astype(np.float32)*trees).astype(np.uint16)
    spec=Path(manifest['source_specification']);specification=json.loads(spec.read_text())
    if digest(spec.parent/'valid_tissue.geojson')!=specification['valid_tissue_sha256']:raise ValueError('Tissue changed')
    tissue=rasterize([shape(json.loads((spec.parent/'valid_tissue.geojson').read_text()))],p.shape)
    foreground=(votes==trees)&tissue;regions,foci=clean_regions(foreground,tissue)
    for name,array in [('pixel_classes',foreground),('boundary_regions',regions),('small_foci',foci)]:
        Image.fromarray(array.astype('uint8')).save(target/f'{name}.png')
    pixel_area=(manifest['source_pixel_size_um']*manifest['prediction_downsample'])**2
    manifest.update(status='unanimous_vote_sensitivity_not_selected_primary',
        foreground_area_um2=float(foreground.sum()*pixel_area),raw_connected_regions=int(ndi.label(foreground)[1]),
        boundary_connected_regions=int(ndi.label(regions)[1]),small_focus_regions=int(ndi.label(foci)[1]),
        boundary_area_um2=float(regions.sum()*pixel_area),
        operating_point=dict(minimum_votes=trees,forest_trees=trees,selection_status='Sensitivity only; not an accepted primary cutoff',
            donor_mean_background_reference_FPR=float(summary.iloc[-1].mean_background_FPR)),
        calibration_result_sha256=digest(out/'calibration_result.json'),
        source_prediction_manifest_sha256=digest(args.prediction_root/'prediction_manifest.json'),
        source_probability_sha256=digest(args.prediction_root/'foreground_probability.npy'),accepted_islets=0)
    manifest['shared_inference_runtime_seconds']=manifest.pop('runtime_seconds')
    manifest['postprocessing_runtime_seconds']=time.monotonic()-started
    manifest['config']['probability_cutoff']=1.
    manifest['raster_sha256']={f.name:digest(f) for f in target.glob('*.png')}
    write_json(target/'prediction_manifest.json',manifest)
    print(json.dumps(calibration),flush=True)


if __name__=='__main__':main()
