#!/usr/bin/env python3
"""Infer local, point-and-box prompted masks; save every failure for review."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from PIL import Image
from scipy import ndimage
from skimage import measure
from shapely.geometry import Polygon, mapping, shape
import torch


def supported_components(mask, points):
    labels, count = ndimage.label(mask)
    xy = np.rint(points).astype(int)
    xy[:, 0] = np.clip(xy[:, 0], 0, mask.shape[1]-1)
    xy[:, 1] = np.clip(xy[:, 1], 0, mask.shape[0]-1)
    support = labels[xy[:, 1], xy[:, 0]]
    keep = np.unique(support[support > 0])
    return np.isin(labels, keep), float((support > 0).mean()), count-len(keep)


def geometry(mask, scale, offset):
    # Contours are evaluated by nesting parity, preserving disconnected components
    # and holes. The original helper retained only the longest contour.
    contours = measure.find_contours(np.pad(mask.astype(float), 1), .5)
    result = Polygon()
    polygons = []
    for contour in contours:
        xy = (contour[:, ::-1]-1)*scale+offset
        polygon = Polygon(xy)
        if not polygon.is_valid: polygon = polygon.buffer(0)
        if not polygon.is_empty: polygons.append(polygon)
    for polygon in sorted(polygons, key=lambda p: p.area, reverse=True):
        result = result.symmetric_difference(polygon)
    return result


def main():
    # Geometry helpers remain importable for validation after the SAM pilot is
    # retired.  The optional model dependency is needed only for inference.
    from mobile_sam import SamPredictor, sam_model_registry

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--checkpoint', required=True, type=Path)
    p.add_argument('--case-ids', nargs='+', default=['78'])
    p.add_argument('--prompt-mode', choices=['points-and-box', 'box-only'], default='box-only')
    a = p.parse_args()
    torch.set_num_threads(4)
    torch.manual_seed(20260917)
    # This installed CUDA version does not promise deterministic cumsum.
    # Preserve that limitation explicitly, rather than claiming strict determinism.
    torch.use_deterministic_algorithms(True, warn_only=True)
    torch.backends.cudnn.benchmark = False
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = sam_model_registry['vit_b'](checkpoint=str(a.checkpoint)).to(device=device)
    model.eval()
    predictor = SamPredictor(model)
    checkpoint_hash = hashlib.sha256(a.checkpoint.read_bytes()).hexdigest()
    for cid in a.case_ids:
        folder = a.root/cid
        output = folder/'sam_predictions.geojson'
        if output.exists(): raise ValueError(f'Preserve previous predictions: {output}')
        spec = json.loads((folder/'prompts.json').read_text())
        tissue = shape(json.loads((folder/'valid_tissue.geojson').read_text()))
        features, audits = [], []
        started = time.time()
        prompts = [p for p in spec['prompts'] if 'crop' in p]
        for number, prompt in enumerate(prompts, 1):
            rgb = np.asarray(Image.open(folder/prompt['crop']).convert('RGB'))
            pts = np.asarray(prompt['points'])
            with torch.inference_mode():
                predictor.set_image(rgb)
                if a.prompt_mode == 'box-only':
                    masks, scores, _ = predictor.predict(box=np.asarray(prompt['box']), multimask_output=False)
                else:
                    masks, scores, _ = predictor.predict(point_coords=pts, point_labels=np.ones(len(pts)),
                        box=np.asarray(prompt['box']), multimask_output=True)
            choices = [supported_components(m, pts) for m in masks]
            # Prioritize seed retention, then the model score. Neither is identity validation.
            best = max(range(len(masks)), key=lambda i: (choices[i][1], float(scores[i]), -i))
            mask, recall, removed = choices[best]
            g = geometry(mask, prompt['downsample'], np.asarray(prompt['offset_px'])).intersection(tissue)
            edge = bool(mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any())
            signal = rgb.astype(float).max(axis=-1)/255
            surround = ndimage.binary_dilation(mask, iterations=max(1, round(10/(spec['pixel_size_um']*prompt['downsample'])))) & ~mask
            inside_mean = float(signal[mask].mean()) if mask.any() else 0.
            surround_mean = float(signal[surround].mean()) if surround.any() else 0.
            record = dict(candidate_id=prompt['candidate_id'], seed_count=prompt['seed_count'],
                          prompted_seed_recall=recall, predicted_mask_quality=float(scores[best]),
                          discarded_unsupported_components=int(removed), crop_edge_contact=edge,
                          mask_index=best, empty_geometry=g.is_empty,
                          mask_hormone_rgb_mean=inside_mean, surround_hormone_rgb_mean=surround_mean,
                          inverse_or_flat_signal_contrast=inside_mean <= surround_mean,
                          valid_geometry=g.is_valid, identity_reviewed=False)
            audits.append(record)
            if not g.is_empty:
                features.append(dict(type='Feature', geometry=mapping(g), properties=record))
            if number % 25 == 0 or number == len(prompts):
                print(f'{cid}: {number}/{len(prompts)} masks; {time.time()-started:.1f}s', flush=True)
        output.write_text(json.dumps(dict(type='FeatureCollection', features=features))+'\n')
        provenance = dict(case_id=cid, model='vit_b', checkpoint=str(a.checkpoint),
                          checkpoint_sha256=checkpoint_hash, torch_version=torch.__version__,
                          device=device, seed=20260917, runtime_seconds=time.time()-started,
                          prompt_mode=a.prompt_mode,
                          strict_determinism=False, deterministic_algorithms='warn_only; CUDA cumsum unsupported',
                          prompts_sha256=hashlib.sha256((folder/'prompts.json').read_bytes()).hexdigest(),
                          output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                          prompt_count=len(prompts), output_count=len(features), records=audits,
                          status='unreviewed boundary candidates')
        (folder/'prediction_audit.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__': main()
