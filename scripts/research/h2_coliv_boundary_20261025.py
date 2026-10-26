#!/usr/bin/env python3
"""Existing-image ColIV boundary coverage. No source masks or app data changed.

Thresholded positive area follows Richardson2023; Otsu follows the supplied
Nature Communications 2026 IBD spatial methods review. The current endocrine
outlines and CD31 capillary objects are read as supplied. A symmetric +/-2 um
boundary band is primary; +/-5 um is a prespecified sensitivity analysis.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/h2-20261025-mpl')
import sys
import argparse
import json
import hashlib
from pathlib import Path
import cv2
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
from shapely.affinity import scale, translate
from skimage.filters import threshold_otsu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h1_features_20261022 as hf
import h1_analysis_20261022 as h1
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform' / 'typing'))
import codex_typing as ct
from redsea_cells import Channels

SRC = Path('data/development/uniform-20donor-20261015')
APP = Path('data/development/uniform-20donor-20261022-v10')
INPUT = Path('output/h1_edge_corrected_20261024/h1_object_features.parquet')
OUT = Path('output/h1_h4_literature_methods_20261025/h2')
PIXEL = hf.PIXEL_UM


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def paint(mask, geom):
    if geom.is_empty:
        return
    if geom.geom_type != 'Polygon':
        for g in geom.geoms:
            paint(mask, g)
        return
    ext = np.round(np.asarray(geom.exterior.coords) * 16).astype(np.int32)
    cv2.fillPoly(mask, [ext], 1, shift=4)
    for ring in geom.interiors:
        r = np.round(np.asarray(ring.coords) * 16).astype(np.int32)
        cv2.fillPoly(mask, [r], 0, shift=4)


def pixel_mask(geom, shape_yx, x0=0, y0=0, step=1):
    mask = np.zeros(shape_yx, np.uint8)
    g = translate(scale(geom, xfact=1 / PIXEL, yfact=1 / PIXEL, origin=(0, 0)),
                  xoff=-x0, yoff=-y0)
    if step != 1:
        g = scale(g, xfact=1 / step, yfact=1 / step, origin=(0, 0))
    paint(mask, g)
    return mask.astype(bool)


def boundary_features(ch, threshold, core, vessels, vtree, tissue=None):
    H, W = ch.shape
    x1, y1, x2, y2 = core.buffer(7).bounds
    x0 = max(0, int(np.floor(x1 / PIXEL))); y0 = max(0, int(np.floor(y1 / PIXEL)))
    xend = min(W, int(np.ceil(x2 / PIXEL)) + 1); yend = min(H, int(np.ceil(y2 / PIXEL)) + 1)
    patch = np.asarray(ch[y0:yend, x0:xend])
    vm = hf.union_local(vtree, vessels, core.buffer(7))
    vesselmask = pixel_mask(vm, patch.shape, x0, y0)
    coremask = pixel_mask(core, patch.shape, x0, y0)
    positive = patch > threshold
    tissue_mask = pixel_mask(tissue.intersection(core.buffer(7)), patch.shape, x0, y0) if tissue is not None else np.ones(patch.shape, bool)
    rec = dict(capillary_core_area_fraction=float((vesselmask & coremask).sum() / max(coremask.sum(), 1)))
    for width in (2, 5):
        bandgeom = core.buffer(width).difference(core.buffer(-width))
        band = pixel_mask(bandgeom, patch.shape, x0, y0) & tissue_mask
        nonvascular = band & ~vesselmask
        n = int(nonvascular.sum())
        rec.update({f'coliv_coverage_{width}um': float(positive[nonvascular].mean()) if n else np.nan,
                    f'coliv_coverage_all_{width}um': float(positive[band].mean()) if band.any() else np.nan,
                    f'boundary_pixels_{width}um': int(band.sum()),
                    f'nonvascular_boundary_pixels_{width}um': n,
                    f'vascular_boundary_fraction_{width}um': float(vesselmask[band].mean()) if band.any() else np.nan,
                    f'boundary_mean_coliv_{width}um': float(patch[nonvascular].mean()) if n else np.nan})
    return rec, (patch, positive, coremask, vesselmask, x0, y0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', nargs='*', type=int)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'by_donor').mkdir(exist_ok=True)
    manifest = ct.manifest(SRC)
    data = h1.eligible(pd.read_parquet(INPUT))
    data['case_id'] = data.case_id.astype(int)
    thresholds = []
    provenance = dict(feature_input=str(INPUT), feature_input_sha256=sha(INPUT), images=[])
    previews = []
    for tag, entry in sorted(manifest.items()):
        case = int(entry['case_id'])
        if args.cases is not None and case not in args.cases:
            continue
        d = data.loc[data.case_id.eq(case)]
        dest = OUT / 'by_donor' / f'{case}.csv'
        meta_dest = OUT / 'by_donor' / f'{case}.json'
        if dest.exists() and meta_dest.exists():
            meta = json.loads(meta_dest.read_text()); thresholds.append(meta['threshold_record']); provenance['images'].append(meta['provenance'])
            print(tag, 'already complete', flush=True)
            continue
        cores, _, tissue, _, gp = hf.geometry(APP, case)
        vp = SRC / 'structures' / 'Capillary' / f'Capillary__{tag}.geojson'
        vg = json.loads(vp.read_text())['features']
        vessels = np.array([shapely.make_valid(scale(shape(f['geometry']), xfact=PIXEL, yfact=PIXEL, origin=(0, 0))) for f in vg], dtype=object)
        vtree = shapely.STRtree(vessels)
        names = ct.channels(SRC, entry['entry_id'])
        img = Channels(entry['image_path'])
        # Memory-map uncompressed stacks; only one channel is decoded otherwise.
        ch = img.arr[names.index('ColIV')] if img.memmap else img(names.index('ColIV'))
        sampled = np.asarray(ch[::8, ::8])
        tm = pixel_mask(tissue, sampled.shape, step=8)
        vals = sampled[tm]
        threshold = float(threshold_otsu(vals))
        threshold_record = dict(case_id=case, tag=tag, n_sampled_tissue_pixels=len(vals), threshold_raw=threshold,
                                tissue_positive_fraction=float((vals > threshold).mean()),
                                raw_p50=float(np.median(vals)), raw_p99=float(np.percentile(vals, 99)))
        thresholds.append(threshold_record)
        recs = []
        preview_id = d.sort_values('core_area_um2').iloc[len(d) // 2].islet_id
        for row in d.itertuples():
            core = cores[row.base_islet_id][1]
            rec, arrays = boundary_features(ch, threshold, core, vessels, vtree)
            recs.append(dict(islet_id=row.islet_id, case_id=case, donor_status=row.donor_status, **rec))
            if row.islet_id == preview_id and tag in ('ND_6516', 'Aab_6450', 'T1D_6533'):
                patch, positive, cm, vm, x0, y0 = arrays
                np.savez_compressed(OUT / f'preview_{case}.npz', image=patch, positive=positive, core=cm, vessel=vm,
                                    threshold=threshold, islet_id=row.islet_id, pixel_um=PIXEL)
        frame = pd.DataFrame(recs)
        frame.to_csv(dest, index=False)
        stat = Path(entry['image_path']).stat()
        prov = dict(case_id=case, source_image=entry['image_path'], bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                    channel='ColIV', channel_index=names.index('ColIV'), source_shape=list(img.shape),
                    current_geometry=str(gp), geometry_sha256=sha(gp), capillary_objects=str(vp), capillary_sha256=sha(vp))
        provenance['images'].append(prov)
        meta_dest.write_text(json.dumps(dict(threshold_record=threshold_record, provenance=prov), indent=2))
        print(tag, len(frame), 'objects', 'threshold', round(threshold, 2), 'mean coverage', round(frame.coliv_coverage_2um.mean(), 3), flush=True)
        del ch
        img.close()
    if args.cases is not None:
        return
    # Tissue-edge sections are retained; measure only observed tissue in their
    # annulus. Cached interior sections are unaffected by this intersection.
    for tag, entry in sorted(manifest.items()):
        case = int(entry['case_id'])
        dp = OUT / 'by_donor' / f'{case}.csv'
        mp = OUT / 'by_donor' / f'{case}.json'
        metadata = json.loads(mp.read_text())
        if metadata.get('boundary_intersected_tissue'):
            continue
        edge = data.loc[data.case_id.eq(case) & data.core_distance_to_tissue_edge_um.le(6)]
        if len(edge):
            cores, _, tissue, _, _ = hf.geometry(APP, case)
            vp = SRC / 'structures' / 'Capillary' / f'Capillary__{tag}.geojson'
            vessels = np.array([shapely.make_valid(scale(shape(f['geometry']), xfact=PIXEL, yfact=PIXEL, origin=(0, 0)))
                               for f in json.loads(vp.read_text())['features']], dtype=object)
            vtree = shapely.STRtree(vessels)
            names = ct.channels(SRC, entry['entry_id']); img = Channels(entry['image_path'])
            ch = img.arr[names.index('ColIV')] if img.memmap else img(names.index('ColIV'))
            frame = pd.read_csv(dp).set_index('islet_id')
            for row in edge.itertuples():
                rec, _ = boundary_features(ch, metadata['threshold_record']['threshold_raw'], cores[row.base_islet_id][1], vessels, vtree, tissue)
                for k, v in rec.items():
                    frame.loc[row.islet_id, k] = v
            frame.reset_index().to_csv(dp, index=False)
            del ch; img.close()
        metadata['boundary_intersected_tissue'] = True
        metadata['tissue_edge_objects_remeasured'] = len(edge)
        mp.write_text(json.dumps(metadata, indent=2))
        print(tag, 'tissue intersection', len(edge), 'objects', flush=True)
    features = pd.concat([pd.read_csv(p) for p in sorted((OUT / 'by_donor').glob('*.csv'))], ignore_index=True)
    assert len(features) == len(data) and features.islet_id.is_unique
    merged = data.merge(features.drop(columns=['case_id', 'donor_status']), on='islet_id', validate='one_to_one')
    merged.to_csv(OUT / 'h2_analysis_objects.csv', index=False)
    pd.DataFrame(thresholds).to_csv(OUT / 'thresholds.csv', index=False)
    provenance['method'] = 'Per-image Otsu on native ColIV tissue pixels sampled every 8 pixels; fraction of pixels positive within symmetric boundary bands intersected with Tissue, excluding current Capillary masks. No morphological closing or historical ColIV masks.'
    (OUT / 'source_provenance.json').write_text(json.dumps(provenance, indent=2))
    plot_previews()


def plot_previews():
    paths = sorted(OUT.glob('preview_*.npz'))
    fig, axes = plt.subplots(len(paths), 2, figsize=(9, 4 * len(paths)), squeeze=False)
    for (left, right), path in zip(axes, paths):
        a = np.load(path); im = a['image']; vmax = np.percentile(im, 99.5)
        for ax in (left, right):
            ax.imshow(im, cmap='gray', vmin=0, vmax=vmax)
            ax.contour(a['core'], levels=[.5], colors=['#ffcc00'], linewidths=.6)
            bar=10/float(a['pixel_um']); by=im.shape[0]-5
            ax.plot([5,5+bar],[by,by],color='white',lw=3)
            ax.text(5,by-3,'10 µm',color='white',fontsize=8)
            ax.axis('off')
        rgba = np.zeros((*im.shape, 4)); rgba[a['positive']] = [.85, .1, .8, .5]; right.imshow(rgba)
        right.contour(a['vessel'], levels=[.5], colors=['#00ccff'], linewidths=.7)
        left.set_title(f"{a['islet_id']} · raw ColIV")
        right.set_title('Threshold-positive ColIV; CD31 vessels in cyan')
    fig.tight_layout(); fig.savefig(OUT / 'h2_image_examples.png', dpi=180); plt.close(fig)


if __name__ == '__main__':
    main()
