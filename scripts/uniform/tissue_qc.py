#!/usr/bin/env python3
"""Check each image's Tissue outline, cells and endocrine objects against the image (owner standing rule, 2026-10-26).

For every image of a uniform run, measured against the image itself and against the original 15 images:
  - Tissue area; the proportion of the image's nucleated tissue (redraw_tissue.py's DAPI rule, no hole filling)
    inside the outline; interior holes under 500 um2;
  - cells per mm2 of Tissue and per mm2 of nucleated tissue, and the proportion of cells whose Parent is
    "Annotation (Tissue)", the value every later step keys tissue cells on (cells/<image>.parquet);
  - Islet and EndocrineCluster objects (islets/<image>.geojson) within 1 um of the Tissue outline.
A value outside the original 15's range is flagged. Renders, for the images named in --render:
  - <image>_tissue.png: DAPI at 1.51 um/px, Tissue outline yellow, area outside darkened;
  - <image>_cells.png: four 300 um patches at 0.755 um/px with every detected cell's centroid marked, chosen from
    where the image's nucleated tissue is densest in cells per area the least (the likeliest misses) and at random.

    python scripts/uniform/tissue_qc.py --out $OUT --report-dir <dir> [--render Aab_6090 ...]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import cv2
import numpy as np
import pyarrow.parquet as pq
import tifffile
from PIL import Image, ImageDraw, ImageFont
from shapely.geometry import shape
from shapely.prepared import prep
from skimage import filters

import redraw_tissue as RT

HERE = Path(__file__).resolve().parent
PATCH_UM = 300


def tissue_geometry(out: Path, tag: str):
    raw = json.loads((out / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    geoms = [shape(f['geometry']) for f in feats
             if ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue']
    if len(geoms) != 1:
        raise SystemExit(f'{tag}: expected one Tissue feature, found {len(geoms)}')
    return geoms[0]


def rasterize(geom, shape_, down: int) -> np.ndarray:
    m = np.zeros(shape_, np.uint8)
    for p in getattr(geom, 'geoms', [geom]):
        cv2.fillPoly(m, [np.round(np.asarray(p.exterior.coords) / down).astype(np.int32)], 1)
        for h in p.interiors:
            cv2.fillPoly(m, [np.round(np.asarray(h.coords) / down).astype(np.int32)], 0)
    return m.astype(bool)


def nucleated(dapi: np.ndarray) -> np.ndarray:
    d = np.log1p(filters.gaussian(dapi.astype(float), 2, preserve_range=True))
    return RT.above(d, d[d > 0], 0.6)


def measure(out: Path, tag: str) -> dict:
    path = RT.IMAGES / f'{tag}.ome.tiff'
    dapi, _, size = RT.channels(path)
    geom = tissue_geometry(out, tag)
    tis = rasterize(geom, dapi.shape, RT.DOWNSAMPLE)
    nuc = nucleated(dapi)
    um2 = RT.PIXEL_UM ** 2
    holes = [len([h for h in p.interiors if abs(cv2.contourArea(np.asarray(h.coords, np.float32))) * um2 < 500])
             for p in getattr(geom, 'geoms', [geom])]
    rec = {'image': tag, 'tissue_mm2': geom.area * um2 / 1e6, 'nucleated_mm2': nuc.sum() * RT.PX2 / 1e6,
           'nucleated_inside_pct': 100 * (nuc & tis).sum() / max(nuc.sum(), 1), 'holes_under_500um2': sum(holes)}
    cells_path = out / 'cells' / f'{tag}.parquet'
    if cells_path.is_file():
        n = pq.read_metadata(cells_path).num_rows
        parent = pq.read_table(cells_path, columns=['Parent']).column('Parent').to_pandas()
        rec.update(cells=n, cells_per_mm2_tissue=n / rec['tissue_mm2'], cells_per_mm2_nucleated=n / rec['nucleated_mm2'],
                   # the pipeline keys tissue cells on this exact Parent (islet_cell_test.TISSUE_PARENT)
                   cells_tissue_parent_pct=100 * float((parent == 'Annotation (Tissue)').mean()))
    islets = out / 'islets' / f'{tag}.geojson'
    if islets.is_file():
        raw = json.loads(islets.read_text())
        objs = [shape(f['geometry']) for f in (raw['features'] if isinstance(raw, dict) else raw)
                if ((f.get('properties') or {}).get('classification') or {}).get('name') in ('Islet', 'EndocrineCluster')]
        edge = prep(geom.boundary.buffer(1.0 / RT.PIXEL_UM))
        touching = sum(1 for o in objs if edge.intersects(o))
        rec.update(objects=len(objs), objects_touching_outline=touching,
                   objects_touching_pct=100 * touching / max(len(objs), 1))
    return rec


def gray(d):
    v = np.log1p(d.astype(np.float32))
    lo, hi = np.percentile(v, [1, 99.7])
    return (np.clip((v - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)


def label(img: np.ndarray, text: str, size: int = 30) -> Image.Image:
    im = Image.fromarray(img)
    out = Image.new('RGB', (im.width, im.height + size + 16))
    ImageDraw.Draw(out).text((10, 6), text, fill=(255, 255, 255), font=ImageFont.load_default(size=size))
    out.paste(im, (0, size + 16))
    return out


def render(out: Path, tag: str, report: Path) -> None:
    path = RT.IMAGES / f'{tag}.ome.tiff'
    geom = tissue_geometry(out, tag)
    with tifffile.TiffFile(path) as t:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        d2 = t.series[0].levels[2].pages[names.index('DAPI')].asarray()
        d1 = t.series[0].levels[1].pages[names.index('DAPI')].asarray()
    g = gray(d2)
    tis = rasterize(geom, g.shape, 4)
    rgb = np.repeat(g[..., None], 3, 2).astype(np.float32)
    rgb[~tis] *= 0.3
    edge = cv2.morphologyEx(tis.astype(np.uint8), cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)).astype(bool)
    rgb[edge] = (255, 225, 0)
    label(rgb.astype(np.uint8), f'{tag}: Tissue outline (yellow), outside darkened; DAPI 1.51 um/px').save(
        report / f'{tag}_tissue.png', optimize=True)

    cells = pq.read_table(out / 'cells' / f'{tag}.parquet', columns=['Centroid X µm', 'Centroid Y µm']).to_pandas()
    cx = (cells['Centroid X µm'].to_numpy() / RT.PIXEL_UM / 2).astype(int)
    cy = (cells['Centroid Y µm'].to_numpy() / RT.PIXEL_UM / 2).astype(int)
    g1 = gray(d1)
    p = int(PATCH_UM / (RT.PIXEL_UM * 2))
    tis1 = rasterize(geom, g1.shape, 2)
    H, W = g1.shape[0] // p, g1.shape[1] // p
    density = np.zeros((H, W))
    np.add.at(density, (np.clip(cy // p, 0, H - 1), np.clip(cx // p, 0, W - 1)), 1)
    cover = tis1[:H * p, :W * p].reshape(H, p, W, p).mean((1, 3))
    nuc = gray(d1)[:H * p, :W * p].reshape(H, p, W, p).mean((1, 3))
    ok = np.argwhere(cover > 0.9)
    per = density[ok[:, 0], ok[:, 1]] / np.maximum(nuc[ok[:, 0], ok[:, 1]], 1)
    picks = [tuple(ok[i]) for i in np.argsort(per)[:2]]
    rng = np.random.default_rng(0)
    picks += [tuple(ok[i]) for i in rng.choice(len(ok), 2, replace=False)]
    tiles = []
    for r, c in picks:
        y0, x0 = r * p, c * p
        rgb = np.repeat(g1[y0:y0 + p, x0:x0 + p, None], 3, 2).copy()
        sel = (cx >= x0) & (cx < x0 + p) & (cy >= y0) & (cy < y0 + p)
        for x, y in zip(cx[sel] - x0, cy[sel] - y0):
            cv2.circle(rgb, (int(x), int(y)), 3, (255, 40, 200), -1)
        tiles.append(np.asarray(label(rgb, f'x={x0 * 2}-{(x0 + p) * 2} y={y0 * 2}-{(y0 + p) * 2} px: {sel.sum()} cells', 22)))
    pair = lambda a, b: np.hstack([a, np.zeros((a.shape[0], 10, 3), np.uint8), b])
    grid = np.vstack([pair(tiles[0], tiles[1]), np.zeros((10, tiles[0].shape[1] * 2 + 10, 3), np.uint8), pair(tiles[2], tiles[3])])
    label(grid, f'{tag}: detected cell centroids (magenta) on DAPI, 0.755 um/px; top: fewest cells per DAPI signal; bottom: random').save(
        report / f'{tag}_cells.png', optimize=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--report-dir', type=Path, required=True)
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    ap.add_argument('--render', nargs='*', default=[])
    a = ap.parse_args()
    a.report_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(a.cohort.open()))
    recs = []
    for r in rows:
        rec = measure(a.out, r['image_name'])
        rec['cohort'] = r['cohort']
        recs.append(rec)
        print(json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in rec.items()}), flush=True)
    keys = ['nucleated_inside_pct', 'holes_under_500um2', 'cells_per_mm2_tissue', 'cells_per_mm2_nucleated',
            'cells_tissue_parent_pct', 'objects_touching_pct']
    ref = [x for x in recs if x['cohort'] == 'original']
    for x in recs:
        x['flags'] = ';'.join(k for k in keys if k in x and all(k in y for y in ref) and
                              not min(y[k] for y in ref) <= x[k] <= max(y[k] for y in ref))
    fields = list(dict.fromkeys(k for x in recs for k in x))
    with open(a.report_dir / 'tissue_qc.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(recs)
    for tag in a.render:
        render(a.out, tag, a.report_dir)
        print(f'{tag}: rendered', flush=True)
    flagged = [(x['image'], x['flags']) for x in recs if x['flags']]
    print('outside the original 15 range: ' + (', '.join(f'{t} ({f})' for t, f in flagged) if flagged else 'none'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
