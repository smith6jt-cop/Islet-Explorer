#!/usr/bin/env python3
"""Check all 20 content-rule Tissue outlines against the image before downstream use.

Per image: the outline in use (old) and the content-rule redraw (new), at pyramid level 3.
Quantities: areas; area only in the new (added) and only in the old (lost); the nuclei the rule sees (DAPI) inside
each; pieces and perimeter; the run's detected cells and endocrine objects (islets, clusters, singles) that lie in the
lost area (they would fall outside Tissue). The five additions (their outline in use is already this rule) are listed
for the cohort comparison.

Writes output/tissue_redraw_check_20261029/check.csv and per image <image>_whole.png (DAPI blue, stroma green; old
outline white, new outline magenta; added area tinted magenta, lost area red) and <image>_tiles.png (six 500 um tiles:
the four largest added regions and the two largest lost regions, read at level 1). Individual tiles and their physical
coordinates are also saved. Cell, nucleus and endocrine containment is checked in exact vector geometry, never from
downsampled masks. Raster comparisons describe signal coverage and display only.

    python scripts/research/tissue_redraw_check_20261029.py --new-dir <dir of redrawn outlines>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import shapely
from scipy import ndimage

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform'))
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
sys.path.insert(0, str(REPO / 'scripts' / 'research'))
import codex_typing as CT  # noqa: E402
import redraw_tissue as RT  # noqa: E402
import tissue_outline_consistency_20261029 as TC  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'tissue_redraw_check_20261029'
ADDITIONS = ('Aab_6090', 'Aab_6558', 'ND_6232', 'ND_6539', 'T1D_6520')


def geometry_from(path: Path):
    raw = json.loads(path.read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    return shapely.union_all([shapely.make_valid(shapely.geometry.shape(f['geometry'])) for f in feats
                              if ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue'])


def mask_from(path: Path, shape, downsample: int = RT.DOWNSAMPLE) -> np.ndarray:
    raw = json.loads(path.read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    m = np.zeros(shape, np.uint8)
    for f in feats:
        if ((f.get('properties') or {}).get('classification') or {}).get('name') != 'Tissue':
            continue
        g = f['geometry']
        for poly in (g['coordinates'] if g['type'] == 'MultiPolygon' else [g['coordinates']]):
            rings = [np.round(np.asarray(r, float) / downsample * 16).astype(np.int32) for r in poly if len(r) >= 3]
            if rings:
                # Rasterize each polygon separately before union. A hole in one
                # polygon must not erase a disjoint island drawn earlier inside it.
                x0, y0 = np.maximum(rings[0].min(axis=0) // 16 - 1, 0)
                x1, y1 = np.minimum((rings[0].max(axis=0) + 15) // 16 + 2,
                                    [shape[1], shape[0]])
                if x1 <= x0 or y1 <= y0:
                    continue
                local = np.zeros((y1 - y0, x1 - x0), np.uint8)
                shifted = [r - [x0 * 16, y0 * 16] for r in rings]
                cv2.fillPoly(local, [shifted[0]], 1, lineType=cv2.LINE_8, shift=4)
                if len(rings) > 1:
                    cv2.fillPoly(local, shifted[1:], 0, lineType=cv2.LINE_8, shift=4)
                m[y0:y1, x0:x1] |= local
    return m.astype(bool)


def perimeter_mm(m: np.ndarray) -> float:
    cs, _ = cv2.findContours(m.astype(np.uint8), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    return sum(cv2.arcLength(c, True) for c in cs) * RT.PIXEL_UM * RT.DOWNSAMPLE / 1000


def one(tag: str, new_dir: Path, run: Path = RUN, out_dir: Path = OUT, compare_dir: Path | None = None) -> dict:
    dapi, stroma, _ = RT.channels(RT.IMAGES / f'{tag}.ome.tiff')
    nuclei, strom, _, _ = TC.rule_parts(dapi, stroma)
    old_path = run / 'inputs' / 'tissue' / f'{tag}.geojson'
    old = mask_from(old_path, dapi.shape)
    newp = new_dir / f'{tag}.geojson'
    if not newp.exists():
        raise FileNotFoundError(newp)
    new = mask_from(newp, dapi.shape)
    added, lost = new & ~old, old & ~new
    px = RT.PX2 / 1e6
    cell_frame = pd.read_parquet(run / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    xy = cell_frame[['Centroid X µm', 'Centroid Y µm']].to_numpy()
    old_geom, new_geom = geometry_from(old_path), geometry_from(newp)
    shapely.prepare(old_geom)
    shapely.prepare(new_geom)
    points = shapely.points(xy / RT.PIXEL_UM)
    in_old, in_new = shapely.covers(old_geom, points), shapely.covers(new_geom, points)
    missing = {'cell': 0, 'nucleus': 0, 'endocrine': 0}
    maximum_missing_area = 0.0
    minimum_endocrine_fraction = 1.0
    footprint_diagnostics = {}
    for batch, kinds in RT.footprint_batches(run, tag, diagnostics=footprint_diagnostics):
        covered = shapely.covers(new_geom, batch)
        for i in np.flatnonzero(~covered):
            difference = batch[i].difference(new_geom).area
            maximum_missing_area = max(maximum_missing_area, difference * RT.PIXEL_UM ** 2)
            # GEOS may report uncovered coincident boundaries with zero lost area.
            # Retain strict covers counts alongside the exact missing-area measure.
            missing[kinds[i]] += 1
            if kinds[i] == 'endocrine':
                minimum_endocrine_fraction = min(minimum_endocrine_fraction, 1 - difference / batch[i].area)
    pieces_old, pieces_new = ndimage.label(old)[1], ndimage.label(new)[1]
    rec = {'image': tag, 'cohort': 'addition' if tag in ADDITIONS else 'original',
           'old_mm2': old.sum() * px, 'new_mm2': new.sum() * px, 'added_mm2': added.sum() * px, 'lost_mm2': lost.sum() * px,
           'added_pct_of_old': 100 * added.sum() / old.sum(), 'lost_pct_of_old': 100 * lost.sum() / old.sum(),
           'nuclei_in_old_pct': 100 * (nuclei & old).sum() / nuclei.sum(),
           'nuclei_in_new_pct': 100 * (nuclei & new).sum() / nuclei.sum(),
           'added_with_nuclei_pct_of_added': 100 * (added & nuclei).sum() / max(added.sum(), 1),
           'cells_in_lost_area': int((in_old & ~in_new).sum()), 'cells_detected': len(xy),
           'cells_outside_old': int((~in_old).sum()), 'cells_outside_new': int((~in_new).sum()),
           'cell_footprints_not_covered': missing['cell'], 'nucleus_footprints_not_covered': missing['nucleus'],
           'endocrine_footprints_not_covered': missing['endocrine'],
           'maximum_missing_footprint_area_um2': maximum_missing_area,
           'minimum_endocrine_fraction_inside': minimum_endocrine_fraction,
           'empty_source_footprints': len(footprint_diagnostics.get('empty_source_footprints', [])),
           'pieces_old': pieces_old, 'pieces_new': pieces_new,
           'perimeter_old_mm': perimeter_mm(old), 'perimeter_new_mm': perimeter_mm(new),
           'cells_per_mm2_old': len(xy) / (old.sum() * px),
           'content_in_old_pct': 100 * (old & (nuclei | strom)).sum() / old.sum(),
           'content_in_new_pct': 100 * (new & (nuclei | strom)).sum() / new.sum(),
           'lost_with_content_pct_of_lost': 100 * (lost & (nuclei | strom)).sum() / max(lost.sum(), 1),
           'exact_new_area_mm2': new_geom.area * RT.PIXEL_UM ** 2 / 1e6,
           'exact_new_piece_count': len(RT.polygon_parts(new_geom)),
           'exact_new_hole_count': sum(len(p.interiors) for p in RT.polygon_parts(new_geom)),
           'exact_new_perimeter_mm': new_geom.length * RT.PIXEL_UM / 1000,
           'cells_per_mm2_new': len(xy) / (new_geom.area * RT.PIXEL_UM ** 2 / 1e6),
           'whole_pyramid_level': RT.LEVEL, 'tile_pyramid_level': 1,
           'new_outline_sha256': hashlib.sha256(newp.read_bytes()).hexdigest(),
           'old_outline_sha256': hashlib.sha256(old_path.read_bytes()).hexdigest(),
           'new_outline_valid': bool(new_geom.is_valid)}
    provenance_path = new_dir / f'{tag}.preservation.json'
    if provenance_path.exists():
        provenance = json.loads(provenance_path.read_text())
        rec['content_rule_version'] = provenance.get('content_rule', {}).get('rule_version', 'unrecorded')
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f'{tag}_source_geometry_notes.json').write_text(json.dumps(footprint_diagnostics, indent=2) + '\n')
    exceptions = cell_frame.loc[~in_new].copy()
    exceptions['outside_old_tissue'] = ~in_old[~in_new]
    exceptions.to_csv(out_dir / f'{tag}_centroid_exceptions.csv', index=False)

    def u8(a):
        lo, hi = np.percentile(a, [1, 99.5])
        return (np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)
    rgb = np.dstack([np.zeros(dapi.shape, np.uint8), u8(stroma), u8(np.log1p(dapi.astype(float)))])
    over = rgb.copy()
    over[added] = (0.5 * over[added] + 0.5 * np.array([255, 0, 255])).astype(np.uint8)
    over[lost] = (0.5 * over[lost] + 0.5 * np.array([255, 0, 0])).astype(np.uint8)
    k = np.ones((3, 3), np.uint8)
    over[cv2.morphologyEx(old.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 255, 255)
    over[cv2.morphologyEx(new.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 0, 255)
    cv2.imwrite(str(out_dir / f'{tag}_whole.png'), np.hstack([rgb, over])[:, :, ::-1])
    # Load genuine level-1 data: no enlarged low-resolution image is presented as a zoom.
    hi_dapi, hi_stroma, _ = RT.channels(RT.IMAGES / f'{tag}.ome.tiff', level_index=1)
    hi_rgb = np.dstack([np.zeros(hi_dapi.shape, np.uint8), u8(hi_stroma), u8(np.log1p(hi_dapi.astype(float)))])
    hi_old = mask_from(old_path, hi_dapi.shape, downsample=2)
    hi_new = mask_from(newp, hi_dapi.shape, downsample=2)
    del hi_dapi, hi_stroma
    tiles = []
    tile_rows = []
    for region, n in ((added, 4), (lost, 2)):
        lab, m = ndimage.label(region)
        if not m:
            continue
        sizes = np.bincount(lab.ravel())[1:]
        for i in np.argsort(sizes)[::-1][:n]:
            ys, xs = np.nonzero(lab == i + 1)
            cy, cx = int(ys.mean()), int(xs.mean())
            width = round(500 / (RT.PIXEL_UM * 2))
            cy, cx = cy * 4, cx * 4
            y0, x0 = max(0, min(cy - width // 2, hi_rgb.shape[0] - width)), max(0, min(cx - width // 2, hi_rgb.shape[1] - width))
            raw_tile = hi_rgb[y0:y0 + width, x0:x0 + width]
            o, nmask = hi_old[y0:y0 + width, x0:x0 + width], hi_new[y0:y0 + width, x0:x0 + width]
            overlay = raw_tile.copy()
            overlay[nmask & ~o] = (0.5 * overlay[nmask & ~o] + [127.5, 0, 127.5]).astype(np.uint8)
            overlay[o & ~nmask] = (0.5 * overlay[o & ~nmask] + [127.5, 0, 0]).astype(np.uint8)
            overlay[cv2.morphologyEx(o.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 255, 255)
            overlay[cv2.morphologyEx(nmask.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 0, 255)
            t = np.hstack([raw_tile, overlay])
            label = f'{tag} | {"added" if region is added else "removed"} {len(tiles)+1} | 500 um | raw / outlines'
            t = cv2.copyMakeBorder(t, 28, 0, 0, 0, cv2.BORDER_CONSTANT)
            cv2.putText(t, label, (8, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1, cv2.LINE_AA)
            tiles.append(t)
            cv2.imwrite(str(out_dir / f'{tag}_tile{len(tiles):02d}.png'), t[:, :, ::-1])
            tile_rows.append({'tile': len(tiles), 'region': 'added' if region is added else 'removed',
                              'x_um': x0 * RT.PIXEL_UM * 2, 'y_um': y0 * RT.PIXEL_UM * 2, 'width_um': width * RT.PIXEL_UM * 2,
                              'pyramid_level': 1})
    if tiles:
        cv2.imwrite(str(out_dir / f'{tag}_tiles.png'), np.vstack(tiles)[:, :, ::-1])
    (out_dir / f'{tag}_tiles.json').write_text(json.dumps(tile_rows, indent=2) + '\n')
    if compare_dir is not None:
        previous_path = compare_dir / f'{tag}.geojson'
        previous_geom = geometry_from(previous_path)
        delta = {'previous_outline_sha256': hashlib.sha256(previous_path.read_bytes()).hexdigest(),
                 'new_outline_sha256': rec['new_outline_sha256'],
                 'equal_geometry': previous_geom.equals(new_geom),
                 'added_um2': new_geom.difference(previous_geom).area * RT.PIXEL_UM ** 2,
                 'removed_um2': previous_geom.difference(new_geom).area * RT.PIXEL_UM ** 2,
                 'tiles': []}
        previous = mask_from(previous_path, dapi.shape)
        changed, count = ndimage.label(previous != new)
        sizes = np.bincount(changed.ravel())[1:]
        hi_previous = mask_from(previous_path, hi_rgb.shape[:2], downsample=2)
        for index in np.argsort(sizes)[::-1][:min(count, 3)]:
            yy, xx = np.nonzero(changed == index + 1)
            width = round(500 / (RT.PIXEL_UM * 2))
            y0 = max(0, min(int(yy.mean()) * 4 - width // 2, hi_rgb.shape[0] - width))
            x0 = max(0, min(int(xx.mean()) * 4 - width // 2, hi_rgb.shape[1] - width))
            raw_tile = hi_rgb[y0:y0 + width, x0:x0 + width]
            prior = hi_previous[y0:y0 + width, x0:x0 + width]
            revised = hi_new[y0:y0 + width, x0:x0 + width]
            overlay = raw_tile.copy()
            overlay[prior & ~revised] = (overlay[prior & ~revised] * .5 + [127.5, 0, 0]).astype(np.uint8)
            overlay[revised & ~prior] = (overlay[revised & ~prior] * .5 + [127.5, 0, 127.5]).astype(np.uint8)
            overlay[cv2.morphologyEx(prior.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 255, 255)
            overlay[cv2.morphologyEx(revised.astype(np.uint8), cv2.MORPH_GRADIENT, k).astype(bool)] = (255, 0, 255)
            filename = f'{tag}_conversion_change{len(delta["tiles"]) + 1:02d}.png'
            panel = cv2.copyMakeBorder(np.hstack([raw_tile, overlay]), 28, 0, 0, 0, cv2.BORDER_CONSTANT)
            cv2.putText(panel, f'{tag} | candidate comparison | prior white / new magenta / removed red',
                        (8, 20), cv2.FONT_HERSHEY_SIMPLEX, .5, (255, 255, 255), 1, cv2.LINE_AA)
            cv2.imwrite(str(out_dir / filename), panel[:, :, ::-1])
            delta['tiles'].append({'file': filename, 'x_um': x0 * RT.PIXEL_UM * 2,
                                   'y_um': y0 * RT.PIXEL_UM * 2, 'width_um': width * RT.PIXEL_UM * 2,
                                   'pyramid_level': 1})
        (out_dir / f'{tag}_conversion_comparison.json').write_text(json.dumps(delta, indent=2) + '\n')
    (out_dir / f'{tag}_check.json').write_text(json.dumps(rec, indent=2) + '\n')
    print(tag, {k2: round(v, 3) if isinstance(v, float) else v for k2, v in rec.items() if k2 not in ('image',)},
          flush=True)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--new-dir', type=Path, required=True)
    ap.add_argument('--run', type=Path, default=RUN)
    ap.add_argument('--out-dir', type=Path, default=OUT)
    ap.add_argument('--images', nargs='+')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--compare-dir', type=Path, help='Previous candidate for exact comparison and changed-region tiles')
    a = ap.parse_args()
    tags = a.images or [r['image_name'] for r in CT.cohort_rows()]
    with ProcessPoolExecutor(a.workers) as ex:
        rows = list(ex.map(one, tags, repeat(a.new_dir), repeat(a.run), repeat(a.out_dir), repeat(a.compare_dir)))
    # Keep completed records when reviewing the cohort in bounded batches.
    d = pd.DataFrame([json.loads(p.read_text()) for p in sorted(a.out_dir.glob('*_check.json'))])
    d.to_csv(a.out_dir / 'check.csv', index=False)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
