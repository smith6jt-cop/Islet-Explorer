#!/usr/bin/env python3
"""Collagen IV area objects, including acinar mesh and vessel collagen (owner 2026-10-29).

Current default `--method all-tissue`: local log1p nuclear-reference median plus one lower-side robust scale.
The reference is the union of this image's current nucleus polygons, eroded two native pixels and intersected
with current Tissue. It uses no phenotype labels. In each 200 um block the scale is the median lower-side
absolute deviation divided by 0.6744897501960817. This is an empirical mesh-retaining cutoff, not a claimed
Gaussian specificity or independent false-positive estimate. Nuclear interiors can contain spillover or
segmentation error; positive-reference residuals and unsupported blocks are reported for review. Missing or
degenerate block references borrow the nearest valid block, then a 3x3 log-grid median and raw bilinear
interpolation form the map. An image with no valid reference block fails. Donors are never pooled.
Objects are clipped to exact Tissue and edge-touching components are retained.
Per-image levels, geometry/code hashes, area diagnostics and whole-image/native-resolution zoom renders are saved.

Historical `--method above-mesh` below is retained only to reproduce the superseded comparison:

Owner 2026-10-29: "ColIV needs to be converted to objects like vessels, lymphatics, and nerves." There is no original
classifier for ColIV to fit a threshold against, so the image sets its own, by the typing's known-negative method
(negative_p99) applied locally: in each BLOCK_UM square block, the 99th percentile of the smoothed channel over the pixels
of Exocrine cells (typing/identity) at least REFERENCE_CLEARANCE_UM from every Capillary, Lymphatic and Nerve object
and every endocrine object (their basement membranes lie within about 1 um of the objects; the typing's 20 um
ISOLATION_UM left 0-11 % of blocks with enough reference pixels on five images checked, 93-100 % at 5 um); a block with fewer reference pixels than BLOCK_MIN_REFERENCE of its area takes the nearest
block's value; the block values are smoothed (3 x 3 median, log scale) and interpolated bilinearly to a threshold per
pixel. Every acinus is outlined by a dim collagen IV basement membrane; an object is collagen IV brighter than the mesh
around it. Local since the first image check (2026-10-29, all 40 renders): one threshold per image traced the acinar
mesh wherever the mesh was brighter than the image's 99th percentile (ND_6548, Aab_6538, ND_6539, the bright lobules of
ND_6516 and T1D_6533, a brighter acquisition tile of Aab_6090) and missed capillaries where it was dimmer (T1D_6550,
ND_6516's upper right); those objects are kept in structures/superseded_collagen_iv_global_p99/. Agent determination
2026-10-29, compared on ND_0112, ND_6232 and T1D_6533 against the background peak plus 1, 2 or 3 widths (3 widths drew
0.2-2 % of tissue and missed most vessels; 1 width outlined the acinar mesh; 2 widths reached 44-93 % of capillary
pixels), a Gaussian-mixture or Otsu split (42-66 % of tissue, one tissue-wide network) and the triangle method
(2-95 %). Drawing as structure_thresholds.make: the despeckled channel (redsea_cells.despeckle; owner 2026-10-28)
Gaussian-smoothed at SIGMA px (the 15's best Dice for Nerve and Capillary), >= threshold inside Tissue, 8-connected
components as separate objects, at least MIN_AREA_UM2 (Capillary's). No object is dropped for touching the image edge
(owner 2026-10-29: "the edge rule is only supposed to drop objects that are significantly obscurred by being cut off on
the edge e.g. an small islet cut off into a fraction"; the objects' only uses are density in tissue and endocrine objects
and distances of cells to them, which a cut does not distort). The first draw dropped every component touching the
10 um edge ring, as structure_thresholds.make does, which removed whole artery walls (16,726 um2 on Aab_6450).

Writes structures/Collagen IV/Collagen IV__<image>.geojson, structures/collagen_iv.json (threshold and checks per
image) and output/structure_check_20261029/collagen_iv/ renders (whole image and zoomed tiles: ColIV with DAPI and
the object outlines) for the image check every gating geometry needs (CLAUDE.md, 2026-10-26).

    python scripts/uniform/typing/collagen_structures.py --run data/development/uniform-20donor-20261029-v19-work --workers 2 --resume
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import shapely
from scipy.ndimage import distance_transform_edt, gaussian_filter, median_filter
from shapely.geometry import mapping

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import redsea_cells as RS  # noqa: E402
import structure_thresholds as ST  # noqa: E402
import levels as LV  # noqa: E402

CLASS = 'Collagen IV'
MARKER = 'ColIV'
SIGMA = 2.0
MIN_AREA_UM2 = 4.0
METHOD_VERSION = 'local_nuclear_log_median_plus_one_lower_sigma/1'
NUCLEAR_EROSION_PX = 2
NUCLEAR_MIN_PIXELS = 1000
NUCLEAR_MIN_FRACTION = 0.005
LOWER_MAD_NORMALIZER = 0.6744897501960817
REFERENCE_QUANTILE = 99.0           # negative_p99
REFERENCE_STRUCTURES = ('Capillary', 'Lymphatic', 'Nerve')
REFERENCE_CLEARANCE_UM = 5.0
BLOCK_UM = 200.0
BLOCK_MIN_REFERENCE = 0.05          # proportion of a block's pixels that must be reference pixels for its own value
PREVIOUS = 'superseded_collagen_iv_global_p99'
# spots the first image check flagged, in its _whole.png coordinates (x in the outlined right half unless < width)
CHECKPOINTS = {'Aab_6090': [(335, 1085)], 'Aab_6450': [(3415, 1440)], 'Aab_6505': [(2125, 285), (3215, 50)],
               'Aab_6538': [(3450, 925)], 'Aab_6558': [(2880, 940), (3820, 920)], 'ND_0112': [(3910, 190), (2650, 1210)],
               'ND_6356': [(2240, 175), (3560, 1147)], 'ND_6516': [(2450, 1175), (3650, 475)], 'ND_6539': [(2250, 775)],
               'ND_6548': [(3150, 975)], 'T1D_6520': [(3760, 900)], 'T1D_6533': [(3725, 1000), (2580, 1322)],
               'T1D_6534': [(3925, 1365)], 'T1D_6550': [(3850, 920)], 'T1D_6563': [(3745, 230)]}
REPO = Path(__file__).resolve().parents[3]
RENDERS = REPO / 'output' / 'structure_check_20261029' / 'collagen_iv'


def load(run: Path, tag: str, marker: str) -> np.ndarray:
    entry = CT.manifest(run)[tag]
    names = CT.channels(run, entry['entry_id'])
    img = RS.Channels(entry['image_path'])
    try:
        return RS.despeckle(img.native(names.index(marker))).astype(np.float32)
    finally:
        img.close()


def stretch(a: np.ndarray, region: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(a[region], [1, 99.5]) if region.any() else (0, 1)
    return (np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def outline(mask: np.ndarray) -> np.ndarray:
    return cv2.morphologyEx(mask.astype(np.uint8), cv2.MORPH_GRADIENT, np.ones((3, 3), np.uint8)).astype(bool)


def render(tag: str, col: np.ndarray, dapi: np.ndarray, tissue: np.ndarray, binary: np.ndarray, capillary: np.ndarray,
           previous: np.ndarray | None, renders: Path = RENDERS):
    """Whole image (downsampled: raw | outlined) and zoomed tiles (raw | outlined | previous objects): ColIV green, DAPI
    blue, Tissue outline cyan, Capillary objects' outline yellow, the Collagen IV objects' outline magenta on top.
    Tiles: four at seeded Tissue positions, then one per CHECKPOINTS spot of the image,
    then additional seeded Tissue positions to reach at least six without changing earlier tile order."""
    renders.mkdir(parents=True, exist_ok=True)
    g, b = stretch(col, tissue), stretch(dapi, tissue)
    rgb = np.dstack([np.zeros_like(g), g, b])
    t_edge, cap = outline(tissue), outline(capillary)

    def draw(objs):
        o = rgb.copy()
        o[t_edge] = (0, 255, 255)
        o[cap] = (255, 255, 0)
        o[outline(objs)] = (255, 0, 255)
        return o
    over = draw(binary)
    before = draw(previous) if previous is not None else np.zeros_like(rgb)
    H, W = col.shape
    f = 2000 / max(H, W)
    small = cv2.resize(np.hstack([rgb, over]), (int(2 * W * f), int(H * f)), interpolation=cv2.INTER_AREA)
    cv2.imwrite(str(renders / f'{tag}_whole.png'), small[:, :, ::-1])
    ys, xs = np.nonzero(tissue[::50, ::50])
    rng = np.random.default_rng(0)
    seeded_indices = rng.choice(len(ys), size=min(4, len(ys)), replace=False)
    centers = [(ys[i] * 50, xs[i] * 50) for i in seeded_indices]
    for xw, yw in CHECKPOINTS.get(tag, []):
        x = (xw - W * f) / f if xw >= W * f else xw / f
        centers.append((int(yw / f), int(x)))
    if len(centers) < 6:
        remaining = np.setdiff1d(np.arange(len(ys)), seeded_indices)
        extra = rng.choice(remaining, size=min(6 - len(centers), len(remaining)), replace=False)
        centers.extend((ys[i] * 50, xs[i] * 50) for i in extra)
    tiles = []
    tile_positions = []
    for cy, cx in centers:
        y, x = int(min(max(cy - 400, 0), H - 800)), int(min(max(cx - 400, 0), W - 800))
        tile = np.hstack([rgb[y:y + 800, x:x + 800], over[y:y + 800, x:x + 800], before[y:y + 800, x:x + 800]])
        tiles.append(tile)
        cv2.imwrite(str(renders / f'{tag}_tile{len(tiles):02}.png'), tile[:, :, ::-1])
        tile_positions.append({'tile':len(tiles), 'x_px':x, 'y_px':y, 'pixel_um':CT.PIXEL_UM,
                               'panels':['raw','current','previous']})
    if tiles:
        sheet = np.vstack(tiles)
        cv2.imwrite(str(renders / f'{tag}_tiles.png'), cv2.resize(sheet, None, fx=0.75, fy=0.75,
                                                                  interpolation=cv2.INTER_AREA)[:, :, ::-1])
    (renders / f'{tag}_tiles.json').write_text(json.dumps(tile_positions, indent=2)+'\n')


def threshold_map(sm: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, dict]:
    """A threshold per pixel: each block's REFERENCE_QUANTILE over its reference pixels, blocks without enough of them
    from the nearest block that has, smoothed and interpolated (module docstring)."""
    H, W = sm.shape
    b = int(round(BLOCK_UM / CT.PIXEL_UM))
    ny, nx = -(-H // b), -(-W // b)
    grid = np.full((ny, nx), np.nan)
    for i in range(ny):
        for j in range(nx):
            r = ref[i * b:(i + 1) * b, j * b:(j + 1) * b]
            if r.sum() >= BLOCK_MIN_REFERENCE * r.size:
                grid[i, j] = np.percentile(sm[i * b:(i + 1) * b, j * b:(j + 1) * b][r], REFERENCE_QUANTILE)
    own = np.isfinite(grid)
    if not own.any():
        raise ValueError('no block with enough reference pixels')
    near = distance_transform_edt(~own, return_distances=False, return_indices=True)
    filled = grid[near[0], near[1]]
    smooth = np.expm1(median_filter(np.log1p(filled), size=3, mode='nearest'))
    # scaled to whole blocks and cropped, so each block's value is centered on the block it was measured in (resizing
    # straight to the image size spread ny rows over H and shifted them by up to ~5 % of H; first re-check, 2026-10-29)
    tmap = cv2.resize(smooth.astype(np.float32), (nx * b, ny * b), interpolation=cv2.INTER_LINEAR)[:H, :W]
    stats = {'blocks': int(grid.size), 'blocks_with_reference_pct': float(100 * own.mean()),
             'threshold_raw_p10': float(np.percentile(grid[own], 10)),
             'threshold_raw_median': float(np.median(grid[own])),
             'threshold_raw_p90': float(np.percentile(grid[own], 90)),
             'threshold_raw_whole_image': float(np.percentile(sm[ref], REFERENCE_QUANTILE))}
    return tmap, stats


def nucleus_reference(path: Path, shp: tuple, tissue: np.ndarray) -> tuple[np.ndarray, dict]:
    """Raster union of current nucleus polygons, with no classification/measurement reads.

    Empty/zero-area source nuclei are counted and skipped; source objects are never changed.
    Each polygon is drawn separately before union, so its holes cannot erase another nucleus.
    """
    raw = json.loads(path.read_text())
    feats = raw['features'] if isinstance(raw, dict) else raw
    ref = np.zeros(shp, np.uint8)
    used = absent = degenerate = 0
    for feat in feats:
        geom = feat.get('nucleusGeometry')
        if not geom or geom.get('type') not in ('Polygon', 'MultiPolygon'):
            absent += 1
            continue
        polys = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
        good = False
        for poly in polys:
            if not poly or len(poly[0]) < 3:
                continue
            rings = [np.asarray(r, np.float64) for r in poly if len(r) >= 3]
            if any(r.ndim != 2 or r.shape[1] != 2 or not np.isfinite(r).all() for r in rings):
                raise ValueError(f'Invalid nonfinite nucleus coordinates in {path}')
            if cv2.contourArea(rings[0].astype(np.float32)) <= 0:
                continue
            x0, y0 = np.maximum(np.floor(rings[0].min(axis=0)).astype(int) - 1, 0)
            x1, y1 = np.minimum(np.ceil(rings[0].max(axis=0)).astype(int) + 2, [shp[1], shp[0]])
            if x1 <= x0 or y1 <= y0:
                continue
            tile = np.zeros((y1-y0, x1-x0), np.uint8)
            coords = [np.round((r-[x0,y0])*16).astype(np.int32) for r in rings]
            cv2.fillPoly(tile, [coords[0]], 1, shift=4)
            if len(coords) > 1:
                cv2.fillPoly(tile, coords[1:], 0, shift=4)
            ref[y0:y1,x0:x1] |= tile
            good = True
        used += int(good)
        degenerate += int(not good)
    before = int(ref.sum())
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*NUCLEAR_EROSION_PX+1,)*2)
    ref = cv2.erode(ref, kernel, borderType=cv2.BORDER_CONSTANT, borderValue=0).astype(bool)
    ref &= tissue
    return ref, {'nuclei_used': used, 'missing_nucleus': absent, 'degenerate_nucleus_skipped': degenerate,
                 'nuclear_union_pixels_before_erosion': before, 'reference_pixels': int(ref.sum()),
                 'reference_tissue_pct': float(100*ref.sum()/max(tissue.sum(),1)),
                 'erosion_native_px': NUCLEAR_EROSION_PX, 'erosion_um': NUCLEAR_EROSION_PX*CT.PIXEL_UM,
                 'erosion': 'union of current nucleus polygons; OpenCV elliptical erosion; Tissue intersection'}


def nuclear_threshold_map(sm: np.ndarray, ref: np.ndarray, tissue: np.ndarray) -> tuple[np.ndarray, dict]:
    """The calibrated one-lower-side-scale rule; unsupported blocks never get an invented intensity floor."""
    if sm.shape != ref.shape or sm.shape != tissue.shape or not tissue.any():
        raise ValueError('ColIV arrays must have matching shapes and nonempty Tissue')
    H, W = sm.shape
    b = int(round(BLOCK_UM/CT.PIXEL_UM))
    ny, nx = -(-H//b), -(-W//b)
    grid = np.full((ny,nx), np.nan)
    counts = np.zeros((ny,nx), np.int64)
    in_tissue = np.zeros((ny,nx), bool)
    diagnostics = []
    zero_scale = 0
    for i in range(ny):
        for j in range(nx):
            sl = np.s_[i*b:(i+1)*b,j*b:(j+1)*b]
            r = ref[sl] & tissue[sl] & np.isfinite(sm[sl]) & (sm[sl] >= 0)
            counts[i,j] = int(r.sum())
            in_tissue[i,j] = bool(tissue[sl].any())
            rec = {'row':i, 'column':j, 'reference_pixels':int(counts[i,j]), 'in_tissue':bool(in_tissue[i,j])}
            if counts[i,j] < max(NUCLEAR_MIN_PIXELS, NUCLEAR_MIN_FRACTION*r.size):
                rec['status'] = 'insufficient_reference'
            else:
                values = np.log1p(sm[sl][r])
                center = float(np.median(values))
                scale = float(np.median(center-values[values<=center])/LOWER_MAD_NORMALIZER)
                if not np.isfinite(scale) or scale <= 0:
                    rec['status'] = 'degenerate_lower_scale'
                    zero_scale += 1
                else:
                    level = center+scale
                    if not np.isfinite(np.expm1(level)) or np.expm1(level) > np.finfo(np.float32).max:
                        raise ValueError('ColIV local threshold overflows; reference distribution is unsupported')
                    grid[i,j] = level
                    rec.update(status='own_reference', median_log=center, lower_scale_log=scale, threshold_log=level)
            diagnostics.append(rec)
    own = np.isfinite(grid)
    if not own.any():
        raise ValueError('ColIV: no block with sufficient nuclear reference and nondegenerate lower-side scale')
    distances, near = distance_transform_edt(~own, return_distances=True, return_indices=True)
    smooth = median_filter(grid[near[0],near[1]], size=3, mode='nearest')
    # Full block sizes preserve their original center positions even for partial image-edge blocks.
    tmap = cv2.resize(np.expm1(smooth).astype(np.float32), (nx*b,ny*b), interpolation=cv2.INTER_LINEAR)[:H,:W]
    for rec in diagnostics:
        i,j = rec['row'],rec['column']
        rec.update(nearest_valid_row=int(near[0,i,j]), nearest_valid_column=int(near[1,i,j]),
                   nearest_valid_distance_blocks=float(distances[i,j]))
    valid_ref = ref & tissue & np.isfinite(sm) & (sm >= 0)
    residual = float(100*np.count_nonzero((sm >= tmap)&valid_ref)/max(valid_ref.sum(),1))
    unsupported = int((in_tissue & ~own).sum())
    max_distance = float(distances[in_tissue].max())
    flags = ['reference_residual_is_not_independent_specificity', 'weak_signal_background_not_independently_identifiable']
    if residual > 25:
        flags.append('high_reference_positive_residual_gt_25pct')
    if unsupported:
        flags.append('unsupported_tissue_blocks_interpolated')
    if max_distance > 1:
        flags.append('reference_borrowing_beyond_one_block')
    if zero_scale:
        flags.append('degenerate_reference_blocks')
    return tmap, {'method':METHOD_VERSION, 'reference':'current_image_eroded_nucleus_union',
                  'mesh_is_target_signal':True, 'specificity_claim':False,
                  'block_um':BLOCK_UM, 'block_native_px':b, 'grid_shape':[ny,nx],
                  'minimum_reference_pixels':NUCLEAR_MIN_PIXELS, 'minimum_reference_block_fraction':NUCLEAR_MIN_FRACTION,
                  'minimum_reference_rule':'maximum of pixel and fractional requirements',
                  'lower_scale_divisor':LOWER_MAD_NORMALIZER, 'lower_scale_multiplier':1.0,
                  'tissue_blocks':int(in_tissue.sum()), 'own_reference_blocks':int(own.sum()),
                  'own_reference_tissue_blocks_pct':float(100*own.sum()/in_tissue.sum()),
                  'unsupported_tissue_blocks':unsupported, 'degenerate_scale_blocks':zero_scale,
                  'reference_pixel_count_min_median_max_valid_blocks':[int(counts[own].min()),float(np.median(counts[own])),int(counts[own].max())],
                  'nearest_reference_distance_blocks_p90_max_in_tissue':[float(np.percentile(distances[in_tissue],90)),max_distance],
                  'threshold_raw_p10_median_p90_in_tissue':np.percentile(tmap[tissue],[10,50,90]).tolist(),
                  'threshold_reference_positive_pct':residual, 'finite_nonnegative_reference_pixels':int(valid_ref.sum()),
                  'interpolation':'nearest valid block; 3x3 median log-grid; raw bilinear full blocks then crop',
                  'qc_flags':flags, 'block_diagnostics':diagnostics}


def reference(run: Path, tag: str, shp, tissue: np.ndarray) -> np.ndarray:
    """Pixels of the image's Exocrine cells at least REFERENCE_CLEARANCE_UM from every structure and endocrine object."""
    iso = int(round(REFERENCE_CLEARANCE_UM / CT.PIXEL_UM))
    feats = json.loads((run / 'cell_geojson' / f'{tag}.geojson').read_text())
    feats = feats['features'] if isinstance(feats, dict) else feats
    mask, _, oids, _ = RS.rasterize(feats, shp)
    del feats
    ident = pd.read_parquet(run / 'typing' / 'identity' / f'{tag}.parquet', columns=['object_id', 'label'])
    exocrine = set(ident.loc[ident['label'] == 'Exocrine', 'object_id'])
    keep = np.zeros(len(oids) + 1, bool)
    keep[[i + 1 for i, o in enumerate(oids) if o in exocrine]] = True
    near = np.zeros(shp, np.uint8)
    for cls in REFERENCE_STRUCTURES:
        near |= ST.object_mask(run / 'structures' / cls / f'{cls}__{tag}.geojson', shp).astype(np.uint8)
    raw = json.loads((run / 'islets' / f'{tag}.geojson').read_text())
    for f in raw['features'] if isinstance(raw, dict) else raw:
        if f.get('geometry'):
            ST.fill(near, f['geometry'])
    near = cv2.dilate(near, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * iso + 1, 2 * iso + 1))).astype(bool)
    return keep[mask] & ~near & tissue


def background_level(sm: np.ndarray, tissue: np.ndarray) -> tuple[float, dict]:
    """Superseded histogram diagnostic, retained for reproducibility; not called by either CLI method.

    Owner 2026-10-29: acinar mesh is signal, not a negative reference. Use the established
    background-peak estimator on a systematic native-pixel sample, without cell labels or
    cross-donor pooling. The one-width choice is an agent calibration: earlier 1/2/3-width
    trials showed the one-width mask includes the mesh; larger cuts discard it and vessels.
    An unreadable peak is a method failure, never a silent fallback or an empty mask.
    """
    values = sm[tissue & np.isfinite(sm) & (sm > 0)][::8]
    center, width, level = LV.peak_level(np.log1p(values), 1.0, lowest=True)
    if not np.isfinite(level):
        raise ValueError('ColIV background peak/width unreadable; calibrate method before adoption')
    raw = float(np.expm1(level))
    return raw, {'method': 'all_tissue_peak_low_plus_1_width', 'reference': 'image_tissue_pixels',
                 'reference_pixels': int(len(values)), 'background_peak_log': center,
                 'background_width_log': width, 'level_log': level, 'threshold_raw': raw,
                 'mesh_is_target_signal': True}


def clip_polygons(geoms: list, tissue_geometry, minimum_area: float) -> list:
    """Exact tissue clipping, retaining only polygonal positive-area components."""
    shapely.prepare(tissue_geometry)
    pending = [g if tissue_geometry.covers(g) else g.intersection(tissue_geometry) for g in geoms]
    result = []
    while pending:
        g = pending.pop()
        if g.geom_type == 'Polygon':
            if g.area >= minimum_area:
                result.append(g)
        elif hasattr(g, 'geoms'):
            pending.extend(g.geoms)
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


LOADED_CODE_SHA256 = {str(Path(p).resolve()):sha256(Path(p))
                      for p in (__file__,ST.__file__,RS.__file__,CT.__file__)}


def loaded_code_fingerprint() -> dict:
    current = {path:sha256(Path(path)) for path in LOADED_CODE_SHA256}
    if current != LOADED_CODE_SHA256:
        raise ValueError('ColIV code changed after module import; restart this process before computing or validating')
    return current


def source_paths(run: Path, tag: str, tissue_dir=None, cell_dir=None) -> tuple[Path, Path]:
    return ((Path(tissue_dir) if tissue_dir else run/'inputs'/'tissue')/f'{tag}.geojson',
            (Path(cell_dir) if cell_dir else run/'cell_geojson')/f'{tag}.geojson')


def provenance(run: Path, tag: str, method='all-tissue', tissue_dir=None, cell_dir=None) -> dict:
    code = loaded_code_fingerprint()
    tissue_path, cells_path = source_paths(run, tag, tissue_dir, cell_dir)
    entry = CT.manifest(run)[tag]
    files = [tissue_path, cells_path, run/'project_manifest.json',
             run/'qupath'/'data'/str(entry['entry_id'])/'server.json',
             run/'structures'/'Capillary'/f'Capillary__{tag}.geojson']
    previous = run/'structures'/PREVIOUS/f'{CLASS}__{tag}.geojson'
    if previous.exists():
        files.append(previous)
    if method == 'above-mesh':
        files += [run/'typing'/'identity'/f'{tag}.parquet',run/'islets'/f'{tag}.geojson']
        files += [run/'structures'/cls/f'{cls}__{tag}.geojson' for cls in REFERENCE_STRUCTURES]
    image_path = Path(entry['image_path']).resolve()
    st = image_path.stat()
    return {'schema':'uniform.collagen_iv/2', 'image':tag, 'selected_method':method,
            'method_version':METHOD_VERSION if method == 'all-tissue' else 'historical_above_mesh/1',
            'tissue_path':str(tissue_path.resolve()), 'nuclei_path':str(cells_path.resolve()),
            'files':{str(p.resolve()):sha256(p) for p in files},
            'raw_image':{'path':str(image_path),'size':st.st_size,'mtime_ns':st.st_mtime_ns},
            'raw_image_freshness':'resolved path, size and nanosecond mtime; image bytes are not rehashed',
            'code':code,
            'parameters':{'pixel_um':CT.PIXEL_UM,'sigma_px':SIGMA,'minimum_area_um2':MIN_AREA_UM2,
                          'block_um':BLOCK_UM,'nuclear_erosion_px':NUCLEAR_EROSION_PX,
                          'nuclear_min_pixels':NUCLEAR_MIN_PIXELS,'nuclear_min_fraction':NUCLEAR_MIN_FRACTION,
                          'lower_scale_divisor':LOWER_MAD_NORMALIZER,'despeckle':RS.DESPECKLE}}


def checkpoint_path(run: Path, tag: str) -> Path:
    return run/'structures'/CLASS/f'{CLASS}__{tag}.json'


def atomic_json(path: Path, value, indent=1) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name+'.part')
    part.write_text(json.dumps(value, indent=indent, allow_nan=False)+'\n')
    os.replace(part,path)


def validate_record(run: Path, tag: str, record: dict, method='all-tissue', tissue_dir=None, cell_dir=None,
                    check_renders=False) -> dict:
    if record.get('provenance') != provenance(run, tag, method, tissue_dir, cell_dir):
        raise ValueError(f'Stale Collagen IV inputs or method: {tag}; regenerate mask')
    output = run/'structures'/CLASS/f'{CLASS}__{tag}.geojson'
    if not output.is_file() or record.get('output_sha256') != sha256(output):
        raise ValueError(f'Stale Collagen IV output: {tag}')
    if check_renders:
        for path,digest in record.get('renders_sha256',{}).items():
            if not Path(path).is_file() or sha256(Path(path)) != digest:
                raise ValueError(f'Stale Collagen IV review render: {path}')
    return record


def validate_current(run: Path, images: list[str] | None = None) -> dict:
    """Require current run Tissue/nuclei, not an isolated precomputation's alternative paths."""
    run = Path(run)
    path = run/'structures'/'collagen_iv.json'
    if not path.is_file():
        raise ValueError('Missing Collagen IV manifest; generate all-tissue masks first')
    record = json.loads(path.read_text())
    if record.get('selected_method') != 'all-tissue' or record.get('method_version') != METHOD_VERSION:
        raise ValueError('Stale Collagen IV method; regenerate all-tissue masks')
    tags = images if images is not None else list(record.get('images',{}))
    if not tags:
        raise ValueError('Empty Collagen IV image manifest')
    for tag in tags:
        if tag not in record.get('images',{}):
            raise ValueError(f'Missing Collagen IV image: {tag}')
        validate_record(run,tag,record['images'][tag])
    return record


def one(run: Path, tag: str, renders: Path = RENDERS, method: str = 'all-tissue',
        tissue_dir: Path | None = None, cell_dir: Path | None = None) -> dict:
    started = time.monotonic()
    initial = provenance(run,tag,method,tissue_dir,cell_dir)
    tissue_path, cells_path = source_paths(run,tag,tissue_dir,cell_dir)
    col = load(run, tag, MARKER)
    shp = col.shape
    raw_tissue = json.loads(tissue_path.read_text())
    tissue_features = raw_tissue['features'] if isinstance(raw_tissue,dict) else raw_tissue
    tissue_features = [f for f in tissue_features if
                       ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue']
    tissue_geometry = shapely.union_all([shapely.geometry.shape(f['geometry']) for f in tissue_features])
    if tissue_geometry.is_empty or not tissue_geometry.is_valid or tissue_geometry.area <= 0:
        raise ValueError(f'{tag}: missing, invalid or empty Tissue')
    tissue = np.zeros(shp,np.uint8)
    ST.fill(tissue,mapping(tissue_geometry))
    tissue = tissue.astype(bool)
    sm = gaussian_filter(col, SIGMA)
    if method == 'all-tissue':
        ref, rstats = nucleus_reference(cells_path,shp,tissue)
        tmap, tstats = nuclear_threshold_map(sm,ref,tissue)
        tstats.update(rstats)
    else:
        ref = reference(run, tag, shp, tissue)
        if ref.sum() < 10000:
            raise ValueError(f'{tag}: {int(ref.sum())} reference pixels')
        tmap, tstats = threshold_map(sm, ref)
        tstats.update(method='historical_above_mesh', reference_pixels=int(ref.sum()))
    binary = (sm >= tmap) & tissue
    threshold_tissue_pct = float(100*binary.sum()/tissue.sum())
    del sm, tmap
    geoms = ST.polygons(binary, MIN_AREA_UM2 / CT.PIXEL_UM ** 2)
    geoms = clip_polygons(geoms, tissue_geometry, MIN_AREA_UM2 / CT.PIXEL_UM ** 2)
    feats = [{'type': 'Feature', 'geometry': mapping(g),
              'properties': {'objectType': 'annotation', 'classification': {'name': CLASS}}} for g in geoms]
    dest = run / 'structures' / CLASS / f'{CLASS}__{tag}.geojson'
    if not dest.parent.resolve().is_relative_to(run.resolve()):
        raise ValueError(f'Refusing Collagen IV output through a symlink outside run: {dest.parent}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    drawn = np.zeros(shp, np.uint8)
    for f in feats:
        ST.fill(drawn, f['geometry'])
    drawn = drawn.astype(bool)
    capillary = ST.object_mask(run / 'structures' / 'Capillary' / f'Capillary__{tag}.geojson', shp)
    near = cv2.dilate(drawn.astype(np.uint8), np.ones((7, 7), np.uint8)).astype(bool)   # within ~1 um
    areas = np.array([g.area for g in geoms]) * CT.PIXEL_UM ** 2
    rec = {**tstats, 'objects': len(feats), 'tissue_pct': float(100 * (drawn & tissue).sum() / tissue.sum()),
           'threshold_raster_tissue_pct':threshold_tissue_pct,
           'filtered_reference_positive_pct':float(100*np.count_nonzero(drawn & ref)/max(ref.sum(),1)),
           'exact_object_area_um2':float(areas.sum()), 'exact_tissue_area_um2':float(tissue_geometry.area*CT.PIXEL_UM**2),
           'object_area_um2_median': float(np.median(areas)) if len(areas) else None,
           'object_area_um2_p99': float(np.percentile(areas, 99)) if len(areas) else None,
           'largest_object_pct_of_objects': float(100 * areas.max() / areas.sum()) if len(areas) else None,
           'capillary_pixels_within_1um_pct': float(100 * (capillary & near).sum() / max(capillary.sum(), 1)),
           'objects_per_mm2': float(len(feats) / (tissue.sum() * CT.PIXEL_UM ** 2 / 1e6))}
    old = run / 'structures' / PREVIOUS / f'{CLASS}__{tag}.geojson'
    rec['tissue_sha256'] = sha256(tissue_path)
    rec['nuclei_sha256'] = sha256(cells_path)
    rec['code_sha256'] = sha256(Path(__file__))
    rec['edge_touching_retained'] = True
    render(tag, col, load(run, tag, 'DAPI'), tissue, drawn, capillary,
           ST.object_mask(old, shp) if old.exists() else None, renders)
    if provenance(run,tag,method,tissue_dir,cell_dir) != initial:
        raise ValueError(f'{tag}: Collagen IV inputs changed during processing; output not published')
    atomic_json(dest, {'type':'FeatureCollection','features':feats},indent=None)
    rec['output_sha256'] = sha256(dest)
    rec['provenance'] = initial
    rec['renders_sha256'] = {str(p.resolve()):sha256(p) for p in sorted(renders.glob(f'{tag}_*.png'))}
    rec['renders_sha256'][str((renders/f'{tag}_tiles.json').resolve())] = sha256(renders/f'{tag}_tiles.json')
    rec['elapsed_seconds'] = round(time.monotonic()-started,2)
    rec['review_status'] = 'requires_whole_and_native_tile_review; no independent accuracy claim'
    atomic_json(checkpoint_path(run,tag),rec)
    print(tag,json.dumps({k:rec.get(k) for k in ('objects','tissue_pct','threshold_reference_positive_pct',
          'filtered_reference_positive_pct','unsupported_tissue_blocks','qc_flags','elapsed_seconds')}),flush=True)
    return rec


def run_or_resume(run, tag, renders, method, tissue_dir, cell_dir, resume):
    path = checkpoint_path(run,tag)
    if resume and path.is_file():
        try:
            record = validate_record(run,tag,json.loads(path.read_text()),method,tissue_dir,cell_dir,True)
            if str((renders/f'{tag}_whole.png').resolve()) not in record.get('renders_sha256',{}):
                raise ValueError('render directory changed')
            print(f'{tag}: fresh checkpoint retained',flush=True)
            return record
        except (ValueError,OSError) as exc:
            print(f'{tag}: rebuilding checkpoint ({exc})',flush=True)
    return one(run,tag,renders,method,tissue_dir,cell_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--images', nargs='*', default=[])
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--render-dir', type=Path, default=RENDERS)
    ap.add_argument('--method', choices=['all-tissue', 'above-mesh'], default='all-tissue',
                    help='all tissue collagen (current); above-mesh only reproduces historical output')
    ap.add_argument('--tissue-dir',type=Path,help='alternative accepted Tissue directory for isolated calibration')
    ap.add_argument('--cell-dir',type=Path,help='alternative current nucleus GeoJSON directory for isolated calibration')
    ap.add_argument('--resume',action='store_true',help='reuse only checkpoints with matching inputs, code, output and renders')
    ap.add_argument('--verify-only',action='store_true',help='validate selected image checkpoints without computing')
    a = ap.parse_args()
    tags = a.images or [r['image_name'] for r in CT.cohort_rows()]
    if len(tags) != len(set(tags)) or a.workers < 1:
        ap.error('images must be unique and workers must be positive')
    if a.verify_only:
        for tag in tags:
            record=json.loads(checkpoint_path(a.run,tag).read_text())
            validate_record(a.run,tag,record,a.method,a.tissue_dir,a.cell_dir,True)
        print(f'{len(tags)} Collagen IV checkpoints are fresh')
        return 0
    path = a.run/'structures'/'collagen_iv.json'
    prior = json.loads(path.read_text()) if path.exists() else {}
    images = {}
    for tag,record in prior.get('images',{}).items():
        if tag in tags:
            continue
        try:
            validate_record(a.run,tag,record,a.method,a.tissue_dir,a.cell_dir)
            images[tag]=record
        except (ValueError,OSError):
            pass  # never mix stale previous-code images into a partial current manifest
    def publish():
        atomic_json(path,{'selected_method':a.method,
                         'method_version':METHOD_VERSION if a.method=='all-tissue' else 'historical_above_mesh/1',
                         'sigma_px':SIGMA,'min_area_um2':MIN_AREA_UM2,'images':images})
    with ProcessPoolExecutor(a.workers) as ex:
        for tag,record in zip(tags,ex.map(run_or_resume,repeat(a.run),tags,repeat(a.render_dir),repeat(a.method),
                                         repeat(a.tissue_dir),repeat(a.cell_dir),repeat(a.resume))):
            images[tag]=record
            publish()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
