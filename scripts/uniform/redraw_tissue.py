#!/usr/bin/env python3
"""Redraw the Tissue outline of the five added donors' images (owner, 2026-10-27).

The additions' Tissue outlines came from the `Tissue` pixel classifier of the six-donor review project and were reused
unchanged; ND_6232's covered 38 % of the tissue visible in DAPI, and ND_6539 / Aab_6558 held thousands of cell-sized
holes (CLAUDE.md, 2026-10-26 standing rule). The redraw, approved by the owner "if all holes are filled":

  - pyramid level 3 (3.0192 um/px, the original Tissue classifier's resolution);
  - nuclei: DAPI, Gaussian sigma 2 px, log; Otsu on the nonzero pixels, threshold moved to 0.6 of the way from the
    5th percentile to Otsu's;
  - stroma: ColIV, VIM and SMA, each log-scaled to its own 1st-99th percentiles in the image, averaged, Gaussian
    sigma 2 px; Otsu threshold at 1.2 of the way from the 5th percentile;
  - Tissue = nuclei OR stroma; pieces under 80 um2 removed (the original 15 outlines' smallest piece is 82 um2); every
    hole filled (the original 15 keep holes of 1,000 um2 and more; the owner's rule for the redraw fills all);
  - clipped to the image less a 5 um border, as the original Tissue was (`edgeMarginUm` and `clipObjectsToInner` in
    the original project's 1_Panc_Objects.groovy, "to fully cover the bright rim").

Without the hole filling the same rule reproduces the original 15 outlines at Dice 0.958-0.991.

Owner 2026-10-29, after the image check of the hole-filled redraw: one content rule for all 20 images ("Content rule,
all 20"; `--content-rule`). Filling every hole made inclusion depend on shape, not content: nucleus-free gaps between
lobules became Tissue when enclosed and stayed out when they reached the image edge (T1D_6533: 40 % of its dark septal
space added, 58 % left out; Aab_6450: 4.4 % of the added pixels held nuclei or stroma). The content rule: nuclei or
stroma, holes under 1,000 um2 filled, pieces under 80 um2 removed, no other hole filled, plus the footprint of every
detected cell and endocrine object of the run (`--keep-from`, dilated by one level-3 pixel), so nothing detected leaves
Tissue; empty gaps and lumens stay out wherever they are. The five additions are redrawn the same way.
After conversion to polygons, union any missing parts of the ORIGINAL cell, nucleus and endocrine geometries.
The 5 um crop applies to newly inferred tissue only; an existing footprint at the image edge is preserved exactly.

2026-10-29 continuation, image-checked refinement: absolute DAPI/stroma thresholds alone admitted smooth lumen haze
(T1D_6533 and T1D_6550). The common content rule additionally requires contrast above a 51 um square opening:
log-DAPI residual > 0.6 times its positive-residual Otsu threshold, or stromal residual > its positive-residual Otsu
threshold. The brightest 1% of the stromal composite remains eligible without contrast, preserving thick bright
vessel walls. Original footprints and the existing small-hole/piece cleanup stay in force. The same parameters are
used on all 20 images; per-image data-derived thresholds are recorded. Calibration comparisons live under
output/typing_resume_20261029/tissue_strong_rescue_trial/. The 101 um trial failed the diagnosed lumens; contrast
without the strong-stroma safeguard removed bright wall material.

The final cleanup preserves bounded empty components that survived the first 1,000 um2 hole cleanup, before
original-footprint padding. Padding can split or shrink a large lumen below the small-hole threshold; its remaining
unoccupied pixels are therefore restored after the second cleanup. Padded footprints and exact original geometry
remain protected. This corrects the ND_6516 lumen diagnosed in the continuation without changing signal gates or
reopening small holes that were already filled by the content cleanup.

Superseded the same day: the same rule for all 20 images ("Fill every hole, all 20"). The original 15 images' outlines kept
every hole of 1,000 um2 or more, cutting out septal stroma, vessel walls and periductal tissue (0.6-4.0 % of their area
held nuclei, never segmented). Their redraw is merged with the outline in use (`--union-with`), every hole filled:
alone, the rule leaves out 0.15-3.1 % of their area, loose interlobular connective tissue with sparse dim nuclei where
60-712 cells per image were detected and 13 endocrine objects lie (settled, owner 2026-10-28), so nothing detected
leaves Tissue (output/tissue_redraw_check_20261029/).

Writes <out-dir>/<image>.geojson: one Tissue annotation (MultiPolygon, full-resolution pixel coordinates) that
run_uniform.py's tissue_export step copies for images whose cohort.csv `tissue_outline` names it.

    python scripts/uniform/redraw_tissue.py --images Aab_6090 ND_6232 T1D_6520 ND_6539 Aab_6558 \\
        --out-dir scripts/uniform/tissue_outlines
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import cv2
import numpy as np
import shapely
import tifffile
from scipy import ndimage
from shapely.geometry import MultiPolygon, Polygon, box, mapping, shape
from shapely.ops import unary_union
from shapely.validation import make_valid
from skimage import filters, morphology

REPO = Path(__file__).resolve().parents[2]
IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
LEVEL, DOWNSAMPLE = 3, 8
PIXEL_UM = 0.3774
PX2 = (PIXEL_UM * DOWNSAMPLE) ** 2
STROMA = ('ColIV', 'VIM', 'SMA')
MIN_PIECE_UM2 = 80
EDGE_MARGIN_UM = 5.0
CONTENT_CONTRAST_UM = 51.0
DAPI_CONTRAST_FACTOR = 0.6
STRONG_STROMA_PERCENTILE = 99.0
CONTENT_RULE_VERSION = 'content-local-contrast-51um-strong-stroma-protected-holes-v2'
SMALL_HOLE_UM2 = 1000.0


def channels(path: Path, level_index: int = LEVEL) -> tuple[np.ndarray, np.ndarray, tuple[int, int]]:
    with tifffile.TiffFile(path) as t:
        height, width = t.series[0].levels[0].shape[-2:]
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', t.ome_metadata)
        level = t.series[0].levels[level_index]
        downsample = 2 ** level_index
        if level.shape[-2:] != tuple(-(-s // downsample) for s in t.series[0].levels[0].shape[-2:]) and \
                level.shape[-2:] != tuple(s // downsample for s in t.series[0].levels[0].shape[-2:]):
            raise SystemExit(f'{path}: pyramid level {level_index} is not a {downsample}x downsample')
        dapi = level.pages[names.index('DAPI')].asarray()
        stroma = []
        for n in STROMA:
            v = np.log1p(level.pages[names.index(n)].asarray().astype(np.float32))
            lo, hi = np.percentile(v, [1, 99])
            stroma.append(np.clip((v - lo) / (hi - lo + 1e-6), 0, 1))
    return dapi, np.mean(stroma, axis=0), (width, height)


def above(v: np.ndarray, fg: np.ndarray, k: float) -> np.ndarray:
    t, bg = filters.threshold_otsu(fg), np.percentile(fg, 5)
    return v > bg + k * (t - bg)


def positive_otsu(values: np.ndarray) -> float:
    positive = values[values > 0]
    return float(filters.threshold_otsu(positive)) if positive.size else float('inf')


def content_mask(dapi: np.ndarray, stroma: np.ndarray, diagnostics: dict | None = None) -> np.ndarray:
    """Locally contrasting nuclei/stroma or strong stroma, with common physical-scale parameters."""
    d = np.log1p(filters.gaussian(dapi.astype(float), 2, preserve_range=True))
    s = filters.gaussian(stroma, 2)
    nuclei, stromal = above(d, d[d > 0], 0.6), above(s, s, 1.2)
    width = round(CONTENT_CONTRAST_UM / (PIXEL_UM * DOWNSAMPLE))
    width += 1 - width % 2
    d_contrast = d - ndimage.grey_opening(d, size=(width, width))
    s_contrast = s - ndimage.grey_opening(s, size=(width, width))
    d_threshold = DAPI_CONTRAST_FACTOR * positive_otsu(d_contrast)
    s_threshold = positive_otsu(s_contrast)
    strong_threshold = float(np.percentile(s, STRONG_STROMA_PERCENTILE))
    m = (nuclei & (d_contrast > d_threshold)) | (stromal & ((s_contrast > s_threshold) | (s >= strong_threshold)))
    if diagnostics is not None:
        diagnostics.update(rule_version=CONTENT_RULE_VERSION, pyramid_level=LEVEL,
                           contrast_opening_requested_um=CONTENT_CONTRAST_UM, contrast_opening_pixels=width,
                           contrast_opening_actual_um=width * PIXEL_UM * DOWNSAMPLE,
                           dapi_contrast_factor=DAPI_CONTRAST_FACTOR,
                           dapi_contrast_threshold=d_threshold if np.isfinite(d_threshold) else None,
                           stroma_contrast_threshold=s_threshold if np.isfinite(s_threshold) else None,
                           strong_stroma_percentile=STRONG_STROMA_PERCENTILE,
                           strong_stroma_threshold=strong_threshold)
    m = morphology.remove_small_holes(m, int(SMALL_HOLE_UM2 / PX2))
    return morphology.remove_small_objects(m, int(MIN_PIECE_UM2 / PX2))


def add_footprints_preserving_content_holes(content: np.ndarray, footprint: np.ndarray,
                                          diagnostics: dict | None = None) -> np.ndarray:
    """Keep large content-supported voids after padding shrinks or divides them.

    ``content`` has already undergone the first small-hole and small-piece cleanup.
    Only bounded pre-padding holes are protected, and footprint pixels always win.
    """
    threshold = int(SMALL_HOLE_UM2 / PX2)
    protected = ndimage.binary_fill_holes(content) & ~content
    labels, count = ndimage.label(protected)
    sizes = np.bincount(labels.ravel())[1:]
    if np.any(sizes < threshold):
        raise ValueError('Content mask contains a bounded hole below the first cleanup threshold')
    current = morphology.remove_small_holes(content | footprint, threshold)
    result = current & ~(protected & ~footprint)
    restored = current & ~result
    if diagnostics is not None:
        diagnostics.update(small_hole_requested_um2=SMALL_HOLE_UM2,
                           small_hole_threshold_pixels=threshold,
                           small_hole_effective_um2=threshold * PX2,
                           protected_pre_padding_holes=int(count),
                           protected_pre_padding_hole_min_um2=float(sizes.min() * PX2) if count else None,
                           protected_pre_padding_hole_total_um2=float(sizes.sum() * PX2),
                           padded_footprint_pixels_in_protected_holes=int(np.count_nonzero(protected & footprint)),
                           restored_empty_pixels_after_second_fill=int(restored.sum()),
                           restored_empty_area_um2=float(restored.sum() * PX2),
                           protected_hole_connectivity=1,
                           protected_hole_rule='Restore bounded holes surviving first content cleanup after second fill, except padded footprints')
    return result


def footprints(run: Path, tag: str, shape: tuple[int, int]) -> np.ndarray:
    """Every detected cell's and endocrine object's polygon at level 3, dilated by one pixel."""
    m = np.zeros(shape, np.uint8)
    for path in (run / 'cell_geojson' / f'{tag}.geojson', run / 'islets' / f'{tag}.geojson'):
        raw = json.loads(path.read_text())
        for f in raw['features'] if isinstance(raw, dict) else raw:
            g = f.get('geometry') or {}
            for poly in (g.get('coordinates') or []) if g.get('type') == 'MultiPolygon' else [g.get('coordinates') or []]:
                if poly and len(poly[0]) >= 3:
                    cv2.fillPoly(m, [np.round(np.asarray(poly[0], float) / DOWNSAMPLE * 16).astype(np.int32)], 1,
                                 lineType=cv2.LINE_8, shift=4)
    return cv2.dilate(m, np.ones((3, 3), np.uint8)).astype(bool)


def polygon_parts(geom):
    """Only polygonal components, including valid repairs nested in collections."""
    if geom.geom_type == 'Polygon':
        return [geom] if not geom.is_empty else []
    return [part for child in getattr(geom, 'geoms', ()) for part in polygon_parts(child)]


def footprint_batches(run: Path, tag: str, batch_size: int = 10000, diagnostics: dict | None = None):
    """Read original geometry without rasterizing or changing its coordinates."""
    for dirname, fields in (('cell_geojson', ('geometry', 'nucleusGeometry')), ('islets', ('geometry',))):
        raw = json.loads((run / dirname / f'{tag}.geojson').read_text())
        batch = []
        kinds = []
        for feature in raw['features'] if isinstance(raw, dict) else raw:
            for field in fields:
                if feature.get(field):
                    geom = shape(feature[field])
                    kind = 'endocrine' if dirname == 'islets' else 'nucleus' if field == 'nucleusGeometry' else 'cell'
                    if geom.is_empty or geom.area <= 0:
                        # Some source nuclei are the sentinel polygon [(0,0), ...].
                        # They have no footprint: retain the source and record it,
                        # without inventing a tissue island at the image origin.
                        if diagnostics is not None:
                            diagnostics.setdefault('empty_source_footprints', []).append(
                                {'object_id': feature.get('id'), 'kind': kind, 'geometry_type': geom.geom_type})
                        continue
                    if not geom.is_valid:
                        geom = make_valid(geom)
                    batch.append(geom)
                    kinds.append(kind)
                    if len(batch) == batch_size:
                        yield batch, kinds
                        batch, kinds = [], []
        if batch:
            yield batch, kinds


def preserve_exact_footprints(geom, batches):
    """Keep original footprints after conversion and the new-tissue edge crop.

    Only missing portions are unioned. This avoids the many-million-polygon
    union and retains the content rule in unoccupied areas and empty lumens.
    """
    shapely.prepare(geom)
    additions = []
    counts = {'cell': 0, 'nucleus': 0, 'endocrine': 0}
    for batch, kinds in batches:
        covered = shapely.covers(geom, batch)
        for i in np.flatnonzero(~covered):
            additions.append(batch[i])
            counts[kinds[i]] += 1
    old_area = geom.area
    if additions:
        geom = unary_union([geom, *additions])
    geom = MultiPolygon(polygon_parts(geom))
    if not geom.is_valid:
        raise ValueError('Exact-footprint union produced invalid tissue')
    counts['added_area_um2'] = (geom.area - old_area) * PIXEL_UM ** 2
    return geom, counts


def tissue_mask(dapi: np.ndarray, stroma: np.ndarray) -> np.ndarray:
    d = np.log1p(filters.gaussian(dapi.astype(float), 2, preserve_range=True))
    s = filters.gaussian(stroma, 2)
    m = above(d, d[d > 0], 0.6) | above(s, s, 1.2)
    # the order the owner reviewed: holes under 1,000 um2 filled, small pieces removed, then every hole filled
    m = morphology.remove_small_holes(m, int(1000 / PX2))
    m = morphology.remove_small_objects(m, int(MIN_PIECE_UM2 / PX2))
    return ndimage.binary_fill_holes(m)


def to_geometry(mask: np.ndarray, size: tuple[int, int], diagnostics: dict | None = None) -> MultiPolygon:
    """Outer contours traced at 4x (OpenCV traces pixel centres: a quarter-pixel error), corner contacts repaired."""
    scale = DOWNSAMPLE / 4
    cs, hier = cv2.findContours(np.kron(mask.astype(np.uint8), np.ones((4, 4), np.uint8)), cv2.RETR_CCOMP,
                                cv2.CHAIN_APPROX_SIMPLE)
    parts = []
    for i, c in enumerate(cs):
        if hier[0][i][3] != -1 or len(c) < 3:
            continue
        holes, j = [], hier[0][i][2]
        while j != -1:
            if len(cs[j]) >= 3:
                holes.append((cs[j][:, 0, :] + 0.5) * scale)
            j = hier[0][j][0]
        # OpenCV contours can touch at corners. The shell/hole hierarchy is
        # known: structure repair preserves it. Default linework repair may
        # promote a touching hole into filled tissue (observed in T1D_6534).
        parts.append(shapely.make_valid(Polygon((c[:, 0, :] + 0.5) * scale, holes), method='structure'))
    geom = unary_union(parts)
    # Validate the inferred mask before the intentional new-tissue edge crop.
    # Pixel centers must survive vectorization exactly; this detects filled
    # large holes even when every original cell remains safely covered.
    if diagnostics is not None:
        shapely.prepare(geom)
        added = removed = 0
        for y0 in range(0, mask.shape[0], 128):
            yy, xx = np.indices(mask[y0:y0 + 128].shape)
            points = shapely.points((xx.ravel() + 0.5) * DOWNSAMPLE,
                                   (yy.ravel() + y0 + 0.5) * DOWNSAMPLE)
            back = shapely.covers(geom, points).reshape(yy.shape)
            expected = mask[y0:y0 + 128]
            added += int(np.count_nonzero(back & ~expected))
            removed += int(np.count_nonzero(expected & ~back))
        diagnostics.update(repair_method='structure', pre_crop_added_pixel_centers=added,
                           pre_crop_removed_pixel_centers=removed, shapely_version=shapely.__version__,
                           geos_version=shapely.geos_version_string, opencv_version=cv2.__version__)
        if added or removed:
            raise ValueError(f'Raster-to-vector conversion changed pixel centers: +{added}, -{removed}')
    m = EDGE_MARGIN_UM / PIXEL_UM
    geom = geom.intersection(box(m, m, size[0] - m, size[1] - m))
    geom = MultiPolygon([p for g in getattr(geom, 'geoms', [geom]) for p in getattr(g, 'geoms', [g])
                         if p.geom_type == 'Polygon'])
    if not geom.is_valid:
        raise SystemExit('outline geometry is not valid')
    return geom


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--content-rule', action='store_true', help='the content rule (owner 2026-10-29); needs --keep-from')
    ap.add_argument('--keep-from', type=Path, help='the run whose detected cells and endocrine objects stay inside Tissue')
    ap.add_argument('--union-with', type=Path,
                    help='directory of the outlines in use (<image>.geojson); each redraw keeps its image\'s outline inside')
    a = ap.parse_args()
    if a.content_rule and a.union_with:
        ap.error('--content-rule cannot be combined with --union-with, which fills every hole')
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for tag in a.images:
        dapi, stroma, size = channels(IMAGES / f'{tag}.ome.tiff')
        content_diagnostics = {}
        if a.content_rule:
            if not a.keep_from:
                raise SystemExit('--content-rule needs --keep-from <run>')
            content = content_mask(dapi, stroma, diagnostics=content_diagnostics)
            mask = add_footprints_preserving_content_holes(
                content, footprints(a.keep_from, tag, dapi.shape), content_diagnostics)
        else:
            mask = tissue_mask(dapi, stroma)
        conversion_diagnostics = {}
        geom = to_geometry(mask, size, conversion_diagnostics if a.content_rule else None)
        if a.content_rule:
            diagnostics = {}
            geom, preservation = preserve_exact_footprints(geom, footprint_batches(a.keep_from, tag, diagnostics=diagnostics))
            preservation.update(diagnostics)
            preservation['content_rule'] = content_diagnostics
            preservation['mask_conversion'] = conversion_diagnostics
            preservation['source_run'] = str(a.keep_from.resolve())
            preservation['generator_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            print(f'{tag}: exact footprint restoration {preservation}', flush=True)
            (a.out_dir / f'{tag}.preservation.json').write_text(json.dumps(preservation, indent=2) + '\n')
        if a.union_with:
            raw = json.loads((a.union_with / f'{tag}.geojson').read_text())
            current = unary_union([make_valid(shape(f['geometry'])) for f in (raw['features'] if isinstance(raw, dict) else raw)
                                   if ((f.get('properties') or {}).get('classification') or {}).get('name') == 'Tissue'])
            merged = unary_union([geom, current])
            parts = [p for g in getattr(merged, 'geoms', [merged]) for p in getattr(g, 'geoms', [g]) if p.geom_type == 'Polygon']
            geom = unary_union([Polygon(p.exterior) for p in parts])          # every hole filled
            geom = MultiPolygon([p for p in getattr(geom, 'geoms', [geom]) if p.geom_type == 'Polygon'])
            if not geom.is_valid or not geom.covers(current.buffer(-1e-6)):
                raise SystemExit(f'{tag}: merged outline is not valid or does not hold the outline in use')
        holes = sum(len(p.interiors) for p in geom.geoms)
        feature = {'type': 'Feature', 'geometry': mapping(geom),
                   # exactly the original exports' properties: an unnamed annotation, so QuPath's cell `Parent` reads
                   # "Annotation (Tissue)", the value the pipeline keys tissue cells on (islet_cell_test.TISSUE_PARENT)
                   'properties': {'objectType': 'annotation', 'classification': {'name': 'Tissue', 'color': [102, 153, 153]},
                                  'isLocked': True}}
        out = a.out_dir / f'{tag}.geojson'
        out.write_text(json.dumps({'type': 'FeatureCollection', 'features': [feature]}))
        print(f'{tag}: {mask.sum() * PX2 / 1e6:.2f} mm2 (mask), {geom.area * PIXEL_UM ** 2 / 1e6:.2f} mm2 (outline), '
              f'{len(geom.geoms)} pieces, {holes} holes -> {out}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
