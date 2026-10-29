#!/usr/bin/env python3
"""Render new augmentation cells for independent geometry/image review, never auto-accept them.

    python scripts/research/augmentation_geometry_review_20261029.py --run <candidate> \\
        --images Aab_6450 --out output/typing_resume_20261029/augmentation_geometry_review

Requires the full augmentation verify checkpoint. Each image gets an exact vector containment
report, a whole-image view, and up to six spatially diverse, seeded new-cell fields at native
OME pyramid level 0. DAPI is blue, BActin green, ECAD red; these three raw channels are a documented
morphology view, not a recreation of the nine-channel InstanSeg input or a phenotype call.
Each native field includes raw DAPI, the raw composite, and geometry overlays. Whole-image
level 3 is for orientation only. At most one full-resolution marker plane is held at a time.
Existing image review directories are never overwritten. A human/agent must inspect the images
and separately record their verdict; geometric checks alone do not constitute acceptance.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
from itertools import repeat
import json
from pathlib import Path
import re
import sys

import cv2
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts/uniform/typing'))
import codex_typing as CT
import refresh_ready_cells as RF

DEFAULT_IMAGES = REPO / 'data/development/uniform-20donor-20261022-v10/images'
DEFAULT_OUT = REPO / 'output/typing_resume_20261029/augmentation_geometry_review'
CHANNELS = ('DAPI', 'BActin', 'ECAD')
COLORS = {'old_cell': (135, 135, 135), 'old_nucleus': (255, 140, 40),
          'new_cell': (0, 255, 255), 'new_nucleus': (255, 255, 0),
          'baseline_tissue': (255, 255, 255), 'final_tissue': (255, 0, 255)}
AREA_EPS_UM2 = 1e-9  # report strict covers separately; ignore only floating-point area residue


def features(path: Path) -> list[dict]:
    raw = json.loads(path.read_text())
    return raw['features'] if isinstance(raw, dict) else raw


def tissue_geometry(path: Path):
    parts = [shape(f['geometry']) for f in features(path)
             if (f.get('properties', {}).get('classification') or {}).get('name') == 'Tissue']
    if not parts or any(not g.is_valid or g.is_empty for g in parts):
        raise ValueError(f'Missing, empty or invalid Tissue: {path}')
    return shapely.union_all(parts)


def cell_geometries(path: Path):
    ids, cells, nuclei, invalid = [], [], [], []
    for f in features(path):
        oid = f['id']
        cell = shape(f['geometry'])
        raw_nuc = f.get('nucleusGeometry')
        nuc = shape(raw_nuc) if raw_nuc and raw_nuc.get('coordinates') else shapely.GeometryCollection()
        if not cell.is_valid or (nuc.area > 0 and not nuc.is_valid):
            invalid.append(oid)
        # Retain exact valid shapes; record repairs needed for GEOS measurement, never hide them.
        ids.append(oid)
        cells.append(cell if cell.is_valid else shapely.make_valid(cell))
        nuclei.append(nuc if nuc.is_valid else shapely.make_valid(nuc))
    if len(ids) != len(set(ids)) or any(not oid for oid in ids):
        raise ValueError('Null/duplicate GeoJSON object IDs')
    return ids, np.array(cells, dtype=object), np.array(nuclei, dtype=object), invalid


def containment(geoms: np.ndarray, container) -> tuple[np.ndarray, np.ndarray]:
    """Strict covers plus exact positive-area loss; degenerate nuclei are explicitly unassessed."""
    strict = shapely.covers(container, geoms)
    area = shapely.area(geoms)
    lost = np.zeros(len(geoms), float)
    need = np.flatnonzero((area > 0) & ~strict)
    if len(need):
        lost[need] = shapely.area(shapely.difference(geoms[need], container)) * CT.PIXEL_UM ** 2
    return strict, lost


def check_geometry(frame: pd.DataFrame, ids: list[str], cells: np.ndarray, nuclei: np.ndarray,
                   new_ids: set[str], baseline, final, invalid: list[str]) -> tuple[dict, pd.DataFrame]:
    if frame['Object ID'].isna().any() or frame['Object ID'].duplicated().any() or set(frame['Object ID']) != set(ids):
        raise ValueError('Cell table and geometry inventory differ or IDs are duplicated/null')
    if new_ids - set(ids):
        raise ValueError('Accepted new IDs are missing from cells/geometry')
    table = frame.set_index('Object ID').loc[ids]
    exported_xy = table[['Centroid X µm', 'Centroid Y µm']].to_numpy(float)
    if not np.isfinite(exported_xy).all():
        raise ValueError('Nonfinite exported cell centroids')
    # QuPath's TSV export rounds centroids; boundary slivers can move across Tissue when
    # that rounded value is used. Match augmentation's ROI-centroid test exactly instead.
    points = shapely.centroid(cells)
    xy = shapely.get_coordinates(points) * CT.PIXEL_UM
    if xy.shape != exported_xy.shape or not np.isfinite(xy).all():
        raise ValueError('Cell geometry lacks a finite centroid')
    added = final.difference(baseline)
    shapely.prepare(added)
    shapely.prepare(final)
    in_added, in_final = shapely.covers(added, points), shapely.covers(final, points)
    exported_points = shapely.points(exported_xy / CT.PIXEL_UM)
    exported_added = shapely.covers(added, exported_points)
    exported_final = shapely.covers(final, exported_points)
    cell_strict, cell_lost = containment(cells, final)
    nuc_strict, nuc_lost = containment(nuclei, final)
    positive_nuc = shapely.area(nuclei) > 0
    nuc_own = shapely.covers(cells, nuclei)
    nuc_own_lost = np.zeros(len(ids))
    check_own = np.flatnonzero(positive_nuc & ~nuc_own)
    if len(check_own):
        nuc_own_lost[check_own] = shapely.area(shapely.difference(nuclei[check_own], cells[check_own])) * CT.PIXEL_UM ** 2
    is_new = np.array([oid in new_ids for oid in ids])
    report = {'new_cells': int(is_new.sum()), 'old_cells': int((~is_new).sum()),
              'added_tissue_area_um2': added.area * CT.PIXEL_UM ** 2,
              'area_residue_tolerance_um2': AREA_EPS_UM2, 'invalid_source_geometry_ids': invalid}
    for name, select in [('new', is_new), ('old', ~is_new)]:
        report.update({
            f'{name}_centroids_outside_final': int((select & ~in_final).sum()),
            f'{name}_cell_strict_covers_failures': int((select & ~cell_strict).sum()),
            f'{name}_cell_positive_area_outside_final': int((select & (cell_lost > AREA_EPS_UM2)).sum()),
            f'{name}_cell_max_area_outside_final_um2': float(cell_lost[select].max(initial=0)),
            f'{name}_positive_area_nuclei': int((select & positive_nuc).sum()),
            f'{name}_empty_or_zero_area_nuclei': int((select & ~positive_nuc).sum()),
            f'{name}_nucleus_strict_covers_failures': int((select & positive_nuc & ~nuc_strict).sum()),
            f'{name}_nucleus_positive_area_outside_final': int((select & (nuc_lost > AREA_EPS_UM2)).sum()),
            f'{name}_nucleus_max_area_outside_final_um2': float(nuc_lost[select].max(initial=0)),
            f'{name}_nucleus_positive_area_outside_own_cell': int((select & (nuc_own_lost > AREA_EPS_UM2)).sum()),
        })
    report['new_centroids_outside_added'] = int((is_new & ~in_added).sum())
    report['centroid_method'] = 'Exact cell GeoJSON geometry centroid, matching augmentation ROI centroids'
    report['new_exported_rounded_centroids_outside_added'] = int((is_new & ~exported_added).sum())
    report['new_exported_rounded_centroids_outside_final'] = int((is_new & ~exported_final).sum())
    report['maximum_exported_centroid_rounding_distance_um'] = float(np.linalg.norm(xy - exported_xy, axis=1).max(initial=0))
    required_zero = ['new_centroids_outside_added', 'new_centroids_outside_final'] + [
        f'{kind}_{part}_positive_area_outside_final' for kind in ('old', 'new') for part in ('cell', 'nucleus')]
    required_zero.append('new_nucleus_positive_area_outside_own_cell')
    report['geometry_checks_passed'] = not invalid and all(report[k] == 0 for k in required_zero)
    report['visual_review_status'] = 'pending'
    detail = pd.DataFrame({'object_id': ids, 'is_new': is_new, 'centroid_x_um': xy[:, 0], 'centroid_y_um': xy[:, 1],
                           'centroid_inside_added': in_added, 'centroid_inside_final': in_final,
                           'exported_centroid_x_um': exported_xy[:, 0], 'exported_centroid_y_um': exported_xy[:, 1],
                           'exported_centroid_inside_added': exported_added, 'exported_centroid_inside_final': exported_final,
                           'cell_area_um2': shapely.area(cells) * CT.PIXEL_UM ** 2,
                           'cell_area_outside_final_um2': cell_lost,
                           'nucleus_area_um2': shapely.area(nuclei) * CT.PIXEL_UM ** 2,
                           'nucleus_area_outside_final_um2': nuc_lost,
                           'nucleus_area_outside_own_cell_um2': nuc_own_lost})
    return report, detail


def choose_fields(detail: pd.DataFrame, seed: int, count: int = 6) -> pd.DataFrame:
    """Seeded first new cell, then farthest-from-selected centroids for spatial diversity."""
    candidates = detail[detail['is_new']].sort_values('object_id').reset_index(drop=True)
    if candidates.empty:
        return candidates
    xy = candidates[['centroid_x_um', 'centroid_y_um']].to_numpy()
    rng = np.random.default_rng(seed)
    selected = [int(rng.integers(len(candidates)))]
    nearest = np.sum((xy - xy[selected[0]]) ** 2, axis=1)
    for _ in range(1, min(count, len(candidates))):
        nearest[selected] = -1
        choice = int(np.argmax(nearest))
        selected.append(choice)
        nearest = np.minimum(nearest, np.sum((xy - xy[choice]) ** 2, axis=1))
    return candidates.iloc[selected]


def rings(geom):
    if geom.is_empty:
        return
    if geom.geom_type == 'Polygon':
        yield np.asarray(geom.exterior.coords)
        for inner in geom.interiors:
            yield np.asarray(inner.coords)
    elif hasattr(geom, 'geoms'):
        for part in geom.geoms:
            yield from rings(part)


def draw_shape(rgb, geom, color, origin=(0, 0), downsample=1, thickness=1):
    lines = [np.round((ring - origin) / downsample * 16).astype(np.int32) for ring in rings(geom)]
    if lines:
        cv2.polylines(rgb, lines, True, color, thickness, cv2.LINE_AA, shift=4)


def overlay(raw, ids, cells, nuclei, new_ids, baseline, final, origin=(0, 0), downsample=1, selected_id=None):
    rgb = raw.copy()
    for oid, cell, nuc in zip(ids, cells, nuclei):
        prefix = 'new' if oid in new_ids else 'old'
        draw_shape(rgb, cell, COLORS[prefix + '_cell'], origin, downsample)
        draw_shape(rgb, nuc, COLORS[prefix + '_nucleus'], origin, downsample)
    draw_shape(rgb, baseline, COLORS['baseline_tissue'], origin, downsample, 2)
    draw_shape(rgb, final, COLORS['final_tissue'], origin, downsample, 2)
    if selected_id in ids:
        index = ids.index(selected_id)
        draw_shape(rgb, cells[index], COLORS['new_cell'], origin, downsample, 2)
        draw_shape(rgb, nuclei[index], COLORS['new_nucleus'], origin, downsample, 2)
        point = cells[index].centroid
        center = tuple(np.round((np.array([point.x, point.y]) - origin) / downsample).astype(int))
        cv2.drawMarker(rgb, center, (255, 255, 255), cv2.MARKER_CROSS, 13, 1, cv2.LINE_AA)
    return rgb


def with_legend(panels, title):
    rgb = np.hstack(panels)
    canvas = cv2.copyMakeBorder(rgb, 96, 28, 0, 0, cv2.BORDER_CONSTANT)
    cv2.putText(canvas, title, (8, 20), cv2.FONT_HERSHEY_SIMPLEX, .48, (255, 255, 255), 1, cv2.LINE_AA)
    cv2.putText(canvas, 'RAW: DAPI blue | BActin green | ECAD red. Native tiles: DAPI gray / composite / overlays.',
                (8, 41), cv2.FONT_HERSHEY_SIMPLEX, .42, (220, 220, 220), 1, cv2.LINE_AA)
    for i, (key, color) in enumerate(COLORS.items()):
        x, y = 8 + (i % 3) * (canvas.shape[1] // 3), 63 + (i // 3) * 22
        cv2.line(canvas, (x, y - 4), (x + 20, y - 4), color, 2)
        cv2.putText(canvas, key.replace('_', ' '), (x + 25, y), cv2.FONT_HERSHEY_SIMPLEX, .42, color, 1, cv2.LINE_AA)
    cv2.putText(canvas, 'White cross: sampled cell centroid. Geometry statistics do not imply visual acceptance; review pending.',
                (8, canvas.shape[0] - 9), cv2.FONT_HERSHEY_SIMPLEX, .42, (220, 220, 220), 1, cv2.LINE_AA)
    return canvas


def scale(values, limits):
    lo, hi = limits
    return (255 * np.clip((np.log1p(values.astype(np.float32)) - lo) / max(hi - lo, 1e-6), 0, 1)).astype(np.uint8)


def native_views(image_path: Path, selected: pd.DataFrame, tile_um: float, whole_level: int):
    """Read one marker plane at a time, copying only the six native crops before releasing it."""
    with tifffile.TiffFile(image_path) as tf:
        names = re.findall(r'<(?:\w+:)?Channel [^>]*?Name="([^"]*)"', tf.ome_metadata)
        missing = set(CHANNELS) - set(names)
        if missing:
            raise ValueError(f'Missing morphology channels: {sorted(missing)}')
        series = tf.series[0]
        if whole_level >= len(series.levels):
            raise ValueError('Requested whole-image pyramid level does not exist')
        height, width = series.shape[-2:]
        side = min(round(tile_um / CT.PIXEL_UM), height, width)
        level = series.levels[whole_level]
        downsample = 2 ** whole_level
        if any(abs(a * downsample - b) >= downsample for a, b in zip(level.shape[-2:], (height, width))):
            raise ValueError('OME pyramid has an unexpected downsample')
        windows = []
        for row in selected.itertuples():
            x = int(np.clip(round(row.centroid_x_um / CT.PIXEL_UM) - side // 2, 0, width - side))
            y = int(np.clip(round(row.centroid_y_um / CT.PIXEL_UM) - side // 2, 0, height - side))
            windows.append({'selected_object_id': row.object_id, 'x_px': x, 'y_px': y, 'width_px': side,
                            'height_px': side, 'x_um': x * CT.PIXEL_UM, 'y_um': y * CT.PIXEL_UM,
                            'width_um': side * CT.PIXEL_UM, 'pyramid_level': 0})
        whole = np.zeros((*level.shape[-2:], 3), np.uint8)
        crops = [np.zeros((side, side, 3), np.uint8) for _ in windows]
        limits = {}
        for channel, rgb_index in [('DAPI', 2), ('BActin', 1), ('ECAD', 0)]:
            index = names.index(channel)
            low = level.pages[index].asarray()
            limits[channel] = np.percentile(np.log1p(low.astype(np.float32)), [1, 99.5]).tolist()
            whole[:, :, rgb_index] = scale(low, limits[channel])
            del low
            native = series.levels[0].pages[index].asarray()
            for crop, win in zip(crops, windows):
                x, y = win['x_px'], win['y_px']
                crop[:, :, rgb_index] = scale(native[y:y + side, x:x + side], limits[channel])
            del native
        meta = {'source_ome': RF.file_stat(image_path),
                'ome_metadata_sha256': hashlib.sha256(tf.ome_metadata.encode()).hexdigest(),
                'native_shape_yx': [height, width], 'whole_pyramid_level': whole_level,
                'whole_downsample': downsample, 'tile_pyramid_level': 0, 'pixel_um': CT.PIXEL_UM,
                'channels_rgb': {'red': 'ECAD', 'green': 'BActin', 'blue': 'DAPI'},
                'display': f'raw pixels, log1p, per-image per-channel level-{whole_level} percentiles 1..99.5',
                'display_limits_log1p': limits, 'despeckle': False,
                'detection_input_channels_for_context': ['DAPI', 'CD68', 'CD44', 'VIM', 'ECAD', 'LGALS3', 'panCK', 'CD3e', 'BActin']}
    return whole, crops, windows, meta


def render(tag: str, run: Path, out: Path, baseline_run: Path | None = None, seed: int = 29,
           tile_um: float = 160., image_dir: Path = DEFAULT_IMAGES, whole_level: int = 3) -> dict:
    folder = out / tag
    if folder.exists() and any(folder.iterdir()):
        raise ValueError(f'Refusing to overwrite existing review: {folder}')
    snap = RF.ready_inputs(run, tag)
    if snap is None:
        raise ValueError(f'{tag}: full augmentation verification is required before rendering')
    state_path = run / 'augmentation' / tag / 'state.json'
    state = json.loads(state_path.read_text())
    baseline_run = baseline_run or Path(state['inputs']['source_run'])
    if baseline_run.resolve() != Path(state['inputs']['source_run']).resolve():
        raise ValueError('Baseline must match the augmentation source run')
    inputs = {'new_object_ids': run / 'augmentation' / tag / 'new_object_ids.json',
              'baseline_tissue': baseline_run / 'inputs/tissue' / f'{tag}.geojson',
              'final_tissue': run / 'inputs/tissue' / f'{tag}.geojson',
              'baseline_cells': baseline_run / 'cells' / f'{tag}.parquet',
              'cells': run / 'cells' / f'{tag}.parquet', 'cell_geojson': run / 'cell_geojson' / f'{tag}.geojson',
              'renderer_code': Path(__file__)}
    signatures = {name: RF.signature(path) for name, path in inputs.items()}
    declared = state['steps']['augment_detect']['outputs']
    matched = [sha for path, sha in declared.items() if Path(path).resolve() == inputs['new_object_ids'].resolve()]
    if matched != [signatures['new_object_ids']['sha256']]:
        raise ValueError('New-cell ID inventory does not match the augmentation checkpoint')
    listed = json.loads(inputs['new_object_ids'].read_text())
    if not isinstance(listed, list) or len(listed) != len(set(listed)):
        raise ValueError('Accepted new IDs must be a unique list')
    new_ids = set(listed)
    baseline_ids = pd.read_parquet(inputs['baseline_cells'], columns=['Object ID'])['Object ID']
    if baseline_ids.isna().any() or baseline_ids.duplicated().any() or set(baseline_ids) & new_ids:
        raise ValueError('Baseline IDs are invalid or overlap accepted new IDs')
    frame = pd.read_parquet(inputs['cells'], columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    if set(frame['Object ID']) != set(baseline_ids) | new_ids:
        raise ValueError('Current cell inventory is not exactly baseline plus new IDs')
    ids, cells, nuclei, invalid = cell_geometries(inputs['cell_geojson'])
    baseline, final = tissue_geometry(inputs['baseline_tissue']), tissue_geometry(inputs['final_tissue'])
    report, detail = check_geometry(frame, ids, cells, nuclei, new_ids, baseline, final, invalid)
    selected = choose_fields(detail, seed)
    image_path = image_dir / f'{tag}.ome.tiff'
    image_before = RF.file_stat(image_path)
    whole, crops, windows, image_meta = native_views(image_path, selected, tile_um, whole_level)
    folder.mkdir(parents=True, exist_ok=True)
    bounds = shapely.bounds(cells)
    view = overlay(whole, ids, cells, nuclei, new_ids, baseline, final, downsample=2 ** whole_level)
    for i, win in enumerate(windows, 1):
        x, y, side = (round(win[k] / 2 ** whole_level) for k in ('x_px', 'y_px', 'width_px'))
        cv2.rectangle(view, (x, y), (x + side, y + side), (0, 255, 255), 1)
        cv2.putText(view, str(i), (x + 3, y + 13), cv2.FONT_HERSHEY_SIMPLEX, .45, (0, 255, 255), 1, cv2.LINE_AA)
    labeled = with_legend([whole, view], f'{tag} | whole level {whole_level} | raw composite / geometry | {len(new_ids)} new cells')
    cv2.imwrite(str(folder / 'whole.png'), labeled[:, :, ::-1])
    for i, (raw, window) in enumerate(zip(crops, windows), 1):
        x, y, side = window['x_px'], window['y_px'], window['width_px']
        near = np.flatnonzero((bounds[:, 2] >= x) & (bounds[:, 0] <= x + side) &
                              (bounds[:, 3] >= y) & (bounds[:, 1] <= y + side))
        local_ids = [ids[k] for k in near]
        window['visible_new_ids'] = [oid for oid in local_ids if oid in new_ids]
        window['visible_old_ids'] = [oid for oid in local_ids if oid not in new_ids]
        window['tile'] = i
        dapi = np.repeat(raw[:, :, 2:3], 3, axis=2)
        over = overlay(raw, local_ids, cells[near], nuclei[near], new_ids, baseline, final, origin=(x, y),
                       selected_id=window['selected_object_id'])
        title = f'{tag} | tile {i} | level 0, native pixels | {window["width_um"]:.1f} um | center {window["selected_object_id"]}'
        labeled = with_legend([dapi, raw, over], title)
        cv2.imwrite(str(folder / f'tile{i:02d}.png'), labeled[:, :, ::-1])
    detail[detail['is_new']].to_csv(folder / 'new_cells_geometry.csv', index=False)
    # Old exceptions are diagnostic: a concave cell's centroid can lie outside its own preserved shape.
    detail[(~detail['is_new']) & ((~detail['centroid_inside_final']) |
           (detail['cell_area_outside_final_um2'] > AREA_EPS_UM2) |
           (detail['nucleus_area_outside_final_um2'] > AREA_EPS_UM2))].to_csv(folder / 'old_geometry_exceptions.csv', index=False)
    if RF.ready_inputs(run, tag) != snap or RF.file_stat(image_path) != image_before:
        raise ValueError('Inputs changed during rendering; this output directory is incomplete and unreviewable')
    if {name: RF.signature(path) for name, path in inputs.items()} != signatures:
        raise ValueError('Review inputs changed during rendering; output is incomplete and unreviewable')
    report.update(image=tag, baseline_run=str(baseline_run.resolve()), run=str(run.resolve()))
    RF.atomic_json(folder / 'check.json', report)
    metadata = {'image': tag, 'created_at_utc': datetime.now(timezone.utc).isoformat(),
                'visual_review_status': 'pending', 'inputs': signatures, 'augmentation_verified_inputs': snap,
                'image_display': image_meta, 'sampling': {'seed': seed, 'method': 'seeded first new cell, then farthest centroids',
                'requested_tiles': 6, 'actual_tiles': len(windows), 'native_tile_um': tile_um},
                'tiles': windows, 'legend_rgb': COLORS,
                'output_sha256': {p.name: RF.sha256(p) for p in sorted(folder.iterdir()) if p.is_file()}}
    RF.atomic_json(folder / 'metadata.json', metadata)
    print(json.dumps({'image': tag, 'new_cells': len(new_ids), 'tiles': len(windows),
                      'geometry_checks_passed': report['geometry_checks_passed'], 'review': str(folder),
                      'visual_review_status': 'pending'}), flush=True)
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True, type=Path)
    ap.add_argument('--baseline', type=Path, help='default: augmentation state source_run')
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--out', type=Path, default=DEFAULT_OUT)
    ap.add_argument('--image-dir', type=Path, default=DEFAULT_IMAGES)
    ap.add_argument('--seed', type=int, default=29)
    ap.add_argument('--tile-um', type=float, default=160.)
    ap.add_argument('--whole-level', type=int, default=3)
    ap.add_argument('--workers', type=int, choices=(1, 2), default=1)
    args = ap.parse_args()
    cohort = [r['image_name'] for r in CT.cohort_rows()]
    tags = cohort if args.images == ['all'] else args.images
    if set(tags) - set(cohort) or len(set(tags)) != len(tags):
        ap.error('images must be unique cohort names, or all')
    if args.tile_um < 60 or args.whole_level < 0:
        ap.error('tile size must be at least 60 um and whole level nonnegative')
    with ProcessPoolExecutor(args.workers) as pool:
        reports = list(pool.map(render, tags, repeat(args.run.resolve()), repeat(args.out.resolve()),
                                repeat(args.baseline), repeat(args.seed), repeat(args.tile_um),
                                repeat(args.image_dir), repeat(args.whole_level)))
    return 0 if all(r['geometry_checks_passed'] for r in reports) else 1


if __name__ == '__main__':
    raise SystemExit(main())
