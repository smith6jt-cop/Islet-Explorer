#!/usr/bin/env python3
"""Tile sheets for checking the final typing against the image (owner 2026-10-29: "then do what you suggested").

Samples cells (seed 0) from the final typing run's typing/identity, in five groups:
  compartment_tie  cells two compartments tied on depth, settled by brightness and area coverage (the six most
                   frequent pairs, 16 each): green = the placed compartment's canary-positive defining markers,
                   magenta = the runner-up's
  lineage_tie      the same between the lineages of one compartment (the four most frequent pairs, 16 each)
  hema_cd45neg     Hematopoietic CD45- (CD4 alone, CD8a alone, CD20 alone; 16 each): CD45, CD3e and the positive marker
  no_marker        status no_marker_above_background (24): the four markers closest to their levels
  shared_only      status shared_marker_only (24): the positive shared markers and the two defining markers closest
                   to their levels
Each row is one cell: a composite and grayscale panels of the deciding markers and DAPI, 48 um around the centroid,
despeckled as measured (redsea_cells.despeckle). Grayscale: black = 0, mid-gray = the image's level of that marker,
white = twice the level. The cell's outline is white, its neighbors' gray (composite only).

Writes <out>/sheets/<group>_<nn>.png (8 rows each) and <out>/manifest.csv (per row: the cell, its call, the markers
shown with their cell value, level, fold over level, area coverage (typing/coverage) and the coverage of the cell's
interior, 2 px in from its outline).

    python scripts/research/typing_accuracy_tiles_20261029.py --workers 5
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import identity as ID  # noqa: E402
import redsea_cells as RS  # noqa: E402
import registry as RG  # noqa: E402
import area_coverage as AC  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'typing_accuracy_20261029'
CROP_UM = 48.0
ZOOM = 2
ROWS = 8
PANELS = 6                       # composite + up to four markers + DAPI
SEED = 0


def polygon_components(geometry: dict) -> list[list[np.ndarray]]:
    """Retain every exterior and hole for diagnostic outlines and own-cell masks."""
    coordinates = geometry.get('coordinates') or []
    if not coordinates:
        return []
    if geometry['type'] not in ('Polygon', 'MultiPolygon'):
        raise ValueError(f"Unsupported cell geometry: {geometry['type']}")
    polygons = coordinates if geometry['type'] == 'MultiPolygon' else [coordinates]
    result = []
    for polygon in polygons:
        rings = [np.asarray(ring, dtype=float) for ring in polygon]
        if any(ring.ndim != 2 or ring.shape[1] != 2 or len(ring) < 3 or
               not np.isfinite(ring).all() for ring in rings):
            raise ValueError('Invalid cell polygon ring in diagnostic renderer')
        if rings:
            result.append(rings)
    return result


def polygon_bounds(components: list[list[np.ndarray]]) -> tuple[float, float, float, float]:
    exteriors = np.concatenate([rings[0] for rings in components])
    return (*exteriors.min(axis=0), *exteriors.max(axis=0))


def draw_cell_outline(image: np.ndarray, components: list[list[np.ndarray]], origin, color) -> None:
    for rings in components:
        for ring in rings:
            cv2.polylines(image, [np.round(ring - origin).astype(np.int32)], True, color, 1)


def cell_masks(components: list[list[np.ndarray]], shape: tuple[int, int], origin) -> tuple[np.ndarray, np.ndarray]:
    """Rasterize the cell union, keeping holes local, then erode the whole cell by two pixels."""
    own = np.zeros(shape, np.uint8)
    for rings in components:
        part = np.zeros(shape, np.uint8)
        local = [np.round(ring - origin).astype(np.int32) for ring in rings]
        cv2.fillPoly(part, [local[0]], 1)
        if len(local) > 1:
            cv2.fillPoly(part, local[1:], 0)
        own |= part
    interior = cv2.erode(own, np.ones((3, 3), np.uint8), iterations=2).astype(bool)
    return own, interior


def plan_and_markers(run: Path = RUN):
    reg = RG.load()
    markers = ID.used_markers(run, reg)
    return reg, markers, ID.build_plan(reg, markers)


def lineage_runner(row: pd.Series, info: dict, lineage: str) -> str:
    """The other lineage of the compartment with a canary-positive specific marker (else a positive one)."""
    for col in ('canary__', 'pos__'):
        others = [lin for lin, L in info['lineages'].items() if lin != lineage
                  and any(row.get(col + m, 0) for m in L['specific'])]
        if others:
            return others[0]
    return ''


def compartment_runner(row: pd.Series, plan: dict) -> str:
    """The runner-up compartment; a tie settled at plain positivity has none recorded (identity writes the canary
    level's), so it is the other compartment with a positive defining marker."""
    if row['runner_up_compartment']:
        return row['runner_up_compartment']
    for col in ('canary__', 'pos__'):
        others = [c for c, info in plan['comp'].items() if c != row['compartment']
                  and any(row.get(col + m, 0) for m in info['defining'])]
        if others:
            return others[0]
    return ''


def candidates(reg, plan, run: Path = RUN) -> pd.DataFrame:
    rows = []
    for r in CT.cohort_rows():
        tag = r['image_name']
        d = pd.read_parquet(run / 'typing' / 'identity' / f'{tag}.parquet')
        d['image'] = tag
        ct = d['tags'].str.contains('compartment_tie_brightness_coverage', regex=False)
        lt = d['tags'].str.contains('lineage_tie_brightness_coverage', regex=False)
        a = d[ct].copy()
        a['runner_up_compartment'] = a.apply(lambda x: compartment_runner(x, plan), axis=1) if len(a) else ''
        a = a.assign(group='compartment_tie', pair=a['compartment'] + ' over ' + a['runner_up_compartment'])
        b = d[lt].copy()
        b['group'] = 'lineage_tie'
        b['pair'] = b.apply(lambda x: f"{x['lineage']} over "
                                      f"{lineage_runner(x, plan['comp'][x['compartment']], x['lineage'])}", axis=1) \
            if len(b) else pd.Series(dtype=object)
        h = d[d['label'] == 'Hematopoietic CD45-'].copy()
        hm = ['CD3e', 'CD8a', 'CD4', 'CD20', 'CD68', 'CD163']
        h['pair'] = h[[f'pos__{m}' for m in hm]].apply(lambda s: '+'.join(m for m, v in zip(hm, s) if v) or '-', axis=1)
        h = h[h['pair'].isin(['CD4', 'CD8a', 'CD20'])].assign(group='hema_cd45neg')
        n = d[d['status'] == 'no_marker_above_background'].assign(group='no_marker', pair='')
        s = d[d['status'] == 'shared_marker_only'].assign(group='shared_only', pair='')
        rows.append(pd.concat([a, b, h, n, s]))
    return pd.concat(rows, ignore_index=True)


def exclude_reviewed(c: pd.DataFrame, manifests: list[Path]) -> pd.DataFrame:
    excluded = set()
    for path in manifests:
        previous = pd.read_csv(path, dtype=str, keep_default_na=False)
        if not {'image', 'object_id'} <= set(previous):
            raise ValueError(f'{path}: exclusion manifest needs image and object_id columns')
        excluded.update(zip(previous['image'], previous['object_id']))
    return c.loc[[key not in excluded for key in zip(c['image'], c['object_id'])]].copy()


def sample(c: pd.DataFrame, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = []
    seen = set()
    for group, n_pairs, n in (('compartment_tie', 6, 16), ('lineage_tie', 4, 16), ('hema_cd45neg', 3, 16),
                              ('no_marker', 1, 24), ('shared_only', 1, 24)):
        g = c[c['group'] == group]
        for pair in g['pair'].value_counts().index[:n_pairs]:
            p = g[g['pair'] == pair]
            p = p.loc[[key not in seen for key in zip(p['image'], p['object_id'])]]
            chosen = p.iloc[rng.choice(len(p), size=min(n, len(p)), replace=False)]
            seen.update(zip(chosen['image'], chosen['object_id']))
            out.append(chosen)
    return pd.concat(out, ignore_index=True) if out else c.iloc[:0].copy()


def shown_markers(row: pd.Series, plan: dict, markers: list[str], reg) -> tuple[list[str], list[str]]:
    """(green markers, magenta markers): at most two each, the highest fold over level first."""
    def top(ms, k=2, col='canary__'):
        ms = [m for m in ms if row.get(col + m, 0)] or [m for m in ms if row.get('pos__' + m, 0)]
        return sorted(ms, key=lambda m: -row.get('fold__' + m, -9))[:k]
    if row['group'] == 'compartment_tie':
        return top(plan['comp'][row['compartment']]['defining']), top(plan['comp'][row['runner_up_compartment']]['defining'])
    if row['group'] == 'lineage_tie':
        info = plan['comp'][row['compartment']]
        win, lose = row['pair'].split(' over ')
        return top(info['lineages'][win]['specific']), top(info['lineages'].get(lose, {'specific': []})['specific'])
    if row['group'] == 'hema_cd45neg':
        return [row['pair']], ['CD45', 'CD3e'] if row['pair'] != 'CD3e' else ['CD45']
    if row['group'] == 'no_marker':
        return sorted(markers, key=lambda m: -row.get('fold__' + m, -9))[:4], []
    shared = [m for m in markers if row.get('pos__' + m, 0)][:2]
    defining = [m for m in markers if m not in shared and reg.role(m) == 'defining']
    return shared, sorted(defining, key=lambda m: -row.get('fold__' + m, -9))[:2]


def render_image(tag: str, rows: pd.DataFrame, run: Path = RUN) -> list[dict]:
    reg, markers, plan = plan_and_markers(run)
    entry = CT.manifest(run)[tag]
    names = CT.channels(run, entry['entry_id'])
    img = RS.Channels(entry['image_path'])
    H, W = img.shape[1:]
    lv = pd.read_csv(run / 'typing' / 'levels.csv')
    lv = lv[lv['image'] == tag].set_index('marker')
    cov = AC.read_validated(run, tag, rows['object_id'], markers).set_index('object_id')
    rs = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet').set_index('object_id')
    xy = pd.read_parquet(run / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    xy = xy.set_index('Object ID')
    feats = json.loads((run / 'cell_geojson' / f'{tag}.geojson').read_text())
    feats = feats['features'] if isinstance(feats, dict) else feats
    half = int(round(CROP_UM / CT.PIXEL_UM / 2))
    # polygons of every cell, with bounding boxes, to draw the sampled cells and their neighbors
    polys, boxes, ids = [], [], []
    for f in feats:
        components = polygon_components(f.get('geometry') or {})
        if not components:
            continue
        polys.append(components)
        boxes.append(polygon_bounds(components))
        ids.append(f['id'])
    boxes = np.asarray(boxes)
    by_id = {o: i for i, o in enumerate(ids)}
    del feats
    rows = rows.copy()
    shown = [shown_markers(r, plan, markers, reg) for _, r in rows.iterrows()]
    shown = [(mg, g) if r.get('blind_color_swap', False) else (g, mg)
             for (_, r), (g, mg) in zip(rows.iterrows(), shown)]
    rows['green'] = [g for g, _ in shown]
    rows['magenta'] = [mg for _, mg in shown]
    need = sorted({m for gs in rows['green'] for m in gs} | {m for ms in rows['magenta'] for m in ms} | {'DAPI'})
    windows = []
    for _, r in rows.iterrows():
        cx = xy.loc[r['object_id'], 'Centroid X µm'] / CT.PIXEL_UM
        cy = xy.loc[r['object_id'], 'Centroid Y µm'] / CT.PIXEL_UM
        x0, y0 = int(max(0, min(W - 2 * half, cx - half))), int(max(0, min(H - 2 * half, cy - half)))
        windows.append((x0, y0))
    crops = {}                                   # (row number, marker) -> crop; one channel in memory at a time
    for m in need:
        ch = reg.table.loc[m, 'channel'] if m in reg.table.index else m
        ch = ch if ch in names else m
        full = RS.despeckle(img.native(names.index(ch)))
        for k, (x0, y0) in enumerate(windows):
            crops[k, m] = full[y0:y0 + 2 * half, x0:x0 + 2 * half].astype(np.float32)
        del full
    img.close()
    out = []
    for k, (_, r) in enumerate(rows.iterrows()):
        oid = r['object_id']
        x0, y0 = windows[k]
        near = np.flatnonzero((boxes[:, 2] >= x0) & (boxes[:, 0] <= x0 + 2 * half) &
                              (boxes[:, 3] >= y0) & (boxes[:, 1] <= y0 + 2 * half))
        _, interior = cell_masks(polys[by_id[oid]] if oid in by_id else [],
                                 (2 * half, 2 * half), (x0, y0))
        panels, info = [], {}
        for m in list(r['green']) + list(r['magenta']):
            a = crops[k, m]
            level = float(lv.loc[m, 'level_raw']) if m in lv.index else float(np.percentile(a, 99))
            panels.append((m, np.clip(a / (2 * level), 0, 1)))
            meas = lv.loc[m, 'measurement'] if m in lv.index else 'whole_cell'
            col = reg.column(m, meas) if m in reg.table.index else ''
            info[m] = {'value': float(rs.loc[oid, col]) if col in rs.columns else np.nan, 'level': level,
                       'fold': float(r.get('fold__' + m, np.nan)), 'measurement': meas,
                       'coverage': float(cov.loc[oid, 'cov__' + m]) if 'cov__' + m in cov.columns else np.nan,
                       'coverage_interior': float(cov.loc[oid, AC.INTERIOR_PREFIX + m]),
                       'coverage_interior_crop': float((a[interior] >= level).mean()) if interior.any() else np.nan,
                       'n_interior_pixels': int(cov.loc[oid, AC.INTERIOR_COUNT])}
        dapi = crops[k, 'DAPI']
        dapi = np.clip(dapi / max(np.percentile(dapi, 99.5), 1), 0, 1)
        g = np.max([p for m, p in panels if m in r['green']], axis=0) if r['green'] else np.zeros_like(dapi)
        mg = np.max([p for m, p in panels if m in r['magenta']], axis=0) if r['magenta'] else np.zeros_like(dapi)
        comp = np.dstack([mg, np.maximum(g, 0), np.maximum(mg, 0.6 * dapi)])
        comp = (np.clip(comp, 0, 1) * 255).astype(np.uint8)
        for i in near:
            if ids[i] != oid:
                draw_cell_outline(comp, polys[i], (x0, y0), (90, 90, 90))
        tiles = [comp] + [np.dstack([(p * 255).astype(np.uint8)] * 3) for _, p in panels] + \
                [np.dstack([(dapi * 255).astype(np.uint8)] * 3)]
        if oid in by_id:
            for t in tiles:
                draw_cell_outline(t, polys[by_id[oid]], (x0, y0), (255, 255, 255))
        out.append({**{k: r[k] for k in ('group', 'pair', 'image', 'object_id', 'label', 'status', 'tags', 'compartment',
                                        'runner_up_compartment', 'lineage')},
                    'review_id': r.get('review_id', ''), 'blind_color_swap': bool(r.get('blind_color_swap', False)),
                    'green': ';'.join(r['green']), 'magenta': ';'.join(r['magenta']),
                    'panel_names': ['composite'] + [m for m, _ in panels] + ['DAPI'],
                    'tiles': [cv2.resize(t, None, fx=ZOOM, fy=ZOOM, interpolation=cv2.INTER_NEAREST) for t in tiles],
                    'markers': info})
    return out


def sheets(done: list[dict], out: Path = OUT, blind: bool = False) -> pd.DataFrame:
    sheet_dir = out / ('blinded' if blind else 'sheets')
    sheet_dir.mkdir(parents=True, exist_ok=True)
    size = done[0]['tiles'][0].shape[0]
    head = 30
    manifest = []
    for group in (['review'] if blind else dict.fromkeys(d['group'] for d in done)):
        items = list(done) if blind else [d for d in done if d['group'] == group]
        items.sort(key=(lambda d: d['review_id']) if blind else (lambda d: (d['pair'], d['image'], d['object_id'])))
        for s in range(0, len(items), ROWS):
            chunk = items[s:s + ROWS]
            sheet = Image.new('RGB', (PANELS * size, len(chunk) * (size + head)), (0, 0, 0))
            draw = ImageDraw.Draw(sheet)
            name = f'{group}_{s // ROWS + 1:02d}'
            for k, d in enumerate(chunk):
                y = k * (size + head)
                nums = '  '.join(f"{m}: x{np.expm1(v['fold']) + 1 if np.isfinite(v['fold']) else float('nan'):.2f} "
                                 f"cov {v['coverage']:.2f}/{v['coverage_interior']:.2f}" for m, v in d['markers'].items())
                if blind:
                    draw.text((4, y + 2), f"row {k + 1}  {d['review_id']}  green {d['green'] or '-'} | "
                                          f"magenta {d['magenta'] or '-'}", fill=(255, 255, 255))
                else:
                    draw.text((4, y + 2), f"row {k + 1}  {d['image']} {d['object_id'][:8]}  call: {d['label']}  "
                                      f"[{d['pair']}]  green {d['green'] or '-'} | magenta {d['magenta'] or '-'}",
                          fill=(255, 255, 255))
                    draw.text((4, y + 15), nums, fill=(200, 200, 200))
                for j, (t, pn) in enumerate(zip(d['tiles'], d['panel_names'])):
                    col = j if pn in ('composite',) else (j if j < len(d['tiles']) - 1 else PANELS - 1)
                    sheet.paste(Image.fromarray(t), (col * size, y + head))
                    draw.text((col * size + 4, y + head + 4), pn, fill=(255, 255, 0))
                manifest.append({'sheet': name, 'row': k + 1, 'review_id': d['review_id'],
                                 'blind_color_swap': d['blind_color_swap'], **{x: d[x] for x in (
                    'group', 'pair', 'image', 'object_id', 'label', 'status', 'tags', 'compartment',
                    'runner_up_compartment', 'lineage', 'green', 'magenta')},
                    'markers': json.dumps(d['markers'])})
            sheet.save(sheet_dir / f'{name}.png')
    m = pd.DataFrame(manifest)
    key_dir = out / 'key' if blind else out
    key_dir.mkdir(parents=True, exist_ok=True)
    m.to_csv(key_dir / 'manifest.csv', index=False)
    if blind:
        reviewer = m[['sheet', 'row', 'review_id', 'green', 'magenta']].copy()
        reviewer['verdict'] = ''
        reviewer['notes'] = ''
        reviewer.to_csv(sheet_dir / 'manifest.csv', index=False)
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=5)
    ap.add_argument('--run', type=Path, default=RUN)
    ap.add_argument('--out', type=Path, default=OUT)
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--exclude-manifest', type=Path, nargs='+', default=[],
                    help='exclude every (image, object_id) in previous full review manifests')
    ap.add_argument('--blind', action='store_true', help='hide calls and strata, randomize color assignment/order; '
                    'give reviewers only the blinded/ folder, keep key/ separate')
    a = ap.parse_args()
    if a.out.exists() and any(a.out.iterdir()):
        ap.error(f'{a.out} is not empty; choose a fresh output directory to preserve previous reviews')
    reg, markers, plan = plan_and_markers(a.run)
    available = candidates(reg, plan, a.run)
    fresh = exclude_reviewed(available, a.exclude_manifest)
    picked = sample(fresh, a.seed)
    if picked.empty:
        ap.error('no eligible unreviewed cells remain')
    rng = np.random.default_rng(a.seed)
    picked = picked.iloc[rng.permutation(len(picked))].reset_index(drop=True)
    picked['review_id'] = [f'R{i + 1:04d}' for i in range(len(picked))]
    picked['blind_color_swap'] = rng.integers(0, 2, len(picked)).astype(bool) if a.blind else False
    print(picked.groupby(['group', 'pair']).size().to_string(), flush=True)
    done = []
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(render_image, tag, g, a.run) for tag, g in picked.groupby('image')]
        for f in futs:
            done += f.result()
    m = sheets(done, a.out, a.blind)
    (a.out / 'sampling.json').write_text(json.dumps({'run': str(a.run.resolve()), 'seed': a.seed, 'blind': a.blind,
        'candidate_rows_before_exclusion': len(available), 'candidate_rows_after_exclusion': len(fresh),
        'sampled_unique_cells': len(m), 'exclusion_manifests': [
            {'path': str(p.resolve()), 'sha256': AC.sha256(p)} for p in a.exclude_manifest]}, indent=2))
    print(f'{len(m)} cells on {m["sheet"].nunique()} sheets in {a.out}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
