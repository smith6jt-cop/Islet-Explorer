#!/usr/bin/env python3
"""redsea_cells.py -- Senior's geometry QC and compartment-resolved REDSEA, on this panel's images.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/phenocycler/redsea.py` (process_donor, rasterize_mask,
contact_matrix[_gapbridge], compensate_compartments, apply_spillover_context_eligibility) and
`phenocycler/geometry_qc.py` (evaluate_geometry_metrics). The operator and parameters are Senior's
(codex_typing.REDSEA): corrected = clip(data + comp_mode*edge - alpha*(F @ edge), 0), F the donor-normalised
contact matrix, resolved onto Nucleus / Cytoplasm / Membrane / Cell. Channels in codex_typing.NO_SPILLOVER pass
through uncorrected.

CODEX differences: cells and nuclei come from `$OUT/cell_geojson/<image>.geojson` (QuPath's export of the
InstanSeg cells, pixel coordinates, one feature per cell with its nucleusGeometry); the QuPath measurements
`$OUT/cells/<image>.parquet` supply area, solidity and centroids; channels are read from the source image the
QuPath entry names, in the entry's channel order. A cell without a nucleus is excluded from analysis and
from spillover context, as Senior's pipeline config sets it. Of Senior's other QC flags (area under 20 um2 or
over 4x the median, solidity under 0.70, raster/reported area ratio) none excludes a cell here (owner,
2026-10-16); they are written to `qc_flags`.

Despeckle (owner 2026-10-28: "some images look like they need a despeckle filter prior to cell intensity
measurements"): every channel of every image goes through a 3x3 median (ImageJ's Process > Noise > Despeckle) before
the cells are measured. Agent determination on the images (output/typing_final_checks_20261028/despeckle): a
fluorescence peak spreads over several pixels, so its brightest neighbor stands at >= 0.5 of its height (ND_0112:
0-2 % of bright peaks below 0.25), while ND_6232 and Aab_6558 carry single-pixel specks (20-28 % of bright peaks below
0.25), which a 3x3 median removes and smooth signal passes. One filter for every image and channel (owner 2026-10-28:
fix the method, never flag images).

The empty-slide null of every channel is measured from the same despeckled channels (levels.py floors the edge and
gate levels at it): means of 19 x 19-pixel windows (7.2 um, about a cell's measurement area) wholly in the slide area at
least 20 um outside the image's Tissue outline, 20,000 seeded windows, their p50 / p99 / p99.9 (the method of
scripts/research/empty_area_background_20261028.py, which read the raw pixels); main() writes them to
typing/empty_area_nulls.csv.

Before writing, the uncorrected whole-cell mean of every channel (before the despeckle) is compared with QuPath's own
`Cell: <channel>: Mean` over the cells of at least 20 um2 whose raster and polygon areas agree; a channel whose correlation falls below ALIGN_MIN_R stops the image (Senior's
alignment_check: a key or channel mismatch invalidates every corrected value).

    python scripts/uniform/typing/redsea_cells.py --out data/development/uniform-20donor-20261015 --images all
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import cv2
import tifffile
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from shapely.geometry import shape
from skimage.draw import polygon as draw_polygon
from skimage.morphology import binary_dilation, diamond
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402

COMPARTMENTS = ('Nucleus', 'Cytoplasm', 'Membrane', 'Cell')
ALIGN_MIN_R = 0.95
DESPECKLE = {'filter': 'median', 'size': 3, 'implementation': 'cv2.medianBlur', 'channels': 'all'}


EMPTY = {'window_px': 19, 'windows': 20000, 'margin_um': 20.0, 'seed': 0}


def empty_windows(out: Path, tag: str, H: int, W: int) -> tuple[np.ndarray, np.ndarray]:
    """Centers of EMPTY['window_px']-square windows wholly in the slide area at least EMPTY['margin_um'] outside the
    image's Tissue outline (`inputs/tissue/<image>.geojson`, pixels), a seeded sample of EMPTY['windows']."""
    raw = json.loads((out / 'inputs' / 'tissue' / f'{tag}.geojson').read_text())
    geom = shape((raw['features'] if isinstance(raw, dict) else raw)[0]['geometry'])
    mask = np.ones((H, W), np.uint8)
    for poly in getattr(geom, 'geoms', [geom]):
        cv2.fillPoly(mask, [np.round(np.asarray(poly.exterior.coords)).astype(np.int32)], 0)
        for hole in poly.interiors:
            cv2.fillPoly(mask, [np.round(np.asarray(hole.coords)).astype(np.int32)], 1)
    r = int(round(EMPTY['margin_um'] / CT.PIXEL_UM))
    mask = cv2.erode(mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1)))
    w = EMPTY['window_px']
    ok = cv2.erode(mask, np.ones((w, w), np.uint8)).astype(bool)
    e = w // 2 + 1
    ok[:e], ok[-e:], ok[:, :e], ok[:, -e:] = False, False, False, False
    ys, xs = np.nonzero(ok)
    if len(ys):
        pick = np.random.default_rng(EMPTY['seed']).choice(len(ys), min(EMPTY['windows'], len(ys)), replace=False)
        ys, xs = ys[pick], xs[pick]
    return ys, xs


def window_means(a: np.ndarray, ys: np.ndarray, xs: np.ndarray) -> np.ndarray:
    ii = cv2.integral(np.ascontiguousarray(a, dtype=a.dtype.newbyteorder('=')), sdepth=cv2.CV_64F)
    h = EMPTY['window_px'] // 2
    y0, y1, x0, x1 = ys - h, ys + h + 1, xs - h, xs + h + 1
    return (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]) / EMPTY['window_px'] ** 2


def despeckle(a: np.ndarray) -> np.ndarray:
    """ImageJ's Despeckle: the 3x3 median of each pixel (borders replicated), on the channel's own integer values."""
    native = a.dtype.newbyteorder('=')        # OpenCV reads native byte order; some of these TIFFs are big-endian
    return cv2.medianBlur(np.ascontiguousarray(a, dtype=native), DESPECKLE['size'])


def log(tag, msg):
    print(f'[{tag}] {msg}', flush=True)


# --------------------------------------------------------------------------- raster (redsea.py, unchanged)
def draw_geometry(target, geom, label, shape):
    polys = geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]
    for poly in polys:
        ext = np.asarray(poly[0], dtype=np.float64)
        if ext.shape[0] < 3:
            continue
        pr, pc = draw_polygon(ext[:, 1], ext[:, 0], shape=shape)
        target[pr, pc] = label
        for hole in poly[1:]:
            h = np.asarray(hole, dtype=np.float64)
            if h.shape[0] < 3:
                continue
            hr, hc = draw_polygon(h[:, 1], h[:, 0], shape=shape)
            target[hr, hc] = 0


def rasterize(features, shape):
    mask = np.zeros(shape, dtype=np.int32)
    nuc = np.zeros(shape, dtype=np.int32)
    oids, has_nucleus = [], []
    for i, ft in enumerate(features):
        oids.append(ft['id'])
        draw_geometry(mask, ft['geometry'], i + 1, shape)
        ng = ft.get('nucleusGeometry')
        present = bool(ng and ng.get('coordinates'))
        has_nucleus.append(present)
        if present:
            draw_geometry(nuc, ng, i + 1, shape)
    np.putmask(nuc, nuc != mask, 0)       # nucleus a strict subset of its own cell's final pixels
    return mask, nuc, oids, np.asarray(has_nucleus)


def contact_matrix(mask, n):
    H, W = mask.shape
    keys = []
    for dy, dx in [(0, 1), (1, 0), (1, 1), (1, -1)]:
        if dx >= 0:
            a, b = mask[:H - dy, :W - dx], mask[dy:, dx:]
        else:
            a, b = mask[:H - dy, -dx:], mask[dy:, :W + dx]
        v = (a > 0) & (b > 0) & (a != b)
        if not v.any():
            continue
        aa, bb = a[v].astype(np.int64), b[v].astype(np.int64)
        keys.append(np.minimum(aa, bb) * (np.int64(n) + 1) + np.maximum(aa, bb))
    if not keys:
        return sp.csr_matrix((n, n), dtype=np.float64)
    uk, cnt = np.unique(np.concatenate(keys), return_counts=True)
    lo, hi = (uk // (n + 1)).astype(np.int64), (uk % (n + 1)).astype(np.int64)
    M = sp.coo_matrix((cnt.astype(np.float64), (lo - 1, hi - 1)), shape=(n, n)).tocsr()
    return M + M.T


def contact_matrix_gapbridge(mask, n, grow):
    if grow <= 0:
        return contact_matrix(mask, n)
    bg = mask == 0
    dist, idx = distance_transform_edt(bg, return_distances=True, return_indices=True)
    grown = mask[tuple(idx)]
    del idx
    grown[bg & (dist > grow + 0.5)] = 0
    del dist, bg
    return contact_matrix(grown, n)


def compensate_compartments(data, edge, nuc, nuc_edge, sizes, bandcount, nuccount, nucbandcount, M,
                            comp_mode, alpha):
    rowsum = np.asarray(M.sum(1)).ravel()
    inv = np.zeros_like(rowsum)
    nz = rowsum > 0
    inv[nz] = 1.0 / rowsum[nz]
    F = M @ sp.diags(inv)                                  # norm_form "donor"
    sub = np.asarray(alpha, dtype=np.float64) * (F @ edge)
    cyto, cyto_edge = data - nuc, edge - nuc_edge
    cytocount, cytobandcount = sizes - nuccount, bandcount - nucbandcount
    safe_band = np.where(bandcount > 0, bandcount, 1.0)[:, None]
    nuc_share = np.where(bandcount[:, None] > 0, nucbandcount[:, None] / safe_band, 0.0)
    cyto_share = np.where(bandcount[:, None] > 0, cytobandcount[:, None] / safe_band, 0.0)
    sums = {'Nucleus': nuc + comp_mode * nuc_edge - sub * nuc_share,
            'Cytoplasm': cyto + comp_mode * cyto_edge - sub * cyto_share,
            'Membrane': edge + comp_mode * edge - sub,
            'Cell': data + comp_mode * edge - sub}
    counts = {'Nucleus': nuccount, 'Cytoplasm': cytocount, 'Membrane': bandcount, 'Cell': sizes}
    means = {}
    for name in COMPARTMENTS:
        s = np.clip(sums[name], 0, None)
        c = counts[name]
        with np.errstate(invalid='ignore', divide='ignore'):
            m = s / np.where(c > 0, c, np.nan)[:, None]
        m[c <= 0, :] = np.nan
        means[name] = m
    return means, counts, int((~nz).sum())


# --------------------------------------------------------------------------- geometry QC (geometry_qc.py)
def geometry_qc(oids, has_nucleus, cells, sizes_px, nuc_px, cfg=CT.GEOMETRY_QC):
    """Senior's analysis / estimation / spillover-context eligibility, one row per geometry feature."""
    q = cells.set_index('Object ID').reindex(oids)
    x, y = q['Centroid X µm'].to_numpy(float), q['Centroid Y µm'].to_numpy(float)
    area = q['Cell: Area µm^2'].to_numpy(float)
    solidity = q['Cell: Solidity'].to_numpy(float)
    px_area = CT.PIXEL_UM ** 2
    present = np.ones(len(oids), bool)
    missing_nucleus = ~has_nucleus
    valid_centroid = np.isfinite(x) & np.isfinite(y)
    valid_area = np.isfinite(area) & (area > 0)
    valid_raster = sizes_px > 0
    valid_solidity = np.isfinite(solidity) & (solidity >= 0) & (solidity <= 1)
    valid_nuc_raster = nuc_px > 0
    with np.errstate(invalid='ignore', divide='ignore'):
        raster_ratio = sizes_px * px_area / area
        nc_ratio = nuc_px / sizes_px
    valid_scale = np.isfinite(raster_ratio) & (raster_ratio >= cfg['min_raster_area_ratio']) \
        & (raster_ratio <= cfg['max_raster_area_ratio'])

    # duplicates: centroids within 1 um and areas within 15 %; survivor = QC-viable, most raster pixels, id
    cand = np.flatnonzero(present & valid_area & valid_raster & valid_centroid)
    losers = set()
    n_groups = 0
    if len(cand) > 1:
        pairs = cKDTree(np.column_stack([x[cand], y[cand]])).query_pairs(cfg['duplicate_radius_um'],
                                                                           output_type='ndarray')
        parent = {}

        def find(i):
            while parent.get(i, i) != i:
                parent[i] = parent.get(parent[i], parent[i])
                i = parent[i]
            return i
        for a, b in pairs:
            i, j = int(cand[a]), int(cand[b])
            if abs(area[i] - area[j]) / max(area[i], area[j], 1e-12) <= cfg['duplicate_area_tolerance']:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[max(ri, rj)] = min(ri, rj)
        groups = {}
        for i in {int(v) for p in pairs for v in (cand[p[0]], cand[p[1]])}:
            groups.setdefault(find(i), []).append(i)
        viable = (present & ~missing_nucleus & valid_nuc_raster & valid_centroid & valid_area
                  & (area >= cfg['min_cell_area_um2']) & valid_raster & valid_scale & valid_solidity
                  & (solidity >= cfg['min_cell_solidity']))
        for members in groups.values():
            if len(members) < 2:
                continue
            n_groups += 1
            keep = min(members, key=lambda i: (0 if viable[i] else 1, -sizes_px[i], str(oids[i])))
            losers |= {i for i in members if i != keep}
    duplicate = np.zeros(len(oids), bool)
    duplicate[list(losers)] = True

    provisional = present & ~missing_nucleus & valid_area & (area >= cfg['min_cell_area_um2']) \
        & valid_solidity & (solidity >= cfg['min_cell_solidity']) & ~duplicate
    max_area = cfg['max_cell_area_median_multiple'] * float(np.median(area[provisional]))
    flags = {
        'missing_nucleus_geometry': missing_nucleus,
        'invalid_centroid': ~valid_centroid,
        'invalid_cell_area': ~valid_area,
        'cell_area_below_min': valid_area & (area < cfg['min_cell_area_um2']),
        'cell_area_above_max': valid_area & (area > max_area),
        'invalid_solidity': ~valid_solidity,
        'low_solidity': valid_solidity & (solidity < cfg['min_cell_solidity']),
        'raster_empty': ~valid_raster,
        'raster_area_mismatch': valid_area & valid_raster & ~valid_scale,
        'nucleus_raster_empty': has_nucleus & ~valid_nuc_raster,
        'duplicate_detection': duplicate,
    }
    reasons = np.full(len(oids), '', dtype=object)      # what excludes the cell (codex_typing.QC_EXCLUDES)
    noted = np.full(len(oids), '', dtype=object)        # every flag, for the record
    for name, f in flags.items():
        noted[f] = np.where(noted[f] == '', name, noted[f] + ';' + name)
        if name in CT.QC_EXCLUDES:
            reasons[f] = np.where(reasons[f] == '', name, reasons[f] + ';' + name)
    analysis = reasons == ''
    no_cyto = analysis & np.isfinite(nc_ratio) & (nc_ratio >= cfg['no_cytoplasm_nc_ratio'])
    spill = ~missing_nucleus & ~(has_nucleus & ~valid_nuc_raster) & valid_raster & ~duplicate
    frame = pd.DataFrame({'object_id': oids, 'qc_analysis_eligible': analysis,
                          'qc_estimation_eligible': analysis & ~no_cyto,
                          'qc_spillover_context_eligible': spill, 'qc_reasons': reasons, 'qc_flags': noted,
                          'cell_area_px': sizes_px, 'nucleus_area_px': nuc_px})
    summary = {'n': len(oids), 'analysis_eligible': int(analysis.sum()),
               'estimation_eligible': int((analysis & ~no_cyto).sum()),
               'spillover_context_eligible': int(spill.sum()), 'duplicate_groups': n_groups,
               'cell_area_max_um2': max_area,
               'excluding': list(CT.QC_EXCLUDES),
               'reasons': {k: int(v.sum()) for k, v in flags.items()}}
    return frame, summary


# --------------------------------------------------------------------------- image
class Channels:
    """One channel at a time from the source image, in the QuPath entry's channel order."""

    def __init__(self, path):
        self.tf = tifffile.TiffFile(path)
        level = self.tf.series[0].levels[0]
        self.shape = level.shape
        self.memmap = level.dataoffset is not None and len(level.pages) == 1
        self.arr = tifffile.memmap(path, series=0, level=0, mode='r') if self.memmap else None
        self.pages = level.pages

    def __call__(self, k):
        return self.native(k).astype(np.float32)

    def native(self, k):
        """Channel k as stored (the image's own integer type)."""
        if self.memmap:
            return np.asarray(self.arr[k])
        return self.pages[k].asarray()

    def close(self):
        self.arr = None
        self.tf.close()


def process(out: Path, tag: str, force: bool = False, dest_dir: Path | None = None) -> dict:
    dest = (dest_dir or CT.typing_dir(out) / 'redsea') / f'{tag}.parquet'
    rec_path = dest.with_suffix('.json')
    if dest.exists() and rec_path.exists() and not force:
        log(tag, 'done already')
        return json.loads(rec_path.read_text())
    t_all = time.time()
    entry = CT.manifest(out)[tag]
    names = CT.channels(out, entry['entry_id'])
    img = Channels(entry['image_path'])
    if img.shape[0] != len(names):
        raise SystemExit(f'{tag}: image has {img.shape[0]} channels, QuPath entry names {len(names)}')
    H, W = img.shape[1:]

    t0 = time.time()
    features = json.loads((out / 'cell_geojson' / f'{tag}.geojson').read_text())
    features = features['features'] if isinstance(features, dict) else features
    mask, nucm, oids, has_nucleus = rasterize(features, (H, W))
    del features
    n = len(oids)
    log(tag, f'rasterized {n:,} cells on {H}x{W} in {time.time() - t0:.0f}s')

    keep_cols = ['Object ID', 'Centroid X µm', 'Centroid Y µm', 'Cell: Area µm^2', 'Cell: Solidity'] + \
        [f'Cell: {c}: Mean' for c in names]
    cells = pd.read_parquet(out / 'cells' / f'{tag}.parquet', columns=keep_cols)
    missing = set(oids) - set(cells['Object ID'])
    if missing or len(cells) != n:
        raise SystemExit(f'{tag}: {len(missing)} geometry ids not in the measurement table '
                         f'({n:,} features, {len(cells):,} rows)')
    sizes0 = np.bincount(mask.ravel(), minlength=n + 1)[1:].astype(np.float64)
    nuc0 = np.bincount(nucm.ravel(), minlength=n + 1)[1:].astype(np.float64)
    qc, qc_summary = geometry_qc(oids, has_nucleus, cells, sizes0, nuc0)
    excluded = np.flatnonzero(~qc['qc_spillover_context_eligible'].to_numpy()) + 1
    if len(excluded):
        keep_label = np.ones(n + 1, bool)
        keep_label[excluded] = False
        mask = np.where(keep_label[mask], mask, 0)
        nucm = np.where(keep_label[nucm], nucm, 0)
    log(tag, f'QC: {qc_summary["analysis_eligible"]:,} analysis-eligible, {len(excluded):,} out of spillover '
             f'context; {qc_summary["reasons"]}')

    sizes = np.bincount(mask.ravel(), minlength=n + 1)[1:].astype(np.float64)
    inner = find_boundaries(mask, mode='inner', connectivity=2, background=0)
    band = binary_dilation(inner, diamond(CT.REDSEA['edge_radius'])) & (mask > 0)
    del inner
    band_idx = np.flatnonzero(band.ravel())
    del band
    mask_flat = mask.ravel()
    mask_band = mask_flat[band_idx]
    bandcount = np.bincount(mask_band, minlength=n + 1)[1:].astype(np.float64)
    nuc_flat = nucm.ravel()
    nuc_band = nuc_flat[band_idx]
    nuccount = np.bincount(nuc_flat, minlength=n + 1)[1:].astype(np.float64)
    nucbandcount = np.bincount(nuc_band, minlength=n + 1)[1:].astype(np.float64)

    C = len(names)
    data, edge = np.zeros((n, C)), np.zeros((n, C))
    nuc, nuc_edge = np.zeros((n, C)), np.zeros((n, C))
    data_raw = np.zeros((n, C))                       # before the despeckle, for the alignment check only
    ey, ex = empty_windows(out, tag, H, W)
    empty = {}
    t0 = time.time()
    for k in range(C):
        raw = img.native(k)
        data_raw[:, k] = np.bincount(mask_flat, weights=raw.ravel().astype(np.float32), minlength=n + 1)[1:]
        clean = despeckle(raw)
        if len(ey):
            wm, wr = window_means(clean, ey, ex), window_means(raw, ey, ex)
            empty[names[k]] = {'windows': len(wm), 'null_p50': float(np.percentile(wm, 50)),
                               'null_p99': float(np.percentile(wm, 99)), 'null_p999': float(np.percentile(wm, 99.9)),
                               'null_p99_before_despeckle': float(np.percentile(wr, 99))}
        ch = clean.astype(np.float32).ravel()
        del raw, clean
        data[:, k] = np.bincount(mask_flat, weights=ch, minlength=n + 1)[1:]
        edge[:, k] = np.bincount(mask_band, weights=ch[band_idx], minlength=n + 1)[1:]
        nuc[:, k] = np.bincount(nuc_flat, weights=ch, minlength=n + 1)[1:]
        nuc_edge[:, k] = np.bincount(nuc_band, weights=ch[band_idx], minlength=n + 1)[1:]
        del ch
    img.close()
    del nucm, nuc_flat, nuc_band, mask_band, band_idx
    log(tag, f'{C} channels despeckled and summed in {time.time() - t0:.0f}s')

    t0 = time.time()
    M = contact_matrix_gapbridge(mask, n, CT.REDSEA['gap_bridge'])
    del mask, mask_flat
    log(tag, f'contact matrix {M.nnz:,} in {time.time() - t0:.0f}s')

    corrected = np.array([c not in CT.NO_SPILLOVER for c in names])
    alpha = np.where(corrected, CT.REDSEA['alpha'], 0.0)
    comp_mode = np.where(corrected, CT.REDSEA['comp_mode'], 0)
    means, counts, n_isolated = compensate_compartments(data, edge, nuc, nuc_edge, sizes, bandcount,
                                                        nuccount, nucbandcount, M, comp_mode, alpha)

    # alignment check: the uncorrected whole-cell mean (before the despeckle) against QuPath's own, per channel
    with np.errstate(invalid='ignore', divide='ignore'):
        uncorr = data_raw / np.where(sizes > 0, sizes, np.nan)[:, None]
        despeckled = data / np.where(sizes > 0, sizes, np.nan)[:, None]
    q = cells.set_index('Object ID').reindex(oids)
    # The check catches a key or channel mismatch, so it reads the cells whose raster and polygon agree in size:
    # on cells of a few pixels the two pixel sets differ by construction (ND_6232 CD35, r 0.73 over all cells).
    flags = qc['qc_flags'].astype(str)
    ok = (qc['qc_analysis_eligible'].to_numpy()
          & ~flags.str.contains('raster_area_mismatch|cell_area_below_min').to_numpy())
    align = {}
    for k, c in enumerate(names):
        a, b = uncorr[ok, k], q[f'Cell: {c}: Mean'].to_numpy(float)[ok]
        f = np.isfinite(a) & np.isfinite(b)
        r = float(np.corrcoef(a[f], b[f])[0, 1]) if f.sum() > 2 else float('nan')
        align[c] = {'r': r, 'median_ratio': float(np.nanmedian(a[f] / np.where(b[f] > 0, b[f], np.nan)))}
    bad = {c: v for c, v in align.items() if not v['r'] >= ALIGN_MIN_R}
    if bad:
        raise SystemExit(f'{tag}: raster means disagree with QuPath on {bad}')
    # what the despeckle changed: per channel, the median whole-cell mean after over before (analysis cells)
    effect = {c: float(np.nanmedian(despeckled[ok, k] / np.where(uncorr[ok, k] > 0, uncorr[ok, k], np.nan)))
              for k, c in enumerate(names)}

    frame = qc.copy()
    for comp in COMPARTMENTS:
        m = means[comp].astype(np.float32)
        for k, c in enumerate(names):
            frame[f'{comp}__{c}'] = m[:, k]
        frame[f'{comp}__n_px'] = counts[comp].astype(np.float32)
    dest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(dest, index=False)
    rec = {'image': tag, 'image_path': entry['image_path'], 'channels': names, 'n_cells': n,
           'redsea': CT.REDSEA, 'despeckle': {**DESPECKLE, 'median_cell_mean_ratio': effect},
           'empty_area_null': {'method': EMPTY, 'per_channel': empty},
           'uncorrected_channels': [c for c in names if c in CT.NO_SPILLOVER],
           'isolated_cells': n_isolated, 'contacts': int(M.nnz), 'geometry_qc': qc_summary,
           'alignment_check': {'min_r': min(v['r'] for v in align.values()), 'per_channel': align},
           'seconds': round(time.time() - t_all)}
    rec_path.write_text(json.dumps(rec, indent=1))
    log(tag, f'done in {rec["seconds"]}s; alignment r >= {rec["alignment_check"]["min_r"]:.4f}; '
             f'{n_isolated:,} isolated cells')
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--dest', type=Path, default=None, help='default: $OUT/typing/redsea (trials: a scratch folder)')
    a = ap.parse_args()
    if a.dest:
        a.dest.mkdir(parents=True, exist_ok=True)
    tags = [r['image_name'] for r in CT.cohort_rows()] if a.images == ['all'] else a.images
    if a.workers <= 1 or len(tags) == 1:
        for t in tags:
            process(a.out, t, a.force, a.dest)
        write_empty_nulls(a.dest or CT.typing_dir(a.out) / 'redsea')
        return 0
    with ProcessPoolExecutor(a.workers) as pool:
        futures = {t: pool.submit(process, a.out, t, a.force, a.dest) for t in tags}
        failed = []
        for t, f in futures.items():
            try:
                f.result()
            except BaseException as e:   # noqa: BLE001 -- report every image, then fail
                print(f'[{t}] FAILED: {e}', flush=True)
                failed.append(t)
    print(f'{len(tags) - len(failed)} of {len(tags)} images done' + (f'; failed: {failed}' if failed else ''))
    write_empty_nulls(a.dest or CT.typing_dir(a.out) / 'redsea')
    return 1 if failed else 0


def write_empty_nulls(redsea_dir: Path) -> None:
    """typing/empty_area_nulls.csv (levels.empty_nulls), from every image's record in `redsea_dir`."""
    rows = []
    for t in [r['image_name'] for r in CT.cohort_rows()]:
        f = redsea_dir / f'{t}.json'
        if f.exists():
            per = json.loads(f.read_text()).get('empty_area_null', {}).get('per_channel', {})
            rows += [{'image': t, 'marker': c, **v} for c, v in per.items()]
    if rows:
        pd.DataFrame(rows).to_csv(redsea_dir.parent / 'empty_area_nulls.csv', index=False)


if __name__ == '__main__':
    sys.exit(main())
