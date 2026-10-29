#!/usr/bin/env python3
"""Area coverage: per cell and marker, the proportion of the cell's pixels at or above the marker's level.

Owner 2026-10-29: "There are two markers for a T Cell CD8a and both are brighter and cover more cell area than any
myeloid marker so this should not be a conflicting call." A marker a cell expresses covers the cell; a neighbor's signal
that REDSEA leaves behind sits at the contact rim. identity.py uses this measure, with brightness (fold over level), to
settle a canary conflict between two groups of equal depth.

Pixels: the despeckled channel (redsea_cells.despeckle), uncompensated. Level: the image's level of the marker
(typing/levels.csv level_raw, on the marker's measurement). Measured over the measurement's pixels: the whole cell, or
REDSEA's 1-pixel inner boundary band for a marker the level selection reads on the membrane.

Interior coverage (covint__<marker>, owner 2026-10-29 "Interior coverage, clear lead"): the same proportion over the
cell's interior, the pixels at least INTERIOR_PX from its outline, for every marker whatever its measurement. A
neighbor's signal left at the contact rim does not reach the interior; on 127 tie-settled cells checked against the
image (output/typing_accuracy_20261029/), the difference in interior coverage separated "the placed group fills the
cell" from "both or neither fill it" at AUC 0.87, whole-cell coverage at 0.76.

Writes typing/coverage/<image>.parquet (object_id, cov__<marker>, covint__<marker>) for every typing marker with a
finite level.

    python scripts/uniform/typing/area_coverage.py --run data/development/uniform-20donor-20261027 --workers 5
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter, minimum_filter
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import redsea_cells as RS  # noqa: E402
import registry as RG  # noqa: E402

PREFIX = 'cov__'
INTERIOR_PREFIX = 'covint__'
INTERIOR_PX = 2
INTERIOR_COUNT = 'n_interior_pixels'
SCHEMA_VERSION = 1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def input_fingerprints(run: Path, tag: str) -> dict:
    """Content hashes for small/derived inputs; source TIFF identity by path, size and nanosecond mtime.

    TIFFs are many GB, so their filesystem identity is recorded without rereading every channel for validation.
    Relative keys allow a copied run to retain provenance when its inputs are unchanged.
    """
    reg = RG.load()
    entry = CT.manifest(run)[tag]
    sources = {name: run / name for name in (
        'typing/levels.csv', 'typing/levels_selection.csv', f'cell_geojson/{tag}.geojson',
        f'qupath/data/{entry["entry_id"]}/server.json')}
    sources.update({'registry': RG.REGISTRY_DIR / f'{reg.name}.csv',
                    'coverage_code': Path(__file__), 'pixel_processing_code': Path(RS.__file__)})
    image = Path(entry['image_path']).resolve()
    st = image.stat()
    return {'sha256': {name: sha256(path) for name, path in sources.items()},
            'image': {'path': str(image), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns},
            'channels': CT.channels(run, entry['entry_id'])}


def validate_frame(frame: pd.DataFrame, object_ids: pd.Series, markers: list[str], tag: str) -> None:
    required = ['object_id', INTERIOR_COUNT] + [INTERIOR_PREFIX + m for m in markers]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f'{tag}: coverage missing columns: {missing}; rerun area_coverage.py')
    ids = frame['object_id']
    if ids.isna().any() or ids.duplicated().any() or object_ids.isna().any() or object_ids.duplicated().any():
        raise ValueError(f'{tag}: coverage/source object IDs must be non-null and unique')
    missing_ids, extra_ids = set(object_ids) - set(ids), set(ids) - set(object_ids)
    if missing_ids or extra_ids:
        raise ValueError(f'{tag}: coverage object IDs differ from REDSEA: {len(missing_ids)} missing, '
                         f'{len(extra_ids)} extra')
    count = frame[INTERIOR_COUNT].to_numpy(float)
    if not np.isfinite(count).all() or (count < 0).any() or (count != np.floor(count)).any():
        raise ValueError(f'{tag}: invalid interior pixel count')
    values = frame[[INTERIOR_PREFIX + m for m in markers]].to_numpy(float)
    has = count > 0
    if not np.isfinite(values[has]).all() or ((values[has] < 0) | (values[has] > 1)).any():
        raise ValueError(f'{tag}: invalid interior coverage for cells with interior pixels')
    if not np.isnan(values[~has]).all():
        raise ValueError(f'{tag}: cells with no interior must have NaN interior coverage')


def write_coverage(run: Path, tag: str, frame: pd.DataFrame, fingerprints: dict | None = None) -> None:
    markers = [c.removeprefix(INTERIOR_PREFIX) for c in frame if c.startswith(INTERIOR_PREFIX)]
    ids = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet', columns=['object_id'])['object_id']
    validate_frame(frame, ids, markers, tag)
    current = input_fingerprints(run, tag)
    if fingerprints is not None and fingerprints != current:
        raise ValueError(f'{tag}: coverage inputs changed during computation')
    dest = run / 'typing' / 'coverage'
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / f'{tag}.parquet'
    tmp = path.with_suffix('.parquet.tmp')
    frame.to_parquet(tmp, index=False)
    tmp.replace(path)
    metadata = {'schema_version': SCHEMA_VERSION, 'image': tag, 'interior_px': INTERIOR_PX,
                'rows': len(frame), 'markers': markers, 'inputs': current, 'parquet_sha256': sha256(path),
                'cells_without_interior': int((frame[INTERIOR_COUNT] == 0).sum())}
    meta = path.with_suffix('.json')
    temp_meta = meta.with_suffix('.json.tmp')
    temp_meta.write_text(json.dumps(metadata, indent=2))
    temp_meta.replace(meta)


def read_validated(run: Path, tag: str, object_ids: pd.Series, markers: list[str]) -> pd.DataFrame:
    path = run / 'typing' / 'coverage' / f'{tag}.parquet'
    meta = path.with_suffix('.json')
    if not meta.exists():
        raise ValueError(f'{tag}: coverage provenance missing; rerun area_coverage.py')
    metadata = json.loads(meta.read_text())
    if (metadata.get('schema_version') != SCHEMA_VERSION or metadata.get('interior_px') != INTERIOR_PX
            or metadata.get('image') != tag or metadata.get('inputs') != input_fingerprints(run, tag)):
        raise ValueError(f'{tag}: stale coverage inputs or method; rerun area_coverage.py')
    if metadata.get('parquet_sha256') != sha256(path):
        raise ValueError(f'{tag}: coverage Parquet does not match its provenance')
    frame = pd.read_parquet(path)
    ids = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet', columns=['object_id'])['object_id']
    validate_frame(frame, ids, markers, tag)
    if len(frame) != metadata.get('rows') or not set(markers) <= set(metadata.get('markers', [])):
        raise ValueError(f'{tag}: coverage metadata rows or markers do not match')
    if object_ids.isna().any() or object_ids.duplicated().any() or not set(object_ids) <= set(ids):
        raise ValueError(f'{tag}: identity input object IDs must be unique and present in coverage')
    return frame.set_index('object_id').loc[object_ids].reset_index()


def one(run: Path, tag: str) -> str:
    fingerprints = input_fingerprints(run, tag)
    reg = RG.load()
    entry = CT.manifest(run)[tag]
    names = CT.channels(run, entry['entry_id'])
    lv = pd.read_csv(run / 'typing' / 'levels.csv')
    lv = lv[(lv['image'] == tag) & np.isfinite(lv['level_raw'])].set_index('marker')
    img = RS.Channels(entry['image_path'])
    H, W = img.shape[1:]
    feats = json.loads((run / 'cell_geojson' / f'{tag}.geojson').read_text())
    feats = feats['features'] if isinstance(feats, dict) else feats
    mask, _, oids, _ = RS.rasterize(feats, (H, W))
    del feats
    n = len(oids)
    flat = mask.ravel()
    size = np.bincount(flat, minlength=n + 1)[1:].astype(float)
    band = np.flatnonzero(find_boundaries(mask, mode='inner', connectivity=2, background=0).ravel())
    w = 2 * INTERIOR_PX + 1                     # a pixel is interior when its whole w x w window is the same cell
    interior = np.flatnonzero(((minimum_filter(mask, size=w) == mask) & (maximum_filter(mask, size=w) == mask)
                               & (mask > 0)).ravel())
    del mask
    iflat = flat[interior]
    isize = np.bincount(iflat, minlength=n + 1)[1:].astype(float)
    bflat = flat[band]
    bsize = np.bincount(bflat, minlength=n + 1)[1:].astype(float)
    frame = {'object_id': oids, INTERIOR_COUNT: isize.astype(np.int64)}
    for m in reg.typing_markers():
        if m not in lv.index:
            continue
        ch = reg.table.loc[m, 'channel']
        ch = ch if ch in names else m               # the registry's PGP9_5 is PGP9.5 in the images
        if ch not in names:
            raise ValueError(f'{tag}: no image channel for required coverage marker {m}')
        above = RS.despeckle(img.native(names.index(ch))).ravel() >= lv.loc[m, 'level_raw']
        with np.errstate(invalid='ignore', divide='ignore'):
            if lv.loc[m, 'measurement'] == 'membrane':
                frame[PREFIX + m] = np.bincount(bflat, weights=above[band], minlength=n + 1)[1:] / \
                    np.where(bsize > 0, bsize, np.nan)
            else:
                frame[PREFIX + m] = np.bincount(flat, weights=above, minlength=n + 1)[1:] / \
                    np.where(size > 0, size, np.nan)
            frame[INTERIOR_PREFIX + m] = np.bincount(iflat, weights=above[interior], minlength=n + 1)[1:] / \
                np.where(isize > 0, isize, np.nan)
        del above
    img.close()
    write_coverage(run, tag, pd.DataFrame(frame), fingerprints)
    return tag


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--images', nargs='*', default=[])
    ap.add_argument('--workers', type=int, default=5)
    a = ap.parse_args()
    tags = a.images or [r['image_name'] for r in CT.cohort_rows()]
    with ProcessPoolExecutor(a.workers) as ex:
        for t in ex.map(one, [a.run] * len(tags), tags):
            print(t, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
