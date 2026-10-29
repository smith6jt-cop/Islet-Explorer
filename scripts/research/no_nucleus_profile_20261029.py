#!/usr/bin/env python3
"""Marker profile of the detected cells without a nucleus (owner question 2026-10-29), against the nucleated cells.

The typing drops a cell without a nucleus before measuring (owner 2026-10-16), so these cells have no REDSEA value.
Here every cell of every image is measured the same way: the despeckled channel (redsea_cells.despeckle), the whole-cell
mean, with no spillover compensation (a cell without a nucleus is outside REDSEA's spillover context). Positivity is
the typing's per-image level of each marker; markers the typing reads on the membrane are also measured on
the 1-pixel inner boundary band (REDSEA's band), so the cells can be put through identity.py (typing/levels.csv; those levels are on compensated values, so both groups
are read against them the same way and compared with each other, not with the typing's calls).

Writes output/typing_final_checks_20261028/no_nucleus/<image>.parquet (object_id, has_nucleus, area, marker means).

    python scripts/research/no_nucleus_profile_20261029.py --workers 5
"""
from __future__ import annotations

import argparse, json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.segmentation import find_boundaries

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import redsea_cells as RS  # noqa: E402
import registry as RG  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'typing_final_checks_20261028' / 'no_nucleus'
ONLY: set[str] = set()


def one(tag: str) -> str:
    reg = RG.load()
    entry = CT.manifest(RUN)[tag]
    names = CT.channels(RUN, entry['entry_id'])
    img = RS.Channels(entry['image_path'])
    H, W = img.shape[1:]
    feats = json.loads((RUN / 'cell_geojson' / f'{tag}.geojson').read_text())
    feats = feats['features'] if isinstance(feats, dict) else feats
    mask, _, oids, has_nucleus = RS.rasterize(feats, (H, W))
    del feats
    n = len(oids)
    flat = mask.ravel()
    size = np.bincount(flat, minlength=n + 1)[1:].astype(float)
    frame = {'object_id': oids, 'has_nucleus': has_nucleus, 'area_um2': size * CT.PIXEL_UM ** 2}
    band = np.flatnonzero(find_boundaries(mask, mode='inner', connectivity=2, background=0).ravel())
    bflat = flat[band]
    bsize = np.bincount(bflat, minlength=n + 1)[1:].astype(float)
    lv = pd.read_csv(RUN / 'typing' / 'levels.csv')
    membrane = set(lv.loc[(lv['image'] == tag) & (lv['measurement'] == 'membrane'), 'marker'])
    for m in reg.typing_markers():
        ch = reg.table.loc[m, 'channel']
        ch = ch if ch in names else m           # the registry's PGP9_5 is PGP9.5 in the images
        if ch not in names or (ONLY and m not in ONLY):
            continue
        a = RS.despeckle(img.native(names.index(ch))).astype(np.float32).ravel()
        with np.errstate(invalid='ignore', divide='ignore'):
            frame[m] = np.bincount(flat, weights=a, minlength=n + 1)[1:] / np.where(size > 0, size, np.nan)
            if m in membrane:
                frame[f'membrane__{m}'] = np.bincount(bflat, weights=a[band], minlength=n + 1)[1:] / \
                    np.where(bsize > 0, bsize, np.nan)
    img.close()
    OUT.mkdir(parents=True, exist_ok=True)
    new = pd.DataFrame(frame)
    if ONLY and (OUT / f'{tag}.parquet').exists():      # add these markers to an earlier measurement
        old = pd.read_parquet(OUT / f'{tag}.parquet')
        keep = ('object_id', 'has_nucleus', 'area_um2')
        new = old.drop(columns=[c for c in new.columns if c in old.columns and c not in keep]).merge(
            new.drop(columns=['has_nucleus', 'area_um2']), on='object_id', validate='one_to_one')
    new.to_parquet(OUT / f'{tag}.parquet', index=False)
    return tag


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=5)
    ap.add_argument('--only', nargs='*', default=[], help='measure only these markers and add them to the saved tables')
    a = ap.parse_args()
    ONLY.update(a.only)
    tags = [r['image_name'] for r in CT.cohort_rows()]
    with ProcessPoolExecutor(a.workers) as ex:
        for t in ex.map(one, tags):
            print(t, flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
