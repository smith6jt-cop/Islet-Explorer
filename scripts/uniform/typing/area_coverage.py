#!/usr/bin/env python3
"""Area coverage: per cell and marker, the proportion of the cell's pixels at or above the marker's level.

Owner 2026-10-29: "There are two markers for a T Cell CD8a and both are brighter and cover more cell area than any
myeloid marker so this should not be a conflicting call." A marker a cell expresses covers the cell; a neighbor's signal
that REDSEA leaves behind sits at the contact rim. identity.py uses this measure, with brightness (fold over level), to
settle a canary conflict between two groups of equal depth.

Pixels: the despeckled channel (redsea_cells.despeckle), uncompensated. Level: the image's level of the marker
(typing/levels.csv level_raw, on the marker's measurement). Measured over the measurement's pixels: the whole cell, or
REDSEA's 1-pixel inner boundary band for a marker the level selection reads on the membrane.

Writes typing/coverage/<image>.parquet (object_id, cov__<marker>) for every typing marker with a finite level.

    python scripts/uniform/typing/area_coverage.py --run data/development/uniform-20donor-20261027 --workers 5
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.segmentation import find_boundaries

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import redsea_cells as RS  # noqa: E402
import registry as RG  # noqa: E402

PREFIX = 'cov__'


def one(run: Path, tag: str) -> str:
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
    del mask
    bflat = flat[band]
    bsize = np.bincount(bflat, minlength=n + 1)[1:].astype(float)
    frame = {'object_id': oids}
    for m in reg.typing_markers():
        if m not in lv.index:
            continue
        ch = reg.table.loc[m, 'channel']
        ch = ch if ch in names else m               # the registry's PGP9_5 is PGP9.5 in the images
        if ch not in names:
            continue
        above = RS.despeckle(img.native(names.index(ch))).ravel() >= lv.loc[m, 'level_raw']
        with np.errstate(invalid='ignore', divide='ignore'):
            if lv.loc[m, 'measurement'] == 'membrane':
                frame[PREFIX + m] = np.bincount(bflat, weights=above[band], minlength=n + 1)[1:] / \
                    np.where(bsize > 0, bsize, np.nan)
            else:
                frame[PREFIX + m] = np.bincount(flat, weights=above, minlength=n + 1)[1:] / \
                    np.where(size > 0, size, np.nan)
        del above
    img.close()
    out = run / 'typing' / 'coverage'
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(frame).to_parquet(out / f'{tag}.parquet', index=False)
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
