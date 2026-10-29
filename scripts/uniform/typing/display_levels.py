#!/usr/bin/env python3
"""Display-only levels: per-image positivity of the panel markers the typing does not read (registry role 'display').

Owner 2026-10-29, on KRT14, CD44, NaKATPase and BActin: "just use the ecad methods". The app needs every panel marker's
value and positive flag (owner 2026-10-20); these levels give the flag and the value scale's floor and never name a cell.

E-cadherin's methods (levels.py, positivity_class majority_epithelial): the expressing population's lower edge
('edge': its peak minus POSITIVITY_K right-side widths) or the midpoint of a two-component Gaussian mixture ('gate'), each
floored at the image's empty-slide null (typing/empty_area_nulls.csv), on the whole cell or the membrane band. Selection
(agent determination 2026-10-29; E-cadherin's checks are against leukocytes, which carry CD44, NaKATPase and BActin, so
they do not apply): the first option, in the order whole cell / membrane and edge / gate, whose level is finite on every
image and leaves both positive and negative cells on every image.

Writes typing/display_levels.csv (levels.csv's columns) and typing/display_levels_options.csv (every option).

    python scripts/uniform/typing/display_levels.py --run data/development/uniform-20donor-20261027 --workers 10
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import levels as LV  # noqa: E402
import registry as RG  # noqa: E402

COLUMNS = ['image', 'marker', 'measurement', 'method', 'center', 'width', 'level_log', 'level_raw', 'positive_fraction',
           'zero_fraction']


def one(run: Path, tag: str, nulls: dict) -> list[dict]:
    reg = RG.load()
    k = reg.const('POSITIVITY_K')
    markers = reg.markers('display')
    cols = sorted({reg.column(m, meas) for m in markers for meas in reg.candidates(m)})
    d = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet', columns=['qc_analysis_eligible'] + cols)
    d = d[d['qc_analysis_eligible']].reset_index(drop=True)
    rows = []
    for m in markers:
        for meas, method in LV.options(reg, m):
            lv = LV.level_for(d, reg, m, meas, method, k, nulls.get((tag, m), np.nan))
            x = d[reg.column(m, meas)].to_numpy(float)
            rows.append({'image': tag, 'marker': m, 'measurement': meas, 'method': method,
                         **{c: lv.get(c, np.nan) for c in ('center', 'width', 'level_log')},
                         'level_raw': float(np.expm1(lv['level_log'])) if np.isfinite(lv['level_log']) else np.nan,
                         'positive_fraction': float(np.mean(lv['pos'])) if np.isfinite(lv['level_log']) else np.nan,
                         'zero_fraction': float(np.mean(~(x > 0)))})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=10)
    a = ap.parse_args()
    tags = [r['image_name'] for r in CT.cohort_rows()]
    nulls = LV.empty_nulls(a.run)
    with ProcessPoolExecutor(a.workers) as ex:
        opts = pd.DataFrame([r for rs in ex.map(one, repeat(a.run), tags, repeat(nulls)) for r in rs])
    typ = a.run / 'typing'
    opts.to_csv(typ / 'display_levels_options.csv', index=False)
    reg = RG.load()
    chosen = []
    for m in reg.markers('display'):
        for meas, method in LV.options(reg, m):
            o = opts[(opts['marker'] == m) & (opts['measurement'] == meas) & (opts['method'] == method)]
            works = len(o) == len(tags) and np.isfinite(o['level_log']).all() and \
                o['positive_fraction'].between(0, 1, inclusive='neither').all()
            print(f'{m:10s} {meas:10s} {method:5s} positive {o["positive_fraction"].min():.3f}-'
                  f'{o["positive_fraction"].max():.3f}  {"works on every image" if works else "fails"}', flush=True)
            if works:
                chosen.append(o)
                break
        else:
            print(f'STOP: no option works on every image for {m}', flush=True)
            return 2
    pd.concat(chosen)[COLUMNS].to_csv(typ / 'display_levels.csv', index=False)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
