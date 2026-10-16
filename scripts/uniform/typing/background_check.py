#!/usr/bin/env python3
"""background_check.py -- is each marker's ACTUAL background aligned across the 20 images?

Owner, 2026-10-16: "there is a difference between actual background and variation between images."

UniFORM's acceptance bar (align_scale.py, Senior's BAR_MODE_SD) is the spread of each image's histogram MODE after
the shift. The mode is background only when most cells do not carry the marker. For a marker most cells carry
(E-cadherin, panCK, beta-catenin, Na/K-ATPase, vimentin, CK19) the mode is the positive population, and its spread
between images is composition and staining biology, not misalignment. This check measures the background
directly: each marker's level in cells that clearly belong to a lineage that does not carry it, chosen on OTHER
markers only (never on the marker being measured), per image, before and after the shift.

  NEGATIVE_BY[group] -> the selecting markers and their minimum scaled value; a cell is in the negative population
  when any selecting marker is at or above that value. Epithelial cells (E-cadherin or panCK >= 0.5) are the
  negatives for hormone, immune, vascular, stromal and neural markers; immune cells (CD45, CD3e or CD68 >= 0.5) for
  epithelial markers. Hormones also require the cell to be outside every endocrine object. Markers every lineage
  carries (BActin, BCatenin, NaKATPase) have no negative population and are reported as not assessable.

Background spread = SD over images of log(median of the negative population), the quantity Senior's bar is stated
on (0.10). Written to `$OUT/typing/background_check.csv`.

    python scripts/uniform/typing/background_check.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402

EPITHELIAL = (('ECAD', 'panCK'), 0.5)
IMMUNE = (('CD45', 'CD3e', 'CD68'), 0.5)
GROUPS = {
    'hormone': (('INS', 'GCG', 'SST', 'CHGA', 'Ghrelin', 'PGP9.5'), EPITHELIAL, True),
    'immune': (('CD3e', 'CD4', 'CD8a', 'CD20', 'CD68', 'CD163', 'CD45', 'HLADR', 'CD35'), EPITHELIAL, False),
    'vascular_stromal_neural': (('CD31', 'CD34', 'PDPN', 'SMA', 'VIM', 'ColIV', 'B3TUBB', 'GAP43', 'CD44', 'LGALS3',
                                 'Ki67'), EPITHELIAL, False),
    'epithelial': (('ECAD', 'panCK', 'CK19', 'KRT14', 'KRT8-18', 'epCAM'), IMMUNE, False),
}
NOT_ASSESSABLE = ('BActin', 'BCatenin', 'NaKATPase')
BAR = 0.10


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    typ = CT.typing_dir(a.out)
    scale = CL.load_scale(a.out)
    shifts = json.loads((typ / 'shifts.json').read_text())['shifts']
    tags = [r['image_name'] for r in CT.cohort_rows()]
    rows = []
    med = {}
    for t in tags:
        al = pd.read_parquet(typ / 'aligned' / f'{t}.parquet')
        ctx = pd.read_parquet(typ / 'context' / f'{t}.parquet', columns=['object_id', 'in_islet'])
        al = al.merge(ctx, on='object_id', how='left')
        al = al[al['qc_analysis_eligible']]
        for group, (markers, (sel, level), outside_objects) in GROUPS.items():
            sel = [m for m in sel if m in al.columns]
            s = pd.DataFrame(CL.apply_scale(al[sel].to_numpy(), sel, scale), columns=sel)
            neg = (s >= level).any(axis=1).to_numpy()
            if outside_objects:
                neg &= ~al['in_islet'].fillna(False).to_numpy(bool)
            for m in markers:
                if m not in al.columns:
                    continue
                v = al.loc[neg, m].to_numpy(float)
                v = v[np.isfinite(v) & (v > 0)]
                factor = shifts[m]['scale'][t] if isinstance(shifts[m]['scale'].get(t), float) else 1.0
                med[(m, t)] = (float(np.median(v)) if len(v) else np.nan, factor, int(len(v)))
    for group, (markers, _, _) in GROUPS.items():
        for m in markers:
            vals = [(t, *med[(m, t)]) for t in tags if (m, t) in med]
            after = np.array([v for _, v, _, _ in vals])
            before = np.array([v / f for _, v, f, _ in vals])
            ok = np.isfinite(after) & (after > 0)
            sd_after = float(np.log(after[ok]).std())
            sd_before = float(np.log(before[ok]).std())
            rows.append({'marker': m, 'group': group, 'images': int(ok.sum()),
                         'negative_cells_median': int(np.median([n for *_, n in vals])),
                         'background_sd_log_before': round(sd_before, 3), 'background_sd_log_after': round(sd_after, 3),
                         'background_fold_after': round(float(np.exp(np.log(after[ok]).max() - np.log(after[ok]).min())), 2),
                         'mode_sd_log_after': round(shifts[m].get('mode_sd_after', np.nan), 3),
                         'background_aligned': sd_after <= BAR,
                         'worst_images': ', '.join(t for t, v, _, _ in sorted(
                             [x for x in vals if np.isfinite(x[1])],
                             key=lambda x: abs(np.log(x[1]) - np.median(np.log(after[ok]))))[-2:])})
    for m in NOT_ASSESSABLE:
        rows.append({'marker': m, 'group': 'every lineage', 'images': 0, 'background_aligned': None,
                     'mode_sd_log_after': round(shifts[m].get('mode_sd_after', np.nan), 3)})
    df = pd.DataFrame(rows)
    df.to_csv(typ / 'background_check.csv', index=False)
    print(df.to_string(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
