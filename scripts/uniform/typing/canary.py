#!/usr/bin/env python3
"""canary.py -- each typing marker's canary detection level per image (final typing plan of 2026-10-28, Step 5).

A marker counts toward canary depth (identity.py) only above a stricter level than plain positivity, as the Xenium
canary step counts a gene only above a stricter level than its ordinary detection. The plan set that level at 0-3
background widths above the positivity level, chosen on the 7 ND images against a null (each marker's positive rate in
compartments that cannot carry it, in cells farther than 20 um from any structure object). Two parts of that failed on
these data (agent determination 2026-10-28, under the owner's rule "flag the methods you are using so they will work
with the images and markers"): a background width cannot be read for every marker and image (CD20 on 6 images,
E-cadherin on 2, CD68 on 2: the background is a plateau), and leukocytes farther than 20 um from any structure number
82-1,395 on the ND images, too few for a per-image percentile. The method used instead needs neither:

  canary level = max(the marker's positivity level, the CANARY_NULL_QUANTILE-th percentile of the marker among the
                 image's cells of every compartment that cannot carry it)

the known-negative rule the owner approved for SMA (Barlow et al. eLife 2024), applied to the canary step. Compartments
that can carry a marker: its registry compartment and shared_by. Reference cells (not isolated: a neighbor's
spillover is what the canary level must discount): Epithelial = endocrine-object hormone-positive cells and ordinary-
tissue cells positive for an Epithelial majority marker that are not leukocytes; Hematopoietic = leukocytes (the
Leukocyte marker with a T or Myeloid marker). A marker carried by both has its registry check reference instead (e.g.
T cells). A reference needs KNOWN_NEGATIVE_MIN_ABOVE cells above the percentile, or the run stops.

Writes typing/canary_levels.csv: image, marker, level_log, canary_log, canary_fold_over_level, null_reference,
null_cells, null_rate_at_level (the plain level's rate in the null).

    python scripts/uniform/typing/canary.py --run $OUT --workers 20
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import registry as RG  # noqa: E402
import levels as LV  # noqa: E402
import identity as ID  # noqa: E402

REFERENCE_COMPARTMENTS = ('Epithelial', 'Hematopoietic')


def image_levels(run: Path, tag: str) -> list[dict]:
    reg = RG.load()
    R = LV.roles(reg)
    markers = ID.used_markers(run, reg)
    lv = pd.read_csv(run / 'typing' / 'levels.csv')
    lv = lv[lv['image'] == tag].set_index('marker')
    d = ID.load_values(run, tag, reg, markers, lv)
    n = len(d)
    x = {m: np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None) for m in markers}
    pos = {m: np.log1p(x[m]) > lv.loc[m, 'level_log'] for m in markers}
    in_endo = d['endocrine_object_class'].notna().to_numpy()
    ordinary = d['context'].eq('tissue').to_numpy() & ~in_endo
    hormone = LV.any_of(pos, R['hormones'], n) & in_endo
    leuk = LV.any_of(pos, R['leukocyte'], n) & LV.any_of(pos, R['t'] + R['myeloid'], n)
    refs = {'Epithelial': hormone | (ordinary & LV.any_of(pos, R['epithelial'], n) & ~leuk), 'Hematopoietic': leuk}
    named = {'T cells': LV.any_of(pos, R['t'], n) & LV.any_of(pos, R['leukocyte'], n),
             'endocrine-object hormone-positive cells': hormone}
    q = reg.const('CANARY_NULL_QUANTILE')
    need = reg.const('KNOWN_NEGATIVE_MIN_ABOVE')
    rows = []
    for m in markers:
        carriers = {reg.table.loc[m, 'compartment']} | set(RG.split(reg.table.loc[m, 'shared_by']))
        nulls = [c for c in REFERENCE_COMPARTMENTS if c not in carriers]
        if nulls:
            ref = np.zeros(n, bool)
            for c in nulls:
                ref |= refs[c]
            name = ' + '.join(nulls)
        else:
            name = next((r for r in RG.split(reg.table.loc[m, 'reference']) if r.strip() in named), '').strip()
            if not name:
                raise SystemExit(f'{m}: no reference compartment cannot carry it and no named check reference')
            ref = named[name]
        o = x[m][ref]
        if len(o) * (1 - q / 100) < need:
            raise SystemExit(f'{tag} {m}: {len(o)} reference cells are too few for the {q:g}th percentile')
        level = float(lv.loc[m, 'level_log'])
        canary = max(level, float(np.log1p(np.percentile(o, q))))
        rows.append({'image': tag, 'marker': m, 'level_log': level, 'canary_log': canary,
                     'canary_fold_over_level': float(np.expm1(canary) / np.expm1(level)) if level > 0 else np.nan,
                     'null_reference': name, 'null_cells': int(ref.sum()),
                     'null_rate_at_level': float(pos[m][ref].mean())})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=20)
    a = ap.parse_args()
    run = a.run.resolve()
    tags = [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())]
    with ProcessPoolExecutor(min(a.workers, len(tags))) as ex:
        rows = [r for rs in ex.map(image_levels, [run] * len(tags), tags) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(run / 'typing' / 'canary_levels.csv', index=False)
    status = {r['image_name']: r['donor_status'] for r in RG.csv.DictReader(RG.COHORT.open())}
    nd = df[df['image'].map(status).eq('ND')]
    pd.set_option('display.width', 200)
    print(nd.groupby('marker').agg(null_reference=('null_reference', 'first'),
                                   null_rate_at_level=('null_rate_at_level', 'median'),
                                   canary_fold_over_level=('canary_fold_over_level', 'median')).round(3).to_string())
    return 0


if __name__ == '__main__':
    sys.exit(main())
