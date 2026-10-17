#!/usr/bin/env python3
"""endocrine_checks.py -- what the cells inside the endocrine objects are, before any conclusion is drawn.

Owner, 2026-10-16, on the Excluded share inside the objects: "you can't make this determination yet because you
need to do the checks for other endocrine and the object cells."

Two checks over the analysis-eligible cells whose centroid is inside an accepted endocrine object:

other endocrine  For each cell the classification did NOT call endocrine (Excluded, Immune, structure, Exocrine):
                 does it carry endocrine evidence other than a hormone clearing the image's level --
                   chga       CHGA at or above the image's ordinary-tissue CHGA level (Senior's CD99 reference level)
                   ghrelin    Ghrelin at or above its ordinary-tissue p99 (the four additions that carry it)
                   sub_level  the hormone statistic (classify.hormone_statistic) at or above the ordinary-tissue
                              p95 of that statistic but below the image's hormone level (its tissue p99): by
                              construction 4 % of tissue cells fall in this band, so `sub_level_expected` is
                              the count background alone would put there
                 (Before 2026-10-17 sub_level was "scaled > 0", which two-thirds of ordinary exocrine cells pass.)
object cells     For each object, how many of its cells are endocrine; objects whose cells hold no endocrine call,
                 by object class, and the share of such objects' cells by class.

Written to `$OUT/typing/endocrine_checks.json` and `endocrine_checks_objects.csv`.

    python scripts/uniform/typing/endocrine_checks.py --out data/development/uniform-20donor-20261015
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

ENDOCRINE = ('Beta', 'Alpha', 'Delta', 'Endocrine (unresolved subtype)', 'Endocrine (ambiguous subtype)')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    typ = CT.typing_dir(a.out)
    scale = CL.load_scale(a.out)
    report, objects = {}, []
    for r in CT.cohort_rows():
        t = r['image_name']
        calls = pd.read_parquet(typ / 'calls' / f'{t}.parquet')
        audit = json.loads((typ / 'calls' / f'{t}.json').read_text())
        al = pd.read_parquet(typ / 'aligned' / f'{t}.parquet')
        ctx = pd.read_parquet(typ / 'context' / f'{t}.parquet', columns=['object_id', 'context'])
        markers = [m for m in ('INS', 'GCG', 'SST', 'CHGA', 'Ghrelin') if m in al.columns]
        s = pd.DataFrame(CL.apply_scale(al[markers].to_numpy(), markers, scale), columns=markers)
        s['object_id'] = al['object_id']
        d = calls.merge(s, on='object_id', how='left')
        tissue = ctx.set_index('object_id').reindex(d['object_id'])['context'].eq('tissue').to_numpy()
        inside = d['islet_num'].notna().to_numpy()
        endo = d['cell_class'].isin(ENDOCRINE).to_numpy()
        chga = d['CHGA'].to_numpy() >= audit['levels']['CHGA']
        ghrelin = np.zeros(len(d), bool)
        if 'Ghrelin' in d:
            g = d['Ghrelin'].to_numpy()
            ghrelin = g >= np.nanpercentile(g[tissue], CT.REFERENCE_QUANTILE)
        fold = d['hormone_fold'].to_numpy(dtype=float)
        sub = (fold >= np.nanpercentile(fold[tissue], 95)) & (fold < audit['levels']['hormone'])
        rep = {'cells_in_objects': int(inside.sum()), 'endocrine': int((inside & endo).sum())}
        for cls in sorted(set(d.loc[inside & ~endo, 'cell_class'])):
            m = inside & (d['cell_class'] == cls).to_numpy()
            rep[cls] = {'cells': int(m.sum()), 'chga': int((m & chga).sum()), 'ghrelin': int((m & ghrelin).sum()),
                        'hormone_below_level': int((m & sub & ~chga).sum()),
                        'hormone_below_level_expected': round(0.04 * int(m.sum())),
                        'no_endocrine_evidence': int((m & ~chga & ~ghrelin & ~sub).sum())}
        per = d[inside].groupby('islet_num').agg(object_class=('endocrine_object_class', 'first'),
                                                  cells=('object_id', 'size'),
                                                  endocrine=('cell_class', lambda c: int(c.isin(ENDOCRINE).sum())),
                                                  excluded=('cell_class', lambda c: int((c == 'Excluded').sum())))
        per['image'] = t
        objects.append(per.reset_index())
        rep['objects_without_endocrine_cell'] = {
            k: f'{int(((per.object_class == k) & (per.endocrine == 0)).sum())} of {int((per.object_class == k).sum())}'
            for k in CT.ENDOCRINE_OBJECTS}
        report[t] = rep
        x = rep.get('Excluded', {'cells': 0})
        print(f"{t:10} in objects {rep['cells_in_objects']:6} endocrine {rep['endocrine']:6} | Excluded {x['cells']:5}"
              f" chga {x.get('chga', 0):5} ghrelin {x.get('ghrelin', 0):4} sub-level {x.get('hormone_below_level', 0):5}"
              f" (bg {x.get('hormone_below_level_expected', 0):4})"
              f" none {x.get('no_endocrine_evidence', 0):5} | objects w/o endocrine {rep['objects_without_endocrine_cell']}",
              flush=True)
    (typ / 'endocrine_checks.json').write_text(json.dumps(report, indent=1))
    pd.concat(objects).to_csv(typ / 'endocrine_checks_objects.csv', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
