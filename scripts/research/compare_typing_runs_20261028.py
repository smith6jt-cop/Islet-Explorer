#!/usr/bin/env python3
"""Cell-by-cell comparison of two typing outputs of run uniform-20donor-20261027 (owner's typing review, 2026-10-28).
Measurements only.

Each side is a folder holding `calls/<image>.parquet` and `immune/<image>.parquet` (the run's `typing/` or one of its
`superseded/typing-*` copies). A cell's label is its rule class (T cell, B cell, Macrophage, CD45+ (rule), the
endocrine and Exocrine classes, Excluded) or, when the top-marker step decided it, 'top:' + its class. Reported: the
label table of both sides, the cells that changed, the rule-typed immune cells of side A that changed, and the cells
side A left to the top-marker step that side B types by a rule.

    python scripts/research/compare_typing_runs_20261028.py --a superseded/typing-v16-20261028 --b typing
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
IMMUNE = {'Imm-T-cell': 'T cell', 'Imm-B-cell': 'B cell', 'Imm-macrophage': 'Macrophage',
          'Imm-CD45-leukocyte': 'CD45+ (rule)'}


def labels(folder: Path) -> pd.DataFrame:
    frames = []
    for r in CT.cohort_rows():
        t = r['image_name']
        c = pd.read_parquet(folder / 'calls' / f'{t}.parquet',
                            columns=['object_id', 'cell_class', 'immune_definition', 'resolved_by'])
        im = pd.read_parquet(folder / 'immune' / f'{t}.parquet', columns=['object_id', 'immune_class'])
        c = c.merge(im, on='object_id', how='left')
        use = c['immune_class'].fillna('').where(c['immune_class'].fillna('') != '', c['immune_definition'])
        lab = c['cell_class'].where(c['cell_class'] != 'Immune', use.map(IMMUNE).fillna('Immune (no rule)'))
        lab = lab.where(~c['resolved_by'].str.endswith('top_marker'), 'top:' + lab)
        frames.append(pd.DataFrame({'object_id': c['object_id'], 'image': t, 'label': lab}))
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a', required=True, help='folder under the run, e.g. superseded/typing-v16-20261028')
    ap.add_argument('--b', required=True, help='folder under the run, e.g. typing')
    ap.add_argument('--out', type=Path, help='write the tables here')
    a = ap.parse_args()
    A, B = labels(RUN / a.a), labels(RUN / a.b)
    m = A.merge(B[['object_id', 'label']], on='object_id', suffixes=('_a', '_b'))
    tab = pd.DataFrame({'a': A['label'].value_counts(), 'b': B['label'].value_counts()}).fillna(0).astype(int)
    tab['change'] = tab['b'] - tab['a']
    moved = m[m['label_a'] != m['label_b']]
    trans = moved.groupby(['label_a', 'label_b']).size().sort_values(ascending=False).rename('cells')
    rule_immune = set(IMMUNE.values())
    lost = moved[moved['label_a'].isin(rule_immune)]
    rescued = moved[moved['label_a'].str.startswith('top:') & ~moved['label_b'].str.startswith('top:')]
    pd.set_option('display.width', 200)
    print(f'A = {a.a}, B = {a.b}; cells compared {len(m):,}; changed {len(moved):,}')
    print(tab.sort_values('b', ascending=False).to_string())
    print('\nlargest transitions:\n' + trans.head(20).to_string())
    print(f'\nrule-typed immune cells of A that changed: {len(lost):,}')
    print(lost.groupby(['label_a', 'label_b']).size().sort_values(ascending=False).head(10).to_string())
    print(f'\ncells A left to the top-marker step that B types by a rule: {len(rescued):,}')
    print(rescued.groupby(['label_a', 'label_b']).size().sort_values(ascending=False).head(10).to_string())
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True)
        tab.to_csv(a.out / 'labels.csv')
        trans.to_csv(a.out / 'transitions.csv')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
