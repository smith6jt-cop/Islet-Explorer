#!/usr/bin/env python3
"""Accuracy checks for the per-image positivity levels (owner-approved plan, 2026-10-28). Measurements only.

Compares the run's current typing (v16) with v15's (`superseded/typing-v15-20261027/`), with positivity read the same
way for both: each marker above its image's level (`typing/positivity_levels.csv`).
  1. how many cells each step decides (rules, top marker, excluded), v15 and v16;
  2. label transitions v15 -> v16, overall and per image;
  3. independent consistency checks on the typed cells, per root:
     - T cells (rule) that are CD45+; macrophages (rule) that are CD45+;
     - immune-rule cells that are also E-cadherin+ or panCK+;
     - T cells that are CD20+ or CD68+, B cells that are CD3e+ or CD68+, macrophages that are CD3e+;
  4. reference cells at or above their level (about 1 % by construction), per marker.

    python scripts/research/typing_accuracy_checks_20261028.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'uniform' / 'typing'))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402

RUN = REPO / 'data' / 'development' / 'uniform-20donor-20261027'
OUT = REPO / 'output' / 'typing_accuracy_checks_20261028'
CHECK = ['CD45', 'CD3e', 'CD20', 'CD68', 'CD163', 'ECAD', 'panCK', 'CD8a', 'CD4']


def calls(dirname: Path) -> pd.DataFrame:
    fr = []
    for r in CT.cohort_rows():
        t = r['image_name']
        c = pd.read_parquet(dirname / 'calls' / f'{t}.parquet',
                            columns=['object_id', 'cell_class', 'immune_definition', 'resolved_by', 'excluded_reason'])
        im = pd.read_parquet(dirname / 'immune' / f'{t}.parquet', columns=['object_id', 'immune_class', 'immune_subtype'])
        c = c.merge(im, on='object_id', how='left')
        c['image'] = t
        fr.append(c)
    d = pd.concat(fr, ignore_index=True)
    use = d['immune_class'].fillna('').where(d['immune_class'].fillna('') != '', d['immune_definition'])
    lab = d['cell_class'].astype(str)
    imm = d['cell_class'] == 'Immune'
    lab = lab.where(~imm, use.map({'Imm-T-cell': 'T cell', 'Imm-B-cell': 'B cell', 'Imm-macrophage': 'Macrophage',
                                   'Imm-CD45-leukocyte': 'CD45+ (rule)'}).fillna('Immune (no rule)'))
    d['label'] = lab
    d['step'] = np.select([d['resolved_by'].str.endswith('top_marker'), d['cell_class'] == 'Excluded'],
                          ['top marker', 'excluded'], 'rule')
    return d


def positivity() -> pd.DataFrame:
    typ = CT.typing_dir(RUN)
    scale = CL.load_scale(RUN)
    fr = []
    for r in CT.cohort_rows():
        t = r['image_name']
        al = pd.read_parquet(typ / 'aligned' / f'{t}.parquet', columns=['object_id'] + CHECK)
        e = CL.excess(al, CHECK, CL.load_levels(RUN, t)) > 0
        e.insert(0, 'object_id', al['object_id'])
        fr.append(e)
    return pd.concat(fr, ignore_index=True)


def consistency(d: pd.DataFrame, pos: pd.DataFrame, name: str) -> list[dict]:
    m = d.merge(pos, on='object_id', how='left')
    out = []
    for label, marker, expect in [('T cell', 'CD45', '+'), ('Macrophage', 'CD45', '+'), ('B cell', 'CD45', '+'),
                                  ('T cell', 'CD20', '-'), ('T cell', 'CD68', '-'), ('B cell', 'CD3e', '-'),
                                  ('B cell', 'CD68', '-'), ('Macrophage', 'CD3e', '-'),
                                  ('T cell', 'ECAD', '-'), ('Macrophage', 'ECAD', '-'), ('B cell', 'ECAD', '-'),
                                  ('CD45+ (rule)', 'ECAD', '-'), ('T cell', 'panCK', '-')]:
        for step in ('rule', 'top marker'):
            g = m[(m['label'] == label) & (m['step'] == step)]
            if not len(g):
                continue
            share = float(g[marker].mean())
            out.append({'root': name, 'label': label, 'decided_by': step, 'cells': len(g), 'marker': marker,
                        'expected': expect, 'proportion_positive': round(share, 4)})
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    v15 = calls(RUN / 'superseded' / 'typing-v15-20261027')
    v16 = calls(CT.typing_dir(RUN))
    pos = positivity()
    steps = pd.concat([v15.assign(root='v15'), v16.assign(root='v16')]).groupby(['root', 'step']).size().unstack(0)
    steps.to_csv(OUT / 'cells_by_step.csv')
    labels = pd.concat([v15.assign(root='v15'), v16.assign(root='v16')]).groupby(['root', 'label']).size().unstack(0)
    labels = labels.fillna(0).astype(int).sort_values('v16', ascending=False)
    labels.to_csv(OUT / 'labels_v15_v16.csv')
    m = v15[['object_id', 'image', 'label']].merge(v16[['object_id', 'label']], on='object_id', suffixes=('_v15', '_v16'))
    moved = m[m['label_v15'] != m['label_v16']]
    moved.groupby(['label_v15', 'label_v16']).size().sort_values(ascending=False).rename('cells') \
        .to_csv(OUT / 'transitions_v15_v16.csv')
    moved.groupby('image').size().rename('cells_changed').to_csv(OUT / 'transitions_by_image.csv')
    cons = pd.DataFrame(consistency(v15, pos, 'v15') + consistency(v16, pos, 'v16'))
    cons.to_csv(OUT / 'consistency_checks.csv', index=False)
    lv = pd.read_csv(CT.typing_dir(RUN) / 'positivity_levels.csv')
    lv.to_csv(OUT / 'positivity_levels.csv', index=False)
    print(steps.to_string()); print(labels.to_string())
    print('changed', len(moved)); print(moved.groupby(['label_v15', 'label_v16']).size().sort_values().tail(20).to_string())
    print(cons.pivot_table(index=['label', 'marker', 'expected', 'decided_by'], columns='root',
                           values='proportion_positive').to_string())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
