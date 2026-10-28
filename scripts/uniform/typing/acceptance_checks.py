#!/usr/bin/env python3
"""acceptance_checks.py -- the final typing run's checks before assembly (final typing plan of 2026-10-28, Step 6).

Reads typing/identity/<image>.parquet (identity.py), typing/levels.csv, the registry, and for comparison the Barlow
cluster labels (output/barlow_typing_comparison_20261027/work/barlow_labels.parquet), root v16's cell types and v17's
calls. Writes to typing/checks/:
  per_label.csv      image x label: cells, confirmation rate, canary status mix, label-free co-expression (each
                     defining marker's positive rate in the label, e.g. CD45+ among T cells, CD3e+ among macrophages,
                     E-cadherin+ among hematopoietic cells, cells with two hormones positive among endocrine cells)
  confirmation.csv   image x lineage x confirming marker: its canary rate in the lineage against its rate in the other
                     labels (the plan: "each rate reported against the same marker's rate in unrelated classes")
  barlow.csv         image x label x Barlow label (cells), reported, never matched to
  transitions_v16.csv, transitions_v17.csv  label x earlier label (cells, all images)
  assertions.json    every assertion with the cells that break it (must all be 0)

    python scripts/uniform/typing/acceptance_checks.py --run $OUT --workers 20
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
import registry as RG  # noqa: E402
import identity as ID  # noqa: E402

BARLOW = REPO / 'output' / 'barlow_typing_comparison_20261027' / 'work' / 'barlow_labels.parquet'
V16 = REPO / 'data' / 'development' / 'uniform-20donor-20261027-v16' / 'cell_types'
V17 = 'superseded/typing-v17-20261028'
ALLOWED_MEASUREMENTS = set(RG.COMPARTMENT_COLUMN)


def named_markers(row_named_by: str) -> list[str]:
    """The markers a label was named by; none for the residual lineage ('compartment') and for the established
    endocrine typing ('established'), which the owner kept as it was (2026-10-28) and is reported, not asserted."""
    return [m for m in str(row_named_by).split(';') if m and m not in ('compartment', 'established', 'canary')]


def image_checks(run: Path, tag: str) -> dict:
    reg = RG.load()
    markers = ID.used_markers(run, reg)
    plan = ID.build_plan(reg, markers)
    d = pd.read_parquet(run / 'typing' / 'identity' / f'{tag}.parquet')
    P = {m: d[f'pos__{m}'].to_numpy(bool) for m in markers}
    C = {m: d[f'canary__{m}'].to_numpy(bool) for m in markers}
    lab = d['label'].to_numpy(dtype=object)
    starts = set(d['start_class'])
    out = {'image': tag, 'assert': {}}
    A = out['assert']

    # 1 no label from a marker below its level: a label named by markers has one of them positive
    nb = d['named_by'].to_numpy(dtype=object)
    bad = 0
    for i in np.flatnonzero(nb != ''):
        ms = named_markers(nb[i])
        if ms and not any(P[m][i] for m in ms if m in P):
            bad += 1
    A['label_from_marker_below_level'] = bad
    # residual (named by the compartment): one of the residual lineage's markers positive
    res = np.zeros(len(d), bool)
    for c, info in plan['comp'].items():
        if info['residual']:
            rows = (d['named_by'] == 'compartment').to_numpy() & (d['compartment'] == c).to_numpy()
            ok = np.zeros(len(d), bool)
            for m in info['residual_members']:
                ok |= P[m]
            res |= rows & ~ok
    A['residual_label_without_positive_marker'] = int(res.sum())
    # 2 lineage requirements: broad marker positive for lineages named by it; the registry's requires
    for c, info in plan['comp'].items():
        for lin, L_ in info['lineages'].items():
            rule = ID.label_row(plan, c, lin)
            if not rule or rule.get('named_by') != 'broad':
                continue
            labels = plan['labels']
            own = labels.loc[labels['lineage'].str.split(':').str[0] == lin, 'label'].tolist()
            rows = np.isin(lab, own)
            if rule.get('subtype_rule') == 'profile':
                rows |= np.char.find(lab.astype(str), f"{rule['label']}_") >= 0
            broad_pos = np.zeros(len(d), bool)
            for m in L_['broad']:
                broad_pos |= P[m]
            A[f'{lin}_label_without_broad_marker'] = int((rows & ~broad_pos).sum())
            for m in L_['broad']:
                for req in reg.requires(m):
                    if req in P:
                        A[f'{lin}_label_without_{req}'] = int((rows & ~P[req]).sum())
    # 3 no marker without a typing role names a cell; no marker outside the registry
    role = {m: reg.role(m) for m in reg.table.index}
    by_one = np.asarray([f'{x}+' == y for x, y in zip(nb, lab)])
    A['label_named_by_marker_without_typing_role'] = int(sum(role.get(x) not in RG.TYPING_ROLES
                                                             for x in nb[by_one]))
    A['marker_outside_registry'] = len(set(markers) - set(reg.typing_markers()))
    # 4 no default class: a cell keeps its start class only with a status saying why
    kept = np.isin(lab, list(starts))
    A['start_class_without_status'] = int((kept & (d['status'] == '').to_numpy()).sum())
    A['start_class_with_positive_marker_not_object_defined'] = int(
        (kept & (d['status'] == 'no_marker_above_background').to_numpy() & (d['n_positive'] > 0).to_numpy()).sum())

    # tables
    rows = []
    for L, g in d.groupby('label'):
        r = {'image': tag, 'label': L, 'cells': len(g),
             'confirmation_rate': float(g['confirmed'].mean()) if g['confirmed'].notna().any() else np.nan,
             'confirmable': int(g['confirmed'].notna().sum())}
        for s in ID.STATUSES:
            r[f'compartment_{s}'] = float((g['compartment_status'] == s).mean())
        idx = g.index.to_numpy()
        for m in markers:
            r[f'pos_{m}'] = float(P[m][idx].mean())
        hormones = list(reg.subtypes('Endocrine'))
        r['two_or_more_hormones'] = float(np.column_stack([P[h][idx] for h in hormones if h in P]).sum(axis=1)
                                          .__ge__(2).mean()) if hormones else np.nan
        rows.append(r)
    out['per_label'] = rows
    conf = []
    confirms = {m: RG.split(reg.table.loc[m, 'confirms']) for m in markers} if 'confirms' in reg.table.columns else {}
    for m, lins in confirms.items():
        for lin in lins:
            inside = (d['lineage'] == lin).to_numpy() & ~np.isin(lab, list(starts) + ['Excluded'])
            if not inside.any():
                continue
            conf.append({'image': tag, 'lineage': lin, 'marker': m, 'cells': int(inside.sum()),
                         'canary_rate_in_lineage': float(C[m][inside].mean()),
                         'canary_rate_elsewhere': float(C[m][~inside & ~np.isin(lab, list(starts))].mean())})
    out['confirmation'] = conf
    # comparisons (reported, never matched to)
    ids = d[['object_id', 'label']]
    if BARLOW.exists():
        b = pd.read_parquet(BARLOW, columns=['object_id', 'image', 'barlow_label'])
        b = b[b['image'] == tag]
        out['barlow'] = ids.merge(b, on='object_id', how='left').fillna({'barlow_label': 'absent'}) \
            .groupby(['label', 'barlow_label']).size().rename('cells').reset_index().assign(image=tag).to_dict('records')
    case = tag.split('_')[1].lstrip('0')
    f16 = V16 / f'{case}.parquet'
    if f16.exists():
        v = pd.read_parquet(f16, columns=['object_id', 'phenotype'])
        out['v16'] = ids.merge(v, on='object_id', how='left').fillna({'phenotype': 'Excluded or absent'}) \
            .groupby(['label', 'phenotype']).size().rename('cells').reset_index().to_dict('records')
    calls = run / V17 / 'calls' / f'{tag}.parquet'
    if calls.exists():
        sys.path.insert(0, str(HERE.parent))
        import assemble_uniform as AU  # the v17 labels as v17's own code maps them
        c17 = pd.read_parquet(calls)
        imm = pd.read_parquet(run / V17 / 'immune' / f'{tag}.parquet',
                              columns=['object_id', 'immune_class', 'immune_subtype', 'macrophage_profile'])
        v = pd.DataFrame({'object_id': c17['object_id'], 'v17': AU.app_labels(c17, imm)})
        out['v17'] = ids.merge(v, on='object_id', how='left').fillna({'v17': 'absent'}) \
            .groupby(['label', 'v17']).size().rename('cells').reset_index().to_dict('records')
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=20)
    a = ap.parse_args()
    run = a.run.resolve()
    reg = RG.load()
    tags = [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())]
    dest = run / 'typing' / 'checks'
    dest.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(min(a.workers, len(tags))) as ex:
        res = list(ex.map(image_checks, [run] * len(tags), tags))
    # run-level assertions: one measurement per marker, from an allowed compartment, the same on every image;
    # every typing marker acquired on every image
    lv = pd.read_csv(run / 'typing' / 'levels.csv')
    used = ID.used_markers(run, reg)
    meas = lv[lv['marker'].isin(used)].groupby('marker')['measurement'].nunique()
    run_assert = {'markers_with_two_measurements': int((meas > 1).sum()),
                  'measurement_outside_whole_cell_membrane_nucleus': int((~lv['measurement'].isin(
                      ALLOWED_MEASUREMENTS)).sum()),
                  'typing_marker_not_on_every_image': int(sum(reg.table.loc[m, 'on_every_image'] != 'yes'
                                                               for m in used)),
                  'typing_marker_without_level': int(lv[lv['marker'].isin(used)]['level_log'].isna().sum())}
    assertions = {'run': run_assert, 'images': {r['image']: r['assert'] for r in res}}
    (dest / 'assertions.json').write_text(json.dumps(assertions, indent=1))
    pd.DataFrame([x for r in res for x in r['per_label']]).to_csv(dest / 'per_label.csv', index=False)
    pd.DataFrame([x for r in res for x in r['confirmation']]).to_csv(dest / 'confirmation.csv', index=False)
    if any('barlow' in r for r in res):
        pd.DataFrame([x for r in res for x in r.get('barlow', [])]).to_csv(dest / 'barlow.csv', index=False)
    for k in ('v16', 'v17'):
        t = pd.DataFrame([x for r in res for x in r.get(k, [])])
        if len(t):
            t.groupby([c for c in t.columns if c != 'cells'])['cells'].sum().reset_index() \
                .to_csv(dest / f'transitions_{k}.csv', index=False)
    broken = {k: v for k, v in run_assert.items() if v} | {
        f'{img}:{k}': v for img, a in assertions['images'].items() for k, v in a.items() if v}
    print('assertions broken:', broken if broken else 'none', flush=True)
    return 1 if broken else 0


if __name__ == '__main__':
    sys.exit(main())
