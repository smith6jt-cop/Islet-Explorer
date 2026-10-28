#!/usr/bin/env python3
"""identity.py -- canary evidence and identity for every cell (final typing plan of 2026-10-28: method Steps 3-4, order of
work Step 6).

Inputs, per image: the registry and its labels file (`registry/<panel>__<tissue>.csv`, `..__labels.csv`: marker roles,
compartments, lineages, label rules), every typing marker's positivity level (`typing/levels.csv`, levels.py) and canary
level (`typing/canary_levels.csv`, canary.py), the REDSEA value of each marker's one measurement, and the cell's context.
This module names compartments and lineages, never markers.

1 Start class: the parent object's class, Islet, EndocrineCluster, EndocrineSingle or Tissue (owner 2026-10-28).
2 Canary (Islet-Explorer-Senior Xenium notebook 02 cell 22, adapted as planned): a marker is canary-positive above its
  canary level. Per compartment, depth = the canary-positive markers among its defining markers and the shared markers it
  carries, counted only when a defining marker is canary-positive; groups detected = compartments with a canary-positive
  defining marker; status none / single / resolved (2 groups, top depth >= CANARY_MIN_DEPTH and >= CANARY_MARGIN above
  the runner-up) / lineage_ambiguous (2 otherwise) / doublet_suspected (>= 3). Groups rank by depth x 1000 + the summed
  log fold over level of their canary-positive markers; the margin is taken between compartments (plan: fixes Xenium's
  margin between lineages of one group). The same rule runs between the lineages of the cell's compartment: lineage-
  specific markers detect a lineage; markers on several lineages, the compartment's pan marker and its shared markers add
  depth to every detected lineage. A cell with no canary-positive defining marker is placed the same way by its positive
  markers (plain level) with status none: agent determination, recorded in the ledger (the plan does not cover a
  defining marker positive below its canary level; the owner's rule admits names from positive markers; Xenium types
  cells before its canary step and uses the canary only for conflicts).
3 Identity by the labels file: Endocrine is the established endocrine typing (labels file named_by 'established'; owner
  2026-10-28: "The endocrine typing didn't need to change", "If the singles were already found and verified, don't undo
  that"): its decisions inside the islets, clusters and singles (Beta / Alpha / Delta, CHGA+ hormone-negative, unresolved
  subtype, and its exclusions hormone_unresolved_margin, exocrine_inside_islet_border, endocrine_single_*) are kept as
  they were, and no other cell is named endocrine: a cell positive for an endocrine marker is typed by its other evidence
  and tagged endocrine_marker_positive (endocrine_outside_object outside the objects). Exocrine: an epithelial cell
  without another lineage; inside an endocrine object Excluded (exocrine_inside_islet_border). A lineage named by its broad marker (T, B, Myeloid)
  needs that marker positive; T subtypes by positivity, both positive -> the higher fold over level (owner), neither ->
  the lineage label with t_subtype_status T_subtype_not_shown (interim, owner decision pending); B needs the registry's
  `requires` (CD45 positive, one hematopoietic lineage); Myeloid takes the profile label. A leukocyte without a named
  lineage is the pan label (CD45+) when the pan marker is positive. Every other placed cell, and every cell whose only
  positive markers are shared, is named by its top positive marker ('<marker>+'). A cell with no typing marker above its
  level keeps its start class (status no_marker_above_background). No marker below its level names a cell; markers
  without a typing role are never read.
4 Confirmed: compartment and lineage single or resolved, and a marker that confirms the label's lineage (registry
  `confirms`) and was not used to name it is canary-positive; NaN where the lineage has no such marker.

Writes typing/identity/<image>.parquet (one row per analysis-eligible cell) and <image>.json (counts).

    python scripts/uniform/typing/identity.py --run $OUT --workers 20 [--images ND_0112 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import registry as RG  # noqa: E402

STATUSES = ('none', 'single', 'resolved', 'lineage_ambiguous', 'doublet_suspected')
# the established endocrine typing: the run's typing state before the final run (v17's calls, never assembled; its
# endocrine chain unchanged since v15, decisions 98-99.5 % identical with v16's), and the exclusions it made
ESTABLISHED = 'superseded/typing-v17-20261028'
ESTABLISHED_REASONS = ('hormone_unresolved_margin', 'exocrine_inside_islet_border')
ESTABLISHED_REASON_PREFIX = 'endocrine_single_'
PLACED = ('single', 'resolved')


# -- the plan of one run: compartments, lineages and label rules from the registry ------------------------------------
def labels_table(reg: RG.Registry) -> pd.DataFrame:
    return pd.read_csv(RG.REGISTRY_DIR / f'{reg.name}__labels.csv', keep_default_na=False, dtype=str)


def used_markers(run: Path, reg: RG.Registry) -> list[str]:
    """Typing markers the level selection kept (decision 'use'), in registry order."""
    sel = pd.read_csv(run / 'typing' / 'levels_selection.csv').set_index('marker')
    return [m for m in reg.typing_markers() if m in sel.index and sel.loc[m, 'decision'] == 'use']


def build_plan(reg: RG.Registry, markers: list[str]) -> dict:
    lab = labels_table(reg)
    ix = {m: i for i, m in enumerate(markers)}
    comps = [c for c in reg.compartments() if any(reg.table.loc[m, 'compartment'] == c for m in markers
                                                  if reg.role(m) == 'defining')]
    plan = {'markers': markers, 'ix': ix, 'compartments': comps, 'labels': lab, 'comp': {}}
    for c in comps:
        defining = [m for m in markers if reg.role(m) == 'defining' and reg.table.loc[m, 'compartment'] == c]
        carried = [m for m in markers if c in RG.split(reg.table.loc[m, 'shared_by']) and m not in defining]
        rows = lab[lab['compartment'] == c]
        pan_lin = set(rows.loc[rows['kind'] == 'pan', 'lineage'])
        pan = [m for m in defining if {lin for lin, _ in reg.lineage_entries(m)} <= pan_lin and reg.lineage_entries(m)]
        lineages = {}
        for m in defining:
            if m in pan:
                continue
            for lin, sub in reg.lineage_entries(m):
                L = lineages.setdefault(lin, {'specific': [], 'members': [], 'broad': [], 'subtypes': {}})
                L['members'].append(m)
                if {x for x, _ in reg.lineage_entries(m)} == {lin}:
                    L['specific'].append(m)
                if reg.lineage_entries(m) == [(lin, '')]:
                    L['broad'].append(m)
                if sub:
                    L['subtypes'][sub] = m
        lineages = {k: v for k, v in lineages.items() if v['specific']}
        residual = next(iter(rows.loc[rows['kind'] == 'residual', 'lineage']), '')
        plan['comp'][c] = {'defining': defining, 'carried': carried, 'pan': pan, 'lineages': lineages,
                           'pan_label': next(iter(rows.loc[rows['kind'] == 'pan', 'label']), ''),
                           'residual': residual,
                           'residual_members': [m for m in defining
                                                if residual in {lin for lin, _ in reg.lineage_entries(m)}]}
    return plan


def label_row(plan: dict, compartment: str, lineage: str) -> dict:
    lab = plan['labels']
    r = lab[(lab['compartment'] == compartment) & (lab['lineage'] == lineage)]
    return r.iloc[0].to_dict() if len(r) else {}


# -- canary status (Xenium cell 22) --------------------------------------------------------------------------------------
def group_status(detected: np.ndarray, depth: np.ndarray, score: np.ndarray, margin: float, min_depth: float):
    """Per row: status, top group, runner-up group, their depths and the number of groups detected."""
    n = len(detected)
    rows = np.arange(n)
    ng = detected.sum(axis=1)
    comp = np.where(detected, depth * 1000.0 + score, -np.inf)
    top = comp.argmax(axis=1)
    comp2 = comp.copy()
    comp2[rows, top] = -np.inf
    runner = comp2.argmax(axis=1)
    td = np.where(ng > 0, depth[rows, top], 0)
    rd = np.where(ng > 1, depth[rows, runner], 0)
    status = np.select([ng == 0, ng == 1, ng >= 3, (td - rd >= margin) & (td >= min_depth)],
                       ['none', 'single', 'doublet_suspected', 'resolved'], 'lineage_ambiguous').astype(object)
    return status, top, np.where(ng > 1, runner, -1), td, rd, ng


def groups_matrix(X: np.ndarray, F: np.ndarray, detect: list[list[int]], count: list[list[int]]):
    """detected / depth / score for groups given as column-index lists (detect: markers that detect the group, count:
    markers that add depth once it is detected)."""
    n = len(X)
    det = np.column_stack([X[:, d].any(axis=1) if d else np.zeros(n, bool) for d in detect])
    dep = np.column_stack([X[:, c].sum(axis=1) for c in count]) * det
    sc = np.column_stack([np.where(X[:, c], F[:, c], 0.0).sum(axis=1) for c in count]) * det
    return det, dep, sc


def placed(P, C, F, detect, count, margin, min_depth):
    """Canary status, and the group a cell is placed in: by the canary level, else (status none) by plain positivity."""
    det, dep, sc = groups_matrix(C, F, detect, count)
    status, top, runner, td, rd, ng = group_status(det, dep, sc, margin, min_depth)
    detp, depp, scp = groups_matrix(P, F, detect, count)
    statusp, topp, _, _, _, ngp = group_status(detp, depp, scp, margin, min_depth)
    group = np.where(ng > 0, top, np.where(ngp > 0, topp, -1))
    basis = np.where(ng > 0, 'canary', np.where(ngp > 0, 'plain', '')).astype(object)
    # a cell placed below the canary level keeps its plain-level conflict visible (review 2026-10-28)
    status = np.where(ng > 0, status, np.where(ngp > 1, 'below_canary_' + statusp.astype(str), 'none')).astype(object)
    return {'group': group, 'basis': basis, 'status': status, 'n': ng, 'n_plain': ngp, 'depth': td,
            'runner': runner, 'runner_depth': rd, 'n_detected_at_basis': np.where(ng > 0, ng, ngp)}


# -- one image -------------------------------------------------------------------------------------------------------------
def load_values(run: Path, tag: str, reg: RG.Registry, markers: list[str], lv: pd.DataFrame) -> pd.DataFrame:
    cols = {m: reg.column(m, lv.loc[m, 'measurement']) for m in markers}
    rs = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet',
                         columns=['object_id', 'qc_analysis_eligible'] + sorted(set(cols.values())))
    ctx = pd.read_parquet(run / 'typing' / 'context' / f'{tag}.parquet',
                          columns=['object_id', 'x_um', 'y_um', 'context', 'endocrine_object_class'])
    d = rs[rs['qc_analysis_eligible']].merge(ctx, on='object_id', how='left', validate='one_to_one')
    d.attrs['columns'] = cols
    return d.reset_index(drop=True)


def established(path: Path, tag: str, object_id: pd.Series, reg: RG.Registry, plan: dict, row) -> tuple:
    """The established typing's decisions in a lineage (row of the labels file), mapped to this run's labels: a subtype
    class (e.g. 'Beta') to the subtype's label, the lineage label as it is, its unresolved or ambiguous subtype classes
    ('<Lineage> (...)') to the lineage name, and its exclusions (reasons ESTABLISHED_REASONS) kept as Excluded."""
    calls = pd.read_parquet(path / 'calls' / f'{tag}.parquet', columns=['object_id', 'cell_class', 'excluded_reason'])
    c = pd.DataFrame({'object_id': object_id}).merge(calls, on='object_id', how='left', validate='one_to_one')
    cls, why = c['cell_class'].fillna('').to_numpy(dtype=object), c['excluded_reason'].fillna('').to_numpy(dtype=object)
    lin = row.lineage
    out = np.full(len(c), '', dtype=object)
    for sub in set(reg.subtypes(lin).values()):
        out[cls == sub] = label_row(plan, row.compartment, f'{lin}:{sub}')['label']
    out[cls == row.label] = row.label
    out[np.char.startswith(cls.astype(str), f'{lin} (')] = lin
    excl = (cls == 'Excluded') & np.isin(why, ESTABLISHED_REASONS) | \
        (cls == 'Excluded') & np.char.startswith(why.astype(str), ESTABLISHED_REASON_PREFIX)
    out[excl] = 'Excluded'
    keep = out != ''
    return keep, out, np.where(excl, why, '')


def classify_image(run: Path, tag: str, endocrine_from: str = ESTABLISHED) -> tuple[pd.DataFrame, dict]:
    reg = RG.load()
    markers = used_markers(run, reg)
    state = reg.markers('state')
    lv_all = pd.read_csv(run / 'typing' / 'levels.csv')
    lv = lv_all[lv_all['image'] == tag].set_index('marker')
    cl = pd.read_csv(run / 'typing' / 'canary_levels.csv')
    cl = cl[cl['image'] == tag].set_index('marker')
    missing = [m for m in markers + state if not np.isfinite(lv.loc[m, 'level_log'])]
    missing += [m for m in markers if m not in cl.index or not np.isfinite(cl.loc[m, 'canary_log'])]
    if missing:
        raise SystemExit(f'{tag}: no level or canary level for {sorted(set(missing))}: the method must be fixed first')
    plan = build_plan(reg, markers)
    ix = plan['ix']
    d = load_values(run, tag, reg, markers + state, lv)
    n = len(d)
    V = np.column_stack([np.log1p(np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None))
                         for m in markers])
    L = lv.loc[markers, 'level_log'].to_numpy(float)
    Lc = np.maximum(cl.loc[markers, 'canary_log'].to_numpy(float), L)
    P, C, F = V > L, V > Lc, V - L
    start = d['endocrine_object_class'].fillna('Tissue').to_numpy(dtype=object)
    margin, min_depth = reg.const('CANARY_MARGIN'), reg.const('CANARY_MIN_DEPTH')

    # compartments
    comps = plan['compartments']
    detect = [[ix[m] for m in plan['comp'][c]['defining']] for c in comps]
    count = [[ix[m] for m in plan['comp'][c]['defining'] + plan['comp'][c]['carried']] for c in comps]
    cp = placed(P, C, F, detect, count, margin, min_depth)
    compartment = np.where(cp['group'] >= 0, np.asarray(comps, dtype=object)[np.maximum(cp['group'], 0)], '')
    comp_status = cp['status'].astype(object)

    # lineages inside the cell's compartment
    lineage = np.full(n, '', dtype=object)
    lin_basis = np.full(n, '', dtype=object)
    lin_status = np.full(n, '', dtype=object)
    lin_detected = np.zeros(n, int)
    for c in comps:
        info = plan['comp'][c]
        lins = list(info['lineages'])
        rows = np.flatnonzero(compartment == c)
        if not lins or not len(rows):
            continue
        extra = info['pan'] + info['carried']
        det_cols = [[ix[m] for m in info['lineages'][lin]['specific']] for lin in lins]
        cnt_cols = [[ix[m] for m in dict.fromkeys(info['lineages'][lin]['members'] + extra)] for lin in lins]
        lp = placed(P[rows], C[rows], F[rows], det_cols, cnt_cols, margin, min_depth)
        has = lp['group'] >= 0
        lineage[rows[has]] = np.asarray(lins, dtype=object)[lp['group'][has]]
        lin_basis[rows] = lp['basis']
        lin_status[rows] = np.where(has, lp['status'], '')
        lin_detected[rows] = lp['n_detected_at_basis']

    # identity
    label = start.copy()
    reason = np.full(n, '', dtype=object)
    named_by = np.full(n, '', dtype=object)
    tags = [[] for _ in range(n)]
    t_status = np.full(n, '', dtype=object)
    status = np.full(n, '', dtype=object)
    done = np.zeros(n, bool)
    names = np.asarray(markers, dtype=object)
    comp_basis = cp['basis'].copy()

    # the established typing of a lineage (labels file named_by 'established'; owner 2026-10-28: "The endocrine typing
    # didn't need to change", "If the singles were already found and verified, don't undo that"): its decisions are
    # kept as they were, and this hierarchy names no other cell in that lineage
    lab = plan['labels']
    for r in lab[lab['named_by'] == 'established'].itertuples():
        keep, est_label, est_reason = established(run / endocrine_from, tag, d['object_id'], reg, plan, r)
        for i in np.flatnonzero(keep & (compartment != r.compartment) & (compartment != '')):
            tags[i].append(f'canary_top_{compartment[i]}')
        label[keep], reason[keep] = est_label[keep], est_reason[keep]
        named_by[keep] = 'established'
        compartment[keep], lineage[keep] = r.compartment, r.lineage
        comp_basis[keep], lin_basis[keep], lin_status[keep] = 'established', 'established', ''
        done[keep] = True

    def top_positive(rows: np.ndarray, cols: list[int]) -> np.ndarray:
        """The positive marker with the highest fold over its level, among `cols`, for each row ('' if none)."""
        if not cols or not len(rows):
            return np.full(len(rows), '', dtype=object)
        sub = np.where(P[np.ix_(rows, cols)], F[np.ix_(rows, cols)], -np.inf)
        best = sub.argmax(axis=1)
        return np.where(np.isfinite(sub.max(axis=1)), names[np.asarray(cols)[best]], '')

    def by_marker(rows: np.ndarray, cols: list[int], why: str) -> None:
        tm = top_positive(rows, cols)
        ok = tm != ''
        label[rows[ok]] = np.asarray([f'{m}+' for m in tm[ok]], dtype=object)
        named_by[rows[ok]] = tm[ok]
        done[rows[ok]] = True
        for i in rows[ok]:
            tags[i].append(why) if why else None

    for c in comps:
        info = plan['comp'][c]
        in_c = (compartment == c) & ~done
        defining_cols = [ix[m] for m in info['defining']]
        for lin, L_ in info['lineages'].items():
            rule = label_row(plan, c, lin)
            if not rule or rule.get('kind') == 'residual':       # the residual lineage is labeled below
                continue
            rows = np.flatnonzero(in_c & (lineage == lin))
            if not len(rows):
                continue
            only_in = set(RG.split(rule.get('only_in', '')))
            # a lineage named only by its established typing, or outside its objects, names no cell here: the cell is
            # typed by its other evidence of the compartment (the residual lineage's markers), else named by its top
            # positive marker; tagged either way
            barred = np.ones(len(rows), bool) if rule.get('named_by') == 'established' else \
                (~np.isin(start[rows], list(only_in)) if only_in else np.zeros(len(rows), bool))
            out_rows = rows[barred]
            for i in out_rows:
                tags[i].append(f'{lin.lower()}_marker_positive' if start[i] in only_in else f'{lin.lower()}_outside_object')
            members = [ix[m] for m in info['residual_members']]
            other = P[np.ix_(out_rows, members)].any(axis=1) if members else np.zeros(len(out_rows), bool)
            lineage[out_rows[other]] = info['residual']
            # the lineage's own markers name nothing here (tiles 2026-10-28: 2-6 of 40 INS+/GCG+/SST+ cells showed the
            # hormone on the cell); another positive defining marker names the cell, else it keeps its start class
            rest = out_rows[~other]
            naming = [ix[m] for m in markers if reg.role(m) == 'defining' and m not in L_['specific']]
            by_marker(rest, naming, '')
            kept = rest[~done[rest]]
            done[kept] = True
            status[kept] = f'{lin.lower()}_marker_only'
            rows = rows[~barred]
            if not len(rows):
                continue
            broad = [ix[m] for m in L_['broad']]
            subs = {s: ix[m] for s, m in L_['subtypes'].items()}
            if rule.get('named_by') == 'broad':
                ok = P[np.ix_(rows, broad)].any(axis=1) if broad else np.zeros(len(rows), bool)
            else:                                   # broad_or_subtype
                ok = P[np.ix_(rows, broad + list(subs.values()))].any(axis=1)
            for m in L_['broad']:
                for req in reg.requires(m):
                    if req == 'single_lineage':
                        ok &= lin_detected[rows] == 1
                    elif req in ix:
                        ok &= P[rows, ix[req]]
            for i in rows[~ok]:
                tags[i].append(f'{lin}_not_named')
            rows = rows[ok]
            srule = rule.get('subtype_rule', '')
            if srule == 'exclusive_margin' and subs:
                cols = list(subs.values())
                Fs = F[np.ix_(rows, cols)]
                any_pos = P[np.ix_(rows, cols)].any(axis=1)
                order = np.sort(Fs, axis=1)
                clear = order[:, -1] - order[:, -2] >= np.log(reg.const('HORMONE_MARGIN')) if len(cols) > 1 else \
                    np.ones(len(rows), bool)
                win = np.asarray(list(subs), dtype=object)[Fs.argmax(axis=1)]
                named = any_pos & clear
                for s in subs:
                    sel = rows[named & (win == s)]
                    label[sel] = label_row(plan, c, f'{lin}:{s}')['label']
                    named_by[sel] = markers[subs[s]]
                unresolved = rows[any_pos & ~clear]
                label[unresolved] = 'Excluded'
                reason[unresolved] = 'hormone_unresolved_margin'
                neg = rows[~any_pos]
                label[neg] = rule['label']
                named_by[neg] = ';'.join(L_['broad'])
            elif srule == 'higher' and subs:
                cols = list(subs.values())
                Fs = np.where(P[np.ix_(rows, cols)], F[np.ix_(rows, cols)], -np.inf)
                any_pos = np.isfinite(Fs.max(axis=1))
                win = np.asarray(list(subs), dtype=object)[Fs.argmax(axis=1)]
                for s in subs:
                    sel = rows[any_pos & (win == s)]
                    label[sel] = label_row(plan, c, f'{lin}:{s}')['label']
                    named_by[sel] = markers[subs[s]]
                none_ = rows[~any_pos]
                label[none_] = rule['label']
                named_by[none_] = ';'.join(L_['broad'])
                t_status[none_] = f'{lin}_subtype_not_shown'
                both = rows[P[np.ix_(rows, cols)].sum(axis=1) > 1]
                t_status[both] = 'both_subtype_markers_higher_named'
            elif srule == 'profile':
                head, tail = (RG.split(x.replace(',', ';')) for x in rule['profile'].split('|'))
                sign = lambda m: np.where(P[rows, ix[m]], f'{m}+', f'{m}-') if m in ix else np.full(len(rows), f'{m}?')
                parts = [sign(m) for m in head] + [np.full(len(rows), f"{rule['label']}_")] + [sign(m) for m in tail]
                label[rows] = np.asarray([''.join(p) for p in zip(*parts)], dtype=object)
                named_by[rows] = ';'.join(head)
            else:
                label[rows] = rule['label']
                named_by[rows] = ';'.join(L_['broad'])
            done[rows] = True
        # the compartment's cells without a named lineage
        rows = np.flatnonzero((compartment == c) & ~done)
        if not len(rows):
            continue
        if info['residual']:
            rule = label_row(plan, c, info['residual'])
            excl = np.isin(start[rows], RG.split(rule.get('excluded_in', '')))
            label[rows[~excl]] = rule['label']
            named_by[rows[~excl]] = 'compartment'
            label[rows[excl]] = 'Excluded'
            reason[rows[excl]] = f"{info['residual'].lower()}_inside_islet_border"
            done[rows] = True
            continue
        if info['pan']:
            pan_pos = P[np.ix_(rows, [ix[m] for m in info['pan']])].any(axis=1)
            label[rows[pan_pos]] = info['pan_label']
            named_by[rows[pan_pos]] = ';'.join(info['pan'])
            done[rows[pan_pos]] = True
            rows = rows[~pan_pos]
        by_marker(rows, defining_cols, '')
    # no compartment: a shared marker positive names the cell; else the start class stays
    rest = np.flatnonzero(~done & (compartment == ''))
    by_marker(rest, [ix[m] for m in markers if reg.role(m) == 'shared'], 'shared_marker_only')
    left = ~done
    status[left & ~P.any(axis=1)] = 'no_marker_above_background'
    if (left & P.any(axis=1)).any():
        raise ValueError(f'{tag}: {int((left & P.any(axis=1)).sum())} cells with a positive marker were not named')

    # confirmation
    confirms = {m: RG.split(reg.table.loc[m, 'confirms']) for m in markers} if 'confirms' in reg.table.columns else {}
    confirmed = np.full(n, np.nan)
    confirming = np.full(n, '', dtype=object)
    lab_lineage = lineage.copy()
    for c in comps:
        info = plan['comp'][c]
        if info['residual']:
            lab_lineage[(compartment == c) & (label == label_row(plan, c, info['residual'])['label'])] = info['residual']
    placed_ok = np.isin(comp_status, PLACED) & (np.isin(lin_status, PLACED) | (lin_status == ''))
    for lin in sorted(set(lab_lineage) - {''}):
        cand = [m for m, ls in confirms.items() if lin in ls]
        if not cand:
            continue
        rows = np.flatnonzero(lab_lineage == lin)
        for i in rows:
            used = set(str(named_by[i]).split(';'))
            use = [ix[m] for m in cand if m not in used]
            if not use:
                continue
            hit = C[i, use]
            confirmed[i] = float(placed_ok[i] and hit.any())
            if hit.any():
                confirming[i] = ';'.join(names[np.asarray(use)[hit]])

    out = pd.DataFrame({
        'object_id': d['object_id'].to_numpy(), 'start_class': start, 'label': label, 'excluded_reason': reason,
        'named_by': named_by, 'compartment': compartment, 'compartment_basis': comp_basis,
        'compartment_status': comp_status, 'n_compartments_canary': cp['n'], 'n_compartments_plain': cp['n_plain'],
        'compartment_depth': cp['depth'], 'runner_up_compartment': np.where(cp['runner'] >= 0, np.asarray(
            comps, dtype=object)[np.maximum(cp['runner'], 0)], ''), 'runner_up_depth': cp['runner_depth'],
        'lineage': lineage, 'lineage_basis': lin_basis, 'lineage_status': lin_status,
        'subtype_status': t_status, 'status': status, 'tags': [';'.join(t) for t in tags],
        'confirmed': confirmed, 'confirming_markers': confirming, 'n_positive': P.sum(axis=1).astype(np.int16),
        'n_canary': C.sum(axis=1).astype(np.int16)})
    for m in state:
        x = np.log1p(np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None))
        out[f'{m}_pos'] = (x > lv.loc[m, 'level_log']).astype(np.int8)
    for j, m in enumerate(markers):
        out[f'pos__{m}'] = P[:, j].astype(np.int8)
        out[f'canary__{m}'] = C[:, j].astype(np.int8)
        out[f'fold__{m}'] = F[:, j].astype(np.float32)
    audit = {'image': tag, 'cells': n, 'markers': markers,
             'label': {str(k): int(v) for k, v in out['label'].value_counts().items()},
             'compartment_status': {str(k): int(v) for k, v in out['compartment_status'].value_counts().items()},
             'lineage_status': {str(k): int(v) for k, v in out['lineage_status'].value_counts().items()},
             'status': {str(k): int(v) for k, v in pd.Series(status[status != '']).value_counts().items()},
             'no_marker_above_background': int((status == 'no_marker_above_background').sum()),
             'excluded': {str(k): int(v) for k, v in out.loc[out['label'] == 'Excluded', 'excluded_reason']
                          .value_counts().items()}}
    return out, audit


def process(run: Path, tag: str, dest: Path, endocrine_from: str = ESTABLISHED) -> dict:
    out, audit = classify_image(run, tag, endocrine_from)
    out.to_parquet(dest / f'{tag}.parquet', index=False)
    (dest / f'{tag}.json').write_text(json.dumps(audit, indent=1))
    print(f"{tag}: {len(out)} cells, no marker above background {audit['no_marker_above_background']}", flush=True)
    return audit


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--images', nargs='*', default=None, help='default: every image of cohort.csv')
    ap.add_argument('--workers', type=int, default=20)
    ap.add_argument('--dest', type=Path, default=None, help='default: $RUN/typing/identity (trials: a scratch folder)')
    ap.add_argument('--endocrine-from', default=ESTABLISHED, help='the established endocrine typing, under $RUN')
    a = ap.parse_args()
    run = a.run.resolve()
    tags = a.images or [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())]
    dest = (a.dest or run / 'typing' / 'identity').resolve()
    dest.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(min(a.workers, len(tags))) as ex:
        audits = list(ex.map(process, [run] * len(tags), tags, [dest] * len(tags), [a.endocrine_from] * len(tags)))
    (dest / 'summary.json').write_text(json.dumps(audits, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
