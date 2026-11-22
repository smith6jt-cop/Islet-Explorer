#!/usr/bin/env python3
"""adjust_t_levels.py -- normal-donor T-partition adjustment of a 'higher'-rule lineage's two subtype levels.

Owner, 2026-11-06: "Use the basis of no DN or DP to adjust the thresholds for CD8 and CD4 T Cells in the normal
donors." The literature check (ledger § 2026-11-06) found no support for double-negative or double-positive T-cell
populations in normal human pancreas, so on each normal-donor image the two subtype markers' levels are re-set so the
lineage's cells partition into single-positives, with the residual DN/DP reported, never forced.

Runs AFTER levels.py and BEFORE canary.py, area_coverage.py and identity.py, which all consume the adjusted
typing/levels.csv. Everything is registry/labels-driven (owner 2026-10-28: marker roles are data):

  1. The adjusted lineages: every labels row of kind `lineage` with subtype_rule `higher` and exactly two subtype
     markers (today: Hematopoietic T, CD4/CD8a).
  2. The universe U of an image: positive for each of the broad marker's registry `requires`, negative for every
     other defining marker of the compartment, at the plain levels. The broad marker itself does not gate U — its
     level is the third adjusted threshold (owner: "the thresholding adjustment is done three ways with CD3e").
  3. Three-way mutual fit, iterated from the current levels: T = U above the broad level; within T the higher fold
     over level wins the side (owner: "there are no double positives since the highest marker wins for those" — a
     double positive is never an error and the losing side's values never push a threshold). Each subtype level is
     the optimal separator between its side's T cells and the non-T cells of U (the other side excluded); the broad
     level is the separator between the subtype-positive cells of U (would-be orphans) and the subtype-negative
     ones (would-be DN). Degenerate states keep the last valid iterate.
  4. Floor: an adjusted level never falls below the (1 - tolerated out-of-place rate) quantile of the marker's
     known-negative reference (endocrine-object cells positive for an endocrine-specific marker, and epithelial
     cells without one), the registry's own CD45_AMONG_HORMONE_MAX tolerance. Without it a T-cell-only fit can land
     below tissue background (measured 2026-11-06: 94% of ND_6548 would turn positive).
  5. typing/levels.csv rows of the adjusted (image, marker) pairs are updated in place (level_log, level_raw); the
     originals, floors and fit audit go to typing/t_partition_adjustment.json. Other images are untouched.

The same fit on the current preview inputs is output/typing_revision_20261105/adjust_nd_t_levels.py (keep the two in
step); the preview applies the result through recomputed flags (preview_revision.apply_adjusted_levels).

    python scripts/uniform/typing/adjust_t_levels.py --run $OUT [--images ND_0112 ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import registry as RG
import identity as ID

NORMAL_PREFIX = 'ND_'          # the normal-donor images (owner 2026-11-06: adjust in the normal donors)
OBJECT_CLASSES = ('Islet', 'EndocrineCluster', 'EndocrineSingle')


def separator(x: np.ndarray, positive: np.ndarray) -> float:
    """Optimal 1-D threshold between assigned-positive and assigned-negative values (fewest misclassified; widest
    gap midpoint; lower threshold)."""
    v = np.unique(x)
    cand = np.concatenate([[v[0] - 1e-3], (v[:-1] + v[1:]) / 2, [v[-1] + 1e-3]])
    gaps = np.concatenate([[1.0], np.minimum(np.diff(v), 1.0), [1.0]])
    above = x[None, :] > cand[:, None]
    err = (above & ~positive[None, :]).sum(axis=1) + (~above & positive[None, :]).sum(axis=1)
    k = np.flatnonzero(err == err.min())
    k = k[np.lexsort((cand[k], -gaps[k]))][0]
    return float(cand[k])


def threeway_fit(xg: np.ndarray, xa: np.ndarray, xb: np.ndarray, lg: float, la: float, lb: float,
                 fg: float, fa: float, fb: float) -> dict:
    """Three-way mutual fit over U: T = above the gate (broad) level; sides within T by the higher fold (a double
    positive is not an error); each subtype level separates its side from the non-T cells of U; the gate level
    separates subtype-positive from subtype-negative cells of U. Floored; degenerate states keep the last iterate."""
    tg, ta, tb = lg, la, lb
    iters = 0
    collapsed = False
    for iters in range(1, 16):
        T = xg > tg
        a_side = T & ((xa - ta) > (xb - tb))
        b_side = T & ~a_side
        if not T.any() or not a_side.any() or not b_side.any():
            collapsed = True
            break
        sel_a = a_side | ~T
        na = max(separator(xa[sel_a], a_side[sel_a]), fa)
        sel_b = b_side | ~T
        nb = max(separator(xb[sel_b], b_side[sel_b]), fb)
        sub_pos = (xa > na) | (xb > nb)
        if not sub_pos.any() or sub_pos.all():
            collapsed = True
            break
        ng = max(separator(xg, sub_pos), fg)
        if (ng, na, nb) == (tg, ta, tb):
            break
        tg, ta, tb = ng, na, nb
    T = xg > tg
    a_side = T & ((xa - ta) > (xb - tb))
    dp = int((T & (xa > ta) & (xb > tb)).sum())
    dn = int((T & (xa <= ta) & (xb <= tb)).sum())
    orphan = int((~T & ((xa > ta) | (xb > tb))).sum())
    return {'tg': tg, 'ta': ta, 'tb': tb, 'dp': dp, 'dn': dn, 'orphan_subtype': orphan,
            'n_universe': int(len(xg)), 't_cells': int(T.sum()), 'iterations': iters,
            'collapsed': bool(collapsed), 'a_side': int(a_side.sum()), 'b_side': int((T & ~a_side).sum())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True, type=Path)
    ap.add_argument('--images', nargs='*', default=None)
    a = ap.parse_args()
    run = a.run
    reg = RG.load()
    markers = ID.used_markers(run, reg)
    plan = ID.build_plan(reg, markers)
    lv_all = pd.read_csv(run / 'typing' / 'levels.csv')
    lab = plan['labels']
    targets = []                     # (compartment, lineage, [two subtype markers], gate marker, requires, negatives)
    for r in lab[(lab['kind'] == 'lineage') & (lab['subtype_rule'] == 'higher')].itertuples():
        L_ = plan['comp'].get(r.compartment, {}).get('lineages', {}).get(r.lineage)
        if not L_ or len(L_['subtypes']) != 2 or len(L_['broad']) != 1:
            continue
        req = sorted({q for b in L_['broad'] for q in reg.requires(b) if q in markers})
        negatives = [m for m in plan['comp'][r.compartment]['defining'] if m not in L_['members']]
        targets.append((r.compartment, r.lineage, sorted(L_['subtypes'].values()), L_['broad'][0], req, negatives))
    if not targets:
        raise SystemExit('no higher-rule lineage with two subtype markers in the labels table')
    tol = float(reg.const('CD45_AMONG_HORMONE_MAX'))
    # the floor reference: endocrine-object cells positive for an endocrine-specific marker, and epithelial cells
    # without one (the registry reference recorded for these markers' known-negative levels)
    epi_def = plan['comp'].get('Epithelial', {}).get('defining', [])
    endo_specific = plan['comp'].get('Epithelial', {}).get('lineages', {}).get('Endocrine', {}).get('specific', [])

    tags = a.images or sorted(lv_all['image'].unique())
    tags = [t for t in tags if t.startswith(NORMAL_PREFIX)]
    audit = {'tolerance_out_of_place': tol, 'images': {}}
    lv_new = lv_all.copy()
    for tag in tags:
        lv = lv_all[lv_all['image'] == tag].set_index('marker')
        d = ID.load_values(run, tag, reg, markers, lv)
        V = {m: np.log1p(np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None))
             for m in markers}
        P = {m: V[m] > float(lv.loc[m, 'level_log']) for m in markers}
        in_obj = d['endocrine_object_class'].fillna('Tissue').isin(OBJECT_CLASSES).to_numpy()
        endo_pos = np.any([P[m] for m in endo_specific], axis=0) if endo_specific else np.zeros(len(d), bool)
        epi_pos = np.any([P[m] for m in epi_def], axis=0) if epi_def else np.zeros(len(d), bool)
        ref = (in_obj & endo_pos) | (epi_pos & ~endo_pos)
        for comp, lin, subs, gate, req, negatives in targets:
            ma, mb = subs
            usel = np.ones(len(d), bool)
            for q in req:
                usel &= P[q]
            for m in negatives:
                if m not in subs and m != gate:
                    usel &= ~P[m]
            if usel.sum() < 2 or not ref.any():
                continue
            fg = float(np.quantile(V[gate][ref], 1 - tol))
            fa = float(np.quantile(V[ma][ref], 1 - tol))
            fb = float(np.quantile(V[mb][ref], 1 - tol))
            lg, la, lb = (float(lv.loc[m, 'level_log']) for m in (gate, ma, mb))
            # iterate with leftover DN cells leaving the universe as another population (owner 2026-11-06); in the
            # typing workflow those cells are not named by the lineage and continue through the regular cascade
            active = usel.copy()
            tg, ta, tb = lg, la, lb
            rounds = 0
            dn_removed = 0
            fit = None
            for rounds in range(1, 7):
                fit = threeway_fit(V[gate][active], V[ma][active], V[mb][active], tg, ta, tb, fg, fa, fb)
                same = (fit['tg'], fit['ta'], fit['tb']) == (tg, ta, tb)
                tg, ta, tb = fit['tg'], fit['ta'], fit['tb']
                dn_local = (V[gate][active] > tg) & (V[ma][active] <= ta) & (V[mb][active] <= tb)
                if not dn_local.any() or (same and rounds > 1) or fit['collapsed']:
                    break
                active[np.flatnonzero(active)[dn_local]] = False
                dn_removed += int(dn_local.sum())
            fit['tg'], fit['ta'], fit['tb'] = tg, ta, tb
            told = P[gate] & usel
            before = {'t_cells': int(told.sum()), 'dp': int((P[ma] & P[mb] & told).sum()),
                      'dn': int((~P[ma] & ~P[mb] & told).sum())}
            for m, new, old, floor in ((gate, tg, lg, fg), (ma, ta, la, fa), (mb, tb, lb, fb)):
                sel = (lv_new['image'] == tag) & (lv_new['marker'] == m)
                lv_new.loc[sel, 'level_log'] = new
                if 'level_raw' in lv_new.columns:
                    lv_new.loc[sel, 'level_raw'] = float(np.expm1(new))
                audit['images'].setdefault(tag, {})[m] = {
                    'lineage': f'{comp}/{lin}', 'level_log_old': old, 'level_log_new': new, 'floor_log': floor,
                    'floor_binding': bool(new <= floor)}
            t_final = int((V[gate][usel] > tg).sum())
            dn_final = int(((V[gate][usel] > tg) & (V[ma][usel] <= ta) & (V[mb][usel] <= tb)).sum())
            dp_final = int(((V[gate][usel] > tg) & (V[ma][usel] > ta) & (V[mb][usel] > tb)).sum())
            audit['images'].setdefault(tag, {})[f'{comp}/{lin}'] = {
                'universe': int(usel.sum()), 't_cells': {'before': before['t_cells'], 'after': t_final},
                'before': before, 'after': {'dp': dp_final, 'dn_leftover_to_regular_workflow': dn_final,
                                            'orphan_subtype': fit['orphan_subtype']},
                'rounds': rounds, 'dn_removed_as_population': dn_removed,
                'iterations': fit['iterations'], 'collapse_guard': fit['collapsed'],
                'sides': {ma: fit['a_side'], mb: fit['b_side']}}
            print(f"{tag} {comp}/{lin}: {gate}/{ma}/{mb} DP {before['dp']} -> {dp_final}, DN {before['dn']} -> "
                  f"{dn_final} ({dn_removed} removed over {rounds} rounds), orphans {fit['orphan_subtype']}"
                  f" of U {int(usel.sum())}", flush=True)
    lv_new.to_csv(run / 'typing' / 'levels.csv', index=False)
    (run / 'typing' / 't_partition_adjustment.json').write_text(json.dumps(audit, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
