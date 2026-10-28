#!/usr/bin/env python3
"""immune_subtype.py -- Senior's immune class and subtype, inside the cells already called Immune.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/typing/immune_subtype.py` (evaluate, resolve): for
cells whose cell_class is Immune, level 2 is the rule table's children of Immune and level 3 the children of
Imm-T-cell and Imm-macrophage; a rule fires when every required term is at scaled >= 0.10, one member of each OR
group is, and every forbidden term is below it; a cell whose value is missing on a rule's marker is undecided for
that rule; the most specific firing rule wins (number of terms, then name), and every other match is recorded.
cell_class is never changed.
Positivity is each marker above its image's own level (typing/positivity.py, owner 2026-10-28), not scaled >= 0.10.
A T cell that neither T-CD8 nor T-CD4-helper names (both markers positive, or neither) takes the subtype of the
marker further above its own level (owner, 2026-10-27;
codex_typing.T_SUBTYPE_BY_HIGHER; immune_subtype_resolved_by). A cell classify.py withdrew from the T-cell call on
its whole profile (resolved_by 'tie_profile_*') is resolved without Imm-T-cell.

CODEX departures: the rules on this panel (codex_typing.codex_rules): a rule that needs a marker this panel never
acquired is unavailable on every donor and is listed, never scored as negative. The per-rule term diagnostics
and the reduced-form counts Senior prints are not ported.

    python scripts/uniform/typing/immune_subtype.py --out data/development/uniform-20donor-20261015
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


def rule_markers(rule):
    all_of, any_of, none_of = rule
    return sorted(set(all_of) | set(none_of) | {m for g in any_of for m in g})


def evaluate(scaled: pd.DataFrame, rule, available: set):
    """`scaled` holds each marker's excess over its image's positivity level (classify.excess): > 0 is positive."""
    all_of, any_of, none_of = rule
    needed = set(rule_markers(rule))
    if not needed <= available:
        return None, None, sorted(needed - available)
    block = scaled[sorted(needed)].to_numpy()
    undecided = ~np.isfinite(block).all(axis=1)
    fires = np.ones(len(scaled), bool)
    for m in all_of:
        fires &= scaled[m].to_numpy() > 0
    for g in any_of:
        hit = np.zeros(len(scaled), bool)
        for m in g:
            hit |= scaled[m].to_numpy() > 0
        fires &= hit
    for m in none_of:
        fires &= ~(scaled[m].to_numpy() > 0)
    return fires & ~undecided, undecided, []


def resolve(n, candidates, rules, scaled, available, specificity):
    hits, undecided, skipped = {}, np.zeros(n, bool), {}
    for name in candidates:
        if name not in rules:
            skipped[name] = 'not on this panel'
            continue
        fires, undec, missing = evaluate(scaled, rules[name], available)
        if fires is None:
            skipped[name] = missing
            continue
        hits[name] = fires
        undecided |= undec
    if not hits:
        return np.full(n, None, dtype=object), np.zeros(n, int), np.full(n, '', dtype=object), undecided, skipped
    order = sorted(hits, key=lambda k: (-specificity[k], k))
    stacked = np.vstack([hits[k] for k in order])
    n_match = stacked.sum(axis=0)
    first = np.argmax(stacked, axis=0)
    call = np.where(n_match > 0, np.array(order, dtype=object)[first], None)
    others = np.array([';'.join(k for j, k in enumerate(order) if stacked[j, i] and j != first[i])
                       if n_match[i] > 1 else '' for i in range(n)], dtype=object)
    undecided &= ~(n_match > 0)
    return call, n_match, others, undecided, skipped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', default=['all'])
    a = ap.parse_args()
    typ = CT.typing_dir(a.out)
    tags = [r['image_name'] for r in CT.cohort_rows()] if a.images == ['all'] else a.images
    rules, record = CT.codex_rules()
    specificity = {k: len(r[0]) + sum(len(g) for g in r[1]) + len(r[2]) for k, r in rules.items()}
    # this panel's added rules rank below every Senior rule: they name a cell only when no Senior rule does
    for k, r in CT.ADDED_IMMUNE_RULES.items():
        rules[k] = r
        specificity[k] = -1
    level2 = tuple(CT.IMMUNE_LEVEL2) + tuple(CT.ADDED_IMMUNE_RULES)
    scale = CL.load_scale(a.out)
    (typ / 'immune').mkdir(parents=True, exist_ok=True)
    summary = {}
    for t in tags:
        calls = pd.read_parquet(typ / 'calls' / f'{t}.parquet', columns=['object_id', 'cell_class', 'resolved_by',
                                                                          'immune_definition'])
        immune = calls[calls['cell_class'] == 'Immune'][['object_id']]
        aligned = pd.read_parquet(typ / 'aligned' / f'{t}.parquet')
        markers = [m for m in scale.index if m in aligned.columns
                   and np.isfinite(aligned[m].to_numpy(float)).any()]
        values = immune.merge(aligned[['object_id'] + markers], on='object_id', how='left', validate='one_to_one')
        # every marker as its excess over the image's positivity level (typing/positivity.py; > 0 is positive)
        scaled = CL.excess(values, markers, CL.load_levels(a.out, t)).reset_index(drop=True)
        n = len(values)
        call2, n2, other2, undec2, skip2 = resolve(n, level2, rules, scaled, set(markers), specificity)
        # a T cell withdrawn on its whole profile (classify.retype_tied_t_cells) is never a T cell here either
        withdrawn = values['object_id'].isin(set(calls.loc[calls['resolved_by'].str.startswith('tie_profile_'),
                                                         'object_id'])).to_numpy()
        if withdrawn.any():
            c2, m2, o2, u2, _ = resolve(int(withdrawn.sum()), tuple(k for k in level2 if k != 'Imm-T-cell'), rules,
                                        scaled[withdrawn].reset_index(drop=True), set(markers), specificity)
            call2, n2, other2, undec2 = call2.copy(), n2.copy(), other2.copy(), undec2.copy()
            call2[withdrawn], n2[withdrawn], other2[withdrawn], undec2[withdrawn] = c2, m2, o2, u2
        # a cell placed by its top marker (classify.tissue_chain, owner 2026-10-27) keeps that class here
        placed = calls.loc[calls['resolved_by'].str.endswith('top_marker')]     # also 'tie_profile_top_marker'
        placed = placed.set_index('object_id')['immune_definition']
        by_top = values['object_id'].map(placed).fillna('').to_numpy(dtype=object)
        call2 = np.where(by_top != '', by_top, call2)
        out = pd.DataFrame({'object_id': values['object_id'], 'immune_class': call2, 'immune_class_matches': n2,
                            'immune_class_also': other2, 'immune_class_undecided': undec2,
                            'immune_subtype': None, 'immune_subtype_matches': 0, 'immune_subtype_also': '',
                            'immune_subtype_undecided': False})
        skipped = dict(skip2)
        for parent, kids in CT.IMMUNE_LEVEL3.items():
            mask = (out['immune_class'] == parent).to_numpy()
            if not mask.any():
                continue
            sub, n3, other3, undec3, skip3 = resolve(int(mask.sum()), kids, rules,
                                                     scaled[mask].reset_index(drop=True), set(markers), specificity)
            out.loc[mask, 'immune_subtype'] = sub
            out.loc[mask, 'immune_subtype_matches'] = n3
            out.loc[mask, 'immune_subtype_also'] = other3
            out.loc[mask, 'immune_subtype_undecided'] = undec3
            skipped.update(skip3)
        (s8, m8), (s4, m4) = CT.T_SUBTYPE_BY_HIGHER.items()
        v8, v4 = scaled[m8].to_numpy(), scaled[m4].to_numpy()          # each above its own level
        both = ((out['immune_class'] == 'Imm-T-cell') & out['immune_subtype'].isna()).to_numpy() \
            & np.isfinite(v8) & np.isfinite(v4)
        out['immune_subtype_resolved_by'] = np.where(out['immune_subtype'].notna(), 'rule', '')
        for name, hit in ((s8, both & (v8 > v4)), (s4, both & (v4 > v8))):
            out.loc[hit, 'immune_subtype'] = name
            out.loc[hit, 'immune_subtype_resolved_by'] = f'higher_of_{m8}_{m4}'
        # the macrophage label (codex_typing.MACROPHAGE_PROFILE), on every Imm-macrophage cell
        prof = [m for pair in CT.MACROPHAGE_PROFILE for m in pair]
        pos = {m: scaled[m].to_numpy() > 0 for m in prof}
        out['macrophage_profile'] = [CT.macrophage_label({m: bool(pos[m][i]) for m in prof})
                                     if out['immune_class'].iat[i] == 'Imm-macrophage' else ''
                                     for i in range(n)]
        out.to_parquet(typ / 'immune' / f'{t}.parquet', index=False)
        summary[t] = {'n_immune': n,
                      'immune_class': {str(k): int(v) for k, v in out['immune_class'].fillna('(no rule fired)')
                                       .value_counts().items()},
                      'immune_subtype': {str(k): int(v) for k, v in out['immune_subtype'].dropna()
                                         .value_counts().items()},
                      'subtype_by_higher_marker': int((out['immune_subtype_resolved_by'] == f'higher_of_{m8}_{m4}').sum()),
                      'tied': int((both & (v8 == v4)).sum()),
                      'macrophage_profile': {str(k): int(v) for k, v in out.loc[out['macrophage_profile'] != '',
                                             'macrophage_profile'].value_counts().items()},
                      'unavailable_rules': {k: v for k, v in skipped.items()}}
        print(f'{t}: {summary[t]["immune_class"]} | {summary[t]["immune_subtype"]}', flush=True)
    (typ / 'immune' / 'summary.json').write_text(json.dumps({'rules_on_this_panel': record, 'images': summary},
                                                            indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
