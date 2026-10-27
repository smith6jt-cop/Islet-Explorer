#!/usr/bin/env python3
"""immune_subtype.py -- Senior's immune class and subtype, inside the cells already called Immune.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/typing/immune_subtype.py` (evaluate, resolve): for
cells whose cell_class is Immune, level 2 is the rule table's children of Immune and level 3 the children of
Imm-T-cell and Imm-macrophage; a rule fires when every required term is at scaled >= 0.10, one member of each OR
group is, and every forbidden term is below it; a cell whose value is missing on a rule's marker is undecided for
that rule; the most specific firing rule wins (number of terms, then name), and every other match is recorded.
cell_class is never changed.
A T cell that neither T-CD8 nor T-CD4-helper names (both markers positive, or either below 0.10) takes the subtype
of the higher scaled marker at any level (owner, 2026-10-27; codex_typing.T_SUBTYPE_BY_HIGHER;
immune_subtype_resolved_by); an exact tie keeps none.

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
    all_of, any_of, none_of = rule
    needed = set(rule_markers(rule))
    if not needed <= available:
        return None, None, sorted(needed - available)
    block = scaled[sorted(needed)].to_numpy()
    undecided = ~np.isfinite(block).all(axis=1)
    fires = np.ones(len(scaled), bool)
    for m in all_of:
        fires &= scaled[m].to_numpy() >= CT.IMMUNE_REQUIRE_MIN
    for g in any_of:
        hit = np.zeros(len(scaled), bool)
        for m in g:
            hit |= scaled[m].to_numpy() >= CT.IMMUNE_REQUIRE_MIN
        fires &= hit
    for m in none_of:
        fires &= scaled[m].to_numpy() < CT.IMMUNE_REQUIRE_MIN
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
        calls = pd.read_parquet(typ / 'calls' / f'{t}.parquet', columns=['object_id', 'cell_class'])
        immune = calls[calls['cell_class'] == 'Immune'][['object_id']]
        aligned = pd.read_parquet(typ / 'aligned' / f'{t}.parquet')
        markers = [m for m in scale.index if m in aligned.columns
                   and np.isfinite(aligned[m].to_numpy(float)).any()]
        values = immune.merge(aligned[['object_id'] + markers], on='object_id', how='left', validate='one_to_one')
        scaled = pd.DataFrame(CL.apply_scale(values[markers].to_numpy(), markers, scale), columns=markers)
        n = len(values)
        call2, n2, other2, undec2, skip2 = resolve(n, level2, rules, scaled, set(markers), specificity)
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
        v8, v4 = scaled[m8].to_numpy(), scaled[m4].to_numpy()
        both = ((out['immune_class'] == 'Imm-T-cell') & out['immune_subtype'].isna()).to_numpy() \
            & np.isfinite(v8) & np.isfinite(v4)
        out['immune_subtype_resolved_by'] = np.where(out['immune_subtype'].notna(), 'rule', '')
        for name, hit in ((s8, both & (v8 > v4)), (s4, both & (v4 > v8))):
            out.loc[hit, 'immune_subtype'] = name
            out.loc[hit, 'immune_subtype_resolved_by'] = f'higher_of_{m8}_{m4}'
        out.to_parquet(typ / 'immune' / f'{t}.parquet', index=False)
        summary[t] = {'n_immune': n,
                      'immune_class': {str(k): int(v) for k, v in out['immune_class'].fillna('(no rule fired)')
                                       .value_counts().items()},
                      'immune_subtype': {str(k): int(v) for k, v in out['immune_subtype'].dropna()
                                         .value_counts().items()},
                      'subtype_by_higher_marker': int((out['immune_subtype_resolved_by'] == f'higher_of_{m8}_{m4}').sum()),
                      'tied': int((both & (v8 == v4)).sum()),
                      'unavailable_rules': {k: v for k, v in skipped.items()}}
        print(f'{t}: {summary[t]["immune_class"]} | {summary[t]["immune_subtype"]}', flush=True)
    (typ / 'immune' / 'summary.json').write_text(json.dumps({'rules_on_this_panel': record, 'images': summary},
                                                            indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
