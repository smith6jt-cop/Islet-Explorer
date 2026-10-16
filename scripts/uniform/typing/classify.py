#!/usr/bin/env python3
"""classify.py -- Senior's subpanel classification on this panel: context from geometry, one question per context.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/typing/subpanel.py` (classify and everything it
calls) and `subpanel_io.py`, as run on 2026-10-08 with `--membership geometry --emit-detail --refit` and
per-donor levels. One image per donor here, so "per donor" is "per image".

  levels      per donor: the p99 of max(scaled subpanel markers) over the donor's ORDINARY TISSUE cells (in no
              islet object and no structure object) -- for the hormones, for CHGA (Senior's CD99) and for each
              structure.
  1 STRUCTURE asked of every cell, islet included, weakest class first so Blood vessel > Nerve > Lymphatic:
              centroid inside an object of the class, MORE THAN HALF the cell's area inside the class's objects,
              the class's subpanel clearing the donor's level, and the marker that drew the objects being the
              cell's top expressed identity marker (correction C, on the aligned values before scaling).
              Otherwise the cell keeps its own identity and carries passed_through_<class>.
  2 ENDOCRINE inside an endocrine object only. A first pass admits a cell whose highest scaled hormone clears
              the donor's hormone level; the islet sizes it leaves decide each cell's stratum (1 / 2-3 / 4-10 /
              >10 admitted cells). Then per donor x stratum a CHGA level is fitted as the 1st percentile of that
              stratum's hormone-confirmed cells (complete hormones, top >= 0.35 and >= 2x the runner-up; at least
              500 of them; refused at the floor); in the >10 stratum nucleus elongation may also ADD. Where no
              level can be fitted the hormone question of the first pass stands. A complete, unique positive
              hormone winner names Beta / Alpha / Delta; incomplete or zero comparisons are
              "Endocrine (unresolved subtype)", positive ties "Endocrine (ambiguous subtype)".
  3 ORDINARY  every cell left: the immune definitions in name order (first firing wins; scaled >= 0.10, an
              unmeasured term is unknown, never negative), then the supporting tier (bottom-tier immune >= 0.8715
              while the exocrine markers stay below it), then Exocrine as the residual.
  4 FINISH    Exocrine inside an endocrine object is Excluded (exocrine_inside_islet_border, or no_marker_values
              when no exocrine marker was measured). An EndocrineSingle's cell is named by its unique highest
              positive scaled hormone, with no level test, or Excluded.

Departures (codex_typing.py): CHGA in CD99's place; the subpanels and rules on this panel's markers, with the
removed terms recorded in the audit; the Capillary objects (drawn on CD31) are Senior's Blood vessel; no hormone
fragments exist any more (step 1b), so no hormone is masked.

    python scripts/uniform/typing/classify.py --out data/development/uniform-20donor-20261015 --images all
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
import codex_typing as CT  # noqa: E402

STRUCTURES = CT.STRUCTURES
CROSSED = {'Blood vessel': 'crossed_frac_vessel', 'Nerve': 'crossed_frac_nerve', 'Lymphatic': 'crossed_frac_lymphatic'}
ORDINARY_ENDOCRINE_NAMES = ('Endocrine (unresolved subtype)', 'Endocrine (ambiguous subtype)')
SHAPE_FEATURE = 'nucleus_elongation'


# ----------------------------------------------------------------------------------------------- inputs
def load_scale(out: Path) -> pd.DataFrame:
    return pd.read_parquet(CT.typing_dir(out) / 'scale.parquet')


def apply_scale(values, markers, scale):
    lo = scale.loc[markers, 'floor'].to_numpy(float)
    span = scale.loc[markers, 'span'].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    return np.clip((np.asarray(values, dtype=np.float32) - lo) / span, 0.0, 1.0)


def scaled_markers(aligned: pd.DataFrame, wanted: set, scale: pd.DataFrame) -> pd.DataFrame:
    """Scaled values by name; an absent marker is a column of NaN. The unscaled values ride along as
    attrs['expressed'] for correction C."""
    present = sorted(m for m in wanted if m in aligned.columns and m in scale.index)
    values = pd.DataFrame(apply_scale(aligned[present].to_numpy(), present, scale), columns=present,
                          index=aligned.index)
    for m in sorted(wanted - set(present)):
        values[m] = np.nan
    raw = pd.DataFrame(aligned[present].to_numpy(dtype=float), columns=present, index=aligned.index)
    for m in sorted(wanted - set(present)):
        raw[m] = np.nan
    values.attrs['expressed'] = raw
    return values


def reference_level(pooled, quantile=CT.REFERENCE_QUANTILE) -> float:
    level = float(np.nanpercentile(pooled, quantile))
    if not np.isfinite(level):
        raise ValueError(f'reference level is not finite ({level!r}) -- refusing to classify')
    return level


def endocrine_level(values, quantile=CT.ENDOCRINE_QUANTILE, min_cells=CT.MIN_REFERENCE_CELLS, floor=0.0):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < min_cells:
        return None
    level = float(np.percentile(finite, quantile))
    if not np.isfinite(level):
        raise ValueError(f'endocrine level is not finite ({level!r})')
    return None if level <= floor else level


def size_stratum(islet_cells) -> np.ndarray:
    n = np.asarray(islet_cells, dtype=float)
    out = np.full(n.shape, '', dtype=object)
    known = np.isfinite(n)
    out[known] = CT.SIZE_STRATA[0]
    out[known & (n >= 2) & (n <= 3)] = CT.SIZE_STRATA[1]
    out[known & (n >= 4) & (n <= CT.LARGE_ISLET_MIN_CELLS)] = CT.SIZE_STRATA[2]
    out[known & (n > CT.LARGE_ISLET_MIN_CELLS)] = CT.SIZE_STRATA[3]
    return out


def fit_levels(ctx: pd.DataFrame, values: pd.DataFrame) -> dict:
    """subpanel.fit_levels for one donor (per-donor mode)."""
    tissue = ctx['context'].to_numpy() == 'tissue'
    levels = {}
    for name, markers in [('hormone', list(CT.HORMONES)), (CT.LEVEL_MARKER, [CT.LEVEL_MARKER])] + \
                         [(n, list(CT.SUBPANEL[n])) for n in STRUCTURES]:
        present = [m for m in markers if values[m].notna().any()]
        if not present:
            raise ValueError(f'{name}: no marker measured -- refusing to classify')
        levels[name] = reference_level(np.nanmax(values.loc[tissue, present].to_numpy(), axis=1))
    return levels


# -------------------------------------------------------------------------------------------- the rules
def rule_evidence(values: pd.DataFrame, rule: tuple):
    all_of, any_of, none_of = rule
    n = len(values)
    contradicted = np.zeros(n, bool)
    unknown = np.zeros(n, bool)
    for m in all_of:
        v = values[m].to_numpy(dtype=float) if m in values else np.full(n, np.nan)
        measured = np.isfinite(v)
        contradicted |= measured & (v < CT.IMMUNE_REQUIRE_MIN)
        unknown |= ~measured
    for group in any_of:
        if not group:
            contradicted[:] = True
            continue
        members = [m for m in group if m in values.columns]
        if not members:
            unknown[:] = True
            continue
        v = values[members].to_numpy(dtype=float)
        measured = np.isfinite(v)
        witness = (measured & (v >= CT.IMMUNE_REQUIRE_MIN)).any(axis=1)
        complete = measured.all(axis=1) & (len(members) == len(group))
        contradicted |= ~witness & complete
        unknown |= ~witness & ~complete
    for m in none_of:
        v = values[m].to_numpy(dtype=float) if m in values else np.full(n, np.nan)
        measured = np.isfinite(v)
        contradicted |= measured & (v >= CT.IMMUNE_REQUIRE_MIN)
        unknown |= ~measured
    return ~contradicted & ~unknown, ~contradicted & unknown


def immune_rules() -> dict:
    rules, _ = CT.codex_rules()
    return {k: v for k, v in rules.items() if k.startswith(('Imm-', 'T-'))}


def identity_argmax(values: pd.DataFrame) -> np.ndarray:
    """top_expressed_marker: the identity marker this cell EXPRESSES most, on the aligned values before scaling;
    '' where no identity marker was measured (correction C reads this column)."""
    raw = values.attrs['expressed']
    cols = [m for m in CT.IDENTITY_MARKERS if m in raw.columns]
    M = raw[cols].to_numpy(dtype=float)
    none = np.isnan(M).all(axis=1)
    top = np.asarray(cols, dtype=object)[np.argmax(np.where(np.isnan(M), -np.inf, M), axis=1)]
    return np.where(none, '', top)


def hormone_evidence(values: pd.DataFrame) -> dict:
    m = values.reindex(columns=list(CT.HORMONES)).to_numpy(dtype=float)
    finite = np.isfinite(m)
    filled = np.where(finite, m, -np.inf)
    order = np.sort(filled, axis=1)
    top, runner_up = order[:, -1], order[:, -2]
    which = np.asarray(CT.HORMONES, dtype=object)[np.argmax(filled, axis=1)]
    complete = finite.all(axis=1)
    top_count = ((filled == top[:, None]) & finite).sum(axis=1)
    top = np.where(top < 0, np.nan, top)
    runner_up = np.where(runner_up < 0, 0.0, runner_up)
    confirmed = complete & (top >= CT.CONFIRMED_MIN) & (top >= CT.CONFIRMED_MARGIN * np.maximum(runner_up, 1e-6))
    return {'top': top, 'runner_up': runner_up, 'which': which, 'measured': finite.any(axis=1),
            'complete': complete, 'confirmed': confirmed, 'top_count': top_count}


def fit_islet_levels(donor_level, confirmed, stratum, chga, shape) -> dict:
    out = {'donor_reference_level': donor_level}
    for name in CT.SIZE_STRATA:
        pool = confirmed & (stratum == name)
        entry = {'n_confirmed': int(pool.sum()), CT.LEVEL_MARKER: None, SHAPE_FEATURE: None}
        if pool.any():
            entry[CT.LEVEL_MARKER] = endocrine_level(chga[pool])
            if name == CT.SIZE_STRATA[-1]:
                entry[SHAPE_FEATURE] = endocrine_level(shape[pool], floor=CT.SHAPE_FLOOR)
        out[name] = entry
    return out


def endocrine_decision(open_islet, stratum, chga, elongation, islet_levels, measured) -> dict:
    n = len(open_islet)
    open_islet = np.asarray(open_islet) & np.asarray(measured)
    level = np.array([islet_levels.get(s, {}).get(CT.LEVEL_MARKER) for s in stratum], dtype=object)
    have_level = np.array([v is not None for v in level]) & open_islet
    chga_ok = np.zeros(n, bool)
    if have_level.any():
        chga_ok[have_level] = chga[have_level] >= np.array(list(level[have_level]), dtype=float)
    large = np.asarray(stratum) == CT.SIZE_STRATA[-1]
    shape_level = islet_levels.get(CT.SIZE_STRATA[-1], {}).get(SHAPE_FEATURE)
    shape_ok = (np.zeros(n, bool) if shape_level is None
                else large & have_level & (np.asarray(elongation, dtype=float) >= shape_level))
    return {'have_level': have_level, 'chga_ok': chga_ok, 'shape_ok': shape_ok,
            'endocrine': have_level & (chga_ok | shape_ok), 'large': large}


def structure_calls(ctx, values, levels, top_marker):
    n = len(ctx)
    cell_class = np.full(n, None, dtype=object)
    resolved_by = np.full(n, '', dtype=object)
    passed_through, audit = {}, {}
    for name in reversed(STRUCTURES):
        markers = [m for m in CT.SUBPANEL[name] if m in values.columns]
        inside = ctx[f'in_{name}'].to_numpy()
        clears = np.nanmax(values[markers].to_numpy(), axis=1) >= levels[name]
        owns = np.isin(top_marker, CT.CLASSIFIER_MARKER[name])
        frac = ctx[CROSSED[name]].to_numpy(dtype=float)
        known = ~np.isnan(frac)
        majority = known & (frac > CT.MAJORITY_AREA)
        candidate = inside & majority
        hit = candidate & clears & owns
        cell_class[hit] = name
        resolved_by[hit] = 'structure'
        passed_through[name] = np.where(inside & ~hit, frac, np.nan)
        audit[name] = {'inside': int(inside.sum()), 'inside_and_majority': int(candidate.sum()),
                       'majority_withholds': int((inside & known & ~majority).sum()),
                       'inside_and_clears': int((candidate & clears).sum()), 'called': int(hit.sum()),
                       'top_marker_withholds': int((candidate & clears & ~owns).sum()),
                       'top_marker_withholds_in_islet': int((candidate & clears & ~owns
                                                             & ctx['in_islet'].to_numpy()).sum())}
    return cell_class, resolved_by, passed_through, audit


def tissue_chain(values, rules, cell_class, definition, resolved_by, rule_audit=None, unavailable_count=None):
    n = len(cell_class)
    unclaimed = cell_class == None  # noqa: E711
    for phenotype, rule in sorted(rules.items()):
        satisfied, unavailable = rule_evidence(values, rule)
        hit = unclaimed & satisfied
        if unavailable_count is not None:
            unavailable_count += unavailable.astype(unavailable_count.dtype)
        if rule_audit is not None:
            rule_audit[phenotype] = {'n_firing': int(satisfied.sum()), 'n_called': int(hit.sum()),
                                     'n_unavailable': int(unavailable.sum())}
        cell_class[hit] = 'Immune'
        definition[hit] = phenotype
        resolved_by[hit] = 'definition'
        unclaimed &= ~hit
    # this panel's added rules, after every Senior definition (codex_typing.ADDED_IMMUNE_RULES)
    for phenotype, rule in CT.ADDED_IMMUNE_RULES.items():
        satisfied, unavailable = rule_evidence(values, rule)
        hit = unclaimed & satisfied
        if rule_audit is not None:
            rule_audit[phenotype] = {'n_firing': int(satisfied.sum()), 'n_called': int(hit.sum()),
                                     'n_unavailable': int(unavailable.sum()), 'added_for_this_panel': True}
        cell_class[hit] = 'Immune'
        definition[hit] = phenotype
        resolved_by[hit] = 'definition_added'
        unclaimed &= ~hit
    tier = [m for m in CT.BOTTOM_TIER_IMMUNE if m in values.columns]
    exocrine = [m for m in CT.SUBPANEL['exocrine'] if m in values.columns]
    if tier:
        clears = np.nanmax(values[tier].to_numpy(), axis=1) >= CT.SUPPORTING_LEVEL
        speaks_exocrine = (np.nanmax(values[exocrine].to_numpy(), axis=1) >= CT.SUPPORTING_LEVEL
                           if exocrine else np.zeros(n, bool))
        hit = unclaimed & clears & ~speaks_exocrine
        cell_class[hit] = 'Immune'
        resolved_by[hit] = 'supporting'
        unclaimed &= ~hit
    cell_class[unclaimed] = 'Exocrine'
    resolved_by[unclaimed] = 'fallthrough'


def finish(in_islet, cell_class, exocrine_measured):
    reason = np.full(len(cell_class), '', dtype=object)
    ruled_out = in_islet & (cell_class == 'Exocrine')
    cell_class[ruled_out] = 'Excluded'
    reason[ruled_out] = CT.ISLET_EXOCRINE_REASON
    reason[ruled_out & ~np.asarray(exocrine_measured)] = CT.NO_MEASUREMENT_REASON
    return reason


def run_chain(values, rules, structure_class, structure_resolved, in_islet, endocrine, endocrine_name,
              endocrine_route, exocrine_measured, rule_audit=None, unavailable_count=None):
    cell_class = structure_class.copy()
    resolved_by = structure_resolved.copy()
    definition = np.full(len(cell_class), '', dtype=object)
    hit = endocrine & (cell_class == None)  # noqa: E711
    cell_class[hit] = endocrine_name[hit]
    resolved_by[hit] = endocrine_route[hit]
    tissue_chain(values, rules, cell_class, definition, resolved_by, rule_audit, unavailable_count)
    reason = finish(in_islet, cell_class, exocrine_measured)
    return cell_class, definition, resolved_by, reason


def direct_single_identity(values):
    markers = [m for m in CT.HORMONES if m in values]
    names = np.full(len(values), 'Excluded', dtype=object)
    status = np.full(len(values), 'no_usable_hormone', dtype=object)
    x = values[markers].to_numpy(float)
    finite = np.isfinite(x)
    safe = np.where(finite, x, -np.inf)
    top = safe.max(axis=1)
    usable = finite.any(axis=1) & (top > 0)
    tied = ((safe == top[:, None]) & finite).sum(axis=1) > 1
    unique = usable & ~tied
    winner = safe.argmax(axis=1)
    names[unique] = np.asarray([CT.HORMONE_NAMES[m] for m in markers], dtype=object)[winner[unique]]
    status[unique] = 'unique_highest_hormone'
    status[usable & tied] = 'tied_highest_hormones'
    return names, status


def ordinary_endocrine_identity(evidence):
    top, complete, top_count, which = (evidence[k] for k in ('top', 'complete', 'top_count', 'which'))
    names = np.full(len(top), ORDINARY_ENDOCRINE_NAMES[0], dtype=object)
    status = np.full(len(top), 'incomplete_hormone_comparison', dtype=object)
    positive = complete & np.isfinite(top) & (top > 0)
    tied = positive & (top_count > 1)
    unique = positive & (top_count == 1)
    names[tied] = ORDINARY_ENDOCRINE_NAMES[1]
    names[unique] = np.asarray([CT.HORMONE_NAMES[m] for m in which[unique]], dtype=object)
    status[complete & ~positive] = 'no_usable_hormone'
    status[tied] = 'tied_highest_hormones'
    status[unique] = 'unique_highest_hormone'
    return names, status


# --------------------------------------------------------------------------------------------- classify
def classify(tag: str, ctx: pd.DataFrame, aligned: pd.DataFrame, scale: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rules = immune_rules()
    wanted = set(CT.HORMONES) | {CT.LEVEL_MARKER} | set(CT.IDENTITY_MARKERS) | set(CT.BOTTOM_TIER_IMMUNE)
    for markers in CT.SUBPANEL.values():
        wanted |= set(markers)
    for all_of, any_of, none_of in list(rules.values()) + list(CT.ADDED_IMMUNE_RULES.values()):
        wanted |= set(all_of) | set(none_of) | {m for g in any_of for m in g}
    values = scaled_markers(aligned, wanted, scale)
    n = len(ctx)
    in_islet = ctx['in_islet'].to_numpy()
    A = {'image': tag, 'n_cells': n, 'n_in_islet': int(in_islet.sum())}
    levels = fit_levels(ctx, values)
    A['levels'] = levels

    top_marker = identity_argmax(values)
    A['top_expressed_marker'] = {str(k): int(v) for k, v in pd.Series(top_marker).value_counts().items()}

    structure_class, structure_resolved, passed_through, A['structures'] = structure_calls(ctx, values, levels,
                                                                                          top_marker)
    structural = structure_class != None  # noqa: E711

    horm = hormone_evidence(values)
    exo = [m for m in CT.SUBPANEL['exocrine'] if m in values.columns]
    exo_measured = values[exo].notna().any(axis=1).to_numpy()
    base_endocrine = in_islet & horm['measured'] & (horm['top'] >= levels['hormone'])
    endocrine_name, endocrine_status = ordinary_endocrine_identity(horm)
    endocrine_names = set(CT.HORMONE_NAMES.values()) | set(ORDINARY_ENDOCRINE_NAMES)
    base_class, _, _, _ = run_chain(values, rules, structure_class, structure_resolved, in_islet, base_endocrine,
                                    endocrine_name, np.full(n, 'hormone', dtype=object), exo_measured)

    admitted = in_islet & ~structural & (base_class != 'Excluded')
    islet_num = ctx['islet_num'].to_numpy()
    sizes = pd.Series(islet_num[admitted]).value_counts()
    islet_cells = pd.Series(islet_num).map(sizes).to_numpy(dtype=float)
    islet_cells = np.where(in_islet, np.nan_to_num(islet_cells, nan=0.0), np.nan)
    stratum = size_stratum(islet_cells)
    A['n_admitted'] = int(admitted.sum())
    A['stratum_cells'] = {s: int((stratum == s).sum()) for s in CT.SIZE_STRATA}

    chga = values[CT.LEVEL_MARKER].to_numpy(dtype=float)
    elongation = ctx[SHAPE_FEATURE].to_numpy(dtype=float)
    confirmed_pool = admitted & horm['confirmed']
    islet_levels = fit_islet_levels(levels.get(CT.LEVEL_MARKER), confirmed_pool, stratum, chga, elongation)
    A['islet_levels'] = islet_levels
    open_islet = in_islet & ~structural
    D = endocrine_decision(open_islet, stratum, chga, elongation, islet_levels, horm['measured'])
    fallback = open_islet & ~D['have_level']
    endocrine = D['endocrine'] | (fallback & base_endocrine)
    route = np.full(n, '', dtype=object)
    route[endocrine & D['chga_ok']] = 'chga'
    route[endocrine & ~D['chga_ok'] & D['shape_ok']] = 'elongation'
    route[endocrine & fallback] = 'hormone'

    rule_audit = {}
    unavailable_count = np.zeros(n, dtype=np.int16)
    cell_class, definition, resolved_by, reason = run_chain(
        values, rules, structure_class, structure_resolved, in_islet, endocrine, endocrine_name, route,
        exo_measured, rule_audit, unavailable_count)
    A['immune_rule_evidence'] = rule_audit

    single = ctx['endocrine_object_class'].eq('EndocrineSingle').to_numpy()
    single_name, single_status = direct_single_identity(values)
    cell_class[single] = single_name[single]
    resolved_by[single] = 'single_hormone'
    route[single] = 'single_hormone'
    endocrine[single] = True
    definition[single] = ''
    reason[single] = ''
    single_excluded = single & (single_status != 'unique_highest_hormone')
    cell_class[single_excluded] = 'Excluded'
    resolved_by[single_excluded] = 'excluded_single_hormone'
    route[single_excluded] = ''
    endocrine[single_excluded] = False
    reason[single_excluded] = np.asarray(['endocrine_single_' + str(s) for s in single_status[single_excluded]],
                                         dtype=object)
    A['direct_endocrine_singles'] = {str(k): int(v) for k, v in pd.Series(single_status[single]).value_counts().items()}
    ordinary = ~single & np.isin(cell_class, list(endocrine_names))
    subtype_status = np.where(ordinary, endocrine_status, '')
    subtype_status[single] = single_status[single]

    A['n_endocrine_base'] = int((base_endocrine & open_islet).sum())
    A['n_endocrine'] = int(endocrine.sum())
    A['n_by_chga'] = int((endocrine & D['chga_ok']).sum())
    A['n_by_elongation_only'] = int((endocrine & ~D['chga_ok'] & D['shape_ok']).sum())
    A['n_by_hormone_fallback'] = int((endocrine & fallback).sum())
    A['n_immune_to_endocrine'] = int(((base_class == 'Immune') & np.isin(cell_class, list(endocrine_names))).sum())
    A['cell_class'] = {str(k): int(v) for k, v in pd.Series(cell_class).value_counts().items()}
    A['cell_class_in_islet'] = {str(k): int(v) for k, v in pd.Series(cell_class[in_islet]).value_counts().items()}
    A['immune_definition'] = {str(k): int(v) for k, v in pd.Series(definition[definition != '']).value_counts().items()}

    frame = pd.DataFrame({
        'object_id': ctx['object_id'].to_numpy(), 'cell_class': cell_class.astype(str),
        'context': ctx['context'].to_numpy(), 'object_context': ctx['object_context'].to_numpy(),
        'excluded_reason': reason.astype(str), 'islet_num': ctx['islet_num'].to_numpy(),
        'endocrine_object_class': ctx['endocrine_object_class'].to_numpy(),
        'endocrine_object_id': ctx['endocrine_object_id'].to_numpy(),
        'islet_cells': islet_cells, 'islet_stratum': stratum, 'islet_route': route,
        'immune_definition': definition.astype(str), 'resolved_by': resolved_by.astype(str),
        'top_expressed_marker': top_marker, 'endocrine_subtype_status': subtype_status,
        'immune_rules_unavailable_count': unavailable_count})
    for name in STRUCTURES:
        frame[f'passed_through_{name.replace(" ", "_")}'] = passed_through[name]
    return frame, A


def process(out: Path, tag: str) -> dict:
    typ = CT.typing_dir(out)
    scale = load_scale(out)
    ctx = pd.read_parquet(typ / 'context' / f'{tag}.parquet')
    aligned = pd.read_parquet(typ / 'aligned' / f'{tag}.parquet')
    if not ctx['object_id'].equals(aligned['object_id']):
        aligned = ctx[['object_id']].merge(aligned, on='object_id', how='left', validate='one_to_one')
    keep = ctx['qc_analysis_eligible'].to_numpy(bool)
    ctx = ctx[keep].reset_index(drop=True)
    aligned = aligned[keep].reset_index(drop=True)
    frame, audit = classify(tag, ctx, aligned, scale)
    audit['n_not_analysis_eligible'] = int((~keep).sum())
    dest = typ / 'calls' / f'{tag}.parquet'
    dest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(dest, index=False)
    dest.with_suffix('.json').write_text(json.dumps(audit, indent=1, default=float))
    print(f'{tag}: {audit["cell_class"]}', flush=True)
    return audit


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    tags = [r['image_name'] for r in CT.cohort_rows()] if a.images == ['all'] else a.images
    _, record = CT.codex_rules()
    (CT.typing_dir(a.out) / 'immune_rules_on_this_panel.json').write_text(json.dumps(record, indent=1))
    with ProcessPoolExecutor(a.workers) as pool:
        list(pool.map(process, [a.out] * len(tags), tags))
    return 0


if __name__ == '__main__':
    sys.exit(main())
