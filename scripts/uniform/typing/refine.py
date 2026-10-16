#!/usr/bin/env python3
"""refine.py -- Senior's refinement pass, applied once: a class is kept only where its own defining marker is.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/typing/refine.py` (representation_params, to_Z,
neighbours, concordance, fit_bounds, _demote, route, immune_leaf, refine_donor, load_cells):

  representation  cohort-wide per marker: arcsinh(x / (median / 5)), standardised (seeded 200,000 cells per
                  donor), over the rule table's markers on this panel (REFINE_MARKERS)
  neighbours      K = 30 nearest in that space, within one donor (pynndescent, seed 0)
  concordance     the share of a cell's neighbours carrying its label; >= 0.40 is "clear"
  bound           per donor x islet stratum x class: the 5th percentile of max(OWN markers of the class) over its
                  clear members, at least 50 of them (else none, and the group is left alone)
  demotion        below the bound: Exocrine -> "Unresolved (no lineage evidence)"; any other class -> Immune if
                  the cell meets an immune leaf rule (scaled >= 0.10), else Exocrine

`refined_class` is written BESIDE `cell_class`; nothing else is changed. Excluded cells are outside the pass.

CODEX departures: the representation uses the 23 of Senior's 54 rule-table markers this panel carries (CHGA for
CD99); OWN_MARKERS on this panel (codex_typing.OWN_MARKERS); Senior's CD8 exception (unshifted CD8 in refinement,
from its own CD8 reagent change) has no counterpart here and is not ported.

    python scripts/uniform/typing/refine.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402
import classify as CL  # noqa: E402

LOW = 0.40
K = 30
BOUND_PCT = 5.0
SEED = 0
MIN_CLEAR = 50
N_PER_DONOR = 200_000
UNRESOLVED = 'Unresolved (no lineage evidence)'
UNTOUCHED = ('Excluded',)
STRUCTURAL = set(CT.STRUCTURES)
# Senior's derive_markers() (54 markers, 2026-10-16) intersected with this panel and mapped to its names.
SENIOR_REFINE_MARKERS = ['B3TUBB', 'Bcl_2', 'Beta_actin', 'CD107a', 'CD11b', 'CD11c', 'CD163', 'CD20', 'CD206',
                         'CD209', 'CD31', 'CD34', 'CD39', 'CD3e', 'CD4', 'CD44', 'CD56', 'CD66', 'CD68', 'CD79a',
                         'CD8', 'CD99', 'Caveolin', 'E_cadherin', 'EpCAM', 'FOXP3', 'GCG', 'Granzyme_B', 'HLA_A',
                         'HLA_DR', 'ICOS', 'IDO1', 'INS', 'Iba1', 'Ker8_18', 'Keratin_5', 'Ki67', 'LAG3', 'M2Gal3',
                         'MPO', 'PCNA', 'PD_1', 'PD_L1', 'Pan_Cytokeratin', 'Podoplanin', 'SMA', 'SOX2', 'SST',
                         'TCF_1', 'TOX', 'TP63', 'VISTA', 'Vimentin', 'iNOS']
REFINE_MARKERS = [CT.SENIOR_TO_CODEX[m] for m in SENIOR_REFINE_MARKERS
                  if m in CT.SENIOR_TO_CODEX and CT.SENIOR_TO_CODEX[m] not in CT.NOT_ON_EVERY_IMAGE]


class RefinementAlreadyApplied(RuntimeError):
    pass


def load_cells(out: Path, tag: str) -> pd.DataFrame:
    df = pd.read_parquet(CT.typing_dir(out) / 'calls' / f'{tag}.parquet',
                         columns=['object_id', 'cell_class', 'islet_num'])
    core = df['islet_num'].notna() & ~df['cell_class'].isin(STRUCTURAL) & (df['cell_class'] != 'Excluded')
    size = df.loc[core].groupby('islet_num')['object_id'].size()
    n = df['islet_num'].map(size).where(df['islet_num'].notna())
    df['islet_cells'] = n.to_numpy(dtype=float)
    df['size_stratum'] = CL.size_stratum(df['islet_cells'].to_numpy())
    return df


def donor_values(out: Path, tag: str, markers: list[str]) -> tuple[pd.DataFrame, list[str]]:
    f = pd.read_parquet(CT.typing_dir(out) / 'aligned' / f'{tag}.parquet')
    have = [m for m in markers if m in f.columns]
    return f[['object_id', 'qc_analysis_eligible'] + have], have


def representation_params(out: Path, tags: list[str], markers: list[str]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    pool = {m: [] for m in markers}
    for d in tags:
        t, have = donor_values(out, d, markers)
        t = t[t['qc_analysis_eligible']]
        if len(t) > N_PER_DONOR:
            t = t.iloc[rng.choice(len(t), N_PER_DONOR, replace=False)]
        for m in have:
            pool[m].append(t[m].to_numpy('float64'))
    rows = {}
    for m, chunks in pool.items():
        if not chunks:
            continue
        v = np.concatenate(chunks)
        cof = np.nanmedian(v) / 5.0
        if not np.isfinite(cof) or cof <= 0:
            cof = 1.0
        z = np.arcsinh(v / cof)
        rows[m] = {'cofactor': float(cof), 'centre': float(np.nanmean(z)), 'spread': float(np.nanstd(z)),
                   'n': int(np.isfinite(v).sum())}
    out_ = pd.DataFrame(rows).T
    out_.index.name = 'marker'
    return out_


def to_Z(values, markers, params):
    cof = params.loc[markers, 'cofactor'].to_numpy(float)
    centre = params.loc[markers, 'centre'].to_numpy(float)
    spread = params.loc[markers, 'spread'].to_numpy(float)
    spread = np.where(np.isfinite(spread) & (spread > 0), spread, 1.0)
    Z = (np.arcsinh(np.asarray(values, dtype='float64') / cof) - centre) / spread
    if not np.isfinite(Z).all():
        raise ValueError('non-finite value in the representation')
    return Z.astype('float32')


def neighbours(Z, k=K, seed=SEED):
    from pynndescent import NNDescent
    idx = NNDescent(Z, n_neighbors=k + 1, metric='euclidean', random_state=seed, n_jobs=-1)
    nbr, _ = idx.neighbor_graph
    return nbr[:, 1:]


def concordance(nbr, labels):
    return (labels[nbr] == labels[:, None]).mean(axis=1)


def own_marker_max(Z, markers, cls):
    cols = [markers.index(m) for m in CT.OWN_MARKERS[cls] if m in markers]
    return Z[:, cols].max(axis=1) if cols else None


def fit_bounds(Z, markers, labels, strata, conc):
    out = {}
    for cls in CT.OWN_MARKERS:
        v = own_marker_max(Z, markers, cls)
        if v is None:
            continue
        in_class = labels == cls
        for st in (pd.unique(strata[in_class]) if in_class.any() else []):
            m = in_class & (strata == st)
            clear = m & (conc >= LOW)
            bound = float(np.percentile(v[clear], BOUND_PCT)) if int(clear.sum()) >= MIN_CLEAR else None
            out[(cls, str(st))] = {'n': int(m.sum()), 'n_clear': int(clear.sum()),
                                   'bound_sd': None if bound is None else round(bound, 4)}
    return out


def demote(Z, markers, labels, strata, bounds):
    demoted = np.zeros(len(labels), bool)
    why = np.full(len(labels), '', dtype=object)
    for cls in CT.OWN_MARKERS:
        v = own_marker_max(Z, markers, cls)
        in_class = labels == cls
        if v is None:
            why[in_class] = 'no own marker on this donor'
            continue
        for st in (pd.unique(strata[in_class]) if in_class.any() else []):
            m = in_class & (strata == st)
            rec = bounds.get((cls, str(st)))
            if rec is None or rec['bound_sd'] is None:
                why[m] = f"bound unavailable ({rec['n_clear'] if rec else 0} clear < {MIN_CLEAR})"
                continue
            below = m & (v < rec['bound_sd'])
            demoted |= below
            why[below] = f"below {cls} bound {rec['bound_sd']:+.3f} SD in stratum {st or 'none'}"
    return demoted, why


def route(labels, demoted, immune_leaf):
    new = labels.copy()
    exo = demoted & (labels == 'Exocrine')
    oth = demoted & (labels != 'Exocrine')
    new[exo] = UNRESOLVED
    new[oth & immune_leaf] = 'Immune'
    new[oth & ~immune_leaf] = 'Exocrine'
    return new


def immune_leaf(aligned: pd.DataFrame, scale: pd.DataFrame) -> np.ndarray:
    """Does each cell meet AT LEAST ONE immune leaf (Imm-/T-/Mac-) at scaled >= 0.10 (score.meets_requirements:
    NaN fails, an unmeasured term fails the rule)."""
    rules, _ = CT.codex_rules()
    have = [m for m in scale.index if m in aligned.columns]
    Xs = pd.DataFrame(CL.apply_scale(aligned[have].to_numpy(), have, scale), columns=have)
    hot, cold = Xs >= CT.IMMUNE_REQUIRE_MIN, Xs < CT.IMMUNE_REQUIRE_MIN
    ok = np.zeros(len(Xs), bool)
    for name, (all_of, any_of, none_of) in rules.items():
        if not name.startswith(('Imm-', 'T-', 'Mac-')):
            continue
        if any(m not in Xs for m in list(all_of) + list(none_of) + [x for g in any_of for x in g]):
            continue
        col = np.ones(len(Xs), bool)
        for m in all_of:
            col &= hot[m].to_numpy()
        for g in any_of:
            col &= hot[list(g)].to_numpy().any(axis=1)
        for m in none_of:
            col &= cold[m].to_numpy()
        ok |= col
    return ok


def refine_donor(out: Path, tag: str, params: pd.DataFrame, scale: pd.DataFrame, markers: list[str]):
    cells = load_cells(out, tag)
    if (cells['cell_class'] == UNRESOLVED).any():
        raise RefinementAlreadyApplied(tag)
    values, have = donor_values(out, tag, markers)
    frame = cells.merge(values, on='object_id', how='left', validate='one_to_one')
    aligned = pd.read_parquet(CT.typing_dir(out) / 'aligned' / f'{tag}.parquet')
    aligned = cells[['object_id']].merge(aligned, on='object_id', how='left', validate='one_to_one')
    measurable = (frame['qc_analysis_eligible'].fillna(False).astype(bool)
                  & ~frame['cell_class'].isin(UNTOUCHED)
                  & np.isfinite(frame[have].to_numpy('float64')).all(axis=1)).to_numpy()
    reason = np.full(len(frame), '', dtype=object)
    reason[~measurable] = 'not measured (excluded, ineligible or no aligned value)'
    for c in UNTOUCHED:
        reason[(frame['cell_class'] == c).to_numpy()] = f'{c}: outside the pass'
    Z = to_Z(frame.loc[measurable, have].to_numpy('float32'), have, params)
    labels = frame.loc[measurable, 'cell_class'].to_numpy().astype(object)
    strata = frame.loc[measurable, 'size_stratum'].to_numpy().astype(object)
    leaf = immune_leaf(aligned, scale)
    t0 = time.time()
    nbr = neighbours(Z)
    before = concordance(nbr, labels)
    bounds = fit_bounds(Z, have, labels, strata, before)
    demoted, why = demote(Z, have, labels, strata, bounds)
    new = route(labels, demoted, leaf[measurable])
    after = concordance(nbr, new)
    res = pd.DataFrame({'object_id': frame['object_id'].to_numpy(), 'cell_class': frame['cell_class'].to_numpy(),
                        'refined_class': frame['cell_class'].to_numpy().astype(object),
                        'size_stratum': frame['size_stratum'].to_numpy(), 'islet_cells': frame['islet_cells'],
                        'immune_leaf': leaf, 'concordance': np.nan, 'concordance_before': np.nan,
                        'demoted': False, 'refine_note': reason, 'refine_pass': np.where(measurable, 1, 0)})
    res.loc[measurable, 'refined_class'] = new
    res.loc[measurable, 'concordance'] = after
    res.loc[measurable, 'concordance_before'] = before
    res.loc[measurable, 'demoted'] = demoted
    res.loc[measurable, 'refine_note'] = why
    moved = pd.Series(list(zip(labels[demoted], new[demoted]))).value_counts()
    log = {'image': tag, 'markers': have, 'n_measurable': int(measurable.sum()), 'n_demoted': int(demoted.sum()),
           'seconds_graph': round(time.time() - t0),
           'bounds': {f'{c}|{s}': v for (c, s), v in bounds.items()},
           'demoted_by_class': {str(k): int(v) for k, v in pd.Series(labels[demoted]).value_counts().items()},
           'moves': {f'{a} -> {b}': int(c) for (a, b), c in moved.items()},
           'concordance_median_before': float(np.median(before)), 'concordance_median_after': float(np.median(after))}
    return res, log


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', default=['all'])
    a = ap.parse_args()
    typ = CT.typing_dir(a.out)
    all_tags = [r['image_name'] for r in CT.cohort_rows()]
    tags = all_tags if a.images == ['all'] else a.images
    scale = CL.load_scale(a.out)
    pfile = typ / 'refine' / 'representation.parquet'
    pfile.parent.mkdir(parents=True, exist_ok=True)
    if pfile.exists():
        params = pd.read_parquet(pfile)
    else:
        params = representation_params(a.out, all_tags, REFINE_MARKERS)
        params.to_parquet(pfile)
    markers = [m for m in REFINE_MARKERS if m in params.index]
    for t in tags:
        res, log = refine_donor(a.out, t, params, scale, markers)
        res.to_parquet(typ / 'refine' / f'{t}.parquet', index=False)
        (typ / 'refine' / f'{t}.json').write_text(json.dumps(log, indent=1))
        print(f'{t}: {log["n_demoted"]:,} of {log["n_measurable"]:,} demoted; {log["moves"]}', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
