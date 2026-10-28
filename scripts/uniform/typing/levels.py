#!/usr/bin/env python3
"""levels.py -- each marker's positivity level per image, chosen by its own reference (final typing plan, Step 2).

Owner, 2026-10-28 and Senior 10-19/10-26: positivity is binary, from each image's own background, on the values after
REDSEA, one measurement per marker for every cell (whole cell, membrane or nucleus), and "one rule doesn't fit all
situations"; a marker is accepted by its effect on object and cell-type positivity. The registry
(`registry/<panel>__<tissue>.csv`) gives each marker its candidate measurements and its positivity class; this script
computes every candidate level on every image, scores it against the class's reference, and keeps one per marker.

Candidate levels, on log1p of the REDSEA value (exact zeros set aside when reading a peak):
  peak_dense   the densest point of the smoothed histogram + K background widths (its own half-width, this image)
  peak_low     the lowest prominent peak + K widths (owner 10-19: "why would you take the densest point as the
               background peak?")
  tissue_p99   endocrine markers only: the marker's p99 among the image's cells outside every object (Senior's per-donor
               reference level, the earlier hormone level)
  edge / gate  markers most cells carry: the expressing population's lower edge (peak - K right widths) or a
               two-component mixture gate, either floored at the empty-slide null ("the background outside the
               tissue [helps] with the high expressers like ECAD")
  negative_pNN known-negative markers (registry `background_method`, owner decisions of 2026-10-28): the NN-th percentile
               (p99 = 99, p999 = 99.9) of the marker among reference cells that cannot carry it (Barlow et al. eLife
               2024: positive above the 99th percentile of cells known not to express the marker), per image. SMA:
               p99 of endocrine-object hormone-positive cells and leukocytes with no cell positive at SMA's background
               peak within ISOLATION_UM; its check is the Capillary-object enrichment. CD20: p99.9 of epithelial cells
               (E-cadherin or panCK positive, CD45 negative) farther than ISOLATION_UM from any leukocyte; its check is
               reported, not a gate: CD45+ among CD20+ cells, the B-cell calls at plain positivity (CD45+, no other
               hematopoietic lineage), and the expected false calls, 0.1 % of the CD45+ cells that are not B cells.
The marker sets the checks use come from the registry's lineages (hormones = the Endocrine subtypes, leukocytes = the
Leukocyte marker with a T or Myeloid marker, and so on); this module names lineages, never markers.

Checks (constants in `registry/<panel>__<tissue>__constants.csv`, the agent's determination delegated by the owner on
2026-10-28; each counts how often a marker is positive where it cannot be), in dependency order. Every out-of-place rate
is read on reference cells with no cell positive for the marker within ISOLATION_UM (centroids): measured on the 7 ND
images, most cross-lineage positives are a touching neighbor's signal (E-cadherin in leukocytes 17 % overall, 0.9 %
beyond 20 um), which no cutoff can remove and the canary step resolves; the isolated rate reads the cutoff's own errors.
Leukocytes = CD45+ and CD3e+, CD68+ or CD163+ (two markers).
  1 objects     structure and endocrine markers: positive >= OBJECT_ENRICHMENT_MIN times as often inside their objects
                as in ordinary tissue (cells in no object)
  2 immune      every leukocyte-only marker positive in <= CD45_AMONG_HORMONE_MAX of hormone-positive endocrine-object
                cells; CD45+ in >= CD45_AMONG_CD3E_MIN of CD3e+ cells; CD3e+ in <= CD3E_AMONG_MYELOID_MAX of CD68+/CD163+
                cells; among CD20+CD45+ cells CD3e+ and CD68+/CD163+ each <= B_PAIR_OTHER_LINEAGE_MAX
  3 epithelial  markers most cells carry: positive in <= EPITHELIAL_IN_IMMUNE_MAX of CD3e+ or CD68+/CD163+ cells
  4 shared      positive in <= SHARED_OUT_OF_PLACE_MAX of cells of a compartment that cannot carry them
A rate is computed only when its denominator holds MIN_CELLS_FOR_RATE cells ("not evaluable" otherwise, reported).
Selection: on the 7 ND images (owner, Senior 10-01: test classifiers on the non-diabetic controls, because disease
changes the reference populations, e.g. insulitis puts leukocytes among endocrine cells), every check must pass; among
the passing options the highest sensitivity is kept. The same checks are computed and reported for the Aab+ and T1D
images; a failure there is flagged for review. A chosen background-based level more than 10x from its cohort median
stops the run (a known-negative level is judged by its own check; background_based).
A marker with no passing option is left out of typing on all images, unless the registry marks it required (a class
cannot do without it): then the run stops.

Writes typing/levels.csv (chosen level per image and marker), typing/levels_options.csv (every option, image and
check), typing/levels_selection.csv (per marker: chosen option, or dropped / stop, and why).

    python scripts/uniform/typing/levels.py --out $OUT --workers 20 [--images ND_0112 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product, repeat
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
import registry as RG  # noqa: E402

EDGES = np.linspace(0.0, 14.0, 281)            # log1p units, 0.05 wide
CEN = (EDGES[:-1] + EDGES[1:]) / 2
SMOOTH_BINS = 2
PEAK_HEIGHT, PEAK_PROMINENCE = 0.10, 0.05      # of the histogram maximum
HWHM_TO_SIGMA = 1.1774
V17 = 'superseded/typing-v17-20261028'          # holds the empty-slide nulls (raw pixel window p99.9) of v16/v17


# -- level methods ----------------------------------------------------------------------------------------------------
def hist(v: np.ndarray) -> np.ndarray:
    return gaussian_filter1d(np.histogram(v, bins=EDGES)[0].astype(float), SMOOTH_BINS)


def half_width(h: np.ndarray, i: int, side: str) -> float:
    half, j = h[i] / 2, i
    step = -1 if side == 'left' else 1
    while 0 < j < len(h) - 1 and h[j] > half:
        j += step
    if h[j] > half:
        return np.nan
    return abs(CEN[i] - CEN[j]) / HWHM_TO_SIGMA


MAX_WIDTH = 1.5                                 # log units: a wider "peak" is a plateau, not a background peak
AGAINST_ZERO = float(np.log1p(5.0))            # a peak this low may sit against zero; its width is read on the right


def peak_level(v: np.ndarray, k: float, lowest: bool) -> tuple[float, float, float]:
    """(center, width, level) of a background peak. The width is the left half-width (the right one only for a peak
    against zero); a width that cannot be read or exceeds MAX_WIDTH makes the level unusable (NaN), so the option fails
    on that image instead of producing a level far from the background."""
    if len(v) < 1000:
        return np.nan, np.nan, np.nan
    h = hist(v)
    mx = h.max()
    i = int(h.argmax())
    if lowest:
        pk, _ = find_peaks(h, height=PEAK_HEIGHT * mx, prominence=PEAK_PROMINENCE * mx)
        if len(pk):
            i = int(pk.min())
    w = half_width(h, i, 'left')
    if not np.isfinite(w) and CEN[i] <= AGAINST_ZERO:
        w = half_width(h, i, 'right')
    if not np.isfinite(w) or w > MAX_WIDTH:
        return float(CEN[i]), float(w), np.nan
    return float(CEN[i]), float(w), float(CEN[i] + k * w)


def edge_level(v: np.ndarray, k: float) -> tuple[float, float, float]:
    """The expressing population's lower edge: its peak minus k right-side widths."""
    if len(v) < 1000:
        return np.nan, np.nan, np.nan
    h = hist(v)
    i = int(h.argmax())
    w = half_width(h, i, 'right')
    if not np.isfinite(w) or w > MAX_WIDTH:
        return float(CEN[i]), float(w), np.nan
    return float(CEN[i]), float(w), float(CEN[i] - k * w)


def gate_level(v: np.ndarray, seed: int = 0, n: int = 200000) -> tuple[float, float, float]:
    if len(v) < 1000:
        return np.nan, np.nan, np.nan
    if len(v) > n:
        v = np.random.default_rng(seed).choice(v, n, replace=False)
    gm = GaussianMixture(n_components=2, random_state=seed).fit(v.reshape(-1, 1))
    return float(np.min(gm.means_)), float('nan'), float(np.mean(gm.means_))


def negative_quantile(method: str) -> float | None:
    """negative_p99 -> 99.0, negative_p999 -> 99.9; None for every other method."""
    if not method.startswith('negative_p'):
        return None
    digits = method[len('negative_p'):]
    return 100.0 - 10.0 ** (2 - len(digits))


def roles(reg: RG.Registry) -> dict:
    """The marker sets the checks need, from the registry's lineages (owner 2026-10-28: no marker lists in code)."""
    leukocyte, t, myeloid = reg.broad('Leukocyte'), reg.broad('T'), reg.broad('Myeloid')
    core = t + leukocyte + myeloid                     # chosen together against canonical co-expression
    epithelial = [m for m in reg.defining('Epithelial') if reg.positivity_class(m) == 'majority_epithelial']
    return {'hormones': list(reg.subtypes('Endocrine')), 'leukocyte': leukocyte, 't': t, 'myeloid': myeloid,
            'core': core, 'epithelial': epithelial,
            'immune2': [m for m in reg.measured_markers()
                        if reg.positivity_class(m) == 'immune_coexpression' and m not in core]}


def any_of(fixed: dict, markers: list[str], n: int) -> np.ndarray:
    out = np.zeros(n, bool)
    for m in markers:
        out |= fixed[m]
    return out


def methods_for(reg: RG.Registry, m: str) -> list[str]:
    cls = reg.positivity_class(m)
    if cls == 'known_negative':
        return [reg.background_method(m)]
    if cls == 'majority_epithelial':
        return ['edge', 'gate']
    if cls == 'endocrine_object':
        return ['peak_dense', 'peak_low', 'tissue_p99']
    return ['peak_dense', 'peak_low']


def options(reg: RG.Registry, m: str) -> list[tuple[str, str]]:
    return list(product(reg.candidates(m), methods_for(reg, m)))


# -- per image --------------------------------------------------------------------------------------------------------
def empty_nulls(out: Path) -> dict:
    """Raw empty-slide nulls of the markers most cells carry, recovered from v16/v17's levels file
    (center = log1p(null x UniFORM factor))."""
    lv = pd.read_csv(out / V17 / 'positivity_levels.csv')
    shifts = json.loads((out / V17 / 'shifts.json').read_text())['shifts']
    nulls = {}
    for r in lv[lv['source'].isin(['epithelial_edge', 'mixture_gate', 'empty_area'])].itertuples():
        f = shifts.get(r.marker, {}).get('scale', {}).get(r.image)
        f = f if isinstance(f, float) and f > 0 else 1.0
        nulls[(r.image, r.marker)] = float(np.expm1(r.center) / f) if np.isfinite(r.center) else np.nan
    return nulls


def load(out: Path, tag: str, reg: RG.Registry) -> pd.DataFrame:
    cols = sorted({reg.column(m, meas) for m in reg.measured_markers() for meas in reg.candidates(m)})
    rs = pd.read_parquet(out / 'typing' / 'redsea' / f'{tag}.parquet', columns=['object_id', 'qc_analysis_eligible'] + cols)
    ctx = pd.read_parquet(out / 'typing' / 'context' / f'{tag}.parquet',
                          columns=['object_id', 'context', 'endocrine_object_class'])
    xy = pd.read_parquet(out / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    xy.columns = ['object_id', 'x_um', 'y_um']
    d = rs[rs['qc_analysis_eligible']].merge(ctx, on='object_id', how='left').merge(xy, on='object_id', how='left')
    d['in_endocrine'] = d['endocrine_object_class'].notna()
    d['ordinary'] = d['context'].eq('tissue') & ~d['in_endocrine']
    return d.reset_index(drop=True)


def level_for(d: pd.DataFrame, reg: RG.Registry, m: str, meas: str, method: str, k: float, null: float,
              ref: np.ndarray | None = None) -> dict:
    """One option's level on one image; `ref` is the reference-cell mask of a known-negative method."""
    x = d[reg.column(m, meas)].to_numpy(float)
    v = np.log1p(x[np.isfinite(x) & (x > 0)])
    q = negative_quantile(method)
    if method in ('peak_dense', 'peak_low'):
        c, w, lvl = peak_level(v, k, lowest=(method == 'peak_low'))
    elif q is not None:
        o = x[ref & np.isfinite(x)] if ref is not None else np.array([])
        c, w = np.nan, np.nan
        enough = len(o) * (1 - q / 100) >= reg.const('KNOWN_NEGATIVE_MIN_ABOVE')
        lvl = float(np.log1p(np.percentile(o, q))) if enough else np.nan
    elif method == 'tissue_p99':
        o = x[d['ordinary'].to_numpy() & np.isfinite(x)]
        c, w, lvl = np.nan, np.nan, float(np.log1p(np.percentile(o, 99))) if len(o) >= 1000 else np.nan
    else:
        c, w, own = edge_level(v, k) if method == 'edge' else gate_level(v)
        floor = np.log1p(null) if np.isfinite(null) else np.nan
        lvl = float(np.nanmax([own, floor])) if np.isfinite(own) or np.isfinite(floor) else np.nan
    pos = np.log1p(np.clip(np.nan_to_num(x, nan=0.0), 0, None)) > lvl if np.isfinite(lvl) else np.zeros(len(x), bool)
    return {'center': c, 'width': w, 'level_log': lvl, 'level_raw': float(np.expm1(lvl)) if np.isfinite(lvl) else np.nan,
            'positive_fraction': float(pos.mean()), 'zero_fraction': float(np.mean(x <= 0)), 'pos': pos}


def rate(num: np.ndarray, den: np.ndarray, nmin: int) -> float:
    return float(num[den].mean()) if den.sum() >= nmin else np.nan


def isolated_rate(p: np.ndarray, ref: np.ndarray, xy: np.ndarray, cut: float, nmin: int) -> tuple[float, float, int]:
    """(rate in all reference cells, rate in reference cells with no cell positive for the marker outside the reference
    within `cut` um, number of those isolated cells). The isolated rate reads a cutoff's own errors; the difference is
    a touching neighbor's signal, which the canary step resolves."""
    all_rate = rate(p, ref, nmin)
    src = np.flatnonzero(p & ~ref)
    idx = np.flatnonzero(ref)
    if len(idx) == 0:
        return all_rate, np.nan, 0
    if len(src) == 0:
        far = np.ones(len(idx), bool)
    else:
        dist, _ = cKDTree(xy[src]).query(xy[idx], k=1)
        far = dist > cut
    iso = np.zeros(len(p), bool)
    iso[idx[far]] = True
    return all_rate, rate(p, iso, nmin), int(far.sum())


STRUCTURE_CONTEXT = {'Capillary objects': 'Blood vessel', 'Nerve objects': 'Nerve', 'Lymphatic objects': 'Lymphatic'}
# the reference populations a known-negative level is read from (registry `reference` of the known_negative class)
HORMONES_AND_LEUKOCYTES = 'endocrine-object hormone-positive cells; leukocytes'
EPITHELIAL_AWAY_FROM_LEUKOCYTES = 'epithelial cells away from leukocytes'


def isolated_from(src: np.ndarray, ref: np.ndarray, xy: np.ndarray, cut: float) -> np.ndarray:
    """The cells of `ref` with no `src` cell within `cut` um (centroids)."""
    out = np.zeros(len(ref), bool)
    idx = np.flatnonzero(ref)
    s = np.flatnonzero(src)
    if len(idx) == 0:
        return out
    if len(s) == 0:
        out[idx] = True
        return out
    dist, _ = cKDTree(xy[s]).query(xy[idx], k=1)
    out[idx[dist > cut]] = True
    return out


def cell_sets(d: pd.DataFrame, fixed: dict, R: dict) -> tuple[np.ndarray, np.ndarray | None]:
    """(endocrine-object hormone-positive cells, leukocytes) from the chosen positivity; leukocytes need the Leukocyte
    marker and a T or Myeloid marker (two markers), None until those are chosen."""
    n = len(d)
    hormone = any_of(fixed, [h for h in R['hormones'] if h in fixed], n) & d['in_endocrine'].to_numpy()
    if not all(m in fixed for m in R['core']):
        return hormone, None
    return hormone, any_of(fixed, R['leukocyte'], n) & any_of(fixed, R['t'] + R['myeloid'], n)


def negative_reference(d: pd.DataFrame, reg: RG.Registry, m: str, meas: str, fixed: dict, xy: np.ndarray, k: float,
                       cut: float, R: dict) -> np.ndarray:
    """The isolated reference cells whose percentile is a known-negative marker's level (registry `reference`):
    endocrine-object hormone-positive cells and leukocytes with no cell positive for the marker at its background peak
    within `cut` (isolation as in isolated_rate), or epithelial cells (an Epithelial majority marker positive, the
    Leukocyte marker negative) with no leukocyte within `cut`."""
    name = reg.table.loc[m, 'reference']
    hormone, leuk = cell_sets(d, fixed, R)
    if leuk is None:
        raise SystemExit(f'{m}: the leukocyte markers must be chosen before a known-negative level')
    if name == HORMONES_AND_LEUKOCYTES:
        ref = hormone | leuk
        src = level_for(d, reg, m, meas, 'peak_dense', k, np.nan)['pos'] & ~ref
    elif name == EPITHELIAL_AWAY_FROM_LEUKOCYTES:
        ref = any_of(fixed, R['epithelial'], len(d)) & ~any_of(fixed, R['leukocyte'], len(d))
        src = leuk
    else:
        raise SystemExit(f'{m}: unknown known-negative reference {name!r}')
    return isolated_from(src, ref, xy, cut)


def chosen_positivity(d: pd.DataFrame, reg: RG.Registry, chosen: dict, k: float, nulls: dict, tag: str,
                      xy: np.ndarray, cut: float, R: dict) -> tuple[dict, dict]:
    """Each chosen marker's level and positivity on this image; known-negative markers last (their reference cells
    are defined by the others)."""
    fixed, lvs = {}, {}
    for m in sorted(chosen, key=lambda x: reg.positivity_class(x) == 'known_negative'):
        meas, method = chosen[m]
        ref = negative_reference(d, reg, m, meas, fixed, xy, k, cut, R) if negative_quantile(method) else None
        lvs[m] = level_for(d, reg, m, meas, method, k, nulls.get((tag, m), np.nan), ref)
        fixed[m] = lvs[m]['pos']
    return fixed, lvs


def lineage_report(reg: RG.Registry, m: str, p: np.ndarray, fixed: dict, q: float, nmin: int) -> dict:
    """A reported known-negative marker whose lineage needs other markers (registry `requires`; CD20: CD45 and one
    hematopoietic lineage): positives carrying the required markers, calls at plain positivity, and the expected false
    calls -- by construction (1 - q) of the required-positive cells that are not calls exceed the level."""
    n = len(p)
    need = reg.requires(m)
    req = [r for r in need if r in reg.table.index]
    required = np.ones(n, bool)
    for r in req:
        required &= fixed[r]
    comp = reg.table.loc[m, 'compartment']
    skip = {lin for x in [m] + req for lin, _ in reg.lineage_entries(x)}
    lineages = dict.fromkeys(lin for x in reg.defining(comp) for lin, _ in reg.lineage_entries(x))
    others = [x for lin in lineages if lin not in skip for x in reg.broad(lin)]
    other = any_of(fixed, others, n)
    calls = p & required & ~other if 'single_lineage' in need else p & required
    false = (1 - q / 100) * float((required & ~calls).sum())
    return {'positive_cells': int(p.sum()), 'required_among_positive': rate(required, p, nmin),
            'other_lineage_among_positive_required': rate(other, p & required, nmin),
            'calls_plain_level': int(calls.sum()), 'expected_false_calls': false,
            'expected_false_proportion': false / calls.sum() if calls.sum() else np.nan}


def image_pass(out: Path, tag: str, reg_name: tuple, stage: str, chosen: dict, nulls: dict) -> list[dict]:
    """Every option of the stage's markers on one image, with its checks; `chosen` holds earlier stages' choices.
    Out-of-place rates are read on isolated reference cells (ISOLATION_UM); the all-cell rate is reported beside."""
    reg = RG.load(*reg_name)
    R = roles(reg)
    k = reg.const('POSITIVITY_K')
    nmin = int(reg.const('MIN_CELLS_FOR_RATE'))
    cut = reg.const('ISOLATION_UM')
    d = load(out, tag, reg)
    n = len(d)
    xy = d[['x_um', 'y_um']].to_numpy(float)
    fixed, chosen_levels = chosen_positivity(d, reg, chosen, k, nulls, tag, xy, cut, R)
    inside_endo = d['in_endocrine'].to_numpy()
    ordinary = d['ordinary'].to_numpy()
    hormone_cells, leuk = cell_sets(d, fixed, R)
    rows = []

    def row(m, meas, method, lv, **checks):
        r = {'image': tag, 'stage': stage, 'marker': m, 'measurement': meas, 'method': method,
             **{k2: v for k2, v in lv.items() if k2 != 'pos'}, **checks}
        if not np.isfinite(lv['level_log']):          # no level on this image: the option is unusable here
            r['passes'] = False
        rows.append(r)

    def le(v, bar):
        return (v <= bar) if np.isfinite(v) else np.nan

    def enrichment(p, inside):
        ri, ro = rate(p, inside, nmin), rate(p, ordinary, nmin)
        if not (np.isfinite(ri) and np.isfinite(ro)):
            enr = np.nan
        else:
            enr = ri / ro if ro > 0 else np.inf
        return ri, ro, enr, ((enr >= reg.const('OBJECT_ENRICHMENT_MIN')) if not np.isnan(enr) else np.nan)

    if stage == 'objects':
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) in ('structure_object', 'endocrine_object')]:
            ref = reg.table.loc[m, 'reference']
            inside = inside_endo if reg.positivity_class(m) == 'endocrine_object' else \
                d['context'].eq(STRUCTURE_CONTEXT[ref]).to_numpy()
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, np.nan)
                ri, ro, enr, ok = enrichment(lv['pos'], inside)
                row(m, meas, method, lv, sensitivity=ri, out_of_place=ro, enrichment=enr, passes=ok)
    elif stage == 'immune':
        core = R['core']
        lv = {m: {opt: level_for(d, reg, m, *opt, k, np.nan) for opt in options(reg, m)} for m in core}
        for combo in product(*(lv[m].keys() for m in core)):
            pos = {m: lv[m][o]['pos'] for m, o in zip(core, combo)}
            t, l, mye = any_of(pos, R['t'], n), any_of(pos, R['leukocyte'], n), any_of(pos, R['myeloid'], n)
            c45 = rate(l, t, nmin)
            t_all, t_iso, _ = isolated_rate(t, mye, xy, cut, nmin)
            checks = {'leukocyte_among_t': c45, 't_among_myeloid_all': t_all, 't_among_myeloid': t_iso}
            for m in core:
                r_all, r_iso, _ = isolated_rate(pos[m], hormone_cells, xy, cut, nmin)
                checks[f'{m}_in_hormone_all'], checks[f'{m}_in_hormone'] = r_all, r_iso
            ok = [(c45 >= reg.const('CD45_AMONG_CD3E_MIN')) if np.isfinite(c45) else True,
                  le(t_iso, reg.const('CD3E_AMONG_MYELOID_MAX'))] + \
                 [le(checks[f'{m}_in_hormone'], reg.const('CD45_AMONG_HORMONE_MAX')) for m in core]
            ok = [True if (isinstance(x, float) and np.isnan(x)) else bool(x) for x in ok]
            ok.append(all(np.isfinite(lv[m][o]['level_log']) for m, o in zip(core, combo)))
            oop = np.nanmax([t_iso] + [checks[f'{m}_in_hormone'] for m in core])
            checks.update({f'level_raw_{m}': lv[m][o]['level_raw'] for m, o in zip(core, combo)})
            rows.append({'image': tag, 'stage': stage, 'marker': '+'.join(core),
                         'measurement': '|'.join(o[0] for o in combo), 'method': '|'.join(o[1] for o in combo),
                         'out_of_place': oop, 'sensitivity': c45, 'passes': all(ok), **checks})
    elif stage == 'immune2':
        t = any_of(fixed, R['t'], n)
        for m in R['immune2']:
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, np.nan)
                p = lv['pos']
                inh_all, inh, _ = isolated_rate(p, hormone_cells, xy, cut, nmin)
                ok = le(inh, reg.const('CD45_AMONG_HORMONE_MAX'))
                ok = True if (isinstance(ok, float) and np.isnan(ok)) else bool(ok)
                row(m, meas, method, lv, out_of_place=inh, sensitivity=rate(p, t, nmin), passes=ok,
                    in_hormone_all=inh_all, in_hormone=inh)
    elif stage == 'epithelial':
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) == 'majority_epithelial']:
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, nulls.get((tag, m), np.nan))
                r_all, r_iso, n_iso = isolated_rate(lv['pos'], leuk, xy, cut, nmin)
                row(m, meas, method, lv, out_of_place=r_iso, out_of_place_all=r_all, isolated_cells=n_iso,
                    sensitivity=rate(lv['pos'], ordinary, nmin), passes=le(r_iso, reg.const('EPITHELIAL_IN_IMMUNE_MAX')))
    elif stage == 'known_negative':
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) == 'known_negative']:
            for meas, method in options(reg, m):
                ref = negative_reference(d, reg, m, meas, fixed, xy, k, cut, R)
                lv = level_for(d, reg, m, meas, method, k, np.nan, ref)
                p = lv['pos']
                inh_all, inh, _ = isolated_rate(p, hormone_cells, xy, cut, nmin)
                extra = {'reference_cells': int(ref.sum()), 'in_hormone_all': inh_all, 'in_hormone': inh}
                check = reg.validation(m)
                if check in STRUCTURE_CONTEXT:
                    ri, ro, enr, ok = enrichment(p, d['context'].eq(STRUCTURE_CONTEXT[check]).to_numpy())
                    row(m, meas, method, lv, sensitivity=ri, out_of_place=ro, enrichment=enr, passes=ok, **extra)
                elif check == 'reported':
                    extra.update(lineage_report(reg, m, p, fixed, negative_quantile(method), nmin))
                    row(m, meas, method, lv, passes=True, **extra)
                else:
                    raise SystemExit(f'{m}: unknown validation {check!r}')
    elif stage == 'shared':
        t_cells = any_of(fixed, R['t'], n) & any_of(fixed, R['leukocyte'], n)
        macro = any_of(fixed, R['myeloid'], n) & any_of(fixed, R['leukocyte'], n)
        refs = {'endocrine-object hormone-positive cells': hormone_cells, 'T cells': t_cells,
                'T cells and macrophages': t_cells | macro,
                'hematopoietic cells; endocrine-object hormone-positive cells': leuk | hormone_cells,
                HORMONES_AND_LEUKOCYTES: leuk | hormone_cells}
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) == 'shared_negative']:
            ref = refs[reg.table.loc[m, 'reference']]
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, np.nan)
                r_all, r_iso, n_iso = isolated_rate(lv['pos'], ref, xy, cut, nmin)
                row(m, meas, method, lv, out_of_place=r_iso, out_of_place_all=r_all, isolated_cells=n_iso,
                    passes=le(r_iso, reg.const('SHARED_OUT_OF_PLACE_MAX')))
    elif stage == 'state':
        for m in reg.markers('state'):
            meas = reg.measurement(m)
            lv = level_for(d, reg, m, meas, 'peak_dense', k, np.nan)
            row(m, meas, 'peak_dense', lv, passes=True)
    elif stage == 'final':                              # the chosen level of every marker kept
        for m, (meas, method) in chosen.items():
            x = d[reg.column(m, meas)].to_numpy(float)
            bc, bw, _ = peak_level(np.log1p(x[np.isfinite(x) & (x > 0)]), k, lowest=False)
            row(m, meas, method, chosen_levels[m], passes=True, background_center=bc, background_width=bw)
    return rows


# -- selection ----------------------------------------------------------------------------------------------------------
def background_based(method: pd.Series) -> np.ndarray:
    """Options whose level is read from the image's background (peak, edge, gate, tissue p99), where a level far from
    the cohort's is a misread background. Measured 2026-10-28: SMA's background peak is 0.7-2.3x its cohort median on
    every image while its known-negative p99 is 2.5-172x that peak (the reference cells' SMA differs between images),
    and the level passes its Capillary-object check on every ND image; the 10x rule would refuse ND_6356 and Aab_6090."""
    return ~method.astype(str).str.split('|').apply(lambda ms: any(negative_quantile(m) for m in ms)).to_numpy()


def select(df: pd.DataFrame, select_on: list[str]) -> pd.DataFrame:
    """Per marker (or marker combination): options passing in every selection image (the ND controls; owner, Senior
    10-01: "use all non-diabetic/controls for testing classifier since the endocrine differences with disease
    interfere"); among them the highest median sensitivity, then the lowest median out-of-place rate (every passing
    option already meets the specificity bar, so the one catching most true cells is kept). Failures on the other
    images are listed for review."""
    df = df.copy()
    for c in ('sensitivity', 'out_of_place'):
        if c not in df.columns:
            df[c] = np.nan
    # a background-based level more than 10x from the same option's cohort median fails that image (a wild level is not
    # a background); a known-negative level follows its reference cells by design and is judged by its own check
    key = ['marker', 'measurement', 'method']
    lcols = [c for c in df.columns if c == 'level_raw' or c.startswith('level_raw_')]
    wild = np.zeros(len(df), bool)
    for c in lcols:
        med = df.groupby(key)[c].transform('median')
        fold = df[c] / med
        wild |= ((fold > 10) | (fold < 0.1)).fillna(False).to_numpy()
    wild &= background_based(df['method'])
    df['wild_level'] = wild
    df.loc[wild, 'passes'] = False
    sel = df[df['image'].isin(select_on)]
    g = sel.assign(fail=sel['passes'].eq(False)).groupby(['marker', 'measurement', 'method'])
    s = g.agg(images=('image', 'nunique'), failing_images=('fail', 'sum'),
              out_of_place=('out_of_place', 'median'), sensitivity=('sensitivity', 'median')).reset_index()
    s['failing'] = g.apply(lambda x: ' '.join(sorted(x.loc[x['passes'].eq(False), 'image'])), include_groups=False).values
    other = df[~df['image'].isin(select_on)]
    flags = other[other['passes'].eq(False)].groupby(['marker', 'measurement', 'method'])['image'].apply(
        lambda x: ' '.join(sorted(x))).rename('flagged_other_images').reset_index()
    s = s.merge(flags, on=['marker', 'measurement', 'method'], how='left')
    s['ok'] = (s['images'] == len(select_on)) & (s['failing_images'] == 0)
    s['sens_key'] = s['sensitivity'].fillna(-1)
    s = s.sort_values(['marker', 'ok', 'sens_key', 'out_of_place'], ascending=[True, False, False, True])
    return s.drop(columns='sens_key')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='*', default=None, help='default: every image of cohort.csv')
    ap.add_argument('--workers', type=int, default=20)
    ap.add_argument('--report', type=Path, default=None, help='write tables here instead of $OUT/typing (trials)')
    a = ap.parse_args()
    out = a.out.resolve()
    reg = RG.load()
    reg_name = tuple(reg.name.split('__'))
    tags = a.images or [r['image_name'] for r in RG.csv.DictReader(RG.COHORT.open())]
    dest = (a.report or out / 'typing').resolve()
    dest.mkdir(parents=True, exist_ok=True)
    nulls = empty_nulls(out)
    status = {r['image_name']: r['donor_status'] for r in RG.csv.DictReader(RG.COHORT.open())}
    select_on = [t for t in tags if status.get(t) == 'ND']
    if not select_on:
        raise SystemExit('no ND image among the images: options are chosen on the ND controls')
    chosen: dict[str, tuple[str, str]] = {}
    decisions, all_rows = [], []
    for stage in ('objects', 'immune', 'immune2', 'epithelial', 'known_negative', 'shared', 'state'):
        with ProcessPoolExecutor(a.workers) as ex:
            rows = [r for rs in ex.map(image_pass, repeat(out), tags, repeat(reg_name), repeat(stage), repeat(dict(chosen)),
                                       repeat(nulls)) for r in rs]
        df = pd.DataFrame(rows)
        all_rows.append(df)
        sel = select(df, select_on)
        for marker, opts in sel.groupby('marker', sort=False):
            best = opts.iloc[0]
            names = marker.split('+') if stage == 'immune' else [marker]
            meas = best['measurement'].split('|') if stage == 'immune' else [best['measurement']]
            meth = best['method'].split('|') if stage == 'immune' else [best['method']]
            for i, m in enumerate(names):
                chosen[m] = (meas[i], meth[i])      # the best option, also when it fails: later stages and the report
                flagged = best.get('flagged_other_images')
                flagged = flagged if isinstance(flagged, str) else ''
                if best['ok']:                      # need it; the decision below says whether it may be used
                    decisions.append({'marker': m, 'decision': 'use', 'measurement': meas[i], 'method': meth[i],
                                      'out_of_place': best['out_of_place'], 'sensitivity': best['sensitivity'],
                                      'why': '', 'flagged_other_images': flagged})
                else:
                    required = reg.table.loc[m, 'required'] == 'yes'
                    decisions.append({'marker': m, 'decision': 'stop' if required else 'drop',
                                      'measurement': meas[i], 'method': meth[i], 'out_of_place': best['out_of_place'],
                                      'sensitivity': best['sensitivity'],
                                      'why': f"no option passes on every ND image; best fails on: {best['failing']}",
                                      'flagged_other_images': flagged})
        print(f'[{stage}]', {m: chosen.get(m) for m in sel['marker'].unique()}, flush=True)
    opts = pd.concat(all_rows, ignore_index=True)
    opts.to_csv(dest / 'levels_options.csv', index=False)
    dec = pd.DataFrame(decisions)
    dec.to_csv(dest / 'levels_selection.csv', index=False)
    with ProcessPoolExecutor(a.workers) as ex:
        final = [r for rs in ex.map(image_pass, repeat(out), tags, repeat(reg_name), repeat('final'), repeat(dict(chosen)),
                                    repeat(nulls)) for r in rs]
    cols = ['image', 'marker', 'measurement', 'method', 'center', 'width', 'level_log', 'level_raw', 'positive_fraction',
            'zero_fraction', 'background_center', 'background_width']
    lv = pd.DataFrame(final)[cols]
    med = lv.groupby('marker')['level_raw'].transform('median')
    lv['fold_from_cohort_median'] = lv['level_raw'] / med
    lv.to_csv(dest / 'levels.csv', index=False)
    used = set(dec.loc[dec['decision'] == 'use', 'marker'])
    wild = lv[lv['marker'].isin(used) & ((lv['fold_from_cohort_median'] > 10) | (lv['fold_from_cohort_median'] < 0.1))
              & background_based(lv['method'])]
    pd.set_option('display.width', 250)
    print(dec.to_string(index=False))
    if len(wild):
        print('LEVELS MORE THAN 10x FROM THEIR COHORT MEDIAN:')
        print(wild.to_string(index=False))
    stops = dec[dec['decision'] == 'stop']
    if len(wild):
        return 3
    if len(stops):
        print(f'STOP: defining markers without a passing option: {", ".join(stops["marker"])}', flush=True)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
