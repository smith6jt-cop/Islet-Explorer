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
               peak within ISOLATION_UM; since the method repair its gate is exocrine epithelium and the Capillary-object
               enrichment is reported (owner decision pending: the owner's decision named the Capillary check). CD20: p99.9 of epithelial cells
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
import identity as ID  # noqa: E402

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


def peak_level(v: np.ndarray, k: float, lowest: bool, side: str = 'left') -> tuple[float, float, float]:
    """(center, width, level) of a background peak. The width is the left half-width (the right one only for a peak
    against zero); side 'narrow' takes the narrower of the two half-widths (REDSEA's exact zeros spread into the left
    side: method repair 2026-10-28, CD68/CD163 and the isolation of known-negative references), side 'fwhm' their mean
    (CD45). A width that cannot be read or exceeds MAX_WIDTH makes the level unusable (NaN), so the option fails on that
    image instead of producing a level far from the background."""
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
    if side in ('narrow', 'fwhm'):
        r = half_width(h, i, 'right')
        both = [x for x in (w, r) if np.isfinite(x)]
        if both:
            w = min(both) if side == 'narrow' else (float(np.mean(both)) if len(both) == 2 else np.nan)
    elif not np.isfinite(w) and CEN[i] <= AGAINST_ZERO:
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
    return {'hormones': list(reg.subtypes('Endocrine')), 'subtype_of': reg.subtypes('Endocrine'),
            'leukocyte': leukocyte, 't': t, 'myeloid': myeloid,
            'core': core, 'epithelial': epithelial,
            'ductal': [m for m in reg.markers('defining') if {lin for lin, _ in reg.lineage_entries(m)} == {'Exocrine'}],
            'immune2': [m for m in reg.measured_markers()
                        if reg.positivity_class(m) == 'immune_coexpression' and m not in core]}


def any_of(fixed: dict, markers: list[str], n: int) -> np.ndarray:
    out = np.zeros(n, bool)
    for m in markers:
        out |= fixed[m]
    return out


def methods_for(reg: RG.Registry, m: str) -> list[str]:
    """The level methods a marker's options are chosen among; one, where the registry fixes it (background_method)."""
    if reg.background_method(m):
        return [reg.background_method(m)]
    cls = reg.positivity_class(m)
    if cls == 'majority_epithelial':
        return ['edge', 'gate']
    if cls == 'endocrine_object':
        return ['peak_dense', 'peak_low', 'tissue_p99']
    return ['peak_dense', 'peak_low']


def options(reg: RG.Registry, m: str) -> list[tuple[str, str]]:
    return list(product(reg.candidates(m), methods_for(reg, m)))


# -- per image --------------------------------------------------------------------------------------------------------
def empty_nulls(out: Path) -> dict:
    """Empty-slide nulls of the markers most cells carry: the EMPTY_NULL_QUANTILE-th percentile of 7.2-um window means
    at least 20 um outside the Tissue outline (raw image units; typing/empty_area_nulls.csv, written from
    scripts/research/empty_area_background_20261028.py). Method repair 2026-10-28: the 99.9th percentile was set by single
    debris or artifact objects off the tissue on ND_6539, ND_0112 and ND_6356 and floored those images' levels."""
    q = RG.load().const('EMPTY_NULL_QUANTILE')
    t = pd.read_csv(out / 'typing' / 'empty_area_nulls.csv')
    col = f'null_p{q:g}'.replace('.', '')
    return {(r.image, r.marker): float(getattr(r, col)) for r in t.itertuples()}


def load(out: Path, tag: str, reg: RG.Registry) -> pd.DataFrame:
    cols = sorted({reg.column(m, meas) for m in reg.measured_markers() for meas in reg.candidates(m)})
    rs = pd.read_parquet(out / 'typing' / 'redsea' / f'{tag}.parquet', columns=['object_id', 'qc_analysis_eligible'] + cols)
    ctx = pd.read_parquet(out / 'typing' / 'context' / f'{tag}.parquet',
                          columns=['object_id', 'context', 'endocrine_object_class'] + list(OVERLAP.values()))
    xy = pd.read_parquet(out / 'cells' / f'{tag}.parquet', columns=['Object ID', 'Centroid X µm', 'Centroid Y µm'])
    xy.columns = ['object_id', 'x_um', 'y_um']
    d = rs[rs['qc_analysis_eligible']].merge(ctx, on='object_id', how='left').merge(xy, on='object_id', how='left')
    d['in_endocrine'] = d['endocrine_object_class'].notna()
    d['ordinary'] = d['context'].eq('tissue') & ~d['in_endocrine']
    # the established endocrine typing (owner 2026-10-28: "The endocrine typing didn't need to change"): the verified
    # endocrine cells are the reference wherever a check needs endocrine cells
    est = out / ID.ESTABLISHED / 'calls' / f'{tag}.parquet'
    if est.exists():
        c = pd.read_parquet(est, columns=['object_id', 'cell_class']).rename(columns={'cell_class': 'established_class'})
        d = d.merge(c, on='object_id', how='left', validate='one_to_one')
    d['established_class'] = d.get('established_class', pd.Series('', index=d.index)).fillna('')
    return d.reset_index(drop=True)


PEAK_SIDE = {'peak_dense': 'left', 'peak_low': 'left', 'peak_narrow': 'narrow', 'peak_fwhm': 'fwhm'}


def level_for(d: pd.DataFrame, reg: RG.Registry, m: str, meas: str, method: str, k: float, null: float,
              ref=None, subset: np.ndarray | None = None) -> dict:
    """One option's level on one image. `ref`: the reference-cell mask of a known-negative method, or a list of masks
    whose percentiles are each read and the highest kept (level_combine 'each': every lineage held to its own
    percentile); `subset`: the cells a background peak is read on (registry background_population)."""
    x = d[reg.column(m, meas)].to_numpy(float)
    xs = x if subset is None else x[subset]
    v = np.log1p(xs[np.isfinite(xs) & (xs > 0)])
    q = negative_quantile(method)
    if method in PEAK_SIDE:
        c, w, lvl = peak_level(v, k, lowest=(method == 'peak_low'), side=PEAK_SIDE[method])
    elif q is not None:
        parts = ref if isinstance(ref, list) else [ref]
        c, w, lvls = np.nan, np.nan, []
        for part in parts:
            o = x[part & np.isfinite(x)] if part is not None else np.array([])
            enough = len(o) * (1 - q / 100) >= reg.const('KNOWN_NEGATIVE_MIN_ABOVE')
            lvls.append(float(np.log1p(np.percentile(o, q))) if enough else np.nan)
        lvl = float(max(lvls)) if lvls and all(np.isfinite(lvls)) else np.nan
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
# registry excludes_overlap: a reference or check cell whose outline overlaps these objects at all is left out
OVERLAP = {'Blood vessel': 'crossed_frac_vessel', 'Nerve': 'crossed_frac_nerve', 'Lymphatic': 'crossed_frac_lymphatic'}
# the reference populations of known-negative levels and their checks (registry level_reference, check_reference); every
# population is defined from registry roles and lineages
HORMONES_AND_LEUKOCYTES = 'endocrine-object hormone-positive cells; leukocytes'
POPULATION_STAGE = {'hormone': 0, 'leukocyte': 0, 'exocrine': 0, 'epithelial_away': 0, 'structure_cells': 0,
                    't_cells': 0, 'acinar_away': 1, 'exocrine_clean': 2, 'leukocyte_clean': 2}


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
    """(endocrine-object hormone-positive cells, leukocytes): the established hormone-named endocrine cells (owner
    2026-10-28: the endocrine typing is settled; by this module's hormone levels where no established typing exists);
    leukocytes need the Leukocyte marker and a T or Myeloid marker (two markers), None until those are chosen."""
    n = len(d)
    est = d['established_class'].isin(set(R['subtype_of'].values())).to_numpy()
    hormone = (est if est.any() else any_of(fixed, [h for h in R['hormones'] if h in fixed], n)) & \
        d['in_endocrine'].to_numpy()
    if not all(m in fixed for m in R['core']):
        return hormone, None
    return hormone, any_of(fixed, R['leukocyte'], n) & any_of(fixed, R['t'] + R['myeloid'], n)


def carrier_specific(reg: RG.Registry, m: str, markers) -> list[str]:
    """Defining markers whose every lineage is one of marker m's carrier lineages (registry carrier_lineages)."""
    lins = set(RG.split(reg.table.loc[m, 'carrier_lineages'])) if m in reg.table.index else set()
    return [x for x in markers if reg.role(x) == 'defining' and reg.lineage_entries(x)
            and {lin for lin, _ in reg.lineage_entries(x)} <= lins]


def populations(d: pd.DataFrame, reg: RG.Registry, m: str, fixed: dict, xy: np.ndarray, cut: float,
                R: dict) -> dict:
    """Named populations for marker m's known-negative level and checks, from the chosen positivity of the others:
      hormone          established hormone-named endocrine cells with no leukocyte within `cut` (itself included):
                       insulitic islets hold T cells that carry hormone spillover (method repair 2026-10-28)
      leukocyte        the Leukocyte marker with a T or Myeloid marker
      exocrine         ordinary-tissue cells positive for an Epithelial majority marker, not leukocyte- or hormone-cells
      epithelial_away  an Epithelial majority marker positive, the Leukocyte marker negative, no leukocyte within `cut`
      acinar_away      exocrine, the Exocrine-specific (ductal) markers negative, no leukocyte within `cut`
      exocrine_clean   ordinary tissue, every Epithelial majority marker positive, no defining marker of another
                       compartment and none specific to m's carrier lineages positive
      leukocyte_clean  leukocytes outside endocrine and Nerve objects with no marker specific to m's carrier lineages
      structure_cells  cells in a structure object that are neither leukocytes nor hormone cells
      t_cells          the T marker with the Leukocyte marker"""
    n = len(d)
    hormone, leuk = cell_sets(d, fixed, R)
    leuk = leuk if leuk is not None else np.zeros(n, bool)
    in_endo, ordinary = d['in_endocrine'].to_numpy(), d['ordinary'].to_numpy()
    L = any_of(fixed, R['leukocyte'], n)
    epi = any_of(fixed, [x for x in R['epithelial'] if x in fixed], n)
    near_leuk = ~isolated_from(leuk, np.ones(n, bool), xy, cut)
    pops = {'hormone': hormone & ~near_leuk, 'leukocyte': leuk, 'exocrine': ordinary & epi & ~L & ~hormone,
            'epithelial_away': epi & ~L & ~near_leuk,
            'structure_cells': d['context'].isin(list(STRUCTURE_CONTEXT.values())).to_numpy() & ~leuk & ~hormone,
            't_cells': any_of(fixed, [x for x in R['t'] if x in fixed], n) & L}
    pops['acinar_away'] = pops['exocrine'] & ~any_of(fixed, [x for x in R['ductal'] if x in fixed], n) & ~near_leuk
    carriers = carrier_specific(reg, m, list(fixed))
    others = [x for x in fixed if reg.role(x) == 'defining' and reg.table.loc[x, 'compartment'] != 'Epithelial']
    all_epi = np.ones(n, bool)
    for x in R['epithelial']:
        all_epi &= fixed[x] if x in fixed else np.zeros(n, bool)
    pops['exocrine_clean'] = ordinary & all_epi & ~any_of(fixed, others + carriers, n)
    nerve = d['context'].eq(STRUCTURE_CONTEXT['Nerve objects']).to_numpy()
    pops['leukocyte_clean'] = leuk & ~any_of(fixed, carriers, n) & ~in_endo & ~nerve
    return pops


def overlap_free(d: pd.DataFrame, reg: RG.Registry, m: str) -> np.ndarray:
    ov = reg.table.loc[m, 'excludes_overlap'] if 'excludes_overlap' in reg.table.columns else ''
    return d[OVERLAP[ov]].fillna(0).to_numpy() == 0 if ov else np.ones(len(d), bool)


def negative_reference(d: pd.DataFrame, reg: RG.Registry, m: str, meas: str, fixed: dict, xy: np.ndarray, k: float,
                       cut: float, R: dict, pops: dict | None = None):
    """The reference cells a known-negative level is read from (registry level_reference: '+'-joined populations;
    excludes_overlap; level_isolation: 'peak_left' / 'peak_narrow' = no cell outside the reference positive at the
    marker's own background-peak level within `cut`, as for SMA, or 'none'). level_combine 'each' returns one mask per
    population (each lineage held to its own percentile), 'union' one mask."""
    pops = pops if pops is not None else populations(d, reg, m, fixed, xy, cut, R)
    names = RG.split(reg.table.loc[m, 'level_reference'].replace('+', ';'))
    clean = overlap_free(d, reg, m)
    parts = [pops[x] & clean for x in names]
    union = np.logical_or.reduce(parts)
    iso = reg.table.loc[m, 'level_isolation']
    if iso in ('peak_left', 'peak_narrow'):
        src = level_for(d, reg, m, meas, {'peak_left': 'peak_dense', 'peak_narrow': 'peak_narrow'}[iso], k,
                        np.nan)['pos'] & ~union
        parts = [isolated_from(src, p_, xy, cut) for p_ in parts]
    elif iso != 'none':
        raise SystemExit(f'{m}: unknown level_isolation {iso!r}')
    return parts if reg.table.loc[m, 'level_combine'] == 'each' else np.logical_or.reduce(parts)


def dependency_rank(reg: RG.Registry, m: str) -> int:
    """Order in which levels are computed: background-peak levels first, a level read on a background population after
    that population's marker, known-negative levels after the markers their reference populations use."""
    if reg.positivity_class(m) != 'known_negative':
        return 1 if reg.table.loc[m, 'background_population'] else 0
    names = RG.split(reg.table.loc[m, 'level_reference'].replace('+', ';'))
    return 2 + max(POPULATION_STAGE[x] for x in names)


def chosen_positivity(d: pd.DataFrame, reg: RG.Registry, chosen: dict, k: float, nulls: dict, tag: str,
                      xy: np.ndarray, cut: float, R: dict) -> tuple[dict, dict]:
    """Each chosen marker's level and positivity on this image, in dependency order (dependency_rank)."""
    fixed, lvs = {}, {}
    for m in sorted(chosen, key=lambda x: dependency_rank(reg, x)):
        meas, method = chosen[m]
        ref = negative_reference(d, reg, m, meas, fixed, xy, k, cut, R) if negative_quantile(method) else None
        bp = reg.table.loc[m, 'background_population']
        subset = any_of(fixed, reg.broad(bp), len(d)) if bp else None
        lvs[m] = level_for(d, reg, m, meas, method, k, nulls.get((tag, m), np.nan), ref, subset)
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

    def ok_all(xs):
        # a gate that cannot be read on this image (fewer than MIN_CELLS_FOR_RATE reference cells) is not a pass:
        # the method is not shown to work there (review of 2026-10-28)
        return all(False if (isinstance(x, float) and np.isnan(x)) else bool(x) for x in xs)

    def enrichment(p, inside, base=None):
        base = ordinary if base is None else base
        ri, ro = rate(p, inside, nmin), rate(p, base, nmin)
        if not (np.isfinite(ri) and np.isfinite(ro)):
            enr = np.nan
        else:
            enr = ri / ro if ro > 0 else np.inf
        return ri, ro, enr, ((enr >= reg.const('OBJECT_ENRICHMENT_MIN')) if not np.isnan(enr) else np.nan)

    if stage == 'objects':
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) in ('structure_object', 'endocrine_object')]:
            ref = reg.table.loc[m, 'reference']
            if reg.positivity_class(m) == 'endocrine_object':
                # the established cells the marker names (a hormone: its subtype; the broad marker: every hormone-named
                # cell), so the check does not presume the cells are many (two T1D images have few beta cells)
                sub = R['subtype_of'].get(m)
                est = d['established_class']
                inside = ((est == sub) if sub else est.isin(set(R['subtype_of'].values()))).to_numpy() & inside_endo
                inside = inside if inside.any() else inside_endo
            else:
                inside = d['context'].eq(STRUCTURE_CONTEXT[ref]).to_numpy()
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, np.nan)
                ri, ro, enr, ok = enrichment(lv['pos'], inside)
                row(m, meas, method, lv, sensitivity=ri, out_of_place=ro, enrichment=enr, passes=ok)
    elif stage == 'immune':
        # one option per marker since the method repair (registry background_method); the background population's
        # marker first (CD3e's background is read on CD45-positive cells)
        core = sorted(R['core'], key=lambda x: dependency_rank(reg, x))
        pos, lv = {}, {}
        for m in core:
            (meas, method), = options(reg, m)
            bp = reg.table.loc[m, 'background_population']
            lv[m] = level_for(d, reg, m, meas, method, k, np.nan,
                              subset=any_of(pos, reg.broad(bp), n) if bp else None)
            pos[m] = lv[m]['pos']
        t, l, mye = any_of(pos, R['t'], n), any_of(pos, R['leukocyte'], n), any_of(pos, R['myeloid'], n)
        mye_all = np.logical_and.reduce([pos[x] for x in R['myeloid']])
        leuk_c = l & (t | mye)
        hormone_free = isolated_from(leuk_c, hormone_cells, xy, cut)       # no leukocyte within cut, itself included
        c45 = rate(l, t, nmin)
        t_all, t_iso, _ = isolated_rate(t, mye_all, xy, cut, nmin)
        checks = {'leukocyte_among_t': c45, 't_among_myeloid_all': t_all, 't_among_myeloid': t_iso,
                  'hormone_reference_cells': int(hormone_free.sum())}
        for m in core:
            r_all, r_iso, _ = isolated_rate(pos[m], hormone_free, xy, cut, nmin)
            checks[f'{m}_in_hormone_all'], checks[f'{m}_in_hormone'] = r_all, r_iso
        ok = [(c45 >= reg.const('CD45_AMONG_CD3E_MIN')) if np.isfinite(c45) else np.nan,
              le(t_iso, reg.const('CD3E_AMONG_MYELOID_MAX'))] + \
             [le(checks[f'{m}_in_hormone'], reg.const('CD45_AMONG_HORMONE_MAX')) for m in core]
        good = ok_all(ok) and all(np.isfinite(lv[m]['level_log']) for m in core)
        oop = np.nanmax([t_iso] + [checks[f'{m}_in_hormone'] for m in core])
        checks.update({f'level_raw_{m}': lv[m]['level_raw'] for m in core})
        order = R['core']
        opts = {m: options(reg, m)[0] for m in order}
        rows.append({'image': tag, 'stage': stage, 'marker': '+'.join(order),
                     'measurement': '|'.join(opts[m][0] for m in order), 'method': '|'.join(opts[m][1] for m in order),
                     'out_of_place': oop, 'sensitivity': c45, 'passes': good, **checks})
    elif stage == 'immune2':
        t = any_of(fixed, R['t'], n)
        hormone_free = isolated_from(leuk, hormone_cells, xy, cut)
        for m in R['immune2']:
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, np.nan)
                p = lv['pos']
                inh_all, inh, _ = isolated_rate(p, hormone_free, xy, cut, nmin)
                ok = le(inh, reg.const('CD45_AMONG_HORMONE_MAX'))
                row(m, meas, method, lv, out_of_place=inh, sensitivity=rate(p, t, nmin), passes=ok_all([ok]),
                    in_hormone_all=inh_all, in_hormone=inh)
    elif stage == 'epithelial':
        pops = populations(d, reg, '', fixed, xy, cut, R)
        for m in [m for m in reg.measured_markers() if reg.positivity_class(m) == 'majority_epithelial']:
            check = reg.table.loc[m, 'check_reference'] or 'leukocyte'
            for meas, method in options(reg, m):
                lv = level_for(d, reg, m, meas, method, k, nulls.get((tag, m), np.nan))
                r_all, r_iso, n_iso = isolated_rate(lv['pos'], pops[check], xy, cut, nmin)
                _, l_iso, _ = isolated_rate(lv['pos'], leuk, xy, cut, nmin)
                _, t_iso, _ = isolated_rate(lv['pos'], pops['t_cells'], xy, cut, nmin)
                row(m, meas, method, lv, out_of_place=r_iso, out_of_place_all=r_all, isolated_cells=n_iso,
                    in_leukocytes=l_iso, in_t_cells=t_iso, sensitivity=rate(lv['pos'], ordinary, nmin),
                    passes=ok_all([le(r_iso, reg.const('EPITHELIAL_IN_IMMUNE_MAX'))]))
    elif stage == 'known_negative':
        local = dict(fixed)                 # earlier known-negative markers feed later references (dependency_rank)
        for m in sorted([m for m in reg.measured_markers() if reg.positivity_class(m) == 'known_negative'],
                        key=lambda x: dependency_rank(reg, x)):
            pops = populations(d, reg, m, local, xy, cut, R)
            for meas, method in options(reg, m):
                ref = negative_reference(d, reg, m, meas, local, xy, k, cut, R, pops)
                lv = level_for(d, reg, m, meas, method, k, np.nan, ref)
                p = lv['pos']
                ref_all = np.logical_or.reduce(ref) if isinstance(ref, list) else ref
                extra = {'reference_cells': int(ref_all.sum())}
                gates = []
                bar_name = reg.table.loc[m, 'check_bar']
                bar = reg.const(bar_name) if bar_name else np.nan
                cr = reg.table.loc[m, 'check_reference']
                if cr:
                    cref = pops[cr] & overlap_free(d, reg, m)
                    r_all, r_iso, n_iso = isolated_rate(p, cref, xy, cut, nmin)
                    extra.update({'out_of_place': r_iso, 'out_of_place_all': r_all, 'isolated_cells': n_iso})
                    gates.append(le(r_iso, bar))
                for spec in RG.split(reg.table.loc[m, 'check_objects']):
                    obj = spec.replace(' reported', '').strip()
                    inside = inside_endo if obj == 'endocrine objects' else d['context'].eq(
                        STRUCTURE_CONTEXT[obj]).to_numpy()
                    ri, ro, enr, ok = enrichment(p, inside, ordinary & ~ref_all)
                    key = obj.split()[0].lower()
                    extra.update({f'sensitivity_{key}': ri, f'enrichment_{key}': enr})
                    if not spec.endswith('reported'):
                        gates.append(ok)
                co = reg.table.loc[m, 'coexpressed_with']
                if co:
                    lacking = rate(~local[co], p, nmin)
                    extra['positive_without_' + co] = lacking
                    gates.append(le(lacking, bar))
                cf = reg.table.loc[m, 'cross_fit']
                if cf:
                    # each half of a 500-um checkerboard predicts the other's false-positive rate in `cf`
                    parity = ((d['x_um'] // 500).astype(int) + (d['y_um'] // 500).astype(int)).to_numpy() % 2
                    worst = 0.0
                    for a in (0, 1):
                        half = [r_ & (parity == a) for r_ in (ref if isinstance(ref, list) else [ref])]
                        la = level_for(d, reg, m, meas, method, k, np.nan, half if isinstance(ref, list) else half[0])
                        held = pops[cf] & (parity != a)
                        _, r_iso, _ = isolated_rate(la['pos'], held, xy, cut, nmin)
                        worst = np.nanmax([worst, r_iso])
                    extra['cross_fit_out_of_place'] = worst
                    gates.append(le(worst, bar))
                if reg.validation(m) == 'reported':
                    extra.update(lineage_report(reg, m, p, local, negative_quantile(method), nmin))
                    gates = [True]
                inh_all, inh, _ = isolated_rate(p, pops['hormone'], xy, cut, nmin)
                row(m, meas, method, lv, in_hormone=inh, passes=ok_all(gates) if gates else True, **extra)
                local[m] = p
    elif stage == 'shared':
        pops = populations(d, reg, '', fixed, xy, cut, R)
        refs = {'endocrine-object hormone-positive cells': pops['hormone'], 'T cells': pops['t_cells'],
                HORMONES_AND_LEUKOCYTES: leuk | pops['hormone']}
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
    fails = other[other['passes'].eq(False)].groupby(['marker', 'measurement', 'method'])['image'].apply(
        lambda x: ' '.join(sorted(x))).rename('failing_other_images').reset_index()
    s = s.merge(fails, on=['marker', 'measurement', 'method'], how='left')
    # owner 2026-10-28: "Don't flag the images or markers, flag the methods you are using so they will work with the
    # images and markers": an option is usable only if it passes on every image
    s['ok'] = (s['images'] == len(select_on)) & (s['failing_images'] == 0) & s['failing_other_images'].isna()
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
                failing = ' '.join(x for x in (best['failing'], best.get('failing_other_images'))
                                   if isinstance(x, str) and x)
                decisions.append({'marker': m, 'decision': 'use' if best['ok'] else 'stop',
                                  'measurement': meas[i], 'method': meth[i], 'out_of_place': best['out_of_place'],
                                  'sensitivity': best['sensitivity'],
                                  'why': '' if best['ok'] else f'method fails on: {failing}'})
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
        print(f'STOP: methods that do not work on every image: {", ".join(stops["marker"])}', flush=True)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
