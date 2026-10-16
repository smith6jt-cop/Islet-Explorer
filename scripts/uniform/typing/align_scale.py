#!/usr/bin/env python3
"""align_scale.py -- Senior's UniFORM alignment and cohort scale, on this panel's REDSEA values.

Port of Islet-Explorer-Senior-v3 `core/phenocycler/scripts/celltype/run_step2.py` (derive, apply_shift),
`uniform_adapter.py` (marker_histograms, background_mode_spread, the pinned UniFORM import) and
`run_step3.build_scale` / `score.marker_scale`:

  1. each marker is read in its own compartment (codex_typing.COMPARTMENT) from `typing/redsea/<image>.parquet`;
  2. UniFORM's automatic_registration aligns each marker's background mode across the 20 images (one image per
     donor) as a rigid log-shift -- one multiplicative factor per (donor, marker); the histograms are
     UniFORM's own (log of values >= 1, common range, 1024 bins) over the analysis-eligible cells;
  3. a marker is left unaligned only when its donor medians step with the acquisition group (ratio >= 2.0,
     usability.acquisition_group_ratio; here the groups are the original 15 and the additions) AND the shift
     leaves the background modes more than 0.10 apart in log (run_step2.REFUSED);
  4. the cohort scale per marker: floor = the aligned background mode (512-bin log1p histogram of positive
     values), top = the cohort p99.9, over a seeded 200,000-cell sample per donor; scaled = clip((x - floor) /
     (top - floor), 0, 1).

CODEX difference: no reference donor is pinned (Senior pinned the reference donors its scale was first measured
on, `data/development/celltype/shifts.json`); UniFORM chooses, and the choice is recorded.

    python scripts/uniform/typing/align_scale.py --out data/development/uniform-20donor-20261015
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import codex_typing as CT  # noqa: E402

UNIFORM_PATH = Path.home() / 'IO60panc2nd/Islet-Explorer-Senior-v3/external/UniFORM'
UNIFORM_COMMIT = '7c33396d12e5a19dc0b36e46a94c3293de0cd814'
NUM_BINS = 1024
BAR_MODE_SD = 0.10
ACQ_RATIO_FLAG = 2.0
REFUSED = 'refused_reagent_change'
SCALE_N_PER_DONOR = 200_000
SCALE_Q = 99.9
SEED = 0
QC_COLUMNS = ['object_id', 'qc_analysis_eligible', 'qc_estimation_eligible', 'qc_reasons']


def import_uniform():
    head = subprocess.run(['git', '-C', str(UNIFORM_PATH), 'rev-parse', 'HEAD'], capture_output=True,
                          text=True, check=True).stdout.strip()
    if head != UNIFORM_COMMIT:
        raise SystemExit(f'UniFORM checkout is at {head}; this port is pinned to {UNIFORM_COMMIT}')
    sys.path.insert(0, str(UNIFORM_PATH))
    from UniFORM import normalization, registration  # noqa: E402
    return registration, normalization


def marker_histograms(values_by_donor: dict, num_bins: int = NUM_BINS):
    logs, mins, maxs = {}, [], []
    for d, v in values_by_donor.items():
        v = np.asarray(v, dtype=float)
        v = v[np.isfinite(v) & (v >= 1.0)]
        lv = np.log(v) if len(v) else np.array([0.0])
        logs[d] = lv
        mins.append(float(lv.min()))
        maxs.append(float(lv.max()))
    lo, hi = min(mins), max(maxs)
    counts = {d: np.histogram(lv, bins=num_bins, range=(lo, hi))[0] for d, lv in logs.items()}
    return counts, {'global_min': lo, 'global_max': hi, 'min_list': mins, 'max_list': maxs}


def background_mode_spread(counts, ranges, shifts=None, num_bins=NUM_BINS):
    lo, hi = ranges['global_min'], ranges['global_max']
    binw = (hi - lo) / (num_bins - 1)
    peaks = {}
    for d, c in counts.items():
        peak = lo + (int(np.argmax(np.asarray(c, float)[1:])) + 1) * binw
        if shifts is not None:
            peak -= shifts[d] * binw
        peaks[d] = peak
    v = np.array(list(peaks.values()))
    return {'peaks': peaks, 'sd_log': float(v.std()), 'fold': float(np.exp(v.max() - v.min()))}


def acquisition_group_ratio(level: pd.Series, group: pd.Series) -> dict:
    df = pd.DataFrame({'level': level, 'group': group}).dropna()
    df = df[df['level'] > 0]
    if df['group'].nunique() < 2:
        return {'acq_between': np.nan, 'acq_within': np.nan, 'acq_ratio': np.nan}
    g = df.groupby('group')['level']
    between = float(g.median().max() / max(g.median().min(), 1e-12))
    within = float(g.apply(lambda s: s.max() / max(s.min(), 1e-12)).max())
    return {'acq_between': between, 'acq_within': within, 'acq_ratio': float(between / max(within, 1e-12))}


def read_measured(out: Path, tags: list[str]) -> tuple[dict, list[str]]:
    """image -> frame of QC columns plus each marker in its own compartment, named by the marker."""
    frames, markers = {}, set()
    for t in tags:
        src = pd.read_parquet(CT.typing_dir(out) / 'redsea' / f'{t}.parquet')
        cols = {f'{CT.COMPARTMENT[m]}__{m}': m for m in CT.COMPARTMENT
                if f'{CT.COMPARTMENT[m]}__{m}' in src.columns and m not in CT.NOT_FOR_TYPING}
        f = src[QC_COLUMNS + list(cols)].rename(columns=cols)
        frames[t] = f
        markers |= set(cols.values())
    order = [m for m in CT.COMPARTMENT if m in markers]
    return frames, order


def derive(frames: dict, markers: list[str], groups: dict, reg, norm) -> dict:
    result = {}
    for m in markers:
        t0 = time.time()
        vals = {}
        unavailable = {}
        for d, f in frames.items():
            if m not in f.columns:
                continue
            v = f.loc[f['qc_analysis_eligible'], m].to_numpy('float64')
            if not np.any(np.isfinite(v) & (v >= 1.0)):
                unavailable[d] = 'no finite analysis-eligible value at or above the log floor (1.0)'
                continue
            vals[d] = v
        have = sorted(vals, key=list(frames).index)
        level = pd.Series({d: float(np.nanmedian(vals[d])) for d in have})
        acq = acquisition_group_ratio(level, pd.Series({d: groups[d] for d in have}))
        flagged = bool(np.isfinite(acq['acq_ratio']) and acq['acq_ratio'] >= ACQ_RATIO_FLAG)
        if len(have) < 2:
            result[m] = {'reference_donor': have[0] if have else 'unavailable', 'scale': {d: 1.0 for d in frames},
                         'provenance': 'too_few_donors_to_align', 'n_donors': len(have), **acq}
            continue
        counts, ranges = marker_histograms(vals)
        before = background_mode_spread(counts, ranges)
        hist = {m: {'counts': [counts[d] for d in have], 'bin_edges': [None] * len(have)}}
        with contextlib.redirect_stdout(io.StringIO()):
            shift_map, refs, _ = reg.automatic_registration(hist, [m], [m], have, None, None, NUM_BINS)
        import matplotlib.pyplot as plt
        plt.close('all')
        bins = dict(zip(have, shift_map[m]))
        after = background_mode_spread(counts, ranges, bins)
        log_shift = dict(zip(have, norm.calculate_shift_in_log_pixels({m: ranges}, [m], NUM_BINS, shift_map)[m]))
        scale = {d: float(np.exp(-log_shift[d])) if d in log_shift else 1.0 for d in frames}
        refused = flagged and after['sd_log'] > BAR_MODE_SD
        if refused:
            scale = {d: 1.0 for d in frames}
        for d in unavailable:
            scale[d] = 'unavailable'
        over_bar = (not refused) and after['sd_log'] > BAR_MODE_SD
        result[m] = {
            'reference_donor': refs[0], 'reference_source': 'uniform_choice',
            'shift_bins': {d: int(bins[d]) for d in have}, 'shift_log': {d: float(log_shift[d]) for d in have},
            'scale': scale,
            'provenance': REFUSED if refused else ('automatic_registration_over_bar' if over_bar
                                                    else 'automatic_registration'),
            'over_bar': bool(over_bar), 'acq_ratio_flagged': flagged, **acq,
            'n_donors': len(have), 'unavailable_donors': unavailable,
            'mode_sd_before': before['sd_log'], 'mode_sd_after': after['sd_log'],
            'mode_fold_before': before['fold'], 'mode_fold_after': after['fold'],
            'global_min': ranges['global_min'], 'global_max': ranges['global_max'],
        }
        tag = 'REFUSED' if refused else ('ok' if after['sd_log'] <= BAR_MODE_SD else 'OVER BAR')
        print(f'  {m:<10} ref={refs[0]:<9} mode fold {before["fold"]:7.2f} -> {after["fold"]:5.2f}   '
              f'sd {before["sd_log"]:.3f} -> {after["sd_log"]:.3f}   acq {acq["acq_ratio"]:.2f}   {tag}   '
              f'{time.time() - t0:.0f}s', flush=True)
    return result


def build_scale(aligned: dict, markers: list[str]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    pool = {m: [] for m in markers}
    for d, f in aligned.items():
        t = f[f['qc_analysis_eligible']]
        if len(t) > SCALE_N_PER_DONOR:
            t = t.iloc[rng.choice(len(t), SCALE_N_PER_DONOR, replace=False)]
        for m in markers:
            if m in t.columns:
                pool[m].append(t[m].to_numpy('float64'))
    rows = []
    for m, chunks in pool.items():
        if not chunks:
            continue
        v = np.concatenate(chunks)
        pos = v[np.isfinite(v) & (v > 0)]
        floor = 0.0
        if len(pos) >= 1000:
            h, e = np.histogram(np.log1p(pos), bins=512)
            floor = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))
        fin = v[np.isfinite(v)]
        top = float(np.percentile(fin, SCALE_Q)) if len(fin) else np.nan
        rows.append({'marker': m, 'floor': floor, 'top': top, 'span': top - floor if np.isfinite(top) else np.nan,
                     'n': int(len(fin))})
    return pd.DataFrame(rows).set_index('marker')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    rows = CT.cohort_rows()
    tags = [r['image_name'] for r in rows]
    groups = {r['image_name']: r['cohort'] for r in rows}
    typ = CT.typing_dir(a.out)
    missing = [t for t in tags if not (typ / 'redsea' / f'{t}.parquet').exists()]
    if missing:
        raise SystemExit(f'no REDSEA values yet for {missing}')
    reg, norm = import_uniform()
    frames, markers = read_measured(a.out, tags)
    print(f'{len(markers)} markers over {len(tags)} images', flush=True)
    shifts = derive(frames, markers, groups, reg, norm)

    (typ / 'aligned').mkdir(parents=True, exist_ok=True)
    aligned = {}
    for d, f in frames.items():
        g = f.copy()
        for m in markers:
            if m in g.columns and isinstance(shifts[m]['scale'].get(d), float):
                g[m] = (g[m].to_numpy('float64') * shifts[m]['scale'][d]).astype(np.float32)
        g.to_parquet(typ / 'aligned' / f'{d}.parquet', index=False)
        aligned[d] = g
    scale = build_scale(aligned, markers)
    scale.to_parquet(typ / 'scale.parquet')
    scale.to_csv(typ / 'scale.csv')
    (typ / 'shifts.json').write_text(json.dumps({
        'uniform_commit': UNIFORM_COMMIT, 'num_bins': NUM_BINS, 'bar_mode_sd': BAR_MODE_SD,
        'acq_ratio_flag': ACQ_RATIO_FLAG, 'acquisition_groups': groups,
        'compartments': {m: CT.COMPARTMENT[m] for m in markers}, 'shifts': shifts}, indent=1, default=float))
    print(scale.round(2).to_string(), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
