#!/usr/bin/env python3
"""Read-only endocrine measurement trace on native pixels and final cell masks.

For every unresolved endocrine-object cell and every currently named endocrine
cell, compare the exact raw raster mean, the 3x3-median raster mean, and stored
REDSEA whole-cell mean. Use the production rasterizer, feature ordering and
stored spillover-context eligibility; verify pixel counts against stored values.
No threshold, geometry or identity is changed. Stage comparisons use the current
corrected-value threshold only as a diagnostic, not as alternate typing calls.

Run on HPG via run_in_site.sh / site_job.sbatch. Output must be a fresh directory.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts/uniform/typing'))
import codex_typing as CT
import redsea_cells as RS

MARKERS = ('INS', 'GCG', 'SST', 'CHGA')
UNRESOLVED = ('Islet', 'EndocrineCluster', 'EndocrineSingle')
NAMED = ('Beta cell', 'Alpha cell', 'Delta cell', 'CHGA+ hormone-negative')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    st = path.stat()
    return {'path': str(path), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns}


def finite_corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
        return None
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def summarize(d):
    out = {'cells': len(d), 'mask_count_mismatches': int((~d.mask_counts_match_stored).sum()),
           'raw_any_positive': int(d.raw_any_positive.sum()),
           'median_any_positive': int(d.median_any_positive.sum()),
           'redsea_any_positive': int(d.redsea_any_positive.sum()),
           'raw_positive_to_median_negative': int((d.raw_any_positive & ~d.median_any_positive).sum()),
           'median_positive_to_redsea_negative': int((d.median_any_positive & ~d.redsea_any_positive).sum()),
           'raw_negative_to_median_positive': int((~d.raw_any_positive & d.median_any_positive).sum()),
           'median_negative_to_redsea_positive': int((~d.median_any_positive & d.redsea_any_positive).sum()),
           'zero_cytoplasm_pixels': int((d.cell_pixels == d.nucleus_pixels).sum()), 'markers': {}}
    for marker in MARKERS:
        raw = d[f'raw_mean__{marker}'].to_numpy(float)
        median = d[f'median_mean__{marker}'].to_numpy(float)
        stored = d[f'redsea_mean__{marker}'].to_numpy(float)
        qupath = d[f'qupath_mean__{marker}'].to_numpy(float)
        out['markers'][marker] = {
            'raw_positive': int(d[f'raw_positive__{marker}'].sum()),
            'median_positive': int(d[f'median_positive__{marker}'].sum()),
            'redsea_positive': int(d[f'redsea_positive__{marker}'].sum()),
            'median_minus_raw_median': float(np.nanmedian(median - raw)),
            'redsea_minus_median_median': float(np.nanmedian(stored - median)),
            'raw_raster_vs_qupath_r': finite_corr(raw, qupath),
            'raw_raster_minus_qupath_median': float(np.nanmedian(raw - qupath)),
            'raw_raster_minus_qupath_max_abs': float(np.nanmax(np.abs(raw - qupath))),
        }
    return out


def trace_image(run, out, tag):
    started = time.time()
    identity = pd.read_parquet(run / 'typing/identity' / f'{tag}.parquet',
                               columns=['object_id', 'label', 'status'] + ['pos__' + m for m in MARKERS])
    selected = identity[identity.label.isin(UNRESOLVED + NAMED)].copy()
    if selected.empty:
        raise ValueError(f'{tag}: no endocrine or unresolved cells')
    stored = pd.read_parquet(run / 'typing/redsea' / f'{tag}.parquet', columns=[
        'object_id', 'qc_analysis_eligible', 'qc_spillover_context_eligible', 'qc_flags',
        'cell_area_px', 'nucleus_area_px', 'Cell__n_px', 'Nucleus__n_px',
    ] + ['Cell__' + m for m in MARKERS]).set_index('object_id', verify_integrity=True)
    q = pd.read_parquet(run / 'cells' / f'{tag}.parquet',
                        columns=['Object ID', 'Cell: Area µm^2', 'Nucleus: Area µm^2']
                        + [f'Cell: {m}: Mean' for m in MARKERS]).set_index('Object ID', verify_integrity=True)
    lv = pd.read_csv(run / 'typing/levels.csv')
    lv = lv[lv.image.eq(tag)].set_index('marker', verify_integrity=True)
    if not lv.loc[list(MARKERS), 'measurement'].eq('whole_cell').all():
        raise ValueError(f'{tag}: endocrine positivity is not whole-cell')
    entry = CT.manifest(run)[tag]
    channels = CT.channels(run, entry['entry_id'])
    img = RS.Channels(entry['image_path'])
    try:
        if img.shape[0] != len(channels):
            raise ValueError(f'{tag}: source image / channel-name mismatch')
        features = json.loads((run / 'cell_geojson' / f'{tag}.geojson').read_text())
        features = features['features'] if isinstance(features, dict) else features
        mask, nuc, ids, _ = RS.rasterize(features, img.shape[1:])
        del features
        if len(ids) != len(set(ids)) or set(ids) != set(stored.index) or set(ids) != set(q.index):
            raise ValueError(f'{tag}: geometry / REDSEA / QuPath identity mismatch')
        stored = stored.reindex(ids)
        n = len(ids)
        counts0 = np.bincount(mask.ravel(), minlength=n + 1)[1:]
        nuc0 = np.bincount(nuc.ravel(), minlength=n + 1)[1:]
        initial_mismatch = (counts0 != stored.cell_area_px.to_numpy()) | (nuc0 != stored.nucleus_area_px.to_numpy())
        eligible = stored.qc_spillover_context_eligible.to_numpy(bool)
        keep = np.r_[True, eligible]
        # Exactly the mask filtering in RS.process; excluded labels become slide background.
        np.putmask(mask, ~keep[mask], 0)
        np.putmask(nuc, ~keep[nuc], 0)
        counts = np.bincount(mask.ravel(), minlength=n + 1)[1:]
        nuc_counts = np.bincount(nuc.ravel(), minlength=n + 1)[1:]
        final_mismatch = (counts != stored['Cell__n_px'].to_numpy()) | (nuc_counts != stored['Nucleus__n_px'].to_numpy())
        del nuc, nuc0, counts0
        positions = pd.Series(np.arange(n), index=ids).loc[selected.object_id].to_numpy()
        s = stored.iloc[positions]
        if not s.qc_analysis_eligible.all() or not s.qc_spillover_context_eligible.all():
            raise ValueError(f'{tag}: a selected identity is not measurement-eligible')
        d = selected.reset_index(drop=True)
        d['image'] = tag
        d['population'] = np.where(d.label.isin(UNRESOLVED), 'unresolved', 'named_endocrine')
        d['cell_pixels'] = counts[positions]
        d['nucleus_pixels'] = nuc_counts[positions]
        d['stored_cell_pixels'] = s['Cell__n_px'].to_numpy()
        d['stored_nucleus_pixels'] = s['Nucleus__n_px'].to_numpy()
        d['mask_counts_match_stored'] = ~(initial_mismatch | final_mismatch)[positions]
        d['qc_flags'] = s.qc_flags.to_numpy()
        sq = q.reindex(d.object_id)
        d['qupath_cell_area_um2'] = sq['Cell: Area µm^2'].to_numpy()
        d['qupath_nucleus_area_um2'] = sq['Nucleus: Area µm^2'].to_numpy()
        flat = mask.ravel()
        denom = np.where(counts[positions] > 0, counts[positions], np.nan)
        for marker in MARKERS:
            raw = img.native(channels.index(marker))
            raw_sums = np.bincount(flat, weights=raw.ravel().astype(np.float32), minlength=n + 1)[1:]
            clean = RS.despeckle(raw)
            clean_sums = np.bincount(flat, weights=clean.ravel().astype(np.float32), minlength=n + 1)[1:]
            d[f'raw_mean__{marker}'] = raw_sums[positions] / denom
            d[f'median_mean__{marker}'] = clean_sums[positions] / denom
            d[f'redsea_mean__{marker}'] = s['Cell__' + marker].to_numpy()
            d[f'qupath_mean__{marker}'] = sq[f'Cell: {marker}: Mean'].to_numpy()
            level = float(lv.loc[marker, 'level_log'])
            d[f'level_raw__{marker}'] = float(lv.loc[marker, 'level_raw'])
            for stage in ('raw', 'median', 'redsea'):
                d[f'{stage}_positive__{marker}'] = np.log1p(d[f'{stage}_mean__{marker}'].to_numpy(float)) > level
            if not np.array_equal(d[f'redsea_positive__{marker}'], d['pos__' + marker].astype(bool)):
                raise ValueError(f'{tag}/{marker}: stored positivity disagrees with stored corrected values')
            del raw, clean, raw_sums, clean_sums
        del mask, flat
    finally:
        img.close()
    for stage in ('raw', 'median', 'redsea'):
        d[stage + '_any_positive'] = d[[f'{stage}_positive__{m}' for m in MARKERS]].any(axis=1)
    result = {
        'image': tag, 'status': 'pass' if not (initial_mismatch | final_mismatch).any() else 'mask_count_mismatch',
        'all_geometry_cells': n, 'all_initial_count_mismatches': int(initial_mismatch.sum()),
        'all_final_count_mismatches': int(final_mismatch.sum()),
        'by_population': {key: summarize(g) for key, g in d.groupby('population')},
        'by_label': {key: summarize(g) for key, g in d.groupby('label')},
        'inputs': {name: file_identity(run / name) for name in (
            f'cell_geojson/{tag}.geojson', f'cells/{tag}.parquet', f'typing/redsea/{tag}.parquet',
            f'typing/identity/{tag}.parquet')},
        'source_image': file_identity(entry['image_path']), 'seconds': round(time.time() - started, 1),
    }
    d.to_parquet(out / f'{tag}.parquet', index=False)
    (out / f'{tag}.json').write_text(json.dumps(result, indent=2) + '\n')
    print(tag, result['status'], len(d), 'cells', result['seconds'], 'seconds', flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--images', nargs='+', help='optional image subset; default all identity files')
    a = p.parse_args()
    if a.out.exists():
        p.error('output already exists; choose a fresh directory')
    tags = sorted(path.stem for path in (a.run / 'typing/identity').glob('*.parquet'))
    if a.images:
        if set(a.images) - set(tags):
            p.error('requested image has no identity file')
        tags = sorted(set(a.images))
    if not tags or a.workers < 1:
        p.error('identity files and a positive worker count are required')
    a.out.mkdir(parents=True)
    with ProcessPoolExecutor(a.workers) as pool:
        results = list(pool.map(trace_image, [a.run] * len(tags), [a.out] * len(tags), tags))
    report = {'status': 'pass' if all(r['status'] == 'pass' for r in results) else 'mask_count_mismatch',
              'note': __doc__, 'run': str(a.run.resolve()), 'images': results,
              'code_sha256': {str(path): sha(path) for path in (Path(__file__), Path(RS.__file__), Path(CT.__file__))},
              'levels_sha256': sha(a.run / 'typing/levels.csv')}
    (a.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    return 0 if report['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
