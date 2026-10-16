#!/usr/bin/env python3
"""Compare an image's new cells with its previous segmentation, by mutual nearest centroid.

For review only: the earlier segmentations (InstanSeg 0.1.0 for the original 15 donors, 0.1.1 on 22 channels for
the additions) are references, not targets. Old centroids come from the native exports of the source projects
(`native_centroids.tsv`, full-resolution pixels); new cells from run_uniform.py's cells Parquet (micrometres).

    python scripts/uniform/compare_segmentation.py --new OUT/cells/Aab_6505.parquet \\
        --old data/development/all-donor-20260921-restored/native_cells/6505/native_centroids.tsv --out report.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--new', type=Path, required=True)
    ap.add_argument('--old', type=Path, required=True)
    ap.add_argument('--pixel-um', type=float, default=0.3774)
    ap.add_argument('--radius-um', type=float, default=3.0)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    new = pq.read_table(a.new, columns=['Centroid X µm', 'Centroid Y µm', 'Cell: Area µm^2', 'Nucleus: Area µm^2']).to_pandas()
    old = pd.read_csv(a.old, sep='\t')
    p_new = new[['Centroid X µm', 'Centroid Y µm']].to_numpy() / a.pixel_um
    p_old = old[['x_source_px', 'y_source_px']].to_numpy()
    area_new = new['Cell: Area µm^2'].to_numpy()
    area_old = old['cell_area_px2'].to_numpy() * a.pixel_um ** 2
    r_px = a.radius_um / a.pixel_um
    d_on, i_on = cKDTree(p_new).query(p_old, distance_upper_bound=r_px)   # each old cell's nearest new cell
    d_no, i_no = cKDTree(p_old).query(p_new, distance_upper_bound=r_px)   # each new cell's nearest old cell
    ok_old = np.isfinite(d_on)
    mutual = np.zeros(len(p_old), dtype=bool)
    mutual[ok_old] = i_no[i_on[ok_old]] == np.flatnonzero(ok_old)
    j = i_on[mutual]
    ratio = area_new[j] / area_old[mutual]
    matched_new = np.zeros(len(p_new), dtype=bool)
    matched_new[j] = True
    report = {
        'new_cells': int(len(p_new)), 'old_cells': int(len(p_old)),
        'difference_percent': round(100 * (len(p_new) - len(p_old)) / len(p_old), 2),
        'radius_um': a.radius_um, 'mutual_matches': int(mutual.sum()),
        'old_matched_percent': round(100 * mutual.mean(), 2), 'new_matched_percent': round(100 * matched_new.mean(), 2),
        'matched_centroid_distance_um_median': round(float(np.median(d_on[mutual]) * a.pixel_um), 3),
        'matched_area_ratio_new_over_old': {q: round(float(np.percentile(ratio, p)), 3) for q, p in
                                            (('p10', 10), ('median', 50), ('p90', 90))},
        'cell_area_um2_median': {'new': round(float(np.nanmedian(area_new)), 1), 'old': round(float(np.nanmedian(area_old)), 1),
                                 'new_unmatched': round(float(np.nanmedian(area_new[~matched_new])), 1) if (~matched_new).any() else None,
                                 'old_unmatched': round(float(np.nanmedian(area_old[~mutual])), 1) if (~mutual).any() else None},
        'new_cells_without_nucleus': int(new['Nucleus: Area µm^2'].isna().sum()),
        'new_file': str(a.new), 'old_file': str(a.old),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
