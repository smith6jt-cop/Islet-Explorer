#!/usr/bin/env python3
"""One review table for images run by run_uniform.py: Tissue, cells against the previous segmentation, the
classifier fit, the islet objects, and where the cells the earlier typing called endocrine ended up: inside script
1's outlines, and inside the objects step 1b leaves (only the cells that fit whole). Reads provenance and outputs and
reports only; nothing it computes enters the pipeline.

The earlier labels are references, not truth: the frozen 2025 typing for the original 15 donors and the accepted
reconstruction rules (2026-10-03) for the additions, both on the earlier segmentations. An earlier cell is matched
to a new cell by mutual nearest centroid within 3 µm.

    python scripts/uniform/summarize_review.py --out data/development/uniform-20donor-20261015 --images all
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OLD_CENTROIDS = REPO / 'data/development/all-donor-20260921-restored/native_cells/{case}/native_centroids.tsv'
FROZEN_H5AD = REPO / 'single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad'
ADDITION_H5AD = REPO / 'data/development/protein-expanded-20261015-v3/additions/{case}.h5ad'
PIXEL_UM = 0.3774
MATCH_UM = 3.0
ENDOCRINE = ('Beta cell', 'Alpha cell', 'Delta cell', 'Endocrine')
ISLET_PARENT = r'^(Islet|Cluster|Single)_\d+$'


def load(path: Path) -> dict:
    return json.loads(path.read_text()) if path.is_file() else {}


def _obs(f, name: str) -> np.ndarray:
    import h5py  # noqa: PLC0415
    g = f['obs'][name]
    if isinstance(g, h5py.Group):
        return np.asarray(g['categories'].asstr()[:], dtype=object)[g['codes'][:]]
    return g.asstr()[:] if g.dtype.kind in 'OS' else g[:]


_FROZEN: pd.DataFrame | None = None


def earlier_labels(case: str, cohort: str) -> pd.DataFrame:
    """object_id, phenotype, parent of the earlier segmentation's cells for one donor."""
    global _FROZEN
    import h5py  # noqa: PLC0415
    if cohort == 'original':
        if _FROZEN is None:
            with h5py.File(FROZEN_H5AD, 'r') as f:
                _FROZEN = pd.DataFrame({'imageid': _obs(f, 'imageid').astype(str), 'object_id': _obs(f, 'Object ID'),
                                        'phenotype': _obs(f, 'phenotype'), 'parent': _obs(f, 'Parent')})
        return _FROZEN[_FROZEN.imageid == case].drop(columns='imageid')
    with h5py.File(str(ADDITION_H5AD).format(case=case), 'r') as f:
        return pd.DataFrame({'object_id': _obs(f, 'Object ID'), 'phenotype': _obs(f, 'phenotype'), 'parent': _obs(f, 'Parent')})


def in_outlines(x_px: np.ndarray, y_px: np.ndarray, geojson: Path) -> np.ndarray:
    """Centroids inside any outline of a QuPath GeoJSON export (image pixels); all False without the file."""
    import shapely  # noqa: PLC0415
    from shapely import STRtree  # noqa: PLC0415
    from shapely.geometry import shape  # noqa: PLC0415
    inside = np.zeros(len(x_px), dtype=bool)
    if not geojson.is_file():
        return inside
    data = json.loads(geojson.read_text())
    geoms = [shape(f['geometry']) for f in (data['features'] if isinstance(data, dict) else data)]
    if geoms:
        hit, _ = STRtree(geoms).query(shapely.points(x_px, y_px), predicate='within')
        inside[np.unique(hit)] = True
    return inside


def where_endocrine_went(final_cells: Path, case: str, cohort: str, script1_outlines: Path | None = None) -> dict:
    """Share of the earlier endocrine cells (by label) whose matched new cell lies in an Islet / EndocrineCluster /
    EndocrineSingle object, and of the earlier islet-core cells; with `script1_outlines`, also the share inside
    script 1's outlines, before step 1b keeps only the cells that fit whole."""
    import pyarrow.parquet as pq  # noqa: PLC0415
    from scipy.spatial import cKDTree  # noqa: PLC0415
    old = pd.read_csv(str(OLD_CENTROIDS).format(case=case), sep='\t').merge(earlier_labels(case, cohort), on='object_id')
    new = pq.read_table(final_cells, columns=['Parent', 'Centroid X µm', 'Centroid Y µm']).to_pandas()
    p_new = new[['Centroid X µm', 'Centroid Y µm']].to_numpy() / PIXEL_UM
    p_old = old[['x_source_px', 'y_source_px']].to_numpy()
    r = MATCH_UM / PIXEL_UM
    d_on, i_on = cKDTree(p_new).query(p_old, distance_upper_bound=r)
    d_no, i_no = cKDTree(p_old).query(p_new, distance_upper_bound=r)
    ok = np.isfinite(d_on)
    mutual = np.zeros(len(old), dtype=bool)
    mutual[ok] = i_no[i_on[ok]] == np.flatnonzero(ok)
    in_object = new['Parent'].astype(str).str.match(ISLET_PARENT).to_numpy()
    in_script1 = in_outlines(p_new[:, 0], p_new[:, 1], script1_outlines) if script1_outlines else None
    j = i_on[mutual]
    old = old[mutual].assign(inside=in_object[j])
    if in_script1 is not None:
        old = old.assign(in_script1=in_script1[j])
    pct = lambda s: round(100 * float(s.mean()), 1) if len(s) else None  # noqa: E731
    out = {'earlier_label_source': 'frozen 2025 typing' if cohort == 'original' else 'accepted reconstruction rules 2026-10-03',
           'earlier_cells_matched': int(len(old))}
    for ph in ENDOCRINE:
        s = old[old.phenotype == ph]
        out[ph] = {'n': int(len(s)), 'inside_pct': pct(s.inside)}
        if in_script1 is not None:
            out[ph]['in_script1_outlines_pct'] = pct(s.in_script1)
    endo = old[old.phenotype.isin(ENDOCRINE)]
    out['endocrine_all'] = {'n': int(len(endo)), 'inside_pct': pct(endo.inside)}
    core = old[old.parent.astype(str).str.match(r'^Islet_\d+$')]
    out['earlier_islet_core_cells'] = {'n': int(len(core)), 'islets': int(core.parent.nunique()), 'inside_pct': pct(core.inside)}
    out['non_endocrine_inside'] = int((~old.phenotype.isin(ENDOCRINE) & old.inside).sum())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--images', nargs='+', required=True, help='image names or "all" (images without final cells are skipped)')
    ap.add_argument('--cohort', type=Path, default=HERE / 'cohort.csv')
    a = ap.parse_args()
    out = a.out.resolve()
    rows = {r['image_name']: r for r in csv.DictReader(a.cohort.open())}
    tags = list(rows) if a.images == ['all'] else a.images
    table, flat = [], []
    for tag in tags:
        prov = out / 'provenance' / tag
        if not (out / 'final_cells' / f'{tag}.parquet').is_file():
            print(f'{tag}: no final cells yet, skipped', file=sys.stderr)
            continue
        exp = (load(prov / 'tissue_export.json').get('step_result') or {}).get('details', {})
        imp = (load(prov / 'tissue_import.json').get('step_result') or {}).get('details', {})
        det = load(prov / '0.json')
        fit = load(out / 'fit' / f'{tag}.json')
        one_b = (load(prov / '1b.json').get('step_result') or {}).get('details', {})
        final = (load(prov / '8.json').get('counts') or {}).get('annotations_by_class', {})
        case, cohort = rows[tag]['case_id'], rows[tag]['cohort']
        cmp_path = out / 'review' / f'{tag}.segmentation_comparison.json'
        old = Path(str(OLD_CENTROIDS).format(case=case))
        cells = out / 'cells' / f'{tag}.parquet'
        if not cmp_path.is_file() and cells.is_file() and old.is_file():
            subprocess.run([sys.executable, str(HERE / 'compare_segmentation.py'), '--new', str(cells), '--old', str(old),
                            '--out', str(cmp_path)], check=True, capture_output=True)
        cmp = load(cmp_path)
        went = where_endocrine_went(out / 'final_cells' / f'{tag}.parquet', case, cohort, out / 'script1_outlines' / f'{tag}.geojson')
        objects = {k: final.get(k, 0) for k in ('Islet', 'EndocrineCluster', 'EndocrineSingle', 'Hormone Fragment')}
        entry = {
            'image': tag, 'case_id': case, 'status': rows[tag]['donor_status'], 'cohort': cohort,
            'tissue_mm2': round(imp.get('tissue_area_mm2', float('nan')), 3),
            'cells_new': (det.get('counts') or {}).get('detections'), 'cells_previous': cmp.get('old_cells'),
            'previous_cells_matched_pct': cmp.get('old_matched_percent'), 'new_cells_matched_pct': cmp.get('new_matched_percent'),
            'matched_area_ratio_median': (cmp.get('matched_area_ratio_new_over_old') or {}).get('median'),
            'fit_thresholds': fit.get('thresholds'), 'fit_floor_applied': (fit.get('cohort_scale') or {}).get('floor_applied'),
            'fit_endocrine_cells': (fit.get('cells') or {}).get('endocrine'),
            'fit_confirmed_seeds': (fit.get('cells') or {}).get('confirmed'),
            'script1_outlines': one_b.get('islets_from_script_1'), 'objects': objects,
            'whole_cells_in_objects': one_b.get('whole_cells'),
            'source_islet_objects': {k: v for k, v in exp.get('all_annotation_classes', {}).items() if 'Islet' in k or 'Endocrine' in k},
            'earlier_endocrine': went,
        }
        table.append(entry)
        flat.append({'image': tag, 'status': entry['status'], 'cohort': cohort, 'cells_new': entry['cells_new'],
                     **{f'objects_{k}': v for k, v in objects.items()},
                     **{f'thr_{k}': (round(v) if v else None) for k, v in (entry['fit_thresholds'] or {}).items()},
                     'floor_applied': ','.join(entry['fit_floor_applied'] or []),
                     **{f'{ph.split()[0].lower()}_inside_pct': went[ph]['inside_pct'] for ph in ENDOCRINE},
                     **{f'{ph.split()[0].lower()}_in_script1_pct': went[ph].get('in_script1_outlines_pct') for ph in ENDOCRINE},
                     **{f'{ph.split()[0].lower()}_n': went[ph]['n'] for ph in ENDOCRINE},
                     'earlier_islet_core_inside_pct': went['earlier_islet_core_cells']['inside_pct'],
                     'earlier_islets': went['earlier_islet_core_cells']['islets']})
        print(f"{tag}: objects {objects}; earlier beta/alpha/delta inside the objects "
              f"{went['Beta cell']['inside_pct']} / {went['Alpha cell']['inside_pct']} / {went['Delta cell']['inside_pct']} %, "
              f"inside script 1's outlines {went['Beta cell'].get('in_script1_outlines_pct')} / "
              f"{went['Alpha cell'].get('in_script1_outlines_pct')} / {went['Delta cell'].get('in_script1_outlines_pct')} %", flush=True)
    (out / 'review').mkdir(exist_ok=True)
    (out / 'review' / 'summary.json').write_text(json.dumps(table, indent=2) + '\n')
    pd.DataFrame(flat).to_csv(out / 'review' / 'summary.csv', index=False)
    print(f"wrote {out / 'review' / 'summary.json'} and summary.csv ({len(table)} images)")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
