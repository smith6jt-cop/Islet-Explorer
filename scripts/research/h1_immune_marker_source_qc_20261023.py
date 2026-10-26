#!/usr/bin/env python3
"""Reproduce marker coverage checks and the source explanation for cytoplasm NA."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import h1_features_20261022 as hf
from h1_immune_marker_features_20261023 import MARKERS, IMMUNE_MARKERS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feature-root', type=Path, required=True)
    parser.add_argument('--uniform-source-root', type=Path, required=True)
    args = parser.parse_args()
    root, source = args.feature_root.resolve(strict=True), args.uniform_source_root.resolve(strict=True)
    d = pd.read_csv(root / 'object_marker_counts.csv')
    denominators = pd.read_csv(root / 'population_denominators.csv')
    profiles = pd.read_parquet(root / 'marker_profiles.parquet')
    assert d.islet_id.nunique() == 4763 and d.case_id.nunique() == 20 and 78 not in set(d.case_id)
    assert len(d) == 4763 * 2 * (31 + 2 * 12) and len(denominators) == 4763 * 2 * 8
    assert set(d.loc[d.threshold.eq(.1), 'marker']) == set(MARKERS)
    assert set(d.loc[d.threshold.eq(.05), 'marker']) == set(IMMUNE_MARKERS)
    assert set(d.loc[d.threshold.eq(.20), 'marker']) == set(IMMUNE_MARKERS)
    assert d.n_immune.le(d.n_all).all()
    assert (d.n_immune_positive.le(d.n_immune_measured) | d.n_immune_positive.isna()).all()
    assert (d.n_all_positive.le(d.n_all_measured) | d.n_all_positive.isna()).all()
    assert d.loc[d.n_immune.eq(0), 'n_immune_positive'].eq(0).all()
    assert d.loc[d.n_immune.gt(0) & d.n_immune_measured.eq(0), 'n_immune_positive'].isna().all()
    assert profiles.n_cells.gt(0).all() and profiles.n_measured.le(profiles.n_cells).all()
    for population in ['all', 'immune']:
        p = d.loc[d.marker.isin(IMMUNE_MARKERS)].pivot(index=['islet_id', 'compartment', 'marker'],
                                                      columns='threshold', values=f'n_{population}_positive')
        assert p[.05].isna().eq(p[.10].isna()).all() and p[.20].isna().eq(p[.10].isna()).all()
        assert ((p[.05] >= p[.10]) | p[.10].isna()).all()
        assert ((p[.10] >= p[.20]) | p[.10].isna()).all()
    x = d.loc[d.threshold.eq(.1)].copy()
    x['all_missing'] = x.n_all - x.n_all_measured
    x['immune_missing'] = x.n_immune - x.n_immune_measured
    x['immune_unmeasured_pool'] = x.n_immune.gt(0) & x.n_immune_measured.eq(0)
    x['immune_partly_measured_pool'] = x.n_immune_measured.gt(0) & x.n_immune_measured.lt(x.n_immune)
    coverage = x.groupby(['case_id', 'donor_status', 'compartment', 'marker'], observed=True).agg(
        objects=('islet_id', 'size'), all_cell_memberships=('n_all', 'sum'),
        all_measured_memberships=('n_all_measured', 'sum'), all_missing_memberships=('all_missing', 'sum'),
        immune_cell_memberships=('n_immune', 'sum'), immune_measured_memberships=('n_immune_measured', 'sum'),
        immune_missing_memberships=('immune_missing', 'sum'),
        immune_fully_unmeasured_objects=('immune_unmeasured_pool', 'sum'),
        immune_partly_measured_objects=('immune_partly_measured_pool', 'sum')).reset_index()
    coverage.to_csv(root / 'marker_coverage_by_donor.csv', index=False)
    ki67_path = root / 'ki67' / 'ki67_immune_cells.parquet'
    immune = pd.read_parquet(ki67_path, columns=['case_id', 'object_id', 'phenotype', 'image_name',
                              'Cytoplasm__n_px', 'qc_flags', 'paired_core', 'membership_multiplicity'])
    assert not immune.duplicated(['case_id', 'object_id']).any()
    frames, hashes = [], {str(ki67_path): hf.sha256(ki67_path)}
    for case, g in immune.groupby('case_id'):
        assert case != 78 and g.image_name.nunique() == 1
        path = source / 'typing' / 'aligned' / (g.image_name.iloc[0] + '.parquet')
        signals = pd.read_parquet(path, columns=['object_id', 'LGALS3', 'VIM'])
        frame = g.merge(signals, on='object_id', how='left', validate='one_to_one', indicator=True)
        assert len(frame) == len(g)
        assert frame._merge.eq('both').all()
        frame = frame.drop(columns='_merge')
        frame['LGALS3_missing'] = ~np.isfinite(frame.LGALS3)
        frame['VIM_missing'] = ~np.isfinite(frame.VIM)
        frames.append(frame)
        hashes[str(path)] = hf.sha256(path)
    immune = pd.concat(frames, ignore_index=True)
    missing = immune.loc[immune.LGALS3_missing].copy()
    assert immune.LGALS3_missing.eq(immune.VIM_missing).all()
    assert missing['Cytoplasm__n_px'].eq(0).all()
    missing.to_parquet(root / 'cytoplasm_missing_immune_cells.parquet', index=False)
    by_donor = immune.groupby('case_id').agg(immune_unique_cells=('object_id', 'size'),
        LGALS3_missing_unique=('LGALS3_missing', 'sum'), VIM_missing_unique=('VIM_missing', 'sum'),
        no_cytoplasm_pixels=('Cytoplasm__n_px', lambda x: (x == 0).sum())).reset_index()
    by_donor.to_csv(root / 'cytoplasm_missingness_by_donor.csv', index=False)
    local = x.loc[x.marker.eq('LGALS3')]
    report = dict(universe='all current curated immune cells in served donor file, persistent IDs joined to alignedLGALS3/VIM and sourcegeometryQC',
                  immune_unique_cells=len(immune), LGALS3_missing_unique=len(missing),
                  missing_with_zero_cytoplasm_pixels=int(missing['Cytoplasm__n_px'].eq(0).sum()),
                  missing_with_nonzero_cytoplasm_pixels=int(missing['Cytoplasm__n_px'].gt(0).sum()),
                  missing_qc_flags=missing.qc_flags.value_counts().to_dict(), LGALS3_and_VIM_missing_agree=True,
                  paired_core_missing_unique=int(missing.paired_core.sum()),
                  paired_missing_immune_memberships=int(local.immune_missing.sum()),
                  paired_fully_unmeasured_immune_pools=int(local.immune_unmeasured_pool.sum()),
                  paired_partly_measured_immune_pools=int(local.immune_partly_measured_pool.sum()),
                  source_reason='cytoplasmic measurement undefined when no cytoplasm raster pixels remain; eligible cells and membrane/nuclear measurements retained; no imputation',
                  source_sha256=hashes, script_sha256=hf.sha256(Path(__file__)))
    hf.write_json(root / 'cytoplasm_missingness_audit.json', report)
    hf.write_json(root / 'marker_validation_summary.json', dict(objects=4763, donors=20,
        compact_rows=len(d), population_denominator_rows=len(denominators), primary_nonempty_profile_rows=len(profiles),
        measurement_aware_count_checks_pass=True, threshold_monotonicity_pass=True,
        source_identity_reason_for_cytoplasm_missing_verified=True,
        source_sha256={str(root / n): hf.sha256(root / n) for n in ['object_marker_counts.csv', 'population_denominators.csv', 'marker_profiles.parquet']},
        script_sha256=hf.sha256(Path(__file__))))
    print(f'Coverage/count/threshold/source-QC checks passed; {len(missing)} unique immune cells lack cytoplasm measurements.')


if __name__ == '__main__':
    main()
