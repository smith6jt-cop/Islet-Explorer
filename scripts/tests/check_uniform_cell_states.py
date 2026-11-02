#!/usr/bin/env python3
"""Exercise state preservation with existing real typing and one real donor export.

No scientific fixtures are generated and no input is modified. The fresh output
is a diagnostic test fragment, not a complete app root or accepted candidate.
"""
import argparse
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import assemble_uniform as assembly


def rejected(callback, text):
    try:
        callback()
    except ValueError as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError(f'Expected rejection: {text}')


def validate_assembled(args):
    """Read every completed cell export; never invoke geometry or assembly."""
    root = args.assembled_root.resolve(strict=True)
    reg = assembly.RG.load()
    fields = [f'{marker}_pos' for marker in reg.markers('state')]
    rows = assembly.CT.cohort_rows()
    levels = pd.read_csv(args.run / 'typing/levels.csv')
    manifest = json.loads((root / 'uniform_release.json').read_text())
    assert manifest['donors'] == len(rows) == 20
    expected_states, reports = {}, {}
    distance_columns = ['case_id', 'phenotype', 'cell_region', *fields]
    distance = pd.read_csv(root / 'cell_distance.csv', usecols=distance_columns,
                           dtype={'case_id': str}, keep_default_na=False)
    assert set(distance.case_id) == {str(int(row['case_id'])) for row in rows}
    expected_membership = {f"{int(row['case_id'])}.parquet" for row in rows}
    assert {p.name for p in (root / 'cell_types').glob('*.parquet')} == expected_membership
    assert {p.stem for p in (root / 'donors').glob('*.csv')} == {Path(p).stem for p in expected_membership}
    lookup = pd.read_csv(root / 'islet_spatial_lookup.csv', dtype={'case_id': str})
    # The lookup retains unpaired geometry objects. Only paired analysis rows in
    # the assembled H5AD have per-object CSVs (e.g. frozen v20: 4,853 vs 4,848).
    adata = ad.read_h5ad(root / 'islet_explorer.h5ad', backed='r')
    assert list(adata.var_names) == assembly.MARKERS and len(adata.var_names) == 31
    object_ids = adata.obs['islet_id'].astype(str)
    assert not object_ids.isna().any() and not object_ids.duplicated().any()
    analysis_objects = set(object_ids)
    adata.file.close()
    assert analysis_objects <= set(lookup.combined_islet_id)
    assert not lookup.combined_islet_id.duplicated().any()
    expected_objects = {value + '.csv' for value in analysis_objects}
    assert {p.name for p in (root / 'cells').glob('*.csv')} == expected_objects
    total_object_rows = 0
    for row in rows:
        tag, case = row['image_name'], str(int(row['case_id']))
        lv = levels[levels.image == tag].set_index('marker')
        expected_states[tag] = {'cell_states': {
            f'{marker}_pos': {'marker': marker, 'measurement': str(lv.loc[marker, 'measurement']),
                             'redsea_column': reg.column(marker, str(lv.loc[marker, 'measurement'])),
                             'level_log': float(lv.loc[marker, 'level_log'])}
            for marker in reg.markers('state')}}
        source = pd.read_parquet(args.run / 'typing/identity' / f'{tag}.parquet',
                                 columns=['object_id', 'label', *fields])
        source = source[source.label != 'Excluded'].set_index('object_id')
        assert not source.index.has_duplicates and not source.index.isna().any()
        member = pd.read_parquet(root / 'cell_types' / f'{case}.parquet')
        assert not member.object_id.duplicated().any() and not member.object_id.isna().any()
        assert set(member.object_id) == set(source.index)
        reference = source.loc[member.object_id]
        np.testing.assert_array_equal(member.phenotype, reference.label)
        for field in fields:
            assert pd.api.types.is_numeric_dtype(member[field]) and member[field].isin([0, 1]).all()
            np.testing.assert_array_equal(member[field], reference[field])

        key = ['X_centroid', 'Y_centroid', 'phenotype']
        coordinate_reference = member[['x_um', 'y_um', 'phenotype', *fields]].rename(
            columns={'x_um': 'X_centroid', 'y_um': 'Y_centroid'}).drop_duplicates()
        assert not coordinate_reference.duplicated(key).any(), f'{tag}: ambiguous coordinate/phenotype state key'
        def check_coordinates(frame):
            expected = frame[key].merge(coordinate_reference, on=key, how='left',
                                        validate='many_to_one', indicator=True)
            assert expected['_merge'].eq('both').all(), f'{tag}: export contains unknown cell coordinate/phenotype'
            for field in fields:
                assert pd.api.types.is_numeric_dtype(frame[field]) and frame[field].isin([0, 1]).all()
                np.testing.assert_array_equal(frame[field], expected[field])

        donor = pd.read_csv(root / 'donors' / f'{case}.csv', float_precision='round_trip', keep_default_na=False)
        assert len(donor) == len(member)
        check_coordinates(donor)
        # Match the complete donor key multiset as well as values, so duplicating
        # one coordinate cannot silently replace another cell with the same state.
        full_key = [*key, 'cell_region', 'islet_name', *fields]
        expected_donor = member.rename(columns={'x_um': 'X_centroid', 'y_um': 'Y_centroid'})[full_key]
        pd.testing.assert_frame_equal(
            donor[full_key].sort_values(full_key).reset_index(drop=True),
            expected_donor.sort_values(full_key).reset_index(drop=True), check_dtype=False)
        per_object = lookup[(lookup.case_id.astype(int) == int(case)) &
                            lookup.combined_islet_id.isin(analysis_objects)]
        object_rows = 0
        for name in per_object.combined_islet_id:
            frame = pd.read_csv(root / 'cells' / f'{name}.csv',
                                usecols=[*key, 'cell_region', *fields], float_precision='round_trip')
            assert frame.cell_region.isin(['core', 'peri']).all()
            check_coordinates(frame)
            object_rows += len(frame)
        # cell_distance has neither UUID nor coordinates. Assembly writes each
        # donor in the exact cell_types UUID order, then concatenates cohort order.
        # Verify this documented row association, including phenotype and region.
        actual_distance = distance.loc[distance.case_id == case, ['phenotype', 'cell_region', *fields]]
        pd.testing.assert_frame_equal(actual_distance.reset_index(drop=True),
                                      member[['phenotype', 'cell_region', *fields]].reset_index(drop=True),
                                      check_dtype=False)
        reports[tag] = {'cells': len(member), 'positive': {field: int(member[field].sum()) for field in fields},
                        'object_files': len(per_object), 'object_rows': object_rows,
                        'uuid_flags_and_labels_exact': True, 'coordinate_export_flags_exact': True,
                        'distance_row_associations_exact': True}
        total_object_rows += object_rows
    assert manifest['cell_states'] == assembly.cell_state_metadata(expected_states)
    report = {'result': 'pass', 'scope': 'All-20 completed cell-state exports; no assembly or classification run',
              'scientific_accuracy_established': False, 'run': str(args.run.resolve()),
              'assembled_root': str(root), 'images': reports,
              'cells': sum(record['cells'] for record in reports.values()),
              'object_files': len(expected_objects), 'object_rows': total_object_rows,
              'fields': fields, 'measurement_and_level_metadata_exact': True,
              'embedding_marker_count': 31,
              'distance_row_contract': 'Within each donor, cell_types UUID order; phenotype/region/state matched exactly'}
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('result', 'cells', 'object_files', 'object_rows')}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--donor-source', type=Path)
    parser.add_argument('--assembled-root', type=Path,
                        help='validate every existing export; never rerun build_image')
    parser.add_argument('--image', default='ND_0112')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--validate-only', action='store_true',
                        help='check existing state rows without generating a donor assembly (workstation-safe)')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('use a fresh output directory')
    if args.assembled_root and args.validate_only:
        parser.error('--assembled-root and --validate-only select different validation scopes')
    if not args.assembled_root and not args.validate_only and args.donor_source is None:
        parser.error('--donor-source is required to generate a donor assembly fragment')
    args.out.mkdir(parents=True)
    if args.assembled_root:
        validate_assembled(args)
        return
    reg = assembly.RG.load()
    rows = assembly.CT.cohort_rows()
    levels = pd.read_csv(args.run / 'typing/levels.csv')
    state_fields = [f'{marker}_pos' for marker in reg.markers('state')]
    per_image, results = {}, {}
    selected = None
    for row in rows:
        tag = row['image_name']
        ident = pd.read_parquet(args.run / 'typing/identity' / f'{tag}.parquet',
                                columns=['object_id', *state_fields])
        lv = levels[levels.image == tag].set_index('marker')
        columns = [reg.column(marker, lv.loc[marker, 'measurement']) for marker in reg.markers('state')]
        rs = pd.read_parquet(args.run / 'typing/redsea' / f'{tag}.parquet',
                            columns=['object_id', 'qc_analysis_eligible', *columns])
        rs = rs[rs.qc_analysis_eligible].reset_index(drop=True)
        values, provenance = assembly.cell_states(ident, rs, lv, reg, tag)
        # Changing file row order must not change a UUID's state.
        shuffled, _ = assembly.cell_states(ident, rs.iloc[::-1], lv, reg, tag)
        pd.testing.assert_frame_equal(values, shuffled)
        per_image[tag] = {'cells': len(values), 'positive': {field: int(values[field].sum()) for field in state_fields}}
        results[tag] = {'cell_states': provenance}
        if tag == args.image:
            selected = (row, ident, rs, lv)
    assert selected is not None, args.image
    row, ident, rs, lv = selected
    flag = state_fields[0]
    rejected(lambda: assembly.cell_states(ident.drop(columns=flag), rs, lv, reg, args.image), 'missing or nonbinary')
    altered = ident.copy()
    altered.loc[altered.index[0], flag] = 1 - altered.loc[altered.index[0], flag]
    rejected(lambda: assembly.cell_states(altered, rs, lv, reg, args.image), 'differs from identity')
    altered.loc[altered.index[0], flag] = 2
    rejected(lambda: assembly.cell_states(altered, rs, lv, reg, args.image), 'missing or nonbinary')
    rejected(lambda: assembly.cell_states(ident, rs.iloc[:-1], lv, reg, args.image), 'UUIDs differ')
    rejected(lambda: assembly.cell_states(pd.concat([ident, ident.iloc[:1]]), rs, lv, reg, args.image), 'invalid identity UUIDs')

    report = {'result': 'pass', 'purpose': 'Real-data state preservation, not scientific acceptance',
              'run': str(args.run.resolve()), 'image': args.image, 'all_images': per_image,
              'state_contract': assembly.cell_state_metadata(results),
              'uuid_reordering_invariant': True,
              'corruption_checks': ['missing flag', 'flipped flag', 'nonbinary flag', 'missing UUID', 'duplicate UUID']}
    if args.validate_only:
        report['scope'] = 'Existing state validation only; assembly/export integration not run'
        (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'result': 'pass', 'images': len(per_image),
                          'cells': sum(record['cells'] for record in per_image.values()), 'assembly_run': False}))
        return

    output = args.out / 'assembly-fragment'
    for directory in ('cells', 'donors', 'cell_types', 'gson'):
        (output / directory).mkdir(parents=True)
    source = ad.read_h5ad(args.donor_source / 'islet_explorer.h5ad', backed='r')
    donors = source.obs.drop_duplicates('imageid')[
        ['imageid', 'age', 'gender', 'donor_status', 'GADA', 'ZnT8A', 'IA2A', 'mIAA']].copy()
    source.file.close()
    donors = assembly.apply_donor_metadata(donors, pd.read_csv(assembly.HERE / 'donor_metadata.csv',
                                                               dtype={'case_id': str}))
    groups = assembly.label_groups(args.run, rows)
    result = assembly.build_image(args.run, output, args.image, row['case_id'], donors, groups)
    result['distance'].to_csv(output / 'cell_distance.csv', index=False)
    metadata = assembly.cell_state_metadata({args.image: result})
    (output / 'uniform_release.json').write_text(json.dumps({
        'purpose': 'Real-data cell-state export test fragment; not an app release',
        'cell_states': metadata}, indent=2) + '\n')
    case = str(int(row['case_id']))
    membership = pd.read_parquet(output / 'cell_types' / f'{case}.parquet')
    original = pd.read_parquet(args.run / 'typing/identity' / f'{args.image}.parquet').set_index('object_id')
    for field in state_fields:
        np.testing.assert_array_equal(membership[field], original.loc[membership.object_id, field])
    np.testing.assert_array_equal(membership.phenotype, original.loc[membership.object_id, 'label'])
    assert list(result['adata'].var_names) == assembly.MARKERS
    assert len(assembly.MARKERS) == 31 and 'Ki67' not in assembly.MARKERS
    donor = pd.read_csv(output / 'donors' / f'{case}.csv')
    for field in state_fields:
        np.testing.assert_array_equal(donor[field], membership[field])
        np.testing.assert_array_equal(result['distance'][field], membership[field])
    coordinates = membership[['x_um', 'y_um', *state_fields]].rename(
        columns={'x_um': 'X_centroid', 'y_um': 'Y_centroid'})
    assert not coordinates.duplicated(['X_centroid', 'Y_centroid']).any()
    object_files = sorted((output / 'cells').glob('*.csv'))
    assert object_files
    object_rows = 0
    for file in object_files:
        frame = pd.read_csv(file, float_precision='round_trip')
        expected = frame[['X_centroid', 'Y_centroid']].merge(
            coordinates, on=['X_centroid', 'Y_centroid'], how='left', validate='many_to_one')
        for field in state_fields:
            np.testing.assert_array_equal(frame[field], expected[field])
        object_rows += len(frame)
    report.update(assembly_fragment=str(output.resolve()), donor_rows=len(donor),
                  object_files=len(object_files), object_rows=object_rows,
                  identity_labels_unchanged=True, embedding_marker_count=len(assembly.MARKERS))
    (args.out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('result', 'image', 'donor_rows', 'object_files', 'object_rows')}))


if __name__ == '__main__':
    main()
