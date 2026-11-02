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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--donor-source', type=Path, required=True)
    parser.add_argument('--image', default='ND_0112')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--validate-only', action='store_true',
                        help='check existing state rows without generating a donor assembly (workstation-safe)')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('use a fresh output directory')
    args.out.mkdir(parents=True)
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
