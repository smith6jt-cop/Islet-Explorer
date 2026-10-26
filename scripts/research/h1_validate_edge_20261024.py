#!/usr/bin/env python3
"""Independently reconcile corrected H1 identities, memberships and exclusions."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from h1_features_20261022 import IMMUNE, LABELS, sha256, write_json
from h1_analysis_20261022 import eligible, fit_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, default=Path('output/h1_20261022'))
    parser.add_argument('--snapshot-only', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    if args.snapshot_only:
        out.mkdir(parents=True, exist_ok=True)
        path = out / 'historical_checksums.json'
        if not path.exists():
            write_json(path, {str(p.resolve()): sha256(p)
                for root in [args.baseline, Path('output/h1_immune_size_20261023')]
                for p in root.rglob('*') if p.is_file()})
        return
    audit = json.loads((out / 'feature_audit.json').read_text())
    app = Path(audit['data_root'])
    df = pd.read_parquet(out / 'h1_object_features.parquet')
    old = pd.read_csv(args.baseline / 'h1_object_features.csv', dtype={'case_id': str})
    primary, previous = eligible(df), eligible(old)
    assert set(primary.islet_id).issubset(previous.islet_id)
    assert len(df) == len(old) == 4763 and primary.case_id.nunique() == 20
    assert audit['count_discrepancies'] == audit['core_geometry_identity_discrepancies'] == 0
    census, flows = [], []
    for case, features in df.groupby('case_id'):
        cells = pd.read_parquet(out / 'analysis_cells_by_donor' / f'{case}.parquet')
        source = pd.read_parquet(app / 'cell_types' / f'{case}.parquet')
        assert len(cells) == len(source) and not cells.object_id.duplicated().any()
        assert set(cells.object_id) == set(source.object_id)
        cols = ['phenotype', 'x_um', 'y_um']
        pd.testing.assert_frame_equal(cells.set_index('object_id')[cols].sort_index(),
                                      source.set_index('object_id')[cols].sort_index(), check_dtype=False)
        assert cells.phenotype.isin(LABELS).all() and cells.phenotype.notna().all()
        relation = pd.read_parquet(out / 'object_memberships_by_donor' / f'{case}.parquet')
        assert not relation.duplicated(['case_id', 'object_id', 'islet_id']).any()
        assert set(relation.object_id).issubset(cells.loc[cells.edge_analysis_eligible, 'object_id'])
        exclusive = relation[relation.core_complete & relation.exclusive_owner]
        assert not exclusive.object_id.duplicated().any()
        primary_ids = set(cells.loc[cells.primary_compartment.isin(['core', 'peri']), 'object_id'])
        assert primary_ids == set(exclusive.object_id)
        assert not (relation.loc[~relation.core_complete, 'exclusive_owner']).any()
        grouped = relation.assign(immune=relation.phenotype.isin(IMMUNE), beta=relation.phenotype.eq('Beta cell'))
        for region in ['core', 'peri']:
            group = grouped[grouped.cell_region.eq(region)].groupby('islet_id').agg(
                n=('object_id', 'size'), immune=('immune', 'sum'), beta=('beta', 'sum'))
            for col, target in [('n', f'n_cells_{region}'), ('immune', f'n_immune_{region}')]:
                assert np.array_equal(group[col].reindex(features.islet_id, fill_value=0).to_numpy(), features[target])
            if region == 'core':
                assert np.array_equal(group.beta.reindex(features.islet_id, fill_value=0).to_numpy(), features.n_beta_core)
        exc_counts = exclusive[exclusive.cell_region.eq('peri')].assign(immune=lambda x:x.phenotype.isin(IMMUNE)).groupby('islet_id').immune.sum()
        assert np.array_equal(exc_counts.reindex(features.islet_id, fill_value=0).to_numpy(), features.n_immune_peri_exclusive)
        for (compartment, phenotype), n in cells.groupby(['primary_compartment', 'phenotype'], observed=True).size().items():
            census.append(dict(case_id=int(case), compartment=compartment, phenotype=phenotype, count=int(n)))
        base = previous[previous.case_id.eq(str(case))]
        current = primary[primary.case_id.eq(str(case))]
        subset = features[features.islet_id.isin(base.islet_id)]
        core_ok = subset.core_complete
        peri_ok = core_ok & subset.neighborhood_complete_20
        field_ok = peri_ok & subset.neighborhood_complete_100
        flows.append(dict(case_id=int(case), donor_status=features.donor_status.iloc[0],
            source_cells=len(cells), source_paired_objects=len(features),
            baseline_primary=len(base), core_excluded=int((~core_ok).sum()),
            additional_peri_incomplete=int((core_ok & ~peri_ok).sum()),
            additional_field_incomplete=int((peri_ok & ~field_ok).sum()),
            corrected_primary=len(current), corrected_core_immune=int(current.n_immune_core.sum()),
            corrected_core_positive=int(current.has_core_immune.sum())))
    pd.DataFrame(census).to_csv(out / 'exclusive_cell_census.csv', index=False)
    pd.DataFrame(flows).to_csv(out / 'eligibility_flow_by_donor.csv', index=False)
    objects = pd.read_csv(out / 'eligibility/object_eligibility.csv')
    objects['size_bin_um2'] = pd.cut(objects.core_area_um2, [0,100,1000,10000,np.inf], right=False).astype(str)
    objects.groupby(['object_class','size_bin_um2'], observed=True).agg(
        source_objects=('core_complete','size'), complete=('core_complete','sum'),
        neighborhood100_complete=('neighborhood_complete_100','sum')).to_csv(out / 'exclusions_by_class_size.csv')
    rows, _, _ = fit_model(previous, 'Baseline reproduction')
    saved = json.loads((args.baseline / 'analysis_summary.json').read_text())['primary']
    for key in ['rr','lower','upper','p_cluster_t']:
        assert np.isclose(rows[0][key], saved[key], rtol=1e-10, atol=1e-12), key
    sources = dict(audit['source_sha256'])
    for p in (out / 'eligibility/audits').glob('*.json'):
        a = json.loads(p.read_text())
        assert a['version'] == 'h1-edge-20261024-v3'
        sources.update(a['source_sha256'])
    for path, digest in sources.items():
        assert sha256(path) == digest, f'Input changed: {path}'
    historical = json.loads((out / 'historical_checksums.json').read_text())
    for path, digest in historical.items():
        assert sha256(path) == digest, f'Historical result changed: {path}'
    # Replay one predetermined canary using saved thresholds, never refitting.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
    from islet_classifier_fit import read_planes, smooth, classifier_mask
    source_run = Path(json.loads((out / 'eligibility/edge_manifest.json').read_text())['source_run'])
    tag = 'Aab_6505'
    config = json.loads((source_run / 'fit' / f'{tag}.json').read_text())
    planes_path = source_run / 'planes' / f'{tag}.planes.tif'
    planes = read_planes(planes_path)
    smoothed = {c: smooth(v, config['sigma_um'] / config['resolution_um']) for c, v in planes.items()}
    replay = classifier_mask(smoothed, config['thresholds'],
        **{k: config[k] for k in ['seed_px', 'close_px', 'open_px', 'margin_px']})
    mask_path = source_run / 'fit' / f'{tag}.mask.npy'
    different = int(np.count_nonzero(replay != np.load(mask_path)))
    assert different == 0
    write_json(out / 'preclip_replay_canary.json', dict(image=tag, different_pixels=different,
        thresholds_refitted=False, planes_sha256=sha256(planes_path), saved_mask_sha256=sha256(mask_path)))
    totals = pd.DataFrame(flows).select_dtypes('number').sum().to_dict()
    totals.pop('case_id')
    write_json(out / 'validation_summary.json', dict(status='pass', totals=totals,
        distinct_physical_cells=int(pd.DataFrame(census)['count'].sum()),
        donor_count=20, original_counts_reproduced=True, phenotype_coordinate_identity_preserved=True,
        source_hashes_verified=len(sources), historical_hashes_verified=len(historical),
        baseline_primary_reproduced=rows[0], corrected_primary_objects=len(primary),
        checks=['UUID uniqueness and source phenotype/coordinate equality', 'one compartment per physical cell',
                'one exclusive owner per eligible core/peri cell', 'unique cells within every focal object union',
                'core and peri Beta/immune/total counts independently reproduced from membership relations',
                'source H5AD and all paired drilldown counts reproduced before exclusions',
                'complete core and acquisition neighborhoods required for H1',
                'old primary estimate reproduced and historical/source hashes unchanged']))
    print(json.dumps({'status':'pass', 'totals':totals}, indent=2))


if __name__ == '__main__':
    main()
