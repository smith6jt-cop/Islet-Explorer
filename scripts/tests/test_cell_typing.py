"""Scientific contracts for disease-blind cell typing and biological annotation."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/islet-typing-test-numba')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-typing-test-mpl')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
from common import MARKERS, sha256, write_json
from cell_typing import (TypingConfig, balanced_sample, fit_transform, transform,
    consensus_assignments, read_annotations, catalog_from_annotations,
    discover, apply, attach_positivity, merge_refinement, configure_runtime,
    ENDOCRINE_EVIDENCE, REQUIRED_REFINEMENT_MARKERS, merge_supplemental,
    segmentation_quality_mask, attach_quality_status, cell_eligibility)
from cell_typing_validation import validation_variants, stability_metrics, align_training_lineages
from merge_typing_refinements import merge as merge_typing_deliveries


def cells(cid, n=64, seed=1):
    rng = np.random.default_rng(seed)
    value = rng.lognormal(0, .2, (n, len(MARKERS)))
    value[:n // 2, :8] += 20
    value[n // 2:, 8:20] += 20
    result = pd.DataFrame(value, columns=MARKERS)
    result.insert(0, 'object_id', [f'cell-{i:04d}' for i in range(n)])
    result.insert(0, 'case_id', str(cid))
    result['x_um'] = np.arange(n, dtype=float)
    result['y_um'] = np.arange(n, dtype=float) + 1
    result['Cell Area'] = 25.
    return result


def projections(labels, supported=None):
    supported = [True] * len(labels) if supported is None else supported
    return pd.DataFrame({'cluster_id': labels, 'projection_supported': supported,
                         'neighbor_agreement': [1.] * len(labels)})


def annotation(base_type='Beta cell', decision='accepted', state=''):
    return {'base_type': base_type, 'lineage': 'Endocrine', 'decision': decision,
            'phenotype_state': state, 'endocrine_evidence': ['INS'], 'review_note': 'reviewed',
            'evidence': []}


class CellTypingContracts(unittest.TestCase):
    def test_two_parent_pools_agreeing_on_biology_are_not_a_conflict(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        mm = {(f, 0): {**annotation(''), 'annotation_level': 'lineage',
              'lineage': 'Stromal' if f in ['a', 'b'] else 'Vascular'} for f in pp}
        parent = consensus_assignments(keys, pp, mm)
        mapping = {(f, 0): {**annotation('Lymphatic endothelial cell'), 'lineage': 'Vascular'} for f in pp}
        a = consensus_assignments(keys, pp, mapping)
        b = a.copy()
        a.loc[0, 'confidence'] = .8
        a.at[0, 'phenotype_states'] = ['Lymphatic endothelial cell: VIM high']
        b.at[0, 'phenotype_states'] = ['Lymphatic endothelial cell: CD44 high']
        first = merge_refinement(merge_refinement(parent, a, 'Stromal'), b, 'Vascular')
        second = merge_refinement(merge_refinement(parent, b, 'Vascular'), a, 'Stromal')
        self.assertEqual(first.classification_status.iloc[0], 'validated')
        self.assertEqual(first.lineage.iloc[0], 'Vascular')
        self.assertEqual(first.confidence.iloc[0], 1.)
        self.assertEqual(first.phenotype_states.iloc[0], second.phenotype_states.iloc[0])
        self.assertEqual(len(first.phenotype_states.iloc[0]), 2)

    def test_separate_refinement_deliveries_preserve_conflicts_and_parent_binding(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            parent_root = root / 'broad'
            (parent_root / 'typing/cells').mkdir(parents=True)
            keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
            pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
            mm = {}
            for fit in pp:
                row = annotation('')
                row.update(annotation_level='lineage', lineage='Endocrine' if fit in ['a', 'b'] else 'Neural')
                mm[(fit, 0)] = row
            parent = consensus_assignments(keys, pp, mm)
            parent.to_parquet(parent_root / 'typing/cells/a.parquet', index=False)
            write_json(parent_root / 'typing/catalog.json', catalog_from_annotations({}))
            write_json(parent_root / 'typing/application.json', {'status': 'complete'})
            children = []
            for lineage, base in [('Endocrine', 'Beta cell'), ('Neural', 'Neural-associated cell')]:
                child_root = root / lineage
                (child_root / 'typing/cells').mkdir(parents=True)
                model = child_root / 'model'
                model.mkdir()
                write_json(model / 'specification.json', {'config': {'refinement': lineage}, 'identity': lineage})
                mapping = {(fit, 0): {**annotation(base), 'lineage': lineage} for fit in pp}
                child = consensus_assignments(keys, pp, mapping)
                child = merge_refinement(parent, child, lineage)
                child.to_parquet(child_root / 'typing/cells/a.parquet', index=False)
                write_json(child_root / 'typing/catalog.json', catalog_from_annotations(mapping))
                write_json(child_root / 'typing/application.json', {'status': 'complete', 'model_root': str(model),
                    'parent_typing_hashes': {'a': sha256(parent_root / 'typing/cells/a.parquet')}})
                children.append(child_root)
            first = merge_typing_deliveries(parent_root, children, root / 'merged')
            second = merge_typing_deliveries(parent_root, list(reversed(children)), root / 'reverse')
            self.assertEqual(first['identity'], second['identity'])
            output = pd.read_parquet(root / 'merged/typing/cells/a.parquet')
            self.assertEqual(output.classification_status.iloc[0], 'ambiguous_lineage_refinement')
            self.assertTrue(pd.isna(output.is_endocrine.iloc[0]))
            self.assertEqual(set(output.refinement_model_identities.iloc[0]), {'Endocrine', 'Neural'})
            parent.loc[0, 'confidence'] = .01
            parent.to_parquet(parent_root / 'typing/cells/a.parquet', index=False)
            with self.assertRaisesRegex(ValueError, 'another parent'):
                merge_typing_deliveries(parent_root, children, root / 'bad')

    def test_hierarchical_broad_alignment_uses_training_donors_and_not_final_identities(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            reference = root / 'reference/typing/cells'
            reference.mkdir(parents=True)
            train = pd.DataFrame({'case_id': ['a'] * 12, 'object_id': [str(i) for i in range(12)],
                'candidate_lineages': ['["Endocrine"]'] * 12,
                'base_type': ['Deliberately irrelevant final identity'] * 12})
            train.to_parquet(reference / 'a.parquet', index=False)
            # A corrupt held-out file proves the alignment does not even read it.
            (reference / 'held_out.parquet').write_text('must not be read')
            members = train[['case_id', 'object_id']].copy()
            members['cluster_id'] = 0
            model = root / 'model/arcsinh_k15'
            model.mkdir(parents=True)
            members.to_parquet(model / 'memberships.parquet', index=False)
            spec = {'training_inputs': {'a': 'hash'}, 'config': {'transformations': ['arcsinh'], 'neighbor_counts': [15]}}
            mapping, audit = align_training_lineages(root / 'model', root / 'reference', spec)
            self.assertEqual(mapping[('arcsinh_k15', 0)]['lineage'], 'Endocrine')
            self.assertEqual(mapping[('arcsinh_k15', 0)]['annotation_level'], 'lineage')
            self.assertEqual(mapping[('arcsinh_k15', 0)]['base_type'], '')
            self.assertEqual(audit.n_training_members.iloc[0], 12)

    def test_technical_exclusion_uses_source_pixels_and_never_disease_or_islet_labels(self):
        a = cells('a', 6)
        a['x_source_px'] = [0, 1, 2, 3, 4, 5]
        a['y_source_px'] = 2.
        polygon = {'type': 'Polygon', 'coordinates': [[[1, 1], [3, 1], [3, 3], [1, 3], [1, 1]]]}
        regions = ({'case_id': 'a', 'geometry': polygon, 'evidence': [{'path': 'test', 'sha256': 'test'}]},)
        e = cell_eligibility(a, 5., regions)
        self.assertEqual(e.technical_region_excluded.tolist(), [False, True, True, True, False, False])
        sampled, _ = balanced_sample({'a': a}, 10, technical_exclusion_regions=regions)
        self.assertEqual(sampled.x_source_px.tolist(), [0, 4, 5])
        a['donor_status'] = 'changed'
        a['current_islet_id'] = 'unrelated'
        self.assertEqual(e.typing_eligible.tolist(), cell_eligibility(a, 5., regions).typing_eligible.tolist())

    def test_subcellular_fragments_never_enter_landmarks_or_endocrine_counts(self):
        a = cells('a', 10)
        a.loc[:2, 'Cell Area'] = [.1, 4.99, np.nan]
        sampled, _ = balanced_sample({'a': a}, 10)
        self.assertEqual(len(sampled), 7)
        pp = {f: projections([0] * len(a)) for f in ['a', 'b', 'c', 'd']}
        typed = consensus_assignments(a[['case_id', 'object_id']], pp, {(f, 0): annotation() for f in pp})
        typed = attach_quality_status(typed, a)
        self.assertEqual(typed.classification_status.iloc[:3].tolist(), ['segmentation_quality_excluded'] * 3)
        self.assertTrue(typed.is_endocrine.iloc[:3].isna().all())
        self.assertTrue(typed.base_type.iloc[:3].isna().all())
        self.assertTrue(typed.is_endocrine.iloc[3:].all())

    def test_intermediate_endocrine_lineage_never_assigns_final_hormone_type(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        row = annotation(base_type='')
        row['annotation_level'] = 'lineage'
        typed = consensus_assignments(keys, pp, {(f, 0): row for f in pp})
        self.assertEqual(typed.classification_status.iloc[0], 'lineage_supported')
        self.assertEqual(typed.lineage.iloc[0], 'Endocrine')
        self.assertTrue(pd.isna(typed.base_type.iloc[0]))
        self.assertTrue(pd.isna(typed.is_endocrine.iloc[0]))
        child = consensus_assignments(keys, pp, {(f, 0): annotation() for f in pp})
        refined = merge_refinement(typed, child, 'Endocrine')
        self.assertEqual(refined.classification_status.iloc[0], 'validated')
        self.assertTrue(refined.is_endocrine.iloc[0])

    def test_landmark_sampling_balances_donors_and_ignores_order_labels(self):
        a, b = cells('a', 80), cells('b', 51)
        a['donor_status'], b['donor_status'] = 'T1D', 'ND'
        first, report = balanced_sample({'a': a, 'b': b}, 60)
        self.assertEqual(first.groupby('case_id').size().to_dict(), {'a': 51, 'b': 51})
        a['donor_status'], b['donor_status'] = 'ND', 'T1D'
        second, _ = balanced_sample({'b': b.sample(frac=1, random_state=4),
                                     'a': a.sample(frac=1, random_state=5)}, 60)
        self.assertEqual(first[['case_id', 'object_id']].values.tolist(), second[['case_id', 'object_id']].values.tolist())

    def test_all_required_markers_affect_transformed_features(self):
        sample = cells('a')
        for method in ['log1p', 'arcsinh']:
            x, specification = fit_transform(sample[MARKERS].values, MARKERS, method)
            self.assertEqual(specification['markers'], MARKERS)
            for marker in REQUIRED_REFINEMENT_MARKERS:
                altered = sample[MARKERS].values.copy()
                altered[0, MARKERS.index(marker)] *= 2
                y = transform(altered, specification)
                self.assertNotEqual(x[0, MARKERS.index(marker)], y[0, MARKERS.index(marker)])
                np.testing.assert_allclose(x[1:], y[1:])

    def test_scaling_uses_training_only_and_sparse_variation_not_dropped(self):
        raw = np.array([[0, 2, 1], [0, 2, 2], [0, 2, 3], [0, 2, 4], [20, 2, 5]], float)
        x, spec = fit_transform(raw, ['sparse', 'constant', 'varying'], 'log1p')
        self.assertEqual(spec['dropped_constant'], ['constant'])
        self.assertEqual(spec['sd_fallback_markers'], ['sparse'])
        before = json.dumps(spec, sort_keys=True)
        projected = transform(np.array([[100000, 2, 100000]], float), spec)
        self.assertGreater(projected.max(), x.max())  # no tail clipping
        self.assertEqual(before, json.dumps(spec, sort_keys=True))

    def test_unreviewed_clusters_cannot_produce_validated_cells(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        mm = {(f, 0): annotation(decision='pending') for f in pp}
        result = consensus_assignments(keys, pp, mm)
        self.assertEqual(result.classification_status.iloc[0], 'unreviewed')
        self.assertTrue(pd.isna(result.is_endocrine.iloc[0]))
        self.assertTrue(pd.isna(result.base_type.iloc[0]))

    def test_three_supported_votes_and_separate_state_consensus(self):
        keys = pd.DataFrame({'case_id': ['a', 'a'], 'object_id': ['1', '2']})
        pp = {f: projections([0, 0], [True, f in ['a', 'b']]) for f in ['a', 'b', 'c', 'd']}
        mm = {(f, 0): annotation(state='CHGA high' if f in ['a', 'b'] else '') for f in pp}
        result = consensus_assignments(keys, pp, mm)
        self.assertEqual(result.classification_status.tolist(), ['validated', 'lineage_candidate'])
        self.assertEqual(result.phenotype_state.iloc[0], '')
        self.assertTrue(result.is_endocrine.iloc[0])
        self.assertTrue(pd.isna(result.is_endocrine.iloc[1]))

    def test_overlapping_states_have_independent_three_vote_consensus(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a','b','c','d']}
        mm = {(f,0): annotation() for f in pp}
        for f in pp:
            mm[(f,0)]['phenotype_states'] = ['CD44 high','LGALS3 high'] + (['VIM high'] if f in ['a','b'] else [])
        result = consensus_assignments(keys,pp,mm)
        self.assertEqual(result.phenotype_states.iloc[0], ['Beta cell: CD44 high','Beta cell: LGALS3 high'])
        self.assertNotIn('VIM',result.phenotype_state.iloc[0])

    def test_competing_lineage_refinements_remain_ambiguous_in_either_order(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a','b','c','d']}
        mapping = {}
        for f in pp:
            row=annotation('');row['annotation_level']='lineage';row['lineage']='Endocrine' if f in ['a','b'] else 'Neural';mapping[(f,0)]=row
        parent=consensus_assignments(keys,pp,mapping)
        self.assertEqual(parent.classification_status.iloc[0],'lineage_candidate')
        self.assertTrue(pd.isna(parent.is_endocrine.iloc[0]))
        endocrine=consensus_assignments(keys,pp,{(f,0):annotation() for f in pp})
        neural_row=annotation('Neural');neural_row['lineage']='Neural'
        neural=consensus_assignments(keys,pp,{(f,0):neural_row for f in pp})
        for first,second,ln1,ln2 in [(endocrine,neural,'Endocrine','Neural'),(neural,endocrine,'Neural','Endocrine')]:
            result=merge_refinement(merge_refinement(parent,first,ln1),second,ln2)
            self.assertEqual(result.classification_status.iloc[0],'ambiguous_lineage_refinement')
            self.assertTrue(pd.isna(result.base_type.iloc[0]))
            self.assertTrue(pd.isna(result.is_endocrine.iloc[0]))

    def test_chga_unresolved_membership_requires_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'mapping.json'
            evidence = Path(tmp) / 'review.txt'
            evidence.write_text('Reviewed CHGA-positive hormone-low image references.')
            row = annotation('Hormone-unresolved endocrine cell')
            row.update(fit_id='log1p_k15', cluster_id=0,
                       evidence=[{'path': 'review.txt', 'sha256': sha256(evidence)}])
            spec = {'identity': 'test', 'config': {'transformations': ['log1p', 'arcsinh'], 'neighbor_counts': [15, 30]}}
            write_json(path, {'model_identity': 'test', 'clusters': [row]})
            with self.assertRaisesRegex(ValueError, 'endocrine evidence'):
                read_annotations(path, spec)
            row['endocrine_evidence'] = ['CHGA']
            write_json(path, {'model_identity': 'test', 'clusters': [row]})
            mapping = read_annotations(path, spec)
            self.assertIn('Hormone-unresolved endocrine cell', catalog_from_annotations(mapping)['base_type_names'])
            evidence.write_text('Changed evidence')
            with self.assertRaisesRegex(ValueError, 'missing or changed'):
                read_annotations(path, spec)

    def test_unknown_positivity_stays_nullable_not_false(self):
        result = attach_positivity(pd.DataFrame({'case_id': ['a'], 'object_id': ['1']}))
        self.assertTrue(result['positive__CD45'].isna().all())
        self.assertEqual(str(result['positive__CD45'].dtype), 'boolean')

    def test_failed_refinement_preserves_validated_parent(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        parent = consensus_assignments(keys, pp, {(f, 0): annotation() for f in pp})
        child = consensus_assignments(keys, pp, {})
        updated = merge_refinement(parent, child, 'Endocrine')
        self.assertEqual(updated.base_type.iloc[0], 'Beta cell')
        self.assertTrue(updated.is_endocrine.iloc[0])
        self.assertEqual(updated.refinement_status.iloc[0], 'unreviewed')

    def test_catalog_states_match_cell_state_and_are_parent_specific(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        mapping = {(f, 0): annotation(state='NaKATPase high') for f in pp}
        typed = consensus_assignments(keys, pp, mapping)
        catalog = catalog_from_annotations(mapping)
        self.assertEqual(typed.phenotype_state.iloc[0], 'Beta cell: NaKATPase high')
        self.assertEqual(catalog['states'][0]['name'], typed.phenotype_state.iloc[0])
        self.assertEqual(catalog['states'][0]['parent'], 'Beta cell')

    def test_supplemental_cannot_override_supported_identity(self):
        keys = pd.DataFrame({'case_id': ['a'], 'object_id': ['1']})
        pp = {f: projections([0]) for f in ['a', 'b', 'c', 'd']}
        mapping = {(f, 0): annotation() for f in pp}
        supported = consensus_assignments(keys, pp, mapping)
        with self.assertRaisesRegex(ValueError, 'only update unsupported'):
            merge_supplemental(supported, supported)
        unsupported = supported.copy()
        unsupported['projection_supported'] = False
        unsupported['classification_status'] = 'unsupported'
        resolved = merge_supplemental(unsupported, supported)
        self.assertTrue(resolved.supplemental_resolved.iloc[0])
        self.assertEqual(resolved.classification_status.iloc[0], 'validated')
        with self.assertRaisesRegex(ValueError, 'capped at 5000'):
            TypingConfig(supplemental=True).validate()

    def test_validation_reconstructs_every_donor_and_nine_marker_ablations(self):
        donors = [str(i) for i in range(21)]
        variants = validation_variants(donors, TypingConfig())
        lodo = [v for v in variants if v[0].startswith('lodo_')]
        self.assertEqual(len(lodo), 21)
        self.assertTrue(all(v[2] == v[3] and len(v[2]) == 1 for v in lodo))
        self.assertEqual({v[1].omitted_markers[0] for v in variants if v[0].startswith('ablate_')},
                         set(REQUIRED_REFINEMENT_MARKERS))

    def test_stability_metrics_report_support_loss_separately(self):
        keys = pd.DataFrame({'case_id': ['a', 'a'], 'object_id': ['1', '2']})
        pp = {f: projections([0, 0]) for f in ['a', 'b', 'c', 'd']}
        mapping = {(f, 0): annotation() for f in pp}
        baseline = consensus_assignments(keys, pp, mapping)
        reconstructed = baseline.copy()
        reconstructed['classification_status'] = ['diagnostic_label_transfer', 'unsupported']
        reconstructed['projection_supported'] = [True, False]
        report = stability_metrics(baseline, reconstructed)
        self.assertEqual(report['n_lost_support'], 1)
        self.assertEqual(report['base_label_agreement'], 1.)
        self.assertEqual(report['supported_coverage'], .5)

    def test_complete_small_discovery_apply_and_checkpoint_integrity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, output = root / 'input', root / 'output'
            (source / 'corrected').mkdir(parents=True)
            pd.DataFrame({'case_id': ['a', 'b'], 'donor_status': ['ND', 'T1D']}).to_csv(source / 'cohort.csv', index=False)
            for i, cid in enumerate(['a', 'b']):
                f = cells(cid, 40, i + 2)
                f.loc[0, MARKERS] = np.nan
                f.to_parquet(source / 'corrected' / f'{cid}-cells.parquet', index=False)
            configure_runtime(root / 'cache', 1)
            config = TypingConfig(cells_per_donor=25, neighbor_backend='exact')
            spec = discover(source, output, config)
            mapping_file = output / 'typing/model/cluster_annotations.template.json'
            # Pending template application must run without fabricating cell types.
            summary = apply(source, output, mapping_file)
            self.assertEqual(sum(v['n_cells'] for v in summary if v['classification_status'] == 'unmeasured'), 2)
            self.assertEqual(sum(v['n_cells'] for v in summary if v['classification_status'] == 'validated'), 0)
            f = pd.read_parquet(output / 'typing/cells/a.parquet')
            self.assertEqual(len(f), 40)
            self.assertTrue(f.is_endocrine.isna().all())
            resumed = discover(source, output, config)
            self.assertEqual(spec['identity'], resumed['identity'])
            # Disease-label edits must not change model identity or force a refit.
            pd.DataFrame({'case_id': ['a', 'b'], 'donor_status': ['T1D', 'ND']}).to_csv(source / 'cohort.csv', index=False)
            changed = discover(source, output, config)
            self.assertEqual(spec['identity'], changed['identity'])
            with self.assertRaisesRegex(ValueError, 'different inputs/configuration'):
                discover(source, output, TypingConfig(cells_per_donor=24, neighbor_backend='exact'))


if __name__ == '__main__':
    unittest.main()
