#!/usr/bin/env python3
"""Reconstruct cell typing for donor, seed, resolution, graph and marker sensitivity.

A fresh cluster-to-label alignment uses only reference cells in training donors.
The resulting labels are *stability diagnostics*, never newly reviewed production
identities. Comparing them to the full-cohort classifier measures reproducibility,
not independent biological accuracy. Image-reference validation is separate.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
from multiprocessing import get_context
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from common import MARKERS, sha256, write_json
from cell_typing import (TypingConfig, REQUIRED_REFINEMENT_MARKERS, KEYS, atomic_parquet,
                        discover, read_cells, project_frame, consensus_assignments,
                        configure_runtime, safe_name, canonical_hash, attach_quality_status,
                        parent_lineage_candidates, restrict_to_parent, merge_refinement,
                        catalog_from_annotations)


def validation_variants(donors, config):
    variants = []
    for cid in sorted(donors):
        variants.append((f'lodo_{safe_name(cid)}', config, [str(cid)], [str(cid)]))
    for seed in [config.seed + 1, config.seed + 2]:
        variants.append((f'seed_{seed}', replace(config, seed=seed), [], list(donors)))
    for resolution in [.5, 1.5]:
        variants.append((f'resolution_{resolution}', replace(config, resolution=resolution), [], list(donors)))
    for counts in [(10, 30), (15, 45)]:
        variants.append((f'neighbors_{counts[0]}_{counts[1]}', replace(config, neighbor_counts=counts), [], list(donors)))
    for marker in REQUIRED_REFINEMENT_MARKERS:
        variants.append((f'ablate_{safe_name(marker)}', replace(config, omitted_markers=(marker,)), [], list(donors)))
    for minimum in [1.0, 10.0]:
        variants.append((f'cell_area_minimum_{minimum:g}', replace(config, minimum_cell_area_um2=minimum), [], list(donors)))
    return variants


def align_training_clusters(model_root, reference_root, specification, minimum_reference_cells=10,
                            minimum_purity=.80, minimum_reference_coverage=.50):
    """Align only from training-donor validated reference identities; save purity."""
    reference = []
    for cid in specification['training_inputs']:
        frame = pd.read_parquet(Path(reference_root) / 'typing' / 'cells' / f'{cid}.parquet')
        for key in KEYS:
            frame[key] = frame[key].astype(str)
        frame = frame.loc[frame.classification_status.eq('validated') & frame.projection_supported.fillna(False)]
        if specification['config'].get('refinement'):
            frame = frame.loc[parent_lineage_candidates(frame, specification['config']['refinement'])]
        if 'phenotype_states' not in frame:
            frame['phenotype_states'] = frame.phenotype_state.map(lambda s: [s] if isinstance(s, str) and s else [])
        reference.append(frame[KEYS + ['base_type', 'lineage', 'phenotype_state', 'phenotype_states']])
    reference = pd.concat(reference, ignore_index=True)
    mapping, audit = {}, []
    for method in specification['config']['transformations']:
        for k in specification['config']['neighbor_counts']:
            fit_id = f'{method}_k{k}'
            memberships = pd.read_parquet(Path(model_root) / fit_id / 'memberships.parquet')
            joined = memberships.merge(reference, on=KEYS, how='left', validate='one_to_one')
            for cluster, group in joined.groupby('cluster_id'):
                named = group.loc[group.base_type.notna()]
                counts = named.groupby(['base_type', 'lineage']).size().sort_values(ascending=False)
                winner = None if counts.empty else counts.index[0]
                n = 0 if counts.empty else int(counts.iloc[0])
                purity = n / len(named) if len(named) else 0.
                coverage = len(named) / len(group)
                supported = n >= minimum_reference_cells and purity >= minimum_purity and coverage >= minimum_reference_coverage
                states = []
                if supported:
                    selected = named.loc[named.base_type.eq(winner[0]) & named.lineage.eq(winner[1])]
                    state_counts = selected.phenotype_states.explode().dropna().value_counts()
                    states = [str(state).removeprefix(str(winner[0]) + ': ') for state, count in state_counts.items()
                              if state and count / len(selected) >= minimum_purity]
                    mapping[(fit_id, int(cluster))] = {
                        'decision': 'accepted', 'base_type': winner[0], 'lineage': winner[1],
                        'phenotype_states': sorted(states), 'alignment_only': True}
                audit.append({'fit_id': fit_id, 'cluster_id': int(cluster), 'n_training_members': len(group),
                    'n_named_reference': len(named), 'dominant_reference_type': None if winner is None else winner[0],
                    'purity': purity, 'reference_coverage': coverage, 'alignment_supported': supported})
    return mapping, pd.DataFrame(audit)


def align_training_lineages(model_root, broad_reference_root, specification,
                            minimum_reference_cells=10, minimum_purity=.80,
                            minimum_reference_coverage=.50):
    """Align broad candidate pools using training donors only, never final bases.

    Ambiguous reference candidate pools do not nominate a unique training lineage.
    The result is a diagnostic transfer of reviewed broad nominations, not a new
    biological validation or an endocrine-cell denominator.
    """
    reference = []
    for cid in specification['training_inputs']:
        frame = pd.read_parquet(Path(broad_reference_root) / 'typing/cells' / f'{cid}.parquet')
        candidates = frame.candidate_lineages.map(lambda s: json.loads(s) if isinstance(s, str) and s else [])
        frame = frame[KEYS].copy()
        frame['candidate_lineage'] = candidates.map(lambda values: values[0] if len(values) == 1 else None)
        reference.append(frame)
    reference = pd.concat(reference, ignore_index=True)
    mapping, audit = {}, []
    for method in specification['config']['transformations']:
        for k in specification['config']['neighbor_counts']:
            fit_id = f'{method}_k{k}'
            memberships = pd.read_parquet(Path(model_root) / fit_id / 'memberships.parquet')
            joined = memberships.merge(reference, on=KEYS, how='left', validate='one_to_one')
            for cluster, group in joined.groupby('cluster_id'):
                named = group.candidate_lineage.dropna()
                counts = named.value_counts()
                winner = None if counts.empty else str(counts.index[0])
                n = 0 if counts.empty else int(counts.iloc[0])
                purity = n / len(named) if len(named) else 0.
                coverage = len(named) / len(group)
                accepted = n >= minimum_reference_cells and purity >= minimum_purity and coverage >= minimum_reference_coverage
                if accepted:
                    mapping[(fit_id, int(cluster))] = {'decision': 'accepted', 'annotation_level': 'lineage',
                        'base_type': '', 'lineage': winner, 'alignment_only': True}
                audit.append({'fit_id': fit_id, 'cluster_id': int(cluster), 'n_training_members': len(group),
                    'n_nominated_reference': len(named), 'dominant_lineage': winner, 'purity': purity,
                    'reference_coverage': coverage, 'alignment_supported': accepted})
    return mapping, pd.DataFrame(audit)


_PROJECTION_CACHE = {}


def _project_broad_donor(task):
    input_root, root, model_root, spec, mapping, cid = task
    configure_runtime(root / 'projection_cache', 1)
    cells = read_cells(input_root / 'corrected' / f'{cid}-cells.parquet', cid)
    pp, _ = project_frame(cells, model_root, spec, spec['config']['batch_size'], _PROJECTION_CACHE)
    predicted = attach_quality_status(consensus_assignments(cells[KEYS], pp, mapping), cells,
        spec['config']['minimum_cell_area_um2'], spec['config'].get('technical_exclusion_regions', ()))
    predicted['independent_biological_validation'] = False
    atomic_parquet(predicted, root / 'diagnostic_broad_parent/typing/cells' / f'{cid}.parquet')
    diagnostic = cells[KEYS].copy()
    for fit_id, frame in pp.items():
        for column in frame:
            diagnostic[f'{fit_id}__{column}'] = frame[column].to_numpy()
    atomic_parquet(diagnostic, root / 'broad_projection_diagnostics' / f'{cid}.parquet')
    return {'case_id': cid, 'n_cells': len(predicted)}


def prepare_hierarchy_models(input_root, root, broad_reference_root, variant_spec,
                             variant, excluded, lineages, workers=1, fit_workers=1):
    """Precompute fold-local nominations and unsupervised refiners, without final labels."""
    input_root, root, broad_reference_root = map(Path, [input_root, root, broad_reference_root])
    model_root = root / 'typing/model'
    mapping, alignment = align_training_lineages(model_root, broad_reference_root, variant_spec)
    atomic_parquet(alignment, root / 'training_broad_lineage_alignment.parquet')
    parent_root = root / 'diagnostic_broad_parent'
    parent_identity = canonical_hash({'model': variant_spec['identity'],
        'broad_reference_training_hashes': {cid: sha256(broad_reference_root / 'typing/cells' / f'{cid}.parquet')
                                           for cid in variant_spec['training_inputs']},
        'policy': 'training_only_broad_candidate_transfer_v1'})
    checkpoint = root / 'broad_candidate_projection_complete.json'
    if checkpoint.exists():
        previous = json.loads(checkpoint.read_text())
        if previous['identity'] != parent_identity:
            raise ValueError('Fold-local candidate parent belongs to another broad reference')
        if any(sha256(parent_root / 'typing/cells' / f'{cid}.parquet') != digest for cid, digest in previous['parent_hashes'].items()):
            raise ValueError('Fold-local candidate checkpoint changed')
    else:
        tasks = [(input_root, root, model_root, variant_spec, mapping, cid) for cid in variant_spec['all_cohort_donors']]
        if workers > 1:
            with ProcessPoolExecutor(max_workers=min(workers, 8), mp_context=get_context('spawn')) as pool:
                list(pool.map(_project_broad_donor, tasks))
        else:
            for task in tasks:
                _project_broad_donor(task)
            _PROJECTION_CACHE.clear()
        write_json(checkpoint, {'identity': parent_identity, 'parent_hashes': {
            cid: sha256(parent_root / 'typing/cells' / f'{cid}.parquet') for cid in variant_spec['all_cohort_donors']}})
    write_json(parent_root / 'typing/catalog.json', catalog_from_annotations({}))
    write_json(parent_root / 'typing/application.json', {'status': 'diagnostic_only',
        'model_identity': variant_spec['identity'], 'biological_validation': False})
    refinement_records = []
    for lineage in sorted(lineages):
        refinement_root = root / 'hierarchy' / safe_name(lineage)
        refinement_model = refinement_root / 'typing/model'
        config = replace(variant, refinement=lineage)
        spec = discover(input_root, refinement_root, config, workers, fit_workers,
            model_root=refinement_model, parent_typing_root=parent_root, exclude_donors=excluded)
        refinement_records.append({'lineage': lineage, 'root': str(refinement_root.resolve()),
            'model_identity': spec['identity']})
    write_json(root / 'hierarchy_discovery_complete.json', {'parent_identity': parent_identity,
        'refinements': refinement_records, 'final_biological_labels_used': False,
        'excluded_donors': list(excluded), 'omitted_markers': list(variant.omitted_markers)})
    return parent_root, refinement_records


def reconstruct_hierarchy(input_root, root, reference_root, broad_reference_root,
                          variant_spec, variant, excluded, project_donors, lineages,
                          workers=1, fit_workers=1):
    """Reconstruct preprocessing, broad pools and all refiners within a fold."""
    parent_root, _ = prepare_hierarchy_models(input_root, root, broad_reference_root,
        variant_spec, variant, excluded, lineages, workers, fit_workers)
    model_root = root / 'typing/model'
    partitions = []
    for cid in project_donors:
        diagnostic = pd.read_parquet(root / 'broad_projection_diagnostics' / f'{cid}.parquet')
        pp = {}
        for method in variant.transformations:
            for k in variant.neighbor_counts:
                fit_id = f'{method}_k{k}'
                columns = [column for column in diagnostic if column.startswith(fit_id + '__')]
                pp[fit_id] = diagnostic[KEYS + columns].rename(columns={column: column.split('__', 1)[1] for column in columns})
        partitions.extend(partition_metrics(reference_root, model_root, pp, cid))
    predictions = {cid: pd.read_parquet(parent_root / 'typing/cells' / f'{cid}.parquet') for cid in project_donors}
    refinement_records = []
    for lineage in sorted(lineages):
        refinement_root = root / 'hierarchy' / safe_name(lineage)
        refinement_model = refinement_root / 'typing/model'
        spec = json.loads((refinement_model / 'specification.json').read_text())
        labels, audit = align_training_clusters(refinement_model, reference_root, spec)
        atomic_parquet(audit, refinement_root / 'training_cluster_alignment.parquet')
        model_cache = {}
        for cid in project_donors:
            cells = restrict_to_parent(read_cells(input_root / 'corrected' / f'{cid}-cells.parquet', cid), parent_root, lineage, False)
            pp, _ = project_frame(cells, refinement_model, spec, variant.batch_size, model_cache)
            predicted = attach_quality_status(consensus_assignments(cells[KEYS], pp, labels), cells,
                variant.minimum_cell_area_um2, variant.technical_exclusion_regions)
            predicted['independent_biological_validation'] = False
            atomic_parquet(predicted, refinement_root / 'diagnostic_predictions' / f'{cid}.parquet')
            predictions[cid] = merge_refinement(predictions[cid], predicted, lineage)
        del model_cache
        refinement_records.append({'lineage': lineage, 'model_identity': spec['identity'],
            'training_alignment_sha256': sha256(refinement_root / 'training_cluster_alignment.parquet')})
    metrics = []
    for cid, predicted in predictions.items():
        predicted.loc[predicted.classification_status.eq('validated'), 'classification_status'] = 'diagnostic_label_transfer'
        predicted['independent_biological_validation'] = False
        atomic_parquet(predicted, root / 'diagnostic_predictions' / f'{cid}.parquet')
        reference = pd.read_parquet(reference_root / 'typing/cells' / f'{cid}.parquet')
        metrics.append({'case_id': cid, **stability_metrics(reference, predicted)})
    return metrics, partitions, refinement_records


def stability_metrics(reference, reconstructed):
    """Count changed and unavailable predictions explicitly; do not hide attrition."""
    merged = reference.merge(reconstructed, on=KEYS, suffixes=('_reference', '_new'), validate='one_to_one')
    reference_valid = merged.classification_status_reference.eq('validated') & merged.projection_supported_reference.fillna(False)
    new_valid = merged.classification_status_new.eq('diagnostic_label_transfer') & merged.projection_supported_new.fillna(False)
    shared = reference_valid & new_valid
    agreement = merged.loc[shared, 'base_type_reference'].eq(merged.loc[shared, 'base_type_new'])
    state_reference = merged.loc[shared, 'phenotype_state_reference'].fillna('')
    state_new = merged.loc[shared, 'phenotype_state_new'].fillna('')
    return {'n_cells': len(merged), 'n_reference_validated': int(reference_valid.sum()),
            'n_supported_diagnostic': int(new_valid.sum()), 'n_shared_supported': int(shared.sum()),
            'n_changed_base_type': int((~agreement).sum()),
            'n_changed_state': int(state_reference.ne(state_new).sum()),
            'n_lost_support': int((reference_valid & ~new_valid).sum()),
            'n_gained_support': int((~reference_valid & new_valid).sum()),
            'base_label_agreement': float(agreement.mean()) if len(agreement) else None,
            'supported_coverage': float(new_valid.mean()) if len(merged) else None,
            'mean_confidence_change': float((merged.loc[shared, 'confidence_new'] - merged.loc[shared, 'confidence_reference']).mean()) if shared.any() else None}


def partition_metrics(reference_root, variant_model, projections, cid):
    path = Path(reference_root) / 'typing' / 'model' / 'projections' / f'{cid}.parquet'
    if not path.is_file():
        return []
    baseline = pd.read_parquet(path)
    rows = []
    for fit_id, projected in projections.items():
        column = fit_id + '__cluster_id'
        support = fit_id + '__projection_supported'
        if column not in baseline:
            continue
        # project_frame preserves the canonical corrected-cell key ordering.
        left = baseline[KEYS + [column, support]]
        right = projected[KEYS + ['cluster_id', 'projection_supported']]
        both = left.merge(right, on=KEYS, validate='one_to_one')
        valid = both[support] & both.projection_supported & both[column].ge(0) & both.cluster_id.ge(0)
        if valid.sum() < 3:
            continue
        rows.append({'case_id': cid, 'fit_id': fit_id, 'n_shared_supported': int(valid.sum()),
            'adjusted_rand': float(adjusted_rand_score(both.loc[valid, column], both.loc[valid, 'cluster_id'])),
            'normalized_mutual_information': float(normalized_mutual_info_score(both.loc[valid, column], both.loc[valid, 'cluster_id']))})
    return rows


def run_validation(input_root, output_root, reference_root, config=None, workers=1, fit_workers=1, only=()):
    input_root, output_root, reference_root = map(Path, [input_root, output_root, reference_root])
    if output_root.resolve() == reference_root.resolve():
        raise ValueError('Validation outputs must be separate from production typing')
    spec = json.loads((reference_root / 'typing' / 'model' / 'specification.json').read_text())
    config = config or TypingConfig(**spec['config'])
    if config.refinement:
        raise ValueError('Broad typing validation must precede independently audited lineage validation')
    donors = pd.read_csv(input_root / 'cohort.csv', dtype=str).case_id.tolist()
    variants = validation_variants(donors, config)
    if only:
        missing = set(only) - {v[0] for v in variants}
        if missing:
            raise ValueError(f'Unknown validation variants: {sorted(missing)}')
        variants = [v for v in variants if v[0] in set(only)]
    reference_hashes = {cid: sha256(reference_root / 'typing' / 'cells' / f'{cid}.parquet') for cid in donors}
    application_path = reference_root / 'typing/application.json'
    application = json.loads(application_path.read_text()) if application_path.exists() else {}
    hierarchy = application.get('method') == 'reviewed_hierarchical_ensemble'
    broad_reference_root = Path(application['parent_root']) if hierarchy else None
    lineages = [row['lineage'] for row in application.get('refinements', [])]
    identity = canonical_hash({'reference': reference_hashes, 'config': spec['config'],
        'reference_application_sha256': sha256(application_path) if application_path.exists() else None,
        'reconstruction': 'complete_hierarchy_v1' if hierarchy else 'broad_only_v1'})
    output_root.mkdir(parents=True, exist_ok=True)
    complete_path = output_root / 'validation_manifest.json'
    if complete_path.exists() and json.loads(complete_path.read_text())['identity'] != identity:
        raise ValueError('Reference typing changed; start a new validation output root')
    completed, failed = [], []
    for name, variant, excluded, project_donors in variants:
        root = output_root / name
        checkpoint = root / 'validation_complete.json'
        if checkpoint.exists():
            record = json.loads(checkpoint.read_text())
            if record['reference_identity'] != identity:
                raise ValueError('Validation checkpoint belongs to a different reference')
            completed.append(record)
            continue
        try:
            variant_spec = discover(input_root, root, variant, workers, fit_workers, exclude_donors=excluded)
            model_root = root / 'typing' / 'model'
            if hierarchy:
                metrics, partitions, refinement_records = reconstruct_hierarchy(input_root, root,
                    reference_root, broad_reference_root, variant_spec, variant, excluded,
                    project_donors, lineages, workers, fit_workers)
                record = {'name': name, 'reference_identity': identity, 'model_identity': variant_spec['identity'],
                    'excluded_donors': excluded, 'omitted_markers': list(variant.omitted_markers),
                    'metrics': metrics, 'partition_stability': partitions, 'refinements': refinement_records,
                    'reconstruction_scope': 'fold_local_broad_candidates_and_all_lineage_refinements',
                    'interpretation': 'Full hierarchical reconstruction stability with training-donor-only label alignment; not independent biological accuracy.',
                    'training_label_alignment': {'minimum_cells': 10, 'minimum_purity': .8, 'minimum_reference_coverage': .5},
                    'status': 'complete'}
                write_json(checkpoint, record)
                completed.append(record)
                write_json(complete_path, {'identity': identity, 'requested': [v[0] for v in variants],
                    'completed': [v['name'] for v in completed], 'failed': failed, 'status': 'running'})
                continue
            mapping, alignment = align_training_clusters(model_root, reference_root, variant_spec)
            atomic_parquet(alignment, root / 'training_cluster_alignment.parquet')
            metrics, partitions, model_cache = [], [], {}
            for cid in project_donors:
                cells = read_cells(input_root / 'corrected' / f'{cid}-cells.parquet', cid)
                pp, measured = project_frame(cells, model_root, variant_spec, variant.batch_size, model_cache)
                predicted = consensus_assignments(cells[KEYS], pp, mapping)
                predicted = attach_quality_status(predicted, cells, variant.minimum_cell_area_um2, variant.technical_exclusion_regions)
                predicted.loc[predicted.classification_status.eq('validated'), 'classification_status'] = 'diagnostic_label_transfer'
                predicted['independent_biological_validation'] = False
                atomic_parquet(predicted, root / 'diagnostic_predictions' / f'{cid}.parquet')
                reference = pd.read_parquet(reference_root / 'typing' / 'cells' / f'{cid}.parquet')
                metrics.append({'case_id': cid, **stability_metrics(reference, predicted)})
                for fit_id in pp:
                    pp[fit_id] = pd.concat([cells[KEYS], pp[fit_id]], axis=1)
                partitions.extend(partition_metrics(reference_root, model_root, pp, cid))
            record = {'name': name, 'reference_identity': identity, 'model_identity': variant_spec['identity'],
                'excluded_donors': excluded, 'omitted_markers': list(variant.omitted_markers),
                'metrics': metrics, 'partition_stability': partitions,
                'interpretation': 'Reconstruction stability against reviewed full-cohort identities; not independent biological accuracy.',
                'training_label_alignment': {'minimum_cells': 10, 'minimum_purity': .8, 'minimum_reference_coverage': .5},
                'status': 'complete'}
            del model_cache
            write_json(checkpoint, record)
            completed.append(record)
        except Exception as error:
            record = {'name': name, 'status': 'failed', 'error_type': type(error).__name__, 'error': str(error)}
            write_json(root / 'validation_failure.json', record)
            failed.append(record)
        write_json(complete_path, {'identity': identity, 'requested': [v[0] for v in variants],
                   'completed': [v['name'] for v in completed], 'failed': failed, 'status': 'running'})
    rows = [{'variant': r['name'], **m} for r in completed for m in r['metrics']]
    atomic_parquet(pd.DataFrame(rows), output_root / 'cell_typing_stability.parquet')
    report = {'identity': identity, 'requested': [v[0] for v in variants],
              'completed': [v['name'] for v in completed], 'failed': failed,
              'status': 'complete' if not failed else 'incomplete',
              'reference_hashes': reference_hashes, 'production_labels_modified': False,
              'independent_biological_validation': False}
    write_json(complete_path, report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    p.add_argument('--reference-root', required=True, type=Path)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--fit-workers', type=int, default=1)
    p.add_argument('--only', action='append', default=[])
    args = p.parse_args()
    configure_runtime(args.output_root / 'cache', args.workers)
    result = run_validation(args.input_root, args.output_root, args.reference_root,
                            workers=args.workers, fit_workers=args.fit_workers, only=args.only)
    print(json.dumps(result, indent=2))
    if result['failed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
