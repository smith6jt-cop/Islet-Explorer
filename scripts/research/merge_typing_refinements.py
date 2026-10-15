#!/usr/bin/env python3
"""Merge independently reviewed lineage refinements of one immutable parent.

Every child application is bound to the same parent files. Conflicting accepted
lineage identities become unavailable, irrespective of child input order.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
from common import MARKERS, sha256, write_json
from cell_typing import KEYS, atomic_parquet, canonical_hash, merge_refinement, parent_lineage_candidates, attach_positivity


def merge(parent_root, refinement_roots, output_root, positivity_root=None):
    parent_root, output_root = Path(parent_root), Path(output_root)
    parent_application = parent_root / 'typing/application.json'
    parent_files = sorted((parent_root / 'typing/cells').glob('*.parquet'))
    if not parent_files:
        raise ValueError('Explicit parent contains no typed cells')
    parent_hashes = {p.stem: sha256(p) for p in parent_files}
    positivity_hashes = ({p.stem: sha256(Path(positivity_root) / p.name) for p in parent_files}
                        if positivity_root is not None else None)
    children = []
    for root in map(Path, refinement_roots):
        application_file = root / 'typing/application.json'
        application = json.loads(application_file.read_text())
        if application.get('status') != 'complete' or application.get('parent_typing_hashes') != parent_hashes:
            raise ValueError('Refinement is incomplete or belongs to another parent')
        spec = json.loads((Path(application['model_root']) / 'specification.json').read_text())
        lineage = spec['config'].get('refinement')
        if not lineage:
            raise ValueError('Only explicit lineage refinements may be merged')
        if any(item['lineage'] == lineage for item in children):
            raise ValueError('Only one final reviewed refinement per lineage is permitted')
        children.append({'root': root, 'lineage': lineage, 'application': application,
                         'application_sha256': sha256(application_file), 'model_identity': spec['identity']})
    children.sort(key=lambda row: row['lineage'])
    identity = canonical_hash({'parent': parent_hashes, 'children': [
        {k: row[k] for k in ['lineage', 'application_sha256', 'model_identity']} for row in children],
        'policy': 'conflicting_supported_lineages_unavailable_v1', 'positivity_hashes': positivity_hashes})
    types, states = {}, {}
    catalog = json.loads((parent_root / 'typing/catalog.json').read_text())
    for root in [parent_root] + [row['root'] for row in children]:
        current = json.loads((root / 'typing/catalog.json').read_text())
        for record in current.get('base_types', []):
            name = record['name']
            if name in types and types[name] != record:
                raise ValueError(f'Conflicting base catalog records: {name}')
            types[name] = record
        for record in current.get('states', []):
            name = record['name']
            if name in states and states[name] != record:
                raise ValueError(f'Conflicting state catalog records: {name}')
            states[name] = record
    summary = []
    for path in parent_files:
        cid = path.stem
        parent = pd.read_parquet(path)
        original = parent.copy()
        identities = {tuple(row): [] for row in parent[KEYS].itertuples(index=False, name=None)}
        routes = {tuple(row): [] for row in parent[KEYS].itertuples(index=False, name=None)}
        for child in children:
            frame = pd.read_parquet(child['root'] / 'typing/cells' / path.name)
            if not pd.MultiIndex.from_frame(frame[KEYS]).equals(pd.MultiIndex.from_frame(original[KEYS])):
                raise ValueError('Child cell keys or order differ from immutable parent')
            eligible = parent_lineage_candidates(original, child['lineage'])
            selected = frame.loc[eligible].copy()
            # Nonselected child rows are copies of the original parent, never votes.
            parent = merge_refinement(parent, selected, child['lineage'])
            accepted = selected.classification_status.eq('validated') & selected.projection_supported.fillna(False)
            for row in selected.loc[accepted, KEYS].itertuples(index=False, name=None):
                identities[tuple(row)].append(child['model_identity'])
                routes[tuple(row)].append(child['lineage'])
        parent['refinement_model_identities'] = [identities[tuple(row)] for row in parent[KEYS].itertuples(index=False, name=None)]
        parent['refinement_parent_pools'] = [routes[tuple(row)] for row in parent[KEYS].itertuples(index=False, name=None)]
        parent['typing_run_identity'] = identity
        if positivity_root is not None:
            parent = parent.drop(columns=[f'positive__{marker}' for marker in MARKERS], errors='ignore')
            parent = attach_positivity(parent, Path(positivity_root) / path.name)
            unavailable = parent.classification_status.isin(['segmentation_quality_excluded', 'unmeasured'])
            parent.loc[unavailable, [f'positive__{marker}' for marker in MARKERS]] = pd.NA
        atomic_parquet(parent, output_root / 'typing/cells' / path.name)
        for status, n in parent.classification_status.value_counts().items():
            summary.append({'case_id': cid, 'classification_status': status, 'n_cells': int(n),
                'n_endocrine': int(parent.loc[parent.classification_status.eq(status), 'is_endocrine'].fillna(False).sum())})
    catalog['base_types'] = [types[k] for k in sorted(types)]
    catalog['base_type_names'] = sorted(types)
    catalog['states'] = [states[k] for k in sorted(states)]
    write_json(output_root / 'typing/catalog.json', catalog)
    result = {'schema_version': 2, 'identity': identity, 'status': 'complete', 'method': 'reviewed_hierarchical_ensemble',
        'parent_root': str(parent_root.resolve()), 'parent_application_sha256': sha256(parent_application),
        'parent_typing_hashes': parent_hashes, 'positivity_hashes': positivity_hashes, 'refinements': [
            {'lineage': row['lineage'], 'root': str(row['root'].resolve()), 'model_identity': row['model_identity'],
             'application_sha256': row['application_sha256']} for row in children],
        'summary': summary, 'validated_cells': sum(row['n_cells'] for row in summary if row['classification_status'] == 'validated'),
        'conflict_policy': 'Conflicting accepted lineage refinements have no final identity or endocrine denominator membership.'}
    write_json(output_root / 'typing/application.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-root', type=Path, required=True)
    parser.add_argument('--refinement-root', type=Path, action='append', required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--positivity-root', type=Path)
    args = parser.parse_args()
    print(json.dumps(merge(args.parent_root, args.refinement_root, args.output_root, args.positivity_root), indent=2))


if __name__ == '__main__':
    main()
