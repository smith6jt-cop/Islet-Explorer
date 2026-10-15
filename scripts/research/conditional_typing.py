"""Shared fold-local reconstruction of the actually reviewed typing hierarchy."""
from dataclasses import replace
import json
from pathlib import Path
from common import sha256
from cell_typing import TypingConfig,discover
from cell_typing_validation import reconstruct_hierarchy


def hierarchy_reference(reference,excluded):
    reference=Path(reference);application_path=reference/'typing/application.json'
    application=json.loads(application_path.read_text())
    if application.get('method')!='reviewed_hierarchical_ensemble':
        raise ValueError('Supporting v4 inference requires the reviewed broad-plus-refinement hierarchy, not a broad-only surrogate')
    broad=Path(application['parent_root'])
    broad_application_path=broad/'typing/application.json'
    broad_application=json.loads(broad_application_path.read_text())
    broad_specification_path=Path(broad_application['model_root'])/'specification.json'
    specification=json.loads(broad_specification_path.read_text())
    training=[cid for cid in specification['all_cohort_donors'] if cid not in excluded]
    fingerprint=dict(application_sha256=sha256(application_path),
        broad_model_sha256=sha256(broad_specification_path),broad_application_sha256=sha256(broad_application_path),
        broad_reference={cid:sha256(broad/'typing/cells'/f'{cid}.parquet') for cid in training},
        final_reference={cid:sha256(reference/'typing/cells'/f'{cid}.parquet') for cid in training},
        hierarchy_code_sha256=sha256(Path(__file__).with_name('cell_typing_validation.py')),
        conditional_code_sha256=sha256(__file__))
    return application,broad,specification,fingerprint


def reconstruct(root,directory,reference,excluded,omitted,project_donors,excluded_types=(),workers=1,fit_workers=1,endocrine_only=False):
    application,broad,source,_=hierarchy_reference(reference,excluded)
    config=replace(TypingConfig(**source['config']),omitted_markers=tuple(sorted(omitted)),refinement=None,supplemental=False)
    specification=discover(root,directory,config,workers,fit_workers,exclude_donors=excluded)
    if set(excluded)&set(specification['training_inputs']):raise ValueError('Held-out donor entered conditional typing fit')
    if set(omitted)&set(specification['features']):raise ValueError('Omitted outcome marker entered conditional typing fit')
    catalog=json.loads((reference/'typing/catalog.json').read_text())
    allowed={entry['lineage'] for entry in catalog['base_types'] if entry['name'] not in excluded_types}
    lineages=[entry['lineage'] for entry in application['refinements'] if entry['lineage'] in allowed]
    if endocrine_only:lineages=[lineage for lineage in lineages if lineage.lower()=='endocrine']
    if not lineages:raise ValueError('No independent supported typing lineage remains for this outcome module')
    reconstruct_hierarchy(root,directory,reference,broad,specification,config,excluded,project_donors,lineages,workers,fit_workers)
    return specification,config,directory/'diagnostic_predictions'
