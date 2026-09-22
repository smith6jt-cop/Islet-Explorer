#!/usr/bin/env python3
"""Verify and inventory the six-donor review delivery without accepting new donors."""
import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import shutil

import cv2
import numpy as np
import pandas as pd
import shapely
import tifffile

from prepare_endocrine_candidates import digest, write_json


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('rollout-root', 'project-root', 'gallery-root', 'service-record'):
        parser.add_argument('--' + arg, required=True, type=Path)
    args = parser.parse_args()
    root = args.rollout_root.resolve(strict=True)
    project = args.project_root.resolve(strict=True)
    gallery = args.gallery_root.resolve(strict=True)
    final = project / 'delivery_manifest.json'
    if final.exists():
        raise ValueError('Preserve delivered review and verification records')
    staging = read(project / 'staging_manifest.json')
    approval = read(project / 'boundary_approval.json')
    promotion = read(root / 'exports/78/promotion_manifest.json')
    regions = read(root / 'regions/78/verification.json')
    loaded = read(root / 'tests/import_gate/verification.json')
    rendered = read(gallery / 'review_manifest.json')
    if loaded['status'] != 'passed' or loaded['image_count'] != 6:
        raise ValueError('Fresh native six-image load verification required')
    if regions['status'] != 'approved_geometry_regions_verified' or rendered['approved_cases'] != ['78']:
        raise ValueError('Approved region/galllery verification incomplete')
    if digest(approval['source_review_data']) != approval['source_review_data_sha256']:
        raise ValueError('Owner-reviewed source changed')
    if digest(approval['source_review_project']) != approval['source_review_project_sha256']:
        raise ValueError('Owner-reviewed project changed')
    documents = read(project / 'project.qpproj')['images']
    rows, prediction_times, boundary_times, checks = [], {}, {}, []
    for donor in staging['donors']:
        cid = donor['case_id']
        native = Path(donor['review_data'])
        source = Path(donor['source_data'])
        if source.samefile(native) or digest(source) != donor['data_sha256']:
            raise ValueError('Source changed or native review is not independent: ' + cid)
        if cid == '78':
            evidence = promotion
            expected_hash = evidence['output_native_sha256']
            approved, pending = evidence['accepted_islets'], 0
            if evidence['unchanged_cell_content_sha256_before'] != evidence['unchanged_cell_content_sha256_after']:
                raise ValueError('Approved donor raw cell content changed')
        else:
            evidence = read(project / f'verification/{cid}-pixel_review_verification.json')
            if evidence['status'] != 'native_review_roundtrip_passed' or evidence['accepted_islets'] != 0:
                raise ValueError('Missing donor review verification or unsupported approval')
            if evidence['cell_record_sha256_before'] != evidence['cell_record_sha256_after']:
                raise ValueError('Unreviewed donor cell content changed')
            expected_hash = evidence['review_data_sha256']
            approved, pending = 0, evidence['new_boundary_candidates']
            candidate = root / 'candidates' / cid
            manifest = read(candidate / 'prediction_manifest.json')
            geometry = read(candidate / 'independent_geometry_check.json')
            if geometry['status'] != 'passed' or digest(candidate / 'review_candidates.geojson') != geometry['review_candidates_sha256']:
                raise ValueError('Independent geometry verification differs')
            if manifest['operating_point']['minimum_votes'] != 160 or manifest['operating_point']['forest_trees'] != 160:
                raise ValueError('Fixed operating point changed')
            if manifest['source_image'] != donor['source_image']:
                raise ValueError('Donor source image differs')
            if manifest['model_sha256'] != rendered['summaries'][0]['model_sha256']:
                raise ValueError('Classifier differs across donors')
            prediction_times[cid] = manifest['runtime_seconds']
            boundary_times[cid] = manifest['boundary_runtime_seconds']
            checks.append(geometry)
        if digest(native) != expected_hash or evidence['cells'] != donor['cell_count']:
            raise ValueError('Delivered native objects changed')
        entry = next(e for e in documents if int(e['entryID']) == int(donor['entry_id']))
        if cid == '78' and not entry['metadata']['Delivery status'].startswith('HPAP boundaries accepted'):
            raise ValueError('Stale HPAP acceptance metadata')
        if cid != '78' and not entry['metadata']['Delivery status'].startswith('Awaiting donor-specific visual review'):
            raise ValueError('Stale other-donor review metadata')
        summary = read(native.parent / 'summary.json')['hierarchy']['annotationClassificationCounts']
        if int(summary.get('Islet', 0)) != approved or int(summary.get('Endocrine pixel candidate', 0)) != pending:
            raise ValueError('Native class counts differ')
        rows.append(dict(case_id=cid, donor_status=donor['donor_status'], entry_id=int(donor['entry_id']),
            cell_count=donor['cell_count'], current_boundaries=approved+pending,
            accepted_boundaries=approved, unreviewed_candidates=pending,
            review_data=str(native), review_data_sha256=expected_hash, source_data_sha256=donor['data_sha256'],
            source_image=donor['source_image'], boundary_status='owner_approved' if approved else 'awaiting_visual_review',
            cell_membership_rebuilt=bool(approved), regional_measurements_rebuilt=bool(approved)))
    if len(rows) != 6 or len({r['source_image'] for r in rows}) != 6:
        raise ValueError('Expected six distinct acquisitions')
    shutil.copy2(project / 'donors.csv', project / 'source_donors.csv')
    pd.DataFrame(rows).to_csv(project / 'donors.csv', index=False)
    readme = f'''# Endocrine boundary review — 18 September 2026

Open **File → Project → Open project** in QuPath 0.7 and select:

```
{project / 'project.qpproj'}
```

Double-click a donor. Its existing cell/nucleus objects, raw measurements and
annotations load with the image; no import or segmentation rerun is needed.

- HPAP-078: 1,431 cyan `Islet` objects are the exact owner-approved shapes.
  Native core cell parents, 31 marker summaries, signed distances and 20-µm
  context regions are rebuilt. Grey objects preserve the previous segmentation.
- 6090, 6232, 6520, 6539 and 6558: 2,030 cyan `Endocrine pixel candidate`
  objects await individual image review. Their existing cell parents and regional
  measurements still refer to the grey previous segmentation.
- All 943,068 original cells remain. Boundary approval does not certify phenotype
  gates. Save review edits here; earlier projects remain unchanged.

The [offline gallery](http://127.0.0.1:3840/{gallery.name}/) contains whole-field
and size-stratified examples with PNG/PDF plots. [donors.csv](donors.csv) gives
current counts; `source_donors.csv` preserves the original staging inventory.
`boundary_approval.json` records the owner's decision and source hashes.

Rebuilt HPAP tables are in `{root / 'regions/78'}`. Core/peri measurements include
all 31 markers and 10/20/30-µm bands, with missing measurements and overlapping
footprints flagged. The measurement export preserves native coordinate precision.
The current app still uses its earlier candidate artifacts.

The earlier owner-reviewed project is `{Path(approval['source_review_project'])}`.
Reopen it to return to the prior review; source data restoration is unnecessary.
'''
    (project / 'README.md').write_text(readme)
    service = read(args.service_record)
    service['url'] = f'http://127.0.0.1:3840/{gallery.name}/'
    service['existing_service_reused'] = True
    write_json(root / 'review_http.json', service)
    storage = {label: dict(path=str(folder), bytes=sum(p.stat().st_size for p in folder.rglob('*') if p.is_file()),
                          files=sum(p.is_file() for p in folder.rglob('*')))
               for label, folder in [('native_project', project), ('analysis_and_diagnostics', root), ('gallery', gallery)]}
    for record in storage.values():
        record['gib'] = record['bytes'] / 1024**3
    versions = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                    shapely=shapely.__version__, opencv=cv2.__version__, tifffile=tifffile.__version__, qupath='0.7.0')
    source_dir = Path(__file__).resolve().parent
    script_names = ['stage_adaptive_rollout.py', 'promote_reviewed_islets.groovy',
        'refine_endocrine_pixel_classifier.py', 'vectorize_endocrine_pixel_regions.groovy',
        'build_size_adaptive_islets.groovy', 'size_adaptive_islet_geometry.groovy',
        'import_endocrine_pixel_review.groovy', 'verify_size_adaptive_output.py',
        'export_exact_reviewed_geometry.groovy', 'rebuild_reviewed_islet_regions.py',
        'render_endocrine_rollout.py', 'finalize_endocrine_rollout.py']
    record = dict(version=1, recorded_at_utc=datetime.now(timezone.utc).isoformat(),
        status='native_six_donor_review_ready_hpap_boundaries_approved', donors=rows,
        raw_cells=sum(r['cell_count'] for r in rows), accepted_boundaries=sum(r['accepted_boundaries'] for r in rows),
        awaiting_review_candidates=sum(r['unreviewed_candidates'] for r in rows),
        reviewed_source_native_unchanged=True, reviewed_source_project_unchanged=True,
        original_processing_native_files_unchanged=True, phenotype_gate_approval=False,
        app_updated=False, offline_regional_verification=str(root / 'regions/78/verification.json'),
        regional_verification_sha256=digest(root / 'regions/78/verification.json'),
        native_project_load=loaded, independent_geometry_checks=checks, software_versions=versions,
        measured_storage=storage, storage_excludes_external_source_images=True,
        prediction_seconds_by_donor=prediction_times, adaptive_boundary_seconds_by_donor=boundary_times,
        hp_ap_regional_builder_seconds=read(root / 'regions/78/analysis/78-complete.json')['seconds'],
        configured_java_heap_gib=dict(promotion=24, native_import_each=24, vectorization_each=16, adaptive_each=12),
        peak_rss_measured=False,
        script_sha256={name: digest(source_dir / name) for name in script_names},
        artifacts={str(p.relative_to(project)): digest(p) for p in sorted(project.rglob('*'))
                   if p.is_file() and p.name not in ['.import-metadata.lock']},
        gallery_url=service['url'], http_service_record=str(root / 'review_http.json'))
    write_json(final, record)
    print(json.dumps(dict(status=record['status'], raw_cells=record['raw_cells'],
        accepted_boundaries=record['accepted_boundaries'], awaiting_review_candidates=record['awaiting_review_candidates'],
        storage=storage, gallery_url=record['gallery_url']), indent=2), flush=True)


if __name__ == '__main__':
    main()
