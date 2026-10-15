#!/usr/bin/env python3
"""Validate and label the independent native review delivery; never accept islets."""
import argparse
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from PIL import Image
import shapely
from shapely.geometry import shape
import tifffile

from prepare_endocrine_candidates import digest, write_json, refuse_withdrawn_candidate
from convert_images import canonical_channel


def verify_images(candidate, project, inventory):
    sources = {str(d['case_id']): d for d in json.loads(inventory.read_text())['cases']}
    donors = json.loads((project/'review_manifest.json').read_text())['donors']
    target = candidate/'source-image-verification.json'
    previous = json.loads(target.read_text()) if target.exists() else None
    if previous and previous['source_inventory_sha256'] != digest(inventory):
        raise ValueError('Source inventory changed')
    records = []
    for donor in donors:
        cid = donor['case_id']
        spec = json.loads((candidate/cid/'prompts.json').read_text())
        path = Path(spec['source_image']).resolve(strict=True)
        expected = sources[cid]
        if path != Path(expected['source']).resolve(strict=True):
            raise ValueError('Candidate image differs from selected reconciled acquisition')
        stat = path.stat()
        report = json.loads(Path(expected['conversion_report']['path']).read_text())
        if stat.st_size != report['source_bytes'] or stat.st_mtime_ns != report['source_mtime_ns']:
            raise ValueError('Source image fingerprint changed since reconciled acquisition verification')
        if previous:
            record = next(r for r in previous['images'] if r['case_id'] == cid)
            if record['bytes'] != stat.st_size or record['mtime_ns'] != stat.st_mtime_ns or record['previously_verified_sha256'] != expected['source_sha256']:
                raise ValueError('Image changed since this candidate source verification')
        else:
            # Verify the exact pixels used for this candidate. Full TIFF hashes are
            # linked to the prior image inventory, not misreported as newly recomputed.
            with tifffile.TiffFile(path) as tf:
                channels = [canonical_channel(s) for s in tf.imagej_metadata['Labels']]
            pixels = tifffile.memmap(path, mode='r')
            if list(pixels.shape) != spec['image_shape_cyx']:
                raise ValueError('Source image dimensions changed')
            count = 0
            for prompt in spec['prompts']:
                if 'crop' not in prompt: continue
                rendered = np.asarray(Image.open(candidate/cid/prompt['crop']).convert('RGB'))
                x, y = prompt['offset_px']; step = prompt['downsample']
                raw = np.stack([pixels[channels.index(m), y:y+rendered.shape[0]*step:step,
                                       x:x+rendered.shape[1]*step:step] for m in spec['rgb_channels']], axis=-1)
                low, high = np.asarray(spec['display_min']), np.asarray(spec['display_max'])
                reproduced = np.uint8(np.clip((raw-low)/(high-low), 0, 1)*255)
                if not np.array_equal(rendered, reproduced):
                    raise ValueError(f'SAM crop no longer matches raw pixels: {cid} {prompt["candidate_id"]}')
                count += 1
            if path.stat().st_mtime_ns != stat.st_mtime_ns:
                raise ValueError('Source changed during pixel verification')
            record = dict(case_id=cid, path=str(path), bytes=stat.st_size,
                          mtime_ns=stat.st_mtime_ns, previously_verified_sha256=expected['source_sha256'],
                          full_image_sha256_recomputed=False, verified_sam_crops=count,
                          every_sam_crop_matches_current_raw_image=True,
                          checked_at_utc=datetime.now(timezone.utc).isoformat())
            print(f'{cid}: {count} SAM crops match current raw pixels; source path/size/mtime match reconciled acquisition', flush=True)
        records.append(record)
    if len({r['previously_verified_sha256'] for r in records}) != len(records):
        raise ValueError('Duplicate source images')
    if not previous:
        write_json(target, dict(status='passed', source_inventory=str(inventory),
                               source_inventory_sha256=digest(inventory), images=records,
                               verification_scope='Fresh equality check of every SAM crop against raw image pixels; source fingerprints against prior full-image hashes',
                               full_image_sha256_recomputed=False))


def package(candidate, project):
    refuse_withdrawn_candidate(candidate)
    refuse_withdrawn_candidate(project)
    destination = project/'regenerated_review_manifest.json'
    if destination.exists(): raise ValueError('Preserve delivered review state')
    original = json.loads((project/'review_manifest.json').read_text())
    verified = json.loads((project/'regeneration-verification/native-summary.json').read_text())
    if verified['status'] != 'passed': raise ValueError('Native verification is incomplete')
    expected = {d['case_id'] for d in original['donors']}
    checks = {r['case_id']: r for r in verified['donors']}
    if set(checks) != expected: raise ValueError('Native verification cohort differs')
    rows, classes, artifacts = [], {}, []
    for donor in original['donors']:
        cid = donor['case_id']
        folder, proof = candidate/cid, checks[cid]
        if digest(donor['review_data']) != proof['new_data_sha256']:
            raise ValueError('Preserve changed review data')
        if digest(donor['source_data']) != donor['data_sha256']:
            raise ValueError('Source processing objects changed')
        summary = json.loads((folder/'candidate_summary.json').read_text())
        spec = json.loads((folder/'prompts.json').read_text())
        fs = json.loads((folder/'review_candidates.geojson').read_text())['features']
        if digest(folder/'review_candidates.geojson') != proof['candidate_overlay_sha256']:
            raise ValueError('Native objects no longer match source overlay')
        if len({f['properties']['name'] for f in fs}) != len(fs):
            raise ValueError('Duplicate candidate names')
        geometries = np.array([shape(f['geometry']) for f in fs], dtype=object)
        if not shapely.is_valid(geometries).all(): raise ValueError('Invalid review geometry')
        boundaries = [g for g in geometries if g.geom_type in ['Polygon', 'MultiPolygon']]
        union = shapely.union_all(boundaries)
        if abs(sum(g.area for g in boundaries)-union.area) > 1e-5:
            raise ValueError('Overlapping candidate areas remain')
        tissue = shape(json.loads((folder/'valid_tissue.geojson').read_text()))
        if union.difference(tissue).area > 1e-5: raise ValueError('New boundaries extend outside valid tissue')
        seeds = pd.read_csv(folder/'seed_cells.csv')
        xy = seeds[['X_centroid', 'Y_centroid']].to_numpy()/spec['pixel_size_um']
        pairs = shapely.STRtree(geometries).query(shapely.points(xy), predicate='dwithin', distance=1e-5)
        if len(np.unique(pairs[0])) != len(seeds): raise ValueError('Endocrine seed coverage lost in export')
        for feature in fs:
            cls = feature['properties']['classification']
            classes[cls['name']] = cls['color']
        row = dict(case_id=cid, cell_count=int(proof['cells']), previous_islets=summary['old_islet_count'],
                   new_boundaries=summary['new_boundary_count'],
                   endocrine_supported_candidates=summary['endocrine_supported_candidates'],
                   boundaries_needing_extra_review=summary['boundaries_needing_extra_review'],
                   point_focus_annotations=summary['separate_point_focus_annotations'],
                   hormone_seed_cells=summary['seed_cell_count'], accepted_islets=0,
                   native_data_sha256=proof['new_data_sha256'],
                   all_seeds_represented=True, no_overlapping_boundaries=True,
                   all_boundaries_within_valid_tissue=True)
        rows.append(row)
        for path in folder.iterdir():
            if path.is_file(): artifacts.append(dict(path=str(path), bytes=path.stat().st_size, sha256=digest(path)))
    # The copy manifest remains intact as starting evidence; the new manifest is authoritative for delivery.
    document = json.loads((project/'project.qpproj').read_text())
    for entry in document['images']:
        donor = next(d for d in original['donors'] if d['entry_id'] == entry['entryID'])
        entry.setdefault('metadata', {}).update({
            'Boundary review': 'New endocrine-seeded candidates; unreviewed',
            'Previous segmentation': 'Grey outlines, cell parents, old distances and regional measurements retained',
            'Review guide': str(project/'README.md'),
        })
    write_json(project/'project.qpproj', document)
    palette = json.loads((project/'classifiers/classes.json').read_text())
    classes.update({'Previous segmentation: Islet': [150,150,150], 'Previous segmentation: IsletExpanded': [150,150,150]})
    palette['pathClasses'].extend(dict(name=name, color=(255<<24 | rgb[0]<<16 | rgb[1]<<8 | rgb[2])-(1<<32))
        for name, rgb in classes.items() if name not in {c['name'] for c in palette['pathClasses']})
    write_json(project/'classifiers/classes.json', palette)
    pd.DataFrame(rows).to_csv(project/'regenerated_donors.csv', index=False)
    negative_names = ['Islet_130', 'Islet_206', 'Islet_1148', 'Islet_124', 'Islet_161', 'Islet_179']
    negatives = pd.read_csv(candidate/'78/previous_islet_comparison.csv').set_index('previous_islet').loc[negative_names]
    negative_overlap = negatives.fraction_covered_by_new_boundaries.to_dict()
    negative_sentence = ('All six checked hormone-negative HPAP-078 examples have zero overlap with the new boundaries.'
                         if all(v == 0 for v in negative_overlap.values()) else
                         'Some previous hormone-negative HPAP-078 examples overlap a new boundary; inspect the comparison table.')
    manifest = dict(status='ready_for_human_review', created_at_utc=datetime.now(timezone.utc).isoformat(),
        project=str(project/'project.qpproj'), project_sha256=digest(project/'project.qpproj'),
        candidate_root=str(candidate), donor_count=len(rows), donors=rows,
        total_cells=sum(r['cell_count'] for r in rows), total_boundaries=sum(r['new_boundaries'] for r in rows),
        analytical_acceptance=False, threshold_validation=False, accepted_islets=0,
        current_app_or_analysis_modified=False, prior_review_project_modified=False,
        native_object_bytes=sum(Path(d['review_data']).stat().st_size for d in original['donors']),
        source_copy_manifest=str(project/'review_manifest.json'), artifacts=artifacts,
        source_image_verification=str(candidate/'source-image-verification.json'),
        hpap_negative_control_overlap_fractions=negative_overlap,
        support_coverage_definition='The historical prompted_seed_recall field describes coverage of at most 32 local signal-support locations; box-only inference does not pass those points to SAM',
        source_code=[dict(path=str(path), sha256=digest(path)) for path in sorted(Path(__file__).parent.glob('*endocrine*')) if path.is_file()],
        software=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, shapely=shapely.__version__),
        scientific_limitations=[
            'Thresholds are operational, not biologically validated phenotype gates.',
            'This seed rule can miss hormone-low or hormone-negative endocrine structures; previous outlines remain available.',
            'Local SAM masks can fragment or overextend; all new boundaries require visual identity and boundary review.',
            'The 20-um seed grouping can merge nearby endocrine structures; review grouping before islet-count or size inference.',
            'Single-cell and two-cell foci remain point annotations, not accepted islet boundaries.',
            'Cell hierarchy parents, distances, old compartments and old regional measurements use previous geometry.',
            'The installed CUDA cumsum operation lacks a strict determinism guarantee; that limitation is recorded.',
            'Geometry-dependent downstream measurements and trajectory require a separate rebuild after review.',
        ])
    write_json(destination, manifest)
    table = '\n'.join(f'| {r["case_id"]} | {r["new_boundaries"]} | {r["endocrine_supported_candidates"]} | '
                      f'{r["boundaries_needing_extra_review"]} | {r["point_focus_annotations"]} |' for r in rows)
    (project/'README.md').write_text(f'''# Regenerated islet candidates — review delivery

Open `project.qpproj` in QuPath 0.7. All six donors and {manifest['total_cells']:,} original cell objects are included.
The {manifest['total_boundaries']:,} new boundaries are **candidates for your review**, not accepted islets.

- Green: endocrine-supported boundary candidates (at least three threshold-positive cell seeds).
- Yellow: sparse support, partial masks, crop-edge contact or poor hormone-signal contrast; prioritize for review.
- Blue: endocrine seed points, including single-cell and two-cell foci; these are not islet outlines.
- Grey: previous Islet/IsletExpanded outlines. Hide these classes for a clear view of new boundaries.

| Donor | New boundaries | Green candidates | Extra boundary review | Separate point foci |
|---|---:|---:|---:|---:|
{table}

The prior review project and any edits there are preserved. Each native save was read again and checked:
every cell UUID, cell/nuclear outline, classification, parent ID and measurement is unchanged;
every original annotation geometry and every imported candidate coordinate is preserved.
New boundary seed/total-cell counts are computed by spatial intersection, independently of the old hierarchy.
Cell parents, old distance measurements and previous regional measurements still refer to the old boundaries.

For raw INS/GCG/SST/CHGA plots, open [the review gallery]({candidate}/index.html).
That gallery shows previous outlines in orange and new outlines in green.
{negative_sentence}

Review hormone signal, CHGA and tissue context before accepting candidates. Hormone thresholds have not
been biologically validated, so weak or hormone-negative structures must also be assessed using the retained
previous outlines. No new islet has been accepted automatically. Record decisions in each donor's
`candidate_review.csv` (in the candidate folder) and save edits only in this project.

`regenerated_review_manifest.json` and `regeneration-verification/` describe this delivery.
`review_manifest.json` describes the initial independent copies, before regeneration overlays were added.
The app, trajectory, original processing objects and the frozen release have not been replaced.
''')
    gallery_rows = []
    for row in rows:
        cid = row['case_id']
        figures = ''.join(f'<figure><a href="{cid}/{name}.pdf"><img loading="lazy" src="{cid}/{name}.png" alt="{html.escape(name)} for {cid}"></a><figcaption>{html.escape(name.replace("_", " "))}</figcaption></figure>'
                          for name in ['old_positive_controls', 'old_negative_controls', 'new_candidates_size_range'])
        gallery_rows.append(f'<h2>Donor {cid}</h2><p>{row["new_boundaries"]} boundaries; '
                            f'{row["boundaries_needing_extra_review"]} need extra boundary review. '
                            f'<a href="{cid}/candidate_review.csv">Review decisions CSV</a></p>{figures}')
    (candidate/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<title>Regenerated islet candidate review</title><style>body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:1rem}img{max-width:100%}figure{margin:2rem 0}</style>'
        '<h1>Regenerated islet candidates</h1><p>Ready for human review. No biological acceptance or new trajectory fit. '
        'Orange outlines are previous segmentation; green outlines are newly generated candidates. '
        'Fixed per-image display ranges; click figures for PDF.</p>'+'\n'.join(gallery_rows)+'</html>')
    print(json.dumps(dict(status=manifest['status'], donors=len(rows), cells=manifest['total_cells'],
                         new_boundaries=manifest['total_boundaries'], project=str(project/'project.qpproj'))), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate-root', required=True, type=Path)
    p.add_argument('--project-root', required=True, type=Path)
    p.add_argument('--source-inventory', required=True, type=Path)
    p.add_argument('--verify-source-images-only', action='store_true')
    a = p.parse_args()
    candidate, project = a.candidate_root.resolve(strict=True), a.project_root.resolve(strict=True)
    refuse_withdrawn_candidate(candidate)
    refuse_withdrawn_candidate(project)
    verify_images(candidate, project, a.source_inventory.resolve(strict=True))
    if not a.verify_source_images_only: package(candidate, project)


if __name__ == '__main__': main()
