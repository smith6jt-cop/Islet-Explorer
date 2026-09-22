#!/usr/bin/env python3
"""Save size-adaptive boundary comparisons and audits outside the application."""
import argparse
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import shapely
from shapely.geometry import shape
import tifffile

from audit_endocrine_boundaries import raw_crop
from prepare_endocrine_candidates import digest, write_json
from review_hpap_boundary_strategy import geometry_lines


def load(path, classification=None):
    features = json.loads(path.read_text())['features']
    return [shape(f['geometry']) for f in features
            if classification is None or f['properties'].get('classification', {}).get('name') == classification]


def save(fig, out, name):
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=150)
    plt.close(fig)


def overlay(ax, geometries, tree, bounds, offset, color='cyan', scale=1.):
    for i in tree.query(shapely.box(*bounds), predicate='intersects'):
        geom = geometries[i]
        if scale != 1:
            geom = shapely.transform(geom, lambda xy: xy / scale)
        geometry_lines(ax, geom, np.asarray(offset), color)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['adaptive-root', 'refinement-root', 'sam-root', 'selection', 'output-root']:
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    out = args.output_root.resolve()
    out.mkdir(parents=True, exist_ok=False)
    candidates = args.adaptive_root / 'candidates'
    variants = [args.refinement_root / 'sensitivity/78-unanimous',
                candidates / '78-unanimous', candidates / '78-majority']
    manifests = [json.loads((root / 'prediction_manifest.json').read_text()) for root in variants]
    if len({r['model_sha256'] for r in manifests}) != 1:
        raise ValueError('Comparison must hold the fitted model fixed')
    geoms = [load(root / 'review_candidates.geojson') for root in variants]
    trees = [shapely.STRtree(g) for g in geoms]
    spec = json.loads((args.sam_root / '78/prompts.json').read_text())
    pixels = tifffile.memmap(spec['source_image'], mode='r')
    previous = {f['properties']['name']: shape(f['geometry']) for f in
                json.loads((args.sam_root / '78/review_candidates.geojson').read_text())['features']}
    selection = pd.read_csv(args.selection)
    pages, comparisons = [], []
    for start in range(0, len(selection), 3):
        rows = selection.iloc[start:start+3]
        fig, axes = plt.subplots(len(rows), 4, figsize=(15, 3.5 * len(rows)), squeeze=False)
        for row, entry in enumerate(rows.itertuples()):
            old = previous[entry.candidate_id]
            raw, _, bounds = raw_crop(pixels, spec, old, pad_um=20)
            x0, y0, x1, y1 = bounds
            rgb = np.clip((raw-spec['display_min']) /
                          (np.asarray(spec['display_max'])-spec['display_min']), 0, 1)
            for col in range(4):
                ax = axes[row, col]
                ax.imshow(rgb)
                if col:
                    overlay(ax, geoms[col-1], trees[col-1], bounds, (x0, y0))
                ax.set_xlim(0, x1-x0)
                ax.set_ylim(y1-y0, 0)
                ax.set_xticks([])
                ax.set_yticks([])
                if row == 0:
                    ax.set_title(['Raw INS/GCG/SST', 'Previous uniform cleanup\n160/160 votes',
                                  'Size-adaptive boundaries\n160/160 votes',
                                  'Size-adaptive boundaries\n80/160 votes'][col], fontsize=10)
                if col == 0:
                    ax.set_ylabel(entry.candidate_id, fontsize=9)
            comparisons.append(dict(reference_location=entry.candidate_id,
                old_area_um2=old.area*spec['pixel_size_um']**2,
                **{name:len(tree.query(old, predicate='intersects')) for name, tree in
                   zip(['previous_regions', 'adaptive_unanimous_regions', 'adaptive_majority_regions'], trees)},
                review_status='unreviewed'))
        fig.suptitle('HPAP-078: fixed comparison locations; classification is unchanged within each cutoff\n'
                     'Boundaries need image review. Filled interiors are not measurements of hormone-positive area.', fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .94))
        name = f'hpap_comparison_{start//3+1:02d}'
        save(fig, out, name)
        pages.append(name)
        print(f'HPAP comparison {start+len(rows)}/{len(selection)}', flush=True)
    pd.DataFrame(comparisons).to_csv(out / 'hpap_comparison_metrics.csv', index=False)
    with gzip.open(spec['source_geometry'], 'rt') as stream:
        old_native = {f['properties'].get('name'): shape(f['geometry']) for f in json.load(stream)['features']
                      if f['properties'].get('classification', {}).get('name') == 'Islet'}
    unions = [shapely.union_all(g) for g in geoms]
    negatives = [dict(previous_islet=name, **{key:old_native[name].intersection(g).area / old_native[name].area
                 for key, g in zip(['previous_fraction', 'adaptive_unanimous_fraction', 'adaptive_majority_fraction'], unions)})
                 for name in ['Islet_130', 'Islet_206', 'Islet_1148', 'Islet_124', 'Islet_161', 'Islet_179']]
    pd.DataFrame(negatives).to_csv(out / 'previous_negative_controls.csv', index=False)
    # Full field reveals new grouping and any islands outside the fixed examples.
    raw = np.stack([pixels[spec['source_channel_names'].index(m), ::8, ::8] for m in ['INS','GCG','SST']], -1)
    rgb = np.clip((raw-spec['display_min'])/(np.asarray(spec['display_max'])-spec['display_min']), 0, 1)
    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    bounds=(0, 0, pixels.shape[2], pixels.shape[1])
    for col, ax in enumerate(axes):
        ax.imshow(rgb)
        if col:
            overlay(ax, geoms[col], trees[col], bounds, (0, 0), scale=8)
        ax.set_axis_off()
        ax.set_title(['Raw hormones', 'Adaptive 160/160 votes', 'Adaptive 80/160 votes'][col])
    fig.tight_layout()
    save(fig, out, 'hpap_whole_image')
    reference_check = json.loads((args.adaptive_root / 'reference_geometry/verification.json').read_text())
    if reference_check['status'] != 'passed' or reference_check['donor_count'] != 15:
        raise ValueError('Incomplete original-donor geometry check')
    if {m['boundary_method_reference']['helper_sha256'] for m in manifests[1:]} != {reference_check['geometry_helper_sha256']}:
        raise ValueError('Different boundary rules used in reference and target checks')
    ref_pages, summaries = [], []
    for record in reference_check['donors']:
        cid = record['case_id']
        source = Path(record['source_reference'])
        if digest(source) != record['source_reference_sha256'] or digest(Path(record['output'])) != record['output_sha256']:
            raise ValueError('Reference comparison geometry changed')
        original = load(source, 'Islet')
        revised = load(Path(record['output']))
        old_tree, new_tree = shapely.STRtree(original), shapely.STRtree(revised)
        root = args.refinement_root / 'samples_v2' / cid
        raw = np.load(root / 'endocrine_overview.npy')
        mask = np.asarray(Image.open(root / 'source_islet_regions.png')) > 0
        high = np.maximum(np.quantile(raw[mask], .995, axis=0), 1.)
        low = np.quantile(raw.reshape(-1, 3), .05, axis=0)
        rgb = np.clip((raw-low)/np.maximum(high-low, 1), 0, 1)
        bins = [(0, 2000, 'Small'), (2000, 10000, 'Medium'), (10000, np.inf, 'Large')]
        selected=[]
        for low_area, high_area, label in bins:
            subset=sorted([g for g in original if low_area <= g.area*.3774**2 < high_area], key=lambda g:g.area)
            selected.append((subset[len(subset)//2] if subset else None, label))
        fig, axes=plt.subplots(2, 3, figsize=(11, 7))
        for col, (geom, label) in enumerate(selected):
            if geom is None:
                for ax in axes[:, col]:
                    ax.text(.5, .5, f'No {label.lower()} reference', ha='center', va='center', transform=ax.transAxes)
                    ax.set_axis_off()
                continue
            x0,y0,x1,y1=geom.bounds
            bounds=(max(0, int(x0/4)-16), max(0, int(y0/4)-16), min(raw.shape[1], int(x1/4)+17), min(raw.shape[0], int(y1/4)+17))
            a,b,c,d=bounds
            source_bounds=tuple(x*4 for x in bounds)
            for row, (gs, tree, color) in enumerate([(original, old_tree, 'yellow'), (revised, new_tree, 'cyan')]):
                ax=axes[row,col]
                ax.imshow(rgb[b:d,a:c])
                overlay(ax, gs, tree, source_bounds, (a,b), color, scale=4)
                ax.set_xlim(0,c-a); ax.set_ylim(d-b,0);ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(f'{label}: {geom.area*.3774**2:,.0f} µm²')
                if col==0:ax.set_ylabel(['Source reference','Size-adaptive geometry'][row])
        fig.suptitle(f'Original donor {cid}: size-stratified geometry check\nExisting reference shapes as input; no new classification or biological accuracy estimate.', fontsize=11)
        fig.tight_layout(rect=(0,0,1,.92)); name=f'reference_{cid}';save(fig,out,name);ref_pages.append(name)
        summaries.append(dict(case_id=cid, input_kind='preserved reference geometry', **record['summary']))
        print(f'Reference comparisons {cid}', flush=True)
    for key, manifest in zip(['78-unanimous','78-majority'],manifests[1:]):
        summaries.append(dict(case_id=key,input_kind='fixed endocrine classifier foreground',**manifest['boundary_summary']))
        pd.DataFrame(json.loads((candidates/key/'parent_audit.json').read_text())).to_csv(out/f'{key}_parent_audit.csv',index=False)
    for row in summaries:
        for key in ['parents_by_size','output_regions_by_parent_size']:
            row[key]=json.dumps(row[key],sort_keys=True)
    pd.DataFrame(summaries).to_csv(out/'geometry_summary.csv',index=False)
    failure=previous['New_78_0242']
    target_features=json.loads((candidates/'78-unanimous/review_candidates.geojson').read_text())['features']
    matching=[f['properties']['name'] for f in target_features if shape(f['geometry']).intersects(failure)]
    result=dict(status='size_adaptive_boundary_review',HPAP_current_failure_candidates=matching,
        HPAP_unanimous_candidate_count=len(geoms[1]),HPAP_majority_candidate_count=len(geoms[2]),
        original_donor_geometry_checks=15,fixed_hpap_locations=len(selection),
        negative_controls=negatives,accepted_islets=0,other_new_donors_regenerated=False,
        model_sha256=manifests[1]['model_sha256'],geometry_helper_sha256=reference_check['geometry_helper_sha256'])
    write_json(out/'review_manifest.json',result)
    def page(name):return f'<a href="{name}.pdf"><img loading="lazy" src="{name}.png"></a>'
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Size-adaptive endocrine boundaries</title>'
        '<style>body{font:16px system-ui;max-width:1500px;margin:2rem auto;padding:1rem}img{max-width:100%}a{color:#075ac0}</style>'
        '<h1>Size-adaptive endocrine boundaries</h1>'
        '<p>Rules from IO60panc2nd/scripts/1_panc_objects.groovy: erosion of 1, 2 or 3 µm for regions below 2,000, '
        '2,000–9,999 or at least 10,000 µm²; 85% restoration; 50 µm² fragments; enclosed-hole filling; 5 µm image-edge exclusion. '
        'This cohort retains its 0.3774 µm/pixel calibration. Tissue exclusions remain empty; colliding restoration areas are excluded from both outlines.</p>'
        '<p><strong>Review candidates:</strong> the endocrine classifier and its two previous cutoffs are unchanged. '
        'The IO60 Islets2 classifier requires CD99, absent from this panel, so its model was not copied. '
        'Morphology does not validate weak staining or the identity of small endocrine foci. Filled outlines are not hormone-positive area measurements.</p>'
        f'<p>HPAP-078: {len(geoms[1]):,} conservative candidates and {len(geoms[2]):,} majority-cutoff candidates. '
        'The original 15 donors below test geometry handling on preserved reference shapes, without rerunning classification. Other new donors have not been regenerated.</p>'
        '<p>In QuPath, choose <strong>File → Project → Open project</strong> and open '
        '<code>/data/panc_CODEX/Islet-Explorer/data/reviews/qupath-size-adaptive-islets-20260917/project.qpproj</code>. '
        'Search for <strong>AdaptiveIslet_78_00853</strong> to inspect the reported failure location.</p>'
        '<p><a href="../qupath-size-adaptive-islets-20260917/README.md">QuPath opening instructions</a> · '
        '<a href="geometry_summary.csv">Geometry summary</a> · <a href="78-unanimous_parent_audit.csv">Conservative region audit</a> · '
        '<a href="78-majority_parent_audit.csv">Majority region audit</a> · <a href="previous_negative_controls.csv">Previous negative examples</a> · '
        '<a href="review_manifest.json">Review record</a></p>'
        '<h2>Full HPAP-078 image</h2>'+page('hpap_whole_image')+
        ''.join(f'<h2>HPAP comparison {i+1}</h2>'+page(name) for i,name in enumerate(pages))+
        '<h2>Original donor size checks</h2>'+''.join(page(name) for name in ref_pages)+'</html>')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
