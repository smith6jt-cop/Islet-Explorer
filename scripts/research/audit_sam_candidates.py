#!/usr/bin/env python3
"""Trace similarly sized candidate masks to SAM prompts, without changing objects."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from common import sha256, write_json
from islet_identity_review import load_islets, raw_contact_sheet


def geometry_key(geometry):
    # QuPath removes collinear vertices when importing polygons. Those vertex
    # changes must not be mistaken for a change in location or mask shape.
    return shapely.to_wkb(shapely.normalize(geometry.simplify(0)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--identity-review', required=True, type=Path)
    parser.add_argument('--review-project', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    args = parser.parse_args()
    root = args.input_root.resolve(strict=True)
    out = args.output_root.resolve()
    if out.exists():
        raise ValueError('Preserve existing review evidence; choose a fresh output folder')
    out.mkdir(parents=True)
    sources = pd.read_csv(root/'source_manifest.csv', dtype=str).set_index('case_id')
    review = json.loads((args.review_project/'review_manifest.json').read_text())
    summaries = []
    for donor in review['donors']:
        cid = donor['case_id']
        exported = Path(donor['source_project']).parent.parent/'exports'/f'{int(cid):04d}'
        sam_path, prompt_path = exported/'sam_islets.geojson', exported/'sam_prompts.json'
        sam = json.loads(sam_path.read_text())
        prompt_data = json.loads(prompt_path.read_text())
        downsample = prompt_data['downsample']
        prompts = {p['index']: p for p in prompt_data['boxes']}
        masks = {}
        for feature in sam['features']:
            geom = shapely.geometry.shape(feature['geometry'])
            masks.setdefault(geometry_key(geom), []).append(feature['properties'])
        annotations = load_islets(sources.loc[cid, 'geometry'])
        data = pd.read_csv(args.identity_review/cid/'islet_review.csv', dtype={'case_id': str})
        metrics = []
        for row in data.itertuples(index=False):
            geom = shapely.geometry.shape(annotations[row.islet_key]['geometry'])
            matching = masks.get(geometry_key(geom), [])
            record = dict(islet_id=row.islet_id, area_px2=geom.area,
                          width_px=geom.bounds[2]-geom.bounds[0], height_px=geom.bounds[3]-geom.bounds[1],
                          circularity=4*np.pi*geom.area/geom.length**2,
                          sam_geometry_matches=len(matching), sam_quality_corrected=np.nan,
                          source_prompt_index=None, prompt_width_px=np.nan, prompt_height_px=np.nan,
                          mask_prompt_box_iou=np.nan)
            if len(matching) == 1:
                prop = matching[0]
                prompt = prompts[prop['source_prompt_index']]
                box = shapely.box(*(np.asarray(prompt['bbox'])*downsample))
                record.update(sam_quality_corrected=prop['sam_quality'], source_prompt_index=prompt['index'],
                              source_candidate_uuid=prompt['source_object_id'],
                              prompt_width_px=box.bounds[2]-box.bounds[0], prompt_height_px=box.bounds[3]-box.bounds[1],
                              mask_prompt_box_iou=geom.intersection(box).area/geom.union(box).area)
            metrics.append(record)
        data = data.merge(pd.DataFrame(metrics), on='islet_id', validate='one_to_one')
        data['equivalent_diameter_um'] = 2*np.sqrt(data.area_um2/np.pi)
        data['bbox_width_um'] = data.width_px*.3774
        data['bbox_height_um'] = data.height_px*.3774
        folder = out/cid
        folder.mkdir()
        data.to_csv(folder/'candidate_shape_and_prompt_audit.csv', index=False)
        zero = data.loc[data.hormone_positive_union == 0]
        repeated = zero.groupby(['width_px', 'height_px']).agg(
            candidates=('islet_id', 'size'), islet_names=('islet_key', lambda v: ';'.join(sorted(v))))
        repeated.sort_values(['candidates', 'width_px', 'height_px'], ascending=[False, True, True]).to_csv(folder/'repeated_bounding_boxes.csv')
        summary = dict(case_id=cid, candidates=len(data), sam_masks=len(sam['features']), prompts=len(prompts),
                       uniquely_matching_sam_geometries=int(data.sam_geometry_matches.eq(1).sum()),
                       ambiguous_geometry_matches=int(data.sam_geometry_matches.gt(1).sum()),
                       unmatched_candidates=int(data.sam_geometry_matches.eq(0).sum()),
                       zero_hormone_candidates=len(zero),
                       median_zero_hormone_width_um=float(zero.bbox_width_um.median()),
                       median_zero_hormone_circularity=float(zero.circularity.median()),
                       source_geometry_sha256=sha256(sources.loc[cid, 'geometry']),
                       sam_sha256=sha256(sam_path), prompt_sha256=sha256(prompt_path))
        summaries.append(summary)
        if cid == '78':
            # Select two examples from each of the three most frequent small
            # bounding boxes, sorted by ID. Selection does not use disease or
            # SAM confidence and does not make biological identity decisions.
            groups = repeated.reset_index().sort_values(
                ['candidates', 'width_px', 'height_px'], ascending=[False, True, True]).head(3)
            examples = pd.concat([zero.loc[(zero.width_px == g.width_px) &
                                          (zero.height_px == g.height_px)].sort_values('islet_key').head(2)
                                  for g in groups.itertuples(index=False)])
            examples.to_csv(folder/'small_candidate_examples.csv', index=False)
            raw_contact_sheet(donor['source_image'], annotations, examples.to_dict('records'), folder/'small_candidate_raw_channels')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 2, figsize=(11, 4))
            for label, subset, color in [('No hormone-positive cells', zero, '#c76724'),
                                          ('At least one hormone-positive cell', data.loc[data.hormone_positive_union > 0], '#496e83')]:
                axes[0].hist(subset.equivalent_diameter_um, bins=np.arange(0, 126, 3), histtype='step', label=label, color=color)
                axes[1].scatter(subset.prompt_width_px*.3774, subset.bbox_width_um, s=9, alpha=.35, color=color)
            axes[0].set(xlabel='Equivalent circular diameter (µm)', ylabel='Candidate count')
            axes[0].legend(fontsize=8)
            axes[1].set(xlabel='SAM prompt width (µm)', ylabel='Final mask width (µm)')
            fig.suptitle('HPAP-078 candidate shapes: exported threshold flags are review priorities')
            fig.tight_layout()
            for suffix in ['png', 'pdf']:
                fig.savefig(folder/f'candidate_sizes.{suffix}', dpi=160)
            plt.close(fig)
        print(cid, summary['uniquely_matching_sam_geometries'], 'of', len(data), 'uniquely match SAM', flush=True)
    pd.DataFrame(summaries).to_csv(out/'sam_geometry_alignment_corrected.csv', index=False)
    write_json(out/'audit.json', dict(donors=summaries, source_candidate=str(root),
        previous_identity_review=str(args.identity_review.resolve()), objects_modified=False,
        annotation_decisions_applied=False, script_sha256=sha256(__file__),
        correction='The first identity audit compared serialized polygon vertices. QuPath removes collinear vertices; canonical comparison now simplifies with tolerance zero before normalization. Use this table for SAM alignment, not the earlier sam_geometry_matches column.',
        finding='Many small, round candidate masks lack cells above the exported hormone thresholds. The classifier supplies boxes and SAM returns masks; the historical import did not validate islet identity.',
        interpretation='Repeated dimensions partly reflect the 6-pixel SAM grid. Similar size, roundness, SAM quality and threshold flags cannot alone determine biological identity.'))
    (out/'README.md').write_text('''# Candidate mask review

HPAP-078 small-candidate examples are in `78/small_candidate_raw_channels.pdf`
and PNG, with the original full-resolution channels and saved mask outlines.
The six names and measurements are in `78/small_candidate_examples.csv`.
`78/repeated_bounding_boxes.csv` lists all repeated sizes among the hormone-negative
export flags. `candidate_shape_and_prompt_audit.csv` traces each mask back to its
classifier candidate UUID and SAM box where geometry matches uniquely.

The historical importer accepted masks based on minimal geometry checks, without
confirming islet identity. Roundness and repeated dimensions are diagnostic clues,
not an identity classifier. The 6-pixel sampling grid quantizes outline dimensions.
No objects were deleted or reclassified, no current QuPath review data were changed,
and no absent-hormone rule was used to remove potential depleted islets.

This supplement corrects the earlier SAM alignment diagnostic: QuPath removes
collinear polygon vertices. `sam_geometry_alignment_corrected.csv` compares
normalized geometries after zero-tolerance simplification. The earlier strict
vertex comparison understated matches; it was not evidence of image misalignment.

Review INS/GCG/SST together with CHGA, DAPI and tissue context; record confirmed
islet, non-islet or uncertain decisions in the parent review decision table.
The previous trajectory remains provisional. Corrected region eligibility needs
a fresh fit; edited geometry also needs new membership and compartment measurements.
''')


if __name__ == '__main__':
    main()
