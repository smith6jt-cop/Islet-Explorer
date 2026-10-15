#!/usr/bin/env python3
"""Audit endocrine support of candidate islets; flags are not biological decisions."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
import tifffile
from common import sha256, write_json


def load_islets(path):
    with gzip.open(path, 'rt') as stream:
        doc = json.load(stream)
    return {f['properties']['name']: f for f in doc['features']
            if f.get('properties', {}).get('classification', {}).get('name') == 'Islet'}


def raw_contact_sheet(image, annotations, rows, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with tifffile.TiffFile(image) as tf:
        labels = [s.removesuffix('.tif').removesuffix('.tiff') for s in tf.imagej_metadata['Labels']]
    pixels = tifffile.memmap(image, mode='r')
    markers = ['DAPI', 'INS', 'GCG', 'SST', 'CHGA']
    levels = {m: float(np.quantile(pixels[labels.index(m), ::16, ::16], .998)) for m in markers}
    fig, axes = plt.subplots(len(rows), len(markers), figsize=(14, 2.7*len(rows)), squeeze=False)
    crop_records = []
    for rownum, row in enumerate(rows):
        geom = shapely.geometry.shape(annotations[row['islet_key']]['geometry'])
        cx, cy = geom.centroid.x, geom.centroid.y
        half = max(110, int(max(geom.bounds[2]-geom.bounds[0], geom.bounds[3]-geom.bounds[1])/2)+40)
        x0, y0 = max(0, int(cx)-half), max(0, int(cy)-half)
        x1, y1 = min(pixels.shape[2], int(cx)+half), min(pixels.shape[1], int(cy)+half)
        for col, marker in enumerate(markers):
            ax = axes[rownum, col]
            ax.imshow(pixels[labels.index(marker), y0:y1, x0:x1], cmap='gray', vmin=0, vmax=levels[marker])
            for part in shapely.get_parts(geom):
                xy = np.asarray(part.exterior.coords)
                ax.plot(xy[:, 0]-x0, xy[:, 1]-y0, color='#ff7f00', linewidth=.65)
            ax.set_xticks([]); ax.set_yticks([])
            if rownum == 0: ax.set_title(marker)
            if col == 0:
                ax.set_ylabel(f"{row['islet_key']}\n{row['hormone_positive_union']} hormone+ cells", fontsize=8)
        crop_records.append(dict(islet_key=row['islet_key'], x_px=x0, y_px=y0, width_px=x1-x0, height_px=y1-y0))
    fig.suptitle('HPAP-078: raw channels with saved islet outlines\nFixed per-channel display levels across rows; orange outlines are candidate regions', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .95))
    for suffix in ['png', 'pdf']: fig.savefig(destination.with_suffix('.'+suffix), dpi=160)
    plt.close(fig)
    write_json(destination.with_suffix('.json'), dict(source_image=str(image), display_max=levels,
        display_method='99.8th percentile sampled every 16 source pixels; no per-crop contrast fitting', crops=crop_records))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    p.add_argument('--review-project', required=True, type=Path)
    a = p.parse_args(); root = a.input_root.resolve(strict=True); out = a.output_root.resolve()
    if out.exists(): raise ValueError('Preserve existing review evidence; use a new output folder')
    out.mkdir(parents=True)
    cohort = pd.read_csv(root/'source_manifest.csv', dtype=str)
    features = pd.read_csv(root/'trajectory/features.csv', dtype={'case_id': str})
    review = json.loads((a.review_project/'review_manifest.json').read_text())
    review_donors = {d['case_id']: d for d in review['donors']}
    tables, summaries = [], []
    for donor in cohort.itertuples(index=False):
        cid = donor.case_id
        source = root/'regional'/f'islet_composition_{int(cid):04d}.tsv'
        comp = pd.read_csv(source, sep='\t')
        comp['islet_id'] = cid+'_'+comp.name
        comp['hormone_positive_union'] = comp.cells_total-comp.Triple_Neg
        data = features.loc[features.case_id == cid, ['islet_id', 'case_id', 'islet_key', 'area_um2', 'measured_cells', 'eligible', 'centroid_x_um', 'centroid_y_um']].merge(
            comp[['islet_id', 'cells_total', 'hormone_positive_union', 'Ins_any', 'Glu_any', 'Stt_any']], on='islet_id', validate='one_to_one')
        data['review_flag'] = np.where(data.hormone_positive_union == 0, 'no_hormone_positive_cells',
                                     np.where(data.hormone_positive_union < 5, 'one_to_four_hormone_positive_cells', 'hormone_positive_cells_present'))
        data['identity_decision'] = 'unreviewed'
        annotations = load_islets(donor.geometry)
        if set(annotations) != set(data.islet_key): raise ValueError(f'Geometry identity mismatch: {cid}')
        data['sam_quality'] = np.nan
        if cid in review_donors:
            source_project = Path(review_donors[cid]['source_project']).parent.parent
            tag = '0078' if cid == '78' else cid
            sam = json.loads((source_project/'exports'/tag/'sam_islets.geojson').read_text())
            masks = {}
            for f in sam['features']:
                key = shapely.to_wkb(shapely.normalize(shapely.geometry.shape(f['geometry']).simplify(0)))
                masks.setdefault(key, []).append(f['properties'])
            for i, row in data.iterrows():
                key = shapely.to_wkb(shapely.normalize(shapely.geometry.shape(annotations[row.islet_key]['geometry']).simplify(0)))
                if key in masks: data.loc[i, 'sam_quality'] = min(v['sam_quality'] for v in masks[key])
        folder = out/cid; folder.mkdir()
        data.to_csv(folder/'islet_review.csv', index=False)
        flagged = []
        for row in data.loc[data.hormone_positive_union == 0].itertuples(index=False):
            original = annotations[row.islet_key]
            flagged.append(dict(type='Feature', geometry=original['geometry'], properties=dict(
                objectType='annotation', name='Review only: '+row.islet_key,
                classification=dict(name='QC: endocrine support pending', color=[255, 128, 0]),
                measurements={'Source islet area um2': row.area_um2, 'Source hormone-positive cells': 0},
                source_islet_id=row.islet_id, review_decision='pending; zero threshold-positive cells does not establish non-islet identity')))
        write_json(folder/'flagged_islets.geojson', dict(type='FeatureCollection', features=flagged))
        if cid == '78':
            empty = data.loc[data.hormone_positive_union == 0].sort_values(['area_um2', 'islet_key'], ascending=[False, True]).head(4)
            positive = data.sort_values(['hormone_positive_union', 'islet_key'], ascending=[False, True]).head(1)
            raw_contact_sheet(review_donors[cid]['source_image'], annotations,
                              pd.concat([empty, positive]).to_dict('records'), folder/'raw_channel_review')
        summaries.append(dict(case_id=cid, islets=len(data), no_hormone_positive_cells=int((data.hormone_positive_union == 0).sum()),
            trajectory_eligible=int(data.eligible.sum()), eligible_without_hormone_positive_cells=int((data.eligible & (data.hormone_positive_union == 0)).sum()),
            source_composition_sha256=sha256(source), source_geometry_sha256=sha256(donor.geometry),
            sam_geometry_matches=int(data.sam_quality.notna().sum()), decision='requires islet-identity review'))
        tables.append(data); print(cid, summaries[-1]['no_hormone_positive_cells'], 'flagged for review', flush=True)
    combined = pd.concat(tables)
    combined.to_csv(out/'islet_review_all_donors.csv', index=False)
    decisions = combined[['islet_id', 'case_id', 'islet_key', 'review_flag']].copy()
    decisions['decision'] = 'unreviewed'
    decisions['geometry_sha256'] = decisions.case_id.map({r['case_id']: r['source_geometry_sha256'] for r in summaries})
    decisions['review_note'] = ''
    decisions.to_csv(out/'islet_review_decisions.csv', index=False)
    pd.DataFrame(summaries).to_csv(out/'donor_summary.csv', index=False)
    write_json(out/'review.json', dict(status='identity_review_required', donors=summaries,
        identified_gap='SAM import checks geometry size but does not confirm endocrine/islet identity; trajectory eligibility checks geometry and measured-cell count.',
        automatic_exclusions_applied=False, review_project_modified=False,
        sam_geometry_matching_method='normalized polygons after zero-tolerance simplification of collinear vertices',
        interpretation='Flags prioritize review; absence of threshold-positive cells alone does not prove that a candidate is not an islet.',
        trajectory_action='Current coordinates are provisional. Confirmed region/eligibility corrections require a new offline fit and validation, preserving the previous fit.'))


if __name__ == '__main__': main()
