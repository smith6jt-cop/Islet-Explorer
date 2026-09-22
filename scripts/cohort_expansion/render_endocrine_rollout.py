#!/usr/bin/env python3
"""Offline native-boundary review gallery with deterministic size-stratified crops."""
import argparse
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape
import tifffile

from audit_endocrine_boundaries import raw_crop
from prepare_endocrine_candidates import digest, write_json
from render_size_adaptive_review import overlay, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('rollout-root', 'approved-prediction-root', 'display-root', 'project-root', 'output-root'):
        parser.add_argument('--' + arg, type=Path, required=True)
    args = parser.parse_args()
    rollout = args.rollout_root.resolve(strict=True)
    project = args.project_root.resolve(strict=True)
    out = args.output_root.resolve()
    out.mkdir(parents=True, exist_ok=False)
    staging = json.loads((project / 'staging_manifest.json').read_text())
    approval = json.loads((project / 'boundary_approval.json').read_text())
    promotion = json.loads((rollout / 'exports/78/promotion_manifest.json').read_text())
    summaries, selections, pages = [], [], []
    for donor in staging['donors']:
        cid = donor['case_id']
        approved = cid == approval['case_id']
        candidate_root = args.approved_prediction_root if approved else rollout / 'candidates' / cid
        manifest = json.loads((candidate_root / 'prediction_manifest.json').read_text())
        geometry_path = rollout / 'exports/78/accepted_islets.geojson' if approved else candidate_root / 'review_candidates.geojson'
        features = json.loads(geometry_path.read_text())['features']
        geoms = [shape(f['geometry']) for f in features]
        names = [f['properties']['name'] for f in features]
        tree = shapely.STRtree(geoms)
        spec = json.loads((args.display_root / cid / 'prompts.json').read_text())
        source_spec = Path(manifest['source_specification'])
        if digest(source_spec) != manifest['source_specification_sha256']:
            raise ValueError('Prediction source specification changed')
        fingerprint = json.loads(source_spec.read_text())
        # Old review files supply only the saved display ranges and raw channel
        # order. Their discarded masks are never read by this renderer.
        image = Path(manifest['source_image']).resolve(strict=True)
        stat = image.stat()
        if (image != Path(spec['source_image']).resolve(strict=True)
                or str(image) != str(Path(fingerprint['source_image']).resolve(strict=True))
                or stat.st_size != fingerprint['source_bytes'] or stat.st_mtime_ns != fingerprint['source_mtime_ns']):
            raise ValueError('Raw display image identity differs for ' + cid)
        if spec['rgb_channels'] != ['INS', 'GCG', 'SST']:
            raise ValueError('Unexpected display channel order')
        pixels = tifffile.memmap(image, mode='r')
        low, high = np.asarray(spec['display_min']), np.asarray(spec['display_max'])
        px = spec['pixel_size_um']
        areas = np.asarray([g.area * px**2 for g in geoms])
        status = 'Owner-approved boundaries' if approved else 'Awaiting donor-specific visual review'
        summary = dict(case_id=cid, boundary_status=status, boundaries=len(geoms), raw_cells=donor['cell_count'],
            small=int((areas < 2000).sum()), medium=int(((areas >= 2000) & (areas < 10000)).sum()),
            large=int((areas >= 10000).sum()), section_area_um2=float(areas.sum()),
            model_sha256=manifest['model_sha256'], geometry_sha256=digest(geometry_path),
            minimum_votes=160, forest_trees=160, image=str(image),
            new_cell_membership_available=approved, new_regional_measurements_available=approved)
        if approved and len(geoms) != promotion['accepted_islets']:
            raise ValueError('Approved count differs')
        summaries.append(summary)
        # All crops use fixed quantiles of section area within each size bin.
        chosen = []
        for start, end, label in [(0, 2000, 'Small'), (2000, 10000, 'Medium'), (10000, np.inf, 'Large')]:
            ids = sorted(np.flatnonzero((areas >= start) & (areas < end)), key=lambda i: (areas[i], names[i]))
            for index in sorted(set(int((len(ids) - 1) * q) for q in (.25, .75))) if ids else []:
                chosen.append((ids[index], label))
        fig, axes = plt.subplots(len(chosen), 2, figsize=(10, 3.2 * len(chosen)), squeeze=False)
        for row, (index, label) in enumerate(chosen):
            raw, _, bounds = raw_crop(pixels, spec, geoms[index], pad_um=15)
            x0, y0, x1, y1 = bounds
            rgb = np.clip((raw - low) / np.maximum(high - low, 1), 0, 1)
            for col, ax in enumerate(axes[row]):
                ax.imshow(rgb)
                if col:
                    overlay(ax, geoms, tree, bounds, (x0, y0))
                ax.set_xlim(0, x1 - x0); ax.set_ylim(y1 - y0, 0)
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title(('Raw INS/GCG/SST' if col == 0 else names[index]) + f' — {label}, {areas[index]:,.0f} µm²', fontsize=10)
            selections.append(dict(case_id=cid, name=names[index], size_stratum=label,
                                   area_um2=areas[index], source_pixel_bounds=json.dumps(bounds)))
        fig.suptitle(f'Donor {cid}: {status}\nINS red · GCG green · SST blue · boundaries cyan', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .97))
        save(fig, out, f'{cid}_examples')
        raw = np.stack([pixels[spec['source_channel_names'].index(m), ::8, ::8] for m in spec['rgb_channels']], -1)
        rgb = np.clip((raw - low) / np.maximum(high - low, 1), 0, 1)
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        bounds = (0, 0, pixels.shape[2], pixels.shape[1])
        for col, ax in enumerate(axes):
            ax.imshow(rgb)
            if col:
                overlay(ax, geoms, tree, bounds, (0, 0), scale=8)
            ax.set_axis_off(); ax.set_title('Raw endocrine signal' if not col else f'{len(geoms):,} boundaries')
        fig.suptitle(f'Donor {cid}: {status}')
        fig.tight_layout()
        save(fig, out, f'{cid}_whole_image')
        pages.append(f'<h2 id="donor-{cid}">Donor {cid}: {html.escape(status)}</h2>'
            f'<p>{len(geoms):,} boundaries · {donor["cell_count"]:,} preserved cells</p>'
            f'<a href="{cid}_whole_image.pdf"><img loading="lazy" src="{cid}_whole_image.png" alt="Donor {cid} whole image"></a>'
            f'<details><summary>Size-stratified examples and PDF</summary><a href="{cid}_examples.pdf">'
            f'<img loading="lazy" src="{cid}_examples.png" alt="Donor {cid} size examples"></a></details>')
        del pixels, raw, rgb
        print(f'Rendered donor {cid}: {len(geoms)} boundaries', flush=True)
    pd.DataFrame(summaries).to_csv(out / 'donor_boundary_summary.csv', index=False)
    pd.DataFrame(selections).to_csv(out / 'example_locations.csv', index=False)
    regional_root = rollout / 'regions/78'
    verification = json.loads((regional_root / 'verification.json').read_text())
    regional = pd.read_parquet(regional_root / 'analysis/78.parquet')
    core = regional.query("compartment == 'core' and band_um == 20").copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].scatter(core.area_um2, core.measured_cells, s=7, alpha=.5)
    axes[0].set(xscale='log', yscale='symlog', xlabel='Section area (µm²)', ylabel='Measured core cells')
    axes[0].axhline(10, color='darkorange', ls='--', label='Trajectory minimum'); axes[0].legend()
    classes = ['Nerve', 'Capillary', 'Lymphatic']
    for i, comp in enumerate(['core', 'peri']):
        frame = regional.query('band_um == 20 and compartment == @comp')
        axes[1].bar(np.arange(3) + i*.35, [frame[c.lower()+'_area_fraction'].median() for c in classes], .35, label=comp)
    axes[1].set_xticks(np.arange(3) + .175, classes); axes[1].set_ylabel('Median positive area / valid area'); axes[1].legend()
    fig.suptitle('HPAP-078 rebuilt measurements — descriptive, not phenotype validation')
    fig.tight_layout(); save(fig, out, '78_rebuilt_measurements')
    marker_cols = [c for c in core if c.startswith('median_log1p_')]
    core['section_area_quartile'] = pd.qcut(core.area_um2, 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
    matrix = core.groupby('section_area_quartile', observed=True)[marker_cols].median()
    fig, ax = plt.subplots(figsize=(16, 4))
    im = ax.imshow(matrix.to_numpy(), aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(marker_cols)), [c.removeprefix('median_log1p_') for c in marker_cols], rotation=60, ha='right')
    ax.set_yticks(range(len(matrix)), matrix.index); ax.set_ylabel('Section-area quartile')
    ax.set_title('HPAP-078: all 31 markers; median of per-islet core summaries')
    fig.colorbar(im, ax=ax, label='Median log1p raw intensity'); fig.tight_layout()
    save(fig, out, '78_marker_panel')
    matrix.to_csv(out / '78_marker_panel.csv')
    write_json(out / 'review_manifest.json', dict(status='offline_boundary_review_gallery',
        approved_cases=['78'], awaiting_visual_review=[s['case_id'] for s in summaries if s['case_id'] != '78'],
        qupath_project=str(project / 'project.qpproj'), approved_regional_verification=verification,
        summaries=summaries, source_images_modified=False, application_modified=False))
    table = pd.DataFrame(summaries)[['case_id', 'boundary_status', 'boundaries', 'raw_cells']].to_html(index=False, border=0)
    (out / 'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Endocrine boundary rollout</title>'
        '<style>body{font:16px system-ui;max-width:1450px;margin:2rem auto;padding:1rem}img{max-width:100%}'
        'a{color:#075ac0}td,th{padding:.5rem;text-align:left;border-bottom:1px solid #ddd}summary{cursor:pointer;padding:1rem}</style>'
        '<h1>Endocrine boundary rollout · 18 September 2026</h1>'
        '<p>HPAP-078 preserves the exact 1,431 boundaries approved by the owner. Its native cell parents, '
        'signed distances and regional measurements have been rebuilt. The other five donors use the same '
        'fixed 160/160-vote classifier and IO60 size-adaptive rules and await their own image review.</p>'
        '<p>Open QuPath → File → Project → Open project:</p><pre>' + html.escape(str(project / 'project.qpproj')) + '</pre>'
        '<p><a href="../' + html.escape(project.name) + '/README.md">QuPath instructions</a> · '
        '<a href="donor_boundary_summary.csv">Donor summary CSV</a> · <a href="example_locations.csv">Example coordinates</a> · '
        '<a href="review_manifest.json">Measurement and provenance record</a></p>' + table +
        '<p>Cyan marks current outlines. Grey objects in QuPath preserve previous segmentation. Filled outlines '
        'are islet geometry; hormone-positive tissue burden requires a separate pixel-area measurement. '
        'Boundary approval does not certify cell phenotype gates.</p>' + ''.join(pages) +
        '<h2>HPAP-078 rebuilt measurements</h2><p>Descriptive plots and exports remain outside the application.</p>'
        '<a href="78_rebuilt_measurements.pdf"><img src="78_rebuilt_measurements.png" alt="Rebuilt measurements"></a>'
        '<a href="78_marker_panel.pdf"><img src="78_marker_panel.png" alt="All 31 core marker summaries"></a></html>')


if __name__ == '__main__':
    main()
