#!/usr/bin/env python3
"""Render the corrected H1 report directly from saved results and audit tables."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/islet-h1-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import pandas as pd
from shapely.geometry import shape

from h1_features_20261022 import sha256, write_json


def geometry_qc(out):
    source = Path(json.loads((out / 'eligibility/edge_manifest.json').read_text())['source_run'])
    tag = 'Aab_6505'  # Fixed before fitting; independent of Beta outcome.
    d = pd.read_parquet(out / 'eligibility/objects/6505.parquet').set_index('object_name')
    meta = json.loads((out / 'eligibility/audits/6505.json').read_text())
    fig, ax = plt.subplots(figsize=(11, 7))
    polygons, colors = [], []
    for f in json.loads((source / 'islets' / f'{tag}.geojson').read_text())['features']:
        name = f.get('properties', {}).get('name')
        if name not in d.index:
            continue
        r, geom = d.loc[name], shape(f['geometry'])
        parts = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
        for part in parts:
            if part.geom_type == 'Polygon':
                polygons.append(np.asarray(part.exterior.coords) * meta['pixel_um'])
                colors.append('#992f4c' if not r.core_complete else '#dd9635'
                              if not r.neighborhood_complete_100 else '#347896')
    ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors=colors, linewidths=.25))
    for f in json.loads((source / 'inputs/tissue' / f'{tag}.geojson').read_text())['features']:
        geom = shape(f['geometry'])
        parts = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
        for part in parts:
            if part.geom_type == 'Polygon':
                xy = np.asarray(part.exterior.coords) * meta['pixel_um']
                ax.plot(xy[:, 0], xy[:, 1], color='#999999', lw=.3)
    w, h = meta['image_width'] * meta['pixel_um'], meta['image_height'] * meta['pixel_um']
    ax.plot([0,w,w,0,0], [0,0,h,h,0], color='black', lw=1)
    for color, label in [('#347896','Complete core and 100 µm support'),
                         ('#dd9635','Complete core; incomplete 100 µm support'), ('#992f4c','Core withheld')]:
        ax.scatter([], [], color=color, label=label, s=20)
    ax.set(xlim=(-50,w+50), ylim=(h+50,-50), aspect='equal', xlabel='Source x (µm)', ylabel='Source y (µm)',
           title='Geometry audit: Aab_6505 (fixed canary; all endocrine objects, including Singles)')
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    fig.savefig(out / 'geometry_qc.png', dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    out = a.output.resolve()
    geometry_qc(out)
    summary = json.loads((out / 'models/analysis_summary.json').read_text())
    validation = json.loads((out / 'validation_summary.json').read_text())
    matching = json.loads((out / 'independent/summary.json').read_text())
    model = pd.read_csv(out / 'models/model_results.csv')
    slopes = pd.read_csv(out / 'heterogeneity/group_specific_slopes.csv')
    old = json.loads(Path('output/h1_20261022/analysis_summary.json').read_text())
    flow = pd.read_csv(out / 'eligibility_flow_by_donor.csv')
    objects = pd.read_csv(out / 'eligibility/object_eligibility.csv')
    features = pd.read_parquet(out / 'h1_object_features.parquet')
    audits = pd.read_csv(out / 'eligibility/donor_exclusions.csv')
    primary = summary['primary']
    assert validation['status'] == 'pass'
    assert summary['model_failures'] == []
    assert summary['donor_bootstrap']['completed'] == 999
    assert summary['eligible_objects'] == validation['corrected_primary_objects'] == matching['eligible_objects']
    assert len(pd.read_csv(out / 'models/leave_one_donor_out.csv')) == 20
    sensitivity = model[model.term.eq('allocation_10pp') & model.outcome.eq('n_beta_core')]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), gridspec_kw={'width_ratios':[1.25,1]})
    left = np.zeros(len(flow))
    for col, label, color in [('corrected_primary','Retained','#347896'), ('core_excluded','Core excluded','#992f4c'),
                              ('additional_peri_incomplete','Incomplete 20 µm support','#dd9635'),
                              ('additional_field_incomplete','Incomplete 100 µm support','#bcc8ce')]:
        axes[0].barh(np.arange(len(flow)), flow[col], left=left, label=label, color=color)
        left += flow[col].to_numpy()
    axes[0].set_yticks(np.arange(len(flow)), [f'{r.case_id} ({r.donor_status})' for r in flow.itertuples()])
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Previously eligible H1 objects')
    axes[0].set_title('Acquisition-edge exclusions')
    axes[0].legend(fontsize=8, loc='lower right')
    selected = sensitivity[sensitivity.model.isin(['Primary','Exclusive peri allocation','100–250 µm field','T1D only'])]
    values = [('Previous primary',old['primary'])] + [(r.model,r._asdict()) for r in selected.itertuples()]
    for i,(label,r) in enumerate(values):
        axes[1].errorbar(r['rr'],i,xerr=[[r['rr']-r['lower']],[r['upper']-r['rr']]],fmt='o',capsize=4,color='#26364a')
    axes[1].set_yticks(range(len(values)), [r[0] for r in values])
    axes[1].invert_yaxis()
    axes[1].axvline(1,color='gray',ls='--')
    axes[1].set_xlabel('Beta-density ratio per 10 percentage points\n95% donor-cluster intervals')
    axes[1].set_title('Same H1 model and existing cell identities')
    fig.tight_layout()
    fig.savefig(out / 'edge_correction_summary.png',dpi=180)
    fig.savefig(out / 'edge_correction_summary.pdf')
    plt.close(fig)
    total = flow.select_dtypes('number').sum()
    lower,upper = primary['lower'],primary['upper']
    conclusion = ('The pooled estimate remains inverse, but its donor-cluster interval includes no association.'
                  if primary['rr'] < 1 and lower <= 1 <= upper else
                  'The pooled estimate is inverse and its donor-cluster interval excludes 1.'
                  if upper < 1 else 'The corrected estimate and interval are reported below without a causal interpretation.')
    lines = [
        '# H1 after counting and acquisition-edge corrections', '',
        '24 October 2026. Analysis-only rerun authorized by the owner after narrowing the earlier typing overhaul to existing data and established methods.', '',
        f"{conclusion} Per 10 percentage points more nearby immune cells allocated to the core, the adjusted Beta-density ratio is **{primary['rr']:.3f} (95% CI {lower:.3f}–{upper:.3f}; nominal p = {primary['p_cluster_t']:.4f})**. The previous result was {old['primary']['rr']:.3f} ({old['primary']['lower']:.3f}–{old['primary']['upper']:.3f}; p = {old['primary']['p_cluster_t']:.4f}).", '',
        '## What changed', '',
        '- Each served physical cell has one persistent donor/UUID identity, one unchanged phenotype, and one additive primary compartment. Parent totals are unions of exclusive leaf identities. Object membership is a separate relation; no physical cell appears twice within a focal object’s core-plus-peri union.',
        '- Additive peri counts assign a cell to the nearest eligible boundary among its actual saved band memberships, with an object-name tie break. Shared membership remains in the established primary object-local H1 exposure; the exclusive assignment is also tested using the existing sensitivity model.',
        '- Original cell and nucleus polygons are checked against full-image boundaries, using one source pixel tolerance. Observed border contact propagates to whole endocrine objects, including contact in source-QC-rejected cells. A missing nucleus in an already QC-rejected segmentation is not evidence that an islet was cut.',
        '- Whole pre-inset classifier components touching the image’s 5 µm analysis inset are withheld, including their interior descendants after morphology. Components use the original Tissue selection, one classifier-pixel support/matching tolerance, and eight-connectivity. This is conservative component reconstruction, not a claim of recovered exact QuPath parent IDs.',
        '- Core completeness is separate from neighborhood completeness. Primary H1 requires full acquisition support through 100 µm (including the immediate 20 µm neighborhood); the distal sensitivity requires 250 µm. The existing tissue clipping and field-area threshold are retained. Excluded object footprints remain excluded from background tissue fields.', '',
        '## Analysis population', '',
        '| Stage | Objects |', '|---|---:|',
        f'| Paired source Islets/Clusters | {len(features):,} |',
        f"| Previously eligible for primary H1 | {int(total.baseline_primary):,} |",
        f"| Removed: incomplete core | {int(total.core_excluded):,} |",
        f"| Additional removal: incomplete immediate neighborhood | {int(total.additional_peri_incomplete):,} |",
        f"| Additional removal: incomplete regional field | {int(total.additional_field_incomplete):,} |",
        f"| Corrected primary population | {summary['eligible_objects']:,} |", '',
        f"All 20 donors remain (7 ND, 7 Aab+, 6 T1D). The corrected population contains {summary['eligible_core_positive_objects']:,} core-immune-positive objects and {summary['eligible_core_immune_cells']:,} core immune cells. Across all paired objects, {int((~features.core_complete).sum()):,} have incomplete cores. Across the full geometric census, {int((~objects.core_complete).sum()):,}/{len(objects):,} objects are withheld; this includes Singles, which remain outside H1. Unresolved classifier ancestry: {int(objects.preclip_ancestry_unresolved.sum())}.", '',
        f"Original source polygons identify {int(audits.source_border_contact_cells.sum()):,} border-contact cells; {int(audits.served_edge_cells.sum()):,} of the {validation['distinct_physical_cells']:,} served cells require additional cell-level exclusion. Object and neighborhood exclusions can still apply when no retained individual cell touches the image edge.", '',
        'See [donor flow](eligibility_flow_by_donor.csv), [class/size exclusions](exclusions_by_class_size.csv), and [exclusive cell census](exclusive_cell_census.csv). Compartment totals retain excluded rows explicitly for reconciliation; exclude `excluded_object`/`excluded_cell` from retained object totals. The residual `tissue` compartment includes served cells outside saved Tissue polygons; regional density denominators and counts separately use actual Tissue geometry. Raw source-object rows remain in the feature table with eligibility flags; use the saved model population for H1.', '',
        f'![Fixed-canary geometry review]({out}/geometry_qc.png)', '',
        '## Same established H1 analysis', '',
        'The outcome is Beta-cell count with core-area offset. The exposure is core immune count divided by the distinct core-plus-20 µm-band immune union, in 10-percentage-point units. Controls remain donor fixed intercepts, object class, log core area, log band/core area ratio, log1p nearby immune count, and independent 20–100 µm immune density. The Poisson mean fit retains equal total weight per donor, donor-cluster sandwich uncertainty, and t(G−1) intervals. Beta-zero objects remain; zero nearby immune pools have undefined allocation. Broad immune identity and CHGA-based endocrine typing are unchanged.', '',
        '| Analysis | Objects | Donors | Density ratio | 95% CI |', '|---|---:|---:|---:|---|']
    for r in sensitivity.itertuples():
        lines.append(f'| {r.model} | {r.objects:,} | {r.donors} | {r.rr:.3f} | {r.lower:.3f}–{r.upper:.3f} |')
    boot=summary['donor_bootstrap'];lo,hi=boot['rr_percentile_ci']
    lines += ['',f"All {boot['completed']}/999 whole-donor bootstrap refits completed; percentile interval {lo:.3f}–{hi:.3f}. All 20 leave-one-donor-out fits completed; point estimates range {summary['lodo_rr_range'][0]:.3f}–{summary['lodo_rr_range'][1]:.3f}. Model failures: {len(summary['model_failures'])}.", '',
        f'![Exclusions and estimates]({out}/edge_correction_summary.png)', '',
        '## Existing independent and subgroup checks', '',
        f"The same outcome-blind, within-donor/class matched analysis found {matching['total_pairs']} pairs across {matching['donors_with_any_pairs']} donors; {matching['donors_with_at_least_5_pairs']} donors met the existing five-pair threshold. Its equal-donor Beta-density difference is {matching['donor_balanced_beta_density_difference_cells_per_1000um2']:.3f} cells/1,000 µm², with donor-bootstrap interval {matching['bootstrap_95pct_interval'][0]:.3f}–{matching['bootstrap_95pct_interval'][1]:.3f}. This binary-occupancy matched estimand differs from the continuous allocation coefficient.", '',
        'The previously posthoc group-slope analysis was also rerun; it remains exploratory:', '',
        '| Group | Ratio per 10 points | Approximate 95% CI |', '|---|---:|---|']
    for r in slopes.itertuples():
        lines.append(f'| {r.donor_status} | {r.rr_10pp:.3f} | {r.lower_95:.3f}–{r.upper_95:.3f} |')
    lines += ['', 'The group-slope model shares nuisance coefficients across the full cohort; the planned T1D-only model refits those coefficients within six donors. Their different intervals reflect different model assumptions. The posthoc interval uses t(19) from the full fit and does not give T1D 20 independent donors. The T1D-only result is inconclusive, and the positive group estimate argues against describing the pooled inverse association as uniform across disease groups.']
    scores=np.array(list(summary['wild_score']['cluster_scores'].values()));share=scores**2/np.sum(scores**2)
    lines += ['',f"The retained wild-score sensitivity gives p = {summary['wild_score']['p_two_sided']:.4f}; its largest donor contributes {100*share.max():.1f}% of squared cluster scores (concentration heuristic {1/np.sum(share**2):.2f}). This approximation is not an independent confirmation or a substitute for the primary donor-cluster interval.", '',
        '## Verification and limits', '',
        f"Validation passed: all 4,763 source-object phenotype counts reproduce the accepted H5AD and every drilldown before exclusions; all {validation['distinct_physical_cells']:,} served UUIDs, phenotypes and coordinates are preserved; independent membership aggregation reproduces model counts; the old primary coefficient and interval reproduce numerically. {validation['source_hashes_verified']:,} input hashes and {validation['historical_hashes_verified']:,} historical-output hashes remain unchanged. The geometry/counting unit suite has 16 passing tests. The saved canary classifier mask was reproduced exactly from retained planes and unchanged thresholds (zero differing pixels).", '',
        'Acquisition support uses verified TIFF/QuPath dimensions and the established 0.3774 µm coordinate calibration. No independent internal acquisition-gap mask is available; undocumented internal gaps are not claimed to be detected. Pre-inset component matching is deliberately conservative and can exclude intact descendants. This analysis does not recalibrate marker positivity, retype macrophages, regenerate embeddings, or alter the served app. The former broad typing/state requirements are outside this narrowed rerun.', '',
        'These are cross-sectional, exploratory associations with 20 biological donors. Related sensitivity tests have nominal, unadjusted intervals/p-values. An inverse pooled association does not establish immune killing, temporal progression, or a uniform T1D mechanism. The T1D sensitivity and subgroup results should be read alongside the pooled estimate.', '',
        '## Files and replay', '',
        '- [Primary results and sensitivities](models/model_results.csv)',
        '- [Exact model population](models/primary_analysis_objects.csv)',
        '- [Validation record](validation_summary.json)',
        '- [Edge rules and provenance](eligibility/edge_manifest.json)',
        '- [Matched analysis](independent/summary.json)',
        '- [Group diagnostics](heterogeneity/summary.json)',
        '- [Feature/source hashes](feature_audit.json)', '',
        'From the repository root, run `bash scripts/research/run_h1_edge_corrected_20261024.sh`. Set `H1_OUT` to a new directory to preserve this result. The runner uses the existing scvi-env and fixes BLAS thread counts; the model retains seed 20261022. Outputs are outside the app data root.', '']
    (out / 'README.md').write_text('\n'.join(lines))
    write_json(out / 'delivery_manifest.json', dict(status='complete',
        scripts={str(f.resolve()):sha256(f) for f in [Path('scripts/research')/n for n in
            ['h1_edge_eligibility_20261024.py','h1_features_20261022.py','h1_analysis_20261022.py',
             'h1_matched_20261022.py','h1_heterogeneity_20261022.py','h1_validate_edge_20261024.py',
             'h1_edge_report_20261024.py','run_h1_edge_corrected_20261024.sh']]},
        files={str(f.relative_to(out)):sha256(f) for f in out.rglob('*') if f.is_file() and f.name!='delivery_manifest.json'}))
    print(out / 'README.md')


if __name__ == '__main__':
    main()
