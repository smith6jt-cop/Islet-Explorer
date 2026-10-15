#!/usr/bin/env python3
"""Package existing hypothesis results as offline tables, plots and an HTML index."""
import argparse
import html
import json
import re
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from common import sha256, write_json


def effect_plot(row, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 2.8))
    if np.isfinite(row.ci_low) and np.isfinite(row.ci_high):
        ax.plot([row.ci_low, row.ci_high], [0, 0], color='#256b7d', linewidth=3)
    ax.scatter([row.effect], [0], color='#173f4a', zorder=3)
    ax.axvline(0, color='gray', linestyle='--', linewidth=.8)
    ax.set(title=row.id+' — effect and reported confidence interval', xlabel=row.unit, yticks=[])
    ax.text(.02, .05, f'n={row.n_donors} donors; {row.n_islets} islets; BH q={row.q_value:.4g}', transform=ax.transAxes, fontsize=9)
    fig.tight_layout()
    for suffix in ['png', 'pdf']: fig.savefig(output.with_suffix('.'+suffix), dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--source-report-root', required=True, type=Path)
    p.add_argument('--output-root', required=True, type=Path)
    p.add_argument('--identity-review', required=True, type=Path)
    a = p.parse_args(); root = a.input_root.resolve(strict=True); source = a.source_report_root.resolve(strict=True); out = a.output_root.resolve()
    if out.exists() or out.is_relative_to(root): raise ValueError('Use a fresh report folder outside the app runtime')
    out.mkdir(parents=True)
    shutil.copytree(source, out/'source_tables')
    results = pd.read_csv(source/'results.csv')
    if set(results.id) != {f'H{i}' for i in range(1, 11)}: raise ValueError('Expected the ten registered hypotheses')
    if results.loc[results.status == 'unavailable', 'p_value'].notna().any(): raise ValueError('Unavailable tests cannot have p-values')
    shutil.copy2(source/'results.csv', out/'results.csv')
    shutil.copy2(source/'registry.json', out/'registry.json')
    inclusion = root/'evidence/hypothesis_donor_inclusion.csv'
    if inclusion.is_file(): shutil.copy2(inclusion, out/'donor_inclusion.csv')
    shutil.copy2(a.identity_review/'donor_summary.csv', out/'islet_identity_review.csv')
    cohort = pd.read_csv(root/'source_manifest.csv', dtype=str).set_index('case_id')
    sections = []
    warning = ('These are the previously computed estimates, pending scientific review. '
               'Islet-identity review has identified suspect annotations, including HPAP-078 regions without hormone-positive cells. '
               'Correcting confirmed region errors requires rebuilding dependent analyses. No tests were refitted for this report export.')
    for row in results.itertuples(index=False):
        folder = out/row.id; folder.mkdir()
        for path in source.iterdir():
            if path.is_file() and re.match(r'^'+row.id+r'(?:[.-]|$)', path.name): shutil.copy2(path, folder/path.name)
        results.loc[results.id == row.id].to_csv(folder/'result.csv', index=False)
        body = f'<h2>{html.escape(row.id+": "+row.title)}</h2><p>{html.escape(str(row.status))}</p>'
        if np.isfinite(row.effect):
            effect_plot(row, folder/'effect_interval')
            body += f'<img src="{row.id}/effect_interval.png" alt="{row.id} effect and confidence interval"><p><a href="{row.id}/effect_interval.pdf">Vector PDF</a></p>'
        else:
            body += '<p>'+html.escape(str(row.reason))+'</p>'
        representatives = folder/f'{row.id}-representative-islets.csv'
        if representatives.exists():
            examples = pd.read_csv(representatives, dtype={'case_id': str})
            examples['source_image'] = examples.case_id.map(cohort.image)
            examples.to_csv(folder/'representative_image_references.csv', index=False)
        links = [f'<li><a href="{row.id}/{html.escape(p.name)}">{html.escape(p.name)}</a></li>' for p in sorted(folder.iterdir()) if p.is_file()]
        sections.append(body+'<ul>'+''.join(links)+'</ul>')
        (folder/'README.md').write_text(f'# {row.id}: {row.title}\n\nStatus: {row.status}\n\n'+warning+'\n\n'+(str(row.reason) if pd.notna(row.reason) else '')+'\n')
    followups = out/'followups'; followups.mkdir()
    for path in source.glob('*followup*'): shutil.copy2(path, followups/path.name)
    content = ('<!doctype html><html lang="en"><meta charset="utf-8"><title>All-donor research reports</title>'
        '<style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem}img{max-width:100%}td,th{padding:.4rem;border-bottom:1px solid #ddd}table{border-collapse:collapse}aside{background:#fff0d5;padding:1rem}</style>'
        '<h1>All-donor research reports</h1><aside>'+html.escape(warning)+'</aside>'
        '<p>All 21 donors are in the dataset. The primary contrast is T1D versus ND (14 donors). '
        'Pending annotation review excludes ND 6539, leaving 13 donors; H9 additionally excludes T1D 6520 because no islets meet its endocrine-cell sample requirement. '
        'Aab-positive contrasts are secondary and have not been run.</p>'
        '<p><a href="results.csv">Combined results CSV</a> · <a href="registry.json">Hypothesis registry</a> · '
        '<a href="donor_inclusion.csv">Donor inclusion audit</a> · <a href="islet_identity_review.csv">Islet-identity review counts</a></p>'
        +results[['id','status','effect','ci_low','ci_high','p_value','q_value','n_donors','n_islets']].to_html(index=False, na_rep='Unavailable', escape=True)
        +''.join(sections)+'<h2>Exploratory follow-ups</h2><p><a href="followups/followups.csv">Follow-up results and availability</a></p>'
        '<p>Complete source tables, fold diagnostics, permutations and sensitivity results are in source_tables/.</p></html>')
    (out/'index.html').write_text(content)
    (out/'README.md').write_text('# Offline research reports\n\nOpen index.html to browse the results, plots and downloads.\n\n'+warning+'\n\nThe application does not serve or require this folder. The source tables and estimates were copied without refitting.\n')
    files = [dict(path=str(p.relative_to(out)), bytes=p.stat().st_size, sha256=sha256(p)) for p in sorted(out.rglob('*')) if p.is_file()]
    write_json(out/'analysis_manifest.json', dict(version=1, source_report_root=str(source), source_cohort_sha256=sha256(root/'cohort.csv'),
        source_results_sha256=sha256(source/'results.csv'), results_refitted=False, scientific_acceptance=False,
        islet_identity_review=str(a.identity_review.resolve()), primary_tests=10, files=files,
        runtime_dependency=False, total_bytes=sum(f['bytes'] for f in files)))
    print(out/'index.html', flush=True)


if __name__ == '__main__': main()
