#!/usr/bin/env python3
"""Maintain factual H5 narrative/legend text without fitting or testing models."""
from pathlib import Path
import hashlib,json,re

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/h5_h10_original_plan_20261026/h5'


def refresh_narratives():
    path=OUT/'README.md';text=path.read_text();parts=text.split('\n\n')
    parts[1]=('The donor-level conditional disease comparison has Kruskal–Wallis p=0.914550; '
        'all three pairwise BH-adjusted p values are 1. Measurements, cohort counts, model estimability, '
        'and technical diagnostics are recorded below. Branch decisions are unassigned (user judgment).')
    text='\n\n'.join(parts)
    text=text.replace('Median-INS, raw-INS, RedSEA-corrected-INS, and stronger hormone-margin variants also give no significant disease-by-penetration interaction.',
        'Median-INS, raw-INS, RedSEA-corrected-INS, and stronger hormone-margin variant statistics are recorded in model_anova.csv and disease_slope_comparisons.csv.')
    text=text.replace('A brief qualification is that retained β-cell signal is conditional on the established INS-associated phenotype and sparse T1D β-cell retention; a follow-up using an INS-independent β-cell marker would strengthen inference about cell identity versus hormone state. This limitation does not prevent reporting the present computational result.',
        'The retained β-cell INS endpoint uses the current typed Beta-cell population; INS is one of the markers used by the established typing rules.')
    note=('## Original plan conditions\n\nThe original H5 if–then conditions are preserved verbatim in '
        '`original_if_then_conditions.json`. The branch decision in `if_then_decision.json` is unassigned (user judgment).')
    if '## Original plan conditions' not in text:text=text.rstrip()+'\n\n'+note+'\n'
    path.write_text(text)
    path=OUT/'insulin_quantity_and_state_legend.txt';text=path.read_text()
    text=text.replace(' The lower mean INS in T1D does not reach significance after the three pairwise comparisons are adjusted.','')
    path.write_text(text)
    path=OUT/'immune_localization_beta_insulin_legend.txt';text=path.read_text()
    text=text.replace('Nineteen donors support a full-rank regression','Nineteen donors have a full-rank regression')
    text=text.replace('It supports 18 donors','It has 18 estimable donors')
    text=text.replace(' The results provide no evidence that this conditional association differs among disease states; they do not establish identical relationships or a negative within-β-cell insulin response.','')
    path.write_text(text)
    plan_path=ROOT/'docs/T1D_research_hypotheses_20261022.md'
    section=plan_path.read_text().split('### H5 ',1)[1].split('\n### ',1)[0]
    conditions=[dict(condition=title,text=value.strip()) for title,value in re.findall(r'\*\*(If[^*]+)\*\*\s*(.*?)(?=\n\n|\Z)',section,flags=re.S)]
    (OUT/'original_if_then_conditions.json').write_text(json.dumps(dict(source=str(plan_path),hypothesis='H5',original_conditions=conditions),indent=2)+'\n')
    (OUT/'if_then_decision.json').write_text(json.dumps(dict(hypothesis='H5',branch_decision='unassigned (user judgment)'),indent=2)+'\n')
    manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json']
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Updated H5 factual narrative/legends and preserved original conditions; no tests or models run.')


if __name__=='__main__':refresh_narratives()
