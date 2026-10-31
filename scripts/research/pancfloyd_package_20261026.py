#!/usr/bin/env python3
"""Package the 20-donor figure recreation, crosswalk and reproducible inputs."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import csv
import hashlib
import json
import re
import shutil
from xml.etree import ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/pancfloyd_20donor_20261026'
REF=Path('/home/smith6jt/Downloads/PancFloyd_Figures_v3.pptx')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    files=sorted(p for d in ['context','slides1_4','neighborhoods','trajectory','structural_cells'] for p in (OUT/d/'figures').glob('*.png'))
    figures={}
    for p in files:
        m=re.match(r'slide(\d+)_image(\d+)',p.name)
        if m:
            key=tuple(map(int,m.groups()))
            if key in figures: raise ValueError(f'Duplicate figure mapping {key}')
            figures[key]=p
    alias={}
    recovered={(3,1),(3,2),(3,3),(3,4),(4,4),(4,5),(9,1),(9,2)}
    required=[(2,i) for i in range(1,7)]+[(3,i) for i in range(1,9)]+[(4,i) for i in range(3,7)]
    required += [(5,i) for i in range(1,8)]+[(6,i) for i in range(1,9)]+[(7,2)]+[(8,i) for i in range(1,4)]
    required += [(9,i) for i in range(1,5)]+[(10,i) for i in range(1,5)]+[(11,i) for i in range(1,6)]
    assert all(alias.get(k,k) in figures for k in required),'A requested quantitative panel has no current counterpart'
    rows=[]
    with ZipFile(REF) as z:
        for s in range(1,12):
            rel=ET.fromstring(z.read(f'ppt/slides/_rels/slide{s}.xml.rels'))
            for r in rel:
                target=r.get('Target','')
                m=re.search(r'image-\d+-(\d+)\.png$',target)
                if not m: continue
                image=int(m.group(1));key=(s,image);counterpart=alias.get(key,key)
                current=figures.get(counterpart)
                status='recreated';change='Current 20-donor dataset and corrected acquisition-edge eligibility.'
                if key in recovered:
                    status='recovered_structure_cell_rule'
                    change='Existing centroid, majority cell-area overlap, highest-marker and reference-level rule replayed on accepted edge-eligible cells. Correct union-based overlap and fractional passthrough retained. Recovered identities are analysis outputs; served app labels are unchanged.'
                elif key in [(3,7),(3,8)]:
                    status='amended_endpoint';change='SMA-positive area measured with the existing Otsu and minimum-component-area workflow replaces SMA-positive cell fraction. No arteriole identity is assigned.'
                elif key==(5,3):
                    status='amended_endpoint';change='Accepted beta-cell fraction among core endocrine cells replaces continuous INS intensity; no smooth curve is drawn.'
                elif key==(5,1):
                    change='Current saved endocrine-composition UMAP/Leiden0.5, with the same coordinates and edge-eligible display cohort across panels.'
                elif key in [(2,3),(2,4),(8,2)]:
                    change='Same trapezoidal AUC calculation over the shared observed25–225 µm diameter-midpoint range; full-range app integrals retained separately.'
                elif key in [(1,2),(1,3),(1,4),(1,9)]:
                    change='Original physical example field redrawn using current donor6533 annotations; current object name Islet_79 replaces historical Islet_199.'
                elif key==(1,5):
                    status='updated_native_labels';change='Cohort labels updated to ND7, Aab+7, T1D6.'
                elif key==(7,1):
                    status='updated_native_diagram';change='App architecture and cohort/data counts updated as editable PowerPoint text and shapes.'
                elif current is None:
                    status='retained_source_illustration';change='Source schematic, microscopy image, arrow or scale-bar artwork retained unchanged. No new quantitative conclusion or donor identity assigned.'
                rows.append(dict(source_slide=s,source_image=image,source_media=target,
                    status=status,current_panel=f'slide{counterpart[0]:02d}_image{counterpart[1]:02d}' if current else '',
                    current_figure=str(current.relative_to(OUT)) if current else '',change=change))
    with (OUT/'panel_crosswalk.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    legend_files=sorted(p for part in ['context','slides1_4','neighborhoods','trajectory','structural_cells'] for p in (OUT/part).rglob('*')
                        if p.is_file() and (p.name.endswith('_legend.txt') or p.name.endswith('.legend.txt') or p.name in ['legends.md','figure_legends.md']))
    legends=['# Figure legends','', 'All current group comparisons use ND, Aab+ and T1D. The donor cohort is 7, 7, 6. Exact panel populations and measurements are given below. These legends report methods and numerical tests without biological interpretation.','']
    for p in legend_files:
        legends+=['## '+str(p.relative_to(OUT)), '',p.read_text().strip(),'']
    (OUT/'Figure_legends.md').write_text('\n'.join(legends))
    (OUT/'README.md').write_text('''# PancFloyd figures: current 20-donor dataset

The PowerPoint and PDF recreate the supplied figure deck with the accepted uniform-20donor-20261022-v10 data: 7 ND, 7 Aab+, 6 T1D donors. The original source presentation is unchanged. Donor 78 remains excluded. All numerical source panels have a recreated panel or an explicitly recorded current counterpart in `panel_crosswalk.csv`.

## Files

- `PancFloyd_Figures_20donors.pptx`: revised 12-slide presentation, preserving the 11 figure slides and source credit slide. Scientific plots remain picture objects, as in the supplied deck. Editable plotting data and vector figures are supplied separately.
- `PancFloyd_Figures_20donors.pdf`: presentation preview.
- `Figure_legends.md`: separate manuscript-style methods legends and exact numerical test details.
- `prism/`: native GraphPad Prism data file and CSV/TSV tables. These contain data, not embedded Prism analyses or graph layouts.
- `context/`, `slides1_4/`, `neighborhoods/`, `trajectory/`, `structural_cells/`: individual PNG/PDF/SVG figures, exact object/donor/bin values, eligibility records and provenance.
- `source/`: scripts used for plotting, area correction, deck assembly and export.

## Measurement updates

Vascular, neural and lymphatic cell panels use the existing structure-based rule: centroid within the object, majority cell-area overlap, highest expressed defining marker and the existing marker-reference check. Senior's union-based overlap and fractional passthrough definitions are retained. The current CODEX run computes the rule but disables retention of these cell identities before serving the app. This package recovers them offline on accepted edge-eligible cells, keeps original app identities alongside the recovered calls, and leaves served data unchanged. Cell fractions/densities and structure-area curves remain separate panels.

The SMA panels measure SMA-positive area using the previously established image-threshold and component-filter method, because the three-structure cell rule contains no SMA/arteriole class. They do not label that area as arterioles. The quantitative INS trajectory panel is replaced by accepted beta-cell fraction. Other endocrine fractions use accepted categorical cell labels.

Structure areas use the union of overlapping same-class annotations before intersection with core/peri regions. This corrects double counting in the saved summed-area summaries without changing annotations or rebuilding the served pipeline. Main endocrine AUC bars use the same 25–225 µm midpoint range across groups. Original full-available-range integrals remain as supporting data. Old significance stars are removed and all displayed stars are generated from current tests, defined in the legends.

The current composition embedding and DPT are reused. DPT fits all 1,215 eligible app objects; the reference trajectory panels show the 1,126 objects with complete cores and 20 µm neighborhoods, retaining all 20 donors, without reranking or refitting after display filtering. Cell counts, morphology summaries and distance profiles have their own appropriate populations, listed in their legends. Donor-level tests retain donors as biological replicates. Object-level curves and pooled distance profiles retain their explicitly stated descriptive weighting.

## Hypothesis work

The parallel app-based H1–H10 revision is saved in `docs/app_based_hypothesis_plan_20261026.md` in the repository and copied here to `reports/`. It uses the same app composition states, pseudotime and physical neighborhoods. This figure recreation does not assign biological conclusions or original-plan if–then branches and does not claim to have run amended H5/H7 hypothesis models.
''')
    (OUT/'reports').mkdir(exist_ok=True)
    shutil.copy2(ROOT/'docs/app_based_hypothesis_plan_20261026.md',OUT/'reports/app_based_hypothesis_plan_20261026.md')
    (OUT/'source').mkdir(exist_ok=True)
    for p in (ROOT/'scripts/research').glob('pancfloyd_*_20261026.*'):
        shutil.copy2(p,OUT/'source'/p.name)
    for name in ['build.mjs','finalize.mjs','finalize_verified.mjs']:
        shutil.copy2(OUT/'.build'/name,OUT/'source'/name)
    sources={'reference':{'path':str(REF),'sha256':sha(REF)},'plot_images':len(files),
             'required_source_quantitative_slots':len(required),'mapped_quantitative_slots':len(required),
             'group_donors':{'ND':7,'Aab+':7,'T1D':6},'interpretation_added':False,
             'served_pipeline_rebuilt':False,'panel_crosswalk_rows':len(rows)}
    (OUT/'delivery_manifest.json').write_text(json.dumps(sources,indent=2))
    deliverable_files=sorted(p for p in OUT.rglob('*') if p.is_file() and '.build' not in p.parts and p.name!='SHA256SUMS.txt')
    (OUT/'SHA256SUMS.txt').write_text('\n'.join(f'{sha(p)}  {p.relative_to(OUT)}' for p in deliverable_files)+'\n')
    archive=ROOT/'output/PancFloyd_Figures_20donors.zip'
    with ZipFile(archive,'w',ZIP_DEFLATED,compresslevel=6) as z:
        for p in [*deliverable_files,OUT/'SHA256SUMS.txt']:
            z.write(p,arcname=Path('PancFloyd_Figures_20donors')/p.relative_to(OUT))
    with ZipFile(archive) as z: assert z.testzip() is None
    print(json.dumps({**sources,'archive':str(archive),'archive_bytes':archive.stat().st_size},indent=2))
if __name__=='__main__': main()
