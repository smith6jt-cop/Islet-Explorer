#!/usr/bin/env python3
"""Assemble measured outputs and manuscript legends without interpretation."""
from pathlib import Path
import json,re,hashlib,shutil,html,io,ast
from xml.sax.saxutils import escape
import numpy as np
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.pagesizes import A4,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PILImage
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/app_hypotheses_20261026'
PARTS=['shared','h1_h4','h5','h6_h8','h9_h10']
TITLE={'H1':'Immune localization and Beta-cell density','H2':'ColIV coverage and immune localization','H3':'Vessel-crossing immune localization','H4':'LGALS3 source allocation and Beta density','H5':'Beta-cell presence (amended)','H6':'Excess endocrine adjacency and Beta density','H7':'Nerve-vessel proximity and Beta density (amended)','H8':'Immune field and small Beta-bearing Clusters','H9':'Coverage association within shared Beta range','H10':'Fixed spatial score and C-peptide'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fmt(x):
    if pd.isna(x):return 'Unavailable'
    return f'{float(x):.5g}'

def figure_title(r):
    stem=r['stem'];match=re.match(r'^(h\d+)_(.*)',stem)
    if match:
        subject=TITLE[match.group(1).upper()].replace(' (amended)','')
        detail=match.group(2).replace('states_states','state profiles').replace('_',' ')
        return f'{subject} · {detail}'
    special={
        'app_representation_overview':'Current app composition map, clusters and pseudotime',
        'direct_disease_comparisons':'Disease comparisons of adjusted associations',
        'conditional_associations':'Conditional associations by disease state',
        'crossing_compartment_comparison':'Immune density near vessel crossings and in reference tissue',
        'immune_neighborhood_profiles':'Core and peri immune neighborhoods',
        'cell_quantity':'Beta-cell presence and endocrine-cell density',
        'cluster_occupancy':'Composition-cluster occupancy',
        'endocrine_states_states':'Endocrine-cell measurements across composition states and pseudotime',
        'endocrine_states_clusters':'Endocrine-cell measurements by composition cluster',
    }
    return special.get(stem,stem.replace('_',' ').capitalize())

def primary():
    sources=[OUT/p/'primary_comparisons.csv' for p in ['h1_h4','h5','h6_h8','h9_h10']]
    sources=[p if p.exists() else p.with_name('primary_tests.csv') for p in sources]
    for p in sources:assert p.exists(),p
    t=pd.concat([pd.read_csv(p) for p in sources],ignore_index=True)
    assert len(t)==10 and t.hypothesis.nunique()==10
    t['order']=t.hypothesis.str.extract(r'(\d+)').astype(int);t=t.sort_values('order').drop(columns='order')
    valid=t.p_value.notna();ix=t.loc[valid,'p_value'].sort_values().index;p=t.loc[ix,'p_value'].to_numpy()
    q=np.minimum.accumulate((p*10/np.arange(1,len(p)+1))[::-1])[::-1]
    t['p_bh_ten_slots']=np.nan;t.loc[ix,'p_bh_ten_slots']=np.minimum(1,q)
    t['family']='Current direct disease comparisons H1-H8; shared-range H9; nested clinical H10. Ten slots retained.'
    t.to_csv(OUT/'primary_statistics.csv',index=False)
    return t,sources

def figures():
    records=[]
    for part in PARTS:
        for p in sorted((OUT/part).rglob('figures/*.png')):
            leg=p.with_name(p.stem+'_legend.txt')
            if not leg.exists():leg=p.with_suffix('.legend.txt')
            if not leg.exists():raise ValueError(f'Missing legend: {p}')
            assert p.with_suffix('.pdf').exists() and p.with_suffix('.svg').exists(),p
            records.append({'part':part,'stem':p.stem,'png':str(p.relative_to(OUT)),'pdf':str(p.with_suffix('.pdf').relative_to(OUT)),'svg':str(p.with_suffix('.svg').relative_to(OUT)),'legend':str(leg.relative_to(OUT))})
    # Put main measured comparisons before their representation companions.
    def order(r):
        priorities={'app_representation_overview':-10,'direct_disease_comparisons':-9,'conditional_associations':-8,'cell_quantity':-9,'cluster_occupancy':-8,'coverage_association_and_beta_quantity':-9,'clinical_held_out_prediction':-8}
        return PARTS.index(r['part']),priorities.get(r['stem'],0),r['png']
    records=sorted(records,key=order)
    for i,r in enumerate(records,1):r['figure']=i
    pd.DataFrame(records).to_csv(OUT/'figure_index.csv',index=False)
    return records

def paragraph(text,style):
    return Paragraph(escape(str(text)).replace('\n','<br/>'),style)

def build_pdf(t,records):
    pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBold','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
    styles=getSampleStyleSheet()
    for s in styles.byName.values():s.fontName='DejaVu'
    styles.add(ParagraphStyle('Main',fontName='DejaVuBold',fontSize=20,leading=25,spaceAfter=12))
    styles.add(ParagraphStyle('BodySmall',fontName='DejaVu',fontSize=9,leading=12,spaceAfter=7))
    styles.add(ParagraphStyle('Caption',fontName='DejaVu',fontSize=8,leading=10.5,spaceAfter=6))
    styles.add(ParagraphStyle('Cell',fontName='DejaVu',fontSize=7.7,leading=10))
    styles.add(ParagraphStyle('Section',fontName='DejaVuBold',fontSize=13,leading=17,spaceAfter=10))
    width,height=landscape(A4);mw=width-64;mh=height-72
    path=OUT/'App_Hypotheses_Figures_Statistics.pdf'
    doc=SimpleDocTemplate(str(path),pagesize=(width,height),rightMargin=32,leftMargin=32,topMargin=32,bottomMargin=40,title='App-based hypotheses: figures and statistics',author='Islet Explorer analysis')
    story=[paragraph('App-based hypotheses: figures and statistics',styles['Main']),paragraph('26 October 2026 · 20 donors: 7 ND, 7 Aab+, 6 T1D',styles['BodySmall']),paragraph('Corrected cell counts, acquisition-edge exclusions, structure-qualified cell measurements and union areas are linked to the app’s saved composition UMAP, Leiden 0.5 clusters, fixed pseudotime and physical neighborhoods. Donors are the biological replicates. Endpoint-specific cohorts remain distinct from the objects eligible for the embedding.',styles['BodySmall']),paragraph('This report presents measured values, estimates, confidence intervals and tests. It assigns no biological interpretations or original-plan support/exclusion branches. H5 and H7 use the specified categorical/count-based amendments. Quantitative CODEX INS is not an endpoint.',styles['BodySmall']),Spacer(1,5),paragraph('Primary comparison statistics',styles['Section'])]
    data=[[paragraph(x,styles['Cell']) for x in ['Hypothesis','Comparison','Donors','Objects / tiles','p','BH q']]]
    for r in t.itertuples():
        nd='—' if pd.isna(r.n_donors) else str(int(r.n_donors));no='—' if pd.isna(r.n_objects) else str(int(r.n_objects))
        data.append([paragraph(x,styles['Cell']) for x in [r.hypothesis,TITLE[r.hypothesis],nd,no,fmt(r.p_value),fmt(r.p_bh_ten_slots)]])
    table=Table(data,colWidths=[60,350,55,75,110,110],repeatRows=1,hAlign='LEFT')
    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eef4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#778899')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f7f9fb')]),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
    story += [table,Spacer(1,8),paragraph('BH adjustment retains ten declared slots, including unavailable tests. These are the current direct disease comparisons and approved H9/H10 tests; they are not a relabeling of the original ten conditional-direction tests. Asterisks use each figure’s stated comparison family; the ten-primary-test adjustment is listed here. H8 reports tissue tiles; H10 uses donor clinical observations. Within-disease estimates and supporting comparison families are reported separately in legends and CSV files. H9 retains an unavailable primary p when its prespecified shared-range model has no residual degrees of freedom.',styles['Caption']),PageBreak()]
    for r in records:
        story.append(paragraph(f"Figure {r['figure']} · {figure_title(r)}",styles['Section']))
        imagepath=OUT/r['png'];im=PILImage.open(imagepath);w,h=im.size;k=min(mw/w,(mh-55)/h)
        story.append(Image(str(imagepath),width=w*k,height=h*k));story.append(PageBreak())
    story.append(paragraph('Figure legends',styles['Main']))
    for r in records:
        text=(OUT/r['legend']).read_text()
        story.append(paragraph(f"Figure {r['figure']} · {figure_title(r)}",styles['Section']))
        # Keep CSV-derived exact statistics intact and wrap long lines rather than
        # truncating them or assigning conclusions.
        for chunk in re.split(r'\n\s*\n',text):
            if chunk.strip():story.append(paragraph(chunk.strip(),styles['Caption']))
        story.append(Spacer(1,12))
    def footer(canvas,doc):
        canvas.setFont('DejaVu',7);canvas.setFillColor(colors.HexColor('#555555'));canvas.drawString(32,22,'Islet Explorer · measured outputs and statistics · 26 October 2026');canvas.drawRightString(width-32,22,str(doc.page))
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    return path,len(PdfReader(path).pages)

def build_html(t,records):
    lines=['<!doctype html><html><head><meta charset="utf-8"><title>App hypotheses: figures and statistics</title><style>body{font:16px system-ui;max-width:1400px;margin:32px auto;padding:0 24px;color:#202832}h1{font-size:28px}h2{margin-top:42px}img{width:100%;height:auto;border:1px solid #ddd}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}th{background:#eef2f6}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:13px/1.6 system-ui}nav{position:sticky;top:0;background:white;padding:10px;z-index:2}a{color:#245a8b}details{padding:12px;background:#f7f9fb}</style></head><body><h1>App-based hypotheses: figures and statistics</h1><p>20 donors: 7 ND, 7 Aab+, 6 T1D. No biological interpretation or branch assignment.</p><nav><a href="App_Hypotheses_Figures_Statistics.pdf">Full PDF</a> · <a href="primary_statistics.csv">Primary statistics CSV</a> · <a href="prism/App_Hypotheses_plot_data.pzfx">Prism data</a></nav>',t[['hypothesis','comparison','test','n_donors','n_objects','p_value','p_bh_ten_slots','status']].to_html(index=False,escape=True,na_rep='Unavailable')]
    for r in records:
        lines += [f'<h2 id="figure{r["figure"]}">Figure {r["figure"]}: {html.escape(figure_title(r))}</h2>',f'<img src="{r["png"]}" alt="{html.escape(figure_title(r))}" loading="lazy">',f'<p><a href="{r["pdf"]}">PDF</a> · <a href="{r["svg"]}">SVG</a></p>',f'<details><summary>Legend, methods and exact statistics</summary><pre>{html.escape((OUT/r["legend"]).read_text())}</pre></details>']
    lines+=['</body></html>'];(OUT/'index.html').write_text('\n'.join(lines))

def audit():
    p=ROOT/'docs/T1D_research_hypotheses_20261022.md';text=p.read_text();parts=re.split(r'(?m)^### (H\d+) ',text);records=[]
    for i in range(1,len(parts),2):
        h,body=parts[i:i+2];fields={k:v.strip() for k,v in re.findall(r'\*\*([^*]+)\.\*\*\s*(.*?)(?=\n\n|\Z)',body,re.S) if k=='If supported' or k.startswith('If a meaningful') or k=='Same-data follow-up'}
        records.append({'hypothesis':h,'original_conditions':fields,'branch_judgment':'not assigned','amendment':'Original quantitative INS endpoint superseded; categorical cell-quantity amendment used' if h=='H5' else 'Observed proximity with Beta-density endpoint replaces original spatial-null/INS specification' if h=='H7' else 'INS-response component omitted; retained Beta-quantity component used' if h=='H6' else 'No tissue-INS follow-up' if h=='H10' else 'Retained endpoint; current disease-comparison presentation and approved methods'})
    assert len(records)==10;(OUT/'original_conditional_rules.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n')
    return p

def main():
    t,sources=primary();r=figures();orig=audit();pdf,pages=build_pdf(t,r);build_html(t,r)
    (OUT/'README.md').write_text('# App-based H1–H10 figures and statistics\n\nOpen `index.html` for the figure gallery or `App_Hypotheses_Figures_Statistics.pdf` for all figures and manuscript-style legends. `primary_statistics.csv` gives the ten current primary comparisons and BH adjustment. Within-disease estimates, sensitivities, direct pairwise comparisons and all exact plotted observations remain in each hypothesis folder and the Prism export.\n\nThe original app embedding, cluster labels and fixed runtime DPT snapshot are reused; corrected cell and geometry measurements join by persistent IDs. The common representation display contains 1,126 objects from 20 donors. Each hypothesis retains its own model cohort, including eligible objects outside the embedding. No live data pipeline was rebuilt. Original accepted app phenotypes and recovered structural qualification remain separate fields; structural sidecars are not added to endocrine/immune totals.\n\nH5 tests the specified categorical Beta-presence/composition/state quantities. H7 uses log1p Beta-cell density and includes zero-Beta objects when proximity is defined. These amendments do not test within-cell INS abundance. No quantitative INS endpoint or biological interpretation is included. Original if–then conditions are preserved without branch assignments. H9’s shared-range primary test remains unavailable when the specified model has zero residual degrees of freedom; the supporting all-range fit is labeled separately.\n\nAll standalone figures have PNG, PDF and SVG versions, with separate legends. Exact statistical and plot-data files retain donor/object IDs, units, eligibility and missing values. Prism contains data tables; fitted R/Python analyses and graph layouts are not embedded.\n')
    (OUT/'source').mkdir(exist_ok=True)
    scriptroot=ROOT/'scripts/research'
    pending=list(scriptroot.glob('app_*_20261026*'));copied=set()
    while pending:
        p=pending.pop()
        if p in copied or not p.is_file():continue
        shutil.copy2(p,OUT/'source'/p.name);copied.add(p)
        if p.suffix=='.py':
            for node in ast.walk(ast.parse(p.read_text())):
                modules=[a.name for a in node.names] if isinstance(node,ast.Import) else [node.module] if isinstance(node,ast.ImportFrom) and node.module else []
                for module in modules:
                    dependency=scriptroot/(module.split('.')[0]+'.py')
                    if dependency.exists():pending.append(dependency)
    (OUT/'reports').mkdir(exist_ok=True);shutil.copy2(ROOT/'docs/app_based_hypothesis_plan_20261026.md',OUT/'reports/app_based_hypothesis_plan_20261026.md')
    manifest={'figures':len(r),'pdf_pages':pages,'primary_slots':10,'available_primary_p':int(t.p_value.notna().sum()),'branch_assignments':0,'biological_interpretations':0,'pdf_sha256':sha(pdf),'primary_sources':{str(p.relative_to(ROOT)):sha(p) for p in sources},'original_plan_sha256':sha(orig)}
    (OUT/'report_validation.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
