#!/usr/bin/env python3
"""Export plotted observations as GraphPad-schema-validated Prism data and CSV."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import shutil
import zipfile
import numpy as np
import pandas as pd
from lxml import etree as ET

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/prism_export_20261026'
OLD=ROOT/'output/h1_h4_literature_methods_20261025'
NEW=ROOT/'output/h1_h4_disease_comparisons_20261026'
REST=ROOT/'output/h5_h10_original_plan_20261026'
AUDIT=ROOT/'output/h1_h10_conditional_audit_20261026'
GROUPS=['ND','Aab+','T1D']
NS='http://graphpad.com/prism/Prism.htm'
TABLES=[]
FIGURES=[]
INS_NOTICE=OUT/'INS_MEASUREMENT_NOTICE.md'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def node(parent,name,text=None,**attributes):
    e=ET.SubElement(parent,ET.QName(NS,name),**{k:str(v) for k,v in attributes.items()})
    if text is not None:e.text=str(text)
    return e


def read(path):return pd.read_csv(path,dtype={'case_id':str})


def register(key,title,frame,numeric,description,source,paired=False,x=None,ci=False,table_type=None):
    frame=frame.copy()
    affected=(
        (key.startswith('H5_') and key not in ['H5_beta_density','H5_hormone_negative_cell_density']) or
        (key.startswith('H6_') and any(x in key for x in ['INS','analysis_variants','composition_size'])) or
        (key.startswith('H7_') and any(x in key for x in ['coefficients','conditional','contrasts','object_scatter'])) or
        key.startswith('H10_tissue_INS'))
    if INS_NOTICE.exists() and affected:
        description='HISTORICAL INS-BASED ANALYSIS: superseded as a current quantitative endpoint following the owner\'s acquisition clarification; see INS_MEASUREMENT_NOTICE.md. '+description
    frame.to_csv(OUT/'csv'/f'{key}.csv',index=False,float_format='%.17g')
    frame.to_csv(OUT/'txt'/f'{key}.txt',index=False,sep='\t',float_format='%.17g')
    TABLES.append(dict(key=key,title=title,frame=frame,numeric=numeric,description=description,
                       source=str(Path(source).relative_to(ROOT)),paired=paired,x=x,ci=ci,table_type=table_type))


def column_groups(key,title,source,value,description,subset=None,factor=1):
    d=read(source).sort_values(['donor_status','case_id'])
    if subset is not None:d=d.loc[subset(d)]
    assert not d.case_id.duplicated().any()
    z=pd.DataFrame({'Donor_ID':d.case_id.to_numpy()})
    for g in GROUPS:z[g]=np.where(d.donor_status.eq(g),d[value]*factor,np.nan)
    register(key,title,z,GROUPS,description,source)


def paired_groups(key,title,source,value,category,labels,factor,description,subset=None):
    d=read(source)
    if subset is not None:d=d.loc[subset(d)]
    for g in GROUPS:
        q=d.loc[d.donor_status.eq(g)].pivot(index='case_id',columns=category,values=value)
        q=q.reindex(columns=list(labels)).rename(columns=labels)*factor
        q.index.name='Donor_ID'
        z=q.reset_index()
        register(key+'_'+g.replace('+','pos'),title+' — '+g,z,list(labels.values()),description,source,paired=True)


def xy_groups(key,title,source,x,y,description,groups=GROUPS):
    d=read(source)
    for g in groups:
        q=d.loc[d.donor_status.eq(g)].copy()
        q.insert(0,'Row_ID',q.case_id if 'case_id' in q else np.arange(1,len(q)+1).astype(str))
        q=q[['Row_ID',x,y]]
        register(key+'_'+g.replace('+','pos'),title+' — '+g,q,[y],description,source,x=x)


def figure(key,image,legend,legend_text=None):
    image=Path(image);legend=Path(legend)
    if not image.exists():return
    for f in image.parent.glob(image.stem+'.*'):
        if f.suffix in ['.png','.pdf','.svg']:shutil.copy2(f,OUT/'plots'/(key+f.suffix))
    text=legend_text if legend_text is not None else legend.read_text() if legend.exists() else ''
    (OUT/'legends'/f'{key}.txt').write_text(text)
    FIGURES.append(dict(key=key,description=text,source=str(image.relative_to(ROOT))))


def build_tables():
    column_groups('H1_immune_beta_donor_coefficients','Immune localization and beta-cell density',
                  NEW/'h1/donor_coefficients.csv','standardized_beta',
                  'One adjusted standardized regression coefficient per donor. Independent disease groups; rows must not be paired across disease columns.')
    column_groups('H2_coliv_immune_donor_coefficients','Boundary ColIV coverage and core immune density',
                  NEW/'h2/primary_donor_coefficients.csv','standardized_beta',
                  'One adjusted standardized boundary-coverage coefficient per donor. Independent disease groups.')
    column_groups('H3_capillary_density_difference','Immune density near capillary crossings',
                  NEW/'h3/h3_primary_donor_values.csv','difference_density_mm2',
                  'One near-crossing minus reference immune-density difference per donor, cells/mm2; primary 10-µm compartments.')
    d=read(NEW/'h4/disease_slopes.csv')
    d=d.loc[d.model.eq('primary')].set_index('donor_status').loc[GROUPS].reset_index()
    z=pd.DataFrame({'Disease_state':GROUPS,'X':range(1,4),'Estimate':d['source_z.trend'],
                    'Upper_95CI':d['upper.CL'],'Lower_95CI':d['lower.CL']})
    register('H4_lgals3_beta_model_estimates','LGALS3 source allocation and beta-cell density',z,
             ['Estimate','Upper_95CI','Lower_95CI'],
             'Model estimates and confidence limits, not individual donor observations. XY table uses actual upper/lower 95% limits, not SD or SEM.',
             NEW/'h4/disease_slopes.csv',x='X',ci=True)
    s=OLD/'h1/paired_donor_measurements.csv'
    paired_groups('H1_beta_presence_paired','Endocrine objects with beta cells (%)',s,'beta_presence','core_immune',
                  {False:'Peri_only',True:'Core_present'},100,
                  'Matched donor rows. Percentage of endocrine objects with at least one typed Beta cell. Peri only = immune cells in the 20-µm peri band and none in core. Core present = at least one core immune cell.')
    paired_groups('H1_beta_density_paired','Donor mean core beta-cell density (cells/mm2)',s,'beta_density','core_immune',
                  {False:'Peri_only',True:'Core_present'},1,
                  'Matched donor rows. Arithmetic mean of individual object Beta densities; Beta-absent objects contribute zero. Density is not per whole tissue section.')
    # Useful supporting donor measurements are included even when their comparisons are not significant.
    paired_groups('H2_coliv_coverage_paired','Median boundary ColIV coverage (%)',OLD/'h2/h2_donor_coverage.csv',
                  'median_coverage','core_immune_present',{False:'Peri_only',True:'Core_present'},100,
                  'Matched donor rows. Median nonvascular ColIV-positive coverage in the 2-µm boundary band, grouped by immune location.')
    paired_groups('H4_source_by_beta_presence','Donor mean myeloid LGALS3 share (%)',
                  OLD/'h4/donor_beta_presence_source_summary.csv','mean_myeloid_fraction','beta_present',
                  {False:'Beta_absent',True:'Beta_present'},100,
                  'Matched donor rows. Mean myeloid share of measured macrophage-plus-epithelial LGALS3 signal within each Beta-presence category.')
    h3=read(OLD/'h3/h3_donor_densities.csv')
    for g in GROUPS:
        q=h3.loc[h3.width_um.eq(10)&h3.donor_status.eq(g),['case_id','near_density_mm2','reference_density_mm2']].rename(columns={'case_id':'Donor_ID'})
        register('H3_near_reference_'+g.replace('+','pos'),'Capillary crossing immune density — '+g,q,
                 ['near_density_mm2','reference_density_mm2'],
                 'Matched donor rows. Counts and areas summed within each donor before density calculation. Units cells/mm2.',
                 OLD/'h3/h3_donor_densities.csv',paired=True)
    figure('H1_disease_comparison',NEW/'h1/immune_beta_association_by_disease.png',NEW/'h1/figure_legend.txt')
    figure('H2_disease_comparison',NEW/'h2/coliv_disease_comparison.png',NEW/'h2/figure_legend.md')
    figure('H3_disease_comparison',NEW/'h3/h3_density_difference_by_disease.png',NEW/'h3/figure_legend.txt')
    figure('H4_disease_comparison',NEW/'h4/source_association_by_disease.png',NEW/'h4/figure_legend.md')
    for metric in ['beta_presence','beta_density']:
        figure('H1_'+metric+'_paired',OLD/f'h1/h1_{metric}_paired.png',OLD/f'h1/h1_{metric}_paired_legend.txt')
    h5=REST/'h5/plot_data_donor_coefficients.csv'
    if h5.exists():
        column_groups('H5_immune_INS_donor_coefficients','Immune localization and beta-cell insulin',h5,'standardized_beta',
                      'One adjusted standardized regression coefficient per estimable donor. Within-Beta mean INS is the outcome. Nineteen donors: ND 7, Aab+ 7, T1D 5. Donor 6563 cannot support the adjusted regression and has no coefficient, not a zero.')
        h5m=REST/'h5/plot_data_donor_ins_measurements.csv'
        for value,title,definition in [
            ('beta_density','Beta-cell abundance','Typed Beta cells divided by summed eligible endocrine-core area, cells/mm2; not whole tissue density.'),
            ('mean_beta_INS','Retained beta-cell insulin','Mean scaled INS across typed Beta cells per donor, weighted by cell count. A staining measurement, not insulin secretion.'),
            ('hormone_negative_cell_density','CHGA-positive hormone-negative cells','CHGA-positive hormone-negative cells divided by summed eligible endocrine-core area, cells/mm2.')]:
            column_groups('H5_'+value,title,h5m,value,'One observation per donor; all 20 donors. '+definition)
        src=REST/'h5/plot_data_insulin_decomposition.csv';d=read(src)
        ids=d.case_id.drop_duplicates().tolist();components=d.component.drop_duplicates().tolist()
        z=d.pivot(index='case_id',columns='component',values='signal_density').reindex(ids)[components].reset_index().rename(columns={'case_id':'Donor_ID'})
        register('H5_insulin_decomposition','Sources of core insulin signal',z,components,
                 'Grouped table: rows are donors, columns are additive signal components for stacked bars. Units are sum of scaled cell means/mm2 endocrine-core area, not protein mass. Source residual cells were excluded from accepted cell typing and are shown separately. Column components are not independent sample groups.',src,table_type='TwoWay')
        for key in ['immune_localization_beta_insulin','insulin_quantity_and_state','insulin_signal_decomposition']:
            figure('H5_'+key,REST/f'h5/{key}.png',REST/f'h5/{key}_legend.txt')
    h6=REST/'h6/plot_data_observed_adjacency.csv'
    if h6.exists():
        column_groups('H6_observed_adjacency','Observed beta–alpha/delta adjacency',h6,'mean_observed_combined_fraction',
                      'Descriptive donor mean percentage of object-level Beta–Alpha plus Beta–Delta graph edges among all endocrine graph edges (20-µm centroid graph); source fractions multiplied by 100 to match plotted percentages. Nineteen donors have eligible objects; donor 6563 remains a missing value. Observed adjacency; permutation-based and conditional model data are supplied in separate tables.',factor=100)
        figure('H6_observed_adjacency',REST/'h6/observed_beta_alpha_delta_adjacency.png',REST/'h6/observed_beta_alpha_delta_adjacency_legend.txt')
    h7=REST/'h7/donor_observed_proximity.csv'
    if h7.exists():
        column_groups('H7_observed_proximity','Observed nerve–vessel proximity near Beta-containing objects',h7,'mean_proximity20_percent',
                      'Descriptive donor mean percentage of peri B3TUBB-positive area within 20 µm of CD31 structures. All 20 donors; 946 Beta-positive objects with measured nerve area/proximity. Undefined proximity where nerve area is absent is excluded, not assigned zero. No statistical or spatial-null test is included.')
        figure('H7_observed_proximity',REST/'h7/observed_proximity_by_disease.png',REST/'h7/figure_legend.md')
    h9=REST/'h9/prism_donor_points.csv'
    if h9.exists():
        column_groups('H9_coverage_donor_coefficients','Boundary coverage and core immune density',h9,
                      'coverage_core_immune_slope','Donor regression coefficients reused from H2. ND is a descriptive reference; the H9 adjusted comparison is Aab+ versus T1D.')
        xy_groups('H9_coverage_beta_quantity','Boundary coverage association and sampled beta-cell quantity',h9,
                  'tissue_beta_density_mm2','coverage_core_immune_slope',
                  'One donor per XY row. X: unique edge-eligible tissue Beta cells/mm2. Y: adjusted standardized coverage coefficient. Aab+ and T1D are independent disease groups.',groups=['Aab+','T1D'])
        xy_groups('H9_fitted_lines','Fitted coverage association versus sampled beta-cell quantity',
                  REST/'h9/prism_adjusted_regression_lines.csv','tissue_beta_density_mm2','fitted_slope',
                  'Fitted ANCOVA regression line coordinates, not additional observations. Do not use these rows as replicates.',groups=['Aab+','T1D'])
        figure('H9_coverage_and_beta_quantity',REST/'h9/coverage_association_and_beta_quantity.png',REST/'h9/figure_legend.txt')
    h8=REST/'h8/donor_field_coefficients.csv'
    if h8.exists():
        inventory=read(h8)[['model','grid_um','field']].drop_duplicates()
        for _,r in inventory.iterrows():
            column_groups('H8_'+r['model']+'_'+r['field'],'Small beta-bearing clusters and '+r['field']+' context — '+r['model'],
                          h8,'standardized_beta',
                          f"One adjusted standardized donor regression coefficient per disease column; {int(r['grid_um'])}-µm tissue grid. Zero-object tiles are retained. Field: {r['field']}; specification: {r['model']}. Predictor fields are standardized within donor; outcome density across the selected grid/specification. Independent donor groups.",
                          subset=lambda d,r=r:d.model.eq(r['model'])&d.field.eq(r['field']))
        src=REST/'h8/donor_object_class_densities.csv'
        for endpoint in read(src).endpoint.drop_duplicates():
            column_groups('H8_density_'+endpoint,'Endocrine-object tissue density — '+endpoint,src,'density_mm2',
                          'One observed object count per donor divided by valid sampled Tissue area, objects/mm2. All 20 donors; zeros retained. Object class/size definition: '+endpoint,
                          subset=lambda d,e=endpoint:d.endpoint.eq(e))
        legends=REST/'h8/figure_legends.md'
        if legends.exists():
            legend_text=legends.read_text()
            common=legend_text.split('\nprimary_200um,')[0]
            for p in sorted((REST/'h8').glob('*_field_comparisons.png')):
                model=p.stem.removesuffix('_field_comparisons')
                statistics='\n'.join(line for line in legend_text.splitlines() if line.startswith(model+','))
                figure('H8_'+p.stem,p,legends,common+'\n\nSpecification: '+model+'\n'+statistics+'\n')
            statistics='\n'.join(line for line in legend_text.splitlines() if ' density: Kruskal' in line)
            figure('H8_object_class_densities',REST/'h8/object_class_densities.png',legends,common+'\n\n'+statistics+'\n')
    h10=REST/'h10/prism_clinical_input_points.csv'
    if h10.exists():
        d=read(h10)
        for g in GROUPS:
            q=d.loc[d.donor_status.eq(g)].copy()
            z=pd.DataFrame({'Donor_ID':q.case_id,'tissue_beta_density_mm2':q.tissue_beta_density_mm2,
                            'Measured_C_peptide_ng_mL':q.C_peptide_ng_mL.where(~q.C_peptide_censored),
                            'Below_limit_bound_ng_mL':q.C_peptide_ng_mL.where(q.C_peptide_censored)})
            register('H10_clinical_inputs_'+g.replace('+','pos'),'C-peptide and sampled beta-cell quantity — '+g,z,
                     ['Measured_C_peptide_ng_mL','Below_limit_bound_ng_mL'],
                     'Descriptive XY data, with measured values separated from below-limit bounds. Donor 6550 is <0.02 ng/mL, not an exact 0.02 concentration. Draw its bound with a downward triangle/arrow. Model and prediction data are supplied in separate tables.',h10,x='tissue_beta_density_mm2')
        figure('H10_clinical_inputs',REST/'h10/clinical_quantity_overview.png',REST/'h10/clinical_quantity_figure_legend.txt')


def ci_table(key,title,source,label,estimate,lower,upper,description,subset=None):
    d=read(source)
    if subset is not None:d=d.loc[subset(d)]
    d=d.reset_index(drop=True)
    z=pd.DataFrame({'Label':d[label].astype(str),'X':np.arange(1,len(d)+1),
                    'Estimate':d[estimate],'Upper_95CI':d[upper],'Lower_95CI':d[lower]})
    register(key,title,z,['Estimate','Upper_95CI','Lower_95CI'],
             description+' Values are estimates with actual upper/lower 95% confidence limits, not replicate observations.',source,x='X',ci=True)


def build_approved_inference_tables():
    h6=REST/'h6'
    src=h6/'plot_data_permutation_donor_summary.csv'
    column_groups('H6_excess_adjacency','Observed minus permuted endocrine adjacency',src,'excess_percentage_points',
                  'Donor means of observed combined Beta–Alpha/Delta edge fraction minus the fixed-composition label-permutation mean, multiplied by 100. Units percentage-point difference, not a regression effect per ten percentage points. Donor 6563 has no eligible graph and remains missing.')
    for group in GROUPS:
        d=read(src);d=d[d.donor_status.eq(group)]
        z=d[['case_id','observed_percent','null_percent']].rename(columns={'case_id':'Donor_ID'})
        register('H6_observed_permuted_'+group.replace('+','pos'),'Observed and permuted adjacency — '+group,z,
                 ['observed_percent','null_percent'],'Paired donor means of observed and label-permutation expected adjacency (%). Graphs preserve cell composition and geometry.',src,paired=True)
    column_groups('H6_donor_INS_coefficients','Endocrine adjacency and beta-cell insulin',h6/'plot_data_primary_donor_coefficients.csv',
                  'standardized_beta','Fully adjusted independent donor coefficients; 15 estimable donors, including one T1D donor. Missing/nonestimable donors are retained in source eligibility tables.')
    ci_table('H6_mixed_INS_estimates','Adjacency–insulin model estimates',h6/'plot_data_primary_mixed_model_slopes.csv',
             'donor_status','contact_z.trend','lower.CL','upper.CL',
             'Separate random-intercept sensitivity with unchanged fixed covariates. Assumes no residual donor slope heterogeneity within states; primary donor/random-slope diagnostics are retained.')
    ci_table('H6_analysis_variants','Conditional adjacency associations across variants',h6/'plot_data_sensitivity_estimates.csv',
             'label','standardized_beta','lower','upper','Each row is a separately specified analysis; exact model choices, donor/object counts and tests are in the source table. The horizontal coefficient plot can be reconstructed from these estimates.')
    src=h6/'plot_data_matched_composition_size.csv';d=read(src)
    for group in GROUPS:
        z=d.loc[d.donor_status.eq(group),['case_id','nonpositive_INS','positive_INS']].rename(columns={'case_id':'Donor_ID'})
        register('H6_composition_size_'+group.replace('+','pos'),'Matched composition/size insulin measurements — '+group,z,
                 ['nonpositive_INS','positive_INS'],'Paired donor measurements aggregated within the recorded matched composition/size strata. Columns are nonpositive versus positive excess adjacency; INS is scaled within-Beta cell mean.',src,paired=True)
    for stem in ['observed_vs_permuted_adjacency','excess_adjacency_by_disease','adjacency_insulin_association',
                 'adjacency_analysis_variants','matched_composition_size_insulin']:
        figure('H6_'+stem,h6/f'{stem}.png',h6/f'{stem}_legend.txt')

    h7=REST/'h7'
    for model in read(h7/'conditional_plot_data.csv').model.drop_duplicates():
        column_groups('H7_'+model+'_donor_coefficients','Nerve–vessel proximity and insulin — '+model,
                      h7/'individual_donor_coefficients.csv','standardized_beta',
                      'Fully adjusted independent donor coefficients. Nonestimable coefficients remain blank, with donor identifiers retained. Exact specifications and estimability are in source tables.',
                      subset=lambda d,m=model:d.model.eq(m))
        ci_table('H7_'+model+'_conditional','Adjusted proximity–insulin estimates — '+model,
                 h7/'conditional_plot_data.csv','donor_status','estimate','lower','upper',
                 'State-specific conditional estimates using the source-designated model or mean of estimable independent donor coefficients; method is preserved in the source table.',subset=lambda d,m=model:d.model.eq(m))
    ci_table('H7_random_intercept_conditional','Proximity–insulin random-intercept sensitivity',
             h7/'prism_csv/random_intercept_conditional_estimates.csv','donor_status','estimate','lower','upper',
             'Supporting model with the same fixed covariates and donor random intercept; assumes no residual within-state donor slope heterogeneity.')
    ci_table('H7_random_intercept_contrasts','Proximity–insulin disease contrasts in sensitivity',
             h7/'prism_csv/random_intercept_disease_contrasts.csv','contrast','estimate','lower','upper','Disease contrasts from the separate random-intercept sensitivity.')
    for value,title in [('nerve_percent','Peri nerve area (%)'),('cd31_percent','Peri CD31 area (%)'),('mean_core_area_um2','Mean endocrine-core area (µm²)')]:
        column_groups('H7_'+value,title,h7/'donor_abundance_geometry.csv',value,'Descriptive donor companion measurement; one row per donor.')
    src=h7/'prism_csv/observed_proximity_INS_objects.csv'
    for donor,d in read(src).groupby('case_id',sort=False):
        z=d[['islet_id','nerve_near_cd31_20um','ins_mean_beta_core']].rename(columns={'islet_id':'Object_ID'})
        register('H7_object_scatter_'+str(donor),'Observed proximity and insulin — donor '+str(donor),z,['ins_mean_beta_core'],
                 'Raw object scatter shown in the donor facets. X is observed nerve-area fraction within 20 µm of CD31; Y is scaled mean INS in typed Beta cells. Object points are not independent donor replicates.',src,x='nerve_near_cd31_20um')
    legends=(h7/'inference_figure_legends.md').read_text()
    for stem in read(h7/'plot_manifest.csv').figure.drop_duplicates():
        if stem=='observed_proximity_by_disease.png':continue
        heading='## '+stem
        specific=legends.split(heading,1)[1].split('\n## ',1)[0] if heading in legends else legends
        common=legends.split('\n## ',1)[0]
        figure('H7_'+Path(stem).stem,h7/stem,h7/'inference_figure_legends.md',common+'\n\n'+heading+'\n'+specific)

    h10=REST/'h10'
    for stem in ['clinical_held_out_prediction','clinical_disease_conditioned_prediction',
                 'tissue_INS_held_out_prediction','tissue_INS_disease_conditioned_prediction']:
        src=h10/f'prism_{stem}.csv';d=read(src)
        for group in GROUPS:
            z=d.loc[d.donor_status.eq(group),['case_id','observed','baseline_prediction','augmented_prediction']].rename(columns={'case_id':'Donor_ID'})
            register('H10_'+stem+'_'+group.replace('+','pos'),'Held-out prediction — '+stem+' — '+group,z,
                     ['baseline_prediction','augmented_prediction'],
                     'Observed versus donor-held-out predictions from baseline and fixed-score models. C-peptide values are natural logarithms; tissue_INS values are scaled mean INS. Each point is one donor, and both models predict the same held-out donor.',src,x='observed')
        figure('H10_'+stem,h10/f'{stem}.png',h10/f'{stem}_legend.txt')
    src=h10/'prism_T1D_individual_clinical_predictions.csv';d=read(src)
    z=pd.DataFrame({'Donor_ID':d.case_id,'Measured_C_peptide_ng_mL':d.C_peptide_ng_mL.where(~d.C_peptide_censored),
                    'Below_limit_bound_ng_mL':d.C_peptide_ng_mL.where(d.C_peptide_censored),
                    'Baseline_predicted_median':d.baseline_predicted_median_ng_mL,
                    'With_score_predicted_median':d.augmented_predicted_median_ng_mL})
    register('H10_individual_T1D','Individual T1D observations and predictions',z,list(z.columns[1:]),
             'Rows are individual donors; columns are observed values, separately flagged below-limit bounds, and exponentiated predictions. The censored bound is not a measured concentration and has no squared error. Donor 6550 is predicted from all 19 measured donors; each measured donor is held out from its own model.',src,table_type='TwoWay')
    ci_table('H10_censored_score','Censored-model spatial-score coefficients',h10/'prism_censored_model_score_association.csv',
             'analysis','estimate','lower_95','upper_95','Censored-normal likelihood using all 20 donors, with the below-0.02 donor contributing log probability below the bound. These are full-sample coefficients.')
    for stem in ['T1D_individual_clinical_predictions','censored_model_score_association']:
        figure('H10_'+stem,h10/f'{stem}.png',h10/f'{stem}_legend.txt')

    src=AUDIT/'association_audit/plot_data_conditional_means.csv'
    for h in ['H1','H2','H5']:
        ci_table(h+'_conditional_mean_audit','Supporting state-specific conditional means — '+h,src,
                 'donor_status','mean_beta','lower_95','upper_95',
                 'Supporting mean of independent adjusted donor coefficients, using one-sample t intervals. The existing primary disease comparisons remain separate. Exact raw/BH p values and donor signs are in the source table.',subset=lambda d,h=h:d.hypothesis.eq(h))
    figure('conditional_association_means',AUDIT/'association_audit/conditional_association_means.png',AUDIT/'association_audit/figure_legend.txt')


def prism_file():
    root=ET.Element(ET.QName(NS,'GraphPadPrismFile'),nsmap={None:NS},PrismXMLVersion='5.00',LowestVersionThatCanOpen='8.00')
    created=node(root,'Created')
    node(created,'OriginalVersion',CreatedByProgram='Islet Explorer data export',CreatedByVersion='2026.10.26',
         DateTime=datetime.now(timezone.utc).isoformat(),MinVersion='8.00')
    seq=node(root,'InfoSequence');node(seq,'Ref',ID='Info0',Selected='1')
    info=node(root,'Info',ID='Info0');node(info,'Title','Islet Explorer plot data')
    notice=('INS MEASUREMENT AMENDMENT: The owner states these images were not acquired for quantitative INS comparison. H5 INS analyses, H6/H7 INS-response models and H10 tissue-INS follow-ups are historical; see INS_MEASUREMENT_NOTICE.md in the accompanying package. No replacement endpoint results are included. ' if INS_NOTICE.exists() else '')
    node(info,'Notes',notice+'Data tables only. Figure PNG/PDF/SVG files and full legends are supplied separately. Blank cells are absent observations, not zero. Column disease comparisons are unpaired. Peri/core and Beta-absent/present tables are paired by Donor_ID. Model estimates with confidence limits are not replicate observations. No Prism analyses or graph layouts are embedded.')
    seq=node(root,'TableSequence')
    for i,spec in enumerate(TABLES):node(seq,'Ref',ID=f'Table{i}',Selected='1' if i==0 else '0')
    expected=[];expected_x=[]
    for i,spec in enumerate(TABLES):
        ci=spec['ci'];f=spec['frame']
        table=node(root,'Table',ID=f'Table{i}',XFormat='numbers' if spec['x'] else 'none',
                   YFormat='upper-lower-limits' if ci else 'replicates',Replicates='1',TableType=spec['table_type'] or ('XY' if spec['x'] else 'OneWay'))
        node(table,'Title',spec['title'])
        meta=node(table,'Info');node(meta,'Title','Definitions');node(meta,'Notes',spec['description'])
        rows=node(node(table,'RowTitlesColumn'),'Subcolumn')
        for value in f.iloc[:,0]:node(rows,'d',str(value))
        if spec['x']:
            xcol=node(table,'XColumn',Decimals='10',Subcolumns='1');node(xcol,'Title',spec['x']);sub=node(xcol,'Subcolumn')
            for value in f[spec['x']]:node(sub,'d',format(float(value),'.17g'))
            expected_x.append([float(value) for value in f[spec['x']]])
        cols=[spec['numeric']] if ci else [[c] for c in spec['numeric']]
        for names in cols:
            ycol=node(table,'YColumn',Decimals='10',Subcolumns=len(names));node(ycol,'Title','Adjusted coefficient' if ci else names[0])
            for name in names:
                sub=node(ycol,'Subcolumn');node(sub,'Title',name)
                values=[]
                for v in f[name]:
                    text='' if pd.isna(v) else format(float(v),'.17g')
                    node(sub,'d',text);values.append(None if pd.isna(v) else float(v))
                expected.append(values)
    path=OUT/'Islet_Explorer_plot_data.pzfx'
    ET.ElementTree(root).write(str(path),xml_declaration=True,encoding='UTF-8',pretty_print=True)
    # GraphPad's published XSD has overlapping formatted-text/wildcard branches
    # in cell type d, which libxml2 rejects as non-deterministic. For our plain
    # text/numeric cells, remove only that optional extension wildcard. This is
    # a stricter subset of the published content model; keep the original file.
    schema_doc=ET.parse(str(OUT/'support/Prism8XMLSchema.xsd'))
    xs={'xs':'http://www.w3.org/2001/XMLSchema'}
    extension=schema_doc.xpath('//xs:complexType[@name="d"]/xs:complexContent/xs:extension',namespaces=xs)[0]
    wildcard=extension.find('{http://www.w3.org/2001/XMLSchema}sequence')
    assert wildcard.get('minOccurs')=='0'
    extension.remove(wildcard)
    schema_doc.write(str(OUT/'support/Prism8XMLSchema_plain_cells.xsd'),encoding='UTF-8',xml_declaration=True)
    schema=ET.XMLSchema(schema_doc)
    parsed=ET.parse(str(path));schema.assertValid(parsed)
    ns={'p':NS}
    actual=[]
    for col in parsed.xpath('//p:Table/p:YColumn/p:Subcolumn',namespaces=ns):
        actual.append([None if n.text in [None,''] else float(n.text) for n in col.findall(f'{{{NS}}}d')])
    assert actual==expected,'Numeric roundtrip mismatch'
    actual_x=[[float(n.text) for n in col.findall(f'{{{NS}}}d')]
              for col in parsed.xpath('//p:Table/p:XColumn/p:Subcolumn',namespaces=ns)]
    assert actual_x==expected_x,'X-coordinate roundtrip mismatch'
    for table,spec in zip(parsed.xpath('//p:Table',namespaces=ns),TABLES):
        labels=table.xpath('./p:RowTitlesColumn/p:Subcolumn/p:d/text()',namespaces=ns)
        assert labels==spec['frame'].iloc[:,0].astype(str).tolist(),'Row identity mismatch'
    return {'graphpad_schema_plain_cell_subset_valid':True,
            'schema_note':'Published Prism 8 XSD has a non-deterministic optional wildcard in cell type d; validation removes only that wildcard, producing a stricter plain-cell subset.',
            'numeric_roundtrip_exact':True,'row_identity_roundtrip_exact':True,'tables':len(TABLES),
            'numeric_values':sum(v is not None for c in actual for v in c)+sum(len(c) for c in actual_x),
            'prism_application_open_test':False}


def main():
    for name in ['csv','txt','plots','legends','source_data','statistics','reports']: (OUT/name).mkdir(parents=True,exist_ok=True)
    build_tables()
    build_approved_inference_tables()
    assert len({s['key'] for s in TABLES})==len(TABLES),'Duplicate table key'
    assert len({s['key'] for s in FIGURES})==len(FIGURES),'Duplicate figure key'
    assert all(f['description'].strip() for f in FIGURES),'A delivered figure lacks its separate legend'
    validation=prism_file()
    index=[]
    for s in TABLES:
        index.append({k:v for k,v in s.items() if k!='frame'})
    (OUT/'table_dictionary.json').write_text(json.dumps(index,indent=2,ensure_ascii=False)+'\n')
    pd.DataFrame([{k:v for k,v in s.items() if k not in ['numeric']} for s in index]).to_csv(OUT/'table_dictionary.csv',index=False)
    # Copy source tables and exact statistical outputs separately from graph-ready data.
    sources=[OLD/'h1/analysis_objects.csv',OLD/'h2/h2_analysis_objects.csv',
             OLD/'h4/object_source_features.csv',NEW/'h1/donor_coefficients.csv',NEW/'h2/primary_donor_coefficients.csv']
    sources += [p for p in [REST/'h9/prism_donor_points.csv',REST/'h9/prism_adjusted_regression_lines.csv',
                            REST/'h10/prism_clinical_input_points.csv',REST/'h5/plot_data_donor_coefficients.csv',
                            REST/'h5/plot_data_donor_ins_measurements.csv',REST/'h5/plot_data_insulin_decomposition.csv',
                            REST/'h5/donor_regression_eligibility.csv',REST/'h5/object_ins_features.csv',
                            REST/'h6/plot_data_observed_adjacency.csv',
                            REST/'h7/donor_observed_proximity.csv',REST/'h7/matched_object_plot_data.csv',
                            REST/'h7/object_neurovascular_features.csv',
                            REST/'h8/donor_field_coefficients.csv',REST/'h8/donor_object_class_densities.csv',
                            REST/'h8/donor_census.csv',REST/'h8/full_object_census.csv',REST/'h8/all_tissue_tiles.csv'] if p.exists()]
    sources += [ROOT/s['source'] for s in TABLES]
    sources += [p for h in ['h6','h7','h10'] for p in (REST/h).glob('*.csv')]
    sources += list((AUDIT/'association_audit').glob('*.csv'))
    sources = list(dict.fromkeys(sources))
    for p in sources:shutil.copy2(p,OUT/'source_data'/(p.parent.name+'_'+p.name))
    for h,paths in {
        'H1':['donor_disease_global_test.csv','donor_pairwise_disease_tests.csv'],
        'H2':['primary_disease_overall_test.csv','primary_disease_pairwise_tests.csv'],
        'H3':['h3_global_test.csv','h3_pairwise_disease_tests.csv'],
        'H4':['model_anova.csv','disease_slope_comparisons.csv'],
    }.items():
        for name in paths:shutil.copy2(NEW/h.lower()/name,OUT/'statistics'/(h+'_'+name))
    shutil.copy2(NEW/'primary_disease_comparison_tests.csv',OUT/'statistics/primary_disease_comparison_tests.csv')
    if (REST/'h9/adjusted_group_contrasts.csv').exists():
        shutil.copy2(REST/'h9/adjusted_group_contrasts.csv',OUT/'statistics/H9_adjusted_group_contrasts.csv')
    for p in (REST/'h5').glob('*test*.csv'):
        shutil.copy2(p,OUT/'statistics'/('H5_'+p.name))
    for p in (REST/'h8').glob('*test*.csv'):
        shutil.copy2(p,OUT/'statistics'/('H8_'+p.name))
    for h in ['h6','h7','h10']:
        for p in (REST/h).glob('*.csv'):
            if any(word in p.stem for word in ['test','comparison','contrast','coefficient','slopes','anova','eligibility','conditional_means','direction','model_metrics','error_by_group']):
                shutil.copy2(p,OUT/'statistics'/(h.upper()+'_'+p.name))
    for p in AUDIT.glob('*.csv'):shutil.copy2(p,OUT/'statistics'/p.name)
    for h,name in [('h5','README.md'),('h5','analysis_specification.json'),('h5','validation.json'),
                   ('h5','original_if_then_conditions.json'),('h5','if_then_decision.json'),
                   ('h6','README.md'),('h7','README.md'),('h8','README.md'),('h9','RESULTS.txt'),
                   ('h6','original_plan_execution.csv'),('h6','validation.json'),
                   ('h7','original_plan_branch_coverage.csv'),('h7','inference_validation.json'),
                   ('h9','conditional_branch.json'),('h9','support_assessment.json'),
                   ('h10','README.md'),('h10','approved_analysis_specification.json'),('h10','validation.json'),
                   ('h10','original_if_then_conditions.json'),
                   ('h10','if_then_interpretation.json')]:
        p=REST/h/name
        if p.exists():shutil.copy2(p,OUT/'reports'/(h.upper()+'_'+name))
    for h in ['H5','H10']:
        for name in ['README.md','if_then_decision.json','if_then_interpretation.json']:
            p=OUT/'reports'/(h+'_'+name)
            if p.exists():
                content=p.read_text()
                for ref in ['original_if_then_conditions.json','if_then_decision.json','if_then_interpretation.json']:
                    content=content.replace(ref,h+'_'+ref)
                p.write_text(content)
    stale=OUT/'reports/H10_conventional_model_specification.txt'
    if stale.exists():stale.unlink()  # Superseded generated proposal copy, not a source analysis file.
    shutil.copy2(AUDIT/'README.md',OUT/'reports/plan_execution_audit.md')
    shutil.copy2(AUDIT/'original_plan_conditions.json',OUT/'reports/original_plan_conditions.json')
    shutil.copy2(AUDIT/'validation.json',OUT/'reports/plan_audit_validation.json')
    measurement_review=ROOT/'docs/ins_measurement_review_20261026.md'
    if measurement_review.exists():shutil.copy2(measurement_review,OUT/'reports/INS_method_review.md')
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    status=read(AUDIT/'current_test_inventory.csv').rename(columns={'execution_status':'status'})
    status.to_csv(OUT/'analysis_status.csv',index=False)
    lines=['# Islet Explorer: Prism plot-data export','',
        'Open `Islet_Explorer_plot_data.pzfx` in GraphPad Prism. It contains native data tables with donor row titles. '
        'CSV and tab-delimited copies are in `csv/` and `txt/`. All current primary H1–H4 plots are included, regardless of significance, '
        'along with the revised paired H1 figures and additional paired donor measurements. H5–H10 approved analyses, model sensitivities, '
        'individual T1D clinical observations and supporting conditional-association plots are included. '
        'Reports provide numerical results and plan execution without assigning biological interpretation.','',
        'The file contains data tables, not saved Prism graph layouts or rerunnable Prism versions of the fitted R models. '
        'Plots and manuscript legends are supplied in `plots/` and `legends/`. Exact test results are in `statistics/`; analysis specifications, '
        'technical diagnostics and the original-plan execution audit are in `reports/`. The audit distinguishes completed steps, approved adaptations and follow-ups not run.','',
        '## Current analysis status','',
        '| Hypothesis | Current test | Unadjusted p | Status |','|---|---|---:|---|']
    if INS_NOTICE.exists():
        lines[2:2]=['**Measurement amendment:** H5 INS analyses, H6/H7 INS-response models and H10 tissue-INS follow-ups are historical following the owner\'s acquisition clarification. '
                    'Read [INS_MEASUREMENT_NOTICE.md](INS_MEASUREMENT_NOTICE.md) before using these tables. The numerical inventory below records the previous specification; revised endpoint models have not been substituted.','']
    for _,s in status.iterrows():
        p='—' if pd.isna(s.p) else f'{s.p:.4g}'
        lines.append(f"| {s.hypothesis} | {s.relationship} | {p} | {s.status} |")
    lines += ['','The inventory follows the current disease-comparison presentation and the approved H10 nested test; it is not the original ten conditional-direction tests. '
        'The separate plan audit records those distinctions. `statistics/current_test_inventory.csv` includes BH adjustment across ten declared slots with unavailable H6/H9 primary p values left missing. '
        'H6 original main-association and random-intercept sensitivity outputs, and H9 supplementary all-range models, remain separately identified. '
        'Previously declared within-figure pairwise and sensitivity families are retained in their legends. No biological branch decisions are assigned.','',
        '## Import and measurement notes','',
        '- Disease-comparison columns are independent samples. Do not pair their rows across disease groups.',
        '- Paired tables preserve one donor per row across both conditions.',
        '- Empty cells remain missing; they are not zeros. Genuine zero measurements are retained.',
        '- H4 main-figure values are fitted estimates with upper/lower 95% confidence limits. They are not donor replicates or standard deviations.',
        '- Percentages are exported on the 0–100 scale. Density units and coefficient definitions are stated in each table description.',
        '- The data file passed a stricter plain-cell subset of GraphPad’s published XML schema and exact numeric roundtrip checks. The published schema has an overlapping optional wildcard that prevents compilation in libxml2; only that optional branch was removed for validation. It has not been opened in a running Prism application in this environment.','',
        'GraphPad documents [text-file importing](https://www.graphpad.com/guides/prism/latest/user-guide/using_transferring_data_via_txt_files.htm) '
        'and its [Prism XML schema](https://www.graphpad.com/support/faq/prism-xml-style-sheet-and-schema-that-define-the-pzfx-format/).','',
        '## Table index','', '| Table | Layout | Meaning |','|---|---|---|']
    for s in index:
        layout='XY with confidence limits' if s['ci'] else 'XY' if s['x'] else 'Grouped components' if s['table_type']=='TwoWay' else 'Paired columns' if s['paired'] else 'Independent columns'
        lines.append(f"| {s['key']} | {layout} | {s['description']} |")
    lines+=['','## Figures and legends','']
    for f in FIGURES:lines += [f"![{f['key']}](plots/{f['key']}.png)",'',f['description'],'']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n')
    manifest={'created_date':'2026-10-26','validation':validation,'source_files':{str(p.relative_to(ROOT)):sha(p) for p in sources},
              'export_files':{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'}}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    archive=ROOT/'output/Islet_Explorer_Prism_export_20261026.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.rglob('*'):
            if p.is_file() and 'support' not in p.relative_to(OUT).parts:z.write(p,Path(OUT.name)/p.relative_to(OUT))
    print(json.dumps(validation));print(archive)


if __name__=='__main__':main()
