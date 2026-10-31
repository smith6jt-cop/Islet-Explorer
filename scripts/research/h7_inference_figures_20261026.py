#!/usr/bin/env python3
"""Conventional H7 figures, exact plot data and original-plan interpretation."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h7-inference-mpl')
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
OUT=Path('output/h5_h10_original_plan_20261026/h7')
GROUPS=['ND','Aab+','T1D']; LABELS=['ND','AAb+','T1D']; COLORS=['#4477AA','#EE7733','#228833']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
dc=pd.read_csv(OUT/'individual_donor_coefficients.csv');means=pd.read_csv(OUT/'fallback_conditional_means.csv')
mixed=pd.read_csv(OUT/'conditional_slopes.csv');fp=pd.read_csv(OUT/'fallback_disease_pairs.csv')
mp=pd.read_csv(OUT/'disease_slope_comparisons.csv');elig=pd.read_csv(OUT/'model_donor_eligibility.csv')
spec=pd.read_csv(OUT/'inference_model_specifications.csv')
names={'primary_20um':'Nerve–vessel proximity and β-cell INS','distance_10um':'Nerve–vessel proximity at 10 µm and β-cell INS',
 'distance_5um':'Nerve–vessel proximity at 5 µm and β-cell INS','overlap':'Nerve–vessel overlap and β-cell INS',
 'GAP43_20um':'GAP43-supported proximity and β-cell INS','raw_INS_20um':'Nerve–vessel proximity and raw β-cell INS',
 'corrected_INS_20um':'Nerve–vessel proximity and corrected β-cell INS'}
manifest=[];summary=[];donor_plot=[];conditional_plot=[]
def save(fig,stem):
    for ext in ['png','pdf']:fig.savefig(OUT/f'{stem}.{ext}',dpi=250,bbox_inches='tight')
    plt.close(fig)
def dotbox(ax,table,y):
    vv=[table.loc[table.donor_status.eq(g),y].dropna().to_numpy() for g in GROUPS]
    bp=ax.boxplot(vv,positions=[1,2,3],widths=.5,showfliers=False,patch_artist=True,
      medianprops={'color':'#222','linewidth':1.3},whiskerprops={'linewidth':1},capprops={'linewidth':1})
    for j,(v,c) in enumerate(zip(vv,COLORS),1):
        bp['boxes'][j-1].set(facecolor=c,edgecolor=c,alpha=.18)
        ax.scatter(j+np.random.default_rng(20261026+j).uniform(-.08,.08,len(v)),v,c=c,s=36,edgecolor='white',linewidth=.4,zorder=3)
    ax.set_xticks([1,2,3],LABELS);ax.set_xlabel('Disease state')
def groupci(ax,table):
    for i,(g,c) in enumerate(zip(GROUPS,COLORS),1):
        r=table.loc[table.donor_status.eq(g)].iloc[0]
        ax.errorbar(i,r.estimate,yerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',c=c,capsize=4,ms=7)
    ax.axhline(0,color='#777',lw=.8,ls='--');ax.set_xticks([1,2,3],LABELS);ax.set_xlabel('Disease state')
    ax.set_ylabel('Adjusted standardized regression β')
def wide(table,y):
    cols={}
    for g in GROUPS:
        z=table.loc[table.donor_status.eq(g)].sort_values('case_id').reset_index(drop=True)
        cols[g+'_donor_id']=z.case_id;cols[g+'_value']=z[y]
    return pd.DataFrame(cols)
def tablemd(d):return d.to_markdown(index=False,floatfmt='.5g')
for model in spec.model:
    t=dc.loc[dc.model.eq(model)].copy();est=t.loc[t.estimable].copy()
    is_mixed=model in set(mixed.model)
    if is_mixed:
        s=mixed.loc[mixed.model.eq(model)].rename(columns={'exposure_z.trend':'estimate','lower.CL':'lower','upper.CL':'upper','p.value':'p'}).copy()
        source='conditional_slopes.csv';testsource='disease_slope_comparisons.csv'
        s['method']='Donor random intercept and slope mixed model'
        p=mp.loc[mp.model.eq(model)]
    else:
        s=means.loc[means.model.eq(model)].rename(columns={'mean':'estimate'}).copy()
        source='fallback_conditional_means.csv';testsource='fallback_disease_pairs.csv'
        s['method']='Mean of fully adjusted, estimable donor coefficients'
        p=fp.loc[fp.model.eq(model)]
    conditional_plot.append(s[['model','donor_status','estimate','lower','upper','p','p_BH','method']])
    donor_plot.append(t)
    fig,axs=plt.subplots(1,2,figsize=(9,4.2))
    dotbox(axs[0],est,'standardized_beta');axs[0].axhline(0,color='#777',lw=.8,ls='--')
    axs[0].set_ylabel('Adjusted standardized regression β');axs[0].set_title('Individual donor estimates')
    groupci(axs[1],s);axs[1].set_title('Conditional association, 95% CI')
    fig.suptitle(names[model],fontsize=13);fig.tight_layout(rect=(0,0,1,.94))
    stem=model+'_ins_associations';save(fig,stem)
    fn=f'prism_csv/{model}_adjusted_donor_coefficients.csv';wide(t,'standardized_beta').to_csv(OUT/fn,index=False)
    fn2=f'prism_csv/{model}_conditional_estimates.csv';s[['donor_status','estimate','lower','upper','p','p_BH']].to_csv(OUT/fn2,index=False)
    manifest.extend([dict(figure=stem+'.png',panel='individual_donor_estimates',plot_data_csv=str(OUT/fn),long_data_csv='individual_donor_coefficients.csv',filter='model='+model,donor_id='case_id',x='donor_status',y='standardized_beta',statistics_csv=testsource,stars_column='p_BH',biological_replicate='donor',inference='full-covariate donor coefficient; nonestimable donors exported as NA'),
    dict(figure=stem+'.png',panel='conditional_association',plot_data_csv=str(OUT/fn2),long_data_csv=source,filter='model='+model,donor_id='see model_donor_eligibility.csv',x='donor_status',y='estimate; lower; upper',statistics_csv=source,stars_column='',biological_replicate='donor',inference=s.method.iloc[0])])
    summary.append(dict(model=model,primary=model=='primary_20um',method=s.method.iloc[0],disease_tests_BH_min=p.p_BH.min(),conditional_tests_BH_min=s.p_BH.min(),
                        estimable_donors_ND=int((est.donor_status=='ND').sum()),estimable_donors_Aab=int((est.donor_status=='Aab+').sum()),estimable_donors_T1D=int((est.donor_status=='T1D').sum())))
pd.concat(conditional_plot).to_csv(OUT/'conditional_plot_data.csv',index=False)
pd.DataFrame(summary).to_csv(OUT/'inference_summary.csv',index=False)
# Supporting common-within-disease slope regression, retained separately.
ri=pd.read_csv(OUT/'random_intercept_conditional_slopes.csv').rename(columns={'exposure_z.trend':'estimate','lower.CL':'lower','upper.CL':'upper','p.value':'p'})
rp=pd.read_csv(OUT/'random_intercept_disease_comparisons.csv').rename(columns={'lower.CL':'lower','upper.CL':'upper','p.value':'p'})
fig,axs=plt.subplots(1,2,figsize=(9.8,4.3));groupci(axs[0],ri);axs[0].set_title('Conditional association, 95% CI')
for i,r in rp.iterrows():
    axs[1].errorbar(r.estimate,i,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',c='#555',capsize=4)
axs[1].axvline(0,color='#777',lw=.8,ls='--');axs[1].set_yticks(range(3),['ND − AAb+','ND − T1D','AAb+ − T1D']);axs[1].invert_yaxis()
axs[1].set_xlabel('Difference in standardized regression β');axs[1].set_title('Disease differences, 95% CI')
fig.suptitle('Nerve–vessel proximity and β-cell INS',fontsize=13);fig.tight_layout(rect=(0,0,1,.91));save(fig,'random_intercept_ins_associations')
ri.to_csv(OUT/'prism_csv/random_intercept_conditional_estimates.csv',index=False);rp.to_csv(OUT/'prism_csv/random_intercept_disease_contrasts.csv',index=False)
for panel,fn,src,y in [('conditional_association','random_intercept_conditional_estimates.csv','random_intercept_conditional_slopes.csv','estimate; lower; upper'),('disease_comparison','random_intercept_disease_contrasts.csv','random_intercept_disease_comparisons.csv','estimate; lower; upper')]:
    manifest.append(dict(figure='random_intercept_ins_associations.png',panel=panel,plot_data_csv=str(OUT/'prism_csv'/fn),long_data_csv=src,filter='all rows',donor_id='see primary_20um_model_objects.csv',x='donor_status or contrast',y=y,statistics_csv=src,stars_column='',biological_replicate='donor random intercept',inference='Supporting shared-slope sensitivity'))
# Raw object scatter for each donor; these points are not independent replicates.
objs=pd.read_csv(OUT/'primary_20um_model_objects.csv');scatter_cols=['islet_id','case_id','donor_status','nerve_near_cd31_20um','ins_mean_beta_core','nerve_area_fraction','cd31_area_fraction','core_area_um2','core_perimeter_um','n_immune_nearby','core_area_um2','peri_area_um2']
objs.loc[:,list(dict.fromkeys(scatter_cols))].to_csv(OUT/'prism_csv/observed_proximity_INS_objects.csv',index=False)
fig,axs=plt.subplots(4,5,figsize=(12,9),sharex=True,sharey=True)
ordered=objs.assign(order=objs.donor_status.map(dict(zip(GROUPS,range(3))))).sort_values(['order','case_id']).case_id.drop_duplicates()
for ax,ca in zip(axs.flat,ordered):
    z=objs.loc[objs.case_id.eq(ca)];g=z.donor_status.iloc[0];i=GROUPS.index(g)
    ax.scatter(z.nerve_near_cd31_20um,z.ins_mean_beta_core,s=11,alpha=.65,color=COLORS[i]);ax.set_title(f'{LABELS[i]} · {ca}',fontsize=10);ax.set_xlim(-.03,1.03);ax.set_ylim(-.03,1.03)
fig.supxlabel('Observed nerve area within 20 µm of CD31 (fraction)');fig.supylabel('Mean β-cell INS (served scaled signal)')
fig.suptitle('Observed nerve–vessel proximity and β-cell INS',fontsize=14);fig.tight_layout(rect=(.02,.02,1,.96));save(fig,'observed_proximity_INS_by_donor')
manifest.append(dict(figure='observed_proximity_INS_by_donor.png',panel='all donor facets',plot_data_csv=str(OUT/'prism_csv/observed_proximity_INS_objects.csv'),long_data_csv='primary_20um_model_objects.csv',filter='all 946 rows',donor_id='case_id',x='nerve_near_cd31_20um',y='ins_mean_beta_core',statistics_csv='',stars_column='',biological_replicate='donor; raw object points are descriptive',inference='none; no fitted line'))
# Existing outcome-independent median-size image examples; no extraction rerun.
fig,axs=plt.subplots(2,3,figsize=(10,7));image_rows=[]
for j,(ca,lab) in enumerate([(6516,'ND'),(6450,'AAb+'),(6533,'T1D')]):
    z=np.load(OUT/f'example_{ca}.npz');raw=z['raw'].astype(float);lo,hi=np.percentile(raw,[1,99.8]);v=np.clip((raw-lo)/max(hi-lo,1e-12),0,1)
    axs[0,j].imshow(v,cmap='gray');axs[0,j].set_title(f'{lab} · {ca}\nB3TUBB')
    rgb=np.zeros((*raw.shape,3));rgb[z['cd31']]=[.2,.55,.95];rgb[z['nerve']]=[.95,.65,.1];rgb[z['nerve'] & z['cd31']]=[.85,.2,.85]
    axs[1,j].imshow(rgb);axs[1,j].contour(z['peri'],levels=[.5],colors=['white'],linewidths=.5);axs[1,j].set_title(str(z['islet_id']))
    width=20/float(z['pixel_um']);yy=raw.shape[0]-14;axs[1,j].plot([12,12+width],[yy,yy],color='white',lw=3);axs[1,j].text(12,yy-6,'20 µm',color='white',fontsize=8)
    image_rows.append(dict(case_id=ca,islet_id=str(z['islet_id']),archive=f'example_{ca}.npz',raw_channel='B3TUBB',display_low=lo,display_high=hi,pixel_um=float(z['pixel_um']),selection='Prespecified donor and median core area among prepared objects; selected before INS association fitting'))
for ax in axs.flat:ax.set_axis_off()
fig.suptitle('INS-independent nerve masks and CD31 structures',fontsize=14);fig.tight_layout(rect=(0,0,1,.95));save(fig,'nerve_mask_examples')
pd.DataFrame(image_rows).to_csv(OUT/'image_examples.csv',index=False)
manifest.append(dict(figure='nerve_mask_examples.png',panel='B3TUBB and masks',plot_data_csv=str(OUT/'image_examples.csv'),long_data_csv='image_examples.csv',filter='all 3 rows',donor_id='case_id',x='pixel position',y='pixel position',statistics_csv='',stars_column='',biological_replicate='illustrative images',inference='none'))
# Density and geometry companions in exactly the matched set, all 20 donors.
ab=objs.groupby(['case_id','donor_status']).agg(objects=('islet_id','size'),nerve_percent=('nerve_area_fraction',lambda x:100*x.mean()),cd31_percent=('cd31_area_fraction',lambda x:100*x.mean()),mean_core_area_um2=('core_area_um2','mean'),mean_perimeter_um=('core_perimeter_um','mean')).reset_index()
ab.to_csv(OUT/'donor_abundance_geometry.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(10.5,4))
for ax,v,title,yl in zip(axs,['nerve_percent','cd31_percent','mean_core_area_um2'],['Nerve area','CD31 area','Endocrine object area'],['Mean peri area (%)','Mean peri area (%)','Mean core area (µm²)']):
    dotbox(ax,ab,v);ax.set_title(title);ax.set_ylabel(yl)
    fn=f'prism_csv/donor_{v}.csv';wide(ab,v).to_csv(OUT/fn,index=False)
    manifest.append(dict(figure='abundance_geometry_companions.png',panel=v,plot_data_csv=str(OUT/fn),long_data_csv='donor_abundance_geometry.csv',filter='all 20 donors',donor_id='case_id',x='donor_status',y=v,statistics_csv='',stars_column='',biological_replicate='donor',inference='descriptive only'))
fig.suptitle('Structural abundance and geometry near β-cell-containing objects',fontsize=13);fig.tight_layout(rect=(0,0,1,.93));save(fig,'abundance_geometry_companions')
original_manifest=pd.read_csv(OUT/'plot_manifest.csv');original_manifest=original_manifest.loc[original_manifest.figure.eq('observed_proximity_by_disease.png')]
pd.concat([original_manifest,pd.DataFrame(manifest)],ignore_index=True).to_csv(OUT/'plot_manifest.csv',index=False)
branches=pd.DataFrame([
 ('Primary test','Original H7: association of excess-over-null alignment with within-Beta INS, adjusted for separate nerve/vessel abundance, area, perimeter and immune exposure','Completed using the approved observed-proximity substitution; donor-fallback primary and random-intercept sensitivity retained','Not assigned'),
 ('If supported','Original plan condition: a reserve-associated spatial configuration distinguishes endocrine niches despite similar structural density','Conditional coefficients, confidence intervals and tests are reported; this condition has not been assigned','Not assigned'),
 ('If a meaningful effect is excluded','Original plan condition: alignment adds little independent hormone-phenotype information at the measured precision','No minimum effect-size threshold was specified; coefficient intervals are reported; this condition has not been assigned','Not assigned'),
 ('Technical estimability','Full adjustment requires exposure variation, full design rank and positive residual degrees of freedom','Individual coefficients estimable for 7 ND, 7 AAb+ and 2 T1D donors; all 20 donors retained in descriptive and random-intercept data','Not assigned'),
 ('GAP43 corroboration','Original follow-up: compare B3TUBB with GAP43-supported measurements','Completed for 20-µm proximity; 922 matched objects','Not assigned'),
 ('Distance and overlap','Approved sensitivity: 5/10-µm proximity and direct nerve–vessel overlap','Completed; numerical results and plots supplied','Not assigned'),
 ('INS measurement sensitivity','Approved sensitivity: raw cytoplasmic and unscaled spillover-corrected within-Beta mean INS','Completed using existing H5 features','Not assigned'),
 ('Density and geometry','Original test: separate nerve/vessel abundance, area and perimeter adjustment','All required terms retained; donor abundance and geometry tables and figures supplied','Not assigned'),
 ('Core/peri/distal area and fiber length','Original follow-up: compare core/peri/distal fiber area, length and alignment','Peri area and observed proximity completed; core/distal INS-independent fiber measurements and length not run','Not assigned'),
 ('Image omission and inherited masks','Original follow-up: uncertain-image omission and inherited/replacement mask comparison','Source provenance and fixed-rule examples retained; image-omission and inherited-mask comparison not run','Not assigned'),
 ('Other original-plan endpoints','Original plan lists sympathetic identity, neurotransmission, reinnervation and glucagon counterregulation','No tests of these endpoints were run','Not assigned')],columns=['branch','original_condition_or_method','completed_or_not_run','branch_judgment'])
branches.to_csv(OUT/'original_plan_branch_coverage.csv',index=False)
definitions='''All figures use endocrine objects. ND, AAb+ and T1D denote the three disease states. The source set contains 4,179 eligible paired endocrine objects; 1,794 contain measured Beta cells. Proximity and within-Beta INS are jointly measured in 946 objects (509 ND, 383 AAb+, 54 T1D), from all 20 donors. The 848 other Beta-containing objects have no measured B3TUBB-positive peri area; proximity is undefined there and is never assigned zero. No Beta-containing objects in this input fail the 20-µm acquisition-support check. The absent-mask category records zero area above the image threshold.

The INS-independent nerve mask uses a per-image B3TUBB Otsu threshold and connected-component area ≥2 µm², excludes all source endocrine cores, and is assessed in Tissue-supported 20-µm peri regions. INS is never used to construct the nerve mask. The exposure is the fraction of nerve-positive peri area within the stated distance of an existing CD31-supported structure. Overlap uses zero distance; GAP43 corroboration restricts both exposure numerator and denominator to GAP43-positive B3TUBB area. Existing source image, channel, threshold and geometry provenance is retained. Observed proximity replaces the original excess-over-spatial-null endpoint with the user's approval. No spatial-null calculation was run.

For every inference specification, exposure, separate nerve and CD31 area fractions, log core area, log perimeter and log1p nearby immune density are standardized within donor. Nearby immune density is the broad immune count in core plus peri divided by their combined area, in cells/mm². The response is within-Beta mean INS standardized across all objects in that specification; primary INS is the existing served cohort-scaled signal. A standardized regression β expresses outcome SD per within-donor exposure SD conditional on all five covariates. It is not a correlation coefficient. Both geometry covariates and both abundance covariates remain in every fully adjusted model.

The original candidate model includes disease × exposure and random donor intercepts/slopes. The correlated covariance is tried first, followed by independent components only for singularity/convergence. For primary 20-µm proximity the correlation reaches −1 and the independent random-slope variance reaches zero. The prespecified fallback therefore fits the complete seven-parameter regression separately per donor; coefficients require exposure variation, full design rank and positive residual degrees of freedom. Seven ND, seven AAb+ and two T1D donors meet those conditions. T1D donors 6520 and 6563 have constant 20-µm exposure, while 6550 and 6551 have six observations for seven parameters. These donors' descriptive measurements are retained. The two estimable T1D donors are 6533 and 6534. No adjustment term is dropped to obtain a sparse T1D coefficient.

Disease differences in the fallback use Kruskal–Wallis and three two-sided Wilcoxon rank-sum tests, with BH correction across the three pairs per specification. Conditional group mean coefficients use the one-sample t test and 95% Student t interval, with a separate BH family of three disease-specific conditional tests per specification; this assumes approximately normal independent donor coefficients. The T1D group contains two estimable donor coefficients. The direct-overlap sensitivity has a nonsingular random-intercept/slope model; its conditional estimates and disease differences use unchanged lmerTest Satterthwaite inference, with the same separate three-test BH families. Sensitivity specifications are disclosed separately rather than selected by significance. Exact raw and BH-adjusted p values are listed below. No stars are shown. Intervals are pointwise 95% intervals, not multiplicity-adjusted intervals.

The supporting random-intercept-only sensitivity retains exactly the same 946 objects, transformations, outcome, exposure, disease interaction and all five covariates. It assumes no residual donor slope heterogeneity within disease states. Reported standard errors and intervals use that common-within-disease-slope assumption. It does not replace the declared donor-fallback primary. This standard mixed model uses lmerTest/Satterthwaite inference and separate BH families of three conditional slopes and three disease contrasts.

Methods are represented in the supplied literature: Barlow2025 PDF p22 (within-donor predictor standardization, globally standardized response and donor random intercept/slope), Barlow PDF pp8/12 (lmerTest Satterthwaite), Bodenmiller2025 donor mixed regression and donor-level Wilcoxon/FDR, Richardson2023 thresholded neural area and Kruskal–Wallis, and Heide PDF p32 (one-sample t test). The narrow peri nerve/CD31 distance endpoint is an adaptation of thresholded-area and spatial-distance methods. No pipeline rebuild or manual-review gate was introduced.
'''
text=['# Figure legends and complete statistics\n',definitions]
for model in spec.model:
    s=pd.concat(conditional_plot).query('model == @model')
    text.append(f'## {model}_ins_associations.png\n\n{names[model]}. Left: one fully adjusted coefficient per estimable donor; boxes show median and interquartile range with whiskers to values within 1.5 interquartile ranges. Right: group conditional estimate and pointwise 95% confidence interval. The dashed line marks zero association. Nonestimable donor rows remain in the exact CSV with the reason recorded. The exposure and response are specified in inference_model_specifications.csv.\n\nConditional estimates:\n\n'+tablemd(s[['donor_status','estimate','lower','upper','p','p_BH']]))
    if model in set(mixed.model):
        q=mp.loc[mp.model.eq(model)];om=pd.read_csv(OUT/'mixed_model_tests.csv');om=om.loc[om.model.eq(model) & om.term.eq('exposure_z:donor_status')]
        text.append('\nOverall disease-by-exposure test:\n\n'+tablemd(om[['NumDF','DenDF','F value','Pr(>F)']])+'\n\nDisease comparisons:\n\n'+tablemd(q[['contrast','estimate','lower.CL','upper.CL','p.value','p_BH']]))
    else:
        q=fp.loc[fp.model.eq(model)];k=pd.read_csv(OUT/'fallback_disease_tests.csv').query('model == @model').iloc[0]
        text.append(f'\nKruskal–Wallis p={k.p:.6g}; disease pairwise comparisons:\n\n'+tablemd(q[['group1','group2','n1','n2','median1','median2','p','p_BH']]))
text.extend(['## random_intercept_ins_associations.png\n\nSupporting shared-slope sensitivity. Left: conditional standardized proximity slopes by disease and 95% confidence intervals. Right: direct disease slope differences and 95% confidence intervals. All 20 donors and 946 matched objects enter this model. Precision assumes a common adjusted slope within each disease state. Exact raw and BH-adjusted p values are listed below. No stars are shown.\n\n'+tablemd(ri[['donor_status','estimate','lower','upper','p','p_BH']])+'\n\n'+tablemd(rp[['contrast','estimate','lower','upper','p','p_BH']])+'\n\nOverall disease-by-exposure test:\n\n'+tablemd(pd.read_csv(OUT/'random_intercept_model_tests.csv').query("term == 'exposure_z:donor_status'")[["NumDF","DenDF","F value","Pr(>F)"]]),
'## observed_proximity_INS_by_donor.png\n\nEach point is an endocrine object, with one panel per donor. The x-axis is observed 20-µm proximity; the y-axis is served-scaled within-Beta mean INS. These are unadjusted descriptive measurements, not independent biological replicates or fitted correlations. No tests or trend lines are shown. Exact rows are in prism_csv/observed_proximity_INS_objects.csv.\n',
'## nerve_mask_examples.png\n\nTop: original B3TUBB channel with per-image 1st–99.8th percentile display scaling. Bottom: measured B3TUBB mask in orange, CD31 geometry in blue, and overlap in purple; white contours show the assessed peri region. Scale bars are 20 µm. The fixed example donors are ND 6516, AAb+ 6450 and T1D 6533; the median-core-area prepared object within each donor was selected before association fitting. Examples illustrate mask measurements in the source objects, rather than selecting Beta-positive objects for an association figure. The AAb+ and T1D examples contain no measured Beta cells and are therefore absent from the INS regression. The T1D example has no measured nerve-positive area inside the peri contour, so its proximity is undefined. Orange pixels outside the peri contour are visible for context and do not contribute to the measured nerve area. These examples are not a validation sample. Image IDs and display limits are supplied in image_examples.csv.\n',
'## abundance_geometry_companions.png\n\nOne point is the mean among matched objects for one donor, retaining all 20 donors. Panels show peri nerve area fraction, peri CD31 area fraction and core area; boxes show median/interquartile range and 1.5-IQR whiskers. These are descriptive companions for the abundance and geometry terms in the adjusted model. No tests are performed. No additional tests of these companion measurements were run.\n'])
(OUT/'inference_figure_legends.md').write_text('\n\n'.join(text))
primary=means.query("model == 'primary_20um'");primarykw=pd.read_csv(OUT/'fallback_disease_tests.csv').query("model == 'primary_20um'").iloc[0]
readme=f'''# H7: observed neurovascular proximity and within-Beta INS

The approved observed-proximity regressions and prepared sensitivity models have been run. The primary donor-level Kruskal–Wallis p value is {primarykw.p:.6g}. Numerical estimates, intervals and p values are reported below; biological interpretation and original-plan branch judgments are not assigned.

Primary conditional estimates (standardized regression β, mean of estimable donor coefficients):

{tablemd(primary[['donor_status','donors','mean','lower','upper','p_BH']])}

Individual donor coefficients are estimable for seven ND, seven AAb+ and two T1D donors. T1D donors 6520 and 6563 have constant 20-µm exposure; donors 6550 and 6551 have six observations for seven parameters. All six T1D donors remain in the descriptive data and random-intercept sensitivity. The latter uses the same fixed covariates and assumes a common adjusted slope within each disease state.

Random-intercept sensitivity estimates:

{tablemd(ri[['donor_status','estimate','lower','upper','p_BH']])}

The 5/10-µm, direct-overlap, GAP43-supported, raw-INS and corrected-INS specifications are recorded in inference_model_specifications.csv. Their complete coefficient, interval and test tables are reproduced in inference_figure_legends.md. All prepared specifications have corresponding plots and plot CSVs.

The original plan's conditions are identified in original_plan_branch_coverage.csv. Branch judgment is not assigned. Peri area, observed proximity, GAP43 corroboration, distance/overlap and INS-measurement sensitivities were completed. Core/distal fiber measurements, fiber length, inherited-mask comparison and image-omission analyses were not run.

Reproduce with h7_nerve_features_20261026.py (existing completed image extraction), h7_descriptive_20261026.py, h7_models_20261026.R, h7_random_intercept_sensitivity_20261026.R and h7_inference_figures_20261026.py. All model specifications, donor eligibility, raw/scaled predictor tables, R models, complete test tables, figures, legends and Prism CSVs are retained here. plot_manifest.csv maps each panel to its exact data and statistics.
'''
(OUT/'README.md').write_text(readme)
# Independently verify the primary donor coefficients and donor-group tests.
errs=[]
for ca,z in objs.groupby('case_id'):
    row=dc.loc[dc.model.eq('primary_20um') & dc.case_id.eq(ca)].iloc[0]
    if not row.estimable:continue
    cols=['exposure_z','nerve_area_fraction_z','cd31_area_fraction_z','log_core_area_z','log_perimeter_z','log_nearby_immune_density_z']
    X=np.column_stack([np.ones(len(z)),z[cols].to_numpy()]);b=np.linalg.lstsq(X,z.response_z,rcond=None)[0]
    errs.append(abs(b[1]-row.standardized_beta))
v=[dc.loc[dc.model.eq('primary_20um') & dc.estimable & dc.donor_status.eq(g),'standardized_beta'].to_numpy() for g in GROUPS]
check={'max_abs_numpy_vs_R_primary_donor_beta':max(errs),'scipy_KW_p':stats.kruskal(*v).pvalue,'R_KW_p':float(primarykw.p),
       'matched_objects':len(objs),'descriptive_donors':int(objs.case_id.nunique()),'primary_estimable_donors':dict(zip(GROUPS,map(len,v))),
       'no_primary_covariates_dropped':True,'spatial_null_performed':False,'all_prepared_models_plotted':True}
assert max(errs)<1e-10 and abs(check['scipy_KW_p']-check['R_KW_p'])<1e-12
(OUT/'inference_validation.json').write_text(json.dumps(check,indent=2)+'\n')
files=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name!='inference_output_manifest.json':files.append({'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/'inference_output_manifest.json').write_text(json.dumps(files,indent=2)+'\n')
print(json.dumps(check,indent=2))
