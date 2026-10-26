#!/usr/bin/env python3
"""Standard paired donor scatter and adjusted-estimate bar plots for H2."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h2-20261025-mpl')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h1_h4_literature_methods_20261025/h2')
GROUPS=['ND','Aab+','T1D']
COLORS={'Pooled':'#444444','ND':'#377eb8','Aab+':'#df9b22','T1D':'#bf3434'}
d=pd.read_csv(OUT/'h2_analysis_objects.csv')
r=pd.read_csv(OUT/'h2_model_results.csv')
d['core_immune_present']=d.n_immune_core.gt(0)
donors=d.groupby(['case_id','donor_status','core_immune_present'],observed=True).agg(
    objects=('islet_id','size'),median_coverage=('coliv_coverage_2um','median'),
    mean_coverage=('coliv_coverage_2um','mean')).reset_index()
donors.to_csv(OUT/'h2_donor_coverage.csv',index=False)
fig=plt.figure(figsize=(12,8.5))
gs=fig.add_gridspec(2,3,height_ratios=[1,1.15],hspace=.48,wspace=.25)
for j,group in enumerate(GROUPS):
    ax=fig.add_subplot(gs[0,j]); z=donors.loc[donors.donor_status.eq(group)]
    for case, sub in z.groupby('case_id'):
        sub=sub.sort_values('core_immune_present'); x=sub.core_immune_present.astype(int).to_numpy()
        ax.plot(x,sub.median_coverage*100,'o-',color=COLORS[group],alpha=.75,lw=1.2,ms=5)
    ax.set_xticks([0,1],['Around only','Inside core'])
    ax.set_xlim(-.25,1.25);ax.set_ylim(0,20)
    ax.set_title(f'{group}: {z.case_id.nunique()} donors')
    if j==0:ax.set_ylabel('Median boundary ColIV coverage (%)')
    ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,:]);z=r.loc[r.width_um.eq(2)&r.outcome.eq('core_density')&r.model.eq('separate_group')].set_index('group').loc[['Pooled',*GROUPS]]
x=np.arange(len(z))
ax.bar(x,z.estimate,color=[COLORS[g] for g in z.index],alpha=.75,width=.6)
ax.errorbar(x,z.estimate,yerr=np.vstack([z.estimate-z.lower,z.upper-z.estimate]),fmt='none',ecolor='#222',capsize=5)
ax.axhline(0,color='#555',lw=.8);ax.set_xticks(x,z.index)
ax.set_ylabel('Adjusted change in core immune cells/mm²\nper 10 percentage points more ColIV coverage')
ax.set_title('Donor mixed-effects model; bars show estimates and 95% confidence intervals')
ax.spines[['top','right']].set_visible(False)
fig.suptitle('H2: Does greater boundary ColIV coverage predict immune exclusion?',fontsize=15,y=.99)
fig.savefig(OUT/'h2_coliv_results.png',dpi=200,bbox_inches='tight')
fig.savefig(OUT/'h2_coliv_results.pdf',bbox_inches='tight');plt.close(fig)

validity=dict(objects=len(d),donors=int(d.case_id.nunique()),missing_primary=int(d.coliv_coverage_2um.isna().sum()),
              missing_sensitivity=int(d.coliv_coverage_5um.isna().sum()),
              duplicate_islet_ids=int(d.islet_id.duplicated().sum()),
              all_coverage_bounded=bool(d[['coliv_coverage_2um','coliv_coverage_5um']].stack().between(0,1).all()),
              all_nonvascular_pixels_lte_band_pixels=bool((d.nonvascular_boundary_pixels_2um<=d.boundary_pixels_2um).all()),
              image_thresholds=20,source_images_modified=False)
assert validity['objects']==2937 and validity['donors']==20 and not validity['duplicate_islet_ids']
assert validity['all_coverage_bounded'] and validity['all_nonvascular_pixels_lte_band_pixels']
(OUT/'validation.json').write_text(json.dumps(validity,indent=2))
lines=['# H2: boundary ColIV coverage and core immune density','',
'All 2,937 corrected immune-exposed endocrine sections from 20 donors were measured directly in the existing original ColIV images. The current islet/cluster outlines and current CD31 structural objects were retained; the app and data pipeline were not rebuilt.','',
'Otsu thresholds were selected independently for each donor from regularly sampled Tissue pixels, without using immune counts or disease results. Coverage is the ColIV-positive area fraction in a symmetric ±2 µm boundary band after intersecting Tissue and excluding current CD31 object pixels. A ±5 µm band is the sensitivity analysis. No morphological closing or September ColIV masks were used.','',
'The primary linear mixed-effects model uses core immune density, with donor random intercepts and adjustment for peri immune density, core area, perimeter, distance from the Tissue edge, CD31 area fraction, endocrine object class and (in the pooled model) disease group. The parallel model on log1p density checks the influence of large density values. The group-interaction model shares adjustment coefficients, while the separate-group models estimate those coefficients within each group. Reported p-values are unadjusted here so that the combined H1–H4 report can apply its chosen multiple-comparison family.','',
'## Primary estimates','',
'Effect is change in core immune cells/mm² for 10 percentage points more boundary ColIV coverage. Negative estimates support exclusion.','',
'| Group | Estimate | 95% CI | p |','|---|---:|---:|---:|']
for row in z.itertuples():
    lines.append(f'| {row.Index} | {row.estimate:.1f} | {row.lower:.1f} to {row.upper:.1f} | {row.p:.4g} |')
lines+=['',
'H2 is not statistically supported by the primary model. T1D has the hypothesized negative direction, but its confidence interval includes no association. The T1D direction remains negative in the ±5 µm band and log1p-density analyses; neither establishes an association. ND and Aab+ point in the opposite direction. The nominally significant ND result in the ±5 µm log1p-density analysis is positive, so it does not support exclusion and does not alter the primary conclusion.','',
'## Method sources','',
'- `literature/Heide.pdf`, PDF p30: thresholded positive area divided by anatomical area, with geometric intersections; p32–33: donor mixed models, scatter and bar plots.',
'- `literature/Richardson2023.pdf`, PDF p3: image/channel thresholds for object area and ColIV overlap.',
'- `~/Documents/papers/Nature Comm 2026 spatial transcriptomics for IBD research.pdf`, PDF p15: automated Otsu thresholding.',
'- `literature/Bodenmiller2025.pdf`, PDF p28: linear mixed-effects models for cell density using donor random intercepts and lmerTest.','',
'The narrow-boundary application is an adaptation of published thresholded-area measurement to H2. It measures boundary staining coverage; it does not completely reconstruct basement-membrane topology. The paired donor plot is descriptive, while the adjusted model accounts for peri immune exposure and geometry. A targeted follow-up can quantify spatial gaps or repeat the measurement with another thresholding setting.','',
'Tables: `h2_model_results.csv`, `h2_analysis_objects.csv`, `h2_donor_coverage.csv`. Reproducibility: `source_provenance.json`, per-donor threshold records, saved R model objects, and `validation.json`.']
(OUT/'README.md').write_text('\n'.join(lines)+'\n')
print(z[['estimate','lower','upper','p']].to_string())
print(json.dumps(validity))
