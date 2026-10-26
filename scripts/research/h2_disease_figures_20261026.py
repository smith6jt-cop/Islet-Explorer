#!/usr/bin/env python3
"""Disease-state donor comparison plot and a separate complete figure legend."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h2-disease-20261026-mpl')
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h1_h4_disease_comparisons_20261026/h2')
INPUT=Path('output/h1_h4_literature_methods_20261025/h2/h2_analysis_objects.csv')
GROUPS=['ND','Aab+','T1D']
COLORS=['#4477AA','#EE7733','#228833']
d=pd.read_csv(OUT/'primary_donor_coefficients.csv')
pair=pd.read_csv(OUT/'primary_disease_pairwise_tests.csv')
overall=pd.read_csv(OUT/'primary_disease_overall_test.csv').iloc[0]
summary=pd.read_csv(OUT/'primary_disease_summary.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(5.8,4.6))
values=[d.loc[d.donor_status.eq(g)].sort_values('case_id').standardized_beta.to_numpy() for g in GROUPS]
bp=ax.boxplot(values,positions=[1,2,3],widths=.5,showfliers=False,patch_artist=True,
              medianprops={'color':'#222222','linewidth':1.4},
              boxprops={'linewidth':1.1},whiskerprops={'linewidth':1.1},capprops={'linewidth':1.1})
for j,(v,c) in enumerate(zip(values,COLORS),1):
    bp['boxes'][j-1].set(facecolor=c,alpha=.18,edgecolor=c)
    jitter=np.random.default_rng(20261026+j).uniform(-.115,.115,len(v))
    ax.scatter(j+jitter,v,color=c,s=43,edgecolors='white',linewidth=.45,zorder=3)
ax.axhline(0,color='#888888',ls='--',lw=.8,zorder=0)
ax.set_xticks([1,2,3],['ND','Aab+','T1D'])
ax.set_xlabel('Disease state')
ax.set_ylabel('Standardized regression coefficient (β)')
ax.set_title('Boundary ColIV coverage and core immune density',fontsize=11,pad=15)
lo,hi=float(d.standardized_beta.min()),float(d.standardized_beta.max())
span=max(hi-lo,.1);y=hi+.13*span
for row in pair.itertuples():
    if row.p_BH>=.05:
        continue
    a,b=GROUPS.index(row.group1)+1,GROUPS.index(row.group2)+1
    stars='***' if row.p_BH<.001 else '**' if row.p_BH<.01 else '*'
    ax.plot([a,a,b,b],[y,y+.03*span,y+.03*span,y],color='#222222',lw=1)
    ax.text((a+b)/2,y+.035*span,stars,ha='center',va='bottom')
    y+=.15*span
ax.set_ylim(lo-.12*span,max(hi+.16*span,y+.03*span))
fig.tight_layout()
fig.savefig(OUT/'coliv_disease_comparison.png',dpi=300)
fig.savefig(OUT/'coliv_disease_comparison.pdf')
plt.close(fig)

legend=(
'Boundary ColIV coverage and core immune density across disease states. Each point represents one donor '
'(ND, n=7; AAb+, n=7; T1D, n=6); boxes show the median and interquartile range, and whiskers extend to '
'the most extreme values within 1.5 interquartile ranges. The vertical axis is the adjusted standardized '
'linear-regression coefficient (β), not a correlation coefficient. Negative β indicates lower core immune '
'density with greater boundary ColIV coverage. A separate regression was fitted for each donor using all '
'eligible endocrine objects, with adjustment for peri immune density, core area, perimeter, distance from '
'the Tissue edge, CD31 structural area fraction and object class. Coverage and continuous adjustment '
'variables were standardized within donor; core immune density was standardized across the complete '
'dataset. Boundary coverage was the existing ColIV-positive fraction in a ±2 µm band, intersected with '
'Tissue and excluding CD31 object pixels. In total, 2,937 endocrine objects from 20 donors were included. '
'The primary donor random-intercept/random-slope model had a singular covariance; the independent '
'intercept/slope model estimated essentially zero slope variance. Therefore, disease groups were compared '
'using the 20 donor coefficients with a Kruskal–Wallis test, followed by all three two-sided exact Wilcoxon '
'rank-sum tests with Benjamini–Hochberg adjustment. This fallback was triggered by covariance estimation, '
'not statistical significance. '
f'Kruskal–Wallis H={overall.statistic:.3f}, df=2, p={overall.p:.3f}. '
f'Adjusted pairwise p-values: ND versus AAb+, {pair.iloc[0].p_BH:.3f}; ND versus T1D, {pair.iloc[1].p_BH:.3f}; '
f'AAb+ versus T1D, {pair.iloc[2].p_BH:.3f}. No pairwise comparison is significant, so no significance brackets '
'are displayed. The measurement describes boundary staining coverage rather than complete membrane topology.'
)
(OUT/'figure_legend.md').write_text(legend+'\n')
lines=['# ColIV association: direct disease-state comparisons','',
'The association between boundary ColIV coverage and core immune density did not differ detectably among ND, AAb+ and T1D donors.',
'',f'The primary comparison treats the donor as the biological replicate: Kruskal–Wallis H = {overall.statistic:.3f}, df = 2, p = {overall.p:.3f}.',
'','| Comparison | BH-adjusted p |','|---|---:|']
for row in pair.itertuples():lines.append(f'| {row.group1} versus {row.group2} | {row.p_BH:.3f} |')
lines+=['','Median adjusted standardized β: '+', '.join(f'{row.donor_status} {row.median:.3f}' for row in summary.itertuples())+'.',
'','The negative median coefficient in T1D is not evidence that the relationship differs from ND or AAb+. No pooled or separate within-state significance test is used for the primary conclusion.',
'','The intended random-intercept/random-slope model had a singular covariance for raw core immune density; the independent-components alternative estimated essentially zero donor slope variance. Both fits are retained. The primary fallback compares one adjusted regression coefficient per donor with the established donor-level Kruskal–Wallis and Wilcoxon tests. Regression coefficients are not labeled as correlations.',
'','The log1p-density random-slope sensitivity was nonsingular. Its disease-by-coverage interaction was also not significant (F(2, 10.95) = 1.536, p = 0.258; all BH-adjusted pairwise p ≥ 0.335). A wider-band sensitivity agreed. These sensitivity analyses do not replace the primary raw-density endpoint.',
'','Local method precedents: `literature/Barlow2025.pdf`, Methods p22, donor random intercept/slope and within-donor predictor standardization; `literature/Bodenmiller2025.pdf`, donor-level Wilcoxon comparisons and FDR; `literature/Richardson2023.pdf`, Kruskal–Wallis group comparisons. The previously established thresholded-area measurement is reused unchanged.',
'','Plot definitions and all displayed statistics are in `figure_legend.md`. Input source, counts, identifiability checks and checksums are recorded in `validation.json`.']
(OUT/'README.md').write_text('\n'.join(lines)+'\n')
validation=dict(donors=len(d),objects=int(d.objects.sum()),disease_counts=d.donor_status.value_counts().to_dict(),
  all_donor_regressions_full_rank=bool((d.design_rank==d.design_columns).all()),
  all_donor_coefficients_finite=bool(np.isfinite(d.standardized_beta).all()),
  primary_method='Kruskal-Wallis and exact pairwise Wilcoxon rank-sum on one adjusted coefficient per donor; BH over all three pairwise comparisons',
  source=str(INPUT),source_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),new_image_extraction=False)
assert validation['donors']==20 and validation['objects']==2937 and validation['all_donor_regressions_full_rank'] and validation['all_donor_coefficients_finite']
# Independent numerical reproduction of the R donor fits and rank tests.
source=pd.read_csv(INPUT)
y=source.n_immune_core*1e6/source.core_area_um2
ys=(y-y.mean())/y.std(ddof=1)
raw=np.column_stack([source.coliv_coverage_2um,
    np.log1p(source.n_immune_peri*1e6/source.peri_area_um2),np.log(source.core_area_um2),
    np.log(source.core_perimeter_um),np.log1p(source.core_distance_to_tissue_edge_um.clip(lower=0)),
    source.capillary_core_area_fraction])
errors=[]
for case,z in source.groupby('case_id'):
    v=raw[z.index]; xz=(v-v.mean(axis=0))/v.std(axis=0,ddof=1)
    design=np.column_stack([np.ones(len(z)),xz,z.object_class.eq('Islet')])
    beta=np.linalg.lstsq(design,ys[z.index],rcond=None)[0][1]
    errors.append(abs(beta-d.loc[d.case_id.eq(case),'standardized_beta'].iloc[0]))
independent_kw=kruskal(*values)
pair_errors=[]
for row in pair.itertuples():
    pv=mannwhitneyu(d.loc[d.donor_status.eq(row.group1),'standardized_beta'],
                   d.loc[d.donor_status.eq(row.group2),'standardized_beta'],method='exact').pvalue
    pair_errors.append(abs(pv-row.p))
assert max(errors)<1e-10 and abs(independent_kw.pvalue-overall.p)<1e-10 and max(pair_errors)<1e-12
validation.update(independent_numpy_max_beta_difference=float(max(errors)),
                  independent_scipy_kruskal_p=float(independent_kw.pvalue),
                  independent_scipy_max_pairwise_p_difference=float(max(pair_errors)))
(OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print(json.dumps(validation,indent=2))
