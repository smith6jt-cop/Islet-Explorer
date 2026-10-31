#!/usr/bin/env python3
"""Descriptive H7 donor measurements only; no tests or null inference."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h7-descriptive-20261026-mpl')
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path('output/h5_h10_original_plan_20261026/h7')
d=pd.read_parquet(OUT/'object_neurovascular_features.parquet')
d['beta_mean_available']=d.ins_mean_beta_core.notna()
d['matched']=d.beta_mean_available & d.nerve_near_cd31_20um.notna()
rows=[]
for case,g in d.groupby('case_id'):
    z=g.loc[g.matched]
    rows.append(dict(case_id=case,donor_status=g.donor_status.iloc[0],h5_objects=len(g),
      beta_positive_objects=int(g.beta_mean_available.sum()),matched_objects=len(z),
      beta_objects_without_nerve_area=int((g.beta_mean_available & g.neighborhood_complete_20 & g.nerve_area_um2.fillna(0).le(0)).sum()),
      beta_objects_incomplete_20um_support=int((g.beta_mean_available & ~g.neighborhood_complete_20).sum()),
      mean_proximity20_fraction=z.nerve_near_cd31_20um.mean(),median_proximity20_fraction=z.nerve_near_cd31_20um.median(),
      mean_proximity20_percent=100*z.nerve_near_cd31_20um.mean(),
      mean_nerve_area_fraction=z.nerve_area_fraction.mean(),mean_cd31_area_fraction=z.cd31_area_fraction.mean(),
      mean_gap43_supported_nerve_fraction=z.gap43_supported_nerve_fraction.mean()))
donors=pd.DataFrame(rows)
assert len(donors)==20 and donors.matched_objects.gt(0).all()
assert (donors.beta_positive_objects==donors.matched_objects+donors.beta_objects_without_nerve_area+donors.beta_objects_incomplete_20um_support).all()
donors.to_csv(OUT/'donor_observed_proximity.csv',index=False)
d.loc[d.matched,['islet_id','case_id','donor_status','n_beta_core','ins_mean_beta_core','core_area_um2','core_perimeter_um','n_immune_nearby',
 'peri_assessed_area_um2','nerve_area_um2','cd31_area_um2','nerve_area_fraction','cd31_area_fraction','nerve_near_cd31_20um',
 'nerve_near_cd31_10um','nerve_near_cd31_5um','nerve_vessel_overlap','gap43_supported_nerve_fraction']].to_csv(OUT/'matched_object_plot_data.csv',index=False)
groups=['ND','Aab+','T1D'];colors=['#4477AA','#EE7733','#228833']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(6,4.8))
values=[donors.loc[donors.donor_status.eq(g)].sort_values('case_id').mean_proximity20_percent.to_numpy() for g in groups]
bp=ax.boxplot(values,positions=[1,2,3],widths=.5,showfliers=False,patch_artist=True,
 medianprops={'color':'#222','linewidth':1.2},boxprops={'linewidth':1},whiskerprops={'linewidth':1},capprops={'linewidth':1})
cols={}
for i,(g,v,c) in enumerate(zip(groups,values,colors),1):
    bp['boxes'][i-1].set(facecolor=c,edgecolor=c,alpha=.18)
    ax.scatter(i+np.random.default_rng(20261026+i).uniform(-.10,.10,len(v)),v,s=38,color=c,edgecolor='white',linewidth=.4,zorder=3)
    z=donors.loc[donors.donor_status.eq(g)].sort_values('case_id').reset_index(drop=True)
    cols[g+'_donor_id']=z.case_id;cols[g+'_value']=z.mean_proximity20_percent
ax.set_xticks([1,2,3],['ND','AAb+','T1D']);ax.set_xlabel('Disease state')
ax.set_ylabel('Mean nerve area within 20 µm of CD31 (%)')
ax.set_ylim(0,103);ax.set_title('Observed nerve–vessel proximity\nnear β-cell-containing objects',fontsize=12)
fig.tight_layout();fig.savefig(OUT/'observed_proximity_by_disease.png',dpi=300);fig.savefig(OUT/'observed_proximity_by_disease.pdf');plt.close(fig)
(OUT/'prism_csv').mkdir(exist_ok=True);pd.DataFrame(cols).to_csv(OUT/'prism_csv/observed_proximity.csv',index=False)
pd.DataFrame([dict(figure='observed_proximity_by_disease.png',panel='observed_proximity',plot_data_csv=str(OUT/'prism_csv/observed_proximity.csv'),
 long_data_csv='donor_observed_proximity.csv',filter='all20donors',donor_id='case_id',x='donor_status',y='mean_proximity20_percent',
 statistics_csv='',stars_column='',biological_replicate='donor',inference='none_descriptive_only')]).to_csv(OUT/'plot_manifest.csv',index=False)
counts=donors.groupby('donor_status').matched_objects.sum().to_dict()
legend=f'''Observed nerve–vessel proximity near Beta-containing endocrine objects. Each point is one donor (ND n=7, AAb+ n=7, T1D n=6); the value is the mean of that donor's object-level proportions of peri B3TUBB-positive area lying within 20 µm of a CD31-supported structure. Boxes show the median and interquartile range, with whiskers to values within 1.5 interquartile ranges. Matched object counts are ND {counts['ND']}, AAb+ {counts['Aab+']}, T1D {counts['T1D']}; counts for each donor are supplied in donor_observed_proximity.csv. No statistical tests or spatial-null comparisons were performed, and no significance brackets are shown.

The nerve-marker mask was computed directly from original B3TUBB images using per-image Otsu intensity thresholds and exclusion of connected components smaller than 2 µm². Current endocrine cores were excluded from the mask; INS was not used. Measurements used the 20 µm peri region intersected with Tissue and excluding all endocrine cores. CD31 structures are the existing current structural objects. GAP43 overlap is reported as corroboration, without being a required positive gate. Proximity is undefined where no B3TUBB-positive peri area was measured; those objects are counted explicitly rather than assigned zero proximity. Objects lacking complete 20 µm acquisition support are also counted separately. The displayed set requires a measured within-Beta mean INS from H5, but the displayed values do not test its association with proximity. Each donor contributes one plotted value regardless of the number of objects.

These are observed marker-defined proximity measurements. This plot reports observed proximity only; no association test or spatial-null calculation is used in this figure. The user-approved adjusted observed-proximity regression is reported separately in inference_figure_legends.md. Source images, thresholds, UUIDs, geometry hashes and full object measurements are retained. No app pipeline was rebuilt.
'''
(OUT/'figure_legend.md').write_text(legend)
(OUT/'README.md').write_text('# Original H7: prepared observed neurovascular measurements\n\n'+legend)
summary=dict(h5_objects=len(d),proximity_measured_objects=int(d.nerve_near_cd31_20um.notna().sum()),
 beta_proximity_matched_objects=int(d.matched.sum()),matched_donors=len(donors),matched_by_status=counts,
 ins_channel_used_for_nerve_mask=False,statistical_tests_performed=False,spatial_null_performed=False,
 proximity_values_bounded=bool(d.nerve_near_cd31_20um.dropna().between(0,1).all()))
(OUT/'descriptive_validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
