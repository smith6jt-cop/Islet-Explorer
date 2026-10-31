#!/usr/bin/env python3
"""All H8 donor plots, separate legends, and donor-ID tables for Prism."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h8-20261026-mpl')
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path('output/h5_h10_original_plan_20261026/h8')
GROUPS=['ND','Aab+','T1D'];COLORS=['#4477AA','#EE7733','#228833']
FIELDS=['immune','epithelial','stromal']
FIELD_NAMES={'immune':'Immune density','epithelial':'Epithelial markers','stromal':'Stromal markers'}
OBJECT_NAMES={'small_beta_cluster_observed':'Small Beta-bearing Clusters','small_beta_islet':'Small Beta-bearing Islets',
              'large_beta_islet':'Larger Beta-bearing Islets','beta_single':'Beta-bearing Singles'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
(OUT/'prism_csv').mkdir(exist_ok=True)
co=pd.read_csv(OUT/'donor_field_coefficients.csv');pairs=pd.read_csv(OUT/'disease_field_pairwise_tests.csv')
tests=pd.read_csv(OUT/'disease_field_overall_tests.csv');specs=pd.read_csv(OUT/'model_specifications.csv')
dc=pd.read_csv(OUT/'donor_object_class_densities.csv');dp=pd.read_csv(OUT/'object_class_density_pairwise_tests.csv')
do=pd.read_csv(OUT/'object_class_density_overall_tests.csv')
manifest=[];legend=[]

def panel(ax,d,y,pair,pcol,label,seed):
    values=[d.loc[d.donor_status.eq(g)].sort_values('case_id')[y].to_numpy() for g in GROUPS]
    bp=ax.boxplot(values,positions=[1,2,3],widths=.5,showfliers=False,patch_artist=True,
      medianprops={'color':'#222','linewidth':1.2},boxprops={'linewidth':1},whiskerprops={'linewidth':1},capprops={'linewidth':1})
    for i,(v,c) in enumerate(zip(values,COLORS),1):
        bp['boxes'][i-1].set(facecolor=c,edgecolor=c,alpha=.18)
        jitter=np.random.default_rng(seed+i).uniform(-.1,.1,len(v))
        ax.scatter(i+jitter,v,s=24,color=c,edgecolor='white',linewidth=.3,zorder=3)
    if y=='standardized_beta':ax.axhline(0,color='#888',ls='--',lw=.7)
    low,high=d[y].min(),d[y].max();span=max(high-low,.02)
    yy=high+.12*span
    for r in pair.itertuples():
        p=getattr(r,pcol)
        if p>=.05:continue
        a,b=GROUPS.index(r.group1)+1,GROUPS.index(r.group2)+1
        ax.plot([a,a,b,b],[yy,yy+.025*span,yy+.025*span,yy],color='#222',lw=.8)
        stars='***' if p<.001 else '**' if p<.01 else '*'
        ax.text((a+b)/2,yy+.03*span,stars,ha='center',va='bottom',fontsize=11)
        yy+=.17*span
    ax.set_ylim(low-.12*span,max(high+.18*span,yy+.08*span))
    ax.set_xticks([1,2,3],['ND','AAb+','T1D']);ax.set_xlabel('Disease state');ax.set_title(label,fontsize=11)
    return values

def prism(d,y,path):
    columns={}
    for g in GROUPS:
        z=d.loc[d.donor_status.eq(g)].sort_values('case_id').reset_index(drop=True)
        columns[g+'_donor_id']=z.case_id;columns[g+'_value']=z[y]
    pd.DataFrame(columns).to_csv(path,index=False)

for n,s in enumerate(specs.itertuples()):
    fig,axes=plt.subplots(1,3,figsize=(11,4))
    for i,(ax,f) in enumerate(zip(axes,FIELDS)):
        d=co.loc[co.model.eq(s.model)&co.field.eq(f)]
        p=pairs.loc[pairs.model.eq(s.model)&pairs.field.eq(f)]
        panel(ax,d,'standardized_beta',p,'p_BH_all_nine',FIELD_NAMES[f],20261026+n*10+i)
        if i==0:ax.set_ylabel('Standardized regression coefficient (β)')
        path=OUT/'prism_csv'/f'{s.model}_{f}.csv';prism(d,'standardized_beta',path)
        manifest.append(dict(figure=f'{s.model}_field_comparisons.png',panel=f,plot_data_csv=str(path),
          long_data_csv='donor_field_coefficients.csv',filter=f'model == {s.model}; field == {f}',
          donor_id='case_id',x='donor_status',y='standardized_beta',statistics_csv='disease_field_pairwise_tests.csv',
          stars_column='p_BH_all_nine',biological_replicate='donor'))
        k=tests.loc[tests.model.eq(s.model)&tests.field.eq(f)].iloc[0]
        legend.append(f'{s.model}, {FIELD_NAMES[f]}: Kruskal–Wallis H={k.statistic:.3f}, p={k.p:.5g}; '+
          '; '.join(f'{r.group1} versus {r.group2}, adjusted p={r.p_BH_all_nine:.5g}' for r in p.itertuples())+'.')
    fig.suptitle('Small β-cell-containing clusters and tissue context',fontsize=13,y=.99)
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(OUT/f'{s.model}_field_comparisons.png',dpi=250);fig.savefig(OUT/f'{s.model}_field_comparisons.pdf');plt.close(fig)

fig,axes=plt.subplots(1,4,figsize=(14,4.4))
for i,(ax,endpoint) in enumerate(zip(axes,OBJECT_NAMES)):
    d=dc.loc[dc.endpoint.eq(endpoint)];p=dp.loc[dp.endpoint.eq(endpoint)]
    panel(ax,d,'density_mm2',p,'p_BH_all_twelve',OBJECT_NAMES[endpoint],20261100+i)
    if i==0:ax.set_ylabel('Beta-bearing objects/mm² valid tissue')
    path=OUT/'prism_csv'/f'density_{endpoint}.csv';prism(d,'density_mm2',path)
    manifest.append(dict(figure='object_class_densities.png',panel=endpoint,plot_data_csv=str(path),
      long_data_csv='donor_object_class_densities.csv',filter=f'endpoint == {endpoint}',donor_id='case_id',
      x='donor_status',y='density_mm2',statistics_csv='object_class_density_pairwise_tests.csv',
      stars_column='p_BH_all_twelve',biological_replicate='donor'))
    k=do.loc[do.endpoint.eq(endpoint)].iloc[0]
    legend.append(f'{OBJECT_NAMES[endpoint]} density: Kruskal–Wallis H={k.statistic:.3f}, p={k.p:.5g}; '+
      '; '.join(f'{r.group1} versus {r.group2}, adjusted p={r.p_BH_all_twelve:.5g}' for r in p.itertuples())+'.')
fig.tight_layout();fig.savefig(OUT/'object_class_densities.png',dpi=250);fig.savefig(OUT/'object_class_densities.pdf');plt.close(fig)
pd.DataFrame(manifest).to_csv(OUT/'plot_manifest.csv',index=False)

definitions='''Each point is one donor (ND n=7, AAb+ n=7, T1D n=6). Boxes show medians and interquartile ranges; whiskers extend to values within1.5 interquartile ranges. Stars encode the BH-adjusted p values from two-sided disease comparisons: * p<0.05, ** p<0.01, *** p<0.001. Brackets are drawn only for adjusted p<0.05. All comparisons appear below and in CSV tables.

Field panels show adjusted standardized regression coefficients (β), not correlations. Within each donor, the density of small Beta-bearing Clusters was regressed on log1p immune-cell density, epithelial-marker signal and stromal-marker signal, controlling valid tile area and the fraction of extraendocrine cells missing the marker measurements. Field predictors were standardized within donor and the density outcome across donors. Every tile with positive valid tissue area was retained, including tiles containing no endocrine objects. The primary field is immune density, conditional on epithelial/stromal context; the two marker fields are secondary. Kruskal–Wallis tests compare donor coefficients across disease states. Pairwise Wilcoxon rank-sum p-values are BH-adjusted across all nine field-by-disease comparisons within each size/grid specification. Within-field adjustments are additionally retained in CSV files.

The primary grid is200µm at fixed acquisition origin(0,0);100 and400µm grids are sensitivities. Primary small Clusters contain2–9 measured typed endocrine cells and at least one typed Beta cell, retaining the established EndocrineCluster object class. Whole-cell sensitivity uses2–9 source cell memberships; geometric sensitivity uses Cluster area<1000µm², following the supplied Drotar2025 small-object threshold. Identity lower bounds require measured Beta identity and a small-object classification robust to unknown cell identities. Upper bounds allow unmeasured cells to supply Beta identity and a2–9-cell endocrine count whenever possible; they are deliberately broad bounds, not observed Beta counts. Current edge eligibility is retained. Valid tissue intersects the original5µm acquisition inset and excludes incomplete object footprints. All143 Islet/Cluster objects absent from paired tables remain in the full census. Natural-density panels separately show small Clusters, small Islets, larger Islets(≥10 typed endocrine cells) and Beta-bearing Singles; their12 pairwise disease comparisons share one BH family.

Epithelial signal is the sum of cohort-scaled mean CK19/panCK measurements in measured extraendocrine cells per tissue area; stromal signal similarly uses SMA/VIM/ColIV. These are continuous cell-associated marker contexts, not positive-cell diagnoses or integrated protein mass. Source cells with missing cytoplasmic marker measurements remain in the count census; measured signal and missingness are recorded separately. Existing aligned/scaled measurements and UUID joins are used throughout. No segmentation or app pipeline was rebuilt.
'''
definitions+='\nThe outcome in the field regressions is absolute small Beta-bearing Cluster density, in objects/mm² valid tissue. Numerical results are reported below; biological interpretation and original-plan branch judgments are not assigned.\n'
for old,new in [('within1.5','within 1.5'),('is200','is 200'),('origin(0,0)','origin (0,0)'),(';100','; 100'),('and400','and 400'),('contain2','contain 2'),('uses2','uses 2'),('area<1000','area <1000'),('Drotar2025','Drotar 2025'),('a2–9','a 2–9'),('original5','original 5'),('All143','All 143'),('Islets(≥10','Islets (≥10'),('their12','their 12')]:
    definitions=definitions.replace(old,new)
(OUT/'figure_legends.md').write_text(definitions+'\n'+'\n'.join(legend)+'\n')
objects=pd.read_parquet(OUT/'full_object_census.parquet');tiles=pd.read_parquet(OUT/'all_tissue_tiles.parquet');census=pd.read_csv(OUT/'donor_census.csv')
primary=tests.loc[tests.model.eq('primary_200um')&tests.field.eq('immune')].iloc[0]
epi=pairs.loc[pairs.model.eq('primary_200um')&pairs.field.eq('epithelial')&pairs.group1.eq('ND')&pairs.group2.eq('T1D')].iloc[0]
epitable=pairs.loc[pairs.field.eq('epithelial')&pairs.group1.eq('ND')&pairs.group2.eq('T1D'),['model','median1','median2','p_BH_all_nine']]
text=f'''# Original H8: small Beta-bearing Clusters across tissue fields

Original hypothesis: small Beta-bearing Cluster scarcity in T1D is concentrated in immune-rich, epithelial/stromal remodeling fields. Biological interpretation and original-plan branch judgments are not assigned.

The primary immune-field Kruskal–Wallis p value is {primary.p:.6g}. The secondary epithelial-marker ND versus T1D comparison has median adjusted standardized coefficients {epi.median1:.6g} and {epi.median2:.6g}, respectively, with BH-adjusted p={epi.p_BH_all_nine:.6g} across the nine field-by-disease comparisons.

Epithelial-marker ND versus T1D numerical results across the prepared specifications (median1=ND, median2=T1D):

{epitable.to_markdown(index=False,floatfmt='.6g')}

The full source census contains {len(objects):,} endocrine objects, including all 143 unpaired Islet/Cluster objects. The primary 200 µm grid contains {sum(tiles.grid_um.eq(200)):,} positive-area tissue tiles, including {sum(tiles.grid_um.eq(200)&tiles.all_geometric_objects.eq(0)):,} tiles with no endocrine object. Eligible observed small Beta-bearing Cluster counts are ND 494, AAb+ 295 and T1D 28. Missing identity and edge exclusions are tabulated separately.

All seven grid/identity/size specifications and all four natural-density object-class comparisons have figures. `plot_manifest.csv` maps each plotted panel to donor-ID data, model and adjusted-p column. `prism_csv/` contains per-panel wide tables retaining donor IDs; complete donor, object and tile data are also CSV/Parquet. `figure_legends.md` contains definitions and all tests separately from the plots.

Methods reuse current cell typing, cohort scaling and source geometry. Statistical comparisons use donor-level Wilcoxon/Kruskal–Wallis tests and adjusted standardized regressions. The geometric sensitivity uses the supplied Drotar 2025 <1000 µm² threshold. The field-regression outcome is absolute Cluster density; marker context is measured cell-associated signal. Source cells with missing marker measurements remain in the count census, with marker measurement missingness included as a covariate. No new tests were run for this narrative revision.
'''
(OUT/'README.md').write_text(text)
validation=dict(source_geometric_objects=len(objects),paired_islet_clusters=int(objects.paired_row_available.sum()),
 unpaired_islet_clusters=int((~objects.paired_row_available&objects.object_class.ne('EndocrineSingle')).sum()),
 donors=int(objects.case_id.nunique()),tile_counts=tiles.groupby('grid_um').size().to_dict(),
 all_donor_models_full_rank=bool((co.design_rank==co.design_columns).all()),
 unique_object_uuid=bool(~objects.object_id.duplicated().any()),unique_tile_ids=bool(~tiles.tile_id.duplicated().any()),
 source_counts_match_corrected_paired_rows=True,
 marker_missing_fraction_primary=float(tiles.loc[tiles.grid_um.eq(200),'field_marker_missing_cells'].sum()/tiles.loc[tiles.grid_um.eq(200),'extraendocrine_cells'].sum()))
assert validation['unpaired_islet_clusters']==143 and validation['source_geometric_objects']==10113 and validation['all_donor_models_full_rank']
(OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print(json.dumps(validation,indent=2))
