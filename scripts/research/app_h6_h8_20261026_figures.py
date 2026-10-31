#!/usr/bin/env python3
"""Non-INS H6–H8 statistics and linked existing-app state/tissue displays."""
from pathlib import Path
import os,json,hashlib,shutil
os.environ.setdefault('MPLCONFIGDIR','/tmp/app_hypotheses_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
import app_hypotheses_shared_20261026 as S
ROOT=S.ROOT;OUT=ROOT/'output/app_hypotheses_20261026/h6_h8';OLD=ROOT/'output/h5_h10_original_plan_20261026'
DATA=OUT/'data';FIG=OUT/'figures';MODELS=OUT/'models';GROUPS=S.GROUPS;COLORS=S.COLORS
records=[]

def bh(values):
 a=np.asarray(values,float);result=np.full(len(a),np.nan);ix=np.flatnonzero(np.isfinite(a));order=ix[np.argsort(a[ix])]
 if len(order):result[order]=np.minimum(1,np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(order)+1))[::-1])[::-1])
 return result

def add(stem,legend,data):
 (FIG/f'{stem}_legend.txt').write_text(legend+'\n')
 records.append({'stem':stem,'png':str(FIG/f'{stem}.png'),'pdf':str(FIG/f'{stem}.pdf'),'svg':str(FIG/f'{stem}.svg'),'legend':str(FIG/f'{stem}_legend.txt'),'data_csv':[str(Path(p)) for p in data]})

def save(fig,stem,legend,data):S.savefig(fig,FIG,stem);add(stem,legend,data)

def states(d,features,stem,note):
 result=S.state_views(d,features,OUT,stem,note)
 for f in result['figures']:
  data=[DATA/f'{stem}_state_objects.csv',DATA/f'{stem}_pseudotime_profiles.csv',DATA/f'{stem}_donor_pseudotime.csv',DATA/f'{stem}_donor_cluster_values.csv']
  records.append({'stem':f,'png':str(FIG/f'{f}.png'),'pdf':str(FIG/f'{f}.pdf'),'svg':str(FIG/f'{f}.svg'),'legend':str(FIG/f'{f}_legend.txt'),'data_csv':[str(p) for p in data]})

def group_slopes(d,stem,title,legend,data,pairs=None):
 fig,ax=plt.subplots(figsize=(6.5,5.3),layout='constrained')
 for i,g in enumerate(GROUPS):
  r=d[d.donor_status.eq(g)].iloc[0]
  ax.errorbar(i,r.estimate,yerr=np.asarray([[r.estimate-r.lower],[r.upper-r.estimate]]),fmt='o',color=COLORS[g],capsize=5,ms=8,lw=1.7)
 ax.axhline(0,color='#777',ls='--',lw=1);ax.set(xticks=range(3),xticklabels=GROUPS,ylabel='Adjusted standardized regression coefficient (β)',title=title)
 if pairs is not None and len(pairs):
  lo=min(d.lower.min(),0);hi=max(d.upper.max(),0);step=max((hi-lo)*.15,.025);y=hi+step
  for r in pairs.itertuples():
   if not np.isfinite(r.p_bh) or r.p_bh>=.05:continue
   label=r.comparison.replace('(','').replace(')','')
   groups=label.split(' - ')
   if len(groups)!=2 or any(x not in GROUPS for x in groups):continue
   a,b=[GROUPS.index(x) for x in groups];ax.plot([a,a,b,b],[y,y+step*.2,y+step*.2,y],color='black',lw=1)
   ax.text((a+b)/2,y+step*.25,'***'if r.p_bh<.001 else'**'if r.p_bh<.01 else'*',ha='center');y+=step
 save(fig,stem,legend,data)

def tissue_objects(d,background,field,label,stem,note):
 fig,axes=plt.subplots(5,4,figsize=(17,20),layout='constrained');lo=d[field].min();hi=d[field].max()
 donors=background[['case_id','donor_status']].drop_duplicates().sort_values(['donor_status','case_id'])
 for ax,r in zip(axes.flat,donors.itertuples()):
  b=background[background.case_id.eq(r.case_id)];q=d[d.case_id.eq(r.case_id)]
  ax.scatter(b.centroid_x_um,b.centroid_y_um,c='#dddddd',s=4)
  sc=ax.scatter(q.centroid_x_um,q.centroid_y_um,c=q[field],s=18,cmap='viridis',vmin=lo,vmax=hi)
  ax.set(title=f'{r.donor_status} · {r.case_id}',xlabel='x (µm)',ylabel='y (µm)');ax.invert_yaxis();ax.set_aspect('equal',adjustable='box')
 fig.colorbar(sc,ax=list(axes.flat),label=label,shrink=.5)
 p=DATA/f'{stem}_plot_data.csv';d[['islet_id','case_id','donor_status','centroid_x_um','centroid_y_um',field]].to_csv(p,index=False)
 save(fig,stem,note+' Each panel is one donor; colored points are eligible endocrine-object centroids and gray points are the available corrected-object background. Coordinates are physical micrometres, with image y increasing downwards. Color uses one common measurement scale across all donors. Donor identifiers appear in panel labels. No tests are encoded.',[p])

def reused_results():
 for name in ['model_anova.csv','model_coefficients.csv','disease_slopes.csv','disease_slope_comparisons.csv','donor_regression_eligibility.csv','donor_disease_global_tests.csv','donor_disease_pairwise_tests.csv','donor_direction_tests.csv']:
  d=pd.read_csv(OLD/'h6'/name);d=d[d.model.str.startswith('beta_density')]
  d.to_csv(MODELS/f'h6_{name}',index=False)
 for name in ['donor_field_coefficients.csv','disease_field_overall_tests.csv','disease_field_pairwise_tests.csv','donor_object_class_densities.csv','object_class_density_overall_tests.csv','object_class_density_pairwise_tests.csv','model_specifications.csv']:
  shutil.copy2(OLD/'h8'/name,MODELS/f'h8_{name}')
 a=pd.read_csv(MODELS/'h6_model_anova.csv');r=a[(a.model=='beta_density_disease_interaction_correlated')&(a.term=='contact_z:donor_status')].iloc[0]
 primary=[{'hypothesis':'H6','comparison':'ND vs Aab+ vs T1D: excess-adjacency/Beta-density slope','test':'Satterthwaite F test of disease-by-adjacency interaction','p_value':r['Pr(>F)'],'n_donors':19,'n_objects':721,'status':'estimated; existing nonsingular donor random-slope model reused'}]
 g=pd.read_csv(MODELS/'h7_global.csv');r=g[(g.model=='primary_20um')&(g.role!='supporting_sensitivity')].iloc[0]
 primary.append(dict(hypothesis='H7',comparison=r.comparison,test=r.test,p_value=r.p_value,n_donors=r.n_donors,n_objects=r.n_objects,status=r.status))
 d=pd.read_csv(MODELS/'h8_disease_field_overall_tests.csv');r=d[(d.model=='primary_200um')&(d.field=='immune')].iloc[0]
 primary.append(dict(hypothesis='H8',comparison='ND vs Aab+ vs T1D: immune-field/small-Beta-Cluster-density coefficient',test='Kruskal-Wallis on adjusted independent donor field coefficients',p_value=r.p,n_donors=r.donors,n_objects=8362,status='estimated; existing full-tile model reused; n_objects denotes tissue tiles'))
 pd.DataFrame(primary).to_csv(OUT/'primary_tests.csv',index=False)
 pd.DataFrame(primary).to_csv(OUT/'primary_comparisons.csv',index=False)

def six_plots(d,base):
 states(d,[('combined_excess_percentage_points','Excess β–α/δ adjacency (percentage points)'),('beta_density_mm2','β-cell density (cells/mm²)')],'h6_adjacency',
  'The 721 graph-eligible objects from 19 donors retain the original fixed 20 µm centroid graph, endocrine composition and 1999 label permutations per object (seed 20261026). The Beta-density analysis uses no INS measurement. Representation eligibility only limits these descriptive state panels; the complete analysis cohort remains in the model CSV.')
 s=pd.read_csv(MODELS/'h6_disease_slopes.csv').rename(columns={'contact_z.trend':'estimate','lower.CL':'lower','upper.CL':'upper','p.value':'p_value'})
 s['p_bh']=bh(s.p_value);s.to_csv(DATA/'h6_conditional_slopes.csv',index=False)
 p=pd.read_csv(MODELS/'h6_disease_slope_comparisons.csv').rename(columns={'contrast':'comparison','p.value':'p_value'})
 group_slopes(s,'h6_beta_density_slopes','Adjacency and β-cell density',
  'Points and 95% confidence intervals are the retained disease-specific standardized regression slopes for raw Beta-cell density versus combined excess Beta–Alpha/Delta adjacency. The model controls Beta, Delta and other-endocrine fractions, log core area, log1p nearby immune count and log1p core immune count, with correlated donor random intercept and contact slope. Predictors were standardized within donor and outcome across 721 objects/19 donors. This existing fit is nonsingular. Direct disease differences use three BH-adjusted contrasts; within-state slope tests form a separate three-test BH family. Asterisks encode only significant adjusted disease contrasts. Numerical values: '+s[['donor_status','estimate','lower','upper','p_value','p_bh']].to_json(orient='records')+'; pairwise: '+p[['comparison','p_value','p_bh']].to_json(orient='records')+'. The model outcome is cell count per endocrine-core area, not hormone intensity.',[DATA/'h6_conditional_slopes.csv',MODELS/'h6_disease_slope_comparisons.csv'],p)
 donor=d.groupby(['case_id','donor_status']).agg(observed=('observed_adjacency_percent','mean'),permuted=('permuted_adjacency_percent','mean'),objects=('islet_id','size')).reset_index()
 donor.to_csv(DATA/'h6_observed_permuted_donor.csv',index=False)
 fig,axes=plt.subplots(1,3,figsize=(10,4.5),layout='constrained')
 for ax,g in zip(axes,GROUPS):
  z=donor[donor.donor_status.eq(g)]
  for r in z.itertuples():ax.plot([0,1],[r.permuted,r.observed],'-o',color=COLORS[g],alpha=.7,ms=4)
  ax.set(title=g,xticks=[0,1],xticklabels=['Permuted','Observed'],ylabel='Combined adjacency (%)')
 save(fig,'h6_observed_permuted','Each connected pair is one donor: the arithmetic mean of object observed combined Beta–Alpha/Delta edge fractions and of the corresponding 1999-permutation means, shown as percentages. Graph edges, coordinates and cell-type counts are preserved in each permutation. Donor counts are ND 7, Aab+ 7, T1D 5. These are descriptive donor comparisons; original paired-test results remain in the accompanying table. No INS measurement enters this plot.',[DATA/'h6_observed_permuted_donor.csv',MODELS/'h6_observed_vs_permuted_donor_tests.csv'])
 shutil.copy2(OLD/'h6/observed_vs_permuted_donor_tests.csv',MODELS/'h6_observed_vs_permuted_donor_tests.csv')
 # Reuse the pre-existing matching strata; replace only the superseded response with the retained quantity endpoint.
 keys=['islet_id','beta_fraction_stratum','delta_fraction_stratum','log_core_area_stratum','positive_excess']
 assignment=pd.read_csv(OLD/'h6/matched_strata_object_assignments.csv',usecols=keys)
 m=d[['islet_id','case_id','donor_status','beta_density_mm2']].merge(assignment,on='islet_id',validate='one_to_one')
 strata=['case_id','donor_status','beta_fraction_stratum','delta_fraction_stratum','log_core_area_stratum']
 paired=m.groupby(strata+['positive_excess']).beta_density_mm2.mean().unstack().dropna().reset_index().rename(columns={False:'nonpositive_excess',True:'positive_excess'})
 paired.to_csv(DATA/'h6_matched_beta_density_strata.csv',index=False)
 fig,axes=plt.subplots(1,3,figsize=(11,4.5),layout='constrained')
 for ax,g in zip(axes,GROUPS):
  for r in paired[paired.donor_status.eq(g)].itertuples():ax.plot([0,1],[r.nonpositive_excess,r.positive_excess],'-o',color=COLORS[g],alpha=.35,ms=3)
  ax.set(title=g,xticks=[0,1],xticklabels=['Nonpositive excess','Positive excess'],ylabel='Mean β-cell density (cells/mm²)')
 save(fig,'h6_matched_beta_density','Connected values are the mean Beta-cell densities in positive versus nonpositive combined-excess groups within existing donor-specific composition/size strata. Stratum boundaries are the previously fixed pooled tertiles of Beta fraction, Delta fraction and log core area. Only strata containing both signs are connected; exact strata, donor IDs and values are exported. This is a descriptive comparison using the original matching definitions; it does not add a test or claim exact matching/immune adjustment. Core Beta count divided by core area is used; no INS abundance is used.',[DATA/'h6_matched_beta_density_strata.csv'])
 tissue_objects(d,base,'combined_excess_percentage_points','Excess adjacency (percentage points)','h6_tissue_adjacency','The complete 721-object analysis cohort is mapped; all 20 donors are shown, including the donor without an eligible graph.')
 # Deterministic physical graph examples: largest eligible endocrine-cell graph per donor.
 examples=d.sort_values(['case_id','n_endocrine_nodes','islet_id']).groupby('case_id').tail(1)
 fig,axes=plt.subplots(5,4,figsize=(15,18),layout='constrained');example_rows=[];edge_rows=[]
 for ax,r in zip(axes.flat,examples.itertuples()):
  g=np.load(OLD/'h6/graphs_by_donor'/str(r.case_id)/f'{r.islet_id}.npz');xy=g['coordinates'];edges=g['edges'];labels=g['labels']
  ax.add_collection(LineCollection(xy[edges],color='#cccccc',lw=.35,zorder=1))
  ax.scatter(xy[:,0],xy[:,1],c=np.asarray(['#3366AA','#CC6633','#228833','#888888'])[labels],s=8,zorder=2)
  ax.autoscale();ax.invert_yaxis();ax.set_aspect('equal');ax.set(title=r.islet_id,xlabel='x (µm)',ylabel='y (µm)')
  for a,b in edges:edge_rows.append({'islet_id':r.islet_id,'case_id':r.case_id,'cell_id_1':str(g['object_ids'][a]),'cell_id_2':str(g['object_ids'][b]),'x1_um':xy[a,0],'y1_um':xy[a,1],'x2_um':xy[b,0],'y2_um':xy[b,1]})
  for oid,x,y,label in zip(g['object_ids'],xy[:,0],xy[:,1],labels):example_rows.append({'islet_id':r.islet_id,'case_id':r.case_id,'cell_id':str(oid),'x_um':x,'y_um':y,'endocrine_label_code':int(label)})
 for ax in axes.flat[len(examples):]:ax.axis('off')
 pd.DataFrame(example_rows).to_csv(DATA/'h6_contact_graph_cells.csv',index=False)
 pd.DataFrame(edge_rows).to_csv(DATA/'h6_contact_graph_edges.csv',index=False)
 save(fig,'h6_contact_graph_examples','The largest eligible endocrine-cell graph in each of 19 donors is shown, selected by node count with object-ID tie breaking. Coordinates and edges are reused without reconstruction. Blue: Beta; orange: Alpha; green: Delta; gray: other endocrine. Gray lines are the existing 20 µm-centroid-neighborhood graph edges. This is the graph whose labels were permuted; a centroid-neighborhood edge is not a physical membrane contact. Object and cell IDs and coordinates are exported; full graph arrays remain at the cited source path.',[DATA/'h6_contact_graph_cells.csv',DATA/'h6_contact_graph_edges.csv'])

def seven_plots(d,base):
 q=d[d.h7_primary_eligible].copy()
 states(q,[('nerve_near_cd31_20um','Nerve area within 20 µm of CD31 (fraction)'),('beta_density_mm2','β-cell density (cells/mm²)')],'h7_neurovascular',
  'The model-eligible proximity cohort includes zero-Beta objects and uses only existing INS-independent B3TUBB masks, CD31 union masks and corrected typed-cell counts. All 4529 core+20 µm-complete paired objects were assessed; undefined proximity where no nerve area exists is not converted to zero.')
 fields=[('structure_nerve_peri_density_cells_mm2','Recovered nerve-cell density (cells/mm²)'),('structure_nerve_peri_mean_passthrough_fraction_among_measured','Mean measured nerve passthrough fraction'),('nerve_area_fraction','B3TUBB-positive peri-area fraction')]
 states(d[d.core_complete&d.neighborhood_complete_20],fields,'h7_neural_context',
  'Recovered nerve-cell calls, passthrough fractions and B3TUBB-positive area are separate measurements. Cell calls use centroid, majority whole-cell area, marker level and highest-marker criteria with union overlap; app phenotypes remain unchanged. Passthrough means use measured centroid-inside cells without the structure own call. Their missing values remain undefined. These descriptive context fields are not added to the frozen inferential model.')
 s=pd.read_csv(MODELS/'h7_slopes.csv');p=pd.read_csv(MODELS/'h7_pairs_out.csv');g=pd.read_csv(MODELS/'h7_global.csv')
 primary=s[(s.model=='primary_20um')&(s.role!='supporting_sensitivity')];pairs=p[(p.model=='primary_20um')&(p.role!='supporting_sensitivity')]
 primary.to_csv(DATA/'h7_primary_conditional_slopes.csv',index=False)
 group_slopes(primary,'h7_beta_density_slopes','Neurovascular proximity and β-cell density',
  'Points and 95% confidence intervals show disease-specific standardized coefficients for log1p Beta-cell density versus 20 µm nerve–CD31 proximity. Core count/core area supplies cells/mm²; zero-Beta objects remain zero before log1p. The frozen model adjusts separate B3TUBB and CD31 peri-area fractions, log core area, log perimeter and log1p nearby immune density. Outcome is globally standardized and predictors within donor. A correlated donor random-intercept/slope model was the first valid prespecified fit. All 2364 objects from 20 donors enter, including 1325 zero-Beta objects. Three direct disease slope contrasts share BH correction, and three within-state tests form a separate BH family. Asterisks mark significant adjusted disease contrasts. Global test: '+g[g.model.eq('primary_20um')].to_json(orient='records')+'; slopes: '+primary.to_json(orient='records')+'; contrasts: '+pairs.to_json(orient='records')+'. Coefficients are regressions, not correlations.',[DATA/'h7_primary_conditional_slopes.csv',MODELS/'h7_pairs_out.csv',MODELS/'h7_global.csv'],pairs)
 fig,axes=plt.subplots(1,4,figsize=(15,4.5),layout='constrained')
 fig.suptitle('Neurovascular proximity and β-cell density')
 names=['distance_10um','distance_5um','overlap','GAP43_20um'];titles=['10 µm proximity','5 µm proximity','Nerve–vessel overlap','GAP43-supported 20 µm proximity']
 for ax,name,title in zip(axes,names,titles):
  z=s[(s.model==name)&(s.role!='supporting_sensitivity')]
  for i,group in enumerate(GROUPS):
   r=z[z.donor_status.eq(group)].iloc[0];ax.errorbar(i,r.estimate,yerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',color=COLORS[group],capsize=4)
  ax.axhline(0,color='#888',ls='--',lw=1);ax.set(title=title,xticks=range(3),xticklabels=GROUPS,ylabel='Adjusted standardized β')
 save(fig,'h7_proximity_sensitivities','The four predeclared proximity/mask sensitivities use the same log1p Beta-density outcome, five fixed covariates and donor random-intercept/slope family as the primary analysis. Points are disease-specific standardized coefficients with 95% CI. Secondary global disease-interaction tests share one four-test BH family; each specification retains separate three-test pairwise and within-state families. Complete-case cohorts are 2364 objects for 10 µm, 5 µm and overlap, and 2250 objects for GAP43-supported proximity, all 20 donors. No INS measurements are used. Exact statistics and cohort sizes are in the source tables.',[MODELS/'h7_slopes.csv',MODELS/'h7_global.csv',MODELS/'h7_pairs_out.csv'])
 tissue_objects(q,base,'nerve_near_cd31_20um','Nerve area within 20 µm of CD31 (fraction)','h7_tissue_proximity','Defined 20 µm proximity is shown for 2364 objects from 20 donors; absent measured nerve area has undefined proximity and remains gray in the background.')
 paths=[OLD/'h7'/f'example_{case}.npz' for case in [6516,6450,6533]]
 fig,axes=plt.subplots(1,3,figsize=(12,4.5),layout='constrained');examples=[]
 for ax,path in zip(axes,paths):
  f=np.load(path);nerve=f['nerve'].astype(bool);vessel=f['cd31'].astype(bool);peri=f['peri'].astype(bool)
  rgb=np.zeros((*peri.shape,3));rgb[peri]=[.12,.12,.12];rgb[vessel&peri]=[.95,.2,.7];rgb[nerve&peri]=[.1,.8,.4];rgb[nerve&vessel&peri]=[1,1,1]
  px=float(f['pixel_um']);ax.imshow(rgb,extent=[0,peri.shape[1]*px,peri.shape[0]*px,0]);oid=str(f['islet_id']);ax.set(title=oid,xlabel='x within crop (µm)',ylabel='y within crop (µm)')
  examples.append({'islet_id':oid,'source_npz':str(path),'pixel_um':px})
 pd.DataFrame(examples).to_csv(DATA/'h7_mask_examples.csv',index=False)
 save(fig,'h7_neurovascular_masks','Three previously saved physical measurement examples are shown with original object identifiers and calibrated micrometre axes. Dark gray: assessed peri region; green: INS-independent B3TUBB-positive nerve mask; magenta: existing CD31 structure union; white: overlap. The same Otsu/minimum 2 µm² component method was used for all donors. These are binary measurement masks, not quantitative INS displays. The saved examples were selected in the earlier extraction as median-area eligible objects in the listed donors; no selection by the new Beta-density association occurred.',[DATA/'h7_mask_examples.csv'])

def eight_plots(census,tiles):
 q=census[census.density_eligible & census.matched_positive_area_tile].copy()
 features=[('small_beta_cluster_observed_density_mm2','Local small β-bearing Cluster density (objects/mm²)'),('immune_cells_density_mm2','Local immune-cell density (cells/mm²)'),('epithelial_signal_density_mm2','Epithelial marker signal\n(summed scaled signal/mm²)'),('stromal_signal_density_mm2','Stromal marker signal\n(summed scaled signal/mm²)')]
 states(q,features,'h8_field_states','Object centroids are linked to the already defined 200 µm tile. State overlays show the local tile measurement at saved embedded objects; the complete 10113-object census and 8362 positive-area tiles remain in the companion tables, including small nonembedded objects and object-free tiles. Epithelial context is the sum across eligible extraendocrine cells of each cell’s mean cohort-scaled CK19/panCK signal, divided by valid tissue area in mm². Stromal context analogously averages SMA/VIM/ColIV within each cell, sums across eligible extraendocrine cells and divides by valid tissue area. These are summed scaled signal/mm²; neither uses INS.')
 coeff=pd.read_csv(MODELS/'h8_donor_field_coefficients.csv');tests=pd.read_csv(MODELS/'h8_disease_field_overall_tests.csv');pairs=pd.read_csv(MODELS/'h8_disease_field_pairwise_tests.csv')
 primary=coeff[coeff.model.eq('primary_200um')].copy();primary.to_csv(DATA/'h8_primary_field_coefficients.csv',index=False)
 fig,axes=plt.subplots(1,3,figsize=(12,4.7),layout='constrained')
 for ax,field in zip(axes,['immune','epithelial','stromal']):
  for i,group in enumerate(GROUPS):
   z=primary[(primary.field==field)&(primary.donor_status==group)].sort_values('case_id');v=z.standardized_beta.to_numpy();xs=i+np.linspace(-.10,.10,len(v))
   ax.scatter(xs,v,color=COLORS[group],s=25);ax.boxplot(v,positions=[i],widths=.4,showfliers=False,medianprops={'color':'black'},boxprops={'color':COLORS[group]},whiskerprops={'color':COLORS[group]},capprops={'color':COLORS[group]})
  ax.axhline(0,color='#888',ls='--',lw=1);ax.set(title=field.capitalize()+' field',xticks=range(3),xticklabels=GROUPS,ylabel='Adjusted standardized donor coefficient (β)')
  vals=primary.loc[primary.field.eq(field),'standardized_beta'];step=max((vals.max()-vals.min())*.14,.006);y=vals.max()+step
  for r in pairs[(pairs.model=='primary_200um')&(pairs.field==field)&(pairs.p_BH_all_nine<.05)].itertuples():
   a,b=GROUPS.index(r.group1),GROUPS.index(r.group2);ax.plot([a,a,b,b],[y,y+step*.2,y+step*.2,y],color='black',lw=1)
   ax.text((a+b)/2,y+step*.25,'***'if r.p_BH_all_nine<.001 else'**'if r.p_BH_all_nine<.01 else'*',ha='center');y+=step
 fig.suptitle('Local tissue fields and small β-bearing Cluster density')
 save(fig,'h8_field_disease_comparison','Each point is one donor-specific adjusted field coefficient from the retained full 200 µm tile analysis. Boxes show median and quartiles, with 1.5 IQR whiskers. The response is small Beta-bearing Cluster density per valid tissue area; the three fields are log1p immune density, epithelial context and stromal context, mutually adjusted with valid tile area and marker-missing fraction. Predictors are standardized within donor and outcome globally. All 20 donors and 8362 positive-area tiles, including object-free tiles, remain. Primary comparison: immune-field coefficient across disease states by Kruskal–Wallis. Existing three-field omnibus and nine field-by-disease pairwise BH families are retained. Brackets use the nine-comparison BH q values (*q<0.05; **q<0.01; ***q<0.001). Exact global tests: '+tests[tests.model.eq('primary_200um')].to_json(orient='records')+'; pairwise: '+pairs[pairs.model.eq('primary_200um')].to_json(orient='records')+'. No model was refitted for the state linkage.',[DATA/'h8_primary_field_coefficients.csv',MODELS/'h8_disease_field_overall_tests.csv',MODELS/'h8_disease_field_pairwise_tests.csv'])
 sens=coeff[coeff.field.eq('immune')].copy();sens.to_csv(DATA/'h8_immune_field_sensitivities.csv',index=False)
 names=['primary_200um','grid100','grid400','identity_lower','identity_upper','wholecell_size','area_size']
 labels=['200 µm','100 µm','400 µm','Identity lower','Identity upper','Whole-cell size','Area size']
 fig,ax=plt.subplots(figsize=(13,5),layout='constrained')
 for j,g in enumerate(GROUPS):
  for i,name in enumerate(names):
   v=sens[(sens.donor_status==g)&(sens.model==name)].standardized_beta.to_numpy();x=i+(j-1)*.23
   ax.scatter(x+np.linspace(-.04,.04,len(v)),v,s=10,color=COLORS[g],alpha=.65,label=g if i==0 else None)
   ax.plot([x-.08,x+.08],[np.median(v)]*2,color=COLORS[g],lw=2)
 ax.axhline(0,color='#888',ls='--',lw=1);ax.set(xticks=range(7),xticklabels=labels,ylabel='Adjusted immune-field donor coefficient (β)',title='Small β-bearing Cluster density and immune field');ax.legend(frameon=False)
 save(fig,'h8_field_sensitivities','Each point is an existing donor-specific adjusted immune-field coefficient; horizontal marks indicate group medians. The previously computed 100/200/400 µm grids, identity lower/upper bounds, whole-cell count size and geometric-area size definitions are retained unchanged. All 20 donors remain in each specification. Exact tests and existing BH families are in the accompanying tables; this figure does not introduce another test or fit.',[DATA/'h8_immune_field_sensitivities.csv',MODELS/'h8_disease_field_overall_tests.csv',MODELS/'h8_disease_field_pairwise_tests.csv'])
 donors=tiles[['case_id','donor_status']].drop_duplicates().sort_values(['donor_status','case_id'])
 fig,axes=plt.subplots(5,4,figsize=(17,20),layout='constrained');vals=np.log1p(tiles.small_beta_cluster_observed_density_mm2);norm=Normalize(0,vals.max())
 for ax,r in zip(axes.flat,donors.itertuples()):
  z=tiles[tiles.case_id.eq(r.case_id)];patches=[Rectangle((row.x_min_um,row.y_min_um),200,200) for row in z.itertuples()]
  coll=PatchCollection(patches,cmap='viridis',norm=norm,edgecolor='none');coll.set_array(np.log1p(z.small_beta_cluster_observed_density_mm2));ax.add_collection(coll);ax.autoscale();ax.invert_yaxis();ax.set_aspect('equal');ax.set(title=f'{r.donor_status} · {r.case_id}',xlabel='x (µm)',ylabel='y (µm)')
 fig.colorbar(coll,ax=list(axes.flat),label='log1p small β-bearing Cluster density (objects/mm²)',shrink=.5)
 save(fig,'h8_full_tile_census','Every one of the 8362 positive-valid-area 200 µm tiles is shown, including zero-object and zero-endpoint tiles. Color is log1p small Beta-bearing Cluster count per valid tissue area; tile footprint is drawn as its 200 µm grid square, while the numerical denominator is the measured tissue/edge-supported portion within it. Missing tiles remain absent, not zero. Grid origin and eligibility are inherited from the existing full census. Shared color scale and physical micrometre coordinates apply across 20 donor panels. No fitting or tests occur in this map.',[DATA/'h8_tiles_with_states.csv'])
 fig,axes=plt.subplots(5,4,figsize=(17,20),layout='constrained');cmap=plt.get_cmap('tab10')
 for ax,r in zip(axes.flat,donors.itertuples()):
  z=census[census.case_id.eq(r.case_id)];unknown=~z.in_composition_embedding
  ax.scatter(z.loc[unknown,'centroid_x_um'],z.loc[unknown,'centroid_y_um'],s=5,c='#cccccc',label='No saved state')
  known=z[z.in_composition_embedding];ax.scatter(known.centroid_x_um,known.centroid_y_um,c=known['leiden_0.5'],cmap=cmap,vmin=-.5,vmax=9.5,s=12)
  ax.invert_yaxis();ax.set_aspect('equal');ax.set(title=f'{r.donor_status} · {r.case_id}',xlabel='x (µm)',ylabel='y (µm)')
 handles=[plt.Line2D([],[],marker='o',ls='',color=cmap(i),label=str(i)) for i in range(10)]+[plt.Line2D([],[],marker='o',ls='',color='#cccccc',label='No saved state')]
 fig.legend(handles=handles,loc='outside lower center',ncol=11,title='Saved composition cluster')
 save(fig,'h8_census_composition_states','All 10113 geometric endocrine objects, including Singles and small nonembedded Clusters, are shown at their existing centroids. Colors identify the saved composition Leiden 0.5 state where available; gray means no saved composition state, not absence of endocrine cells or exclusion from the census. The state labels are reused without embedding or clustering refits. Coordinates use the same current geometry as the linked tile census. This display does not replace the full-tile primary analysis.',[DATA/'h8_census_with_states.csv'])

def audit_report():
 rows=[
 ('H6','Permutation excess combined adjacency','completed; reused','1999 permutations preserving cell-type counts and saved graph geometry'),
 ('H6','Within-Beta hormone abundance response','superseded','Quantitative INS endpoint not used'),
 ('H6','Beta-density association and direct disease comparison','completed; reused','Existing 721-object 19-donor valid correlated random-slope model'),
 ('H6','Matched composition/size comparison','completed descriptively for Beta density','Existing fixed strata reused; no additional test'),
 ('H6','Separate Alpha/Delta, boundary contacts, contact-by-immune and neurovascular-adjusted Beta-density follow-ups','not run for amended quantity presentation','Existing measurement columns are retained; historical INS-response fits are not presented as Beta-density tests'),
 ('H7','Common INS-independent B3TUBB mask and observed 20 µm proximity','completed and coverage extended','Existing Otsu/component method and cached thresholds;350 additional ROIs'),
 ('H7','Original excess-over-spatial-null alignment','superseded by approved observed-proximity method','No claim of excess alignment over a spatial null'),
 ('H7','Log1p Beta-density response, zero-Beta retained, direct disease interaction','completed','Frozen protocol before fitting; all required abundance/geometry/immune covariates retained'),
 ('H7','5/10µm, overlap and GAP43-supported sensitivities','completed','Same model family and declared multiplicity'),
 ('H7','Core/peri/distal length and alignment; uncertain-image omission','not run in this quantity update','Available peri masks and old source examples retained; no new exclusion cutoffs invented'),
 ('H8','Full 200 µm census and donor field disease comparisons','completed; reused','8362 positive-area tiles; object-free tiles and all 10113 census objects retained'),
 ('H8','100/400µm, identity bounds, whole-cell/area size sensitivities','completed; reused','Previously fixed definitions and tests retained'),
 ('H8','Composition/CHGA-hormone-negative and age-matched/Ki67 follow-ups','partially completed descriptively or not run','Saved state linkage added; age-matched and Ki67 tests not run here'),
 ('H6–H8','Original if-supported/meaningful-effect-excluded branches','judgment not assigned','No numerical meaningful-effect bound was prespecified; no biological interpretation requested')]
 pd.DataFrame(rows,columns=['hypothesis','original_condition_or_followup','execution_status','detail']).to_csv(OUT/'original_plan_execution_audit.csv',index=False)
 text='H6–H8 measurements and linked app displays\n\nThe shared saved composition UMAP/Leiden and app pseudotime are reused without fitting or reranking. Complete measurement/model cohorts and full H8 tissue census remain outside the embedding views. All new exports use accepted cell counts and non-INS structural/context measurements.\n\nH6 retains its existing Beta-density model, permutation measurements and 721-object cohort. H7 uses the protocol recorded before fitting: log1p Beta density, zero-Beta objects retained, observed 20 µm proximity and the existing five abundance/geometry/immune controls. The 350 additional ROIs use cached thresholds and the existing measurement helper. H8 retains its existing 8362-tile primary tests and sensitivities; state linkage does not refit the model.\n\nprimary_tests.csv supplies the three prespecified disease comparisons. Model tables include within-state conditional coefficients and direct contrasts, with separate multiplicity families. Figure legends contain methods, populations and exact statistics. original_plan_execution_audit.csv records completed, superseded and unrun follow-ups without assigning biological branches.\n'
 (OUT/'README.txt').write_text(text)

def main():
 reused_results();six=S.read(DATA/'h6_objects.parquet');seven=S.read(DATA/'h7_objects.parquet');census=S.read(DATA/'h8_census_with_states.parquet');tiles=S.read(DATA/'h8_tiles_with_states.csv');base=S.read(ROOT/'output/app_hypotheses_20261026/shared/object_measurements.parquet')
 six_plots(six,base);print('H6 figures complete',flush=True)
 seven_plots(seven,base);print('H7 figures complete',flush=True)
 eight_plots(census,tiles);print('H8 figures complete',flush=True)
 audit_report()
 (OUT/'figure_manifest.json').write_text(json.dumps({'figures':records,'native_formats':['PNG','PDF','SVG'],'no_quantitative_INS':True,'interpretation_added':False},indent=2)+'\n')
 print('Figures',len(records),flush=True)
if __name__=='__main__':main()
