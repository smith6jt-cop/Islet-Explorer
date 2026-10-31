#!/usr/bin/env python3
"""Amended H5: saved composition-state occupancy and categorical cell quantities."""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations
from app_hypotheses_shared_20261026 import ROOT,OUT,GROUPS,COLORS,read,savefig,state_views,plt,digest

DEST=OUT/'h5'
def prepare():
    (DEST/'data').mkdir(parents=True,exist_ok=True);(DEST/'figures').mkdir(exist_ok=True)
    # Written before any statistical calculation or fitting for this amendment.
    protocol={'hypothesis':'H5 amended','primary_endpoint':'donor percentage of core-complete paired endocrine objects containing at least one accepted Beta cell','primary_test':'Kruskal-Wallis across ND/Aab+/T1D; reserves one slot in current ten-hypothesis disease-comparison family','primary_pairwise':'all three two-sided Wilcoxon rank-sum comparisons; BH within these three','secondary_endpoints':['saved Leiden 0.5 occupancy percentage for each of 10 clusters','equal-object donor means of Beta/Alpha/Delta/CHGA+ hormone-negative endocrine fractions','equal-object donor means of Beta and CHGA+ hormone-negative cell densities per mm2 of core'],'secondary_multiplicity':'BH across 16 omnibus comparisons; separate BH across their 48 pairwise comparisons','cluster_cohort':'saved composition embedding with complete core and 20um neighborhood; occupancy denominator all such objects within donor; absent clusters zero','quantity_cohort':'all paired core-complete objects, including zero-Beta objects; endocrine fractions require a positive endocrine denominator','clinical_replicate':'donor','pairwise_method':'R wilcox.test two-sided; exact if no ties, otherwise asymptotic tie-corrected with continuity correction','plot_families':'donor points/median, common UMAP, 20-bin app pseudotime profiles, cluster heatmaps, physical object maps','unknown_unpaired_census':'retain separately with measured/unknown identities; do not manufacture absent Beta cells','INS_intensity_endpoint':False,'original_H5_replaced_not_equivalent':True,'branch_judgment':'not assigned'}
    (DEST/'analysis_specification_before_fit.json').write_text(json.dumps(protocol,indent=2)+'\n')
    d=read(OUT/'shared/object_measurements.parquet');q=d[d.core_complete].copy();e=d[d.eligible_plotted].copy()
    rows=[];occ=[]
    for (case,g),z in q.groupby(['case_id','donor_status']):
        vals={'beta_presence_percent':z.beta_presence.mean()*100}
        for cell in ['beta','alpha','delta','chga_hormone_negative']:vals[cell+'_endocrine_percent']=z[cell+'_endocrine_fraction'].mean()*100
        for cell in ['beta','chga_hormone_negative']:vals[cell+'_density_mm2']=z[cell+'_density_mm2'].mean()
        for key,value in vals.items():
            n=int(z.n_endocrine_core.gt(0).sum()) if 'endocrine_percent' in key else len(z)
            rows.append({'case_id':case,'donor_status':g,'endpoint':key,'value':value,'objects':n,'family':'primary' if key=='beta_presence_percent' else 'secondary'})
        zz=e[e.case_id.eq(case)]
        for cl in range(10):
            num=int(zz['leiden_0.5'].eq(cl).sum());den=len(zz)
            occ.append({'case_id':case,'donor_status':g,'cluster':cl,'objects_in_cluster':num,'eligible_objects':den,'occupancy_percent':100*num/den if den else np.nan})
            rows.append({'case_id':case,'donor_status':g,'endpoint':f'cluster_{cl}_occupancy_percent','value':100*num/den if den else np.nan,'objects':den,'family':'secondary'})
    pd.DataFrame(rows).to_csv(DEST/'data/donor_measurements.csv',index=False)
    pd.DataFrame(occ).to_csv(DEST/'data/donor_cluster_occupancy.csv',index=False)
    q.to_parquet(DEST/'data/core_complete_objects.parquet',index=False)
    cols=['islet_id','case_id','donor_status','core_area_um2','n_beta_core','n_endocrine_core','n_core_CHGAplus_hormone-negative','beta_presence','beta_density_mm2','chga_hormone_negative_density_mm2','beta_endocrine_fraction','alpha_endocrine_fraction','delta_endocrine_fraction','chga_hormone_negative_endocrine_fraction','in_composition_embedding','eligible_plotted']
    q[cols].to_csv(DEST/'data/quantity_objects.csv',index=False)
    census=read(ROOT/'output/h5_h10_original_plan_20261026/h8/full_object_census.parquet')
    census.to_csv(DEST/'data/full_geometric_census.csv',index=False)
    # Full census companions are exact counts and valid-tissue densities; Singles
    # have count/density outputs only, without states or object response models.
    dc=read(ROOT/'output/h5_h10_original_plan_20261026/h8/donor_census.csv')
    dc.to_csv(DEST/'data/full_census_donor_quantities.csv',index=False)
    print('H5 protocol and donor inputs saved; no statistical fit executed by prepare.')

LABELS={'beta_presence_percent':'Objects containing beta cells (%)','beta_endocrine_percent':'Beta / endocrine cells (%)','alpha_endocrine_percent':'Alpha / endocrine cells (%)','delta_endocrine_percent':'Delta / endocrine cells (%)','chga_hormone_negative_endocrine_percent':'CHGA+ hormone-negative / endocrine cells (%)','beta_density_mm2':'Beta cells / core area (cells/mm²)','chga_hormone_negative_density_mm2':'CHGA+ hormone-negative cells\nper core area (cells/mm²)'}
TITLES={'beta_presence_percent':'Beta-cell presence','beta_endocrine_percent':'Beta-cell composition','alpha_endocrine_percent':'Alpha-cell composition','delta_endocrine_percent':'Delta-cell composition','chga_hormone_negative_endocrine_percent':'CHGA+ hormone-negative composition','beta_density_mm2':'Beta-cell density','chga_hormone_negative_density_mm2':'CHGA+ hormone-negative density'}

def star(p):return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else ''

def figures():
    d=read(OUT/'shared/object_measurements.parquet');v=read(DEST/'data/donor_measurements.csv');tests=read(DEST/'omnibus_tests.csv');pairs=read(DEST/'pairwise_tests.csv')
    def dots(keys,stem,ncol=3):
        fig,axes=plt.subplots(int(np.ceil(len(keys)/ncol)),ncol,figsize=(4.5*ncol,3.6*int(np.ceil(len(keys)/ncol))),squeeze=False,layout='constrained')
        plot=[]
        for ax,key in zip(axes.flat,keys):
            r=v[v.endpoint.eq(key)];label=LABELS.get(key,'Cluster '+key.split('_')[1]+' occupancy (%)')
            for k,g in enumerate(GROUPS):
                z=r[r.donor_status.eq(g)].sort_values('case_id');j=np.linspace(-.11,.11,len(z));ax.scatter(k+j,z.value,c=COLORS[g],s=29,zorder=3)
                ax.plot([k-.2,k+.2],[z.value.median()]*2,c='black',lw=1.4)
                for xx,row in zip(k+j,z.itertuples()):plot.append({'endpoint':key,'case_id':row.case_id,'donor_status':g,'x':xx,'value':row.value})
            t=tests[tests.endpoint.eq(key)].iloc[0];p=t.p_value if key=='beta_presence_percent' else t.p_bh_family
            if np.isfinite(p) and p<.05:ax.text(.5,.96,star(p),transform=ax.transAxes,ha='center',va='top',fontsize=16)
            ax.set(xticks=range(3),xticklabels=GROUPS,ylabel=label,title=TITLES.get(key,'Cluster '+key.split('_')[1]));ax.grid(axis='y',alpha=.15)
            if r.value.min()>=0:ax.set_ylim(bottom=0,top=max(float(r.value.max())*1.18,.01))
        for ax in list(axes.flat)[len(keys):]:ax.set_visible(False)
        savefig(fig,DEST/'figures',stem);pd.DataFrame(plot).to_csv(DEST/'data'/f'{stem}_points.csv',index=False)
        legend='Each point is one donor (ND 7, Aab+ 7, T1D 6); bars show medians. '+json.loads((DEST/'analysis_specification_before_fit.json').read_text())['quantity_cohort']+'. Cluster occupancy instead uses all complete core+20µm embedded objects per donor. Stars denote the omnibus Kruskal–Wallis comparison: raw p for the primary Beta-presence panel and BH-adjusted p across the 16 secondary endpoints for all other panels (*<0.05, **<0.01, ***<0.001). Exact omnibus and pairwise statistics are listed below. No result interpretation is assigned.\n'+tests[tests.endpoint.isin(keys)].to_csv(index=False)+'\nPairwise tests:\n'+pairs[pairs.endpoint.isin(keys)].to_csv(index=False)
        (DEST/'figures'/f'{stem}_legend.txt').write_text(legend)
    dots(['beta_presence_percent','beta_density_mm2','chga_hormone_negative_density_mm2'],'cell_quantity')
    dots(['beta_endocrine_percent','alpha_endocrine_percent','delta_endocrine_percent','chga_hormone_negative_endocrine_percent'],'endocrine_composition',2)
    dots([f'cluster_{cl}_occupancy_percent' for cl in range(10)],'cluster_occupancy',5)
    manifest=[state_views(d[d.eligible_plotted],[(c,LABELS[c.replace('_fraction','_percent')].replace('(%)','(fraction)')) for c in ['beta_endocrine_fraction','chga_hormone_negative_endocrine_fraction']]+[('beta_density_mm2',LABELS['beta_density_mm2'])],DEST,'endocrine_states','The shared core+20µm-complete embedded cohort is displayed; endpoints use corrected accepted-label counts.')]
    primary=tests[tests.family.eq('primary')].copy();primary['hypothesis']='H5';primary['comparison']='Disease difference in donor percentage of Beta-present objects';primary['test']='Kruskal-Wallis';primary['n_donors']=20;primary['n_objects']=int(d.core_complete.sum());primary['status']='computed; amended categorical endpoint'
    primary[['hypothesis','comparison','test','p_value','n_donors','n_objects','status']].to_csv(DEST/'primary_comparisons.csv',index=False)
    # Fixed donor choice: smallest numeric case ID in each disease, independent of outcomes.
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained');maprows=[]
    for ax,g in zip(axes,GROUPS):
        case=min(d.loc[d.donor_status.eq(g),'case_id'].unique(),key=int);z=d[d.case_id.eq(case)&d.core_complete].copy();maprows.append(z[['islet_id','case_id','donor_status','centroid_x_um','centroid_y_um','beta_presence','in_composition_embedding']])
        for status,col,lab in [(0,'#bdbdbd','Beta absent'),(1,COLORS[g],'Beta present')]:
            zz=z[z.beta_presence.eq(status)];ax.scatter(zz.centroid_x_um,zz.centroid_y_um,c=col,s=10,label=lab,rasterized=True)
        ax.set(title=f'{g} · Donor {case}',xlabel='X (µm)',ylabel='Y (µm)');ax.set_aspect('equal');ax.invert_yaxis();ax.legend(frameon=False,fontsize=8)
    savefig(fig,DEST/'figures','beta_presence_tissue_maps');pd.concat(maprows).to_csv(DEST/'data/beta_presence_tissue_maps.csv',index=False)
    (DEST/'figures/beta_presence_tissue_maps_legend.txt').write_text('One fixed donor per disease: smallest numeric donor ID, chosen independently of measurements. Points are complete-core paired endocrine objects at saved physical centroids, including objects outside the composition embedding. Color indicates at least one accepted Beta cell versus none. This is an object-centroid map, not a microscopy overlay. Exact points and eligibility are exported. No tests or interpretation are attached.\n')
    (DEST/'manifest.json').write_text(json.dumps({'state_figures':manifest,'primary':'primary_comparisons.csv','omnibus':'omnibus_tests.csv','pairwise':'pairwise_tests.csv','interpretation_added':False},indent=2)+'\n')

if __name__=='__main__':
    {'prepare':prepare,'figures':figures}[sys.argv[1]]()
