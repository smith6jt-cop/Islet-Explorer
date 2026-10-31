#!/usr/bin/env python3
"""Bind corrected endpoints to existing app representations; no embedding fits."""
from pathlib import Path
import hashlib,json,os,textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/app_hypotheses_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/app_hypotheses_20261026'
GROUPS=['ND','Aab+','T1D']
COLORS={'ND':'#4477AA','Aab+':'#CC6633','T1D':'#228833'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def read(p):
    p=Path(p)
    d=pd.read_parquet(p) if p.suffix=='.parquet' else pd.read_csv(p,dtype={'case_id':str,'islet_id':str})
    if 'case_id' in d:d['case_id']=d.case_id.astype(str)
    return d

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def build():
    out=OUT/'shared';out.mkdir(parents=True,exist_ok=True)
    source=ROOT/'output/h1_edge_corrected_20261024/h1_object_features.parquet'
    rep=ROOT/'output/pancfloyd_20donor_20261026/trajectory/data/object_representation.csv'
    hood=ROOT/'data/development/uniform-20donor-20261022-v10/neighborhood_metrics.csv'
    base=read(source);a=read(rep)
    assert not base.islet_id.duplicated().any() and not a.islet_id.duplicated().any()
    keep=['islet_id','umap_1','umap_2','leiden_0.5','pseudotime','pt_rank','in_composition_embedding','eligible_plotted']
    base=base.merge(a[keep],on='islet_id',how='left',validate='one_to_one')
    for c in ['in_composition_embedding','eligible_plotted']:base[c]=base[c].fillna(False).astype(bool)
    base['leiden_0.5']=pd.to_numeric(base['leiden_0.5']).astype('Int64')
    n=read(hood).drop(columns=['imageid','object_class']);n=n.rename(columns={c:'app_'+c for c in n if c!='islet_id'})
    base=base.merge(n,on='islet_id',how='left',validate='one_to_one')
    inputs=[source,rep,hood]
    for name,group,valuecols in [
      ('structural_cells/object_structure_counts.csv','structure',['n_recovered_cells','count_fraction_percent','density_cells_mm2']),
      ('structural_cells/object_passthrough_summaries.csv','structure',['n_centroid_inside','n_passthrough_measured','mean_passthrough_fraction_among_measured']),
      ('slides1_4/data/structure_union_area_values.csv','class',['area_um2','area_percent','intersecting_annotation_count'])]:
        p=ROOT/'output/pancfloyd_20donor_20261026'/name;inputs.append(p);d=read(p)
        d['field']=d[group].str.lower().str.replace(' ','_',regex=False)+'_'+d.region
        for v in valuecols:
            w=d.pivot(index='islet_id',columns='field',values=v);w.columns=['structure_'+str(c)+'_'+v for c in w]
            base=base.merge(w,on='islet_id',how='left',validate='one_to_one')
    base['beta_density_mm2']=base.n_beta_core*1e6/base.core_area_um2
    base['beta_presence']=base.n_beta_core.gt(0).astype(int)
    for cell,col in [('beta','n_core_Beta_cell'),('alpha','n_core_Alpha_cell'),('delta','n_core_Delta_cell'),('chga_hormone_negative','n_core_CHGAplus_hormone-negative')]:
        base[cell+'_endocrine_fraction']=base[col]/base.n_endocrine_core.replace(0,np.nan)
        base[cell+'_density_mm2']=base[col]*1e6/base.core_area_um2
    base['immune_density_core_mm2']=base.n_immune_core*1e6/base.core_area_um2
    base['immune_density_peri_mm2']=base.n_immune_peri*1e6/base.peri_area_um2.replace(0,np.nan)
    base.to_parquet(out/'object_measurements.parquet',index=False)
    base.to_csv(out/'object_measurements.csv',index=False)
    summary=base.groupby(['case_id','donor_status']).agg(n_paired=('islet_id','size'),n_embedding=('in_composition_embedding','sum'),n_complete_core=('core_complete','sum'),n_complete_core_peri=('eligible_plotted','sum')).reset_index()
    summary.to_csv(out/'donor_population_counts.csv',index=False)
    meta={'sources':{str(p.relative_to(ROOT)):digest(p) for p in inputs},'n_objects':len(base),'n_embedded':int(base.in_composition_embedding.sum()),'n_display_complete':int(base.eligible_plotted.sum()),'n_donors':base.case_id.nunique(),'no_embedding_refit':True,'no_rerank':True,'served_pipeline_changed':False,'neighborhood_names':'app_ prefix retains exact served measures; unprefixed counts/densities use corrected research cell universe','structural_labels':'Recovered structural labels remain separate from original app phenotypes; no double addition to phenotype totals'}
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    assert (len(base),base.case_id.nunique(),int(base.in_composition_embedding.sum()),int(base.eligible_plotted.sum()))==(4763,20,1215,1126)
    print(json.dumps({k:v for k,v in meta.items() if k!='sources'}))

def savefig(fig,out,stem):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:fig.savefig(out/f'{stem}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def state_views(d,features,out,stem,population_note=''):
    """features is a list of (column, axis label); caller supplies its eligible cohort.

    Displays use existing UMAP and rank DPT. Twenty bins show arithmetic object
    means, as in the app heatmap. Donor-bin means and donor-cluster means are
    also exported. No tests are performed here.
    """
    out=Path(out);(out/'data').mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    q=d.loc[d.in_composition_embedding.fillna(False).astype(bool)&d.umap_1.notna()&d.pt_rank.notna()].copy()
    q['pt_bin']=np.minimum((q.pt_rank*20).astype(int),19);q['pt_midpoint']=(q.pt_bin+.5)/20
    fields=[f for f,l in features]
    q[['islet_id','case_id','donor_status','leiden_0.5','umap_1','umap_2','pseudotime','pt_rank','pt_bin',*fields]].to_csv(out/'data'/f'{stem}_state_objects.csv',index=False)
    profiles=[];donorbins=[];clusters=[]
    fig,axes=plt.subplots(len(features),4,figsize=(16,3.6*len(features)),squeeze=False,layout='constrained')
    for j,(field,label) in enumerate(features):
        v=pd.to_numeric(q[field],errors='coerce');finite=v[np.isfinite(v)]
        if not len(finite):continue
        lo,hi=finite.min(),finite.max()
        if lo==hi:hi=lo+1
        for k,g in enumerate(GROUPS):
            z=q[q.donor_status.eq(g)];ok=z[field].notna()
            axes[j,k].scatter(q.umap_1,q.umap_2,s=3,c='#eeeeee',rasterized=True)
            sc=axes[j,k].scatter(z.loc[ok,'umap_1'],z.loc[ok,'umap_2'],c=z.loc[ok,field],s=9,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
            axes[j,k].set_title(g);axes[j,k].set_xlabel('UMAP 1');axes[j,k].set_ylabel('UMAP 2' if k==0 else '')
            axes[j,k].set_aspect('equal',adjustable='box')
        wrapped=textwrap.fill(label,width=43)
        fig.colorbar(sc,ax=list(axes[j,:3]),label=wrapped,shrink=.8)
        b=q.groupby(['donor_status','pt_bin','pt_midpoint'],observed=True)[field].agg(['mean','count']).reset_index();b['feature']=field;profiles.append(b)
        db=q.groupby(['case_id','donor_status','pt_bin','pt_midpoint'],observed=True)[field].agg(['mean','count']).reset_index();db['feature']=field;donorbins.append(db)
        cl=q.groupby(['case_id','donor_status','leiden_0.5'],observed=True)[field].agg(['mean','count']).reset_index();cl['feature']=field;clusters.append(cl)
        ax=axes[j,3]
        for g in GROUPS:
            z=b[b.donor_status.eq(g)&b['count'].gt(0)];ax.plot(z.pt_midpoint,z['mean'],'o-',ms=3,lw=1,color=COLORS[g],label=g)
        ax.set(xlabel='Pseudotime (rank-normalized)',ylabel=wrapped,xlim=(0,1));ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.2)
    savefig(fig,out/'figures',stem+'_states')
    pd.concat(profiles,ignore_index=True).to_csv(out/'data'/f'{stem}_pseudotime_profiles.csv',index=False)
    pd.concat(donorbins,ignore_index=True).to_csv(out/'data'/f'{stem}_donor_pseudotime.csv',index=False)
    cd=pd.concat(clusters,ignore_index=True);cd.to_csv(out/'data'/f'{stem}_donor_cluster_values.csv',index=False)
    fig,axes=plt.subplots(len(features),1,figsize=(11,2.7*len(features)),squeeze=False,layout='constrained')
    for j,(field,label) in enumerate(features):
        z=cd[cd.feature.eq(field)].groupby(['donor_status','leiden_0.5'],observed=True)['mean'].mean().unstack().reindex(index=GROUPS,columns=range(10))
        im=axes[j,0].imshow(z,aspect='auto',cmap='viridis');axes[j,0].set(yticks=range(3),yticklabels=GROUPS,xticks=range(10),xlabel='Saved composition cluster (Leiden 0.5)',title=label);fig.colorbar(im,ax=axes[j,0],label='Mean of donor means')
    savefig(fig,out/'figures',stem+'_clusters')
    counts=q.groupby('donor_status').agg(objects=('islet_id','size'),donors=('case_id','nunique')).to_dict('index')
    legend=f'{population_note} Saved composition UMAP and Leiden resolution 0.5 are reused without refitting. Raw and rank app DPT are retained without reranking. Colored maps show the stated object measurement on a common scale across disease groups; gray points are the matched cohort background. Pseudotime curves are arithmetic object means in 20 fixed rank bins; bins with no measured values are absent. Cluster heatmaps show equal-donor means within each disease/cluster cell, with unobserved cells missing. Individual object values, donor-bin means and donor-cluster means are exported. These representation panels contain no hypothesis tests. Populations: {json.dumps(counts)}. Feature-specific missing values retain their own denominators in the exported count columns. No biological interpretation or conditional-result branch is assigned.'
    (out/'figures'/f'{stem}_states_legend.txt').write_text(legend+'\n');(out/'figures'/f'{stem}_clusters_legend.txt').write_text(legend+'\n')
    return {'stem':stem,'figures':[stem+'_states',stem+'_clusters'],'objects':len(q),'donors':q.case_id.nunique(),'features':features}

if __name__=='__main__':build()
