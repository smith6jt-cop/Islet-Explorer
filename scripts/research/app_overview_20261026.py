#!/usr/bin/env python3
from app_hypotheses_shared_20261026 import *

def main():
    folder=OUT/'shared';(folder/'figures').mkdir(exist_ok=True);(folder/'data').mkdir(exist_ok=True)
    d=read(folder/'object_measurements.parquet');q=d[d.eligible_plotted].copy()
    fig,axes=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
    for g in GROUPS:
        z=q[q.donor_status.eq(g)];axes[0].scatter(z.umap_1,z.umap_2,c=COLORS[g],s=7,label=g,rasterized=True)
    axes[0].legend(frameon=False,markerscale=2,fontsize=9);axes[0].set_title('Disease state')
    colors=plt.get_cmap('tab10')
    for cl in range(10):
        z=q[q['leiden_0.5'].eq(cl)];axes[1].scatter(z.umap_1,z.umap_2,c=[colors(cl)],s=7,label=str(cl),rasterized=True)
    axes[1].legend(frameon=False,ncol=5,markerscale=2,fontsize=8,title='Leiden 0.5',loc='upper center',bbox_to_anchor=(.5,-.17));axes[1].set_title('Saved composition clusters')
    sc=axes[2].scatter(q.umap_1,q.umap_2,c=q.pt_rank,s=7,cmap='viridis',vmin=0,vmax=1,rasterized=True);axes[2].set_title('App pseudotime');fig.colorbar(sc,ax=axes[2],label='Pseudotime (rank-normalized)',shrink=.8)
    for ax in axes:ax.set(xlabel='UMAP 1',ylabel='UMAP 2');ax.set_aspect('equal')
    savefig(fig,folder/'figures','app_representation_overview')
    q[['islet_id','case_id','donor_status','umap_1','umap_2','leiden_0.5','pseudotime','pt_rank','core_complete','neighborhood_complete_20']].to_csv(folder/'data/app_representation_overview.csv',index=False)
    counts=q.groupby('donor_status').agg(objects=('islet_id','size'),donors=('case_id','nunique')).to_dict('index')
    (folder/'figures/app_representation_overview_legend.txt').write_text('The same saved composition UMAP coordinates are shown by disease, saved Leiden resolution 0.5, and the current app DPT rank. The original embedding and runtime DPT fit contain 1,215 objects. The common display contains 1,126 complete-core and complete-20µm-neighborhood objects from 20 donors: '+json.dumps(counts)+'. No embedding, clusters or pseudotime ranks were refitted after display filtering. Endpoint models keep their own eligibility and include eligible objects outside this representation. DPT orders the app core/peri composition profiles and is not a measurement of elapsed time. No statistical test or biological interpretation is assigned to this overview.\n')

if __name__=='__main__':main()
