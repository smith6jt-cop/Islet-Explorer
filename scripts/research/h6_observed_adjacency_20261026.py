#!/usr/bin/env python3
"""Descriptive observed endocrine adjacency; no label permutations or tests."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h6-observed-20261026-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('output/h5_h10_original_plan_20261026/h6')
GROUPS=['ND','Aab+','T1D']
COLORS=['#4477AA','#CC6633','#228833']

def main():
    graph=pd.read_csv(OUT/'graph_inventory.csv')
    features=pd.read_parquet(OUT.parent/'h5/object_ins_features.parquet')
    d=graph.merge(features[['islet_id','ins_mean_beta_core','h5_primary_eligible']],
                  on='islet_id',validate='one_to_one')
    d.to_csv(OUT/'observed_adjacency_object_features.csv',index=False)
    rows=[]
    count_columns=['n_endocrine_nodes','n_beta_nodes','n_alpha_nodes','n_delta_nodes',
                   'n_edges','observed_beta_alpha_edges','observed_beta_delta_edges']
    for (case,state),all_objects in d.groupby(['case_id','donor_status']):
        eligible=all_objects.loc[all_objects.graph_eligible].copy()
        row={'case_id':case,'donor_status':state,'complete_endocrine_objects':len(all_objects),
             'eligible_endocrine_objects':len(eligible),
             'eligible_objects_with_beta_INS':int(eligible.ins_mean_beta_core.notna().sum()),
             'eligible_objects_with_nearby_immune':int(eligible.h5_primary_eligible.sum()),
             'mean_observed_combined_fraction':eligible.observed_combined_fraction.mean()}
        row.update({col:int(eligible[col].sum()) for col in count_columns})
        row['edge_weighted_observed_combined_fraction']=(
            (row['observed_beta_alpha_edges']+row['observed_beta_delta_edges'])/row['n_edges']
            if row['n_edges'] else np.nan)
        rows.append(row)
    donor=pd.DataFrame(rows)
    donor.to_csv(OUT/'plot_data_observed_adjacency.csv',index=False)
    plt.rcParams.update({'font.size':11,'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(6.4,5.1))
    samples=[100*donor.loc[donor.donor_status.eq(g),'mean_observed_combined_fraction'].dropna().to_numpy()
             for g in GROUPS]
    bp=ax.boxplot(samples,positions=range(3),widths=.48,patch_artist=True,showfliers=False,
                  medianprops={'color':'#333333','lw':1.8})
    for i,(vals,col,box) in enumerate(zip(samples,COLORS,bp['boxes'])):
        box.set(facecolor=col,alpha=.18,edgecolor=col)
        ax.scatter(i+np.linspace(-.11,.11,len(vals)),vals,s=42,c=col,edgecolors='white',linewidths=.7,zorder=3)
    ax.set_xticks(range(3),GROUPS)
    ax.set_xlabel('Disease state')
    ax.set_ylabel('β–α/δ edges among endocrine edges (%)')
    ax.set_title('Observed α/δ adjacency to β cells')
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'observed_beta_alpha_delta_adjacency.{ext}',dpi=190)
    plt.close(fig)
    print(donor.to_string(index=False))

if __name__=='__main__':main()
