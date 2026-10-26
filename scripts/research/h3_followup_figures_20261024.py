#!/usr/bin/env python3
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/h234-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from h1_followups_20261024 import OUT,h1


def run():
    d=pd.read_csv(OUT/'h3_donor_effects.csv')
    d=d.loc[d.radius_um.eq(10)].sort_values(['donor_status','case_id'])
    fig,axs=plt.subplots(1,2,figsize=(11,5.5),gridspec_kw={'width_ratios':[1,1.3]})
    ax=axs[0]
    for i,(_,r) in enumerate(d.iterrows()):
        color=h1.COLORS[r.donor_status]
        ax.plot(100*r.excess,i,'o',ms=8,color=color,
                markerfacecolor=color if r.current_mask_accepted else 'white',markeredgewidth=1.5)
    ax.set_yticks(range(len(d)),[f'{r.donor_status} · {r.case_id}' for r in d.itertuples()])
    ax.set_title('Each point is one donor',fontsize=12)
    ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color='#555555',label='Reviewed current mask'),
                       Line2D([],[],marker='o',linestyle='',color='#555555',markerfacecolor='white',
                              label='Reviewed retained alternative')],loc='lower right',fontsize=8)
    ax.set_ylim(-2,len(d)-.4)
    s=pd.read_csv(OUT/'h3_summary.csv')
    ax=axs[1]
    specifications=[('current_reviewed_masks','Pooled','Current masks: 5 donors'),
                    ('all_retained_reviewed_masks','Pooled','All reviewed masks: 10 donors'),
                    ('all_retained_reviewed_masks','T1D','T1D: 3 donors')]
    for i,(selection,group,label) in enumerate(specifications):
        r=s.loc[s.selection.eq(selection)&s.group.eq(group)&s.radius_um.eq(10)].iloc[0]
        color=h1.COLORS['T1D'] if group=='T1D' else '#444444'
        ax.errorbar(100*r.excess,2-i,xerr=[[100*(r.excess-r.lower)],[100*(r.upper-r.excess)]],
                    fmt='o',color=color,capsize=4)
        ax.text(.03,2-i+.22,label,transform=ax.get_yaxis_transform(),fontsize=10)
        ax.text(.03,2-i-.25,f'{100*r.excess:+.1f} pp  [95% CI {100*r.lower:+.1f}, {100*r.upper:+.1f}]',
                transform=ax.get_yaxis_transform(),fontsize=10)
    ax.set_yticks([]);ax.set_ylim(-.6,2.65);ax.set_title('Equal-donor means and donor t intervals',fontsize=12)
    for ax in axs:
        ax.axvline(0,color='#999999',ls='--')
        ax.set_xlabel('Observed − expected contact (percentage points)')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('H3: immune contact within 10 µm of capillary boundary crossings',fontsize=14)
    fig.text(.5,.025,'999 rigid rotations per object preserve clustering and radial position within the same tissue and core.\n'
             'Evaluated: 89 objects / 164 core immune cells; results do not establish absence of a vessel effect.',
             ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.11,1,.95))
    fig.savefig(OUT/'h3_capillary_crossings.png',dpi=180)
    fig.savefig(OUT/'h3_capillary_crossings.pdf')


if __name__=='__main__':run()
