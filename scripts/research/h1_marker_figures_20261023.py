#!/usr/bin/env python3
"""Render audited marker estimates and transparent secondary diagnostics."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

ROOT=Path("output/h1_immune_size_20261023")
MARKERS=["CD45","CD3e","CD8a","CD4","CD20","CD68","CD163","HLADR","CD35","LGALS3","CD44","VIM"]
BINS=["<50 um","50-<100 um",">=100 um"]


def save(fig,name):
    (ROOT/"figures").mkdir(exist_ok=True)
    for ext in ["png","svg"]:
        fig.savefig(ROOT/f"figures/{name}.{ext}",dpi=200,facecolor="white")
    plt.close(fig)


def primary(e):
    d=e[e.family=="primary"]
    fig,axes=plt.subplots(1,3,figsize=(14,8),sharey=True)
    for ax,size in zip(axes,BINS):
        a=d[d.size_bin==size].set_index("marker").reindex(MARKERS)
        for y,(marker,row) in enumerate(a.iterrows()):
            if y%2==0:
                ax.axhspan(y-.46,y+.46,color="#f3f5f7",zorder=0)
            if pd.notna(row.rr):
                ax.errorbar(row.rr,y,xerr=[[row.rr-row.lower],[row.upper-row.rr]],fmt="o",color="#176B89",markersize=5,capsize=2,lw=1.5)
                ax.text(1.02,y,f"q={row.q_bh:.2f}",transform=ax.get_yaxis_transform(),va="center",fontsize=8,color="#24384A")
            else:
                ax.text(.92,y,"Insufficient support",ha="center",va="center",fontsize=8,color="#92999f")
            ax.text(1.02,y+.28,f"{int(row.core_positive_objects)} / {int(row.informative_donors)}",transform=ax.get_yaxis_transform(),va="center",fontsize=7,color="#6b747c")
        ax.axvline(1,color="#92999f",ls="--",lw=1)
        ax.set_xscale("log")
        ax.set_xlim(.60,1.20)
        ax.xaxis.set_major_locator(FixedLocator([.6,.8,1,1.2]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x:.1f}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(size.replace("um","µm"),fontsize=13,pad=17)
        ax.set_xlabel("Beta-density ratio per 10 points",fontsize=10,labelpad=9)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y",length=0)
        ax.text(1.02,1.02,"BH q\ncore+ / D*",transform=ax.transAxes,ha="left",va="bottom",fontsize=8,color="#6b747c")
    axes[0].set_yticks(range(12),MARKERS,fontsize=10)
    axes[0].set_ylim(11.6,-.6)
    fig.suptitle("Marker-positive immune localization and Beta density by size",fontsize=16,x=.08,ha="left",y=.98)
    fig.text(.08,.925,"All 12 immune-associated/shared-state markers; no primary association survives BH correction (minimum q = 0.425).",fontsize=11,color="#425466")
    fig.text(.08,.045,"Marker positivity is measured within curated immune cells at scaled fluorescence ≥0.10; it does not redefine cell identity.\nDots: adjusted ratios; lines: donor-cluster 95% intervals. D*: donors with allocation variation; core+: core-positive objects.\nBH family: 36 marker × size slots. Descriptive-only slots have no p/CI. All groups pooled; no T1D-only stratum meets support gates.",fontsize=9,color="#425466",linespacing=1.6)
    fig.subplots_adjust(left=.09,right=.90,bottom=.19,top=.84,wspace=.60)
    save(fig,"fig4_marker_effects_size")


def secondary(e):
    ms=["CD3e","CD8a","LGALS3"]
    boot=pd.read_csv(ROOT/"marker_models/conditional_bootstrap_summary.csv").set_index("marker")
    leave=pd.read_csv(ROOT/"marker_models/conditional_leave_one_donor_out.csv")
    threshold=pd.read_csv(ROOT/"marker_models/conditional_threshold_point_checks.csv")
    fig,axes=plt.subplots(1,3,figsize=(13.5,5.8),sharey=True)
    for ax,marker in zip(axes,ms):
        a=e[(e.marker==marker)&(e.size_bin=="<50 um")&(e.threshold==.1)&(e.cohort=="All")].set_index("model")
        parallel=a.loc["parallel"]
        cond=a.loc["conditional_broad"]
        for y,row,color in [(0,parallel,"#7B8790"),(1,cond,"#176B89")]:
            ax.errorbar(row.rr,y,xerr=[[row.rr-row.lower],[row.upper-row.rr]],fmt="o",color=color,capsize=3,lw=1.6,ms=6)
        b=boot.loc[marker]
        ax.plot([b.percentile_lower,b.percentile_upper],[2,2],color="#176B89",lw=2)
        ax.scatter([cond.rr],[2],s=38,facecolor="white",edgecolor="#176B89",zorder=3)
        sub=leave[(leave.marker==marker)&(leave.size_bin=="<50 um")]
        ys=3+np.linspace(-.10,.10,len(sub))
        ax.scatter(sub.rr,ys,s=19,color="#6F7D86",alpha=.7)
        t=threshold[(threshold.marker==marker)&(threshold.size_bin=="<50 um")].sort_values("threshold")
        xs=[t.iloc[0].rr,cond.rr,t.iloc[1].rr]
        ax.scatter(xs,[4,4,4],s=35,c=["#D99435","#176B89","#A36DA8"],zorder=3)
        ax.axvline(1,color="#92999f",ls="--",lw=1)
        ax.set_xscale("log")
        ax.set_xlim(.48,2.8)
        ax.xaxis.set_major_locator(FixedLocator([.5,.75,1,1.5,2.5]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f"{x:g}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(f"{marker}\nConditional q={cond.q_bh:.3f}",fontsize=12,pad=13)
        ax.set_xlabel("Beta-density ratio",fontsize=10)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y",length=0)
    axes[0].set_yticks(range(5),["Parallel marker model\n95% cluster interval","Conditional on broad location\n95% cluster interval","Conditional donor bootstrap\n95% percentile interval","Conditional leave-one-donor-out\npoint estimates","Conditional threshold check\npoint estimates"],fontsize=9)
    axes[0].set_ylim(4.6,-.6)
    fig.suptitle("Secondary signals in objects <50 µm: composition, not an entry effect",fontsize=15,y=.98)
    fig.text(.5,.895,"Three conditional-model candidates selected by q<0.05 from a separate 36-slot family; diagnostic checks were added after fitting.",ha="center",fontsize=10,color="#425466")
    fig.text(.48,.14,"Threshold colors:",ha="right",fontsize=9,color="#425466")
    for x,label,color in [(.50,"● 0.05","#D99435"),(.56,"● 0.10","#176B89"),(.62,"● 0.20","#A36DA8")]:
        fig.text(x,.14,label,ha="left",fontsize=9,color=color)
    fig.text(.5,.10,"Successful bootstrap fits /999: CD3e 995; CD8a 992; LGALS3 997. Failed draws remain recorded.",ha="center",fontsize=9,color="#425466")
    fig.text(.5,.025,"Holding broad core allocation and nearby marker+/marker− pools fixed, the contrast exchanges marker+ and marker− immune localization.\nInformation is concentrated in few donors. These pooled, cross-sectional candidates do not establish a T1D-specific mechanism or a size interaction.",ha="center",fontsize=9,color="#425466",linespacing=1.5)
    fig.subplots_adjust(left=.23,right=.98,bottom=.27,top=.76,wspace=.27)
    save(fig,"fig5_secondary_marker_diagnostics")


if __name__=="__main__":
    estimates=pd.read_csv(ROOT/"marker_models/estimates.csv")
    primary(estimates)
    secondary(estimates)
