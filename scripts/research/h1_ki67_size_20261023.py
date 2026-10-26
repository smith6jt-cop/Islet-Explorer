#!/usr/bin/env python3
"""Source nuclear Ki67 localization by size: continuous signal, no positive calls."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("output/h1_immune_size_20261023")
BINS = ["<50 um", "50-<100 um", ">=100 um"]
TYPES = ["IMMUNE_UNION", "CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune"]
STATUS = ["ND", "Aab+", "T1D"]


def main():
    d = pd.read_parquet(ROOT/"ki67/ki67_immune_cells.parquet")
    d["case_id"] = d.case_id.astype(str)
    base = pd.read_csv("output/h1_20261022/h1_object_features.csv", dtype={"case_id": str})
    base["equivalent_diameter_um"] = np.sqrt(4*base.core_area_um2/np.pi)
    base["size_bin"] = pd.cut(base.equivalent_diameter_um, [0,50,100,np.inf], right=False, labels=BINS).astype(str)
    base_columns = ["case_id", "base_islet_id", "islet_id", "size_bin", "core_area_um2", "n_beta_core"]
    d = d[d.in_tissue_geometry & d.qc_analysis_eligible & d.audit_compartment.isin(["paired_core", "paired_peri_exclusive"])].copy()
    d["base_islet_id"] = np.where(d.audit_compartment == "paired_core", d.islet_name, d.exclusive_owner)
    d = d.merge(base[base_columns], on=["case_id", "base_islet_id"], how="left", validate="many_to_one")
    if d.islet_id.isna().any() or d.duplicated(["case_id", "object_id"]).any():
        raise ValueError("Ki67 unique cell to paired object join failed")
    d["compartment"] = np.where(d.audit_compartment == "paired_core", "core", "peri_exclusive")
    d.to_parquet(ROOT/"ki67/ki67_size_joined_cells.parquet", index=False)
    rows = []
    donors = base[["case_id", "donor_status"]].drop_duplicates()
    for donor in donors.itertuples(index=False):
        for size in BINS:
            for phenotype in TYPES:
                for compartment in ["core", "peri_exclusive"]:
                    sub = d[(d.case_id == donor.case_id) & (d.size_bin == size) & (d.compartment == compartment)]
                    if phenotype != "IMMUNE_UNION":
                        sub = sub[sub.phenotype == phenotype]
                    rows.append({"case_id": donor.case_id, "donor_status": donor.donor_status, "size_bin": size, "phenotype": phenotype, "compartment": compartment,
                        "cells": len(sub), "objects": sub.islet_id.nunique(),
                        "mean_donor_relative_ki67_rank": sub.ki67_donor_immune_midrank.mean(),
                        "median_donor_relative_ki67_rank": sub.ki67_donor_immune_midrank.median(),
                        "mean_common_scaled_ki67": sub.ki67_common_scaled.mean(),
                        "mean_nucleus_area_um2": sub.nucleus_area_um2.mean(),
                        "nucleus_area_lt10um2": int((sub.nucleus_area_um2 < 10).sum())})
    counts = pd.DataFrame(rows)
    counts.to_csv(ROOT/"ki67/ki67_donor_type_size_compartment.csv", index=False)
    pairs = counts.pivot(index=["case_id", "donor_status", "size_bin", "phenotype"], columns="compartment", values=["cells", "mean_donor_relative_ki67_rank"])
    pairs.columns = [f"{v}_{c}" for v,c in pairs.columns]
    pairs = pairs.reset_index()
    pairs["core_minus_peri_rank"] = pairs.mean_donor_relative_ki67_rank_core - pairs.mean_donor_relative_ki67_rank_peri_exclusive
    pairs.to_csv(ROOT/"ki67/ki67_paired_donor_size_differences.csv", index=False)
    summary = pairs.groupby(["donor_status", "size_bin", "phenotype"], sort=False).agg(paired_donors=("core_minus_peri_rank", "count"), equal_donor_mean_difference=("core_minus_peri_rank", "mean"), median_difference=("core_minus_peri_rank", "median"), core_cells=("cells_core", "sum"), peri_cells=("cells_peri_exclusive", "sum")).reset_index()
    summary.to_csv(ROOT/"ki67/ki67_size_summary.csv", index=False)
    fig, axes = plt.subplots(1,3,figsize=(11.5,4.6),sharey=True)
    colors = {"ND": "#277DA8", "Aab+": "#E69F00", "T1D": "#AA3377"}
    for ax, size in zip(axes,BINS):
        for i,status in enumerate(STATUS):
            part = pairs[(pairs.phenotype == "IMMUNE_UNION") & (pairs.size_bin == size) & (pairs.donor_status == status)]
            values = part.core_minus_peri_rank.dropna().to_numpy()
            offsets = np.linspace(-.1,.1,len(values)) if len(values)>1 else [0]*len(values)
            ax.scatter(i+np.array(offsets), values, s=42, color=colors[status], alpha=.85, edgecolor="white", lw=.4)
            if len(values):
                ax.plot([i-.18,i+.18], [np.mean(values)]*2, color="black",lw=2)
            ax.text(i,.65,f"n={len(values)}",ha="center",va="top",fontsize=9)
        ax.axhline(0,color="#999999",ls="--",lw=.8)
        ax.set_xticks(range(3),STATUS)
        ax.set_title(size.replace("um", "µm"),fontsize=12)
        ax.set_ylim(-.7,.7)
        ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("Mean donor-relative nuclear Ki67 rank\ncore − exclusive peri",fontsize=10)
    fig.suptitle("Source Ki67: continuous nuclear signal by object size",fontsize=14,y=.99)
    fig.text(.5,.02,"Each dot is a donor with measured immune cells in both compartments; bars are donor means.\nRanks describe relative signal within a donor. They do not identify proliferating cells. Peri ownership is exclusive.",ha="center",fontsize=9)
    fig.tight_layout(rect=[0,.13,1,.94])
    (ROOT/"figures").mkdir(exist_ok=True)
    for ext in ["png","svg"]:
        fig.savefig(ROOT/f"figures/figS3_ki67_size.{ext}",dpi=200,facecolor="white")
    plt.close(fig)
    audit = {"cells_analyzed": len(d), "core_cells": int((d.compartment == "core").sum()), "exclusive_peri_cells": int((d.compartment == "peri_exclusive").sum()), "donors": d.case_id.nunique(), "empty_populations_retained": True,
        "analysis": "Descriptive donor means and paired compartment differences by geometric size and phenotype; no inferential p or positivity threshold", "source": "ID-linked audit nuclear Ki67; within-donor rank reference all curated immune, before clipping", "limits": "No biological positivity calibration; current immune identity refinement includes Ki67; exclusive peri sensitivity differs from overlapping peri primary marker analysis; no across-donor global signal test from ranks"}
    (ROOT/"ki67/size_analysis.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(audit,indent=2))


if __name__ == "__main__":
    main()
