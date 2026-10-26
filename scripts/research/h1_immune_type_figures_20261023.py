#!/usr/bin/env python3
"""Presentation of frozen immune-type H1 size analyses; fits no models."""
from __future__ import annotations

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/islet-h1-type-matplotlib")
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator, PercentFormatter
import numpy as np
import pandas as pd

STATUS = ["ND", "Aab+", "T1D"]
SIZES = ["<50 um", "50-<100 um", ">=100 um"]
SIZE_LABEL = {"<50 um": "<50 µm", "50-<100 um": "50–<100 µm", ">=100 um": "≥100 µm"}
TYPES = ["CD8a T cell", "CD4 T cell", "Other T cell", "B cell", "Macrophage", "Unresolved immune"]
COLOR = {"Broad immune union": "#243447", "CD8a T cell": "#0072B2", "CD4 T cell": "#CC79A7", "Other T cell": "#56B4E9", "B cell": "#E69F00", "Macrophage": "#009E73", "Unresolved immune": "#999999"}


def save(fig, stem, output):
    for suffix in ["png", "svg"]:
        fig.savefig(output/f"{stem}.{suffix}", dpi=240, facecolor="white")
    plt.close(fig)


def summaries(counts):
    source = counts[counts.size_definition == "equivalent_diameter"]
    rows = []
    for size in SIZES:
        for status in STATUS:
            part = source[(source.size_bin == size) & (source.donor_status == status)]
            broad = part[part.immune_type == "Broad immune union"].set_index("case_id")
            positive = broad[broad.type_core_cells > 0].index
            assessable = broad[broad.base_eligible_objects > 0].index
            for immune_type in TYPES:
                typed = part[part.immune_type == immune_type].set_index("case_id")
                share = typed.loc[positive, "type_core_cells"] / broad.loc[positive, "type_core_cells"]
                rows.append({"size_bin": size, "donor_status": status, "immune_type": immune_type,
                    "assessable_donors": len(assessable), "core_event_donors": len(positive),
                    "donor_mean_share_core_positive_donors_only": float(share.mean()),
                    "donor_mean_type_core_density_mm2_including_zero_event_donors": float(typed.loc[assessable, "type_core_density_mm2"].mean()),
                    "core_cells_all_donors": int(typed.type_core_cells.sum()),
                    "broad_core_cells_all_donors": int(broad.type_core_cells.sum()),
                    "broad_core_positive_objects": int(broad.core_positive_objects.sum()),
                    "base_eligible_objects": int(broad.base_eligible_objects.sum()),
                    "donor_mean_broad_core_density_mm2": float(broad.loc[assessable, "type_core_density_mm2"].mean())})
    return pd.DataFrame(rows)


def figure_one(counts, summary, output):
    fig, axes = plt.subplots(2, 3, figsize=(14.6, 10.4), sharey="row")
    fig.subplots_adjust(left=.075, right=.985, top=.74, bottom=.19, wspace=.15, hspace=.34)
    fig.suptitle("Core immune types across endocrine object sizes", x=.075, y=.97, ha="left", fontsize=19, fontweight="bold", color="#172B3A")
    fig.text(.075, .928, "Equal donor weighting · geometric size groups · current 20-donor cohort", fontsize=12.5, color="#526370")
    handles = [Patch(facecolor=COLOR[t], label=t) for t in TYPES]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.068, .89), ncol=3, frameon=False, columnspacing=2.3, fontsize=11)
    rng = np.random.default_rng(20261023)
    for col, size in enumerate(SIZES):
        summ = summary[summary.size_bin == size]
        part = counts[(counts.size_definition == "equivalent_diameter") & (counts.size_bin == size)]
        for row, measure in enumerate(["donor_mean_share_core_positive_donors_only", "donor_mean_type_core_density_mm2_including_zero_event_donors"]):
            ax = axes[row, col]
            bottom = np.zeros(3)
            for immune_type in TYPES:
                values = np.array([summ[(summ.donor_status == s) & (summ.immune_type == immune_type)][measure].iloc[0] for s in STATUS])
                ax.bar(np.arange(3), values, bottom=bottom, width=.62, color=COLOR[immune_type], edgecolor="white", linewidth=.45, zorder=3)
                bottom += np.nan_to_num(values)
            ax.set_xticks(np.arange(3), STATUS)
            ax.tick_params(axis="x", length=0, labelsize=12)
            ax.set_xlim(-.65, 2.65)
            ax.grid(axis="y", color="#D7DDE3", linewidth=.6, zorder=0)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
            if row == 0:
                ax.set_ylim(0, 1.115)
                ax.yaxis.set_major_locator(FixedLocator([0, .25, .5, .75, 1]))
                ax.yaxis.set_major_formatter(PercentFormatter(1))
                ax.set_title(SIZE_LABEL[size], fontsize=15, pad=24, fontweight="bold")
                for j, status in enumerate(STATUS):
                    info = summ[summ.donor_status == status].iloc[0]
                    ax.text(j, 1.035, f"{info.core_event_donors}/{info.assessable_donors} donors", ha="center", va="bottom", fontsize=10, color="#34495A")
            else:
                ax.set_ylim(-12, 650)
                ax.set_yticks([0, 200, 400, 600])
                for j, status in enumerate(STATUS):
                    broad = part[(part.donor_status == status) & (part.immune_type == "Broad immune union")]
                    broad = broad[broad.base_eligible_objects > 0].sort_values("case_id")
                    jitter = rng.uniform(-.18, .18, len(broad))
                    ax.scatter(j+jitter, broad.type_core_density_mm2, facecolor="white", edgecolor="#172B3A", s=30, linewidth=.9, zorder=6)
            if col == 0:
                ax.set_ylabel("Mean donor composition\namong core-positive donors" if row == 0 else "Core immune density (cells/mm²)\nMean donor density and individual donors", fontsize=11.5)
    fig.text(.075, .145, "Composition is undefined for donors with no core immune cells in that size bin; labels show core-positive / assessable donors.", fontsize=10.7, color="#34495A")
    fig.text(.075, .116, "Density includes zero-event donors. Stacked density bars are donor means; hollow points are each donor’s total core immune density.", fontsize=10.7, color="#34495A")
    fig.text(.075, .087, "Within each donor, density = summed core cells / summed endocrine core area. It is not whole-tissue density or a clinical progression measure.", fontsize=10.7, color="#34495A")
    fig.text(.075, .058, "Counts include all base-eligible objects, including those with no nearby immune cells. Islet/Cluster identities are retained in the adjusted models.", fontsize=10.7, color="#34495A")
    save(fig, "fig1_immune_types_size", output)


def figure_two(estimates, output):
    e = estimates[(estimates.size_definition == "equivalent_diameter") & (estimates.cohort == "All")]
    labels = ["Broad immune union"]+TYPES
    y = dict(zip(labels, np.arange(len(labels))[::-1]))
    fig = plt.figure(figsize=(18.4, 9.4))
    grid = fig.add_gridspec(1, 6, width_ratios=[2.05, 1.25]*3, left=.16, right=.987, top=.745, bottom=.26, wspace=.11)
    fig.suptitle("Immune localization and Beta-cell density by size and immune type", x=.03, y=.97, ha="left", fontsize=20, fontweight="bold", color="#172B3A")
    fig.text(.03, .925, "All-cohort isolated fits · density ratio per 10 percentage points of type-specific allocation inside the core", fontsize=12.5, color="#526370")
    fig.text(.03, .875, "No immune subtype meets its primary multiplicity threshold. Unsupported slots are labelled ‘Descriptive only’.", fontsize=12.5, color="#34495A")
    for col, size in enumerate(SIZES):
        ax = fig.add_subplot(grid[0, col*2])
        support = fig.add_subplot(grid[0, col*2+1], sharey=ax)
        ax.set_xscale("log")
        ax.set_xlim(.66, 1.4)
        ax.set_ylim(-.6, 6.6)
        ax.axvline(1, color="#87939E", linestyle="--", linewidth=1, zorder=1)
        ax.set_yticks(np.arange(7))
        if col == 0:
            ax.set_yticklabels(labels[::-1], fontsize=11.5)
        else:
            ax.tick_params(labelleft=False)
        ax.tick_params(axis="y", length=0, pad=10)
        ticks = [.7, .8, 1, 1.2, 1.4]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter(["0.70", "0.80", "1.00", "1.20", "1.40"]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis="x", labelsize=10)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_title(SIZE_LABEL[size], fontsize=15.5, fontweight="bold", pad=24)
        ax.set_xlabel("Beta-density ratio", fontsize=11, labelpad=9)
        support.axis("off")
        support.set_xlim(0, 1)
        support.text(.03, 6.76, "N / core+ / D*", fontsize=10.5, fontweight="bold", color="#34495A")
        for immune_type in labels:
            r = e[(e.size_bin == size) & (e.immune_type == immune_type)].iloc[0]
            yi = y[immune_type]
            if yi % 2 == 0:
                ax.axhspan(yi-.47, yi+.47, color="#EDF1F5", alpha=.6, linewidth=0, zorder=0)
                support.axhspan(yi-.47, yi+.47, color="#EDF1F5", alpha=.6, linewidth=0, zorder=0)
            if r.status == "FITTED_EXPLORATORY":
                color = COLOR[immune_type]
                ax.errorbar(r.rr, yi, xerr=[[r.rr-r.lower], [r.upper-r.rr]], fmt="o", color=color, markersize=7.5, capsize=3.5, elinewidth=1.6, zorder=5)
                support.text(.03, yi+.12, f"{r.objects:,} / {r.core_positive_objects} / {r.informative_donors}", fontsize=10.1, va="center", color="#243447")
                support.text(.03, yi-.19, f"q = {r.q_bh:.3f}", fontsize=9.4, va="center", color="#526370")
            else:
                ax.text(1, yi, "Descriptive only", ha="center", va="center", fontsize=10, color="#8B959E")
                support.text(.03, yi, f"{r.objects:,} / {r.core_positive_objects} / {r.informative_donors}", fontsize=10.1, va="center", color="#8B959E")
    fig.text(.03, .18, "Support: N = eligible type-exposed objects after excluding all-zero-Beta donors; core+ = objects containing this immune type; D* = donors with allocation variation.", fontsize=10.6, color="#34495A")
    fig.text(.03, .145, "Inference requires ≥100 objects, ≥30 core-positive objects, and ≥6 informative donors. Intervals use donor-cluster uncertainty and a t reference with D* − 1 degrees of freedom.", fontsize=10.6, color="#34495A")
    fig.text(.03, .11, "BH families: 18 subtype × size slots; 3 broad-union size slots separately. Untestable slots contribute 1 internally; no p or CI is displayed for them.", fontsize=10.6, color="#34495A")
    fig.text(.03, .075, "No T1D-only geometric-size fit passes the support gate. This panel combines all donor groups; it does not establish a T1D-specific immune mechanism.", fontsize=10.6, color="#34495A")
    fig.text(.03, .04, "Size is equivalent core diameter from area. The endocrine-count sensitivity is separate because endocrine count includes the Beta endpoint. Cells and objects are not independent donors.", fontsize=10.6, color="#34495A")
    save(fig, "fig2_type_effects_size", output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--types-root", type=Path, default=Path("output/h1_immune_size_20261023/types"))
    parser.add_argument("--output", type=Path, default=Path("output/h1_immune_size_20261023/figures"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
    counts = pd.read_csv(args.types_root/"donor_type_size_counts.csv")
    estimates = pd.read_csv(args.types_root/"estimates.csv")
    summary = summaries(counts)
    summary.to_csv(args.types_root/"figure1_donor_composition_density_summary.csv", index=False)
    figure_one(counts, summary, args.output)
    figure_two(estimates, args.output)
    findings = """# Immune-type and object-size findings

The primary geometric-size analyses do not identify an immune subtype whose localization–Beta-density association survives the prespecified 18-slot BH family. Only 5 of the 18 all-cohort subtype slots have enough exposure support to fit; the remaining slots are descriptive. CD4 T cells, Other T cells and B cells are too sparse throughout the geometric size strata for these adjusted tests. The broad-union >=100 micrometer estimate has a nominal interval barely below 1 (RR 0.958, 95% CI 0.918-1.000; nominal p=0.0493), but its separate 3-slot BH q=0.148. No T1D-only geometric-size fit meets the unchanged support gates.

Core immune composition differs across donors and groups. Macrophage and unresolved immune phenotypes account for most core immune cells in the donor-balanced summaries; the broad union is not a cytotoxic T-cell assay. ND cores contain no typed CD8a or CD4 T cells in these eligible data. For objects >=100 micrometers, the mean donor core-immune CD8a fraction is 19.5% in T1D, 3.6% in Aab+, and 0% in ND. These fractions are descriptive identities, not proof of cytotoxicity or antigen-specificity; unresolved immune is an operational residual phenotype.

The 50-<100 micrometer bin has one donor without a core immune event in each group. Their composition is undefined, and those donors are excluded only from the composition average; their zero densities remain in the absolute-density summary. All donors have some core immune events in the other two diameter bins. Mean broad core immune density in ND/Aab+/T1D is respectively 92.4/197.2/259.9 cells/mm2 for <50 micrometers, 14.3/34.2/78.2 for 50-<100 micrometers, and 44.4/69.1/173.3 for >=100 micrometers. These are equal-donor summaries of cell count divided by sampled endocrine core area, not whole-pancreas immune density. Small cores are especially sensitive to single-cell occupancy and boundary placement.

The >=100 micrometer CD8a fit is positive (RR 1.033, CI 0.980-1.089), with only 8 informative donors and 32 core-positive objects; no subtype-specific inference is established (q=1.000). Its positive point estimate persists in the available leave-one-donor-out fits, but concentration in a small subset of donors and phenotype/size selection remain limitations. All six T1D donors are not six equally informative donors for every type and size.

The endocrine-count sensitivity has some smaller nominal p values, including a positive large-endocrine-count T1D broad-immune association. It conditions on a quantity that includes the Beta endpoint and hormone-based typing, so it is not an independent confirmation and cannot replace the geometric primary analysis.

The data contain one acquisition per donor and use current hormone/immune typing plus two-dimensional saved boundaries. Models describe present occupancy, not recruitment, cumulative injury, killing or treatment targets. All-zero-Beta donors are listed and excluded only from their corresponding mean fits; source and descriptive tables retain their observations. Object-associated peri and field counts can overlap and must not be read as unique tissue totals. All 168 size/type/cohort slots remain in estimates.csv, with no fabricated p values for the 136 descriptive-only slots.

Figures are generated without refitting from estimates.csv and donor_type_size_counts.csv. Figure 1 supplies donor-balanced composition and density with individual donor totals. Figure 2 supplies every geometric all-cohort slot, support counts and the declared BH q values. Editable SVG companions accompany the PNG files.
"""
    (args.types_root/"findings.md").write_text(findings)
    audit = {"input_sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.types_root/"donor_type_size_counts.csv", args.types_root/"estimates.csv"]},
             "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "output_sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.glob("fig[12]_*.??*")},
             "software": {"matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__}}
    (args.types_root/"figure_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print("Saved fig1_immune_types_size and fig2_type_effects_size as PNG/SVG, plus findings.md and numeric figure summary.")


if __name__ == "__main__":
    main()
