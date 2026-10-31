#!/usr/bin/env python3
"""Recreate the assigned PancFloyd neighborhood panels from corrected v10 data.

No marker-intensity quantification or phenotype reassignment. The six enrichment
panels retain the donor-level Kruskal-Wallis test used by the reference app.
The app's published-data geometry, typed cells, and neighborhood definitions are
reused; accepted cell and object edge flags are applied explicitly.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from scipy import stats

from h1_features_20261022 import geometry

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/development/uniform-20donor-20261022-v10"
CORRECTED = ROOT / "output/h1_edge_corrected_20261024"
OUT = ROOT / "output/pancfloyd_20donor_20261026/neighborhoods"
GROUPS = ["ND", "Aab+", "T1D"]
COLORS = {"ND": "#4477AA", "Aab+": "#CC6633", "T1D": "#228833"}
TYPES = {
    "immune": ["CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune"],
    "macrophage": ["Macrophage"], "b_cell": ["B cell"],
    "t_cell": ["CD8a Tcell", "CD4 Tcell", "T cell"],
    "cd8_t_cell": ["CD8a Tcell"], "cd4_t_cell": ["CD4 Tcell"],
}
DISPLAY = {"immune": "Immune cells", "macrophage": "Macrophages", "b_cell": "B cells",
           "t_cell": "T cells", "cd8_t_cell": "CD8 T cells", "cd4_t_cell": "CD4 T cells"}
ENRICH_NAMES = {"immune": "slide06_image01_immune_peri_enrichment",
                "macrophage": "slide06_image02_macrophage_peri_enrichment",
                "b_cell": "slide06_image03_b_cell_peri_enrichment",
                "t_cell": "slide06_image04_t_cell_peri_enrichment",
                "cd8_t_cell": "slide11_image01_cd8_t_cell_peri_enrichment",
                "cd4_t_cell": "slide11_image02_cd4_t_cell_peri_enrichment"}
DIST_NAMES = {"immune": "slide06_image05_immune_signed_distance",
              "macrophage": "slide06_image06_macrophage_signed_distance",
              "b_cell": "slide06_image07_b_cell_signed_distance",
              "t_cell": "slide06_image08_t_cell_signed_distance"}
DIST_NORM = {"immune": "count", "macrophage": "proportion", "b_cell": "count", "t_cell": "proportion"}


def token(s):
    return s.replace(" ", "_").replace("+", "plus")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def symlog(x):
    x = np.asarray(x)
    return np.sign(x) * np.log10(1 + np.abs(x) / 25)


def inverse(x):
    x = np.asarray(x)
    return np.sign(x) * 25 * (10 ** np.abs(x) - 1)


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 14, "axes.labelsize": 11,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "svg.fonttype": "none", "savefig.dpi": 240})


def save(fig, name, legend):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / "figures" / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    (OUT / "figures" / f"{name}.legend.txt").write_text(legend.strip() + "\n")


def census(cases):
    rows = []
    for case in cases:
        a = pd.read_parquet(CORRECTED / "analysis_cells_by_donor" / f"{case}.parquet")
        assert not a.object_id.duplicated().any()
        eligible = a.loc[a.edge_analysis_eligible]
        row = {"case_id": case, "n_source_cells": len(a), "n_edge_eligible_cells": len(eligible),
               "n_edge_excluded_cells": len(a) - len(eligible)}
        for key, labels in TYPES.items():
            row[f"n_{key}"] = int(eligible.phenotype.isin(labels).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def make_enrichment(objects, counts):
    all_objects, all_donors, summary, tests = [], [], [], []
    for key, labels in TYPES.items():
        p = counts[f"n_{key}"].sum() / counts.n_edge_eligible_cells.sum()
        x = objects[["islet_id", "case_id", "donor_status", "object_class", "n_cells_peri"]].copy()
        x["target"] = key
        x["observed_peri_count"] = objects[[f"n_peri_{token(s)}" for s in labels]].sum(axis=1)
        x["cohort_cell_frequency"] = p
        x["expected_peri_count"] = p * x.n_cells_peri
        x["peri_enrichment_z"] = (x.observed_peri_count - x.expected_peri_count) / np.sqrt(x.expected_peri_count)
        assert np.isfinite(x.peri_enrichment_z).all()
        d = x.groupby(["case_id", "donor_status"], as_index=False).agg(
            donor_mean_peri_enrichment_z=("peri_enrichment_z", "mean"),
            n_endocrine_objects=("islet_id", "size"))
        d["target"] = key
        d["plot_x"] = d.donor_status.map({g: i for i, g in enumerate(GROUPS)}).astype(float)
        for group in GROUPS:
            ix = d.index[d.donor_status.eq(group)]
            d.loc[ix, "plot_x"] += np.linspace(-.14, .14, len(ix))
        name = ENRICH_NAMES[key]
        x.to_csv(OUT / "data" / f"{name}_objects.csv", index=False)
        d.to_csv(OUT / "data" / f"{name}_plot_data.csv", index=False)
        fig, ax = plt.subplots(figsize=(5.0, 4.2))
        for i, group in enumerate(GROUPS):
            q = d.loc[d.donor_status.eq(group)]
            y = q.donor_mean_peri_enrichment_z
            mean, sem = y.mean(), y.std(ddof=1) / np.sqrt(len(y))
            ax.bar(i, mean, width=.62, color=COLORS[group], alpha=.28,
                   edgecolor=COLORS[group], linewidth=1.2, zorder=2)
            ax.errorbar(i, mean, yerr=sem, color=COLORS[group], capsize=4, linewidth=1.4, zorder=3)
            ax.scatter(q.plot_x, y, color=COLORS[group], s=34, edgecolors="white", linewidths=.5, zorder=4)
            summary.append({"panel": name, "target": key, "donor_status": group,
                            "n_donors": len(q), "mean": mean, "sem": sem,
                            "n_endocrine_objects": int(q.n_endocrine_objects.sum()),
                            "cohort_cell_frequency": p})
        ax.axhline(0, color=".45", linewidth=.8)
        ax.set_xticks(range(3), GROUPS)
        ax.set(xlabel="Disease state", ylabel="Mean peri-enrichment z score", title=DISPLAY[key])
        ax.grid(axis="y", color=".93", linewidth=.7)
        ax.set_axisbelow(True)
        test = stats.kruskal(*[d.loc[d.donor_status.eq(g), "donor_mean_peri_enrichment_z"].to_numpy() for g in GROUPS])
        tests.append({"panel": name, "target": key, "test": "Kruskal-Wallis", "statistic_H": test.statistic,
                      "df": 2, "p_raw": test.pvalue, "multiple_comparison_adjustment": "none"})
        stars = "***" if test.pvalue < .001 else "**" if test.pvalue < .01 else "*" if test.pvalue < .05 else ""
        if stars:
            ax.text(.5, .96, stars, ha="center", va="top", transform=ax.transAxes, fontsize=17)
        fig.tight_layout()
        legend = f"""{DISPLAY[key]} in the 20 µm peri region. Each point is one donor's unweighted mean enrichment z score across that donor's eligible endocrine objects; bars and error bars show the disease-group mean ± SEM. ND n=7, Aab+ n=7, T1D n=6 donors. Eligible objects have complete cores and complete 20 µm peri regions (ND 1,504; Aab+ 1,538; T1D 1,487; 4,529 total). Both Islets and EndocrineClusters are included. The tissue-clipped peri region extends 20 µm outside the core and excludes all endocrine-object cores. Within an object each cell is counted once; peri bands can overlap between objects. The accepted edge-cell exclusions are applied. For each object, z=(observed target peri cells−expected target peri cells)/sqrt(expected target peri cells), with expected=all peri cells × the target-cell frequency among all unique edge-eligible analysed cells in the 20-donor cohort. The cohort target frequency is {p:.12g}. Target phenotypes: {', '.join(labels)}. Zero target counts are retained. No outlier removal or statistical tests are applied. App basis: scripts/uniform/assemble_uniform.py (enrich_z formula) and app/shiny_app/R/mod_spatial_server.R (donor means and mean±SEM display). The plotted mean, SEM, donor values, object values and denominators are exported separately."""
        legend = legend.replace("No outlier removal or statistical tests are applied.", f"No outlier removal is applied. The disease-group comparison uses a Kruskal–Wallis test of the 20 donor means: H(2)={test.statistic:.12g}, p={test.pvalue:.12g} (unadjusted). Asterisks at the top of the panel, if present, denote this single omnibus comparison across all three groups: *p<0.05, **p<0.01, ***p<0.001; they are not a pairwise comparison. No pairwise tests are added.")
        save(fig, name, legend)
        all_objects.append(x)
        all_donors.append(d)
    pd.concat(all_objects).to_csv(OUT / "data/enrichment_all_objects.csv", index=False)
    pd.concat(all_donors).to_csv(OUT / "data/enrichment_all_donors.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "data/enrichment_group_summary.csv", index=False)
    pd.DataFrame(tests).to_csv(OUT / "data/enrichment_kruskal_wallis.csv", index=False)


def make_age(metadata):
    name = "slide07_image02_donor_age"
    d = metadata[["case_id", "donor_status", "age_years_npod"]].copy().sort_values(["donor_status", "case_id"])
    d["plot_x"] = d.donor_status.map({g: i for i, g in enumerate(GROUPS)}).astype(float)
    fig, ax = plt.subplots(figsize=(4.3, 6.3))
    rows = []
    for i, group in enumerate(GROUPS):
        ix = d.index[d.donor_status.eq(group)]
        d.loc[ix, "plot_x"] += np.linspace(-.13, .13, len(ix))
        q = d.loc[ix]
        y = q.age_years_npod
        ax.scatter(q.plot_x, y, color=".12", s=32, zorder=3)
        ax.errorbar(i, y.mean(), yerr=y.std(ddof=1), color="black", capsize=5, linewidth=1.4, zorder=2)
        ax.hlines(y.mean(), i-.22, i+.22, color="black", linewidth=1.8)
        rows.append({"donor_status": group, "n_donors": len(q), "mean_age_years": y.mean(), "sd_age_years": y.std(ddof=1)})
    ax.set_xticks(range(3), GROUPS)
    ax.set(xlabel="Disease state", ylabel="Age (years)", title="Donor age", ylim=(0, 30))
    fig.tight_layout()
    d.to_csv(OUT / "data" / f"{name}_plot_data.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "data" / f"{name}_summary.csv", index=False)
    save(fig, name, "Each point represents one donor's age in years from the current donor_metadata.csv (age_years_npod). Horizontal bars show the arithmetic group mean and error bars show ±1 sample SD. All 20 donors are included: ND n=7, Aab+ n=7, T1D n=6. Horizontal offsets separate points only. No statistical tests are applied.")


def make_object_scatter(objects):
    a = objects[["islet_id", "case_id", "donor_status", "object_class", "n_cells_core", "n_cells_peri", "n_immune_core", "n_immune_peri"]].copy()
    a["t_cells_core_peri"] = objects[[f"n_{r}_{token(l)}" for r in ["core", "peri"] for l in TYPES["t_cell"]]].sum(axis=1)
    a["all_cells_core_peri"] = a.n_cells_core + a.n_cells_peri
    a["immune_cells_core_peri"] = a.n_immune_core + a.n_immune_peri
    a["immune_percent_core_peri"] = 100 * a.immune_cells_core_peri / a.all_cells_core_peri
    donors = sorted(a.case_id.unique())
    cmap = plt.get_cmap("tab20")
    palette = {int(case): matplotlib.colors.to_hex(cmap(i)) for i, case in enumerate(donors)}
    a["donor_color"] = a.case_id.map(palette)
    # Fixed identity-based horizontal jitter, exported exactly.
    a["plot_x"] = a.donor_status.map({g: i for i, g in enumerate(GROUPS)}) + a.islet_id.map(lambda s: (int(hashlib.sha256(s.encode()).hexdigest()[:8], 16) / (2**32-1) - .5)*.62)
    for name, value, title, ylabel in [
        ("slide10_image03_t_cell_count_by_donor", "t_cells_core_peri", "T cells", "T cells in core + peri (count)"),
        ("slide10_image04_immune_percent_by_donor", "immune_percent_core_peri", "Immune cells", "Immune cells in core + peri (%)")]:
        fig, ax = plt.subplots(figsize=(6.3, 4.8))
        for case in donors:
            q = a.loc[a.case_id.eq(case)]
            ax.scatter(q.plot_x, q[value], color=palette[int(case)], s=12, alpha=.65,
                       linewidths=0, label=str(case), rasterized=True)
        ax.set_xticks(range(3), GROUPS)
        ax.set(xlabel="Disease state", ylabel=ylabel, title=title)
        ax.set_ylim(bottom=-max(.01, a[value].max()*.025))
        if value == "immune_percent_core_peri":
            ax.set_ylim(top=103)
        ax.legend(title="Donor", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
                  fontsize=8, title_fontsize=9, markerscale=1.4, ncol=1)
        fig.tight_layout()
        a.to_csv(OUT / "data" / f"{name}_plot_data.csv", index=False)
        definition = "T cells are the mutually exclusive accepted CD8a Tcell, CD4 Tcell and T cell phenotypes. The plotted count is their summed core+peri count. This replaces the historical CD3e+ label with the actual accepted typed-cell definition; no new marker threshold is applied." if value == "t_cells_core_peri" else "The numerator is the core+peri count of CD8a Tcell, CD4 Tcell, T cell, B cell, Macrophage and residual Immune phenotypes. The denominator is every analysed core+peri cell in that object; the ratio is multiplied by 100."
        save(fig, name, f"Each point represents one eligible endocrine object, colored by donor. Both Islets and EndocrineClusters are included; ND 1,504, Aab+ 1,538 and T1D 1,487 objects from 7, 7 and 6 donors, respectively. Cores and 20 µm peri regions must both be complete; accepted edge-cell exclusions are applied. The peri region is tissue-clipped and excludes all endocrine-object cores. Within an object, core and peri cell memberships are disjoint and unique; peri memberships may recur across objects where bands overlap. {definition} Zero counts are retained. Horizontal offsets only separate points. No object-level inference or statistical tests are applied. Exact plotted values, object/donor identity and offsets are provided in the matching CSV.")
    pd.DataFrame([{"case_id": case, "color": color} for case, color in palette.items()]).to_csv(OUT / "data/donor_colors.csv", index=False)


def distance_one(case, status):
    cached = OUT / "data/distance_cells" / f"{case}.parquet"
    if cached.exists():
        return json.loads(cached.with_suffix(".json").read_text())
    a = pd.read_parquet(CORRECTED / "analysis_cells_by_donor" / f"{case}.parquet")
    a = a.loc[a.edge_analysis_eligible & a.phenotype.isin(TYPES["immune"])].copy()
    assert not a.object_id.duplicated().any()
    cores, _, _, _, _ = geometry(SOURCE, case)
    names = [n for n, (cl, _) in cores.items() if cl == "Islet"]
    gs = np.asarray([cores[n][1] for n in names], dtype=object)
    boundaries = shapely.boundary(gs)
    points = shapely.points(a.x_um.to_numpy(), a.y_um.to_numpy())
    nearest = STRtree(boundaries).nearest(points)
    dist = shapely.distance(points, boundaries[nearest])
    # The app defines a negative sign inside any Islet, positive otherwise.
    inside = np.zeros(len(points), dtype=bool)
    intersections = STRtree(gs).query(points, predicate="intersects")
    inside[intersections[0]] = True
    dist[inside] *= -1
    flags = pd.read_parquet(CORRECTED / "eligibility/objects" / f"{case}.parquet").set_index("object_name")
    anchor_names = np.asarray(names, dtype=object)[nearest]
    complete = flags.loc[anchor_names, "core_complete"].to_numpy(bool)
    a["donor_status"] = status
    a["nearest_islet"] = anchor_names
    a["nearest_islet_complete"] = complete
    a["signed_distance_um"] = dist
    a[["object_id", "case_id", "phenotype", "donor_status", "x_um", "y_um", "nearest_islet", "nearest_islet_complete", "signed_distance_um"]].to_parquet(cached, index=False)
    audit = {"case_id": int(case), "donor_status": status, "n_edge_eligible_immune_cells": len(a),
             "n_complete_nearest_islet": int(complete.sum()),
             "n_incomplete_nearest_islet": int((~complete).sum()), "n_islet_anchors": len(names)}
    cached.with_suffix(".json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def make_distances(counts):
    a = pd.concat([pd.read_parquet(p) for p in sorted((OUT / "data/distance_cells").glob("*.parquet"))], ignore_index=True)
    assert not a.duplicated(["case_id", "object_id"]).any()
    a = a.loc[a.nearest_islet_complete].copy()
    denoms = counts.groupby("donor_status").n_edge_eligible_cells.sum()
    totals = []
    for key, name in DIST_NAMES.items():
        q = a.loc[a.phenotype.isin(TYPES[key])]
        lo = min(float(q.signed_distance_um.quantile(.001)), -50.)
        hi = 200.
        tedges = np.linspace(symlog(lo), symlog(hi), 161)
        edges = inverse(tedges)
        widths = np.diff(edges)
        centers = (tedges[:-1] + tedges[1:]) / 2
        rows = []
        fig, ax = plt.subplots(figsize=(6.1, 4.4))
        for group in GROUPS:
            cells = q.loc[q.donor_status.eq(group)]
            dists = cells.signed_distance_um.to_numpy()
            keep = (dists >= lo) & (dists <= hi)
            # Reproduce the app's clamp of an observation exactly at the upper edge.
            bins = np.clip(np.floor((symlog(dists[keep])-tedges[0])/(tedges[-1]-tedges[0])*160).astype(int), 0, 159)
            n = np.bincount(bins, minlength=160)
            y = n / widths
            if DIST_NORM[key] == "proportion":
                y = y / denoms[group]
            ax.fill_between(centers, y, color=COLORS[group], alpha=.13)
            ax.plot(centers, y, color=COLORS[group], linewidth=1.35, label=group)
            for b in range(160):
                rows.append({"donor_status": group, "target": key, "bin_index": b+1,
                             "lower_distance_um": edges[b], "upper_distance_um": edges[b+1],
                             "width_um": widths[b], "center_symlog": centers[b],
                             "center_distance_um": inverse(centers[b]), "target_cell_count": int(n[b]),
                             "all_edge_eligible_group_cells": int(denoms[group]),
                             "normalization": DIST_NORM[key], "plotted_value_per_um": y[b]})
            totals.append({"panel": name, "target": key, "donor_status": group,
                           "n_target_cells_complete_anchor": len(cells), "n_target_cells_plotted": int(keep.sum()),
                           "n_below_plot_range": int((dists < lo).sum()), "n_above_plot_range": int((dists > hi).sum()),
                           "all_edge_eligible_group_cells": int(denoms[group]), "lower_distance_um": lo,
                           "upper_distance_um": hi, "normalization": DIST_NORM[key]})
        ticks = [v for v in [-200, -100, -50, -20, 0, 20, 50, 100, 200] if lo <= v <= hi]
        ax.set_xticks(symlog(ticks), [str(v) for v in ticks])
        ax.axvline(0, color=".35", linewidth=.9, linestyle="--")
        ax.axvline(symlog(20), color=".55", linewidth=.8, linestyle=":")
        ax.set_xlim(symlog(lo), symlog(hi))
        ylabel = "Cells per µm" if DIST_NORM[key] == "count" else "Fraction of all cells per µm"
        ax.set(xlabel="Signed distance to islet boundary (µm)", ylabel=ylabel, title=DISPLAY[key], ylim=(0, None))
        if DIST_NORM[key] == "proportion":
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-3, 3))
            ax.yaxis.set_major_formatter(formatter)
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()
        pd.DataFrame(rows).to_csv(OUT / "data" / f"{name}_plot_data.csv", index=False)
        norm = "The ordinate is the pooled number of target cells in each bin divided by that bin's width in real µm." if DIST_NORM[key] == "count" else "The ordinate is the target-cell count in each bin divided by that bin's width in real µm and by the number of all unique edge-eligible analysed cells in that disease group. The denominator is all cell phenotypes, not just the displayed target phenotype, and is not restricted to the plotted distance window."
        save(fig, name, f"Signed distance of {DISPLAY[key].lower()} from the nearest Islet boundary, using the current app's Islet-only distance option. Negative distances are inside an Islet, zero marks its boundary, and positive distances are outside every Islet; EndocrineClusters do not serve as distance anchors. Target phenotypes: {', '.join(TYPES[key])}. Cell centroids and v10 core outlines are in µm (0.3774 µm/pixel). Every cell has a unique donor+object_id identity; accepted edge-excluded cells are removed. A cell whose nearest Islet has an incomplete core is excluded rather than reassigned to a different Islet. Cells are pooled within ND (7 donors), Aab+ (7) or T1D (6). {norm} The curves connect histogram values, with no kernel smoothing. The 160 bins are uniform in sign(d)×log10(1+|d|/25), matching the app; the x-axis uses this transformation and tick labels are real µm. The lower bound is min(the pooled target distance 0.1th percentile, −50 µm), here {lo:.12g} µm; the upper bound is the current-app default 200 µm. Dashed and dotted vertical lines mark 0 and 20 µm. Curves have no error bands, and no statistical tests are applied. The matching CSV contains every bin count, real width, denominator and ordinate; cell-level distances and donor exclusion counts are also retained. App basis: app/shiny_app/R/spatial_helpers.R (.query_cell_distance_hist_memory) and scripts/uniform/assemble_uniform.py (signed_dist_islet).")
    pd.DataFrame(totals).to_csv(OUT / "data/signed_distance_group_totals.csv", index=False)


def panel_manifest():
    rows = []
    for key, stem in ENRICH_NAMES.items():
        rows.append({"image_id": "_".join(stem.split("_")[:2]), "figure_stem": stem,
                     "source_csvs": [f"data/{stem}_plot_data.csv", f"data/{stem}_objects.csv", "data/enrichment_group_summary.csv", "data/enrichment_kruskal_wallis.csv"],
                     "x": "donor_status (plotted at plot_x)", "y": "donor_mean_peri_enrichment_z", "units": "dimensionless z score",
                     "point_unit": "donor", "aggregation": "unweighted mean of eligible object z scores within donor; group mean ± SEM across donors",
                     "n_donors": 20, "n_objects": 4529, "target": key, "test": "Kruskal-Wallis across 3 disease groups; unadjusted"})
    for key, stem in DIST_NAMES.items():
        rows.append({"image_id": "_".join(stem.split("_")[:2]), "figure_stem": stem,
                     "source_csvs": [f"data/{stem}_plot_data.csv", "data/signed_distance_group_totals.csv", "data/distance_exclusion_audit.csv", "data/unique_cell_census.csv"],
                     "x": "center_symlog (tick labels show center_distance_um)", "y": "plotted_value_per_um",
                     "units": "cells/µm" if DIST_NORM[key] == "count" else "fraction of all disease-group cells/µm",
                     "point_unit": "histogram bin", "aggregation": "unique cells pooled by disease; count divided by real bin width" + (" and all disease-group cells" if DIST_NORM[key] == "proportion" else ""),
                     "n_donors": 20, "n_objects": None, "target": key, "test": "none"})
    for stem, y, units, point, agg, target in [
        ("slide07_image02_donor_age", "age_years_npod", "years", "donor", "group mean ± sample SD", "age"),
        ("slide10_image03_t_cell_count_by_donor", "t_cells_core_peri", "cells", "endocrine object", "sum of mutually exclusive typed T-cell counts in core+peri", "t_cell"),
        ("slide10_image04_immune_percent_by_donor", "immune_percent_core_peri", "%", "endocrine object", "100 × broad immune core+peri cells / all core+peri cells", "immune")]:
        rows.append({"image_id": "_".join(stem.split("_")[:2]), "figure_stem": stem,
                     "source_csvs": [f"data/{stem}_plot_data.csv"], "x": "donor_status (plotted at plot_x)",
                     "y": y, "units": units, "point_unit": point, "aggregation": agg,
                     "n_donors": 20, "n_objects": 4529 if point == "endocrine object" else None,
                     "target": target, "test": "none"})
    (OUT/"panel_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
    flat = [{**r, "source_csvs": ";".join(r["source_csvs"])} for r in rows]
    pd.DataFrame(flat).to_csv(OUT/"panel_manifest.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    for name in ["figures", "data", "data/distance_cells"]:
        (OUT/name).mkdir(parents=True, exist_ok=True)
    style()
    features = pd.read_parquet(CORRECTED / "h1_object_features.parquet")
    cases = sorted(int(c) for c in features.case_id.unique())
    assert len(cases) == 20 and 78 not in cases
    groupmap = features[["case_id", "donor_status"]].drop_duplicates()
    assert len(groupmap) == 20
    meta = pd.read_csv(SOURCE / "donor_metadata.csv").merge(groupmap, on="case_id", validate="one_to_one")
    objects = features.loc[features.core_complete & features.neighborhood_complete_20].copy()
    assert len(objects) == 4529
    assert objects[["n_cells_core", "n_cells_peri"]].gt(0).all().all()
    objects[["islet_id", "case_id", "donor_status", "object_class", "core_complete", "neighborhood_complete_20"]].to_csv(OUT / "data/object_eligibility.csv", index=False)
    counts = census(cases).merge(groupmap, on="case_id", validate="one_to_one")
    counts.to_csv(OUT / "data/unique_cell_census.csv", index=False)
    print("Census complete:", len(objects), "objects;", int(counts.n_edge_eligible_cells.sum()), "edge-eligible unique cells", flush=True)
    make_enrichment(objects, counts)
    make_age(meta)
    make_object_scatter(objects)
    print("Nine enrichment/age/object panels written", flush=True)
    audits = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(distance_one, int(row.case_id), row.donor_status): row.case_id for row in groupmap.itertuples()}
        for future in as_completed(futures):
            result = future.result()
            audits.append(result)
            print("Distance complete", result["case_id"], result["n_edge_eligible_immune_cells"], flush=True)
    pd.DataFrame(audits).sort_values("case_id").to_csv(OUT / "data/distance_exclusion_audit.csv", index=False)
    make_distances(counts)
    panel_manifest()
    sources = [SOURCE/"donor_metadata.csv", CORRECTED/"h1_object_features.parquet", Path(__file__).resolve(),
               ROOT/"scripts/uniform/assemble_uniform.py", ROOT/"app/shiny_app/R/spatial_helpers.R",
               ROOT/"app/shiny_app/R/mod_spatial_server.R"]
    for case in cases:
        sources.extend([CORRECTED/"analysis_cells_by_donor"/f"{case}.parquet",
                        CORRECTED/"eligibility/objects"/f"{case}.parquet", SOURCE/"gson"/f"{case:04d}.geojson.gz"])
    manifest = {"script": str(Path(__file__).resolve()), "data_root": str(SOURCE), "corrected_root": str(CORRECTED),
                "donors": cases, "n_objects": len(objects), "object_counts_by_disease": objects.groupby("donor_status").size().to_dict(),
                "source_sha256": {str(p): sha(p) for p in sources},
                "figures": sorted(p.name for p in (OUT/"figures").glob("*.png")),
                "distance_anchor": "nearest Islet only; nearest incomplete Islet excluded, not reassigned",
                "tests": "donor-level Kruskal-Wallis for six enrichment panels; all other panels descriptive", "quantitative_INS": False}
    (OUT/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (OUT/"README.md").write_text("""# Current 20-donor neighborhood panels

Thirteen panels recreate the assigned scientific plot types in PancFloyd_Figures_v3.pptx using the accepted v10 typed cells, geometry and corrected edge/count tables. Each has PNG, PDF, SVG, a separate methods legend, and exact plotted CSV data. The six enrichment panels retain the donor-level Kruskal–Wallis comparison, with current values in the legends/table and explicitly defined omnibus asterisks. No marker-intensity analyses or biological result interpretation are added. `panel_manifest.json` and `.csv` map every original image ID to plotted sources, x/y fields, units and aggregation.

- Slide 6, images 1–4: donor mean peri-enrichment z scores, group mean ± SEM, for the broad immune union, macrophages, B cells and broad T-cell union.
- Slide 6, images 5–8: current-app Islet-only signed-distance histogram profiles. Immune and B-cell ordinates are pooled count/µm; macrophage and T-cell ordinates are fraction of all disease-group cells/µm, retaining the reference normalization choices. Current-app settings are 160 transformed bins, C=25 µm and 200 µm maximum outside distance.
- Slide 11, images 1–2: CD8 and CD4 T-cell peri-enrichment, as above.
- Slide 7, image 2: donor ages, group mean ± SD.
- Slide 10, images 3–4: object core+peri typed T-cell counts and broad immune percentages, colored by donor.

Object panels use complete cores and complete 20 µm neighborhoods: 4,529 objects (ND 1,504; Aab+ 1,538; T1D 1,487), not the 4,179 subset that additionally requires a complete 100 µm field for other analyses. All 20 donors contribute (7 ND, 7 Aab+, 6 T1D). Peri bands remain object-associated and may overlap between objects. The cohort background and profile denominators count unique edge-eligible cells, never repeated peri memberships. Complete distance anchors and cell exclusions are recorded separately. The historical residual Immune/T cell labels are replaced by the explicit accepted broad unions; the historical CD3e+ panel is labeled as typed T-cell count. No source data or labels are overwritten.

See each `.legend.txt` for its complete definition, and `manifest.json` for source hashes. `data/distance_cells/` retains cell-level identity, anchor and signed distance. These pooled profiles are descriptive and do not treat individual cells as independent donor replicates.
""")
    print("Done:", len(manifest["figures"]), "figures", flush=True)


if __name__ == "__main__":
    main()
