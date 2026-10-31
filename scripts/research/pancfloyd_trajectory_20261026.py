#!/usr/bin/env python3
"""Recreate the requested reference-deck trajectory panels from current app data.

Uses saved UMAP/Leiden without fitting them. Calls the existing app DPT module
once, then the existing R plot helper for descriptive weighted LOESS/jitter.
No hypothesis tests or biological interpretations are added.
"""
from __future__ import annotations
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/pancfloyd_numba")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pancfloyd_mpl")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import hashlib
import importlib.util
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data/development/uniform-20donor-20261022-v10"
OUT = REPO / "output/pancfloyd_20donor_20261026/trajectory"
for p in (OUT / "figures", OUT / "data"):
    p.mkdir(parents=True, exist_ok=True)
STATUS = ["ND", "Aab+", "T1D"]
COLORS = {"ND": "#4477AA", "Aab+": "#CC6633", "T1D": "#228833"}
MARKERS = ["CD3e", "CD8a", "CD68", "CD4", "CD20"]
HEATMARKERS = ["CD68", "CD20", "CD3e", "CD4", "CD8a"]
ENDO = ["Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative", "Endocrine"]
PHENO = ["Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative", "CD8a Tcell", "CD4 Tcell", "B cell", "Macrophage", "Exocrine"]
PHENO_COLORS = ["#D62728", "#FF9900", "#FFD92F", "#9467BD", "#8C2D99", "#386CB0", "#E78AC3", "#5E4FA2", "#999999"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                     "axes.titlesize": 15, "axes.labelsize": 12,
                     "xtick.labelsize": 11, "ytick.labelsize": 11,
                     "legend.fontsize": 11, "svg.fonttype": "none",
                     "pdf.fonttype": 42, "savefig.facecolor": "white"})


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / "figures" / f"{name}.{ext}", dpi=260, bbox_inches="tight")
    plt.close(fig)


def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#B7B7B7")
    ax.spines["bottom"].set_color("#B7B7B7")
    ax.grid(color="#EEEEEE", linewidth=.7)
    ax.set_axisbelow(True)


def main():
    source = ROOT / "adata_ins_root.h5ad"
    a = ad.read_h5ad(source)
    obs = a.obs.copy().reset_index(drop=True)
    obs["islet_id"] = obs["combined_islet_id"].astype(str)
    obs["case_id"] = obs["imageid"].astype(str)
    obs["donor_status"] = obs["donor_status"].astype(str)
    obs["umap_1"] = np.asarray(a.obsm["X_umap"])[:, 0]
    obs["umap_2"] = np.asarray(a.obsm["X_umap"])[:, 1]
    assert obs.case_id.nunique() == 20 and "78" not in set(obs.case_id)
    assert not obs.islet_id.duplicated().any()
    for m in MARKERS:
        obs[m] = np.asarray(a.X[:, a.var_names.get_loc(m)]).reshape(-1)
    obs["n_core_endocrine"] = sum(obs[f"count_{t}"] for t in ENDO)
    obs["beta_endocrine_fraction"] = obs["count_Beta cell"] / obs.n_core_endocrine.replace(0, np.nan)
    edge_path = REPO / "output/h1_edge_corrected_20261024/h1_object_features.parquet"
    edge = pd.read_parquet(edge_path, columns=["islet_id", "core_complete", "neighborhood_complete_20", "edge_exclusion_reason"])
    obs = obs.merge(edge, on="islet_id", how="left", validate="one_to_one")
    assert obs.core_complete.notna().all() and obs.neighborhood_complete_20.notna().all()
    obs["eligible_plotted"] = obs.in_composition_embedding.astype(bool) & obs.core_complete & obs.neighborhood_complete_20
    fitted = obs.loc[obs.in_composition_embedding.astype(bool)].copy()
    assert len(fitted) == 1215 and fitted["leiden_0.5"].nunique() == 10

    # Bind the existing app module to this exact root, bypassing older fallback paths.
    os.environ["ISLET_DATA_ROOT"] = str(ROOT)
    module_path = REPO / "app/shiny_app/python/paga_pseudotime.py"
    spec = importlib.util.spec_from_file_location("app_paga_pseudotime", module_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod._FEATURES_CSV = str(ROOT / "barlow_islet_features.csv")
    bindings = {"h5ad_sha256": digest(source), "barlow_sha256": digest(ROOT / "barlow_islet_features.csv"),
                "dpt_module_sha256": digest(module_path), "min_cells": 10}
    cache = OUT / "data/app_dpt_snapshot.csv"
    cache_meta = OUT / "data/app_dpt_metadata.json"
    if cache.exists() and cache_meta.exists() and json.loads(cache_meta.read_text()).get("source_bindings") == bindings:
        dpt = pd.read_csv(cache, dtype={"islet_id": str})
        dpt_meta = json.loads(cache_meta.read_text())
    else:
        result = mod.compute_pseudotime(min_cells=10, max_cells=None, object_classes=None)
        if "error" in result:
            raise RuntimeError(result["error"])
        dpt = pd.DataFrame(result)
        ctx = mod._context(10, None, None)
        dpt["pt_rank"] = (dpt.pseudotime.rank(method="average") - .5) / dpt.pseudotime.notna().sum()
        root = int(ctx["root"])
        dpt_meta = {"source_bindings": bindings, "data_root": str(ROOT), "objects": len(dpt),
                    "profiles": ctx["A"].n_obs, "features": ctx["features"], "root_profile": root,
                    "root_profile_representative_donor": str(ctx["A"].obs.iloc[root].case_id),
                    "orientation_flipped": bool(ctx.get("flipped", False)),
                    "n_neighbors": mod.N_NEIGHBORS, "metric": mod.METRIC,
                    "n_diffusion_components": mod.N_DIFFMAP_COMPS,
                    "software": {p: version(p) for p in ["scanpy", "anndata", "numpy", "scipy", "scikit-learn", "numba"]}}
        dpt.to_csv(cache, index=False)
        cache_meta.write_text(json.dumps(dpt_meta, indent=2) + "\n")
    assert len(dpt) == 1215 and np.isfinite(dpt.pseudotime).all()
    fitted = fitted.merge(dpt, on="islet_id", how="left", validate="one_to_one")
    assert fitted.pt_rank.notna().all()
    full_pt_range = (float(fitted.pt_rank.min()), float(fitted.pt_rank.max()))
    fields = ["islet_id", "case_id", "donor_status", "object_class", "total_cells_core", "total_cells_peri",
              "umap_1", "umap_2", "leiden_0.5", "pseudotime", "pt_rank", "n_core_endocrine",
              "count_Beta cell", "beta_endocrine_fraction", "in_composition_embedding", "core_complete",
              "neighborhood_complete_20", "eligible_plotted", "edge_exclusion_reason"] + MARKERS
    fitted[fields].to_csv(OUT / "data/object_representation.csv", index=False)
    obs[["islet_id", "case_id", "donor_status", "object_class", "in_composition_embedding", "core_complete",
         "neighborhood_complete_20", "eligible_plotted", "edge_exclusion_reason"]].to_csv(OUT / "data/all_paired_object_eligibility.csv", index=False)
    counts = obs.groupby(["case_id", "donor_status"], observed=True).agg(paired_objects=("islet_id", "size"),
               embedded_objects=("in_composition_embedding", "sum"), plotted_objects=("eligible_plotted", "sum")).reset_index()
    counts.to_csv(OUT / "data/donor_object_counts.csv", index=False)
    fitted = fitted[fitted.eligible_plotted].copy()
    assert len(fitted) == 1126 and fitted.case_id.nunique() == 20
    fitted[fields].to_csv(OUT / "data/plot_object_representation.csv", index=False)

    # Reference slide 5, image 1: the same saved coordinates, colored in two ways.
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    for s in STATUS:
        z = fitted[fitted.donor_status == s]
        axes[0].scatter(z.umap_1, z.umap_2, s=11, alpha=.72, linewidths=0, label=s, color=COLORS[s])
    axes[0].set_title("Endocrine composition by disease")
    axes[0].legend(frameon=False, title="Disease", loc="upper right", markerscale=2)
    labels = sorted(fitted["leiden_0.5"].astype(str).unique(), key=int)
    cluster_colors = dict(zip(labels, plt.get_cmap("tab10").colors))
    for cl in labels:
        z = fitted[fitted["leiden_0.5"].astype(str) == cl]
        axes[1].scatter(z.umap_1, z.umap_2, s=11, alpha=.8, linewidths=0, label=cl, color=cluster_colors[cl])
    axes[1].set_title("Endocrine composition clusters")
    axes[1].legend(frameon=False, title="Leiden 0.5", loc="upper right", ncol=2, markerscale=2)
    for ax in axes:
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.margins(x=.16, y=.1)
    axes[0].legend(frameon=False, title="Disease", loc="upper left", bbox_to_anchor=(.96, 1), markerscale=2)
    axes[1].legend(frameon=False, title="Leiden 0.5", loc="upper left", bbox_to_anchor=(.96, 1), ncol=2, markerscale=2)
    save(fig, "slide05_image01_umap_disease_and_leiden")

    # Exact current cluster-composition definition: average object-level core prop_*.
    comp_rows = []
    comp_objects = fitted[["islet_id", "case_id", "donor_status", "leiden_0.5", "total_cells_core"] + ["prop_"+p for p in PHENO]].copy()
    comp_objects.to_csv(OUT / "data/cluster_composition_objects.csv", index=False)
    fig, ax = plt.subplots(figsize=(10.8, 5.7), constrained_layout=True)
    bottom = np.zeros(len(labels))
    for phenotype, color in zip(PHENO, PHENO_COLORS):
        means = []
        for cl in labels:
            z = fitted[fitted["leiden_0.5"].astype(str) == cl]
            val = float(z["prop_"+phenotype].mean())
            means.append(val)
            comp_rows.append({"cluster": cl, "phenotype": phenotype, "mean_core_fraction": val,
                              "n_objects": len(z), "n_donors": z.case_id.nunique()})
        if max(means) > 0:
            ax.bar(labels, means, bottom=bottom, label=phenotype, color=color, width=.8)
        bottom += np.array(means)
    pd.DataFrame(comp_rows).to_csv(OUT / "data/cluster_composition_means.csv", index=False)
    assert np.all((bottom <= 1 + 1e-10) & (bottom > 0))
    ax.set(xlabel="Leiden cluster (resolution 0.5)", ylabel="Mean fraction of core cells", ylim=(0, 1.02),
           title="Mean phenotype composition by cluster")
    tidy(ax); ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, "slide05_image02_cluster_composition")

    # Marker values come directly from current app .X. INS is not extracted.
    point_rows = []
    for feature in ["Beta endocrine fraction"] + MARKERS:
        z = fitted[["islet_id", "case_id", "donor_status", "total_cells_core", "pt_rank", "pseudotime"]].copy()
        z["feature"] = feature
        z["value"] = fitted.beta_endocrine_fraction if feature == "Beta endocrine fraction" else fitted[feature]
        z["units"] = "Beta cells / typed core endocrine cells" if feature == "Beta endocrine fraction" else "cohort-scaled marker mean across core cells (0–1)"
        point_rows.append(z)
    pd.concat(point_rows, ignore_index=True).to_csv(OUT / "data/trajectory_points.csv", index=False)
    subprocess.run(["Rscript", str(REPO / "scripts/research/pancfloyd_trajectory_20261026.R"), str(REPO), str(OUT)], check=True)
    points = pd.read_csv(OUT / "data/trajectory_plotted_points.csv")
    curves = pd.read_csv(OUT / "data/trajectory_loess_curves.csv")
    mappings = [("Beta endocrine fraction", "slide05_image03_beta_fraction"), ("CD3e", "slide05_image04_cd3e"),
                ("CD8a", "slide05_image05_cd8a"), ("CD68", "slide11_image03_cd68"),
                ("CD4", "slide11_image04_cd4"), ("CD20", "slide11_image05_cd20")]
    for feature, filename in mappings:
        fig, ax = plt.subplots(figsize=(10.8, 6.1), constrained_layout=True)
        beta = feature == "Beta endocrine fraction"
        z = points[points.feature == feature]
        order = np.random.default_rng(42).permutation(len(z))
        zr = z.iloc[order]
        ax.scatter(zr.x_plot, zr.y_plot, s=11, c=zr.donor_status.map(COLORS), alpha=.62, linewidths=0)
        for status in ([] if beta else STATUS):
            cr = curves[(curves.feature == feature) & (curves.donor_status == status)]
            ax.fill_between(cr.x, cr.lower, cr.upper, color=COLORS[status], alpha=.12, linewidth=0)
            ax.plot(cr.x, cr.y, color=COLORS[status], linewidth=2)
        ax.set_xlabel("Pseudotime (rank-normalized)")
        ax.set_ylabel("Beta cells / core endocrine cells" if beta else f"{feature} — scaled mean across core cells (0–1)")
        ax.set_title("Beta-cell fraction along pseudotime" if beta else f"{feature} along pseudotime")
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.set_xlim(-.012, 1.012)
        if beta:
            ax.set_ylim(-.015, 1.015)
            ax.set_yticks([0, .25, .5, .75, 1])
        tidy(ax)
        handles = [Line2D([], [], marker="o", linestyle="None" if beta else "-", color=COLORS[s], label=s, markersize=5) for s in STATUS]
        ax.legend(handles=handles, title="Disease", frameon=False, loc="best")
        save(fig, filename)
        z.to_csv(OUT / "data" / f"{filename}_points.csv", index=False)
        if not beta:
            curves[curves.feature == feature].to_csv(OUT / "data" / f"{filename}_curves.csv", index=False)

    # App heatmap algorithm, including sample SD across 20 bin means and clipping.
    lo, hi = full_pt_range
    normalized = (fitted.pt_rank - lo) / (hi - lo)
    bins = np.minimum((normalized * 20).astype(int), 19)
    heat_rows = []
    for marker in HEATMARKERS:
        means = np.array([fitted.loc[bins == b, marker].mean() for b in range(20)])
        sd = np.nanstd(means, ddof=1)
        if not np.isfinite(sd) or sd == 0:
            sd = 1
        zscore = (means - np.nanmean(means)) / sd
        for b in range(20):
            heat_rows.append({"marker": marker, "bin": b+1, "left": lo+(hi-lo)*b/20,
                              "right": lo+(hi-lo)*(b+1)/20, "midpoint": lo+(hi-lo)*(b+.5)/20,
                              "bin_mean": means[b], "z_unclipped": zscore[b], "z_plot": np.clip(zscore[b], -2.5, 2.5),
                              "n_objects": int((bins == b).sum()),
                              "n_donors": fitted.loc[bins == b, "case_id"].nunique()})
    heat = pd.DataFrame(heat_rows)
    heat.to_csv(OUT / "data/immune_heatmap_values.csv", index=False)
    membership = fitted[["islet_id", "case_id", "donor_status", "pt_rank"]].copy()
    membership["bin"] = bins + 1
    for m in HEATMARKERS:
        membership[m] = fitted[m]
    membership.to_csv(OUT / "data/immune_heatmap_object_membership.csv", index=False)
    fig, ax = plt.subplots(figsize=(11, 4.2), constrained_layout=True)
    mat = heat.pivot(index="marker", columns="bin", values="z_plot").loc[HEATMARKERS].to_numpy()
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto", interpolation="nearest",
                   extent=[lo, hi, len(HEATMARKERS)-.5, -.5])
    ax.set_yticks(range(len(HEATMARKERS)), HEATMARKERS)
    ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xlim(0, 1)
    ax.set_xlabel("Pseudotime (rank-normalized)")
    ax.set_title("Immune marker profiles along pseudotime")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=.03, pad=.025); cb.set_label("Z-score across bin means")
    save(fig, "slide05_image07_immune_heatmap")

    by_status = fitted.groupby("donor_status").agg(objects=("islet_id", "size"), donors=("case_id", "nunique")).to_dict("index")
    mapping_rows = [
        (5, 1, "UMAP by disease and Leiden 0.5", "slide05_image01_umap_disease_and_leiden", "plot_object_representation.csv", "Saved current composition coordinates; ten current clusters"),
        (5, 2, "Cluster composition", "slide05_image02_cluster_composition", "cluster_composition_means.csv;cluster_composition_objects.csv", "Current retained phenotypes and core-cell denominator"),
        (5, 3, "INS trajectory", "slide05_image03_beta_fraction", "slide05_image03_beta_fraction_points.csv", "Replaced quantitative INS with typed Beta fraction of core endocrine cells; no smoothing"),
        (5, 4, "CD3e trajectory", "slide05_image04_cd3e", "slide05_image04_cd3e_points.csv;slide05_image04_cd3e_curves.csv", "Current served marker values"),
        (5, 5, "CD8a trajectory", "slide05_image05_cd8a", "slide05_image05_cd8a_points.csv;slide05_image05_cd8a_curves.csv", "Current served marker values"),
        (5, 7, "Immune heatmap", "slide05_image07_immune_heatmap", "immune_heatmap_values.csv;immune_heatmap_object_membership.csv", "Current app 20-bin heatmap algorithm; reference row order"),
        (11, 3, "CD68 trajectory", "slide11_image03_cd68", "slide11_image03_cd68_points.csv;slide11_image03_cd68_curves.csv", "Current served marker values"),
        (11, 4, "CD4 trajectory", "slide11_image04_cd4", "slide11_image04_cd4_points.csv;slide11_image04_cd4_curves.csv", "Current served marker values"),
        (11, 5, "CD20 trajectory", "slide11_image05_cd20", "slide11_image05_cd20_points.csv;slide11_image05_cd20_curves.csv", "Current served marker values")]
    panelmap = pd.DataFrame(mapping_rows, columns=["reference_slide", "reference_image", "reference_panel", "output_stem", "data_csv", "change_note"])
    panelmap["x_field"] = ["umap_1", "cluster", "x_plot", "x_plot", "x_plot", "midpoint", "x_plot", "x_plot", "x_plot"]
    panelmap["y_field"] = ["umap_2", "mean_core_fraction", "y_plot", "y_plot", "y_plot", "marker", "y_plot", "y_plot", "y_plot"]
    panelmap["units"] = ["UMAP coordinates", "mean phenotype / all core cells", "Beta / core endocrine cells", "served CD3e core mean 0–1", "served CD8a core mean 0–1", "color=z_plot; bin-mean z-score clipped ±2.5", "served CD68 core mean 0–1", "served CD4 core mean 0–1", "served CD20 core mean 0–1"]
    panelmap["n_plotted_objects"] = len(fitted)
    panelmap.to_csv(OUT / "panel_mapping.csv", index=False)
    legends = f"""# Trajectory panel legends and methods

These panels recreate the specified plots from `PancFloyd_Figures_v3.pptx` using the current app root `{ROOT.relative_to(REPO)}`. They show descriptive measurements and the existing app curves; no disease hypothesis tests, biological interpretations or if–then branch judgments are assigned here.

**Common cohort.** The saved composition embedding and the exact app DPT calculation both include 1,215 endocrine objects with at least ten total core cells. All panels display the same 1,126 objects with corrected complete core and complete 20-µm neighborhood, from all 20 donors. Neither coordinates, cluster labels nor pseudotime ranks are refitted after this display filter. Plotted counts by disease are {json.dumps(by_status)}. The 3,548 paired objects outside this representation and 89 embedded objects failing the stated corrected edge eligibility are not plotted. This common display filter does not impose hypothesis-specific immune-exposure or structure-availability requirements. `donor_object_counts.csv` records all donor denominators, `all_paired_object_eligibility.csv` records all 4,763 paired objects, and `object_representation.csv` retains the full 1,215-object fitted representation with `eligible_plotted`. The exact 1,126 displayed rows are in `plot_object_representation.csv`. Individual islets and Clusters are called endocrine objects, not tissue sections.

**Slide 5 image 1, composition UMAP.** Both panels use the same saved `X_umap` coordinates without refitting. The left colors disease; the right colors saved `leiden_0.5`, with ten clusters. Four standardized core-endocrine label shares define the saved representation: Beta, Alpha, Delta and CHGA+ hormone-negative. Its graph uses 30 Euclidean neighbors, seed 0, no donor correction, UMAP min_dist 3/spread 3. The exact plotted coordinates, donor IDs and cluster labels are in `plot_object_representation.csv`. Different current coordinates or cluster counts from the historical reference reflect use of the current app representation.

**Slide 5 image 2, cluster composition.** Each stacked segment is the arithmetic mean, across the objects assigned to that saved cluster, of the current `prop_<phenotype>` field: number of typed core cells of that phenotype divided by all analysed core cells. Objects, rather than cells pooled across the cluster, receive equal weight. Constant-zero and retired phenotype categories are omitted as in the current app; omitted retired categories are not renormalized away, so displayed bars can sum to less than one. The endocrine denominator used to fit the composition embedding differs from this all-core-cell display denominator. `cluster_composition_means.csv` contains plotted means and cluster sample sizes; `cluster_composition_objects.csv` supplies the contributing object values. No structural annotation is represented as a cell phenotype.

**Pseudotime calculation.** The script imports and calls the current `app/shiny_app/python/paga_pseudotime.py` with `min_cells=10`, no upper limit and no class filter, explicitly bound to the current root's `barlow_islet_features.csv`. Stored features are log1p(core or peri phenotype count / max(core endocrine count,1)); retired and constant features are excluded. The 16 active features produce 1,213 distinct rounded profiles. The existing implementation uses 15 cosine neighbors and ten diffusion components, selects an ND-containing root profile nearest the ND diffusion centroid, and applies its existing Beta-ratio orientation. The common displayed rank is `(average rank − 0.5) / 1215`, calculated before marker filtering. Root, orientation, versions and exact input hashes are in `app_dpt_metadata.json`; raw and ranked DPT values are in `app_dpt_snapshot.csv`. This is a cross-sectional profile ordering, not elapsed time.

**Slide 5 image 3, Beta-cell fraction.** This explicitly replaces the reference's quantitative INS panel in accordance with the owner's acquisition instruction. The measured y value is accepted typed Beta-cell count divided by typed core-endocrine count (Beta + Alpha + Delta + CHGA+ hormone-negative + generic Endocrine). No INS mean, sum or intensity transformation is used. Beta composition also contributes to the representation and its orientation; this panel describes that defined state feature. This replacement panel shows the measured fraction directly, with horizontal jitter only and no smooth curve or confidence band. The axis is the 0–1 fraction, not quantitative INS.

**Slide 5 images 4–5 and slide 11 images 3–5.** CD3e, CD8a, CD68, CD4 and CD20 values are read directly from the corresponding current app H5AD `.X` marker column. They are arithmetic means across analysed core cells of the existing cohort-scaled cell measurements, in the served 0–1 scale; they are not marker-positive cell counts, calibrated protein concentrations or within-immune-cell means. The original reference panels' quantitative INS values are not reused.

**Scatter display and curves.** All six panels use the same full-cohort rank axis and disease colors. The five non-INS marker panels use current app descriptive smoothing. R sources the existing app `plot_loess_curve` helper: degree-two LOESS, span 0.75, weight `log1p(total core cells)`, with 80 predictions per disease over its observed range. Shading is the app's pointwise 95% LOESS band (t multiplier with the fitted LOESS degrees of freedom); it is not a donor-replicated disease-comparison confidence interval. All displayed points use seed-42 horizontal jitter of ±1.2% of the pseudotime range. The five marker panels additionally use the app's vertical jitter of ±8% of the feature range, without clipping; Beta fraction values are not vertically jittered. Consequently marker-panel jitter can extend beyond the marker measurement range; `value` and `pt_rank` retain the unjittered measurements. All points have equal display size. Disease curves are drawn from unjittered values. Each panel's points CSV includes donor/object IDs, raw values and plotted x/y; each non-INS marker curves CSV includes predictions and interval bounds. The preserved `diagnostic_beta_fraction_app_loess.csv` records the app smoother calculated during visual checks; it is not a plotted series and must not be used as the replacement panel's curve. No p values, significance stars or regression-correlation claims are added.

**Slide 5 image 7, immune heatmap.** The five marker rows are CD68, CD20, CD3e, CD4 and CD8a, in the reference order. The current app calculation is applied with 20 equal-width bins over the full fitted ranked-pseudotime range, retaining only edge-eligible display objects for bin means, takes each marker's arithmetic object mean in each bin, standardizes those 20 bin means using sample standard deviation, and clips display z values to ±2.5. The last bin includes the maximum endpoint. Every bin and object membership are exported, including unstandardized bin means, unclipped z, displayed z, object counts and donor counts. Color values are standardized marker profiles, not disease-comparison test statistics.

**Reproduction.** Run `/home/smith6jt/miniconda3/envs/scvi-env/bin/python scripts/research/pancfloyd_trajectory_20261026.py`. The DPT snapshot is reused only if its input-H5AD, Barlow-feature and app-module hashes match. No saved app data or original deck is modified. `panel_mapping.csv` maps each output to its original slide/image and exact data tables. PNG, PDF and editable-text SVG versions are saved for every panel; `output_manifest.json` records file hashes.
"""
    (OUT / "figure_legends.md").write_text(legends)
    manifest = {"source_bindings": bindings, "objects": len(fitted), "donors": fitted.case_id.nunique(),
                "figures": 9, "script_sha256": digest(Path(__file__)), "R_helper_sha256": digest(REPO / "scripts/research/pancfloyd_trajectory_20261026.R"), "fitted_objects": 1215, "corrected_edge_source_sha256": digest(edge_path), "by_status": by_status, "files": {str(p.relative_to(OUT)): digest(p)
                for p in sorted(OUT.rglob("*")) if p.is_file() and p.name != "output_manifest.json"}}
    (OUT / "output_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"output": str(OUT), "objects": len(fitted), "donors": fitted.case_id.nunique(), "panels": 9, "dpt": dpt_meta}, indent=2))


if __name__ == "__main__":
    main()
