# Islet Explorer — User Guide

Interactive web application for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Compares endocrine composition, immune cells and spatial microenvironments across ND, Aab+ and T1D donor groups; these cross-sectional observations do not establish disease progression.

**2026-11-02:** this guide covers app features. For the complete data-to-app process, use
[the workflow user guide](workflow_user_guide.md). The
[selected-component v2 build](component_build_v2_20261102.md) preserves the exact
named v19 endocrine cells and applies current methods to the remaining cells.
Production `pubapps-2026-11-components-v1` is live and verified as of
2026-11-03 02:42 UTC (November 2 local), with 15/15 authenticated-access and
49/49 browser checks passing. Consult the [deployment record](container_security_20261102.md)
for the exact image, data and validation evidence. The separate strict-DuckDB
[dev release](evidence/dev_duckdb_deployment_20261103.json) is also live and verified.
The [November 3 display fixes](app_display_fixes_20261103.md) described below are
local changes pending a new dev build; neither current live image contains them.

**Live reviewer URL (authentication required)**:
`https://islet-explorer.rc.ufl.edu/islet-explorer/`

> **Publication profile:** The enabled tabs are Plot, Trajectory, Spatial and Viewer,
> in that order. AI, Statistics, and Inflammation are disabled. The concise
> access and feature instructions are in [reviewer_guide.md](reviewer_guide.md). This longer
> guide retains descriptions of dormant development modules for maintainers; those
> sections are not claims about the reviewer deployment. The detailed measurements and
> limitations behind the current composition DPT/PAGA method are recorded in
> [the trajectory accuracy reference](../app/shiny_app/docs/trajectory.md).

---

## Table of Contents

1. [Overview](#overview)
2. [Cell Phenotypes & Marker Rules](#cell-phenotypes--marker-rules)
3. [Outlier Handling](#outlier-handling)
4. [Downloading Figures](#downloading-figures)
5. [Plot Tab](#plot-tab)
6. [Trajectory Tab](#trajectory-tab)
7. [Viewer Tab](#viewer-tab)
8. [Statistics Tab](#statistics-tab)
9. [Spatial Tab](#spatial-tab)
10. [Single-Cell Drill-Down](#single-cell-drill-down)
11. [AI Assistant](#ai-assistant)
12. [Data Pipeline, Administration & Troubleshooting](#data-pipeline-administration--troubleshooting)

---

## Overview

### Data Source

The selected-component release uses exactly Beta, Alpha, Delta and CHGA+
hormone-negative for endocrine cell totals. Islet, Cluster and Single are
anatomical object classes, not additional cell identities. Unsupported cell
assignments remain explicit as **Unclassified**, retained in overall cell
counts, phenotype controls, plots and downloads. Their grey palette entry
remains visible in legends. Existing Tissue labels are distinct and unchanged.
An object with no named endocrine core cells still has a structural record,
but no available endocrine composition or invented embedding value. Use the
method help and loaded release manifest for the exact generation.

The reviewer release validates a checksummed manifest, loads prepared islet summaries from
read-only Parquet, queries large spatial tables through embedded DuckDB, and reads the
trajectory H5AD and composition feature table as separate immutable artifacts. The
frozen v19 and v20 comparison roots both contain:

- **20 nPOD donors:** 7 ND, 7 Aab+ and 6 T1D. Donor 78 is excluded.
- **3,199,606 analysis-eligible cells**, with the same UUIDs in both versions.
- **4,848 analysed endocrine objects** in the assembled metadata: 1,468 Islets and
  3,380 Clusters. The 5,227 Singles are counted separately, outside this object-analysis
  cohort. These app counts are not the full upstream geometry inventory.
- **1,392 saved endocrine-composition embedding objects**, with at least ten core cells.
- **31 exported marker-value channels**, donor demographics, core/peri composition,
  and spatial structure measurements. Marker availability and measurement definitions
  come from the selected root; a measured channel is not necessarily a typing marker.

These counts describe the two frozen comparison artifacts, not an accepted corrected
release. The root's `uniform_release.json`, label groups and generated release manifest
define each artifact. Displayed counts also depend on the current filters and available
measurements. Phenotype choices are derived from that root rather than a fixed list of
21 labels. See [the component comparison](v19_v20_component_comparison_20261102.md).

**Historical compatibility data:** the earlier 15-donor root had 5,214 trajectory H5AD
observations and 5,382 prepared core-area/composition rows. Its 21 source labels and
19-rule CSV describe that older pipeline only; they are not the current cohort contract.

Local development retains the legacy H5AD/Excel loader for pipeline work. Production fails
closed if a required manifest artifact or hash is missing; it does not silently fall back
to a different scientific dataset.

### Disease Groups

| Group | Color | Description |
|-------|-------|-------------|
| **ND** | Steel Blue (#4477AA) | Non-diabetic control donors |
| **Aab+** | Burnt Umber (#CC6633) | Autoantibody-positive donors |
| **T1D** | Forest Green (#228833) | Type 1 diabetic donors |

Colors are configurable via the **Donor Status Colors** palette selector (Paul Tol default, Bright, Okabe-Ito, Diverging). The selector is available on the Plot sidebar, Trajectory tab, and Spatial tab — all stay synced.

---

## Cell Phenotypes & Marker Rules

Phenotype choices, group totals and marker-profile names follow the loaded root's frozen
identity tables and label groups. The current typing registry is
[`scripts/uniform/typing/registry/`](../scripts/uniform/typing/registry/); an existing frozen
root must be interpreted with its recorded policy and provenance, not silently reinterpreted
using a newer registry. Group choices such as **Endocrine (all)**, **T Cell (all)** and
**Macrophage (all)** are aggregates, not additional cell identities.

### How to read a rule

`+` and `−` describe positivity against the image's recorded marker level. Final v20
uses whole-cell positivity; CD3e, CD68, CD163 and ECAD use separate membrane measurements
when comparing competing identities. Ki67 uses the nucleus upstream, but neither frozen
v19 nor v20 exports that state. The new component build carries `Ki67_pos` through
cell, donor and distance exports with its nuclear measurement/threshold provenance;
it remains separate from identity and the 31-marker embedding panel.

Endocrine names require endocrine-object membership and marker evidence. The final v20
method names a hormone-positive cell by its unique highest scaled INS, GCG or SST value;
CHGA-positive cells without a positive hormone have the explicit **CHGA+ hormone-negative**
label. Pixel-classifier positivity and corrected whole-cell positivity are different
measurements. Retained `Islet`, `EndocrineCluster` or `EndocrineSingle` labels must not be
presented as resolved endocrine types. The unresolved Cluster/Single definitions and
inconsistent endocrine aggregates remain release-blocking discrepancies.
The owner selected v19 endocrine decisions for the combined candidate. Its app and
QuPath builds are complete, with unresolved assignments retained visibly. See
[the build record](component_build_20261102.md) and
[discrepancy trace](endocrine_classifier_discrepancy_20261102.md).

Under the v20 structure-only rule, **Endothelium Blood Vessel**, **Endothelium Lymphatic**
and **Neural** require the corresponding structure membership, more than half of the cell's
area inside it, a positive defining marker and the qualifying highest fold over level.
PDPN/CD34/CD31 positivity alone does not define a current generic Endothelial cell.
Structure area measurements and structure-qualified cell counts are separate endpoints.

### Where the rules appear in the app

- **"?" help icons** (since 2026-10-20, one per card or control panel): a card or Data / View controls panel with something to explain has a "?" beside its title, and its popover explains only that card or panel. Click the icon to open it; click elsewhere or press Escape to close it. Only one is open at a time.
- **"Phenotype rules — marker definitions" button**, at the top of each tab's Data controls help, opens the root's rule reference. Historical fallback descriptions are not authority for newer typing; use the frozen root's manifest and the workflow record when auditing a definition. The Spatial tab's phenotype choices follow its actual data labels.
- There are no inline hints under the phenotype selectors and no rules list in the drill-down panel; both were removed on 2026-10-17, when all help text moved into the "?" popovers.

### Historical rule references

`data/phenotype_rules.csv` retains the older 19-rule hierarchy and its residual categories
for compatibility with historical data. Its 15-marker list, ECAD+/APCs/Structural/Unknown
buckets and marker-only Endothelial route must not be used to explain the final v20
workflow. Current macrophage profiles encode CD68/CD163/LGALS3/CD4 states; conflicts and
unresolved start classes remain visible as recorded data rather than being renamed.

---

## Outlier Handling

Outlier controls apply to the Plot tab. (The Spatial tab's copy was removed on 2026-10-20 with the two
outlier-filtered cards it drove; every Spatial value is now donor-level.) The reviewer Trajectory deliberately has no
outlier exclusion: its cell-count cohort defines and refits the axis, so a second hidden
row-dropping rule would make the trajectory definition ambiguous. The dormant Statistics
module has its own independent sensitivity control when maintainers enable it.

### The Rule

**For each donor status (ND, Aab+, T1D), the Plot tab flags an islet when its robust |z| — the Iglewicz–Hoaglin modified z-score, `0.6745·(value − median)/MAD` — exceeds the shared slider value (default 3, range 0.5–10, step 0.5).** Median and MAD make the rule robust to right skew and two-sided. When MAD = 0, the calculation falls back to classical mean/SD within that group.

- **Dormant Statistics module**: uses its own checkbox and threshold; it never synchronizes with the reviewer-visible Plot control.

### The Control

A checkbox **"Exclude outliers"** plus an **"Outlier z-threshold"** slider lives in the Plot tab's Data controls. **Default: checkbox OFF, slider 3.0.** Trajectory has no corresponding control.

If maintainers enable **Statistics**, its dedicated checkbox and slider operate independently
of Plot. A caption reports how many rows the most recent run excluded.

### What you see when the filter is ON

- Excluded points still render on the Plot panels as **faint gray × markers** ("ghost points") so their locations remain auditable.
- The Plot tab shows an exclusion count and an optional excluded-outlier table with case, group, value and score.
- Plot ghost-marker positions are stably jittered, so toggling the control does not shuffle the cloud.

### What changes when the filter is OFF

- No points are flagged or hidden; ghost markers disappear; counts go to zero.
- The dormant Statistics module, when enabled, shows its own independent "Filter OFF" caption.
- Trend lines, donor-level means, and all statistical tests run on the full data including extreme values.

### Why this matters

Multiple endocrine objects come from each donor. A single extreme object can influence
a summary without representing an independent donor. The count badge and ghost markers
make the effect of the displayed outlier rule inspectable.

---

## Downloading Figures

Every plot in the app can be downloaded as a publication-ready PNG.

### Interactive plots (plotly — most of the app)

Hover over the plot and a small toolbar appears in the top-right of the plot. The **camera icon** ("Download plot as PNG") exports the current view at **2800 × 1800 px** — high enough resolution for slides, posters, and figure panels. Filenames are descriptive (e.g. `statistics_forest_plot.png`, `spatial_distance_kde_All_Immune_Cells.png`). The plotly logo is hidden so it never appears in your exports.

If you have zoomed or panned the plot before clicking the camera, the downloaded PNG reflects that view.

### Static plots (ggplot — large scatters, heatmaps, segmentation)

A few plots are rendered statically because they would overwhelm an interactive viewer (the 177K-cell tissue scatter, the islet segmentation map, the multi-feature pseudotime heatmap, the cell-composition bars). For each of these you'll find a **"Download PNG"** button immediately under the plot:

| Plot | Where | PNG dimensions |
|------|-------|----------------|
| Tissue scatter | Spatial tab | 3600 × 3000 px (12 × 10 in @ 300 dpi) |
| Islet segmentation map | Plot tab + Trajectory tab (when an islet is selected) | 3000 × 2100 px |
| Cell composition bars | Plot tab + Trajectory tab (when an islet is selected) | 3000 × 2100 px |
| Trajectory UMAP (donor status) | Trajectory tab | 3000 × 2100 px |
| Trajectory UMAP (selected feature) | Trajectory tab | 3000 × 2100 px |
| Multi-feature pseudotime heatmap | Trajectory tab | 3600 × variable px (height scales with marker count) |

All PNGs are 300 dpi — sized to drop directly into a manuscript or fit a full-screen slide without resampling. Filenames include a timestamp so repeated downloads don't overwrite.

### CSV downloads

Sidebar "Download CSV" buttons on the Plot, Statistics, and Spatial tabs export the underlying numeric data with the active filters (donor status, min-cells, outlier filter state) applied. The CSV header records the active outlier rule and any zero-value exclusion, so the export is self-documenting for a methods section.

### Font sizes

All plots use a consistent font baseline (16pt body text, 18pt titles, ≥12pt for the densest heatmap labels) so text remains legible whether you read it on screen or at presentation scale.

---

## Plot Tab

The Plot tab is the primary exploration interface. It has a **sidebar** (left) and **main panel** (right).

### Sidebar Controls

- **Mode**: Switch between `Scatter` (feature vs diameter) and `Composition` (cell type proportions)
- **Feature/Composition selector**: Choose what to plot
  - *Scatter mode*: Select from markers (31 proteins) or targets (hormone intensities)
  - *Composition mode*: Three option groups (the previous "Immune Metrics" group of 5 hand-curated ratios was removed in Apr 2026 — use the Spatial tab's per-phenotype enrichment / proportion / distance views instead):
    - **Phenotype Proportions** — the loaded root's phenotype groups and labels. All 20 donors in the uniform comparison roots use the recorded uniform workflow; the older original-versus-additions typing split is historical.
    - **Peri-Islet Proportions** — proportions of those retained phenotypes in the 20 μm expansion zone around the islet
    - **Core+Peri Proportions** — combined phenotype proportions using core and peri cell counts

  QuPath threshold-based INS/GCG/SST measures remain in source artifacts for provenance, but are excluded from Population Measure because the original and added donors do not use the same threshold-counting method.
- **Metric** (Composition / Microenvironment): Picks the y-axis unit. Labels spell out the denominator explicitly so "density" is never ambiguous — see [Density terminology](#density-terminology) below.
- **Region**: Filter by tissue region (core only, peri-islet, or all)
- **Donor Status**: Select disease groups to include (ND, Aab+, T1D)
- **AAb Filter**: Filter by autoantibody positivity (when available)
- **Age Range**: Slider to restrict donor age (H5AD only)
- **Sex**: Checkbox filter (H5AD only)
- **Endocrine objects** (first in the Data controls; the Spatial tab has a copy that stays in step): choose Islets and/or Clusters. Singles (whole-cell endocrine objects, not necessarily named endocrine cell identities) are not analyzed; the Endocrine objects per donor card at the bottom of the Plot tab counts them, and the segmentation map outlines them in orange. This replaced the Min cells/islet filter on 2026-10-17; every tab except Trajectory uses the selected classes (Trajectory fits every object with 10 or more core cells).
- **Diameter Range**: Restrict islet diameter in μm

### Main Panel

The pending November 3 layout update lets Plot selectors use the available card width and
reserves space beside the arrow, so long selected labels can wrap without overlapping it.
The selector values and calculations are unchanged.

- **Scatter (left)**: Interactive Plotly scatter, diameter vs selected feature, with binned summary lines per donor status. Excluded outliers (see [Outlier Handling](#outlier-handling)) appear as faint gray × ghost markers.
- **Distribution (right)**: Violin or box plot of the same feature, grouped by donor status. Has its own **independent** "Exclude zero values" checkbox — toggling it no longer affects the left scatter. Outlier ghost × markers render here too.
- **View controls**: the point controls under *Both plots* (color by, point size, transparency, show individual points, log y-axis) apply to the scatter and the distribution alike; Trend line, log2 x-axis and the summary statistic apply to the scatter, Plot type to the distribution. Individual points are shown by default (since 2026-10-19).
- **Outlier Table**: Lists excluded outliers for the selected feature. Hover the "Show excluded-outlier table" checkbox label in the sidebar for the rule.
- **Click Interaction**: Click any data point to open the segmentation viewer (see [Single-Cell Drill-Down](#single-cell-drill-down))

### Tips

- The two "Exclude zero values" checkboxes (one per panel) are independent — set them differently when you want the violin/box to include zeros but the binned scatter to exclude them, or vice-versa.
- Point positions on both panels use seeded jitter, so the cloud doesn't shuffle when you toggle other controls.
- Use Cell Populations mode with "Peri-Islet Proportions" to visualize the cellular microenvironment around each islet. Use "Core+Peri — Phenotype Proportions" to pool both zones into a single per-islet value (cells from core + peri divided by the matching pooled cell count or pooled annotation area).
- Use Spatial for per-phenotype enrichment, proportion, proximity, and distance views; the retired Plot "Immune Metrics" shortcut group is not part of the reviewer UI.

### Density terminology

The word *density* is reused in several places with different denominators. Each axis label now spells the denominator out, but for reference:

- **Cells per islet (core / peri / core+peri)** — a raw count tallied inside the islet's core, the 20 µm peri-islet zone, or both pooled together. Not divided by area. Surfaced by the Metric = "Cells per islet" choice.
- **Cells per mm² (core / peri / core+peri area)** — phenotype cells divided by the QuPath-measured annotation area in mm². Core uses `core_region_um2` (the islet boundary annotation), peri uses `peri_region_um2` (the `_exp20um` ring annotation), and core+peri uses their sum. Surfaced by Metric = "Cells per mm² …" in Cell Populations mode for all three region variants.
- **{Nerve|Capillary|Lymphatic} area (% of region area)** — the only Microenvironment-mode metric. The *fraction* of the selected region (core, peri-islet, or core+peri) occupied by structure annotation, expressed as a percent: `100 × area_um2 / region_um2`. The y-axis names whichever class and region is selected (e.g. `"Nerve area (% of core area)"`). Counts of these structures are NOT exposed — Nerve / Capillary / Lymphatic are area-valued, and a single nerve may be segmented into several pieces, so object counts depend on segmentation morphology, not biology.
- **Cells per µm of distance (Spatial Card C)** — the signed-distance KDE on Card C is a kernel density rescaled by N, so the y-axis is cells per µm of distance from the islet edge — *not* a per-area density. Useful for comparing distance distributions across disease stages; do not read it as cells/mm².
- Log-scale features (e.g., cell counts) use pseudo-log transformation so zero values remain visible at y=0.

---

## Trajectory Tab

Visualizes one composition-derived diffusion-pseudotime ordering. For the fitted cohort
of endocrine objects with at least ten core cells, the app uses log-normalized core plus 20 µm peri-islet cell-type
composition, a cosine 15-nearest-neighbor graph, and the ND islet nearest the ND centroid
as the root. The implementation uses scanpy DPT; destiny is the explicit axis fallback when
scanpy is unavailable. PAGA/ForceAtlas2 supplies the separate two-dimensional graph card.

### Layout

The Trajectory tab uses a full-width layout (no sidebar). All controls are inline.
The pending November 3 layout update aligns the pseudotime and UMAP cards, then
places the Segmentation map in a full-width row with consistent spacing above and
below. Its PAGA card title is “PAGA force-directed layout”; the method attribution
remains in the help text.

### Components

1. **Pseudotime Scatter** — The selected marker mean or cell-type/endocrine-cell ratio
   against raw or rank-normalized diffusion pseudotime. Points can be colored by donor
   status, donor, or a retained H5AD Leiden assignment.

2. **Pseudotime Strip (Dominant Status × Feature Intensity)** — A single horizontal strip below the UMAP scatter, divided into 25 pseudotime bins. Each bin encodes two variables at once:
   - **Hue** = the *dominant* donor status in that bin, computed from per-status normalized density (count-in-bin divided by total islets for that status). Normalization corrects for cohort imbalance so dominance reflects enrichment relative to each status's pool, not raw count.
   - **Intensity** (saturation toward white) = the bin-mean of the selected feature, robust 5–95th percentile clipped. Floor at 50% saturation so even the dimmest bin remains a clear tint, never white.
   Hover for raw counts per status, normalized density per status, and the feature mean.

3. **Multi-Feature Heatmap** — Z-scored expression of user-selected markers along pseudotime bins. Default markers: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR. Blue-white-red diverging colormap, z-scores clamped to [-2.5, 2.5].

**Y-axis scale on the scatter** — The frozen v19/v20 roots export per-image marker values
scaled as `(x − level) / (p99.9 − level)`, clipped to `[0, 1]`, then averaged over the object's
core cells. Their measurement and level provenance is in `uniform_release.json`; these
are not scVI values or raw fluorescence. The multi-feature heatmap separately standardizes
each marker across pseudotime bins. The owner's acquisition restriction means CODEX INS
intensity is not an approved quantitative between-donor hypothesis endpoint; see
[the measurement review](ins_measurement_review_20261026.md).

4. **UMAP Cards** — The islet-level UMAP (on a uniform root, the endocrine-composition layout the
   Spatial tab shows, for the same ten-cell objects) is shown once by donor status and once by the
   selected Y measure. These coordinates provide a fixed reference and are not the PAGA
   layout.

5. **PAGA Force-Directed Layout** — A separately computed two-dimensional graph initialized
   from PAGA and drawn with ForceAtlas2. Coordinates are unitless; interpret relative
   neighborhood structure, not axis magnitude.

6. **Click Interaction** — Click a pseudotime point to:
   - Highlight the selected islet
   - Open the segmentation viewer panel
   - Jump to the OME-TIFF Viewer tab for that donor

### Controls

- **Y axis**: Switch between a protein marker mean (31-marker panel) and
  `log1p(cells of the selected type / endocrine cells in the islet core)`; ratio counts
  can use core, peri, or core + peri cells.
- **Trend lines**: None, Overall (single LOESS), or By Donor Status (ND/Aab+/T1D separate LOESS curves)
- **Color points by**: Donor Status (default), Donor ID, or Leiden cluster at any of 4 resolutions (0.3 / 0.5 / 0.8 / 1.0)
- **Pseudotime scale**: Rank-normalized (default) spreads the fitted cohort uniformly;
  Raw DPT preserves its original spacing. Rank transformation can change apparent curve
  shape, so use raw DPT when density along the manifold matters.
- **Point size by**: Cell Count (default), Islet Diameter, or Uniform
  - *Cell Count* sizes points by `sqrt(total_cells)` — small islets (3-10 cells) appear as tiny dots, large islets (100+ cells) as large dots. This honestly represents measurement quality.
- **Point transparency**: 0 (fully opaque) to 0.95 (highly transparent). Default 0.4. Moving the slider right increases transparency.
- **Point size**: Base point size slider (0.5-5.0)
- **Fitted objects** (fixed, no control): every endocrine object with 10 or more core
  cells, whatever its class (Islet or Cluster), following Barlow et al.'s
  threshold. The Endocrine objects selector of the other tabs is hidden here.

### Cell-Count-Weighted Trends

Trend lines use `log1p(cell count)` weights. Inspect the cell count in each point's tooltip;
the current fitted cohort includes only objects with at least ten core cells. Weighting and
the cohort definition are analytical choices, not evidence of biological progression.

### Interpretation Limits

Pseudotime is an ordering of compositional similarity, not elapsed time or proof of a
causal ND → Aab+ → T1D transition. Multiple objects from one donor are not independent
donor replicates. Legacy 15-donor trajectory findings do not validate the current 20-donor
candidate or its changed cell identities. Treat marker
trends as exploratory and preserve the displayed cell range, engine, scale, and donor
context in any exported interpretation.

---

## Viewer Tab

Embedded OME-TIFF viewer (Avivator) for full-resolution multiplexed imaging data.

### Usage

1. Select a donor case from the dropdown, or click an islet in Plot/Trajectory to auto-navigate
2. The viewer loads the corresponding OME-TIFF with configurable channel overlays
3. Channel controls allow toggling individual markers and adjusting intensity ranges

### Requirements

- OME-TIFF files must be accessible at the configured URL
- Avivator runs client-side in the browser (WebGL required)

> **macOS note (fixed Jun 2026):** the Viewer previously flickered/reloaded continuously on macOS because the embedded viewer was recreated on every update. It now loads steadily, and switching images updates the existing viewer in place. See the [Administration & Operations Guide](admin_guide.md#troubleshooting) if flickering recurs.

---

## Statistics Tab

**Dormant historical module:** this tab is disabled in the publication profile. The
controls below document its implementation, not a current reviewer feature or acceptance
of its historical analyses.

Rigorous statistical analysis of the currently selected feature from the Plot sidebar. The Statistics tab **shares the Plot sidebar** — changing the feature, filters, or disease groups in the sidebar updates both Plot and Statistics simultaneously.

An **"About this tab"** banner summarizes the module and its donor-aware mixed-effects
model. Historical hard-coded cohort counts are not the current release inventory.

The tab uses a **5-section narrative layout** that walks you through a logical analytical flow, with numbered headings and brief explanations at each step.

### Section 1: Configure Analysis

- Overview banner showing N islets, global p-value, effect size (η²)
- **Run Statistics** button triggers computation; **Download CSV** exports results
- Test type selector with inline explanations:
  - **Parametric** (ANOVA / t-test): assumes roughly normal data within groups
  - **Non-parametric** (Kruskal-Wallis / Wilcoxon): no normality assumption; based on ranks
- Additional controls: significance level (α), min cells/islet, bin width, diameter range. **Outlier removal** is controlled by a dedicated "Exclude outliers" checkbox + "Outlier z-threshold" slider local to the Statistics tab (off by default, slider default 3.0). It is independent of the Plot/Spatial control. A caption reports the most recent exclusion count (or a "Filter OFF" notice). See [Outlier Handling](#outlier-handling).
- **Exclude zero-valued islets** (checkbox, default OFF): filters `value == 0` rows before every Statistics-tab test (mixed-effects, donor-level ANOVA / Kruskal, pairwise Cohen's d, per-bin tests, demographics regressions, covariate-adjusted model). Useful for zero-inflated features — rare phenotype proportions, peri-islet enrichment z-scores, or any metric where many islets register exactly 0. The overview banner shows an orange "Zero-valued excluded (N dropped)" badge when active; the forest-plot title and the Methods Reference paragraph cite the same count. **AUC is read from the Plot tab's binned summary and is not affected by this filter** — toggling re-runs the hypothesis-testing pipeline but leaves the AUC card untouched.

### Section 2: Primary Results

Three equal-height cards showing the core statistical results:

- **Hypothesis Testing** (left): Global test + pairwise comparisons table with Cohen's d effect sizes, 95% CIs, and BH-corrected p-values. The primary "Overall" row uses a mixed-effects model with donor as random intercept. The "Donor-level means" row averages islets within each donor as a conservative cross-check.
- **Effect Size Forest Plot** (center): Visual display of pairwise Cohen's d with confidence intervals. Points colored by significance.
- **Area Under Curve** (right): Trapezoidal AUC by disease group with percentage change relative to ND baseline.

### Section 3: Size-Dependent Patterns

Stratified analysis testing whether group differences depend on islet size (effect modification):

- **Stratified Tests by Islet Diameter**: Heatmap of BH-corrected q-values (ANOVA and Kendall τ) computed within each diameter bin. Identifies which size ranges drive or lack group differences. P-values are corrected across bins to control false discovery rate.
- **Trend Analysis (Kendall τ)**: Measures correlation between disease stage (ND=0, Aab+=1, T1D=2) and the feature value within each size bin. τ > 0 means the feature increases with disease progression; τ < 0 means it decreases. Points colored by significance (red = significant, gray = NS).
- **Bin coverage rule**: The heatmap and the trend plot use the same rule — a bin is tested only if at least 2 donor-status groups each have ≥ 2 donors with finite values after filtering. Bins that fail this rule are silently omitted from both plots so the coverage stays in sync.

### Section 4: Confounders & Deeper Analysis

Full-width **Demographics** card (H5AD only, hidden for Excel fallback):

- **Donor Summary Table**: N islets (primary), N donors, age median/range, % male/female per disease group
- **Age vs Feature** (left): Scatter plot of eligible prepared islet rows after filtering, with overall linear regression line and Pearson correlation. Subtitle notes that islets are correlated within donors.
- **Sex vs Feature** (right): Box plot of feature value by donor status, faceted by sex (Male/Female). Sex-stratified pairwise tests are **BH-corrected across all 6 contrasts** (3 pairs × 2 sexes); the subtitle reports the minimum adjusted q per sex (low-power sex × group cells, n=2-3 donors, are flagged but still corrected jointly).
- **Autoantibody Profile** (Aab+ only, hidden when no AAb data): Per-donor table showing which AAbs are positive, total AAb count, N islets, and mean feature value. Box plot of feature by number of positive autoantibodies (1, 2, 3+) with Kruskal-Wallis test.
- **Covariate-Adjusted Model**: Donor-level linear model adjusting for age and sex, testing whether donor status retains significance.

### Section 5: Methods Reference

Dynamic text (subdued gray background) describing all statistical tests, corrections, and assumptions. Adjusts based on test type and parameters selected in Section 1.

### Using the Statistics Tab

1. Select a feature in the Plot sidebar (e.g., `prop_CD8a Tcell` or `peri_prop_Macrophage`)
2. Switch to the Statistics tab
3. Adjust test parameters in Section 1 (Configure Analysis)
4. Click **Run Statistics** — all sections populate
5. Review the narrative flow: global test → effect sizes → size-dependent patterns → confounders

---

## Spatial Tab

Tissue-wide spatial visualization and peri-islet microenvironment analysis. Combines single-cell tissue scatter plots with Leiden clustering and neighborhood enrichment metrics.

The frozen v19/v20 roots each record 1,392 saved embedding objects from the 20-donor
cohort. Their UMAP and Leiden use endocrine composition for objects with at least ten
core cells; donor status colors the saved coordinates. The earlier 21-donor schema-4
experiment and its 7,493/11,136 object counts are superseded and are not this dataset.
Changing identities requires rebuilding these derived outputs from a coherent candidate;
a diagnostic identity replay cannot be substituted for an assembled root.

### 3-Panel Layout

1. **Controls Sidebar** (left panel)
   - Donor selector (the loaded root's donors, with disease status labels)
   - Color by: Phenotype (root-specific labels) or Leiden cluster
   - Leiden resolution dropdown (0.3, 0.5, 0.8, 1.0) — visible when Leiden coloring selected
   - Region filter: All cells, Core + Peri only, or Core only
   - **Color background cells**: Toggle to show tissue background cells in phenotype colors (dimmed) instead of gray
   - **Phenotype filter**: Show/hide individual phenotypes with checkboxes. Use **All** / **None** links to quickly select or clear all types. Dynamically populated from the selected donor's cell types.
   - Donor status checkboxes
   - Phenotype and Donor Status palette selectors
   - Download CSV button

2. **Tissue Scatter** (center, 800px height)
   - Full tissue-wide scatter plot of ~177K cells per donor
   - Background: tissue cells in light gray by default, or colored by phenotype when "Color background cells" is checked
   - Foreground: core/peri cells colored by phenotype or Leiden cluster
   - Spatial coordinate convention: `coord_fixed()` + `scale_y_reverse()` (microscopy standard)
   - **Zoom**: Drag to select a region, then the view zooms in. Double-click or click **Reset Zoom** to restore the full view.
   - Uses ggplot2 `renderPlot` (not plotly) for performance at >100K points

3. **Leiden clusters** (right): **Show** switches between three views. *Composition* (the default) is a stacked
   bar chart of each cluster's mean core phenotype proportions, clusters in numeric order: what each cluster is. Its
   endocrine bars restate the proportions the clusters are built from. *Donor status* and *Donor* (since 2026-10-20) split
   each cluster's objects by donor status or by donor, with each cluster's donor count above its bar: which
   composition states are disease-associated, and whether a cluster belongs to a few donors.
4. **Islet-Level UMAPs** (full-width row below, since 2026-10-19): the Leiden UMAP and, since 2026-10-20, the same
   layout colored by each object's immune proportion in its 20 µm band, side by side at one size, both interactive (PNG
   from each camera button). The layout by donor status is on the Trajectory tab.
   - Built from each object's endocrine composition (the Beta, Alpha, Delta and CHGA+ hormone-negative proportions of
     its core endocrine cells) for the objects with ten or more core cells, with no donor correction (2026-10-19).
     Smaller objects have no point or cluster in this saved layout. Color patterns are
     descriptive; the layout alone does not establish disease progression or classification accuracy.
   - Hidden with message when Leiden data unavailable in H5AD

The pending November 3 display update places the tissue-scatter legend and both
UMAP legends or continuous color bars below their plots. Dense tissue labels wrap,
and the plot grows vertically to preserve the map area while keeping downloads
below it. Coloring choices, linked selection and exported values are unchanged.

### Data Sources

- **Donor cell tables**: the selected source root's `donors/{imageid}.csv`, or their manifest-backed production Parquets. Columns include X/Y centroids (μm), phenotype, cell_region (core/peri/tissue), and islet_name. The frozen v19/v20 roots contain 20 donors.
- **Leiden clustering**: Stored in `islet_explorer.h5ad` .obs as `leiden_0.3`, `leiden_0.5`, `leiden_0.8`, `leiden_1.0` + UMAP coords `leiden_umap_1`, `leiden_umap_2`.
- **Neighborhood metrics**: the selected root's per-object core/peri composition, enrichment and distance artifacts. Available fields and denominators follow that root and its manifest.

### Neighborhood Analysis Cards

Rebuilt on 2026-10-20 after an owner review that asked each card to answer its own question. Every value is
donor-level: each donor counts once, whatever its number of objects. The cards use the donor statuses and endocrine
objects selected in the tab's Data controls (the immune gradient reads per-cell data, so the object selection does
not apply to it). The section's own control panels and outlier controls were removed; each card carries its own
selector. Cards only appear when neighborhood data is available in the H5AD.

**Insulitis prevalence**: each donor's proportion of objects with at least the chosen number of immune cells (every
immune label, subtyped or not) in the core and the 20 µm band; default 15, after nPOD's 15 CD45+ cells. Bars are
donor means ± SEM with each donor as a point, and a Kruskal-Wallis p on the donor proportions. Small objects rarely reach
the threshold, so uncheck Cluster to count Islets only. On v7, Islets only: donor medians ND 0, Aab+ 0.074, T1D
0.185 (p = 0.01).

**Immune cells and beta-cell content**: for each donor, the mean immune proportion of the band around its objects that
still have beta cells and around its insulin-deficient objects (beta cells under the chosen proportion of the core's
endocrine cells, default 10%), joined by a line. On v7 Islets, in every T1D donor with both kinds the band around
islets with beta cells holds more immune cells (0.17–0.27) than around insulin-deficient islets (0.05–0.13).

**Immune gradient around islets**: for each donor and distance bin from the nearest endocrine object's border
(Core, 0–20 … 200–400, >400 µm), the chosen immune cells' proportion of all cells there divided by the donor's own tissue
proportion (its cells outside every endocrine object and band). 1 is the donor's tissue level. Lines are the donor median
per status with the IQR. Dividing by the donor's own tissue separates gathering at islets from a donor having more
immune cells everywhere. The donor chosen in the Data controls (the tissue scatter's donor) is drawn as a dotted line over
the medians; hover gives its cell counts, and a point resting on fewer than 20 of the chosen cells is faded.

**Phenotype composition & enrichment**: every immune phenotype at once (CD8a and CD4 T cells, B cells,
macrophages; the endocrine and exocrine classes are left out since 2026-10-20). Left, the status summary (median and IQR,
or mean and SEM, of donor means); right, one row per donor grouped by status, so the spread behind each bar shows.
- *Region*: Peri-islet enrichment z (diverging colors), Core proportion, Peri-islet proportion.
- *Baseline* (enrichment z only): **Cohort** compares each band with the phenotype's proportion of all analyzed cells in
  the cohort, so it measures how many are near islets, and a donor with more of a cell type everywhere scores high.
  **Own tissue** compares each band with the donor's own tissue proportion, so it measures gathering at islets. On v7
  the cohort-baseline donor means rank-correlate 0.80–0.91 with the donors' tissue proportions.
- *Drop object z beyond ±5* removes those objects' values before each donor mean.

**Distance-to-islet profile**: the per-cell view of every cell's signed distance to the nearest islet border across
the tissue (symlog axis), as density curves or composition bars, by immune phenotype or by marker positivity (all 31
panel markers from data root v8, each measured in every donor; positive = 0.10 or above on the cohort marker scale).
Both views stop 200 µm outside the border (composition bars 5 µm wide, density in
160 symlog bins), with a dotted line at the 20 µm band edge. In composition bars the Proportion y-axis is the
proportion of all cells at that distance that are of the phenotype or positive for the marker; each phenotype or
marker has its own row and y-axis, and its ND, Aab+ and T1D columns use that one scale, so a rare marker is not
dwarfed by a common one. Until 2026-10-20 the bars divided by the status's total cells, so they followed how many
cells lay at each distance. Every analyzed cell counts, including those without a cell subtype. **Distance to** (data root v9):
*Nearest Islet* (the default) measures each cell to the nearest Islet's border; *Nearest Islet or Cluster* also counts
Clusters, which are small and scattered through the tissue, so most cells beyond 20 µm are then measured to a Cluster
and the profile around Islets blurs (T1D CD3e+ cells at 0–5 µm: 6.1% to the nearest Islet, 5.0% to any object).
Beyond about 30 µm the proportions mostly level off at the tissue level, with a slow fall for immune markers over
hundreds of µm. Every cell
counts once, so donors with more cells weigh more; the immune gradient is the donor-fair version. Its former
donor-averaged KDE view was removed: its table holds only core and band cells, so it stopped at 20 µm.

*Removed on 2026-10-20*: Peri enrichment by status (one column of the composition card), Peri vs core proportion
(immune cells almost never lie in the core and endocrine cells almost never in the band, so the points sat on one
axis), and Distance to nearest cell (measured from the core's center, so it tracked object radius, rank correlation
0.78–0.86; its T1D difference disappeared once radius was subtracted).

### Data Coverage

- Coverage is determined by the selected root, matching UUID/object mappings and the
  displayed filters. Use the manifest and exported row counts rather than historical
  5,214-row totals or old donor-specific counts.
- A disease-group ordering in a plot is an observation, not a validation of typing,
  segmentation or a disease-progression claim.

---

## Single-Cell Drill-Down

Click any islet data point in the Plot or Trajectory tab to inspect individual cells.

### How It Works

1. Click a point in the scatter or UMAP (Trajectory), or in either plot on the Plot tab. Hovering a point first rings the same object in the tab's other plots; **Clear selection** releases a selected object.
2. A segmentation panel appears below the plot with the title showing the islet name and donor info (case ID, disease status, age, sex)

> **macOS note (fixed Jun 2026):** clicking a point previously did nothing on macOS. Drill-down now works on all platforms via a native-click handler. Tip: on the Plot tab, keep **Show individual points** on (the default) so there are per-islet points to click — aggregated summary points are not clickable.

### Panel Layout

- **Left sidebar**: Layer toggles (Single Cells, Peri Boundary, Structures), Color-by dropdown, Phenotype Palette, Show peri-islet cells checkbox. When no cell data is available, shows a boundary legend instead.
- **Center**: Segmentation map with GeoJSON boundaries and optional single-cell overlay
- **Right sidebar**: Cell composition bar chart and cell count table (visible when Single Cells is enabled)

### Single Cells View

- **Color by**: the root's phenotype labels or available protein markers (viridis inferno continuous scale)
- **Show peri-islet**: Toggle checkbox to include/exclude cells from the 20 μm expansion zone
  - Core cells: filled circles
  - Peri-islet cells: open circles (shape 1)
- **Summary panel**: Horizontal bar chart showing phenotype composition of the islet
- **Count table**: Core vs peri cell counts by phenotype

### Data Source

Source-root `cells/` files and the corresponding production spatial tables provide the
selected object's cells. Their inventory is release-specific. Fields include:
- X/Y centroid coordinates (μm, converted to pixel space for GeoJSON overlay)
- Phenotype label from the frozen identity table
- Region (core or peri-islet)
- Cell/Nucleus area
- 31 protein marker expression values

### Availability

- Availability for a selected row depends on its identifier mapping and the released
  cell/geometry artifacts; no fixed historical CSV count establishes completeness.
- If cell data is unavailable for a clicked islet, the panel shows "Boundaries" mode with a message

---

## AI Assistant

**Dormant development module:** AI is disabled in the publication profile and production
has no AI API credential. The following is historical local-module documentation, not a
current reviewer capability or authorization to send scientific data to an AI service.

An embedded AI chat panel powered by the University of Florida Navigator AI Toolkit. Ask questions about your data, plots, statistics, or troubleshooting steps.

### How to Use

1. Type a question in the text area at the bottom of the AI panel (right side of the app)
2. Select a model:
   - **Navigator Fast (gpt-oss-20b)** — Quick responses for simple questions (default)
   - **Navigator Large (gpt-oss-120b)** — More capable reasoning model for complex analysis questions. Shows "Thinking..." while reasoning, then displays the answer.
3. Click **Send** or press Enter
4. Use **New Conversation** to clear chat history and start fresh

### Model Selection Tips

- Use **Fast** for quick lookups: "What markers are in the heatmap?", "How many islets are shown?"
- Use **Large** for analytical questions: "Why does immune infiltration increase in T1D?", "How should I interpret this AUC result?"
- The Large model takes longer to respond because it reasons through the problem before answering

### Token Budget

The app tracks cumulative token usage. If a budget limit is configured (via `OPENAI_TOKEN_BUDGET`), the chat will notify you when the budget is reached. Start a new browser session to reset.

### Requirements

The AI assistant requires:
- `KEY` and `BASE` environment variables configured in `.Renviron`
- The `httr2` R package installed
- Network access to the UF Navigator API endpoint

If credentials are missing, the chat panel displays a configuration message instead of responding.

---


## Data Pipeline, Administration & Troubleshooting

Rebuilding the app data, the full data lineage, deployment/operations, and troubleshooting (data loading, pipeline, deployment, the AI assistant, and the macOS drill-down/Viewer fixes) now live in the **[Administration & Operations Guide](admin_guide.md)**.

Quick pointers:
- **Rebuild app data**: follow the [workflow user guide](workflow_user_guide.md) and selected component provenance. A structure-only update binds its baseline calibration and runs the automatic invariance audit. The old `reaggregate_islets.py` chain is historical compatibility processing, not the current uniform workflow.
- **Local development**: `cd app/shiny_app && Rscript -e 'shiny::runApp(".", port = 7777)'`.
- **App won't load, missing filters, blank panels, or browser-specific issues**: see the admin guide's [Troubleshooting](admin_guide.md#troubleshooting) section.
