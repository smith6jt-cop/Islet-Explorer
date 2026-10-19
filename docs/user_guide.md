# Islet Explorer — User Guide

Interactive web application for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss, immune infiltration, and cellular microenvironment changes across the ND → Aab+ → T1D disease progression.

**Requested reviewer URL (use only after UFRC confirms production access)**:
`https://islet-explorer.rc.ufl.edu/islet-explorer/`

> **Publication profile (August 2026):** Reviewers see Plot, Trajectory, Viewer,
> and Spatial only. AI, Statistics, and Inflammation are disabled. The concise,
> release-current instructions are in [reviewer_guide.md](reviewer_guide.md). This longer
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

The reviewer release validates a checksummed manifest, loads prepared islet summaries from
read-only Parquet, queries large spatial tables through embedded DuckDB, and reads the
trajectory H5AD and composition feature table as separate immutable artifacts. The
underlying cohort contains:

- **15 nPOD donors** (ND, Aab+, T1D); the trajectory H5AD has **5,214 observations**,
  while the prepared core-area and composition tables each have **5,382 rows**
- **31 protein markers** (mean expression per islet)
- **Trajectory source coordinates and marker values**, plus a separately released
  phenotype-composition matrix used to refit diffusion pseudotime in the app
- **21 upstream phenotype labels**, of which 12 are retained in reviewer-facing selectors
- **Donor demographics** (age, sex, autoantibody status)
- **62 neighborhood metrics** (peri-islet composition, immune infiltration, enrichment z-scores, distances)
- **QuPath groovy data** (targets, markers, composition, LGALS3)

The source preparation measured 48,438 target rows, 64,584 marker rows, and 21
paired core/peri feature definitions. These are distinct table counts; filters
and available measurements determine the count in each view. The release
manifest is authoritative for the frozen reviewer data.

Local development retains the legacy H5AD/Excel loader for pipeline work. Production fails
closed if a required manifest artifact or hash is missing; it does not silently fall back
to a different scientific dataset.

### Disease Groups

| Group | Color | Description |
|-------|-------|-------------|
| **ND** | Steel Blue (#4477AA) | Non-diabetic control donors |
| **Aab+** | Burnt Umber (#CC6633) | Autoantibody-positive (pre-T1D) |
| **T1D** | Forest Green (#228833) | Type 1 diabetic donors |

Colors are configurable via the **Donor Status Colors** palette selector (Paul Tol default, Bright, Okabe-Ito, Diverging). The selector is available on the Plot sidebar, Trajectory tab, and Spatial tab — all stay synced.

---

## Cell Phenotypes & Marker Rules

The upstream data carry 21 phenotype labels. The canonical gating CSV contains 19
hierarchical phenotype rules; four additional residual/aggregate definitions are documented
in the app. Nine source labels are retired from reviewer-facing selectors, legends, and
tests, leaving 12 active phenotype choices. The broader rule reference remains visible for
provenance and interpretation of source artifacts.

### How to read a rule

Each phenotype is defined by which markers must be **positive** (`+`) and which must be **negative** (`−`). Rules are hierarchical: a sub-phenotype inherits its parent's constraints. For example, the displayed rule for `CD8 T cell` is the union of its row (`CD8a+`, `CD4−`) and its parent `T cell` row (`CD3e+`, `CD20−`, `CD68−`, `CD163−`, …) and so on. Some categories use **any-of** logic (e.g. `Endothelial` requires at least one of `PDPN+ / CD34+ / CD31+`); the rule text spells that out explicitly when it applies.

### Where the rules appear in the app

- **"Phenotype rules — marker definitions" button** at the top of each enabled tab opens a reference table with the 19 CSV rules and 4 explicitly described residual/aggregate definitions. It is a provenance reference, so it is broader than the 12 labels offered in reviewer selectors.
- **Inline hint under every phenotype selector** — when you change the "Population Measure" (Plot tab) or any of the Spatial Card A / Card C phenotype dropdowns, a small blue box appears underneath showing the active phenotype's rule (e.g. `CD8 T cell` → `CD3e+, CD8a+  —  SMA−, CD20−, CD68−, CD163−, CD4−, …`).
- **Drill-down composition panel** (right rail of the islet viewer on Plot + Trajectory tabs) — under the per-islet cell count table, an expandable "Marker rules for phenotypes in this islet" section lists every phenotype actually present in the selected islet, colour-coded to match the bar chart, with its full rule.

### The 15 markers used by the gating CSV

`INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN, CD34, CD31, B3TUBB`.
The H5AD/scVI panel contains 31 markers; additional markers are available for
visualisation but are not part of the gating CSV. An upstream template once listed
`CD99` and `CD56`, but neither exists in this 31-plex panel, so those phantom columns
were removed from the app's canonical rule file.

### Residual categories

Four source-artifact definitions are not formally gated by the CSV. They remain in the
reference modal for provenance but are retired from reviewer-facing selectors:

- **ECAD+** — epithelial cells (E-cadherin+); residual ductal/exocrine bucket from the upstream single-cell pipeline.
- **APCs** — antigen-presenting cells aggregate (HLADR+ myeloid lineage). Used by neighborhood metrics.
- **Structural** — fallback for fibroblast/stroma-like cells negative for every endocrine and immune marker.
- **Unknown** — cells matching no rule above.

These are described as such inside the in-app modal so users know when a label is a residual bucket rather than a strictly gated phenotype.

---

## Outlier Handling

Outlier controls apply to Plot and Spatial. The reviewer Trajectory deliberately has no
outlier exclusion: its cell-count cohort defines and refits the axis, so a second hidden
row-dropping rule would make the trajectory definition ambiguous. The dormant Statistics
module has its own independent sensitivity control when maintainers enable it.

### The Rule

**For each donor status (ND, Aab+, T1D), Plot and Spatial Card C flag an islet when its robust |z| — the Iglewicz–Hoaglin modified z-score, `0.6745·(value − median)/MAD` — exceeds the shared slider value (default 3, range 0.5–10, step 0.5).** Median and MAD make the rule robust to right skew and two-sided. When MAD = 0, the calculation falls back to classical mean/SD within that group. Spatial Card C applies it to `log(1 + distance)`.

- **Spatial Card A (Peri-Islet Phenotype Enrichment)**: uses **Tukey 1.5×IQR** per donor status, computed on non-zero values only. Many phenotypes are absent from most islets, so a z-score rule would over-flag the zero mode; the IQR variant tolerates that.
- **Dormant Statistics module**: uses its own checkbox and threshold; it never synchronizes with the reviewer-visible Plot/Spatial control.

### The Control

A checkbox **"Exclude outliers"** plus an **"Outlier z-threshold"** slider lives on the Plot sidebar and is mirrored in the Spatial toolbar. The two views stay synchronized. **Default: checkbox OFF, slider 3.0.** Trajectory has no corresponding control.

If maintainers enable **Statistics**, its dedicated checkbox and slider operate independently
of Plot/Spatial. A caption reports how many rows the most recent run excluded.

### What you see when the filter is ON

- Excluded points still render on the Plot panels and Spatial Card C as **faint grey × markers** ("ghost points") so their locations remain auditable.
- Plot, Spatial Card A, and Spatial Card C show an exclusion count using the rule that applies to that view.
- Plot and Spatial Cards A/C provide optional excluded-outlier tables with case, group, value, and score/fence information.
- Plot ghost-marker positions are stably jittered, so toggling the control does not shuffle the cloud.

### What changes when the filter is OFF

- No points are flagged or hidden; ghost markers disappear; counts go to zero.
- The dormant Statistics module, when enabled, shows its own independent "Filter OFF" caption.
- Trend lines, donor-level means, and all statistical tests run on the full data including extreme values.

### Why this matters

With thousands of islet measurements from only 15 donors, a single extreme islet within one donor can dominate the group mean and pull mixed-effects estimates. Keeping the rule visible (count badge + ghost × markers) means you can audit at a glance whether a striking effect is being driven by one or two unusual islets.

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
    - **Phenotype Proportions** — the 12 retained reviewer-facing phenotypes from the 21-label source artifact. The original 15 retain historical labels; five additions use reviewed reconstructed rules, so pooled phenotype comparisons span both methods.
    - **Peri-Islet Proportions** — proportions of those retained phenotypes in the 20 μm expansion zone around the islet
    - **Core+Peri Proportions** — combined phenotype proportions using core and peri cell counts

  QuPath threshold-based INS/GCG/SST measures remain in source artifacts for provenance, but are excluded from Population Measure because the original and added donors do not use the same threshold-counting method.
- **Metric** (Composition / Microenvironment): Picks the y-axis unit. Labels spell out the denominator explicitly so "density" is never ambiguous — see [Density terminology](#density-terminology) below.
- **Region**: Filter by tissue region (core only, peri-islet, or all)
- **Donor Status**: Select disease groups to include (ND, Aab+, T1D)
- **AAb Filter**: Filter by autoantibody positivity (when available)
- **Age Range**: Slider to restrict donor age (H5AD only)
- **Sex**: Checkbox filter (H5AD only)
- **Endocrine objects** (first in the Data controls; the Spatial tab has a copy that stays in step): choose Islets and/or Clusters. Singles (one endocrine cell) are not analysed; the Endocrine objects per donor card at the bottom of the Plot tab counts them, and the segmentation map outlines them in orange. This replaced the Min cells/islet filter on 2026-10-17; every tab except Trajectory uses the selected classes (Trajectory fits every object with 10 or more core cells).
- **Diameter Range**: Restrict islet diameter in μm

### Main Panel

- **Scatter (left)**: Interactive Plotly scatter, diameter vs selected feature, with binned summary lines per donor status. Excluded outliers (see [Outlier Handling](#outlier-handling)) appear as faint grey × ghost markers.
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

Visualizes one composition-derived diffusion-pseudotime ordering. For each selected
cell-count cohort, the app uses log-normalized core plus 20 µm peri-islet cell-type
composition, a cosine 15-nearest-neighbor graph, and the ND islet nearest the ND centroid
as the root. Production uses scanpy DPT; destiny is the explicit axis fallback when
scanpy is unavailable. PAGA/ForceAtlas2 supplies the separate two-dimensional graph card.

### Layout

The Trajectory tab uses a full-width layout (no sidebar). All controls are inline.

### Components

1. **Pseudotime Scatter** — The selected marker mean or cell-type/endocrine-cell ratio
   against raw or rank-normalized diffusion pseudotime. Points can be colored by donor
   status, donor, or a retained H5AD Leiden assignment.

2. **Pseudotime Strip (Dominant Status × Feature Intensity)** — A single horizontal strip below the UMAP scatter, divided into 25 pseudotime bins. Each bin encodes two variables at once:
   - **Hue** = the *dominant* donor status in that bin, computed from per-status normalized density (count-in-bin divided by total islets for that status). Normalization corrects for cohort imbalance so dominance reflects enrichment relative to each status's pool, not raw count.
   - **Intensity** (saturation toward white) = the bin-mean of the selected feature, robust 5–95th percentile clipped. Floor at 50% saturation so even the dimmest bin remains a clear tint, never white.
   Hover for raw counts per status, normalized density per status, and the feature mean.

3. **Multi-Feature Heatmap** — Z-scored expression of user-selected markers along pseudotime bins. Default markers: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR. Blue-white-red diverging colormap, z-scores clamped to [-2.5, 2.5].

**Y-axis scale on the scatter** — Marker values are **scVI-normalized expression** (bounded roughly `0–1`), aggregated as the **mean across the cells in each islet**. They are *not* raw fluorescence intensity and *not* z-scores. Values almost never reach zero because (1) the scVI output for an individual cell rarely sits at exactly 0, and (2) averaging the cells in an islet smooths out any remaining zeros. The **multi-feature heatmap** below the scatter, by contrast, IS z-scored per marker across pseudotime bins (intentional — it standardizes markers for cross-marker comparison).

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

Trend lines are weighted by `log1p(cell count)` so that well-measured islets (many cells) drive the biological signal while noisy small islets (few cells) contribute proportionally less. This is important because 60% of islets have ≤10 cells — their aggregated expression values are inherently noisier.

Hover over any point to see its cell count in the tooltip (e.g., "Cells: 3"). Small islets at unexpected positions (e.g., ND islets at high pseudotime) typically have very few cells.

### Interpretation Limits

Pseudotime is an ordering of compositional similarity, not elapsed time or proof of a
causal ND → Aab+ → T1D transition. Donor identity is a major source of variation in this
15-donor cohort. The audit found no donor-level separation of ND and Aab+ on the legacy
axes, so avoid presenting a clean, evenly spaced three-stage progression. Treat marker
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

Rigorous statistical analysis of the currently selected feature from the Plot sidebar. The Statistics tab **shares the Plot sidebar** — changing the feature, filters, or disease groups in the sidebar updates both Plot and Statistics simultaneously.

An **"About this tab"** banner at the top summarises what the tab does, the donor cohort (N=15), and the primary test (a donor-aware mixed-effects model with donor as a random intercept). Each section heading carries a small **(?)** icon — hover it for a one-line description of that section's contents.

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

- **Hypothesis Testing** (left): Global test + pairwise comparisons table with Cohen's d effect sizes, 95% CIs, and BH-corrected p-values. The primary "Overall" row uses a mixed-effects model with donor as random intercept. The "Donor-level means" row averages islets within each donor (N=15) as a conservative cross-check.
- **Effect Size Forest Plot** (center): Visual display of pairwise Cohen's d with confidence intervals. Points colored by significance.
- **Area Under Curve** (right): Trapezoidal AUC by disease group with percentage change relative to ND baseline.

### Section 3: Size-Dependent Patterns

Stratified analysis testing whether group differences depend on islet size (effect modification):

- **Stratified Tests by Islet Diameter**: Heatmap of BH-corrected q-values (ANOVA and Kendall τ) computed within each diameter bin. Identifies which size ranges drive or lack group differences. P-values are corrected across bins to control false discovery rate.
- **Trend Analysis (Kendall τ)**: Measures correlation between disease stage (ND=0, Aab+=1, T1D=2) and the feature value within each size bin. τ > 0 means the feature increases with disease progression; τ < 0 means it decreases. Points colored by significance (red = significant, grey = NS).
- **Bin coverage rule**: The heatmap and the trend plot use the same rule — a bin is tested only if at least 2 donor-status groups each have ≥ 2 donors with finite values after filtering. Bins that fail this rule are silently omitted from both plots so the coverage stays in sync.

### Section 4: Confounders & Deeper Analysis

Full-width **Demographics** card (H5AD only, hidden for Excel fallback):

- **Donor Summary Table**: N islets (primary), N donors, age median/range, % male/female per disease group
- **Age vs Feature** (left): Scatter plot of eligible prepared islet rows after filtering, with overall linear regression line and Pearson correlation. Subtitle notes that islets are correlated within donors.
- **Sex vs Feature** (right): Box plot of feature value by donor status, faceted by sex (Male/Female). Sex-stratified pairwise tests are **BH-corrected across all 6 contrasts** (3 pairs × 2 sexes); the subtitle reports the minimum adjusted q per sex (low-power sex × group cells, n=2-3 donors, are flagged but still corrected jointly).
- **Autoantibody Profile** (Aab+ only, hidden when no AAb data): Per-donor table showing which AAbs are positive, total AAb count, N islets, and mean feature value. Box plot of feature by number of positive autoantibodies (1, 2, 3+) with Kruskal-Wallis test.
- **Covariate-Adjusted Model**: Donor-level linear model adjusting for age and sex, testing whether donor status retains significance.

### Section 5: Methods Reference

Dynamic text (subdued grey background) describing all statistical tests, corrections, and assumptions. Adjusts based on test type and parameters selected in Section 1.

### Using the Statistics Tab

1. Select a feature in the Plot sidebar (e.g., `prop_CD8a Tcell` or `peri_prop_Macrophage`)
2. Switch to the Statistics tab
3. Adjust test parameters in Section 1 (Configure Analysis)
4. Click **Run Statistics** — all sections populate
5. Review the narrative flow: global test → effect sizes → size-dependent patterns → confounders

---

## Spatial Tab

Tissue-wide spatial visualization and peri-islet microenvironment analysis. Combines single-cell tissue scatter plots with Leiden clustering and neighborhood enrichment metrics.

For the current 21-donor schema-4 candidate, the right panel shows 7,493
feature-supported islets from 11,136 classifier-eligible islets. UMAP and Leiden
are saved offline from typed cell counts per canonical cell inside the
classifier-determined endocrine core; disease status only colors the fixed
coordinates. Eligibility requires only that the classifier core does not touch
an image edge and contains at least one canonical cell. Downstream phenotype
support cannot revoke eligibility. The four Leiden resolutions and cluster
composition remain visible without refitting. References below to 15 donors,
5,214 islets and H5AD raw-marker PCA describe the frozen original compatibility
path.

### 3-Panel Layout

1. **Controls Sidebar** (left panel)
   - Donor selector (15 nPOD donors with disease status labels)
   - Color by: Phenotype (21 cell types) or Leiden cluster
   - Leiden resolution dropdown (0.3, 0.5, 0.8, 1.0) — visible when Leiden coloring selected
   - Region filter: All cells, Core + Peri only, or Core only
   - **Color background cells**: Toggle to show tissue background cells in phenotype colors (dimmed) instead of grey
   - **Phenotype filter**: Show/hide individual phenotypes with checkboxes. Use **All** / **None** links to quickly select or clear all types. Dynamically populated from the selected donor's cell types.
   - Donor status checkboxes
   - Phenotype and Donor Status palette selectors
   - Download CSV button

2. **Tissue Scatter** (center, 800px height)
   - Full tissue-wide scatter plot of ~177K cells per donor
   - Background: tissue cells in light grey by default, or colored by phenotype when "Color background cells" is checked
   - Foreground: core/peri cells colored by phenotype or Leiden cluster
   - Spatial coordinate convention: `coord_fixed()` + `scale_y_reverse()` (microscopy standard)
   - **Zoom**: Drag to select a region, then the view zooms in. Double-click or click **Reset Zoom** to restore the full view.
   - Uses ggplot2 `renderPlot` (not plotly) for performance at >100K points

3. **Cluster Composition** (right): stacked bar chart of mean phenotype composition per Leiden cluster, clusters in
   numeric order.
4. **Islet-Level UMAPs** (full-width row below, since 2026-10-19): the Leiden and Donor Status UMAPs side by side at
   one size, both interactive (PNG from each camera button)
   - Built from each object's endocrine composition (the Beta, Alpha, Delta and CHGA+ hormone-negative shares of
     its core endocrine cells) for the objects with ten or more core cells, with no donor correction (2026-10-19).
     Smaller objects have too few cells to place and have no point or cluster. Disease separates mainly by T1D
     beta-cell loss; ND and Aab+ objects overlap, as Aab+ donors range from ND-like to T1D-like.
   - Hidden with message when Leiden data unavailable in H5AD

### Data Sources

- **Tissue CSVs**: `data/donors/{imageid}.csv` (15 files, ~78 MB total). Columns: X/Y centroids (μm), phenotype, cell_region (core/peri/tissue), islet_name.
- **Leiden clustering**: Stored in `islet_explorer.h5ad` .obs as `leiden_0.3`, `leiden_0.5`, `leiden_0.8`, `leiden_1.0` + UMAP coords `leiden_umap_1`, `leiden_umap_2`.
- **Neighborhood metrics**: ~100 columns merged into H5AD .obs (peri-islet composition, per-phenotype enrichment z-scores and distances, plus aggregate immune metrics).

### Neighborhood Analysis Cards

Below the tissue scatter and Leiden panel, three interactive analysis sections visualize per-phenotype peri-islet neighborhood metrics. The cards expose the 12 retained reviewer-facing phenotypes through dynamic choices rather than the previous hand-curated immune-only lists. Cards only appear when neighborhood data is available in the H5AD (not in Excel fallback mode).

**Global Controls** (toolbar above the cards):
- **Endocrine objects** (Data controls; kept in step with the Plot tab's copy): the neighborhood cards use the selected classes. Selecting Islets alone reduces noise from one- and few-cell objects and lowers NA rates in distance metrics.
- **Point size** and **Opacity**: Adjust scatter plot visualization across all cards.
- **Islet count**: Shows how many islets pass the current filters.

**Card A: Peri-Islet Phenotype Enrichment**
- **Left** (`Phenotype (peri enrichment z)` dropdown): Bar chart of donor-level mean/median enrichment z-score across disease stages for the selected phenotype. Shows Kruskal-Wallis p-value. Donor points overlay with summary bars (mean ± SEM or median ± IQR). Outlier exclusion follows the global toggle (Tukey 1.5×IQR per donor status, zero-aware). A **"Show excluded-outlier table"** checkbox next to the phenotype selector reveals the list of excluded islets directly below the bar chart.
- **Right** (`Phenotype (peri vs core)` dropdown): Scatter of peri-zone proportion vs core proportion for the selected phenotype, with dashed y=x diagonal. Sqrt-scale toggle and per-status trend-line option. Points above the line have higher core proportion than peri.
- *Interpretation*: Pick any phenotype (Macrophage, CD8a Tcell, Beta cell, Endothelial, etc.) to see its peri-zone enrichment trajectory across ND→Aab+→T1D.

**Card B: Phenotype Composition & Enrichment**
- **Left**: Grouped bar chart spanning the 12 retained reviewer-facing phenotypes × 3 disease stages. Toggle median/mean summary; `Clip extreme z > 5` is a **display saturation** that caps the visible axis so the chart isn't dominated by a few hugely enriched phenotypes — it is NOT outlier exclusion (see [Outlier Handling](#outlier-handling) for the actual rule). Error bars show IQR (median) or SEM (mean).
- **Right**: Heatmap — phenotypes (columns) × disease stages (rows). Numeric values annotated on each cell.
- **Region toggle**: Switch between:
  - *Peri-islet (enrichment z)*: Poisson z-scores comparing peri-islet vs tissue-wide. Diverging colorscale (blue = depleted, red = enriched). Default.
  - *Core (proportion)*: Raw cell type proportions in islet core. Sequential colorscale (white → red).
  - *Peri-islet (proportion)*: Raw cell type proportions in peri-islet zone. Sequential colorscale.
- *Interpretation*: z > 0 means that phenotype is enriched near islets relative to the tissue-wide proportion. Useful for spotting structural changes (Endothelial, Lymphatic, Neural) alongside immune infiltration.

**Card C: Phenotype Proximity to Islet**
- **Left** (`Distance to nearest` dropdown): Box plot of minimum distance (μm) from islet core centroid to nearest peri-zone cells of the selected phenotype. Includes "Immune (all)" aggregate plus every individual phenotype with `min_dist_*` data. Shows non-NA counts per group. Outlier exclusion follows the global toggle (z > 3 per donor status applied to `log(1 + distance)`); excluded points still render as faint × markers. A **"Show excluded-outlier table"** checkbox next to the metric selector lists the excluded islets (Case ID, distance, z-score) below the box plot. See [Outlier Handling](#outlier-handling).
- **Right** (`Phenotype` dropdown): KDE of signed distance from islet boundary for individual cells of the selected phenotype. Negative = inside islet (core), positive = outside (peri-islet). Dashed line at zero marks the boundary.
- *Interpretation*: NA values in the box plot indicate islets with zero cells of that phenotype in the peri-islet zone (biological, not a data error — common for rare phenotypes like CD4 Tcell, B cell). Selecting Islets alone reduces NAs.

### Data Coverage

- The historical 5,214-row neighborhood artifact covers the trajectory cohort;
  its coverage must not be inferred for every row in the larger prepared summaries.
- Donor 6533 has 191 islets (fully integrated after Parent annotation fix)
- Biological validation: T1D immune_frac_peri > Aab+ > ND (immune infiltration increases with disease progression)

---

## Single-Cell Drill-Down

Click any islet data point in the Plot or Trajectory tab to inspect individual cells.

### How It Works

1. Click a point in the scatter/UMAP plot
2. A segmentation panel appears below the plot with the title showing the islet name and donor info (case ID, disease status, age, sex)

> **macOS note (fixed Jun 2026):** clicking a point previously did nothing on macOS. Drill-down now works on all platforms via a native-click handler. Tip: on the Plot tab, keep **Show individual points** on (the default) so there are per-islet points to click — aggregated summary points are not clickable.

### Panel Layout

- **Left sidebar**: Layer toggles (Single Cells, Peri Boundary, Structures), Color-by dropdown, Phenotype Palette, Show peri-islet cells checkbox. When no cell data is available, shows a boundary legend instead.
- **Center**: Segmentation map with GeoJSON boundaries and optional single-cell overlay
- **Right sidebar**: Cell composition bar chart and cell count table (visible when Single Cells is enabled)

### Single Cells View

- **Color by**: Phenotype (21 categorical colors) or any of 31 protein markers (viridis inferno continuous scale)
- **Show peri-islet**: Toggle checkbox to include/exclude cells from the 20 μm expansion zone
  - Core cells: filled circles
  - Peri-islet cells: open circles (shape 1)
- **Summary panel**: Horizontal bar chart showing phenotype composition of the islet
- **Count table**: Core vs peri cell counts by phenotype

### Data Source

Per-islet cell CSVs in `data/cells/` (5,214 files, ~203 MB total). Each file contains:
- X/Y centroid coordinates (μm, converted to pixel space for GeoJSON overlay)
- Phenotype label (1 of 21 types)
- Region (core or peri-islet)
- Cell/Nucleus area
- 31 protein marker expression values

### Availability

- The source drill-down collection contains 5,214 per-islet CSV files; availability
  for a selected prepared row depends on its identifier mapping.
- If cell data is unavailable for a clicked islet, the panel shows "Boundaries" mode with a message

---

## AI Assistant

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
- **Rebuild app data**: `conda activate scvi-env`, then run the pipeline scripts (`reaggregate_islets.py` → `compute_neighborhood_metrics.py` → `extract_per_islet_cells.py` → `extract_per_donor_tissue.py` → `build_h5ad_for_app.py`). Produces `data/islet_explorer.h5ad`. Details in the admin guide.
- **Local development**: `cd app/shiny_app && Rscript -e 'shiny::runApp(".", port = 7777)'`.
- **App won't load, missing filters, blank panels, or browser-specific issues**: see the admin guide's [Troubleshooting](admin_guide.md#troubleshooting) section.
