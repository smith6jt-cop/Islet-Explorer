# Islet Explorer — User Guide

Interactive web application for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss, immune infiltration, and cellular microenvironment changes across the ND → Aab+ → T1D disease progression.

**Production URL**: `http://10.15.152.7:8080/islet-explorer/`

---

## Table of Contents

1. [Overview](#overview)
2. [Plot Tab](#plot-tab)
3. [Trajectory Tab](#trajectory-tab)
4. [Viewer Tab](#viewer-tab)
5. [Statistics Tab](#statistics-tab)
6. [Spatial Tab](#spatial-tab)
7. [Single-Cell Drill-Down](#single-cell-drill-down)
8. [Data Pipeline](#data-pipeline)
9. [Deployment & Administration](#deployment--administration)
10. [Troubleshooting](#troubleshooting)

---

## Overview

### Data Source

The app loads from a single enriched H5AD file (`data/islet_explorer.h5ad`, ~48 MB) containing:

- **1,015 islets** from 15 nPOD donors (ND, Aab+, T1D)
- **31 protein markers** (mean expression per islet)
- **Pseudotime trajectory** (DPT with INS root, scVI-corrected)
- **21 phenotype proportions** (from single-cell phenotyping)
- **Donor demographics** (age, gender, autoantibody status)
- **62 neighborhood metrics** (peri-islet composition, immune infiltration, enrichment z-scores, distances)
- **QuPath groovy data** (targets, markers, composition, LGALS3)

If the H5AD is unavailable, the app falls back to `data/master_results.xlsx` with reduced functionality (no phenotype proportions, no demographics, no neighborhood metrics).

### Disease Groups

| Group | Color | Description |
|-------|-------|-------------|
| **ND** | Green (#2ca02c) | Non-diabetic control donors |
| **Aab+** | Yellow (#ffcc00) | Autoantibody-positive (pre-T1D) |
| **T1D** | Purple (#9467bd) | Type 1 diabetic donors |

---

## Plot Tab

The Plot tab is the primary exploration interface. It has a **sidebar** (left) and **main panel** (right).

### Sidebar Controls

- **Mode**: Switch between `Scatter` (feature vs diameter) and `Composition` (cell type proportions)
- **Feature/Composition selector**: Choose what to plot
  - *Scatter mode*: Select from markers (31 proteins) or targets (hormone intensities)
  - *Composition mode*: Four option groups:
    - **Hormone Fractions** — INS, GCG, SST fractions
    - **Cell Type Proportions** — 21 phenotype proportions from single-cell data
    - **Peri-Islet Proportions** — Proportions of each phenotype in the 20 μm expansion zone around the islet
    - **Immune Metrics** — `immune_frac_peri`, `immune_frac_core`, `immune_ratio`, `cd8_to_macro_ratio`, `tcell_density_peri`
- **Region**: Filter by tissue region (core only, peri-islet, or all)
- **Donor Status**: Select disease groups to include (ND, Aab+, T1D)
- **AAb Filter**: Filter by autoantibody positivity (when available)
- **Age Range**: Slider to restrict donor age (H5AD only)
- **Gender**: Checkbox filter (H5AD only)
- **Diameter Range**: Restrict islet diameter in μm

### Main Panel

- **Scatter/Distribution Plot**: Interactive Plotly scatter (diameter vs feature value) or distribution comparison (violin/box by disease stage)
- **Outlier Table**: Lists extreme values for the selected feature
- **Click Interaction**: Click any data point to open the segmentation viewer (see [Single-Cell Drill-Down](#single-cell-drill-down))

### Tips

- Use composition mode with "Peri-Islet Proportions" to visualize the cellular microenvironment around each islet
- The "Immune Metrics" group provides pre-computed ratios that highlight immune infiltration differences between disease stages
- Log-scale features (e.g., cell counts) use pseudo-log transformation so zero values remain visible at y=0

---

## Trajectory Tab

Visualizes the disease progression trajectory computed via diffusion pseudotime (DPT) with INS as the root anchor.

### Layout

The Trajectory tab uses a full-width layout (no sidebar). All controls are inline.

### Components

1. **UMAP Scatter** — Interactive Plotly UMAP colored by:
   - Donor status (default)
   - Selected feature (viridis inferno colormap, scaled to data min/max)
   - Pseudotime (DPT)

2. **Donor Status Heatmap** — Mean marker expression across pseudotime bins, split by disease group. Shows how protein expression changes along the trajectory.

3. **Multi-Feature Heatmap** — Z-scored expression of user-selected markers along pseudotime bins. Default markers: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR. Blue-white-red diverging colormap, z-scores clamped to [-2.5, 2.5].

4. **Click Interaction** — Click any UMAP point to:
   - Highlight the selected islet
   - Open the segmentation viewer panel
   - Jump to the OME-TIFF Viewer tab for that donor

### Key Biological Insights

- **INS vs pseudotime**: Strong negative correlation (r = -0.741) — insulin decreases along the trajectory
- **Disease ordering**: ND (early) → Aab+ (middle) → T1D (late) along pseudotime
- **GCG increase**: Glucagon expression increases as insulin decreases, reflecting alpha cell persistence

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

---

## Statistics Tab

Rigorous statistical analysis of the currently selected feature from the Plot sidebar. The Statistics tab **shares the Plot sidebar** — changing the feature, filters, or disease groups in the sidebar updates both Plot and Statistics simultaneously.

### 7-Card Layout

1. **Overview Banner**
   - N islets included, global p-value, effect size (η²)
   - Inline controls: test type (parametric/non-parametric), significance level (α), outlier removal, bin width, diameter range

2. **Hypothesis Testing**
   - Global test: ANOVA (parametric) or Kruskal-Wallis (non-parametric)
   - Pairwise comparisons table with Cohen's d effect sizes and 95% CIs
   - Forest plot visualization of pairwise differences

3. **Per-Bin Significance Heatmap**
   - Plotly heatmap: x = diameter bins (midpoints), y = test type, z = -log₁₀(p)
   - Reveals which size ranges drive significance

4. **Trend Analysis**
   - Kendall τ per diameter bin
   - Line plot with significance coloring (significant bins highlighted)

5. **Demographics** (H5AD only, hidden for Excel fallback)
   - Age scatter with linear regression
   - Gender-stratified tests
   - Covariate-adjusted model

6. **AUC Analysis**
   - Trapezoidal AUC by disease group (pseudo-log scale handles zeros)
   - Percentage change relative to ND baseline
   - Guards against division by zero when ND AUC = 0

7. **Methods & Interpretation**
   - Dynamic text describing all statistical tests, corrections, and assumptions
   - Adjusts based on the controls selected in card 1

### Using the Statistics Tab

1. Select a feature in the Plot sidebar (e.g., `immune_frac_peri`)
2. Switch to the Statistics tab — analysis runs automatically
3. Adjust test parameters in the Overview banner (card 1)
4. Pairwise results and heatmap update reactively

---

## Spatial Tab

Dedicated analysis of the peri-islet cellular microenvironment (the 20 μm expansion zone around each islet).

### 7-Card Layout

1. **Overview & Controls**
   - Metric category selector (Peri-Islet Composition, Immune Infiltration, Enrichment, Distance)
   - Specific metric dropdown (e.g., `peri_prop_Macrophage`, `immune_frac_peri`)
   - Donor status filter and diameter range

2. **Feature Distribution**
   - Violin + box plot by disease stage (ND/Aab+/T1D)
   - Shows how the selected neighborhood metric varies across disease groups

3. **Core vs Peri Comparison**
   - Paired comparison of the same cell type in core vs peri-islet zones
   - Highlights differential enrichment patterns

4. **Peri-Islet Phenotype Heatmap**
   - Mean proportion of each of 21 phenotypes in the peri-islet zone
   - Columns = disease stages, rows = phenotypes
   - Z-scored for cross-phenotype comparison

5. **Immune Enrichment**
   - Bar chart of Poisson enrichment z-scores by disease stage
   - Positive z = enriched relative to tissue-wide baseline
   - 7 immune cell types: CD8a T cell, Macrophage, CD4 T cell, NK cell, etc.

6. **Pseudotime Correlation**
   - Scatter: neighborhood metric vs DPT pseudotime
   - Reveals how the microenvironment changes along disease progression

7. **Statistical Tests & Methods**
   - Reuses the same test utilities as the Statistics tab (Cohen's d, η², pairwise Wilcoxon)

### Data Coverage

- 949 of 1,015 islets (93.5%) have peri-islet data
- 66 islets lack peri-islet cells → shown as N/A in plots
- Biological validation: T1D immune_frac_peri (0.155) > Aab+ (0.106) > ND (0.069)

---

## Single-Cell Drill-Down

Click any islet data point in the Plot or Trajectory tab to inspect individual cells.

### How It Works

1. Click a point in the scatter/UMAP plot
2. A segmentation panel appears below the plot with two view modes:
   - **Boundaries** — GeoJSON islet boundary overlay (default if no cell data)
   - **Single Cells** — Individual cells plotted in spatial context over the islet boundary

### Single Cells View

- **Color by**: Phenotype (21 categorical colors) or any of 31 protein markers (viridis inferno continuous scale)
- **Show peri-islet**: Toggle checkbox to include/exclude cells from the 20 μm expansion zone
  - Core cells: filled circles
  - Peri-islet cells: open circles (shape 1)
- **Summary panel**: Horizontal bar chart showing phenotype composition of the islet
- **Count table**: Core vs peri cell counts by phenotype

### Data Source

Per-islet cell CSVs in `data/cells/` (949 files, ~111 MB total). Each file contains:
- X/Y centroid coordinates (μm, converted to pixel space for GeoJSON overlay)
- Phenotype label (1 of 21 types)
- Region (core or peri-islet)
- Cell/Nucleus area
- 31 protein marker expression values

### Availability

- 949 of 1,015 islets have cell data
- If cell data is unavailable for a clicked islet, the panel shows "Boundaries" mode with a message

---

## Data Pipeline

### Rebuilding the App Data

The canonical pipeline rebuilds `data/islet_explorer.h5ad` from source:

```bash
# Activate the Python environment
conda activate scvi-env

# Step 1: Compute neighborhood metrics (from 2.6M-cell single-cell H5AD)
python scripts/compute_neighborhood_metrics.py
# → data/neighborhood_metrics.csv (1,015 rows × 62 cols)

# Step 2: Extract per-islet cell CSVs (for drill-down viewer)
python scripts/extract_per_islet_cells.py
# → data/cells/*.csv (949 files, ~111 MB)

# Step 3: Build enriched H5AD (trajectory + groovy + neighborhood + donor metadata)
python scripts/build_h5ad_for_app.py
# → data/islet_explorer.h5ad (~48 MB)
```

### Full Lineage

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, canonical single-cell)
  ↓ islet_analysis/fixed_islet_aggregation.py (core only, min_cells=20)
islets_core_fixed.h5ad  (1,015 islets, proteins + scVI embeddings)
  ↓ notebooks/rebuild_trajectory.ipynb (scVI neighbors, PAGA→UMAP, DPT)
data/adata_ins_root.h5ad  (+ pseudotime + UMAP)
  ↓ scripts/build_h5ad_for_app.py (+ groovy + neighborhood + donor metadata)
data/islet_explorer.h5ad  (complete app data)
```

Branch pipelines:
- `scripts/compute_neighborhood_metrics.py` — peri-islet composition, immune metrics, enrichment z-scores, distances
- `scripts/extract_per_islet_cells.py` — individual cell CSVs for drill-down

### Upstream Data Sources

| Source | Location | Description |
|--------|----------|-------------|
| Single-cell H5AD | `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` | 2.6M cells, scVI batch-corrected, phenotyped |
| Groovy TSV exports | `panc_CODEX/results/groovy_exports/` | QuPath islet measurements (15 donors × 4 types) |
| GeoJSON boundaries | `data/json/*.geojson` | Islet segmentation polygons |
| Spatial lookup | `data/islet_spatial_lookup.csv` | Islet centroid coordinates |

See `data/DATA_PROVENANCE.md` for the complete lineage documentation including validation results.

---

## Deployment & Administration

### Architecture

```
nginx (port 8080)
  ↓ reverse-proxy /islet-explorer/
shiny-server (port 3838)
  ↓ symlink /srv/shiny-server/islet-explorer → app/shiny_app/
Shiny app (R worker process)
```

### Updating the App

Code changes take effect when shiny-server spawns a new R worker. No restart is needed for code-only changes.

If a stale R worker persists with old code:
```bash
# Find and kill stale Shiny R workers
ps aux | grep R | grep shiny
kill <PID>
```

Note: `sudo systemctl restart shiny-server` is not available — use the worker kill approach above.

### Local Development

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
# Then open http://localhost:7777 in browser
```

Port 7777 is for local testing only and is not accessible to end users.

### R Dependencies

```
shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite,
RColorBrewer, scales, anndata, reticulate
```

### Python Dependencies (pipeline only)

```bash
conda activate scvi-env
# scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, pandas, numpy
```

---

## Troubleshooting

### App fails to load

1. Check that `data/islet_explorer.h5ad` exists and is ~48 MB
2. Check that R packages are installed: `library(anndata); library(reticulate)`
3. Check that Python is accessible via reticulate: `reticulate::py_available()`
4. Fall back to Excel: ensure `data/master_results.xlsx` exists (reduced functionality)

### No phenotype or demographic filters

The app is running from Excel fallback. Rebuild the H5AD:
```bash
conda activate scvi-env
python scripts/build_h5ad_for_app.py
```

### No neighborhood metrics in Plot selector

1. Verify `data/neighborhood_metrics.csv` exists (1,015 rows)
2. Rebuild the H5AD: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh the browser

### Single-cell drill-down not available

1. Verify `data/cells/` directory contains ~949 CSV files
2. If missing, regenerate: `python scripts/extract_per_islet_cells.py`

### Segmentation viewer shows "No GeoJSON found"

1. Check that `data/json/` or `data/gson/` has `.geojson` / `.geojson.gz` files
2. Case ID zero-padding: GeoJSON files use 4-digit padded IDs (`0112.geojson`), data uses unpadded (`112`). The app handles this automatically via `sprintf("%04d", ...)` fallback.

### Statistics tab shows no results

The Statistics tab uses data from the Plot sidebar. Ensure:
1. A feature is selected in the Plot sidebar
2. At least 2 disease groups are checked
3. The diameter range includes some islets

### Spatial tab is empty

1. Verify the H5AD has neighborhood columns: look for `peri_prop_*` in the data
2. If absent, run `scripts/compute_neighborhood_metrics.py` and rebuild the H5AD
