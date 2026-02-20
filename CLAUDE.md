# CLAUDE.md - Islet Explorer

This file provides guidance to Claude Code when working with this repository.

## Project Overview

Interactive Shiny app for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss and immune infiltration across ND -> Aab+ -> T1D progression. Features islet-level aggregated views, pseudotime trajectory, spatial neighborhood analysis, and single-cell drill-down.

## App Architecture (Modular - Feb 2026)

The Shiny app lives in `app/shiny_app/` and uses a modular architecture:

```
app/shiny_app/
  app.R                        # Entrypoint (~510 lines): UI layout, CSS, server wiring, root-level outputs
  R/
    00_globals.R               # Package loads, constants, feature flags, paths (h5ad_path, master_path)
    data_loading.R             # load_master_auto(), prep_data(), bin_islet_sizes(), H5AD loading
    utils_stats.R              # summary_stats(), per_bin_anova(), per_bin_kendall(), cohens_d(), eta_squared(), pairwise_wilcox()
    utils_safe_join.R          # safe_left_join(), add_islet_key(), compute_diameter_um()
    segmentation_helpers.R     # GeoJSON cache, build_segmentation_base_plot(), render_islet_segmentation_plot()
    drilldown_helpers.R        # PHENOTYPE_COLORS, load_islet_cells(), render_islet_drilldown_plot(), render_drilldown_summary()
    viewer_helpers.R           # Avivator URL builders, channel config, environment detection
    ai_helpers.R               # OpenAI credential loading, call_openai_chat()
    mod_plot_ui.R / mod_plot_server.R           # Plot tab (scatter, distribution, outliers, segmentation+drilldown panel)
    mod_trajectory_ui.R / mod_trajectory_server.R  # Trajectory tab (UMAP, heatmap, pseudotime, segmentation+drilldown panel)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF, Avivator iframe)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (7-card: hypothesis, heatmap, trend, demographics, AUC)
    mod_spatial_ui.R / mod_spatial_server.R     # Spatial tab (7-card: neighborhood metrics, enrichment, pseudotime correlation)
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  www/                         # Static assets (images, logos)
  app_original.R               # Backup of pre-modularization monolithic app (5,397 lines)
```

### Data Loading (H5AD with Excel Fallback)

The app uses `load_master_auto()` which tries H5AD first, falls back to Excel:
- `h5ad_path` -> `data/islet_explorer.h5ad` (built by `scripts/build_h5ad_for_app.py`)
- `master_path` -> `data/master_results.xlsx` (built by `scripts/build_master_excel.py`)

Both produce the same `list(markers, targets, comp, lgals3, phenotypes, donor_demographics, neighborhood)` structure consumed by `prep_data()`. The `phenotypes`, `donor_demographics`, and `neighborhood` elements are extracted from H5AD `.obs` and are `NULL` when loading from Excel (graceful fallback).

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from H5AD or Excel, consumed by Plot, Trajectory, Statistics, Spatial
- `selected_islet` - reactiveVal, written by Plot and Trajectory click handlers
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer
- `active_tab` - `reactive(input$tabs)`, passed to Plot and Trajectory to prevent duplicate output IDs

### Critical: Root-Level Outputs + Active Tab Guard

The following render outputs are defined in **app.R** at the root level, NOT inside any module:
- `output$islet_segmentation_view` - GeoJSON boundary plot (shared by Plot + Trajectory)
- `output$islet_drilldown_view` - Single-cell spatial overlay plot
- `output$islet_drilldown_summary` - Cell composition bar chart
- `output$islet_drilldown_table` - Core/peri cell counts

Both Plot and Trajectory tabs embed these using non-namespaced `plotOutput("islet_segmentation_view")` etc. **IMPORTANT**: Each module's `segmentation_viewer_panel` renderUI must guard with `if (active_tab() != "Plot") return(NULL)` (or `"Trajectory"` respectively). Without this guard, both renderUIs fire when `selected_islet()` changes, creating duplicate DOM IDs — Shiny can only bind one, causing the other to grey out. Shiny `renderUI` fires based on reactive dependencies, NOT tab visibility.

### Non-Namespaced Inputs for Root-Level Outputs

The segmentation panel `renderUI` in both modules generates non-namespaced inputs:
- `drilldown_view_mode` - "Boundaries" or "Single Cells" toggle
- `drilldown_color_by` - "phenotype" or a marker column name
- `drilldown_show_peri` - Whether to show peri-islet cells

These are read by the root-level `renderPlot` outputs in `app.R`.

### Plotly Click Events in Modules

- Plot module: `source = ns("plot_scatter")` with `event_register("plotly_click")`
- Trajectory module: `source = ns("traj_scatter")` with `event_register("plotly_click")`
- Both use namespaced sources; `event_data()` calls must use matching `ns()` prefix

## Key Data Files

- `data/islet_explorer.h5ad` - **Primary app data** (~48 MB, 1,015 islets, groovy + trajectory + donor + neighborhood metrics)
- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets) -- Excel fallback
- `data/adata_ins_root.h5ad` - Trajectory h5ad (1,015 islets, 31 vars, DPT pseudotime, UMAP)
- `data/neighborhood_metrics.csv` - Per-islet peri-islet metrics (1,015 rows, 62 columns)
- `data/cells/*.csv` - Per-islet single-cell CSVs for drill-down (~949 files, ~111 MB total)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/DATA_PROVENANCE.md` - Documents canonical H5AD lineage and data sources

## Data Pipeline

### Canonical H5AD Lineage

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  |
  |-- islet_analysis/fixed_islet_aggregation.py (core only, min_cells=20)
  |     -> islet_analysis/islets_core_fixed.h5ad  (1,015 islets)
  |         |
  |         |-- notebooks/rebuild_trajectory.ipynb (scVI n=15 cosine, PAGA->UMAP, DPT)
  |         |     -> data/adata_ins_root.h5ad  (+ pseudotime + UMAP)
  |         |
  |         |-- scripts/compute_neighborhood_metrics.py
  |         |     -> data/neighborhood_metrics.csv  (62 cols: peri-islet composition, immune, enrichment, distance)
  |         |
  |         +-- scripts/extract_per_islet_cells.py
  |               -> data/cells/*.csv  (949 files, 37 cols: coords, phenotype, region, 31 markers)
  |
  +-- scripts/build_h5ad_for_app.py (trajectory + groovy + donor + neighborhood)
        -> data/islet_explorer.h5ad  (all app data in one file)
```

### Pipeline Scripts & Notebooks

- `scripts/compute_neighborhood_metrics.py` - Computes per-islet peri-islet metrics from single-cell H5AD
- `scripts/extract_per_islet_cells.py` - Extracts per-islet cell CSVs for drill-down viewer
- `scripts/build_h5ad_for_app.py` - Builds enriched H5AD (trajectory + groovy + neighborhood)
- `scripts/build_master_excel.py` - Builds master_results.xlsx from groovy TSV exports
- `notebooks/scvi_qc_validation.ipynb` - Validates scVI batch correction (silhouette, LISI, UMAP)
- `notebooks/rebuild_trajectory.ipynb` - Regenerates trajectory from fixed pipeline
- `islet_analysis/fixed_islet_aggregation.py` - Core aggregation module (islet-level from single-cell)

### Python Environment

Pipeline scripts use the `scvi-env` conda environment:
```bash
conda activate scvi-env  # scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy
```

## Spatial Neighborhood Analysis (Phase 7, Feb 2026)

### Peri-Islet Metrics (4 categories, 62 columns)

| Category | Columns | Description |
|----------|---------|-------------|
| Peri-islet composition | `peri_prop_*`, `peri_count_*`, `total_cells_peri` | Proportion & count of each of 21 phenotypes in 20um expansion zone |
| Immune infiltration | `immune_frac_peri/core`, `immune_ratio`, `cd8_to_macro_ratio`, `tcell_density_peri` | Immune fractions, ratios, density |
| Enrichment z-scores | `enrich_z_CD8a_Tcell`, `enrich_z_Macrophage`, etc. | Poisson z comparing peri-islet vs tissue-wide proportion |
| Distance metrics | `min_dist_cd8`, `min_dist_Macrophage`, `min_dist_immune_mean` | Min distance from islet centroid to nearest immune cells |

### Data Flow

`compute_neighborhood_metrics.py` reads single-cell H5AD -> `neighborhood_metrics.csv` (1,015 rows) -> merged into `islet_explorer.h5ad` `.obs` by `build_h5ad_for_app.py` (step 4.5) -> extracted in `data_loading.R` `load_master_h5ad()` -> merged into `comp` by `prep_data()` -> available in Plot composition selector (4 option groups) + Statistics tab + Spatial tab.

### Key Details

- 949/1,015 islets have peri-islet data; 66 get NaN (no `_exp20um` cells in single-cell H5AD)
- Column naming sanitizes phenotype names: spaces -> `_`, `+` -> `plus` (e.g., `peri_prop_ECADplus`)
- Immune signal validated: T1D immune_frac_peri (0.155) > Aab+ (0.106) > ND (0.069)
- Plot composition selector now has 4 option groups: Hormone Fractions, Cell Type Proportions, Peri-Islet Proportions, Immune Metrics

### Spatial Tab (mod_spatial_ui/server.R)

Dedicated 7-card layout with inline controls:
1. **Overview + Controls** -- metric category, feature, donor status, diameter range
2. **Feature Distribution** -- violin/box by disease stage (ND/Aab+/T1D)
3. **Core vs Peri Comparison** -- immune fraction paired comparison
4. **Peri-Islet Phenotype Heatmap** -- mean proportion x 21 phenotypes x 3 stages
5. **Immune Enrichment** -- grouped bar chart of enrichment z-scores by stage
6. **Pseudotime Correlation** -- scatter with Kendall tau
7. **Statistical Tests** -- Kruskal-Wallis + pairwise Wilcoxon with Cohen's d

Wired as `spatial_server("spatial", prepared)` -- no sidebar sharing, inline controls only.

## Single-Cell Drill-Down (Phase 8, Feb 2026)

### Per-Islet Cell Data

`extract_per_islet_cells.py` reads single-cell H5AD -> outputs `data/cells/{imageid}_Islet_{N}.csv`:
- 949 files (matching islets with cell data), ~111 MB total
- 37 columns: `X_centroid`, `Y_centroid`, `phenotype`, `cell_region` (core/peri), `Cell Area`, `Nucleus Area`, + 31 protein markers
- File naming matches `combined_islet_id` from `islet_spatial_lookup.csv`

### Segmentation Panel Extension

When cell data exists for a clicked islet, the segmentation panel shows:
- **Boundaries | Single Cells** radio toggle (defaults to "Single Cells")
- **Color-by dropdown**: phenotype (categorical) or any of 31 markers (viridis inferno continuous)
- **Peri-islet toggle**: checkbox to show/hide peri-islet cells
- **Layout**: 8-col plot + 4-col sidebar (composition chart + cell count table)
- Falls back to Boundaries-only when cell data is unavailable for that islet

### Coordinate Alignment

Cell centroids (`X_centroid`, `Y_centroid`) are in micrometers. GeoJSON polygons are in pixels. Conversion: `x_px = X_centroid / PIXEL_SIZE_UM`. This matches the existing `render_islet_segmentation_plot()` coordinate system.

### drilldown_helpers.R

- `PHENOTYPE_COLORS` -- 21 phenotypes with distinctive hex palette (endocrine=warm, immune=cool, structural=neutral)
- `drilldown_available()` -- checks if `data/cells/` exists with CSVs
- `load_islet_cells(imageid, islet_key)` -- CSV loader with `new.env()` caching
- `render_islet_drilldown_plot(info, cells, color_by, show_peri)` -- cells over `build_segmentation_base_plot()`
- `render_drilldown_summary(cells)` -- horizontal bar chart of phenotype composition

### segmentation_helpers.R Refactor

`build_segmentation_base_plot(info)` extracted from `render_islet_segmentation_plot()`. Returns ggplot with GeoJSON polygon layers + coord_sf (WITHOUT crosshairs or title). Reused by both `render_islet_segmentation_plot()` (adds crosshairs + title) and `render_islet_drilldown_plot()` (adds cell scatter layer).

## Statistics Tab (Phase 6, Feb 2026)

### Shared Sidebar Architecture
The Plot sidebar (mode, feature, region, donor status, AAb, age, gender filters) is visible on both the Plot and Statistics tabs. Achieved via `conditionalPanel` condition `"input.tabs == 'Plot' || input.tabs == 'Statistics'"` and matching JS `adjustLayout()` logic in `app.R`. The Statistics module consumes `plot_returns$raw_df` and `plot_returns$summary_df` directly -- no data duplication.

### 7-Card Layout
1. **Overview Banner** -- N islets, global p-value, effect size eta-squared, inline controls
2. **Hypothesis Testing** -- Global test (ANOVA or Kruskal-Wallis) + pairwise table with Cohen's d + forest plot
3. **Per-Bin Significance Heatmap** -- `plot_ly(type="heatmap")` with x=bin midpoints, y=test type, z=-log10(p)
4. **Trend Analysis** -- Kendall tau per bin, line plot with significance coloring
5. **Demographics** -- Conditional on H5AD (hidden for Excel fallback). Age scatter, gender-stratified tests
6. **AUC Analysis** -- Trapezoidal AUC by donor group with percentage change interpretation
7. **Methods & Interpretation** -- Dynamic text describing all tests, corrections, assumptions

### Effect Size Utilities (`utils_stats.R`)
- `cohens_d(x, y)` -- returns `list(d, ci_lo, ci_hi)` (NOT `ci_lower`/`ci_upper`)
- `eta_squared(fit)` -- from `anova()` output
- `pairwise_wilcox(df, group_col, value_col)` -- returns data.frame with columns `group1`, `group2`, `p_value` (NOT `p.adj`/`statistic`)

## Interactive Features (Phase 5, Feb 2026)

### Phenotype Composition Explorer
Plot tab Composition mode exposes 21 cell-type proportions (`prop_*` from H5AD `.obs`) + 21 peri-islet proportions (`peri_prop_*`) + 5 immune metrics alongside 3 hormone fractions. Grouped `selectInput` with 4 option groups.

### Age & Gender Demographic Filters
Age slider and gender checkboxes in Plot sidebar (H5AD only). `renderUI` returns `NULL` when absent.

### Multi-Feature Trajectory Heatmap
Z-scored expression along pseudotime. Clamped [-2.5, 2.5], blue-white-red, dynamic height.

### Segmentation + Drill-Down Pattern
Click a point -> `selected_islet()` updates -> embedded panel renders inline with Boundaries/Single Cells toggle. Close button sets `selected_islet(NULL)`.

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- UMAP "Selected Feature" uses continuous viridis (inferno) colormap scaled to data min/max
- Donor status colors: ND = green (#2ca02c), Aab+ = yellow (#ffcc00), T1D = purple (#9467bd)
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 19 phenotypes, 18 markers
- **Donor metadata**: Comes from `islets_core_fixed.h5ad` obs, NOT from `CODEX_Pancreas_Donors.xlsx` (different cohort)
- **H5AD obs index**: `islets_core_fixed.h5ad` index name is `islet_id` -- same as column, use `reset_index(drop=True)`
- **Case ID zero-padding**: GeoJSON files use `0112.geojson` (4-digit padded), data uses `112` (unpadded). `load_case_geojson()` and click handlers use `sprintf("%04d", ...)` fallback
- **Log-scale with zeros**: Use `scales::pseudo_log_trans(base=10)` instead of `scale_y_log10()` -- zeros map to 0 (visible) instead of -Infinity (dropped)
- **Plot defaults**: Point size = 3.0, transparency = 0.6
- **Peri-islet data guard**: Always check `total_cells_peri > 0` before using peri metrics (66 islets have 0 peri cells)
- **Column name sanitization**: Phenotype names use `_` for spaces, `plus` for `+` in peri-islet columns
- **Single-cell Parent column**: `Islet_N` = core cells, `Islet_N_exp20um` = peri-islet cells
- **CSS overflow for cards with dropdowns**: `selectInput` menus extend below their container. Add `overflow: visible;` to card styles or dropdowns get clipped behind cards.
- **Font size minimums**: Use `h5` (not `h6`) with explicit `font-size: 15px` for panel headings. Legend items minimum 14-15px.

## Deployment Architecture

- **Production URL**: `http://10.15.152.7:8080/islet-explorer/`
- nginx (port 8080) reverse-proxies `/islet-explorer/` -> shiny-server (port 3838)
- shiny-server serves from symlink: `/srv/shiny-server/islet-explorer` -> `app/shiny_app/`
- Changes take effect when shiny-server spawns a fresh R worker (no restart needed; kill stale workers if needed)
- Dev server (`Rscript -e 'shiny::runApp(".", port=7777)'`) is for local testing only

## Running the App

Production (always use this):
```
http://10.15.152.7:8080/islet-explorer/
```

Development only:
```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Running Pipeline Scripts

```bash
conda activate scvi-env

# Phase 7: Compute neighborhood metrics (reads 2.6M-cell H5AD, ~3 min)
python scripts/compute_neighborhood_metrics.py

# Phase 8: Extract per-islet cell CSVs (reads 2.6M-cell H5AD, ~5 min)
python scripts/extract_per_islet_cells.py

# Rebuild enriched H5AD (merges groovy + trajectory + neighborhood)
python scripts/build_h5ad_for_app.py
```

## Analysis Pipeline

See `islet_analysis/CLAUDE.md` for the Python analysis pipeline (phenotyping, aggregation, trajectory).
See `data/DATA_PROVENANCE.md` for full data lineage documentation.
See `docs/user_guide.md` for end-user documentation.

## Skills Registry

The `Skills_Registry/` submodule contains cross-session development knowledge. Custom slash commands:
- `/advise` -- search registry for relevant experiments before starting new work
- `/retrospective` -- save session learnings as a new skill

Skills are registered in `.claude/skills/<name>/SKILL.md` (NOT in Skills_Registry/CLAUDE.md).

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate
