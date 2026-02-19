# CLAUDE.md - Islet Explorer

This file provides guidance to Claude Code when working with this repository.

## Project Overview

Interactive Shiny app for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss and immune infiltration across ND -> Aab+ -> T1D progression.

## App Architecture (Modular - Feb 2026)

The Shiny app lives in `app/shiny_app/` and uses a modular architecture:

```
app/shiny_app/
  app.R                        # Thin entrypoint (~455 lines): UI layout, CSS, server wiring
  R/
    00_globals.R               # Package loads, constants, feature flags, paths (h5ad_path, master_path)
    data_loading.R             # load_master_auto(), prep_data(), bin_islet_sizes(), H5AD loading
    utils_stats.R              # summary_stats(), per_bin_anova(), per_bin_kendall()
    utils_safe_join.R          # safe_left_join(), add_islet_key(), compute_diameter_um()
    segmentation_helpers.R     # GeoJSON cache, load_case_geojson(), render_islet_segmentation_plot()
    viewer_helpers.R           # Avivator URL builders, channel config, environment detection
    ai_helpers.R               # OpenAI credential loading, call_openai_chat()
    mod_plot_ui.R / mod_plot_server.R           # Plot tab (scatter, distribution, outliers)
    mod_trajectory_ui.R / mod_trajectory_server.R  # Trajectory tab (UMAP, heatmap, pseudotime)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF, Avivator iframe)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (ANOVA, pairwise, AUC)
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  www/                         # Static assets (images, logos)
  app_original.R               # Backup of pre-modularization monolithic app (5,397 lines)
```

### Data Loading (H5AD with Excel Fallback)

The app uses `load_master_auto()` which tries H5AD first, falls back to Excel:
- `h5ad_path` → `data/islet_explorer.h5ad` (built by `scripts/build_h5ad_for_app.py`)
- `master_path` → `data/master_results.xlsx` (built by `scripts/build_master_excel.py`)

Both produce the same `list(markers, targets, comp, lgals3, phenotypes, donor_demographics)` structure consumed by `prep_data()`. The `phenotypes` and `donor_demographics` elements are extracted from H5AD `.obs` and are `NULL` when loading from Excel (graceful fallback).

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from H5AD or Excel, consumed by Plot, Trajectory, Statistics
- `selected_islet` - reactiveVal, written by Plot and Trajectory click handlers
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer

### Critical: Segmentation Rendering is Root-Level

The `output$islet_segmentation_view` renderPlot is defined in **app.R** at the root level, NOT inside any module. Both the Plot tab and Trajectory tab use embedded panels that reference this shared output using non-namespaced `plotOutput("islet_segmentation_view")`. This is required for proper output binding since only one tab is visible at a time.

### Plotly Click Events in Modules

- Plot module: `source = ns("plot_scatter")` with `event_register("plotly_click")`
- Trajectory module: `source = ns("traj_scatter")` with `event_register("plotly_click")`
- Both use namespaced sources; `event_data()` calls must use matching `ns()` prefix

## Key Data Files

- `data/islet_explorer.h5ad` - **Primary app data** (47 MB, 1,015 islets, groovy + trajectory + donor metadata)
- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets) — Excel fallback
- `data/adata_ins_root.h5ad` - Trajectory h5ad (1,015 islets, 31 vars, DPT pseudotime, UMAP, rebuilt 2026-02-17)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/DATA_PROVENANCE.md` - Documents canonical H5AD lineage and data sources

## Data Pipeline (Phases 2-3, executed Feb 2026)

Canonical H5AD lineage:
```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells)
  ↓ islet_analysis/fixed_islet_aggregation.py (core only, min_cells=20)
islet_analysis/islets_core_fixed.h5ad  (1,015 islets, proteins in .X, scVI means in .obsm)
  ↓ notebooks/rebuild_trajectory.ipynb (scVI neighbors n=15 cosine, PAGA→UMAP, DPT)
data/adata_ins_root.h5ad  (1,015 islets + pseudotime + UMAP)
  ↓ scripts/build_h5ad_for_app.py (+ groovy TSV data + donor metadata from h5ad obs)
data/islet_explorer.h5ad  (all app data in one file, 47 MB)
```

### Pipeline Scripts & Notebooks

- `notebooks/scvi_qc_validation.ipynb` - Validates scVI batch correction (silhouette, LISI, UMAP)
- `notebooks/rebuild_trajectory.ipynb` - Regenerates trajectory from fixed pipeline
- `scripts/build_h5ad_for_app.py` - Builds enriched H5AD for app from trajectory + groovy exports
- `scripts/build_master_excel.py` - Builds master_results.xlsx from groovy TSV exports
- `islet_analysis/fixed_islet_aggregation.py` - Core aggregation module (islet-level from single-cell)

### Python Environment

Pipeline notebooks use the `scvi-env` conda environment:
```bash
conda activate scvi-env  # scanpy, anndata, scvi-tools, sklearn, scib-metrics
```

## Interactive Features (Phase 5, Feb 2026)

### Phenotype Composition Explorer
The Plot tab Composition mode exposes 21 cell-type proportions (`prop_*` columns from H5AD `.obs`) alongside the original 3 hormone fractions. The selector uses grouped `selectInput` ("Hormone Fractions" / "Cell Type Proportions"). Phenotype values are already 0-1 in `.obs`, scaled to percentage for display. Falls back to hormone-only when loading from Excel.

### Age & Gender Demographic Filters
Age slider and gender checkboxes appear in the Plot sidebar below Donor Status when demographics are available (H5AD path only). Filters apply in `raw_df_base()` after diameter range filtering, before AAb filtering. Both `renderUI` outputs return `NULL` when demographics columns are absent (Excel fallback).

### Multi-Feature Trajectory Heatmap
Below the donor-status heatmap in the Trajectory tab, a multi-row z-scored expression heatmap shows binned marker expression along pseudotime. Users select markers via checkboxes (8 defaults: INS, GCG, SST, CD3e, CD8a, CD68, CD45, HLADR). Z-scores are clamped to [-2.5, 2.5] with blue-white-red diverging colormap. Dynamic height: `max(150, 40 + n_markers * 30)`. Data comes from `tr$adata` (trajectory H5AD), independent of main data path.

### Segmentation Viewer Pattern
Both Plot and Trajectory tabs use **embedded panels** (not modals) for segmentation display. Click a point → `selected_islet()` updates → `output$segmentation_viewer_panel` renders inline with the non-namespaced `plotOutput("islet_segmentation_view")`. Close button sets `selected_islet(NULL)`.

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- UMAP "Selected Feature" uses continuous viridis (inferno) colormap scaled to data min/max
- Donor status colors: ND = green (#2ca02c), Aab+ = yellow (#ffcc00), T1D = purple (#9467bd)
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 19 phenotypes, 18 markers
- **Donor metadata**: Comes from `islets_core_fixed.h5ad` obs, NOT from `CODEX_Pancreas_Donors.xlsx` (different cohort)
- **H5AD obs index**: `islets_core_fixed.h5ad` index name is `islet_id` — same as column, use `reset_index(drop=True)`

## Running the App

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Analysis Pipeline

See `islet_analysis/CLAUDE.md` for the Python analysis pipeline (phenotyping, aggregation, trajectory).
See `data/DATA_PROVENANCE.md` for full data lineage documentation.

## Skills Registry

The `Skills_Registry/` submodule contains cross-session development knowledge. Custom slash commands:
- `/advise` — search registry for relevant experiments before starting new work
- `/retrospective` — save session learnings as a new skill

Skills are registered in `.claude/skills/<name>/SKILL.md` (NOT in Skills_Registry/CLAUDE.md).

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate
