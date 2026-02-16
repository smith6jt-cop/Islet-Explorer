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
    00_globals.R               # Package loads, constants, feature flags (auto-sourced first)
    data_loading.R             # load_master(), prep_data(), bin_islet_sizes()
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

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from master_results.xlsx, consumed by Plot, Trajectory, Statistics
- `selected_islet` - reactiveVal, written by Plot and Trajectory click handlers
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer

### Critical: Segmentation Rendering is Root-Level

The `output$islet_segmentation_view` renderPlot is defined in **app.R** at the root level, NOT inside any module. Both the Plot tab (modal) and Trajectory tab (embedded panel) reference this shared output using non-namespaced `plotOutput("islet_segmentation_view")`. This matches the original monolithic pattern and is required for proper output binding.

### Plotly Click Events in Modules

- Plot module: `source = ns("plot_scatter")` with `event_register("plotly_click")`
- Trajectory module: `source = ns("traj_scatter")` with `event_register("plotly_click")`
- Both use namespaced sources; `event_data()` calls must use matching `ns()` prefix

## Key Data Files

- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets)
- `data/pseudotime_results.rds` - R slingshot trajectory output
- `data/adata_ins_root.h5ad` - Python DPT trajectory output (loaded by Trajectory tab)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/CODEX_Pancreas_Donors.xlsx` - Donor metadata key

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- Expression z-scores in UMAP use diverging colormap (blue-white-red, clamped to +/-3 SD)
- Donor status colors: ND = green (#2ca02c), Aab+ = yellow (#ffcc00), T1D = purple (#9467bd)
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 19 phenotypes, 18 markers

## Running the App

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Analysis Pipeline

See `islet_analysis/CLAUDE.md` for the Python analysis pipeline (phenotyping, aggregation, trajectory).

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate
