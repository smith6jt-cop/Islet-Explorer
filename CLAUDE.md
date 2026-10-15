# CLAUDE.md - Islet Explorer

> **2026-10-03 owner decision:** the accepted reconstruction rules and all six donor boundaries are approved. Use Harmony on protein PCA, retain the uncorrected comparison, and preserve original source records. Exact historical typing/encoder recovery and scVI retraining are no longer prerequisites. See the [current implementation record](docs/protein_expansion_20261003.md); earlier directions below are historical.

This file provides guidance to Claude Code when working with this repository.

## Restoration checkpoint — 2026-09-23

The default app and original pipeline entry points are restored from `000e970`.
The owner requires exact preservation of the original 15 donors and their
islet identities and values. Historical schema-4 records describe the isolated,
untrusted experiment, not the default app; they require revision and evaluation
before reuse under the phase-2 contract. Keep the original data artifacts
untouched and do not select the experimental schema-4 candidate by default.
Candidate-only app helpers/tests are quarantined under
`archive/schema4-before-restore-20260923/`. See the appended verification record
in `app/shiny_app/docs/trajectory.md`.

## Typing reference and revised owner direction — 2026-09-29 (historical)

*Superseded by the 2026-10-03 owner decision at the top of this file; kept as the record of the September 29 direction.*

Read [the typing/trajectory reference record](docs/trajectory_typing_reference_20260919.md)
before further phenotype or trajectory work. Cell typing must precede the
endocrine-normalized composition features used to fit trajectory. The owner now
directs the six additions to use the original phenotype and scVI methods. The
original numeric gates and input scale were recovered from the frozen H5AD;
the checked-in rule CSV fails a 21-label replay, and the original scVI encoder
weights remain missing. Recover and verify the complete rule chain, then either
locate the encoder or explicitly review a joint scVI retrain. The active UI
rules or QuPath labels are not proof of the original typing method. The
September 19 proposal to enhance phenotypes
with CHGA, CD45, VIM, LGALS3, CD44, PGP9.5, GAP43, CD35 and NaKATPase is deferred,
not an acceptance gate for this expansion. Preserve the frozen release.

## Donor-expansion contract — 2026-09-29 (historical)

*Superseded by the 2026-10-03 owner decision at the top of this file; kept as the record of the September 29 direction.*

[Phase 2](docs/donor_expansion_phase2.md) starts from the restored 15-donor app
and raw-only six-donor checkpoint. Preserve all original donor IDs, spatial
geometry, measured values, missingness and phenotype labels; apply the exact
original typing/scVI method and numeric gates to the additions. A common
21-donor trajectory, UMAP and Leiden refit may change those derived outputs for
the original donors without changing their preserved source records. The owner
will visually review donor annotations. A simple development data root selected
with `ISLET_DATA_ROOT` needs no additional candidate-version or optional-loader
gate. As of 2026-09-29 no six-donor integration or common refit had been accepted
or completed; the 2026-10-03 protein expansion has since been built and accepted. The prior schema-4 geometry, typing and trajectory remain
experimental evidence. See the phase-2 document for checkpoint evidence and
remaining scientific dependencies.

## Project Overview

Interactive Shiny app for exploring pancreatic islet CODEX multiplexed imaging
data from nPOD and HPAP donors. Compares endocrine composition and the islet
microenvironment across ND, autoantibody-positive and T1D donors. Section-level
associations and arbitrary trajectory orientation do not establish progression.

## Sub-directory CLAUDE.md files

- `app/shiny_app/CLAUDE.md` — per-tab feature details (Drill-Down, Trajectory, Plot, AI Assistant) + conventions including the May 2026 font-standardisation system
  - `app/shiny_app/docs/spatial.md` — Spatial Neighborhood Analysis deep reference (split out May 2026)
  - `app/shiny_app/docs/statistics.md` — Statistics Tab deep reference (split out May 2026; Jul 2026 correctness pass: min-cells no-op, p-value notation, emmeans d.f. fallback, AAb pseudoreplication, AUC decoupling; **Aug 2026 at the end: why the per-bin pairwise row could never reject — forced `exact=FALSE` on a 5-v-5 rank test plus a BH family of 15 — and the `status × diameter` interaction p, which was inflated ~11 orders of magnitude by a random-intercept-only error term — each companion term is now tested against its own unit of replication (split-plot)**). The tab itself is hidden as of Aug 2026 (`SHOW_STATISTICS_TAB`).
  - `app/shiny_app/docs/trajectory.md` — **Trajectory/pseudotime accuracy reference** (Aug 2026). Ground-truth measurements of the shipped H5AD, what the audit confirmed (donor η² = 2.05× status; ND vs Aab+ never separable), three refuted claims and two non-viable directions that keep being re-proposed, and the axis/CSV/H5AD-integrity defects that were fixed. **Read before quoting any trajectory parameter — `scripts/CLAUDE.md` was stale until Aug 2026.**
- `scripts/CLAUDE.md` — pipeline scripts, batch correction, pseudotime modes
- `islet_analysis/CLAUDE.md` — Python analysis pipeline (phenotyping, aggregation, trajectory). Tracked via a `!` exception; the rest of `islet_analysis/` is ignored.

## App Architecture (Modular - Feb 2026)

The Shiny app lives in `app/shiny_app/` and uses a modular architecture:

```
app/shiny_app/
  app.R                        # Entrypoint: UI layout, CSS, server wiring, root-level outputs
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
    barlow_pseudotime.R        # destiny fallback for the Barlow composition pseudotime
    paga_layout.R              # reticulate bridge to python/paga_pseudotime.py (scanpy PAGA)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF, Avivator iframe)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (5-section narrative)
                                #   HIDDEN since Aug 2026 — gated on SHOW_STATISTICS_TAB
                                #   (00_globals.R, default 0). Set =1 in .Renviron to restore.
    spatial_helpers.R           # Spatial loaders; Parquet/DuckDB production backend + development CSV cache
    mod_spatial_ui.R / mod_spatial_server.R     # Spatial tab (3-panel: sidebar controls, tissue scatter, Leiden panel)
    mod_inflammation_ui.R / mod_inflammation_server.R  # Inflammation tab (edgeR abundance,
                                #   Bodenmiller infiltration score, nPOD donor cohort table)
                                #   HIDDEN since Aug 2026 — gated on SHOW_INFLAMMATION_TAB
                                #   (00_globals.R, default 0). Set =1 in .Renviron to restore.
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  python/
    paga_pseudotime.py         # scanpy DPT + PAGA/ForceAtlas2 layout (Barlow et al. 2024)
  www/                         # Static assets (images, logos)
```

`app.R` is the ONLY entrypoint — it contains no `source()` calls; Shiny auto-sources `R/`.
The superseded `app_original.R` (pre-modularization monolith) and `app_modular.R` (stale
entrypoint snapshot) were moved to `archive/shiny_app/` in Jul 2026: Shiny never loaded
them, but they sat in the deployed directory where 5,400 lines of dead code read as live.
**`www/logo.png` is a symlink to `app/docs/logo.png`** — that file is a live asset, not a
stray doc.

### Data Loading

Production validates the release manifest and loads prepared summary Parquet through
`01_runtime_data_backend.R`; large spatial tables are queried read-only with embedded
DuckDB. It fails closed and never substitutes a development dataset. Local development
uses `load_master_auto()`, which tries H5AD first and falls back to Excel:
- `h5ad_path` -> `data/islet_explorer.h5ad` (built by `scripts/build_h5ad_for_app.py`)
- `master_path` -> `data/master_results.xlsx` (built by `scripts/build_master_excel.py`)

Both produce the same `list(markers, targets, comp, lgals3, phenotypes, donor_demographics, neighborhood)` structure consumed by `prep_data()`. The `phenotypes`, `donor_demographics`, and `neighborhood` elements are extracted from H5AD `.obs` and are `NULL` when loading from Excel (graceful fallback).

### Performance: Reticulate Bridge Optimization (Mar 2026)

**Critical**: `ad$uns[[key]]` crosses the Python-R bridge on EVERY call (~221ms each). With ~62 groovy columns × 3 sheets = ~186 calls, this took ~16s. Fix: cache `uns <- ad$uns` once (0.2s) then index the R list. `reconstruct_groovy_df_from_list(uns, sheet)` replaces per-key bridge crossings.

**Deferred loading**: `annotations.tsv` (72 MB, 576K rows) was loaded eagerly at startup but only serves as a fallback when `islet_spatial_lookup` doesn't have the islet (never happens). Now lazy-loaded via `.seg_lazy` environment on first fallback access.

Total startup: **~22s → ~3.5s**.

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from H5AD or Excel, consumed by Plot, Trajectory, Statistics, Spatial
- `selected_islet` - reactiveVal with `case_id`, `islet_key`, `centroid_x/y`, `donor_status`, `donor_age`, `donor_gender`
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer
- `active_tab` - `reactive(input$tabs)`, passed to Plot and Trajectory to prevent duplicate output IDs
- `active_donor_colors` - reactive from `donor_palette_name()`, synced across 3 non-namespaced inputs: `sidebar_donor_palette` (Plot sidebar), `traj_donor_palette` (Trajectory), `spatial_donor_palette` (Spatial)

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
- **macOS native-click bridge (Jun 2026)**: `plotly_click` is suppressed on macOS, so both scatters also pipe `islet_click_bridge(ns("<plot|traj>_native_click"))` (native mousedown/mouseup → `Shiny.setInputValue`) after `event_register`. The new server observers and the legacy plotly_click handler both call the shared `select_islet_from_key()` (`segmentation_helpers.R`). See `app/shiny_app/CLAUDE.md` § "macOS Native-Click Bridge + Shared Resolver".

## Key Data Files

- `data/islet_explorer.h5ad` - **Primary app data** (~70 MB, 5,214 islets, groovy + trajectory + donor + neighborhood + Leiden)
- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets) -- Excel fallback
- `data/adata_ins_root.h5ad` - Trajectory source (5,214 islets, 31 vars, UMAP, and legacy audited DPT columns)
  Also carries `dpt_pseudotime_small` / `dpt_pseudotime_large` (Aug 2026): DPT recomputed
  **within** each cell-count stratum (<=10 / >10 cells) by `scripts/compute_split_pseudotime.py`.
  Independent diffusion maps — **not on a common scale with each other or with the all-islet
  column**. See `app/shiny_app/docs/trajectory.md` §9-10, which also documents that pseudotime
  is computed on the 10-d scVI latent, NOT on the 31 markers.
- `data/neighborhood_metrics.csv` - Per-islet peri-islet metrics (5,214 rows, 62 columns)
- `data/barlow_islet_features.csv` - Islet composition matrix (5,214 x 42 log-normalised
  core+peri cell-type counts) for the reviewer trajectory; built by
  `scripts/build_barlow_features.py`, consumed live by `app/shiny_app/R/barlow_pseudotime.R`
- `data/cells/*.csv` - Per-islet single-cell CSVs for drill-down (~5,214 files, ~203 MB total)
- `data/donors/*.csv` - Per-donor tissue-wide cell CSVs for Spatial tab scatter (15 files, ~78 MB total)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/donor_metadata.csv` - Per-donor nPOD clinical metadata (15 rows): T1D duration, BMI, HbA1c,
  C-peptide, AAb results, age group. Built by `scripts/build_donor_metadata.py` from an nPOD
  "All-Download" xlsx that is **deliberately gitignored** (741 donors of free-text clinical
  narrative); only the derived CSV is committed. Donor `112` is filed under `HDL112` upstream.
- `data/DATA_PROVENANCE.md` - Documents canonical H5AD lineage and data sources

## Data Pipeline

### Canonical H5AD Lineage

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  |
  |-- scripts/reaggregate_islets.py  (one-stop: aggregate + trajectory + Leiden)
  |     -> islet_analysis/islets_core_fixed.h5ad  (5,214 islets, min_cells=0, require_paired=True)
  |     -> islet_analysis/islets_peri_fixed.h5ad   (5,214 peri-islet)
  |     -> islet_analysis/islets_merged_fixed.h5ad (5,214 paired core+peri)
  |     -> data/adata_ins_root.h5ad  (+ core + combined pseudotimes, UMAP, X_scVI_harmony_pt, X_scVI_combined_harmony_pt)
  |     -> islet_analysis/islets_core_clustered.h5ad  (+ Leiden at 4 resolutions, Harmony θ=2 on imageid)
  |
  |-- scripts/compute_neighborhood_metrics.py
  |     -> data/neighborhood_metrics.csv  (5,214 rows, ~100 cols — per-phenotype enrich_z and min_dist)
  |
  |-- scripts/extract_per_islet_cells.py
  |     -> data/cells/*.csv  (5,214 files, 37 cols: coords, phenotype, region, 31 markers)
  |
  |-- scripts/extract_per_donor_tissue.py
  |     -> data/donors/*.csv  (15 files, 5 cols: X/Y coords, phenotype, cell_region, islet_name)
  |
  +-- scripts/build_h5ad_for_app.py (trajectory + groovy + donor + neighborhood + Leiden)
        -> data/islet_explorer.h5ad  (all app data in one file, incl. leiden_* + leiden_umap_*)
```

Separately, `scripts/build_donor_metadata.py` derives `data/donor_metadata.csv` from an nPOD
"All-Download" spreadsheet (gitignored). `prep_data()` joins it onto `comp` / `targets_all` /
`markers_all`. **The join drops any column the target already has** — the H5AD and the nPOD sheet
both supply `AAb_GADA/IA2A/ZnT8A/mIAA`, and letting dplyr suffix them `.x`/`.y` silently turned the
Plot sidebar's AAb filter into a no-op (every lookup in this app is by exact column name, so a
collision disables a feature rather than erroring). The two sources agree on all 15 donors x 4
assays, so keeping the existing column loses nothing.

`scripts/reaggregate_islets.py` replaces the manual notebook workflow (aggregation + trajectory + Leiden in one script). `islet_analysis/islets_core_clustered.h5ad` provides Leiden clustering (4 resolutions) + UMAP coords, merged by `build_h5ad_for_app.py` step 4.7.

See `scripts/CLAUDE.md` for per-script details, batch correction history, the audited legacy
H5AD pseudotimes, and the current composition-derived reviewer trajectory.

## Spatial / Statistics / Trajectory / AI tabs

Per-tab implementation details — Spatial Neighborhood Analysis (Phase 7), Single-Cell Drill-Down (Phase 8), Statistics Tab pseudoreplication fix (Phase 16), Cell-Count-Weighted Trajectory (Phase 12), Pseudotime Heatmap, AI Assistant (Phase 10) — live in `app/shiny_app/CLAUDE.md`.

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- UMAP "Selected Feature" uses continuous viridis (inferno) colormap scaled to data min/max
- Donor status colors: ND = steel blue (#4477AA), Aab+ = burnt umber (#CC6633), T1D = forest green (#228833) — Paul Tol bright variant, colorblind-safe, centralized as `DONOR_COLORS` in `00_globals.R`
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 19 phenotypes, 15 markers. CSV headers use adata.var_names (CD8a, HLADR, PDPN — not CD8, HLA-DR, Podoplanin). The upstream Rules1 template also listed `CD99` (Endocrine) and `CD56` (Neural), but neither is in this 31-plex panel (`adata.uns['all_markers']` in `single_cell_analysis/CODEX_Panc_scimap_Analysis.ipynb`), so `scimap.tl.phenotype_cells` silently dropped them — effective gates were `Endocrine = any of INS/GCG/SST+` and `Neural = B3TUBB+`. Both phantom columns were removed from the app's CSV (Jul 2026) so the rules helper matches what was applied. This app-facing CSV is not verified as the complete original typing method for the six additions.
- **Donor metadata**: Comes from `islets_core_fixed.h5ad` obs, NOT from `CODEX_Pancreas_Donors.xlsx` (different cohort)
- **User-facing terminology**: Use "Sex" not "Gender" in all UI labels, plot titles, and methods text. The underlying data column remains `gender` for backwards compatibility.
- **H5AD obs index**: `islets_core_fixed.h5ad` index name is `islet_id` -- same as column, use `reset_index(drop=True)`
- **Case ID zero-padding**: GeoJSON files use `0112.geojson` (4-digit padded), data uses `112` (unpadded). `load_case_geojson()` and click handlers use `sprintf("%04d", ...)` fallback
- **Log-scale with zeros**: Use `scales::pseudo_log_trans(base=10)` instead of `scale_y_log10()` -- zeros map to 0 (visible) instead of -Infinity (dropped)
- **Plot defaults**: Point size = 3.0, transparency = 0.6
- **Plot fonts & downloads**: Every plot uses `theme_islet()` (ggplot) and `plotly_islet_fonts()` + `plotly_islet_config()` (plotly) — single source of truth in `app/shiny_app/R/00_globals.R`. Base 16pt / titles 18pt / axes 14-16pt / legends 14-16pt / heatmap rows 12pt min. Plotly modebar PNG downloads ship 2800 × 1800 px (scale=2). ggplot `renderPlot` outputs that should be downloadable are refactored into a `reactive()` consumed by both `renderPlot` and a `downloadHandler` using `ggsave(..., width=10, height=7, dpi=300)` (heatmaps 12 × 8). See `app/shiny_app/CLAUDE.md` § "Plot Fonts & Publication-Quality Downloads".
- **Phenotype marker rules in UI**: `data/phenotype_rules.csv` is parsed once at startup into `PHENOTYPE_RULES_CSV` in `00_globals.R`, with hierarchical composition (sub-phenotypes inherit parent constraints) and an alias map for the four PHENOTYPE_COLORS-vs-CSV naming mismatches (`CD8a Tcell` ↔ `CD8 T cell`, etc.). Surfaced via (1) inline blue hints under phenotype selectors, (2) a global reference modal containing 19 CSV rules plus 4 explicitly described residual definitions, and (3) an expandable rules list in drill-down. Nine source labels are retired from reviewer selectors, leaving 12 active phenotype labels; the broader modal is an audit reference, not the active-choice list. See `app/shiny_app/CLAUDE.md` § "Phenotype Rules — Marker Definitions in the UI".
- **`targets_all` carries per-islet cell counts** (Jul 2026): `prep_data()` joins `total_cells_core` / `total_cells_peri` from `comp` onto `targets_all`. Without them the "Min cells/islet" filter — guarded on those column names in both `mod_plot_server.R::raw_df_base()` and `mod_statistics_server.R::stats_run()` — *silently did nothing* in Microenvironment mode instead of erroring, while the exported CSV still reported the value the user set. Any new filter guarded on a column's presence must surface the skip (see `min_cells_applied`), never swallow it.
- **Peri-islet data guard**: Always check `total_cells_peri > 0` before using peri metrics
- **Column name sanitization**: Phenotype names use `_` for spaces, `plus` for `+` in peri-islet columns
- **Single-cell Parent column**: `Islet_N` = core cells, `Islet_N_exp20um` = peri-islet cells
- **macOS plotly clicks & Viewer iframe (Jun 2026)**: `plotly_click` is suppressed on macOS — drill-down uses a native-click bridge (`islet_click_bridge()` in `00_globals.R`). The Viewer iframe is mounted once and its `src` updated in place via `session$sendCustomMessage` to avoid recreating the WebGL canvas (macOS WebRender flicker). Details + the `reactiveVal`-dedup gotcha in `app/shiny_app/CLAUDE.md`.

App-specific UI conventions (CSS overflow, font minimums, large-scatter rendering, palette syncing, reticulate bridge, synthetic-union rows, macOS click bridge & Viewer iframe persistence) live in `app/shiny_app/CLAUDE.md`.

## Deployment Architecture

- **Requested reviewer URL (pending UFRC confirmation)**: `https://islet-explorer.rc.ufl.edu/islet-explorer/`
- The project requests that the UFRC high-availability proxy terminate HTTPS and Basic Auth and proxy to VM loopback `127.0.0.1:3838`; the final reachable backend address and port must be agreed in the UFRC ticket before deployment.
- Rootless Podman/Quadlet runs one immutable, digest-pinned image; code, scientific data, and images are read-only.
- Releases are checksummed and atomically activated. Never update a live review round from a mutable branch or symlink.
- Dev server (`Rscript -e 'shiny::runApp(".", port=7777)'`) is for local testing only

## Running the App

Reviewer service, once UFRC provisions the requested URL:
```
https://islet-explorer.rc.ufl.edu/islet-explorer/
```

The authoritative production procedure is `deploy/README.md`; provisioning inputs and
unresolved institutional gates are in `PUBAPPS_UFRC_HANDOFF.md`.

Development only:
```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Running Pipeline Scripts

See `scripts/CLAUDE.md` for the full per-script command list. Quick rebuild:

```bash
conda activate scvi-env
python scripts/reaggregate_islets.py          # ~15 min: aggregation + trajectory + Leiden
python scripts/compute_neighborhood_metrics.py # ~3 min
python scripts/build_h5ad_for_app.py          # ~1 min: merge into app H5AD
```

## Analysis Pipeline

See `islet_analysis/CLAUDE.md` for the Python analysis pipeline (phenotyping, aggregation, trajectory).
See `data/DATA_PROVENANCE.md` for full data lineage documentation.
See `docs/user_guide.md` for end-user documentation.

## Handoffs & Archive

- `docs/HANDOFF.md` — historical Jul/Aug 2026 development handoff. Its analysis-audit record remains useful,
  but its mutable nginx deployment notes are superseded by `deploy/README.md` and
  `PUBAPPS_UFRC_HANDOFF.md`.
- `archive/docs/HANDOFF_2026-07-25_shiny-correctness-pass.md` — superseded predecessor; retains the plotly
  in-plot text findings and the chromote selectize gotcha.
- `scripts/codex_phenocycler/docs/HANDOFF.md` — **active** handoff for the in-flight CODEX Phenocycler port. Read it before touching that workstream's gates, floors, or pairs.
- `archive/README.md` — what was retired and why. Nothing there is loaded, sourced, or served; delete freely.

## Cross-Session Knowledge

Cross-session development knowledge for this project lives in Claude's auto-memory (see `MEMORY.md`). The earlier `Skills_Registry/` submodule + `/advise` + `/retrospective` workflow was retired Apr 2026 in favour of the auto-memory mechanism.

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate, httr2, stringr, lmerTest, lme4, emmeans, destiny, data.table

## Python Dependencies (app venv)

The app's reticulate interpreter (`RETICULATE_PYTHON`, default
`~/.local/share/islet-explorer-py`) needs `anndata` for data loading and, since Aug 2026,
`scanpy igraph leidenalg fa2-modified` for the Trajectory tab's PAGA force-directed layout.
Without scanpy the PAGA card degrades to an explicit "unavailable" message and the composition
pseudotime axis falls back to the destiny implementation in `R/barlow_pseudotime.R`.
Without `fa2-modified`, scanpy silently uses Fruchterman-Reingold instead of ForceAtlas2.
