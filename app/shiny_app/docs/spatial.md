# Spatial Neighborhood Analysis — Detailed Reference

Detailed reference for the Spatial tab. Linked from `app/shiny_app/CLAUDE.md`.

> **Source-artifact versus reviewer UI:** the stored artifacts retain 21 upstream
> phenotype labels and systematic metrics for 20 after excluding Unknown. The reviewer
> UI applies the shared nine-label retirement gate and exposes 12 active phenotypes.
> Counts below describe artifact coverage unless a UI count is stated explicitly.

## Spatial Neighborhood Analysis (Phase 7, Feb 2026 — extended Apr 2026)

### Peri-Islet Metrics (~100 columns, all phenotypes systematic)

| Category | Columns | Description |
|----------|---------|-------------|
| Peri-islet composition | `peri_prop_*`, `peri_count_*`, `total_cells_peri` | Proportion & count for each of 21 phenotypes in 20µm expansion zone |
| Aggregate immune metrics (legacy) | `immune_frac_peri/core`, `immune_ratio`, `cd8_to_macro_ratio`, `tcell_density_peri` | Immune fractions, ratios, density. Columns retained for compatibility but no longer surfaced in any UI. |
| Enrichment z-scores (per phenotype) | `enrich_z_<phenotype>` | Poisson z (observed vs tissue-baseline) for **all 20 phenotypes** (excludes "Unknown"). Extended Apr 2026 from a hand-curated 7-immune subset. |
| Distance metrics (per phenotype) | `min_dist_<phenotype>`, `min_dist_immune_mean` | Min distance from islet core centroid to nearest peri cell of each phenotype, plus the aggregate "any immune" mean across the 7 IMMUNE_TYPES. |

### Data Flow

`compute_neighborhood_metrics.py` reads single-cell H5AD -> `neighborhood_metrics.csv` (5,214 rows) -> merged into `islet_explorer.h5ad` `.obs` by `build_h5ad_for_app.py` (step 4.5) -> extracted in `data_loading.R` `load_master_h5ad()` -> merged into `comp` by `prep_data()` -> available in Plot composition selector (4 option groups) + Statistics tab + Spatial tab.

### Key Details

- All 5,214 islets have peri-islet data (100% coverage with min_cells=0 + require_paired=True)
- Column naming sanitizes phenotype names: spaces -> `_`, `+` -> `plus` (e.g., `peri_prop_ECADplus`)
- Immune signal validated: T1D immune_frac_peri (0.155) > Aab+ (0.106) > ND (0.069)
- **Plot composition selector** has 3 option groups: Core Phenotype Proportions, Peri-Islet Proportions, and Core+Peri Phenotype Proportions. QuPath threshold-based INS/GCG/SST measures are excluded because original and added donors use different threshold pipelines. The 5 aggregate ratio metrics were removed (Apr 2026) in favour of phenotype-driven Spatial UI selectors.
- **`IMMUNE_TYPES`** in `compute_neighborhood_metrics.py` (CD8a Tcell, CD4 Tcell, T cell, B cell, Macrophage, APCs, Immune) is intentionally narrow — drives only aggregates where "what counts as immune" needs definition. **`per_type_phenotypes`** (all `obs.phenotype.unique()` minus `PER_TYPE_EXCLUDE = {"Unknown"}`) drives `enrich_z_*` and `min_dist_*` per-phenotype loops.

### Spatial Tab Layout (Phase 9+13+15, mod_spatial_ui/server.R)

Row 1: Data / View control panels, then Tissue Scatter (col-7) + Leiden clusters (col-5; `cluster_view`: Composition, the default, or membership by donor status or donor); a full-width row with the Leiden UMAP and the same layout coloured by band immune share, side by side at one size (2026-10-19, `spatial_umap_plotly()`; the immune colouring replaced a Donor Status UMAP identical to the Trajectory tab's on 2026-10-20). Below, the Neighborhood Analysis heading and object count, then three card rows: (1) Insulitis prevalence | Immune cells and beta-cell content | Immune gradient around islets; (2) one Phenotype composition & enrichment card (status bar chart and donor heatmap); (3) a full-width Distance-to-islet profile card (per-cell view only). Each card's prose lives in its own "?" popover (`card_help()`). The A/B/C labels below are historical names.

**Owner review, 2026-10-20** ("be sure each one serves a unique purpose ... If there is redundancy, suggest a different graph"). Measured on v7 (4,763 objects): Peri enrichment by status was one column of the composition card; Peri vs core proportion was degenerate (core immune share zero in 96–100% of objects, band endocrine share zero in 84–100%); Distance to nearest cell was taken from the core's centre, so it rank-correlated 0.78–0.86 with object radius, and its CD8 status difference (p = 0.006) vanished once radius was subtracted (p = 0.45); Cluster Composition plotted the endocrine shares the clusters are built from (eta² 0.61–0.92 for those, ≤ 0.03 for immune shares) — the owner kept it ("There still needs to be a cluster composition plot"), as it is how a cluster is read, with the membership views beside it; the Donor Status UMAP duplicated the Trajectory tab's; the KDE view's table holds only core and band cells, so it ended at 20 µm. The cohort-baseline enrichment z rank-correlates 0.80–0.91 with each donor's tissue share, so the composition card gained an Own tissue baseline (the owner kept the cohort one). Replacements are donor-level and live in `R/spatial_islet_immune.R` (tests: `tests/testthat/test-spatial-islet-immune.R`); `donor_distance_profile()` aggregates the cell-distance table once per worker (DuckDB on Parquet, `data.table` on the development CSV) without the retired-label filter, so "All immune" counts the unsubtyped T cell / Immune cells. The Spatial outlier controls and their app.R sync went with the outlier box plots. The signed-distance KDE table stays in the release contract (loader, validators) but no card reads it.

1. **Controls Sidebar** (col-2) -- donor selector, color-by (phenotype/Leiden), Leiden resolution, region filter (All/Core+Peri/Core), "Color background cells" checkbox, phenotype filter (checkboxes with All/None links), donor status checkboxes, phenotype + donor palette selectors, download button.
2. **Tissue Scatter** (col-6) -- ggplot2 `renderPlot` (NOT plotly) showing ~177K cells/donor. Background tissue cells in light grey or colored by phenotype (toggle). Foreground (core/peri) colored by phenotype or Leiden cluster. `coord_fixed() + scale_y_reverse()`. Height: 800px. Supports brush-to-zoom.
3. **Leiden Panel** (col-4) -- plotly UMAP of 5,214 islets colored by selected Leiden resolution (0.3/0.5/0.8/1.0) + stacked bar chart of mean phenotype composition per cluster. UMAP uses raw marker PCA visualization coords (same as trajectory). *Frozen-release description.* On a uniform root from v6 (owner, 2026-10-19) the UMAPs and Leiden are the endocrine composition of the objects with ten or more core cells (the Beta, Alpha, Delta and CHGA+ hormone-negative shares of each core's endocrine cells; `scripts/uniform/composition_embedding.py`); smaller objects have no point or cluster, Cluster Composition leaves them out, and the method text comes from `uns['composition_embedding']`.
4. **Neighborhood Analysis Cards** (3 sections, conditionally rendered via `has_neighborhood()` guard). The card selectors are phenotype-driven, dynamically populated from `prepared()$comp`, and passed through the shared retained-phenotype gate.
   - **Row 1 (2026-10-20)**: `insulitis_bars` (`insulitis_prevalence()`, slider `insulitis_threshold`, default 15), `beta_immune` (`immune_by_beta_content()`, slider `beta_cut`, default 0.1) and `immune_gradient` (`immune_gradient()` on `donor_distance_profile()`, select `gradient_phenotype`), each with a CSV download.
   - **Card B: Phenotype Composition & Enrichment** -- `enrich_donor_values()` (one row per donor and phenotype: the mean over its objects) feeds `enrich_summary()` (the bar chart) and the donor heatmap. Region toggle: peri enrichment z, core proportions, peri proportions; `enrich_baseline` (Cohort / Own tissue, `own_tissue_enrich_z()`) for the enrichment z. Band modes skip the endocrine labels. The shared phenotype filter (`enrich_pheno_filter`) governs both plots; an empty selection blanks them. `enrich_clip` drops object z beyond ±5 before the donor mean (it always did; the help used to call it display saturation).
   - **Immune phenotypes only (owner, 2026-10-20)**: Card B and Card C leave out `NON_IMMUNE_LABELS` (Beta, Alpha, Delta, CHGA+ hormone-negative, Exocrine; `R/spatial_islet_immune.R`). Card C's Marker+ list is every `<marker>_pos` column of the root: all 31 panel markers from v8 (`assemble_uniform.py` POS_MARKERS; v7 had 18), coloured by `CELL_DIST_MARKER_COLORS`. Its composition facets share one y range, from the tallest stacked bar of any status.
   - **Profile range (owner, 2026-10-20)**: both views stop at `CD_MAX_OUTSIDE_UM` (200 um) outside the border; composition bars are `CD_COMP_BIN_UM` (5 um) wide and zero-anchored, density uses `CD_HIST_BINS` (160) symlog bins; a dotted line marks `CD_BAND_UM` (20 um), the band edge (`spatial_helpers.R`). The memory and DuckDB queries stay in step (`test-runtime-data-backend.R`).
   - **Card C: Distance-to-islet profile** -- the Senior-ported per-cell view only (~2.65M tissue cells, symlog density or zero-anchored composition; `dist_mode` phenotypes versus marker positivity, `dist_norm` class proportion versus raw count). The picker observer depends on `dist_mode`. Production queries the cell-distance Parquet through embedded DuckDB; CSV is development-only.
   - Cards carry plain `h5` titles; the former `section_heading()` badges were removed in Oct 2026.
   - Distance NA values: caused by zero immune cells in peri-islet zone (20.1% NA for any immune, 88.6% for CD8+ T-cells). Biological, not a bug. Min-cells filter at ≥50 reduces NA from 20% to 6%.

Wired as `spatial_server("spatial", prepared, active_palette, active_donor_colors)` (no outlier arguments since 2026-10-20).

#### Neighborhood Cards Server Architecture
- `nbr_comp()` shared reactive -- filters `prepared()$comp` by `input$groups` (the object classes are applied in `prepared()`). Reused by every neighborhood output except the immune gradient and the per-cell profile, which read per-cell data.
- `phenotype_options()` reactive -- builds retained phenotype records from `prepared()$comp` columns. Each entry: `label`, `safe`, `enrich`, `prop_core`, `peri_prop`, `peri_count`, `min_dist`. Recovers labels from sanitised suffixes, then applies the retired-label gate (`phenotype_is_kept()`). Until Oct 2026 this step was documented but missing, so every neighborhood picker, the Card B filter and both Card B plots listed all 20 source labels. `canonical_phenotype_order()` (tissue-map filter, per-cell distance pickers) also gates its observed extras now; previously retired labels re-entered through that path. The Leiden cluster-composition plot drops retired `prop_*` columns too.
- **No inline marker hints** (Oct 2026): the four `*_pheno_hint` outputs were removed. `phenotype_marker_definitions_ui()` (`00_globals.R`) lists the retained phenotypes' rules in the Spatial Data controls help popover, composed from `protein_phenotype_rules.csv` when the data root ships it, else from the Rules1 reference.
- `phenotype_choices()` reactive -- named vector (labels -> safe-names) behind `phenotype_options()` lookups.
- `enrich_donor_values()` / `enrich_summary()` -- see Card B above. Respect `input$enrich_region`, `enrich_baseline`, `enrich_clip`, `enrich_stat`.
- Every output uses `donor_colors_reactive()` for palette syncing.
- **Plotly categorical axis ordering**: Must set `categoryorder = "array", categoryarray = c("ND", "Aab+", "T1D")` on x-axis. Plotly defaults to alphabetical, which puts "Aab+" before "ND".
- *(Historical)* The May 2026 KDE selector race fix and the later choices-in-renderUI fix applied to the four phenotype dropdowns removed on 2026-10-20.

#### Supporting Files
- `spatial_helpers.R` -- `donor_tissue_available()`, `get_available_donors()`, `load_donor_tissue(imageid)` with `new.env()` caching
- `data/donors/{imageid}.csv` -- 15 per-donor CSVs (5 cols: X_centroid, Y_centroid, phenotype, cell_region, islet_name), ~78 MB total
- `islets_core_clustered.h5ad` -- Leiden clustering source (4 resolutions: leiden_0.3/0.5/0.8/1.0)

#### Tissue Scatter Design
- Uses `ggplot2::renderPlot` NOT plotly -- 177K points would freeze plotly
- Foreground/background layering: tissue cells at `size=0.15, alpha=0.3` in grey; core/peri at `size=0.4, alpha=0.6` in color
- Leiden coloring maps `islet_name -> cluster` via islet-level lookup from `prepared()$comp`
- **Brush-to-zoom**: `plotOutput(brush=brushOpts(...), dblclick=...)` + `reactiveValues(xmin,xmax,ymin,ymax)`. When zoomed, switches from `coord_fixed()` to `coord_cartesian(xlim, ylim)` with `sort()` on ylim for reversed y-axis.
- **Phenotype filter**: `checkboxGroupInput` with All/None action links. Dynamically populated from donor's unique phenotypes.
- Donor 6533 has 191 islets (9,815 core + 13,840 peri cells) — fully integrated after patching Parent annotations.
