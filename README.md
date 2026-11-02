# Islet Explorer

Interactive R Shiny application for exploring pancreatic islet CODEX multiplexed imaging data from nPOD and HPAP donors. Compares ND, autoantibody-positive and T1D donor groups using cross-sectional measurements; autoantibody positivity is not treated as a clinical disease stage.

## Current workflow and data status — 2026-11-02

Start with the [complete workflow user guide](docs/workflow_user_guide.md),
[app user guide](docs/user_guide.md), and [current component comparison](docs/v19_v20_component_comparison_20261102.md).
The owner selected **v19 endocrine cells**, with the remaining components compared against v20.
The Lymphatic-only update had unintentionally recalibrated unrelated markers; the correction separates
calibration from live structure membership and checks every component on the real 20-donor data.
Dev last served the unaccepted v20 candidate; dev and production are now stopped and masked during
[container security remediation](docs/container_security_20261102.md). The in-house app remains stopped.
Released security fixes are active in source `23bd2fd`; the owner deferred the unmerged DuckDB
backport. Its remaining decoder vulnerability blocks release even if the package scan passes.
No corrected mixed candidate has been assembled or promoted. The older cohort and release sections below
describe historical artifacts, not the current selection.

<div align="center">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue.svg)](https://github.com/smith6jt-cop/Islet-Explorer)
[![R Shiny](https://img.shields.io/badge/R-Shiny-blue.svg)](https://shiny.rstudio.com/)

</div>

## Historical all-donor candidate — offline and superseded

The `data/development/all-donor-20260917/` candidate was a 21-donor experiment
(7 ND, 7 Aab+, 7 T1D), including donor 78/HDL078. The owner has since
permanently removed that donor from the app and all analysis. This candidate and
its recorded preview at `http://127.0.0.1:3839/islet-explorer/` are historical;
the preview URL is not a current service endpoint and must not be started or
used as a launch candidate. The local app is offline while the remaining 20
donors await uniform reprocessing. See the [owner directive](docs/owner_directive_20261015.md),
[uniform reprocessing plan](docs/uniform_reprocessing_plan_20261015.md), and
[current UFRC launch-preparation note](docs/ufrc_launch_preparation_20261015.md).

The historical preview included 3,592,722 raw objects and 8,371 islets. Its
outcome-blind trajectory, research outputs, and review limitations remain
documented in the [trajectory revision](docs/trajectory_revision_20260917.md)
and [reproducible offline workflow](scripts/research/README.md). These records
are not scientific or UFRC deployment acceptance.

## Frozen publication-review baseline

This 15-donor artifact is retained as historical evidence and is not the selected
launch candidate while the owner-directed 20-donor reprocessing is pending.

Islet Explorer is a non-mutating R Shiny research application for exploring
pancreatic-islet measurements from **15 nPOD donors** across three disease
stages. The trajectory H5AD contains **5,214 observations**; prepared core-area
and composition tables each contain **5,382 rows**. It combines marker and cell-population plots, composition-based trajectory
analysis, spatial analysis, single-cell segmentation drill-down, and browser-streamed
OME-TIFF imagery. Reviewers can download generated CSV and PNG files, but cannot upload
data or alter the versioned scientific release.

The publication-review configuration is intentionally narrower than the development
tree: it exposes Plot, Trajectory, Viewer, and Spatial; AI, Statistics, and Inflammation
are disabled. Production loads a manifest-validated release, keeps small summaries in
process-wide memory, queries immutable Parquet with embedded DuckDB, and retains images
and direct-lookup drill-down files in file formats suited to their access patterns.

### Key Capabilities

- **4 reviewer-facing tabs**: Plot, Trajectory, Viewer, and Spatial
- **Single-cell drill-down**: Click any islet to view segmentation boundaries and cell composition
- **Composition trajectory**: Barlow-style PAGA/ForceAtlas2 and diffusion pseudotime, refit for the selected islet cell-count range
- **Spatial neighborhood analysis**: Immune infiltration, enrichment z-scores, distance metrics
- **Leiden clustering**: 4 resolution levels with UMAP visualization
- **Configurable palettes**: Bright (default), Paul Tol, Okabe-Ito, Diverging — synced across all tabs
- **Demographic filters**: Age, sex, autoantibody status
- **Frozen reviewer release**: Immutable OCI image plus checksummed scientific and image manifests

### Disease Groups

| Group | Color (default) | Description |
|-------|-----------------|-------------|
| **ND** | Steel Blue (#4477AA) | Non-diabetic control donors |
| **Aab+** | Burnt Umber (#CC6633) | Autoantibody-positive (pre-T1D) |
| **T1D** | Forest Green (#228833) | Type 1 diabetic donors |

## Quick Start

> **Current launch gate:** the local app is offline pending acceptance of the
> uniformly reprocessed 20-donor root. Do not use these generic development
> commands to restart a candidate before that gate is cleared.

```bash
# Install R dependencies, including the local diffusion-pseudotime fallback
INSTALL_BIOC=1 Rscript scripts/install_shiny_deps.R

# Run the application
R -q -e "shiny::runApp('app/shiny_app', launch.browser=TRUE)"
```

For local development, the legacy H5AD/Excel loader remains available. Publication
deployment is documented in [the UFRC handoff](PUBAPPS_UFRC_HANDOFF.md) and
[the deployment runbook](deploy/README.md); do not treat the development fallback as a
production release.

## Installation

### Prerequisites

**R Shiny App (R >= 4.2.0):**

| Package | Purpose |
|---------|---------|
| `shiny`, `shinyjs` | App framework |
| `readxl` | Excel fallback data loading |
| `dplyr`, `tidyr`, `stringr` | Data wrangling |
| `ggplot2`, `plotly` | Visualization |
| `broom`, `jsonlite`, `digest` | Stats, JSON, and release-integrity hashes |
| `scales` | Axis formatting |
| `sf` | GeoJSON/GeoPackage segmentation access |
| `DBI`, `duckdb` | Embedded read-only Parquet queries and release loading |
| `anndata`, `reticulate` | H5AD loading and Python trajectory bridge |
| `destiny` | Composition diffusion pseudotime |
| `httr2` | Dormant AI module only; not used in reviewer production |
| `vitessceR` | Embedded OME-TIFF viewer (optional) |

**Production Python runtime:**
- `anndata`, `scanpy`, `igraph`, `leidenalg`, `fa2-modified`, and their pinned dependencies

The broader offline analysis pipeline also uses `scvi-tools`, scikit-learn,
`scib-metrics`, and SciPy.

### Install

```bash
# Base R dependencies for a minimal/non-trajectory developer setup
Rscript scripts/install_shiny_deps.R

# Complete local app environment, including the destiny trajectory fallback
INSTALL_BIOC=1 Rscript scripts/install_shiny_deps.R

# Optional vitessceR package; the bundled Avivator viewer does not require it
INSTALL_VITESSCER=1 Rscript scripts/install_shiny_deps.R

# Optional: Avivator image viewer static bundle
./scripts/install_avivator.sh
```

## Application Tabs

### Plot Tab

Primary exploration interface with sidebar controls shared with Statistics.

- **Scatter mode**: Feature expression vs islet diameter (31 markers, targets)
- **Composition mode**: Phenotype proportions in the islet core, peri-islet band, and combined core+peri region (12 reviewer-facing phenotypes; 21 labels remain in source data). QuPath threshold-based INS/GCG/SST counts are retained in source data for provenance but are not offered as comparable population measures.
- Phenotype labels for the original 15 donors are preserved; the five added donors use reviewed reconstructed rules. See the [expanded-cohort comparability audit](docs/data_comparability_20261014.md) before interpreting pooled comparisons.
- **Filters**: Donor status, autoantibody (GADA, IA2A, ZnT8A, IAA, mIAA), age range, sex, diameter range
- **Normalization**: None, global z-score, robust per-donor (median/MAD×1.4826)
- **Click interaction**: Click any point to open islet segmentation + single-cell drill-down

### Trajectory Tab

The reviewer view uses one composition-derived trajectory rather than a menu of
non-comparable pseudotimes. It follows the Barlow-style workflow: phenotype-composition
features, scanpy diffusion pseudotime, and a PAGA-initialized ForceAtlas2 graph;
`destiny` supplies the explicit local fallback for the pseudotime axis when scanpy is
unavailable. The scatter fits DPT first; Leiden/PAGA/ForceAtlas2 is deferred until the
reviewer asks to show the PAGA layout. Changing the cell-count range refits the trajectory
on that selected cohort.

- **Pseudotime scatter**: Marker expression or cell-type/endocrine-cell ratio against rank-normalized or raw DPT
- **PAGA layout**: Interactive force-directed embedding colored by status, donor, cluster, or feature
- **Cell-count-aware display**: Point size can reflect islet cell count or diameter
- **Trend and heatmap views**: Overall/by-status trends and binned multi-feature heatmap
- **Click to drill-down**: Opens the selected islet's segmentation and cell view

### Statistics Tab (disabled in reviewer production)

Five-section narrative layout with shared sidebar from Plot tab:

1. **Configure Analysis** — Test type, alpha, outlier handling, bin width, diameter range
2. **Primary Results** — Hypothesis table + forest plot (Cohen's d with CI)
3. **Size-Dependent Patterns** — BH-corrected stratified test heatmap + Kendall τ trend analysis
4. **Confounders & Deeper Analysis** — Demographics + AUC analysis
5. **Methods Reference** — Dynamic statistical methods description

### Spatial Tab

Three-panel layout with neighborhood analysis cards:

- **Tissue scatter** (ggplot2, ~177K cells/donor): Foreground/background layering, brush-to-zoom, phenotype filtering
- **Leiden panel**: UMAP of 5,214 islets colored by Leiden cluster (4 resolutions: 0.3/0.5/0.8/1.0) + stacked bar chart
- **Neighborhood analysis cards**:
  - **Immune Infiltration**: Violin plots + peri vs core scatter
  - **Immune Cell Enrichment**: Grouped bar chart (7 immune types × 3 stages) + heatmap
  - **Immune Proximity**: Distance box plots + enrichment vs distance scatter

### Viewer Tab

Embedded OME-TIFF viewer using Avivator/Viv for multiplex image visualization with configurable channel mappings.

### AI Assistant (disabled in reviewer production)

The source tree retains a UF Navigator chat module for development, but the reviewer
container neither mounts an API key nor initializes the module. Enabling it would change
the approved data-flow and requires a new security/privacy review.

### Inflammation Tab (disabled in reviewer production)

The source tree also retains the cohort-level inflammation module. It is not mounted in
the reviewer UI or server wiring.

### Single-Cell Drill-Down

Available from both Plot and Trajectory tabs when clicking an islet:

- **Segmentation view**: GeoJSON boundary polygons from QuPath
- **Single-cell overlay**: 5,214 per-islet CSV files with 37 columns (coordinates, phenotype, 31 markers)
- **Composition summary**: Horizontal bar chart + cell count table (core/peri breakdown)
- **Color by**: Reviewer-facing phenotype or any available protein marker

## Data

### Development sources

`data/islet_explorer.h5ad` is the enriched analytical source used to construct the
reviewer summaries. It contains a 5,214-observation H5AD axis, 31 protein-marker
features, donor metadata, phenotype proportions, neighborhood metrics, UMAP coordinates,
Leiden clusters, and QuPath-derived measurements. `data/master_results.xlsx` remains a
reduced-feature **development fallback only**.

H5AD observation counts do not describe all embedded QuPath measurements. The
2026-09-08 source preparation produces `core_area` and `comp` with 5,382 rows
each, `targets_all` with 48,438 rows, `markers_all` with 64,584 rows, and
21 `combined_pairs` definitions. See [the readiness evidence index](docs/ufrc_readiness.md)
for these source measurements and the pending source-to-release parity check.

The composition trajectory uses `data/barlow_islet_features.csv` plus the prepared app
tables. `data/adata_ins_root.h5ad` is retained as a versioned scientific/provenance
artifact and for legacy validation; it is not a network database.

### Publication runtime release

`scripts/build_release_data.R` creates a new, immutable release directory and never
overwrites an existing one. It contains:

- `manifest.json`: schema, Git state, sizes, row counts, and SHA-256 checksums.
- `summary/*.parquet`: prepared Plot/Trajectory summaries loaded once per R worker.
- `spatial/cell_distance.parquet`, `spatial/signed_distance_kde.parquet`, and partitioned
  donor-tissue Parquet queried by embedded DuckDB.
- `geometry/*.gpkg`: indexed per-donor segmentation geometry for bounded reads.
- `cells/`: direct-lookup per-islet CSVs, retained because point lookup is faster than
  consolidating this workload into DuckDB.
- `images/`: native OME-TIFFs and sidecars for HTTP byte-range delivery.
- Required H5AD, phenotype, metadata, and lookup artifacts.

Production is fail-closed: `ISLET_RUNTIME_MODE=production` requires a complete release
whose schema and checksums validate. There is no H5AD/Excel or CSV fallback in reviewer
production.

### Supporting Data Files

| File | Description |
|------|-------------|
| `data/adata_ins_root.h5ad` | Trajectory source (5,214 islets, pseudotime, UMAP) |
| `data/cells/*.csv` | Per-islet single-cell CSVs (5,214 files, ~203 MB) |
| `data/donors/*.csv` | Development sources for per-donor tissue Parquet |
| `data/neighborhood_metrics.csv` | Per-islet peri-islet metrics (5,214 rows, 62 columns) |
| `data/islet_spatial_lookup.csv` | Islet centroid coordinates |
| `data/json/*.geojson` | Development source for per-donor GeoPackages |
| `data/gson/*.geojson.gz` | Compressed segmentation source fallback |
| `data/phenotype_rules.csv` | Phenotype gating rules; 12 labels are surfaced in the reviewer UI |

### Data Conventions

**Islet keying**: All joins use (`Case ID`, `Donor Status`, `islet_key`). `islet_key` is derived from region names (e.g., `Islet_200` from `Islet_200_core`).

**Region types**: `islet_core`, `islet_band`, `islet_union` — no other spellings.

**Islet diameter**: `islet_diam_um = 2 * sqrt(core_region_um2 / pi)`

**Region synthesis**: If union rows are missing — Targets: synthesize `count` only. Markers: synthesize `n_cells` and `pos_count` as sums, compute `pos_frac`.

**Factor ordering**: ND < Aab+ < T1D (always).

## Data Pipeline

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  │
  ├── scripts/reaggregate_islets.py          → islet-level aggregation + trajectory + Leiden
  │     → islet_analysis/islets_core_fixed.h5ad      (5,214 islets)
  │     → data/adata_ins_root.h5ad                    (+ pseudotime + UMAP)
  │     → islet_analysis/islets_core_clustered.h5ad   (+ Leiden at 4 resolutions)
  │
  ├── scripts/compute_neighborhood_metrics.py → data/neighborhood_metrics.csv
  ├── scripts/extract_per_islet_cells.py      → data/cells/*.csv
  ├── scripts/extract_per_donor_tissue.py     → data/donors/*.csv
  │
  └── scripts/build_h5ad_for_app.py           → data/islet_explorer.h5ad (all app data)
```

### Running Pipeline Scripts

```bash
conda activate scvi-env

# Full rebuild from single-cell H5AD (~15 min)
python scripts/reaggregate_islets.py

# Compute peri-islet neighborhood metrics (~3 min)
python scripts/compute_neighborhood_metrics.py

# Extract per-islet single-cell CSVs (~5 min)
python scripts/extract_per_islet_cells.py

# Extract per-donor tissue CSVs for Spatial tab (~2 min)
python scripts/extract_per_donor_tissue.py

# Build enriched app H5AD (merges all sources)
python scripts/build_h5ad_for_app.py

# Build Excel fallback (from groovy TSV exports)
python scripts/build_master_excel.py
```

After the scientific sources are frozen and validated, build the reviewer runtime as a
new directory. Use `copy` when the destination itself will be the authoritative frozen
release; `hardlink` is suitable only for a disposable local validation build or when the
source tree is independently immutable and the final transfer creates new files.

```bash
Rscript scripts/build_release_data.R \
  --source-root data \
  --image-root app/shiny_app/www/local_images \
  --output /path/to/islet-data-pubapps-2026-01 \
  --materialize copy

Rscript scripts/verify_release_data.R \
  --root /path/to/islet-data-pubapps-2026-01 \
  --mode full
```

## Project Structure

```
Islet-Explorer/
├── app/shiny_app/                 # Modular R Shiny application
│   ├── app.R                      # Entrypoint: UI layout, CSS, server wiring
│   ├── R/
│   │   ├── 00_globals.R           # Package loads, constants, palettes, paths
│   │   ├── 01_runtime_data_backend.R # Manifest validation + embedded DuckDB
│   │   ├── data_loading.R         # load_master_auto(), prep_data(), H5AD loading
│   │   ├── utils_stats.R          # summary_stats(), cohens_d(), pairwise_wilcox()
│   │   ├── utils_safe_join.R      # safe_left_join(), add_islet_key()
│   │   ├── segmentation_helpers.R # GeoJSON cache, segmentation plot builders
│   │   ├── drilldown_helpers.R    # Single-cell overlay, phenotype colors
│   │   ├── spatial_helpers.R      # Per-donor tissue CSV loader with caching
│   │   ├── viewer_helpers.R       # Avivator URL builders, channel config
│   │   ├── ai_helpers.R           # OpenAI/UF Navigator API, streaming
│   │   ├── mod_plot_*.R           # Plot tab module (UI + server)
│   │   ├── mod_trajectory_*.R     # Trajectory tab module
│   │   ├── mod_statistics_*.R     # Statistics tab module
│   │   ├── mod_spatial_*.R        # Spatial tab module
│   │   ├── mod_viewer_*.R         # Viewer tab module
│   │   └── mod_ai_assistant_*.R   # AI chat panel module
│   └── www/                       # Static assets (logo, avivator/)
├── data/                          # Data files (see Data section)
├── scripts/                       # Pipeline and utility scripts
│   ├── build_release_data.R        # Build a new immutable runtime-data root
│   └── verify_release_data.R       # Full/core release integrity verifier
├── deploy/                        # Container, Quadlet, release and rollback tooling
├── islet_analysis/                # Python analysis notebooks and modules
├── notebooks/                     # Validation notebooks (scVI QC, trajectory)
├── docs/
│   ├── user_guide.md              # End-user documentation
│   └── competitive_landscape.md
└── CLAUDE.md                      # Developer reference (architecture, conventions)
```

## Development

### Running Locally

```bash
# From command line
R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"

# Or use the VS Code task: "Run Shiny app to reproduce issues"
```

### Deployment

The reviewer deployment is an immutable, rootless Podman service managed by a
user-mode Quadlet unit. Shiny listens on container port 3838; Podman publishes it
to the requested PubApps VM loopback address `127.0.0.1:3838`, subject to UFRC
confirming the final reachable backend address and port. UFRC supplies TLS,
whole-site Basic Auth, and the public reverse proxy. Code is baked into the OCI
image. The manifest-backed data root and its `images/` subtree are mounted read-only;
only bounded temp and log directories are writable.

Do not deploy by running `git pull` in an active release. Build and scan a clean tagged
image, build and fully verify a new data release, activate its immutable image/data pair,
and retain the previous pair for rollback. The exact commands and external UFRC steps are
in [deploy/README.md](deploy/README.md) and
[PUBAPPS_UFRC_HANDOFF.md](PUBAPPS_UFRC_HANDOFF.md).

### Testing & Diagnostics

```bash
Rscript scripts/test_shiny_prep.R            # Test data loading
Rscript scripts/diagnose_islets.R            # Islet data quality
Rscript scripts/diagnose_marker_regions.R    # Marker region checks
Rscript scripts/na_audit.R                   # Missing value audit
Rscript app/shiny_app/tests/parse_check.R     # Parse every application R file
Rscript scripts/verify_release_data.R --root /path/to/release --mode full
./scripts/clear_shiny_cache.sh               # Clear Shiny cache
```

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `ISLET_RUNTIME_MODE` | `development` or fail-closed `production` |
| `ISLET_DATA_ROOT` | Manifest-backed release root; production default mount is `/srv/islet-explorer/data` |
| `ISLET_DATA_SCHEMA` | Expected release schema version (currently `1`) |
| `ISLET_VERIFY_RELEASE_HASHES` | `none`, `core`, or `full`; production defaults to `core` at startup |
| `ISLET_DUCKDB_THREADS` | Embedded DuckDB threads (production default `2`) |
| `ISLET_DUCKDB_MEMORY_LIMIT` | DuckDB memory cap (production default `1GB`) |
| `ISLET_DUCKDB_TEMP_DIR` | Writable bounded DuckDB spill directory |
| `LOCAL_IMAGE_ROOT` | OME-TIFF directory; production points to the release's `images/` subtree |
| `RETICULATE_PYTHON` | Python binary for reticulate/anndata |
| `ISLET_ENABLE_AI` | Dormant AI module; must be `FALSE` in reviewer production |
| `SHOW_STATISTICS_TAB` | Dormant Statistics tab; must be `0` in reviewer production |
| `SHOW_INFLAMMATION_TAB` | Dormant Inflammation tab; must be `0` in reviewer production |
| `DEBUG_CREDENTIALS` | Development only: set to `1` for fixed-code credential setup diagnostics |
| `VIEWER_DEBUG` | Development only: set to `1` for Viewer diagnostics; production rejects it |

`KEY`, `OPENAI_API_KEY`, and external AI endpoint variables are forbidden in the
reviewer container. Its entrypoint fails if an AI credential is supplied.

### Working with Images

```bash
./scripts/convert_to_ometiff.sh input.tiff output.ome.tiff
./scripts/convert_to_omezarr.sh input.tiff output.zarr
python scripts/annotate_ome_tiff_channels.py input.ome.tiff channel_names.txt
python scripts/auto_color_ome_tiff_channels.py input.ome.tiff
```

## Contributing

### Conventions

- **Data keys**: Preserve (`Case ID`, `Donor Status`, `islet_key`) three-key system
- **Region labels**: Only `islet_core`, `islet_band`, `islet_union`
- **Factor ordering**: ND < Aab+ < T1D
- **Normalization**: `none`, `global z-score`, `robust per-donor` (median/MAD×1.4826)
- **Statistical models**: Global `value ~ donor_status + islet_diam_um`; pairwise on residuals; BH correction
- **Large scatter plots (>50K points)**: Use `ggplot2::renderPlot()` not plotly
- **Coordinate convention**: `coord_fixed() + scale_y_reverse()` for tissue scatter (microscopy convention)

### Architecture Notes

- The Plot sidebar is shared with Statistics via `conditionalPanel`
- Root-level outputs (`islet_segmentation_view`, `islet_drilldown_view`, etc.) are defined in `app.R` and shared by Plot + Trajectory modules with active-tab guards to prevent duplicate DOM IDs
- Donor palette is synced across 3 non-namespaced `selectInput`s via `observeEvent` + `updateSelectInput`

See [CLAUDE.md](CLAUDE.md) for full architectural documentation.

## Troubleshooting

**"Package not found" errors:**
```bash
Rscript scripts/install_shiny_deps.R
```

**Development: "Cannot find islet_explorer.h5ad" / "Cannot find master_results.xlsx":**
- Run the data pipeline scripts (see [Data Pipeline](#data-pipeline))
- Or ensure files exist at `data/islet_explorer.h5ad` and `data/master_results.xlsx`

**Production: "release manifest is missing/incomplete" or checksum mismatch:**
- Do not enable the development fallback.
- Rebuild with `scripts/build_release_data.R`, run the full verifier, and activate the
  resulting versioned directory. A checksum mismatch means that release must not be used.

**Trajectory/PAGA unavailable:**
- Install `destiny` and use the pinned Python environment documented under `deploy/`.
- Ensure `barlow_islet_features.csv` and the required manifest artifacts are present.

**Image viewer not loading:**
```bash
./scripts/install_avivator.sh
./scripts/verify_static_images.sh
```

**App performance issues:**
```bash
./scripts/clear_shiny_cache.sh
```

## License

This project is licensed under the [MIT License](LICENSE). Copyright (c) 2026 smith6jt-cop.

## Citation

If you use Islet Explorer in your research, please cite this repository:

```
Islet Explorer: Interactive Shiny application for exploring human pancreatic islet
CODEX multiplexed imaging data
https://github.com/smith6jt-cop/Islet-Explorer
```

## Acknowledgments

- Developed for analysis of human pancreatic islet CODEX multiplexed imaging data from nPOD donors
- Built with R Shiny, ggplot2, plotly, and sf
- Avivator/Viv viewer for OME-TIFF visualization
- The dormant AI module supports UF Navigator in development; reviewer production makes no AI/API call
- Statistical analysis with R stats and broom

---

**Version**: 3.0 (publication-review deployment architecture)

**Last Updated**: August 2026
