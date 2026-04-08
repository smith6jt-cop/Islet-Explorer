library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)
library(jsonlite)

# ---- Scaling feature flags ------------------------------------------------
# Toggle these to fall back to the legacy in-memory code paths if something
# regresses. They are read by data_loading.R, spatial_helpers.R, and the
# stats/spatial modules.
USE_DUCKDB <- isTRUE(as.logical(Sys.getenv("ISLET_USE_DUCKDB", "TRUE")))
USE_MIRAI  <- isTRUE(as.logical(Sys.getenv("ISLET_USE_MIRAI",  "TRUE")))

# Pixel size constant for coordinate conversion (micrometers per pixel)
# GeoJSON polygons use pixel coordinates, while islet_spatial_lookup.csv uses micrometers
PIXEL_SIZE_UM <- 0.3774  # micrometers per pixel

# Donor status color palettes
DONOR_COLORS <- c("ND" = "#4477AA", "Aab+" = "#CC6633", "T1D" = "#228833")

DONOR_COLORS_BRIGHT <- c("ND" = "#3366CC", "Aab+" = "#FF6600", "T1D" = "#109618")

DONOR_COLORS_COLORBLIND <- c("ND" = "#0072B2", "Aab+" = "#D55E00", "T1D" = "#009E73")

DONOR_COLORS_CLASSIC <- c("ND" = "#2166AC", "Aab+" = "#B2182B", "T1D" = "#1B7837")

DONOR_PALETTES <- list(
  "Paul Tol (default)" = DONOR_COLORS,
  "Bright"             = DONOR_COLORS_BRIGHT,
  "Okabe-Ito"          = DONOR_COLORS_COLORBLIND,
  "Diverging"          = DONOR_COLORS_CLASSIC
)

get_donor_color_palette <- function(name = "Paul Tol (default)") {
  DONOR_PALETTES[[name]] %||% DONOR_COLORS
}

# Load sf package for GeoJSON/spatial operations (islet segmentation viewer)
SF_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (requireNamespace("sf", quietly = TRUE)) {
    library(sf)
    SF_AVAILABLE <- TRUE
    message("[SF] Package sf loaded for segmentation viewer")
  } else {
    message("[SF] Package 'sf' not found. Islet segmentation viewer will be disabled.")
    message("[SF] Install with: install.packages('sf')")
  }
})

# Check for httr2 package (required for AI assistant)
httr2_available <- requireNamespace("httr2", quietly = TRUE)
if (!httr2_available) {
  message("Package 'httr2' not found. The AI assistant will be disabled until you install it (install.packages('httr2')).")
}

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("Package 'anndata' not found. Install with: BiocManager::install('anndata')")
  }
})

# Load vitessceR for embedded image viewer (GitHub v0.1.0 required)
VITESSCE_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (!requireNamespace("vitessceR", quietly = TRUE)) {
    message("Package 'vitessceR' not found. Install with: remotes::install_github('vitessce/vitessceR')")
    message("Version required: >= 0.1.0")
  } else {
    # Verify we have the GitHub version with required features
    tryCatch({
      pkg_version <- packageVersion("vitessceR")
      if (pkg_version < "0.1.0") {
        warning("vitessceR version ", pkg_version, " detected. Please upgrade to >= 0.1.0 from GitHub")
      }
      test_vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.16", name = "Test")
      message("[VITESSCE] Package verified working (version ", pkg_version, ")")
      VITESSCE_AVAILABLE <<- TRUE
    }, error = function(e) {
      warning("[VITESSCE] Package installed but failed basic test: ", conditionMessage(e))
    })

    # CRITICAL FIX: Ensure all child routes from MultiImageWrapper are registered
    # This fixes segmentation layers not loading properly and appearing at wrong sizes
    if (VITESSCE_AVAILABLE) {
      tryCatch({
        vitessceR::MultiImageWrapper$set("public", "get_routes", function() {
          routes <- list()
          if (length(self$image_wrappers) == 0L) return(routes)
          for (w in self$image_wrappers) {
            wr <- w$get_routes()
            if (wr > 0L) routes <- c(routes, wr)
          }
          routes
        }, overwrite = TRUE)
        message("[VITESSCE] Applied MultiImageWrapper get_routes fix for segmentation layer routing")
      }, error = function(e) {
        warning("[VITESSCE] Failed to apply get_routes fix: ", conditionMessage(e))
      })
    }
  }
})

## for base64 encoding channel_config payload to Avivator
suppressPackageStartupMessages({
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Run scripts/install_shiny_deps.R to install dependencies.")
  }
})

# ---- Runtime diagnostics ----
APP_VERSION <- "modular-v1"
try({
  message("[APP LOAD] Version=", APP_VERSION, " time=", Sys.time())
  message("[APP LOAD] getwd=", tryCatch(getwd(), error = function(e) NA))
  message("[APP LOAD] list.files(.) contains app.R? ", any(grepl('^app.R$', list.files('.'))))
  flag_path <- file.path(dirname(sys.frame(1)$ofile %||% 'app.R'), 'APP_LOADED_FLAG')
  writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
}, silent = TRUE)

# Path constants
master_path <- file.path("..", "..", "data", "master_results.xlsx")
h5ad_path   <- file.path("..", "..", "data", "islet_explorer.h5ad")
parquet_dir <- file.path("..", "..", "data", "parquet")
project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# ---- DuckDB-backed data backend (Phase 1) ---------------------------------
# A single read-only DuckDB connection is created once at app startup and
# shared across all sessions. Parquet files produced by
# `scripts/convert_h5ad_to_parquet.py` are exposed as virtual views so each
# module can write lazy dbplyr pipelines (collect() only when an R-only
# function needs materialised data). When Parquet files are missing or the
# `USE_DUCKDB` flag is FALSE the app falls back to the existing
# H5AD/Excel in-memory loader.
con <- NULL
if (isTRUE(USE_DUCKDB)) {
  have_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
                 requireNamespace("DBI", quietly = TRUE)
  if (!have_duckdb) {
    message("[DUCKDB] Package 'duckdb' not installed; falling back to H5AD loader.")
    USE_DUCKDB <- FALSE
  } else {
    islets_parquet <- file.path(parquet_dir, "islets.parquet")
    tissue_root    <- file.path(parquet_dir, "tissue")
    if (!file.exists(islets_parquet)) {
      message("[DUCKDB] ", islets_parquet, " not found; run scripts/convert_h5ad_to_parquet.py. ",
              "Falling back to H5AD loader for this session.")
      USE_DUCKDB <- FALSE
    } else {
      con <- tryCatch(
        DBI::dbConnect(duckdb::duckdb(), read_only = TRUE),
        error = function(e) { message("[DUCKDB] connect failed: ", conditionMessage(e)); NULL }
      )
      if (is.null(con)) {
        USE_DUCKDB <- FALSE
      } else {
        # Optional spatial extension — used for point-in-bbox click queries on
        # the tissue scatter. Swallow install errors so the app still boots
        # offline (extension lives in DuckDB's CDN by default).
        try(DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;"), silent = TRUE)

        try({
          DBI::dbExecute(con, sprintf(
            "CREATE OR REPLACE VIEW islets AS SELECT * FROM parquet_scan('%s')",
            normalizePath(islets_parquet, winslash = "/", mustWork = TRUE)
          ))
          message("[DUCKDB] Registered view: islets -> ", islets_parquet)
        }, silent = FALSE)

        # Groovy tables — optional; only register if the converter emitted them.
        for (sheet in c("targets", "markers", "composition")) {
          pq <- file.path(parquet_dir, sprintf("uns_%s.parquet", sheet))
          if (file.exists(pq)) {
            try({
              DBI::dbExecute(con, sprintf(
                "CREATE OR REPLACE VIEW uns_%s AS SELECT * FROM parquet_scan('%s')",
                sheet, normalizePath(pq, winslash = "/", mustWork = TRUE)
              ))
            }, silent = TRUE)
          }
        }

        # Partitioned tissue dataset — register only if it exists so the
        # Spatial tab can fall back to CSVs otherwise.
        if (dir.exists(tissue_root)) {
          tissue_glob <- normalizePath(file.path(tissue_root, "**", "*.parquet"),
                                       winslash = "/", mustWork = FALSE)
          try({
            DBI::dbExecute(con, sprintf(
              "CREATE OR REPLACE VIEW tissue AS SELECT * FROM parquet_scan('%s', hive_partitioning=1)",
              tissue_glob
            ))
            message("[DUCKDB] Registered view: tissue -> ", tissue_root)
          }, silent = TRUE)
        }

        # Ensure the connection is closed when the R process exits.
        # (shiny::onStop is only valid inside a running app, so register a
        # reg.finalizer on a harmless environment instead.)
        reg.finalizer(globalenv(), function(e) {
          try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
        }, onexit = TRUE)
      }
    }
  }
}

#' Helper: is the DuckDB backend active for this session?
duckdb_active <- function() isTRUE(USE_DUCKDB) && !is.null(con)

# ---- mirai worker pool for async tasks (Phase 4) --------------------------
# Long-running work (statistical batteries, heatmap resampling) is dispatched
# to a persistent mirai daemon pool so the UI stays responsive. Daemons are
# started once at app load and torn down on exit.
if (isTRUE(USE_MIRAI)) {
  have_mirai <- requireNamespace("mirai", quietly = TRUE)
  if (!have_mirai) {
    message("[MIRAI] Package 'mirai' not installed; async stats disabled.")
    USE_MIRAI <- FALSE
  } else {
    n_workers <- suppressWarnings(as.integer(Sys.getenv("ISLET_MIRAI_WORKERS", "4")))
    if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 4L
    try({
      mirai::daemons(n = n_workers, dispatcher = TRUE)
      message("[MIRAI] Started ", n_workers, " daemon workers")
    }, silent = TRUE)
    reg.finalizer(globalenv(), function(e) {
      try(mirai::daemons(0), silent = TRUE)
    }, onexit = TRUE)
  }
}

# Ensure reticulate/anndata can discover a Python binary when RETICULATE_PYTHON
# isn't set (shiny-server environments often lack PATH).
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  candidate_python <- NULL
  venv_path <- "/home/smith6jt/.local/share/islet-explorer-py/bin/python"
  if (file.exists(venv_path)) candidate_python <- venv_path
  if (is.null(candidate_python) || !nzchar(candidate_python)) {
    candidate_python <- Sys.getenv("PYTHON", Sys.which("python3"))
    if (!nzchar(candidate_python)) candidate_python <- Sys.which("python")
  }
  if (nzchar(candidate_python)) {
    Sys.setenv(RETICULATE_PYTHON = candidate_python)
    message("[PYTHON] RETICULATE_PYTHON=", candidate_python)
  } else {
    message("[PYTHON] No python executable found. Set RETICULATE_PYTHON to a valid interpreter.")
  }
}

# Debug flags
DEBUG_CREDS <- isTRUE(as.logical(Sys.getenv("DEBUG_CREDENTIALS", "0")))
VIEWER_DEBUG_ENABLED <- identical(Sys.getenv("VIEWER_DEBUG", "0"), "1")
