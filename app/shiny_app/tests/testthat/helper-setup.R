# testthat helper — sourced automatically by test_dir() before any test-*.R file.
# Loads the common minimal dependencies and statistical helpers. Backend tests
# source their isolated runtime environments themselves and need no app data.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# base R >= 4.4 ships `%||%`; utils_stats.R's summarize_outlier_filter() uses it.
# Define a fallback on the global environment (where source() places utils_stats.R
# functions) so the suite runs under older R too.
if (!exists("%||%", mode = "function")) {
  assign("%||%", function(a, b) if (is.null(a)) b else a, envir = globalenv())
}

# Locate utils_stats.R regardless of the working directory test_dir() runs in
# (repo root, app/shiny_app, or tests/testthat).
.find_source <- function(rel) {
  cands <- c(rel,
             file.path("app/shiny_app", rel),
             file.path("..", rel),
             file.path("../..", rel))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) {
    stop(sprintf("Could not locate %s from wd=%s", rel, getwd()))
  }
  normalizePath(hit[1])
}

source(.find_source("R/00_audit.R"))
source(.find_source("R/utils_stats.R"))
source(.find_source("R/mod_trajectory_server.R"))
