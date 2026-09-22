#!/usr/bin/env Rscript
# Explicit candidate source/output. Does not touch or borrow frozen summaries.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: build_runtime.R REPOSITORY_ROOT CANDIDATE_ROOT")
repo <- normalizePath(args[1], mustWork = TRUE)
root <- normalizePath(args[2], mustWork = TRUE)
if (!file.exists(file.path(root, "candidate.json"))) stop("Explicit candidate identity required")
schema <- if (file.exists(file.path(root, "candidate.json"))) jsonlite::read_json(file.path(root, "candidate.json"))$schema_version else "3"
Sys.setenv(ISLET_DATA_ROOT = root, ISLET_DATA_SCHEMA = schema, ISLET_RUNTIME_MODE = "development")
setwd(file.path(repo, "app/shiny_app"))
source("R/00_audit.R"); source("R/00_globals.R"); source("R/01_runtime_data_backend.R")
source("R/utils_safe_join.R"); source("R/data_loading.R")
if (identical(as.character(schema), "4")) {
  python <- Sys.getenv("RETICULATE_PYTHON", unset = Sys.which("python3"))
  status <- system2(python, c(shQuote(file.path(repo, "scripts/research/build_runtime_schema4.py")), "--input-root", shQuote(root)))
  if (!identical(status, 0L)) stop("Canonical schema-4 summary materialization failed")
} else {
master <- load_master(file.path(root, "master_results.xlsx"))
master <- lapply(master, function(d) {
  if (!is.null(d) && "Case ID" %in% names(d)) d$`Case ID` <- as.integer(d$`Case ID`)
  d
})
prepared <- prep_data(master)
cohort <- read.csv(file.path(root,"cohort.csv"), colClasses = "character")
con <- islet_duckdb_connection()
for (name in names(ISLET_SUMMARY_ARTIFACTS)) {
  d <- as.data.frame(prepared[[name]])
  if (!ncol(d)) d <- data.frame(`Case ID` = integer(), islet_key = character(), check.names = FALSE)
  if (nrow(d) && "Case ID" %in% names(d) && !setequal(as.character(d$`Case ID`), cohort$case_id)) stop("Mixed cohort in ", name)
  target <- file.path(root, ISLET_SUMMARY_ARTIFACTS[[name]])
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(target)) unlink(target)
  duckdb::duckdb_register(con, "candidate_frame", d)
  DBI::dbExecute(con, paste0("COPY candidate_frame TO ", DBI::dbQuoteString(con, target), " (FORMAT PARQUET, COMPRESSION ZSTD)"))
  duckdb::duckdb_unregister(con, "candidate_frame")
}
}
cohort <- read.csv(file.path(root,"cohort.csv"), colClasses = "character")
dir.create(file.path(root, "geometry"), showWarnings = FALSE)
for (cid in cohort$case_id) {
  key <- sprintf("%04d", as.integer(cid))
  target <- file.path(root, "geometry", paste0(key, ".gpkg"))
  if (file.exists(target)) next
  g <- sf::st_read(paste0("/vsigzip/", file.path(root, "gson", paste0(key, ".geojson.gz"))), quiet = TRUE, options = "ARRAY_AS_STRING=YES")
  sf::st_crs(g) <- NA
  sf::st_write(g, target, layer = "segmentation", layer_options = "SPATIAL_INDEX=YES", quiet = TRUE)
  message("Geometry indexed: ", key)
}
islet_disconnect_backend()
