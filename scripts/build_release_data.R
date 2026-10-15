#!/usr/bin/env Rscript
# Build an immutable, self-contained Islet Explorer runtime-data release.
#
# Usage:
#   Rscript scripts/build_release_data.R \
#     --source-root data \
#     --image-root app/shiny_app/www/local_images \
#     --output /path/to/islet-data-v1 \
#     --materialize copy
#
# The output path must not already exist. "copy" is the production default;
# "hardlink" is for space-saving local tests and falls back to copying across
# filesystems. Do not mutate canonical inputs after making a hard-linked build.
# A dirty worktree is rejected; --allow-dirty true exists only for test builds.

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list(materialize = "copy", schema = "1", allow_dirty = "false")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) stop("Unexpected argument: ", arg, call. = FALSE)
    if (grepl("=", arg, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- bits[[1]]
      value <- paste(bits[-1], collapse = "=")
    } else {
      key <- sub("^--", "", arg)
      i <- i + 1L
      if (i > length(args)) stop("Missing value for --", key, call. = FALSE)
      value <- args[[i]]
    }
    out[[gsub("-", "_", key, fixed = TRUE)]] <- value
    i <- i + 1L
  }
  out
}

script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) return(normalizePath("scripts/build_release_data.R", mustWork = TRUE))
  normalizePath(sub("^--file=", "", arg[[1]]), mustWork = TRUE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing builder packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

tree_files <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files[file.info(files)$isdir %in% FALSE]
}

sha256_tree <- function(root) {
  files <- tree_files(root)
  rel <- substring(files, nchar(normalizePath(root, mustWork = TRUE)) + 2L)
  ord <- order(rel)
  records <- vapply(ord, function(i) paste(rel[[i]], sha256_file(files[[i]]), sep = "  "),
                    character(1))
  digest::digest(paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

copy_one <- function(source, target, mode) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  source_real <- normalizePath(source, mustWork = TRUE)
  ok <- FALSE
  if (identical(mode, "hardlink")) {
    ok <- isTRUE(file.link(source_real, target))
  }
  if (!ok) ok <- isTRUE(file.copy(source_real, target, overwrite = FALSE, copy.mode = TRUE))
  if (!ok) stop("Could not materialize ", source, " -> ", target, call. = FALSE)
  invisible(target)
}

copy_tree <- function(source, target, mode) {
  if (!dir.exists(source)) stop("Required runtime directory is missing: ", source, call. = FALSE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  files <- tree_files(source)
  root <- normalizePath(source, mustWork = TRUE)
  for (file in files) {
    rel <- substring(file, nchar(root) + 2L)
    copy_one(file, file.path(target, rel), mode)
  }
  invisible(target)
}

sql_string <- function(con, value) as.character(DBI::dbQuoteString(con, value))

write_frame_parquet <- function(con, frame, target, view_name) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  duckdb::duckdb_register(con, view_name, as.data.frame(frame))
  on.exit(try(duckdb::duckdb_unregister(con, view_name), silent = TRUE), add = TRUE)
  DBI::dbExecute(con, paste0(
    "COPY (SELECT * FROM ", as.character(DBI::dbQuoteIdentifier(con, view_name)), ") TO ",
    sql_string(con, normalizePath(target, mustWork = FALSE)),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"
  ))
  duckdb::duckdb_unregister(con, view_name)
  on.exit(NULL, add = FALSE)
  invisible(target)
}

csv_to_parquet <- function(con, source, target, character_case_id = FALSE,
                           preserve_blank_islet_name = FALSE) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  types <- if (character_case_id) ", types={'case_id':'VARCHAR'}" else ""
  # fread(), used by the legacy donor-tissue path, keeps empty character
  # islet_name values as strings. Do not change numeric-field null handling.
  blank_strings <- if (preserve_blank_islet_name) ", force_not_null=['islet_name']" else ""
  sql <- paste0(
    "COPY (SELECT * FROM read_csv(", sql_string(con, normalizePath(source, mustWork = TRUE)),
    ", header=true, sample_size=-1", types, blank_strings, ")) TO ",
    sql_string(con, normalizePath(target, mustWork = FALSE)),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"
  )
  DBI::dbExecute(con, sql)
  invisible(target)
}

parquet_rows <- function(con, path) {
  as.numeric(DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n FROM read_parquet(",
    sql_string(con, normalizePath(path, mustWork = TRUE)), ")"
  ))$n[[1]])
}

file_artifact <- function(path, root, row_count = NULL, columns = NULL,
                          startup_verify = TRUE) {
  rel <- substring(normalizePath(path, mustWork = TRUE),
                   nchar(normalizePath(root, mustWork = TRUE)) + 2L)
  out <- list(
    type = "file", relative_path = rel,
    bytes = unname(file.info(path)$size), sha256 = sha256_file(path),
    startup_verify = startup_verify
  )
  if (!is.null(row_count)) out$row_count <- unname(as.numeric(row_count))
  if (!is.null(columns)) out$columns <- as.list(as.character(columns))
  out
}

tree_artifact <- function(path, root, startup_verify = TRUE) {
  files <- tree_files(path)
  rel <- substring(normalizePath(path, mustWork = TRUE),
                   nchar(normalizePath(root, mustWork = TRUE)) + 2L)
  list(
    type = "tree", relative_path = rel, file_count = length(files),
    bytes = unname(sum(file.info(files)$size)), sha256 = sha256_tree(path),
    startup_verify = startup_verify
  )
}

build_release <- function(options) {
  repo_root <- dirname(dirname(script_path()))
  git_commit <- tryCatch(system2("git", c("-C", repo_root, "rev-parse", "HEAD"),
                                 stdout = TRUE, stderr = FALSE)[[1]],
                         error = function(e) NA_character_)
  git_dirty <- tryCatch(length(system2("git", c("-C", repo_root, "status", "--porcelain"),
                                       stdout = TRUE, stderr = FALSE)) > 0,
                        error = function(e) NA)
  allow_dirty <- tolower(options$allow_dirty %||% "false") %in% c("1", "true", "yes")
  if (isTRUE(git_dirty) && !allow_dirty) {
    stop("Refusing to build a release from a dirty Git worktree. Commit/tag the intended ",
         "reviewer release, or use --allow-dirty true only for a non-production test.",
         call. = FALSE)
  }
  source_root <- normalizePath(options$source_root %||% file.path(repo_root, "data"),
                               mustWork = TRUE)
  image_root <- normalizePath(options$image_root %||%
                                file.path(repo_root, "app", "shiny_app", "www", "local_images"),
                              mustWork = TRUE)
  if (is.null(options$output) || !nzchar(options$output)) {
    stop("--output is required and must name a new versioned release directory.", call. = FALSE)
  }
  output <- normalizePath(options$output, mustWork = FALSE)
  if (file.exists(output) || dir.exists(output)) {
    stop("Output already exists; release builds never overwrite: ", output, call. = FALSE)
  }
  mode <- tolower(options$materialize %||% "copy")
  if (!(mode %in% c("hardlink", "copy"))) {
    stop("--materialize must be hardlink or copy", call. = FALSE)
  }
  schema <- as.character(options$schema %||% "1")

  parent <- dirname(output)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output parent: ", parent, call. = FALSE)
  }
  staging <- tempfile(pattern = paste0(".", basename(output), ".building-"), tmpdir = parent)
  if (!dir.create(staging, recursive = FALSE)) stop("Cannot create staging directory", call. = FALSE)
  promoted <- FALSE
  on.exit(if (!promoted && dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)

  # Load the canonical H5AD through the same reconstruction and prep code used
  # by the app. Do not source 00_globals.R: it initializes the complete UI stack.
  data_root <<- source_root
  data_path <<- function(...) file.path(data_root, ...)
  master_path <<- data_path("master_results.xlsx")
  h5ad_path <<- data_path("islet_explorer.h5ad")
  project_root <<- repo_root
  if (!exists("%||%", mode = "function")) {
    assign("%||%", function(x, y) if (is.null(x)) y else x, envir = globalenv())
  }
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readxl)
  })
  source(file.path(repo_root, "app", "shiny_app", "R", "utils_safe_join.R"), local = globalenv())
  source(file.path(repo_root, "app", "shiny_app", "R", "01_runtime_data_backend.R"), local = globalenv())
  source(file.path(repo_root, "app", "shiny_app", "R", "data_loading.R"), local = globalenv())

  master <- load_master_h5ad(h5ad_path)
  if (is.null(master)) stop("Canonical H5AD could not be reconstructed", call. = FALSE)
  prepared <- prep_data(master)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", bigint = "numeric",
                        config = list(threads = "2", preserve_insertion_order = "false"))
  on.exit(if (!is.null(con)) try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  artifacts <- list()

  summary_map <- c(
    core_area = "summary/core_area.parquet",
    targets_all = "summary/targets.parquet",
    markers_all = "summary/markers.parquet",
    comp = "summary/islets.parquet",
    combined_pairs = "summary/combined_pairs.parquet"
  )
  for (name in names(summary_map)) {
    target <- file.path(staging, summary_map[[name]])
    write_frame_parquet(con, prepared[[name]], target, paste0("summary_", name))
    artifacts[[name]] <- file_artifact(target, staging, nrow(prepared[[name]]),
                                       names(prepared[[name]]))
  }

  conversions <- list(
    cell_distance = c("cell_distance.csv", "spatial/cell_distance.parquet"),
    signed_distance_kde = c("signed_distance_kde.csv", "spatial/signed_distance_kde.parquet")
  )
  for (name in names(conversions)) {
    src <- file.path(source_root, conversions[[name]][[1]])
    if (!file.exists(src)) stop("Required source is missing: ", src, call. = FALSE)
    target <- file.path(staging, conversions[[name]][[2]])
    csv_to_parquet(con, src, target, character_case_id = TRUE)
    artifacts[[name]] <- file_artifact(target, staging, parquet_rows(con, target))
  }

  donor_dir <- file.path(source_root, "donors")
  donor_files <- sort(list.files(donor_dir, pattern = "[.]csv$", full.names = TRUE))
  if (!length(donor_files)) stop("No per-donor tissue CSVs found under ", donor_dir, call. = FALSE)
  tissue_rows <- 0
  for (src in donor_files) {
    donor <- sub("[.]csv$", "", basename(src))
    if (!grepl("^[A-Za-z0-9._-]+$", donor)) stop("Unsafe donor filename: ", basename(src))
    target <- file.path(staging, "spatial", "tissue", paste0("case_id=", donor),
                        paste0(donor, ".parquet"))
    csv_to_parquet(con, src, target, preserve_blank_islet_name = TRUE)
    tissue_rows <- tissue_rows + parquet_rows(con, target)
  }
  artifacts$donor_tissue <- tree_artifact(file.path(staging, "spatial", "tissue"), staging)
  artifacts$donor_tissue$row_count <- tissue_rows

  json_dir <- file.path(source_root, "json")
  gson_dir <- file.path(source_root, "gson")
  json_sources <- list.files(json_dir, pattern = "[.]geojson$", full.names = TRUE)
  gson_sources <- list.files(gson_dir, pattern = "[.]geojson[.]gz$", full.names = TRUE)
  donor_ids <- sort(unique(sub("[.]geojson([.]gz)?$", "",
                               basename(c(json_sources, gson_sources)))))
  geometry_sources <- vapply(donor_ids, function(donor) {
    plain <- file.path(json_dir, paste0(donor, ".geojson"))
    zipped <- file.path(gson_dir, paste0(donor, ".geojson.gz"))
    if (file.exists(plain)) plain else zipped
  }, character(1))
  if (!length(geometry_sources)) stop("No segmentation GeoJSON sources found", call. = FALSE)
  geometry_counts <- list()
  for (src in geometry_sources) {
    donor <- sub("[.]geojson([.]gz)?$", "", basename(src))
    target <- file.path(staging, "geometry", paste0(donor, ".gpkg"))
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    geom <- sf::st_read(src, quiet = TRUE, options = "ARRAY_AS_STRING=YES")
    sf::st_crs(geom) <- NA
    sf::st_write(geom, target, layer = "segmentation", driver = "GPKG",
                 layer_options = "SPATIAL_INDEX=YES", quiet = TRUE)
    geometry_counts[[donor]] <- nrow(geom)
  }
  artifacts$geometry <- tree_artifact(file.path(staging, "geometry"), staging)
  artifacts$geometry$feature_counts <- geometry_counts

  passthrough_files <- c(
    "adata_ins_root.h5ad", "barlow_islet_features.csv",
    "islet_spatial_lookup.csv", "phenotype_rules.csv", "donor_metadata.csv"
  )
  for (rel in passthrough_files) {
    src <- file.path(source_root, rel)
    if (!file.exists(src)) stop("Required runtime file is missing: ", src, call. = FALSE)
    target <- file.path(staging, rel)
    copy_one(src, target, mode)
    key <- sub("[.][^.]+$", "", basename(rel))
    artifacts[[key]] <- file_artifact(target, staging)
  }
  copy_tree(file.path(source_root, "cells"), file.path(staging, "cells"), mode)
  artifacts$cells <- tree_artifact(file.path(staging, "cells"), staging)

  # Images are part of the release and are available through data_path("images").
  # Their 100+ GiB tree hash is recorded, but startup does not re-read all TIFFs;
  # the deployment verifier performs that full check before promotion.
  copy_tree(image_root, file.path(staging, "images"), mode)
  artifacts$images <- tree_artifact(file.path(staging, "images"), staging,
                                    startup_verify = FALSE)

  source_files <- c(
    h5ad_path, file.path(source_root, "cell_distance.csv"),
    file.path(source_root, "signed_distance_kde.csv"), donor_files,
    unname(geometry_sources), file.path(source_root, passthrough_files)
  )
  source_files <- unique(source_files[file.exists(source_files)])
  source_hashes <- lapply(source_files, function(path) list(
    path = substring(normalizePath(path), nchar(source_root) + 2L),
    bytes = unname(file.info(path)$size), sha256 = sha256_file(path)
  ))

  manifest <- list(
    schema_version = schema,
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    builder = "scripts/build_release_data.R",
    materialization_mode = mode,
    git_commit = git_commit,
    git_dirty = git_dirty,
    artifacts = artifacts,
    canonical_sources = source_hashes
  )
  jsonlite::write_json(manifest, file.path(staging, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 16)

  DBI::dbDisconnect(con, shutdown = TRUE)
  con <- NULL
  if (!file.rename(staging, output)) stop("Could not promote staging directory to ", output)
  promoted <- TRUE
  message("Release created: ", output)
  message("Schema: ", schema, "; artifacts: ", length(artifacts))
  invisible(output)
}

main <- function() {
  require_packages(c("DBI", "duckdb", "sf", "jsonlite", "digest", "dplyr",
                     "tidyr", "readxl", "anndata"))
  build_release(parse_args(commandArgs(trailingOnly = TRUE)))
}

if (sys.nframe() == 0L) main()
