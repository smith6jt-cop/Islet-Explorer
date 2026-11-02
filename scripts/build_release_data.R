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

# Data-release schemas (app/shiny_app/R/01_runtime_data_backend.R has the runtime side). "1" is the frozen
# 15-donor contract (2026-09). "5" (2026-10-22) is the uniform 20-donor contract: the same summaries and
# spatial tables, plus the uniform root's own files and the H5AD's embedding metadata as JSON. Ids 2-4
# belonged to the retired all-donor experiments and are never reused.
RELEASE_SCHEMAS <- c("1", "5")

#' Files copied as they are, by schema: `required` must exist, `optional` are taken when present. A uniform
#' root (uniform_release.json present) must be released as schema 5: under schema 1 the app would run with
#' the frozen 15-donor behaviour (no uniform typing text, no object-count card, no composition-embedding
#' method text) and never say so.
release_passthrough_files <- function(schema, source_root) {
  schema <- as.character(schema)
  if (!(schema %in% RELEASE_SCHEMAS)) {
    stop("Unsupported --schema ", schema, " (supported: ", paste(RELEASE_SCHEMAS, collapse = ", "), ")",
         call. = FALSE)
  }
  uniform <- file.exists(file.path(source_root, "uniform_release.json"))
  if (identical(schema, "1")) {
    if (uniform) {
      stop("The source root is a uniform 20-donor build (uniform_release.json); release it with --schema 5.",
           call. = FALSE)
    }
    return(list(required = c("adata_ins_root.h5ad", "barlow_islet_features.csv", "islet_spatial_lookup.csv",
                             "phenotype_rules.csv", "donor_metadata.csv"),
                optional = character(0)))
  }
  if (!uniform) {
    stop("--schema 5 needs a uniform root: uniform_release.json is missing under ", source_root, call. = FALSE)
  }
  list(required = c("adata_ins_root.h5ad", "barlow_islet_features.csv", "islet_spatial_lookup.csv",
                    "donor_metadata.csv", "uniform_release.json", "endocrine_objects_by_donor.csv",
                    "protein_phenotype_rules.csv"),
       optional = "phenotype_rules.csv")
}

# The image bundle: every file of the image tree plus, for each OME-TIFF that is a symlink, the XML and
# offsets sidecars beside its target (a uniform root links its TIFFs from other roots and carries no sidecars;
# the Viewer loads faster with the offsets). One row per released file: its path in the release and the
# file it is copied from. A complete tree (the frozen local_images) comes out as it is.
IMAGE_SIDECAR_SUFFIXES <- c(".ome.xml", ".ome.tiff.offsets.json", ".offsets.json")
release_image_files <- function(image_root) {
  root <- normalizePath(image_root, mustWork = TRUE)
  files <- tree_files(root)
  rel <- substring(files, nchar(root) + 2L)
  out <- data.frame(relative_path = rel, source = files, stringsAsFactors = FALSE)
  for (i in grep("[.]ome[.]tiff?$", rel, ignore.case = TRUE)) {
    target <- normalizePath(files[[i]], mustWork = TRUE)
    source_stem <- sub("[.]ome[.]tiff?$", "", basename(target), ignore.case = TRUE)
    release_stem <- sub("[.]ome[.]tiff?$", "", rel[[i]], ignore.case = TRUE)
    for (suffix in IMAGE_SIDECAR_SUFFIXES) {
      sidecar_rel <- paste0(release_stem, suffix)
      sidecar_src <- file.path(dirname(target), paste0(source_stem, suffix))
      if (!(sidecar_rel %in% out$relative_path) && file.exists(sidecar_src)) {
        out <- rbind(out, data.frame(relative_path = sidecar_rel, source = sidecar_src,
                                     stringsAsFactors = FALSE))
      }
    }
  }
  out <- out[order(out$relative_path, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' The digest the materialized image tree will have (sha256_tree() of the copy), computed from the sources.
sha256_image_bundle <- function(image_root) {
  bundle <- release_image_files(image_root)
  records <- vapply(seq_len(nrow(bundle)), function(i) {
    paste(bundle$relative_path[[i]], sha256_file(bundle$source[[i]]), sep = "  ")
  }, character(1))
  digest::digest(paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

#' GDAL data source for a GeoJSON file: the path, or GDAL's gzip virtual file system for a compressed one (a
#' uniform root carries only gson/*.geojson.gz; the app's loader reads them the same way).
release_geometry_dsn <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (grepl("[.]gz$", path)) paste0("/vsigzip/", path) else path
}
# The sources are read whole, and one feature can exceed GDAL's default single-feature limit for GeoJSON
# (OGR_GEOJSON_MAX_OBJ_SIZE): the v20 root's largest Collagen IV polygon (Aab_6505) is about 25 MB of GeoJSON, and
# the read failed with "GeoJSON object too complex/large". Lift the limit for this process (the parity check
# loads this file too).
Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = "0")

#' The H5AD bridge returns one-dimensional arrays; JSON has plain vectors. Applied before writing and before
#' comparing, so the release's metadata equals the source's after one round trip.
release_plain_list <- function(x) {
  if (is.list(x)) return(lapply(x, release_plain_list))
  if (is.array(x) && length(dim(x)) == 1L) return(as.vector(x))
  x
}

#' summary/embedding_metadata.json: the H5AD `uns` records prep_data() returns (composition_embedding,
#' protein_embedding), which production reads back through islet_read_embedding_metadata().
write_embedding_metadata <- function(master, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  meta <- release_plain_list(list(composition_embedding = master$composition_embedding,
                                  protein_embedding = master$protein_embedding))
  # 17 significant digits: every double reads back as the value the H5AD holds.
  jsonlite::write_json(meta, target, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null",
                       digits = I(17))
  invisible(target)
}

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
  # Byte order, as the runtime's .islet_sha256_tree(): the digest must not depend on the host's locale.
  ord <- order(rel, method = "radix")
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
  passthrough <- release_passthrough_files(schema, source_root)

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
  # prep_data() adds the phenotype group aggregates ("All endocrine cells", "All T cells", "All macrophages"; app
  # commit 82ce23e), whose root-aware helpers live in 00_globals.R. Evaluate exactly those definitions, as
  # verify_release_parity.R does, rather than the UI stack.
  aggregate_names <- c("MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN", "is_macrophage_label",
                       "root_phenotype_labels", "root_has_marker_profile_labels", "TOP_MARKER_BASE",
                       "CD45_RULE_LABEL", "CD45_TOP_MARKER_LABEL", "root_has_top_marker_labels",
                       "root_cd45_is_rule_label", "top_marker_labels", "PHENOTYPE_GROUP_AGGREGATES",
                       "phenotype_group", "is_phenotype_aggregate", "phenotype_aggregate_members",
                       "phenotype_aggregates_for", "add_phenotype_aggregate_columns")
  for (expr in parse(file.path(repo_root, "app", "shiny_app", "R", "00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% aggregate_names) eval(expr, globalenv())
  }
  missing <- aggregate_names[!vapply(aggregate_names, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(missing)) stop("App definitions missing from 00_globals.R: ", paste(missing, collapse = ", "), call. = FALSE)

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
  if (identical(schema, "5")) {
    target <- file.path(staging, "summary", "embedding_metadata.json")
    write_embedding_metadata(master, target)
    artifacts$embedding_metadata <- file_artifact(target, staging)
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
    # Each feature keeps its stored type (no promotion of a POLYGON / MULTIPOLYGON mix), as the app reads it.
    geom <- sf::st_read(release_geometry_dsn(src), quiet = TRUE, options = "ARRAY_AS_STRING=YES",
                        promote_to_multi = FALSE)
    sf::st_crs(geom) <- NA
    sf::st_write(geom, target, layer = "segmentation", driver = "GPKG",
                 layer_options = "SPATIAL_INDEX=YES", quiet = TRUE)
    geometry_counts[[donor]] <- nrow(geom)
  }
  artifacts$geometry <- tree_artifact(file.path(staging, "geometry"), staging)
  artifacts$geometry$feature_counts <- geometry_counts

  passthrough_files <- c(
    passthrough$required,
    passthrough$optional[file.exists(file.path(source_root, passthrough$optional))]
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
  image_bundle <- release_image_files(image_root)
  for (i in seq_len(nrow(image_bundle))) {
    copy_one(image_bundle$source[[i]], file.path(staging, "images", image_bundle$relative_path[[i]]), mode)
  }
  artifacts$images <- tree_artifact(file.path(staging, "images"), staging,
                                    startup_verify = FALSE)
  expected_images <- sum(grepl("[.]ome[.]tiff?$", image_bundle$relative_path, ignore.case = TRUE))
  artifacts$images$ome_tiff_count <- expected_images

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
    # The OME-TIFF count in images/, which deploy/generate-release-manifest.sh takes as the reviewed
    # EXPECTED_IMAGES and deploy/verify-release.sh checks against the staged files.
    expected_images = expected_images,
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
