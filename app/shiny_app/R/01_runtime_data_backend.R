# Runtime data backend shared by all Shiny sessions in one R worker.
#
# The publication release is an immutable directory selected by
# ISLET_DATA_ROOT. DuckDB is embedded in the R process: it opens no socket and
# persists no database file; it only scans read-only Parquet artifacts.

.islet_backend <- new.env(parent = emptyenv())

islet_runtime_mode <- function() {
  tolower(trimws(Sys.getenv("ISLET_RUNTIME_MODE", unset = "development")))
}

islet_is_production <- function() identical(islet_runtime_mode(), "production")

islet_release_schema <- function() {
  Sys.getenv("ISLET_DATA_SCHEMA", unset = "1")
}

islet_is_expanded_schema <- function(schema = islet_release_schema()) as.character(schema) %in% c("2", "3")

islet_duckdb_connection <- function() {
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop("The production data backend requires the DBI and duckdb R packages.",
         call. = FALSE)
  }
  if (!is.null(.islet_backend$con) && DBI::dbIsValid(.islet_backend$con)) {
    return(.islet_backend$con)
  }

  threads <- suppressWarnings(as.integer(Sys.getenv("ISLET_DUCKDB_THREADS", "2")))
  if (!is.finite(threads) || threads < 1L) threads <- 2L
  memory_limit <- Sys.getenv("ISLET_DUCKDB_MEMORY_LIMIT", unset = "1GB")
  temp_dir <- Sys.getenv("ISLET_DUCKDB_TEMP_DIR", unset = file.path(tempdir(), "islet-duckdb"))
  if (!dir.exists(temp_dir) && !dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create DuckDB temporary directory: ", temp_dir, call. = FALSE)
  }

  .islet_backend$con <- DBI::dbConnect(
    duckdb::duckdb(), dbdir = ":memory:", bigint = "numeric",
    config = list(
      threads = as.character(threads),
      memory_limit = memory_limit,
      temp_directory = normalizePath(temp_dir, mustWork = TRUE),
      preserve_insertion_order = "false"
    )
  )
  .islet_backend$con
}

islet_disconnect_backend <- function() {
  if (!is.null(.islet_backend$con)) {
    try(DBI::dbDisconnect(.islet_backend$con, shutdown = TRUE), silent = TRUE)
  }
  rm(list = ls(.islet_backend, all.names = TRUE), envir = .islet_backend)
  invisible(TRUE)
}

islet_sql_string <- function(value, con = islet_duckdb_connection()) {
  as.character(DBI::dbQuoteString(con, as.character(value)))
}

islet_sql_identifier <- function(value, con = islet_duckdb_connection()) {
  as.character(DBI::dbQuoteIdentifier(con, as.character(value)))
}

islet_parquet_relation <- function(path, con = islet_duckdb_connection()) {
  if (!file.exists(path)) stop("Required Parquet artifact is missing: ", path, call. = FALSE)
  paste0("read_parquet(", islet_sql_string(normalizePath(path, mustWork = TRUE), con), ")")
}

islet_read_parquet <- function(path, columns = NULL) {
  con <- islet_duckdb_connection()
  select_sql <- if (is.null(columns)) {
    "*"
  } else {
    paste(vapply(columns, islet_sql_identifier, character(1), con = con), collapse = ", ")
  }
  DBI::dbGetQuery(con, paste0("SELECT ", select_sql, " FROM ",
                              islet_parquet_relation(path, con)))
}

islet_parquet_columns <- function(path) {
  con <- islet_duckdb_connection()
  desc <- DBI::dbGetQuery(con, paste0("DESCRIBE SELECT * FROM ",
                                      islet_parquet_relation(path, con)))
  as.character(desc$column_name)
}

islet_manifest_path <- function(root = data_root) file.path(root, "manifest.json")

islet_read_manifest <- function(root = data_root) {
  path <- islet_manifest_path(root)
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to validate the data release.", call. = FALSE)
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) stop("Invalid release manifest ", path, ": ",
                             conditionMessage(e), call. = FALSE)
  )
}

ISLET_SUMMARY_ARTIFACTS <- c(
  core_area = "summary/core_area.parquet",
  targets_all = "summary/targets.parquet",
  markers_all = "summary/markers.parquet",
  comp = "summary/islets.parquet",
  combined_pairs = "summary/combined_pairs.parquet"
)

ISLET_REQUIRED_RELEASE_ARTIFACTS <- c(
  names(ISLET_SUMMARY_ARTIFACTS), "cell_distance", "signed_distance_kde",
  "donor_tissue", "geometry", "adata_ins_root", "barlow_islet_features",
  "islet_spatial_lookup", "phenotype_rules", "donor_metadata", "cells", "images"
)

islet_required_artifacts <- function(schema = islet_release_schema()) {
  if (identical(as.character(schema), "3")) return(c(
    names(ISLET_SUMMARY_ARTIFACTS), "geometry", "islet_spatial_lookup",
    "phenotype_rules", "donor_metadata", "images", "cells", "cohort", "trajectory"
  ))
  if (identical(as.character(schema), "2")) return(c(
    names(ISLET_SUMMARY_ARTIFACTS), "geometry", "islet_spatial_lookup",
    "phenotype_rules", "donor_metadata", "images", "cells", "cohort", "trajectory", "research"
  ))
  ISLET_REQUIRED_RELEASE_ARTIFACTS
}

.islet_tree_files <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files[file.info(files)$isdir %in% FALSE]
}

.islet_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to validate release hashes.", call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.islet_sha256_tree <- function(root) {
  files <- .islet_tree_files(root)
  norm_root <- normalizePath(root, mustWork = TRUE)
  rel <- substring(files, nchar(norm_root) + 2L)
  # Expanded schemas are built by Python: use portable byte ordering.
  # Preserve schema 1's historical tree digest contract for frozen releases.
  ord <- if (islet_is_expanded_schema()) order(rel, method = "radix") else order(rel)
  records <- vapply(ord, function(i) {
    paste(rel[[i]], .islet_sha256_file(files[[i]]), sep = "  ")
  }, character(1))
  digest::digest(paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

.islet_artifact_path <- function(root, relative_path) {
  rel <- as.character(relative_path %||% "")
  parts <- strsplit(gsub("\\\\", "/", rel), "/", fixed = TRUE)[[1]]
  if (!nzchar(rel) || startsWith(rel, "/") || any(parts == "..")) {
    stop("Unsafe artifact path in release manifest: ", rel, call. = FALSE)
  }
  file.path(root, rel)
}

islet_verify_release_hashes <- function(manifest, root = data_root, mode = NULL) {
  if (is.null(mode)) {
    mode <- tolower(Sys.getenv(
      "ISLET_VERIFY_RELEASE_HASHES",
      unset = if (islet_is_production() || islet_is_expanded_schema()) "core" else "none"
    ))
  }
  if (identical(mode, "none")) return(invisible(TRUE))
  if (!(mode %in% c("core", "full"))) {
    stop("ISLET_VERIFY_RELEASE_HASHES must be one of: none, core, full", call. = FALSE)
  }

  for (name in names(manifest$artifacts)) {
    artifact <- manifest$artifacts[[name]]
    verify <- isTRUE(artifact$startup_verify) || identical(mode, "full")
    if (!verify) next
    path <- .islet_artifact_path(root, artifact$relative_path)
    if (!file.exists(path) && !dir.exists(path)) {
      stop("Release artifact is missing: ", artifact$relative_path, call. = FALSE)
    }
    actual <- if (identical(artifact$type, "tree")) {
      .islet_sha256_tree(path)
    } else {
      .islet_sha256_file(path)
    }
    expected <- as.character(artifact$sha256 %||% "")
    if (!nzchar(expected) || !identical(tolower(actual), tolower(expected))) {
      stop("Release SHA-256 mismatch for ", name, " (", artifact$relative_path, ")",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

islet_release_available <- function(root = data_root) {
  file.exists(islet_manifest_path(root)) &&
    all(file.exists(file.path(root, unname(ISLET_SUMMARY_ARTIFACTS))))
}

islet_validate_release <- function(root = data_root) {
  manifest <- islet_read_manifest(root)
  if (is.null(manifest)) stop("Release manifest is missing under ", root, call. = FALSE)

  actual_schema <- as.character(manifest$schema_version %||% "")
  expected_schema <- islet_release_schema()
  if (!identical(actual_schema, expected_schema)) {
    stop("Release schema mismatch: expected ", expected_schema,
         ", found ", actual_schema, call. = FALSE)
  }

  artifacts <- manifest$artifacts %||% list()
  required <- islet_required_artifacts(actual_schema)
  if (islet_is_expanded_schema(actual_schema) && islet_is_production() && !isTRUE(manifest$scientific_acceptance)) {
    stop("Expanded candidate has not passed scientific deployment acceptance.", call. = FALSE)
  }
  missing_entries <- setdiff(required, names(artifacts))
  if (length(missing_entries)) {
    stop("Release manifest lacks required artifacts: ",
         paste(missing_entries, collapse = ", "), call. = FALSE)
  }
  missing_paths <- vapply(required, function(name) {
    path <- .islet_artifact_path(root, artifacts[[name]]$relative_path)
    !file.exists(path) && !dir.exists(path)
  }, logical(1))
  if (any(missing_paths)) {
    stop("Release is incomplete; missing: ",
         paste(required[missing_paths], collapse = ", "),
         call. = FALSE)
  }
  if (identical(actual_schema, "3")) {
    files <- c("counts.csv", "endocrine_counts.csv", "phenotype_availability.csv", "marker_panel.json", "count_manifest.json")
    if (!all(file.exists(file.path(root, "trajectory", files))))
      stop("Schema 3 requires trajectory count, availability and marker-panel artifacts.", call. = FALSE)
    if ("research" %in% names(artifacts)) stop("Schema 3 research reports must be outside the app runtime.", call. = FALSE)
  }
  islet_verify_release_hashes(manifest, root)
  manifest
}

# Shared AnnData object cache. It is safe within a Shiny R worker because the
# event loop is single-threaded and the app treats these objects as immutable.
islet_read_h5ad_shared <- function(path) {
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("The anndata R package is required to read trajectory H5AD data.", call. = FALSE)
  }
  norm <- normalizePath(path, mustWork = TRUE)
  key <- paste0("h5ad::", norm)
  if (!is.null(.islet_backend[[key]])) return(.islet_backend[[key]])
  ad <- anndata::read_h5ad(norm)
  .islet_backend[[key]] <- ad
  ad
}
