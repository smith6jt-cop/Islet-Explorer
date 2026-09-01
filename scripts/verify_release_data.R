#!/usr/bin/env Rscript
# Verify an Islet Explorer release before promotion.
# Usage: Rscript scripts/verify_release_data.R --root /path/to/release [--mode full|core]

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(args) {
  out <- list(mode = "full")
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    i <- i + 1L
    if (i > length(args)) stop("Missing value for --", key, call. = FALSE)
    out[[gsub("-", "_", key, fixed = TRUE)]] <- args[[i]]
    i <- i + 1L
  }
  out
}

script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(if (length(arg)) sub("^--file=", "", arg[[1]])
                else "scripts/verify_release_data.R", mustWork = TRUE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opts$root)) stop("--root is required", call. = FALSE)
  data_root <<- normalizePath(opts$root, mustWork = TRUE)
  data_path <<- function(...) file.path(data_root, ...)
  Sys.setenv(ISLET_VERIFY_RELEASE_HASHES = opts$mode %||% "full")
  repo_root <- dirname(dirname(script_path()))
  source(file.path(repo_root, "app", "shiny_app", "R", "01_runtime_data_backend.R"),
         local = globalenv())
  manifest <- islet_validate_release(data_root)

  con <- islet_duckdb_connection()
  for (name in names(manifest$artifacts)) {
    artifact <- manifest$artifacts[[name]]
    rel <- as.character(artifact$relative_path %||% "")
    expected <- artifact$row_count %||% NULL
    if (is.null(expected) || !grepl("[.]parquet$", rel)) next
    path <- .islet_artifact_path(data_root, rel)
    actual <- as.numeric(DBI::dbGetQuery(con, paste0(
      "SELECT count(*) AS n FROM ", islet_parquet_relation(path, con)
    ))$n[[1]])
    if (!identical(actual, as.numeric(expected))) {
      stop("Row-count mismatch for ", name, ": expected ", expected,
           ", read ", actual, call. = FALSE)
    }
  }
  islet_disconnect_backend()
  message("Verified release schema ", manifest$schema_version, " at ", data_root,
          " (hash mode: ", opts$mode, ")")
}

if (sys.nframe() == 0L) main()
