#!/usr/bin/env Rscript
# Validate real Python assembly-test exports through the actual release CSV to
# Parquet conversion and state parity functions. No historical source is edited.
args <- commandArgs(trailingOnly = TRUE)
identity_mode <- length(args) == 2L && identical(args[[1L]], "--identity")
if (length(args) != 1L && !identity_mode) {
  stop("Usage: check_uniform_cell_states.R REAL_ASSEMBLY_TEST_FRAGMENT | --identity REAL_IDENTITY_PARQUET")
}
source_root <- normalizePath(args[[length(args)]], mustWork = TRUE)
source("scripts/verify_release_parity.R")
builder <- parity_builder(getwd())
runtime <- new.env(parent = globalenv())
sys.source("app/shiny_app/R/01_runtime_data_backend.R", envir = runtime)
release_root <- tempfile("real-cell-state-release-")
dir.create(release_root)
con <- runtime$islet_duckdb_connection()
tryCatch({
  expect_rejection <- function(code, expected) {
    error <- tryCatch({ force(code); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(expected, conditionMessage(error), fixed = TRUE))
  }
  if (identity_mode) {
    # Workstation-safe conversion/parity regression on actual upstream rows.
    # Full donor assembly belongs in the allocated HPG job.
    original <- runtime$islet_read_parquet(source_root, columns = c("object_id", "Ki67_pos"))
    parity_cell_state_values(original, "Ki67_pos")
    csv <- file.path(release_root, "state.csv")
    path <- file.path(release_root, "state.parquet")
    data.table::fwrite(original, csv)
    builder$csv_to_parquet(con, csv, path)
    restored <- runtime$islet_read_parquet(path)
    parity_cell_state_values(restored, "Ki67_pos")
    parity <- parity_compare_frames(original, restored, "object_id")
    altered <- restored
    altered$Ki67_pos[[1L]] <- 1 - altered$Ki67_pos[[1L]]
    expect_rejection(parity_compare_frames(original, altered, "object_id"), "numeric_mismatch")
    altered$Ki67_pos[[1L]] <- 2
    expect_rejection(parity_cell_state_values(altered, "Ki67_pos"), "cell_state_not_binary")
    altered$Ki67_pos <- NULL
    expect_rejection(parity_cell_state_values(altered, "Ki67_pos"), "cell_state_field_missing")
    cat(jsonlite::toJSON(list(result = "pass", source = source_root,
      scope = "Actual identity rows through release CSV-to-Parquet conversion; full assembly not run",
      parity = parity, positive = sum(original$Ki67_pos),
      rejected = c("flipped state", "nonbinary state", "missing state")),
      auto_unbox = TRUE, pretty = TRUE), "\n")
  } else {
    builder$csv_to_parquet(con, file.path(source_root, "cell_distance.csv"),
      file.path(release_root, "spatial", "cell_distance.parquet"), character_case_id = TRUE)
    for (file in list.files(file.path(source_root, "donors"), "[.]csv$", full.names = TRUE)) {
      donor <- sub("[.]csv$", "", basename(file))
      builder$csv_to_parquet(con, file, file.path(release_root, "spatial", "tissue",
        paste0("case_id=", donor), paste0(donor, ".parquet")), preserve_blank_islet_name = TRUE)
    }
    stopifnot(file.copy(file.path(source_root, "cells"), release_root, recursive = TRUE))
    result <- parity_cell_states(source_root, release_root, runtime)
    stopifnot(isTRUE(result$applicable), "Ki67_pos" %in% unlist(result$fields))
    path <- file.path(release_root, "spatial", "cell_distance.parquet")
    rows <- runtime$islet_read_parquet(path)
    original <- rows$Ki67_pos[[1L]]
    rows$Ki67_pos[[1L]] <- 1 - original
    unlink(path)
    builder$write_frame_parquet(con, rows, path, "altered_state_distance")
    expect_rejection(parity_cell_states(source_root, release_root, runtime), "numeric_mismatch")
    rows$Ki67_pos[[1L]] <- 2L
    unlink(path)
    builder$write_frame_parquet(con, rows, path, "nonbinary_state_distance")
    expect_rejection(parity_cell_states(source_root, release_root, runtime), "cell_state_not_binary")
    rows$Ki67_pos[[1L]] <- original
    rows$Ki67_pos <- NULL
    unlink(path)
    builder$write_frame_parquet(con, rows, path, "missing_state_distance")
    expect_rejection(parity_cell_states(source_root, release_root, runtime), "Ki67_pos")
    cat(jsonlite::toJSON(list(result = "pass", source = source_root, parity = result,
      rejected = c("flipped state", "nonbinary state", "missing state")),
      auto_unbox = TRUE, pretty = TRUE), "\n")
  }
}, finally = {
  runtime$islet_disconnect_backend()
  unlink(release_root, recursive = TRUE)
})
