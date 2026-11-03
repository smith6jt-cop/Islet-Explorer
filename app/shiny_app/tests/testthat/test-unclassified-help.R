.unclassified_help_env <- function() {
  env <- new.env(parent = globalenv())
  wanted <- c("UNCLASSIFIED_KEYS", "UNCLASSIFIED_PARENT", "UNCLASSIFIED_RAW_TO_KEY",
              "unclassified_collapse", "ISLET_EXCLUDED_DONORS", "islet_canonical_donor_id",
              "islet_excluded_donors_in", "PHENOTYPE_DISPLAY_NAMES", "MACROPHAGE_LABEL_PATTERN",
              "phenotype_text")
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% wanted) eval(expr, env)
  }
  sys.source(.find_source("R/01_runtime_data_backend.R"), env)
  sys.source(.find_source("R/unclassified_help.R"), env)
  env
}

.unclassified_help_fixture <- function(env) {
  root <- tempfile("unclassified-census-")
  dir.create(root)
  donors <- list(
    `001` = c("Tissue", "Tissue", "Unclassified", "Conflict B/T", "Beta cell"),
    `002` = c("Hematopoietic CD45+", "Hematopoietic CD45-", "T cell", "Doublet suspected", "Conflict B/T"),
    `003` = c("Tissue", "Doublet suspected", "Unclassified", "Unclassified"))
  con <- env$islet_duckdb_connection()
  for (id in names(donors)) {
    path <- file.path(root, "spatial", "tissue", paste0("case_id=", id), "cells.parquet")
    dir.create(dirname(path), recursive = TRUE)
    DBI::dbWriteTable(con, "fixture_cells", data.frame(phenotype = donors[[id]],
                      cell_region = rep(c("core", "peri", "tissue"), length.out = length(donors[[id]]))),
                      overwrite = TRUE)
    DBI::dbExecute(con, paste0("COPY fixture_cells TO ", env$islet_sql_string(path, con), " (FORMAT PARQUET)"))
  }
  metadata <- data.frame(case_id = c("1", "0002", "3"),
                         donor_type = c("No Diabetes", "Autoantibody Positive", "Type 1 Diabetes"))
  utils::write.csv(metadata, file.path(root, "donor_metadata.csv"), row.names = FALSE)
  root
}

test_that("whole-cohort counts join donor status and keep parent distinct from literal Unclassified", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  env <- .unclassified_help_env()
  root <- .unclassified_help_fixture(env)
  on.exit({env$islet_disconnect_backend(); unlink(root, recursive = TRUE)}, add = TRUE)
  result <- env$unclassified_census(root)
  expect_true(result$available)
  expect_equal(result$tissue_cells, 14)
  expect_equal(result$donors, 3)
  expect_identical(result$table$key, c(env$UNCLASSIFIED_PARENT, unname(env$UNCLASSIFIED_KEYS)))
  expect_equal(unlist(result$table[1, c("ND", "Aab+", "T1D", "Total")], use.names = FALSE), c(3, 4, 4, 11))
  expect_equal(result$table$Total[-1], c(3, 3, 3, 2))
  expect_equal(result$table$Total[[1]], sum(result$table$Total[-1]))
  expect_false(any(grepl("^Conflict", result$raw_counts$phenotype)))
  expect_equal(sum(result$raw_counts$cells[result$raw_counts$phenotype == "Unclassified"]), 3)

  html <- as.character(env$unclassified_help_ui(root))
  expect_match(html, "Global cell counts by donor status", fixed = TRUE)
  expect_match(html, "not the current card", fixed = TRUE)
  expect_match(html, "not proof of a physical doublet", fixed = TRUE)
  expect_match(html, "at least three lineage or compartment groups", fixed = TRUE)
  expect_match(html, "that literal source label is only this subclass", fixed = TRUE)
  expect_match(html, "scope=\"col\"", fixed = TRUE)
})

test_that("census cache reuses a snapshot and invalidates when donor metadata changes", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  env <- .unclassified_help_env()
  root <- .unclassified_help_fixture(env)
  on.exit({env$islet_disconnect_backend(); unlink(root, recursive = TRUE)}, add = TRUE)
  reads <- 0L
  original <- env$unclassified_census_read
  env$unclassified_census_read <- function(sources) {reads <<- reads + 1L; original(sources)}
  first <- env$unclassified_census(root)
  expect_identical(env$unclassified_census(root), first)
  expect_equal(reads, 1L)
  metadata_path <- file.path(root, "donor_metadata.csv")
  metadata <- utils::read.csv(metadata_path)
  metadata$donor_type[[1]] <- "T1D"
  utils::write.csv(metadata, metadata_path, row.names = FALSE)
  changed <- env$unclassified_census(root)
  expect_equal(reads, 2L)
  expect_equal(unlist(changed$table[1, c("ND", "Aab+", "T1D", "Total")], use.names = FALSE), c(0, 4, 7, 11))
})

test_that("ambiguous or missing donor evidence never becomes a zero census", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  env <- .unclassified_help_env()
  root <- .unclassified_help_fixture(env)
  on.exit({env$islet_disconnect_backend(); unlink(root, recursive = TRUE)}, add = TRUE)
  path <- file.path(root, "donor_metadata.csv")
  metadata <- utils::read.csv(path)
  duplicate <- rbind(metadata, transform(metadata[1, ], case_id = "0001"))
  utils::write.csv(duplicate, path, row.names = FALSE)
  result <- env$unclassified_census(root)
  expect_false(result$available)
  expect_null(result$table)
  expect_match(result$reason, "unique", fixed = TRUE)
  metadata$donor_type[[1]] <- "Unrecorded status"
  utils::write.csv(metadata, path, row.names = FALSE)
  expect_false(env$unclassified_census(root)$available)
  unlink(path)
  expect_false(env$unclassified_census(root)$available)
  expect_match(as.character(env$unclassified_help_ui(root)), "Global counts are unavailable", fixed = TRUE)
})
