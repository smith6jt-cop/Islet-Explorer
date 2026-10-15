#!/usr/bin/env Rscript
# Run from repository root: Rscript scripts/tests/test_release_parity.R
suppressPackageStartupMessages(library(testthat))
source("scripts/verify_release_parity.R")

test_that("donor Parquet preserves blank identifiers and numeric missingness", {
  # These are required CI dependencies: missing packages must fail, not skip.
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  runtime <- new.env(parent = globalenv())
  sys.source("app/shiny_app/R/01_runtime_data_backend.R", envir = runtime)
  root <- tempfile("donor-csv-parity-")
  dir.create(root)
  source_path <- file.path(root, "source.csv")
  writeLines(c("X_centroid,Y_centroid,phenotype,cell_region,islet_name",
               "1.25,2,Beta,core,Islet_1", ",3,Beta,peri,",
               '4,,Delta,tissue,""', '5,6,Alpha,core,"NA"'), source_path)
  target <- file.path(root, "case_id=0112", "tissue.parquet")
  con <- runtime$islet_duckdb_connection()
  on.exit({ runtime$islet_disconnect_backend(); unlink(root, recursive = TRUE) }, add = TRUE)
  expected <- data.table::fread(source_path, data.table = FALSE)
  builder$csv_to_parquet(con, source_path, target, preserve_blank_islet_name = TRUE)
  actual <- parity_read_stored_parquet(runtime, target)
  expect_setequal(names(actual), names(expected))
  expect_true("case_id" %in% names(runtime$islet_read_parquet(target)))
  expect_identical(actual$islet_name, expected$islet_name)
  expect_identical(sum(actual$islet_name == "", na.rm = TRUE), 2L)
  expect_identical(is.na(actual$X_centroid), is.na(expected$X_centroid))
  expect_identical(is.na(actual$Y_centroid), is.na(expected$Y_centroid))
  result <- parity_compare_frames(expected, actual, c("phenotype", "cell_region", "islet_name"))
  expect_identical(result$rows, 4L)
  expect_equal(result$max_absolute_difference, 0)
})

test_that("row order and duplicate-key order do not change scientific parity", {
  expected <- data.frame(id = c("a", "a", "b"), category = factor(c("x", "y", "x")),
                          n = c(2L, 3L, 4L), value = c(0.1, NA_real_, Inf))
  actual <- expected[c(3, 2, 1), ]
  actual$category <- as.character(actual$category)
  actual$n <- as.numeric(actual$n)
  result <- parity_compare_frames(expected, actual, "id")
  expect_identical(result$rows, 3L)
  expect_equal(result$max_absolute_difference, 0)
  expected$n <- array(expected$n)
  expected$value <- array(expected$value)
  expect_no_error(parity_compare_frames(expected, actual, "id"))
})

test_that("categorical, count, missingness and infinity corruption fail", {
  expected <- data.frame(id = c("a", "b"), category = c("x", "y"), n = c(1, 2), value = c(NA_real_, Inf))
  altered <- expected; altered$category[1] <- "z"
  expect_error(parity_compare_frames(expected, altered, "id"), "categorical_mismatch")
  altered <- expected; altered$n[1] <- altered$n[1] + 1e-13
  expect_error(parity_compare_frames(expected, altered, "id"), "numeric_mismatch")
  altered <- expected; altered$value[1] <- 0
  expect_error(parity_compare_frames(expected, altered, "id"), "missing_value_mismatch")
  altered <- expected; altered$value[1] <- NaN
  expect_error(parity_compare_frames(expected, altered, "id"), "nan_mismatch")
  altered <- expected; altered$value[2] <- -Inf
  expect_error(parity_compare_frames(expected, altered, "id"), "infinity_mismatch")
})

test_that("numeric tolerance is absolute while numerical keys remain exact", {
  expected <- data.frame(id = 1, value = 0.2)
  altered <- expected; altered$value <- altered$value + 5e-13
  expect_no_error(parity_compare_frames(expected, altered, "id"))
  altered$value <- expected$value + 2e-12
  expect_error(parity_compare_frames(expected, altered, "id"), "numeric_mismatch")
  altered <- expected; altered$id <- altered$id + 1e-13
  expect_error(parity_compare_frames(expected, altered, "id"), "numeric_mismatch")
})

test_that("schema and row multiplicity mismatches cannot pass", {
  frame <- data.frame(id = c("a", "a"), value = c(1, 2))
  expect_error(parity_compare_frames(frame, frame[1, ], "id"), "row_count_mismatch")
  expect_error(parity_compare_frames(frame, data.frame(id = c("a", "a")), "id"), "column_mismatch")
  expect_error(parity_compare_frames(frame, frame, "missing"), "stable_key_missing")
})
