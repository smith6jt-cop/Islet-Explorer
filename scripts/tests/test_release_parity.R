#!/usr/bin/env Rscript
# Run from repository root: Rscript scripts/tests/test_release_parity.R
suppressPackageStartupMessages(library(testthat))
source("scripts/verify_release_parity.R")

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
