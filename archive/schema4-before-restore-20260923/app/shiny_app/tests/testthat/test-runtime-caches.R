test_that("immutable trajectory tables are shared without cross-cohort contamination", {
  root <- tempfile(); dir.create(root); dir.create(file.path(root, "trajectory"))
  on.exit(unlink(root, recursive = TRUE))
  env <- new.env(parent = globalenv())
  env$data_path <- function(...) file.path(root, ...)
  sys.source(.find_source("R/mod_trajectory_server.R"), envir = env)
  reads <- 0L
  env$islet_read_parquet <- function(path) {
    reads <<- reads + 1L
    data.frame(case_id = 78L, measured_cells = 20L)
  }
  path <- file.path(root, "trajectory", "coordinates.parquet")
  writeLines("fixture", path)
  first <- env$trajectory_read("coordinates.csv")
  expect_identical(first$case_id, "78")
  first$measured_cells <- 999
  expect_identical(env$trajectory_read("coordinates.csv")$measured_cells, 20L)
  expect_identical(reads, 1L)
  writeLines("changed fixture size", path)
  env$trajectory_read("coordinates.csv")
  expect_identical(reads, 2L)
  other <- tempfile(); dir.create(other); dir.create(file.path(other, "trajectory"))
  on.exit(unlink(other, recursive = TRUE), add = TRUE)
  file.copy(path, file.path(other, "trajectory", "coordinates.parquet"))
  env$data_path <- function(...) file.path(other, ...)
  env$trajectory_read("coordinates.csv")
  expect_identical(reads, 3L)
})

test_that("cell browsing cache obeys entry and byte bounds", {
  env <- new.env(parent = globalenv())
  env$data_path <- function(...) file.path(tempdir(), ...)
  sys.source(.find_source("R/drilldown_helpers.R"), envir = env)
  one <- data.frame(object_id = "cell", INS = 1)
  for (i in 1:5) env$.cache_islet_cells(as.character(i), one, max_entries = 2L)
  expect_setequal(ls(env$drilldown_cache), c("4", "5"))
  expect_identical(length(env$.drilldown_cache_order), 2L)
  env$.cache_islet_cells("6", one, max_bytes = 1.5 * as.numeric(object.size(one)))
  expect_identical(ls(env$drilldown_cache), "6")
  env$.cache_islet_cells("too_large", one, max_bytes = 1)
  expect_false(exists("too_large", envir = env$drilldown_cache, inherits = FALSE))
})

test_that("Python CSV flags preserve approved boundaries and measured cells", {
  root <- tempfile(); dir.create(root); dir.create(file.path(root, "cells")); dir.create(file.path(root, "trajectory"))
  on.exit(unlink(root, recursive = TRUE))
  env <- new.env(parent = globalenv()); env$data_path <- function(...) file.path(root, ...)
  env$islet_safe_path_segment <- function(x) TRUE
  sys.source(.find_source("R/drilldown_helpers.R"), envir = env)
  writeLines(c("object_id,marker_measurements_complete,INS", "measured,True,2", "unmeasured,False,"), file.path(root, "cells", "78_a.csv"))
  cells <- env$load_islet_cells("78", "a")
  expect_identical(cells$marker_measurements_complete, c(TRUE, FALSE))
  sys.source(.find_source("R/mod_trajectory_server.R"), envir = env)
  writeLines(c("case_id,eligible,supported", "78,True,False"), file.path(root, "trajectory", "coordinates.csv"))
  coords <- env$trajectory_read("coordinates.csv")
  expect_identical(coords$eligible, TRUE)
  expect_identical(coords$supported, FALSE)
  env$.cache_islet_cells("empty", cells[FALSE, ])
  expect_equal(nrow(env$drilldown_cache$empty), 0)
})

test_that("overlapping finer states highlight independently and preserve exact names", {
  env <- new.env(parent = globalenv()); env$data_path <- function(...) file.path(tempdir(), ...)
  sys.source(.find_source("R/drilldown_helpers.R"), envir = env)
  cells <- data.frame(phenotype_states = c('["Immune: CD44 high", "Immune: LGALS3 high"]',
    '["Immune: CD44 high"]', '[]', NA_character_), phenotype_state = "ignored display text")
  expect_identical(env$islet_state_membership(cells, "Immune: CD44 high"), c(TRUE,TRUE,FALSE,FALSE))
  expect_identical(env$islet_state_membership(cells, "Immune: LGALS3 high"), c(TRUE,FALSE,FALSE,FALSE))
  expect_identical(env$islet_state_membership(cells, "CD44 high"), rep(FALSE,4))
})
