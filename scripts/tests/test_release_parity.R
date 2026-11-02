#!/usr/bin/env Rscript
# Run from repository root: Rscript scripts/tests/test_release_parity.R
suppressPackageStartupMessages(library(testthat))
source("scripts/verify_release_parity.R")

test_that("releases reject unresolved endocrine calls and inconsistent aggregate definitions", {
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  root <- tempfile("endocrine-group-contract-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  named <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
  app_members <- function(aggregate, labels) intersect(named, labels)
  write_contract <- function(labels, endocrine) {
    jsonlite::write_json(list(labels = labels, label_groups = list(endocrine = endocrine)),
                         file.path(root, "uniform_release.json"), auto_unbox = FALSE)
  }
  # The v20 failure remains a failure even if someone removes the unnamed cells
  # from the denominator without resolving their identities.
  write_contract(c(named, "EndocrineCluster", "EndocrineSingle", "Islet", "Tissue"),
                  c(named, "EndocrineCluster", "EndocrineSingle"))
  expect_error(builder$validate_release_label_groups(root, app_members), "Unresolved endocrine identities")
  write_contract(c(named, "EndocrineSingle"), named)
  expect_error(builder$validate_release_label_groups(root, app_members), "Unresolved endocrine identities")
  # Islet/Tissue are locations, not endocrine identities. Their unresolved cells
  # must not silently enter the endocrine denominator.
  write_contract(c(named, "Islet", "Tissue"), c(named, "Islet"))
  expect_error(builder$validate_release_label_groups(root, app_members), "Unresolved endocrine identities")
  write_contract(c(named, "Islet", "Tissue"), named[-1])
  expect_error(builder$validate_release_label_groups(root, app_members), "Endocrine group mismatch")
  write_contract(c(named, "Islet", "Tissue"), named)
  result <- builder$validate_release_label_groups(root, app_members)
  expect_true(result$applicable)
  expect_setequal(result$endocrine, named)
  write_contract(named, c(named, "Unknown"))
  expect_error(builder$validate_release_label_groups(root, app_members), "Invalid declared endocrine group")
})

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

test_that("the image bundle adds the sidecars beside a linked TIFF and keeps a complete tree as it is", {
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  src <- tempfile("tiff-source-")
  dir.create(src)
  names <- c("ND_0001.ome.tiff", "ND_0001.ome.xml", "ND_0001.ome.tiff.offsets.json", "ND_0001.offsets.json")
  for (name in names) writeLines(name, file.path(src, name))
  linked <- tempfile("image-root-linked-")
  dir.create(linked)
  file.symlink(file.path(src, "ND_0001.ome.tiff"), file.path(linked, "ND_0001.ome.tiff"))
  on.exit(unlink(c(src, linked), recursive = TRUE), add = TRUE)
  bundle <- builder$release_image_files(linked)
  expect_identical(bundle$relative_path, c("ND_0001.offsets.json", "ND_0001.ome.tiff",
                                           "ND_0001.ome.tiff.offsets.json", "ND_0001.ome.xml"))
  expect_identical(normalizePath(bundle$source), normalizePath(file.path(src, bundle$relative_path)))
  # The same files as the source directory itself, with the same digest as its plain tree.
  expect_identical(builder$release_image_files(src)$relative_path, bundle$relative_path)
  expect_identical(builder$sha256_image_bundle(linked), builder$sha256_image_bundle(src))
  expect_identical(builder$sha256_image_bundle(src), builder$sha256_tree(src))
  # A TIFF without sidecars is released alone.
  bare <- tempfile("image-root-bare-")
  dir.create(bare)
  writeLines("x", file.path(bare, "T1D_0002.ome.tiff"))
  expect_identical(builder$release_image_files(bare)$relative_path, "T1D_0002.ome.tiff")
})

test_that("the passthrough list follows the schema and refuses the wrong schema for a root", {
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  frozen <- tempfile("frozen-root-")
  uniform <- tempfile("uniform-root-")
  dir.create(frozen)
  dir.create(uniform)
  writeLines("{}", file.path(uniform, "uniform_release.json"))
  on.exit(unlink(c(frozen, uniform), recursive = TRUE), add = TRUE)
  one <- builder$release_passthrough_files("1", frozen)
  expect_true("phenotype_rules.csv" %in% one$required)
  expect_length(one$optional, 0)
  five <- builder$release_passthrough_files("5", uniform)
  expect_true(all(c("uniform_release.json", "endocrine_objects_by_donor.csv", "protein_phenotype_rules.csv",
                    "donor_metadata.csv", "adata_ins_root.h5ad") %in% five$required))
  expect_identical(five$optional, "phenotype_rules.csv")
  expect_error(builder$release_passthrough_files("1", uniform), "--schema 5")
  expect_error(builder$release_passthrough_files("5", frozen), "uniform root")
  expect_error(builder$release_passthrough_files("4", frozen), "Unsupported --schema 4")
})

test_that("embedding metadata survives the JSON round trip with its arrays flattened and doubles exact", {
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  runtime <- new.env(parent = globalenv())
  sys.source("app/shiny_app/R/01_runtime_data_backend.R", envir = runtime)
  root <- tempfile("embedding-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  master <- list(composition_embedding = list(objects = 1215L, features = array(c("Beta", "Alpha"), 2),
                                              excluded = array("Endocrine", 1), scaling_mean = array(c(0.39927629123456789, 0), 2),
                                              software = list(scanpy = "1.11.5"), shared = TRUE),
                 protein_embedding = NULL)
  target <- file.path(root, "summary", "embedding_metadata.json")
  builder$write_embedding_metadata(master, target)
  manifest <- list(artifacts = list(embedding_metadata = list(relative_path = "summary/embedding_metadata.json")))
  meta <- runtime$islet_read_embedding_metadata(manifest, root)
  expect_identical(meta$composition_embedding$objects, 1215L)
  expect_identical(meta$composition_embedding$features, c("Beta", "Alpha"))
  expect_identical(meta$composition_embedding$excluded, "Endocrine")
  expect_identical(meta$composition_embedding$scaling_mean, c(0.39927629123456789, 0))
  expect_identical(meta$composition_embedding$software$scanpy, "1.11.5")
  expect_true(meta$composition_embedding$shared)
  expect_null(meta$protein_embedding)
})

test_that("compressed GeoJSON is read through GDAL's gzip file system", {
  skip_if_not_installed("sf")
  builder <- new.env(parent = globalenv())
  sys.source("scripts/build_release_data.R", envir = builder)
  root <- tempfile("gz-geometry-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "0112.geojson.gz")
  con <- gzfile(path, "wt")
  writeLines(paste0('{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"Islet_1"},',
                    '"geometry":{"type":"Polygon","coordinates":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}'), con)
  close(con)
  expect_match(builder$release_geometry_dsn(path), "^/vsigzip/")
  geometry <- sf::st_read(builder$release_geometry_dsn(path), quiet = TRUE, options = "ARRAY_AS_STRING=YES")
  expect_equal(nrow(geometry), 1L)
  expect_identical(geometry$name, "Islet_1")
  plain <- file.path(root, "0112.geojson")
  writeLines("{}", plain)
  expect_identical(builder$release_geometry_dsn(plain), normalizePath(plain))
})
