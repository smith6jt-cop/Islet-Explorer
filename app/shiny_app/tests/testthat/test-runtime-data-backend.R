make_backend_env <- function(root) {
  env <- new.env(parent = globalenv())
  env$data_root <- root
  env$data_path <- function(...) file.path(env$data_root, ...)
  env$PHENOTYPE_RETIRED <- c("Retired")
  env$PHENOTYPE_COLORS <- c("Beta cell" = "#000000", "Macrophage" = "#ffffff")
  env$phenotype_keep <- function(x) setdiff(x, env$PHENOTYPE_RETIRED)
  env$phenotype_filter_rows <- function(df) {
    df[!(as.character(df$phenotype) %in% env$PHENOTYPE_RETIRED), , drop = FALSE]
  }
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% c("ISLET_EXCLUDED_DONORS", "islet_canonical_donor_id",
                                       "islet_excluded_donors_in")) {
      eval(expr, envir = env)
    }
  }
  sys.source(.find_source("R/01_runtime_data_backend.R"), envir = env)
  env
}

write_test_parquet <- function(path, frame) {
  con <- DBI::dbConnect(duckdb::duckdb(), bigint = "numeric")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  duckdb::duckdb_register(con, "test_frame", frame)
  on.exit(try(duckdb::duckdb_unregister(con, "test_frame"), silent = TRUE), add = TRUE)
  DBI::dbExecute(con, paste0(
    "COPY test_frame TO ", as.character(DBI::dbQuoteString(con, path)),
    " (FORMAT PARQUET, COMPRESSION ZSTD)"
  ))
}

sort_result <- function(df, keys) {
  rownames(df) <- NULL
  df[do.call(order, df[keys]), , drop = FALSE]
}

test_that("DuckDB cell-distance queries preserve in-memory results", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("duckdb")

  root <- tempfile("islet-backend-")
  dir.create(file.path(root, "spatial"), recursive = TRUE)
  set.seed(42)
  n <- 1200L
  frame <- data.frame(
    case_id = rep(c("0112", "6505", "6533"), length.out = n),
    donor_status = rep(c("ND", "Aab+", "T1D"), each = n / 3),
    cell_region = rep(c("core", "peri", "tissue"), length.out = n),
    phenotype = rep(c("Beta cell", "Macrophage", "Retired", "Unknown"),
                    length.out = n),
    signed_dist = seq(-180, 720, length.out = n) + stats::rnorm(n, sd = 2),
    INS_pos = rep(c(1L, 0L, 0L), length.out = n),
    CD8a_pos = rep(c(0L, 1L), length.out = n),
    stringsAsFactors = FALSE
  )
  parquet <- file.path(root, "spatial", "cell_distance.parquet")
  write_test_parquet(parquet, frame)

  env <- make_backend_env(root)
  sys.source(.find_source("R/spatial_helpers.R"), envir = env)
  env$.cell_distance_cache$data <- env$phenotype_filter_rows(frame)

  memory_hist <- env$.query_cell_distance_hist_memory(
    "marker", "CD8a", c("ND", "Aab+", "T1D"), "prop", 24L
  )
  duck_hist <- env$.query_cell_distance_hist_duckdb(
    "marker", "CD8a", c("ND", "Aab+", "T1D"), "prop", 24L
  )
  expect_equal(duck_hist$lo, memory_hist$lo, tolerance = 1e-12)
  expect_equal(duck_hist$hi, memory_hist$hi, tolerance = 1e-12)
  expect_equal(
    sort_result(duck_hist$hist, c("donor_status", "binx")),
    sort_result(memory_hist$hist, c("donor_status", "binx")),
    tolerance = 1e-12, ignore_attr = TRUE
  )

  memory_comp <- env$.query_cell_distance_composition_memory(
    "phenotype", c("Beta cell", "Macrophage"), c("ND", "Aab+", "T1D"),
    "count"
  )
  duck_comp <- env$.query_cell_distance_composition_duckdb(
    "phenotype", c("Beta cell", "Macrophage"), c("ND", "Aab+", "T1D"),
    "count"
  )
  # Both views stop 200 um outside the border; composition bars are 5 um wide (owner, 2026-10-20).
  expect_equal(memory_hist$hi, 200)
  expect_equal(memory_comp$hi, 200)
  expect_equal(memory_comp$bw, 5)
  expect_true(max(memory_comp$comp$x) < 200)
  expect_equal(duck_comp$lo, memory_comp$lo, tolerance = 1e-12)
  expect_equal(duck_comp$hi, memory_comp$hi, tolerance = 1e-12)
  expect_equal(
    sort_result(duck_comp$comp, c("donor_status", "binx", "item")),
    sort_result(memory_comp$comp, c("donor_status", "binx", "item")),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  env$islet_disconnect_backend()
})

test_that("manifest hash verification detects a modified artifact", {
  skip_if_not_installed("digest")
  root <- tempfile("islet-manifest-")
  dir.create(root)
  path <- file.path(root, "example.txt")
  writeLines("original", path)
  env <- make_backend_env(root)
  manifest <- list(artifacts = list(example = list(
    type = "file", relative_path = "example.txt",
    sha256 = env$.islet_sha256_file(path), startup_verify = TRUE
  )))
  expect_true(env$islet_verify_release_hashes(manifest, root, mode = "core"))
  writeLines("changed", path)
  expect_error(env$islet_verify_release_hashes(manifest, root, mode = "core"),
               "SHA-256 mismatch")
})

test_that("manifest artifact paths cannot escape the release root", {
  root <- tempfile("islet-manifest-path-")
  dir.create(root)
  env <- make_backend_env(root)
  expect_error(env$.islet_artifact_path(root, "../outside.parquet"),
               "Unsafe artifact path")
  expect_error(env$.islet_artifact_path(root, "/absolute/path.parquet"),
               "Unsafe artifact path")
})

test_that("development uses the selected source root even when a manifest exists", {
  root <- tempfile("islet-development-")
  dir.create(root)
  writeLines("historical manifest", file.path(root, "manifest.json"))
  env <- make_backend_env(root)
  sys.source(.find_source("R/data_loading.R"), envir = env)
  env$islet_is_production <- function() FALSE
  env$load_prepared_release <- function(...) stop("development selected release")
  env$load_master_auto <- function(h5ad_path, excel_path) {
    expect_identical(h5ad_path, file.path(root, "islet_explorer.h5ad"))
    expect_identical(excel_path, file.path(root, "master_results.xlsx"))
    list(source = "development")
  }
  env$prep_data <- identity
  env$audit_na <- function(...) invisible(NULL)

  expect_identical(env$load_shared_prepared_data(root)$source, "development")
})

test_that("excluded donor IDs match in every spelling and other donors pass", {
  env <- make_backend_env(tempfile("islet-excluded-"))
  expect_identical(env$ISLET_EXCLUDED_DONORS, "78")
  expect_identical(env$islet_canonical_donor_id(c("78", "0078", "HPAP-078", "HDL078", "78.0", "112", "6232")),
                   c("78", "78", "78", "78", "78", "112", "6232"))
  expect_identical(env$islet_excluded_donors_in(c("112", "0078", NA, "6232", "0078")), "0078")
  expect_length(env$islet_excluded_donors_in(c("112", "780", "1078", "6232")), 0)
})

test_that("a data root that still carries an excluded donor is refused, not filtered", {
  root <- tempfile("islet-excluded-root-")
  dir.create(root)
  env <- make_backend_env(root)
  sys.source(.find_source("R/data_loading.R"), envir = env)
  env$islet_is_production <- function() FALSE
  env$audit_na <- function(...) invisible(NULL)
  env$prep_data <- identity
  env$load_master_auto <- function(h5ad_path, excel_path) {
    list(comp = data.frame(`Case ID` = c("112", "78"), check.names = FALSE),
         targets_all = data.frame(`Case ID` = "112", check.names = FALSE))
  }
  expect_error(env$load_shared_prepared_data(root), "excluded donor\\(s\\) 78")

  env$.islet_backend <- new.env(parent = emptyenv())
  env$load_master_auto <- function(h5ad_path, excel_path) {
    list(comp = data.frame(`Case ID` = c("112", "6232"), check.names = FALSE))
  }
  expect_identical(env$load_shared_prepared_data(root)$comp$`Case ID`, c("112", "6232"))
})

test_that("GeoPackage loading applies an indexed bounding-box filter", {
  skip_if_not_installed("sf")
  root <- tempfile("islet-geometry-")
  dir.create(file.path(root, "geometry"), recursive = TRUE)
  square <- function(x, y) {
    sf::st_polygon(list(matrix(c(x, y, x + 2, y, x + 2, y + 2,
                                 x, y + 2, x, y), ncol = 2, byrow = TRUE)))
  }
  geom <- sf::st_sf(
    id = c("near", "far"),
    classification = c('{"name":"Islet"}', '{"name":"Nerve"}'),
    geometry = sf::st_sfc(square(0, 0), square(100, 100))
  )
  gpkg <- file.path(root, "geometry", "0001.gpkg")
  sf::st_write(geom, gpkg, layer = "segmentation", driver = "GPKG",
               layer_options = "SPATIAL_INDEX=YES", quiet = TRUE)

  env <- make_backend_env(root)
  env$SF_AVAILABLE <- TRUE
  env$PIXEL_SIZE_UM <- 1
  sys.source(.find_source("R/segmentation_helpers.R"), envir = env)
  cropped <- env$load_case_geojson("1", centroid_x_um = 1,
                                   centroid_y_um = 1, buffer_um = 5)
  expect_equal(nrow(cropped), 1L)
  expect_equal(cropped$id, "near")
})

test_that("development segmentation reads compressed GeoJSON for padded donor IDs", {
  skip_if_not_installed("sf")
  root <- tempfile("islet-compressed-geometry-")
  dir.create(file.path(root, "gson"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  con <- gzfile(file.path(root, "gson", "0078.geojson.gz"), "wt")
  writeLines('{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"name":"Islet_1234"},"geometry":{"type":"Polygon","coordinates":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}', con)
  close(con)
  env <- make_backend_env(root)
  env$islet_is_production <- function() FALSE
  env$SF_AVAILABLE <- TRUE
  sys.source(.find_source("R/segmentation_helpers.R"), envir = env)
  geometry <- env$load_case_geojson("78")
  expect_equal(nrow(geometry), 1L)
  expect_identical(geometry$name, "Islet_1234")
  expect_true(is.na(sf::st_crs(geometry)))
})
