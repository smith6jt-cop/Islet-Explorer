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
    signed_dist_islet = seq(-120, 900, length.out = n) + stats::rnorm(n, sd = 3),
    INS_pos = rep(c(1L, 0L, 0L), length.out = n),
    CD8a_pos = rep(c(0L, 1L), length.out = n),
    stringsAsFactors = FALSE
  )
  parquet <- file.path(root, "spatial", "cell_distance.parquet")
  write_test_parquet(parquet, frame)

  env <- make_backend_env(root)
  sys.source(.find_source("R/spatial_helpers.R"), envir = env)
  env$.cell_distance_cache$data <- frame

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
  # Share: each item's cells over every cell in that status and bin, the same in both backends.
  memory_share <- env$.query_cell_distance_composition_memory(
    "marker", c("INS", "CD8a"), c("ND", "Aab+", "T1D"), "prop")
  duck_share <- env$.query_cell_distance_composition_duckdb(
    "marker", c("INS", "CD8a"), c("ND", "Aab+", "T1D"), "prop")
  expect_equal(
    sort_result(duck_share$comp, c("donor_status", "binx", "item")),
    sort_result(memory_share$comp, c("donor_status", "binx", "item")),
    tolerance = 1e-12, ignore_attr = TRUE)
  # The denominator is every cell, retired labels included (they are kept for marker counts).
  inside <- frame[frame$signed_dist >= memory_share$lo & frame$signed_dist < 200, ]
  first <- memory_share$comp[1, ]
  in_bin <- inside[inside$donor_status == first$donor_status & floor(inside$signed_dist / 5) == first$binx, ]
  expect_equal(first$total, nrow(in_bin))
  expect_true(all(memory_share$comp$val <= 1))
  # The distance to the nearest Islet (v9) reads the same in both backends.
  memory_islet <- env$.query_cell_distance_composition_memory(
    "marker", c("INS", "CD8a"), c("ND", "Aab+", "T1D"), "prop", "signed_dist_islet")
  duck_islet <- env$.query_cell_distance_composition_duckdb(
    "marker", c("INS", "CD8a"), c("ND", "Aab+", "T1D"), "prop", "signed_dist_islet")
  expect_equal(duck_islet$lo, memory_islet$lo, tolerance = 1e-12)
  expect_equal(
    sort_result(duck_islet$comp, c("donor_status", "binx", "item")),
    sort_result(memory_islet$comp, c("donor_status", "binx", "item")),
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(memory_islet$comp$val, memory_share$comp$val)))
  hist_islet_mem <- env$.query_cell_distance_hist_memory("marker", "CD8a", c("ND", "Aab+", "T1D"), "prop", 24L,
                                                         "signed_dist_islet")
  hist_islet_duck <- env$.query_cell_distance_hist_duckdb("marker", "CD8a", c("ND", "Aab+", "T1D"), "prop", 24L,
                                                          "signed_dist_islet")
  expect_equal(sort_result(hist_islet_duck$hist, c("donor_status", "binx")),
               sort_result(hist_islet_mem$hist, c("donor_status", "binx")), tolerance = 1e-12, ignore_attr = TRUE)
  expect_error(env$.cd_sql("SELECT signed_dist", "x; DROP TABLE t"))
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

# ---- Data schema 5: the uniform 20-donor release contract (2026-10-22) ----

test_that("schema 5 requires the uniform root's files and the embedding metadata; the frozen schema does not", {
  env <- make_backend_env(tempfile("islet-schema-"))
  expect_identical(env$ISLET_RELEASE_SCHEMAS, c("1", "5"))
  one <- env$islet_required_release_artifacts("1")
  five <- env$islet_required_release_artifacts("5")
  expect_identical(one, env$ISLET_REQUIRED_RELEASE_ARTIFACTS)
  expect_true("phenotype_rules" %in% one)
  expect_false("phenotype_rules" %in% five)
  expect_true(all(setdiff(one, "phenotype_rules") %in% five))
  expect_true(all(c("uniform_release", "endocrine_objects_by_donor", "protein_phenotype_rules",
                    "embedding_metadata") %in% five))
  expect_false(any(c("uniform_release", "embedding_metadata") %in% one))
  # The retired experimental ids are not accepted.
  expect_error(env$islet_required_release_artifacts("4"), "Unsupported data release schema: 4")
})

test_that("a schema-5 release validates under ISLET_DATA_SCHEMA=5 and reads its embedding metadata back", {
  skip_if_not_installed("digest")
  skip_if_not_installed("jsonlite")
  root <- tempfile("islet-schema5-")
  dir.create(file.path(root, "summary"), recursive = TRUE)
  env <- make_backend_env(root)
  artifacts <- list()
  for (name in env$islet_required_release_artifacts("5")) {
    rel <- if (name == "embedding_metadata") env$ISLET_EMBEDDING_METADATA_ARTIFACT else paste0(name, ".txt")
    path <- file.path(root, rel)
    if (name == "embedding_metadata") {
      jsonlite::write_json(list(composition_embedding = list(objects = 1215L, min_core_cells = 10L,
                                                             features = c("Beta", "Alpha"), min_dist = 3),
                                protein_embedding = NULL),
                           path, auto_unbox = TRUE, null = "null")
    } else writeLines(name, path)
    artifacts[[name]] <- list(type = "file", relative_path = rel, sha256 = env$.islet_sha256_file(path),
                              startup_verify = TRUE)
  }
  write_manifest <- function(artifacts) {
    jsonlite::write_json(list(schema_version = "5", artifacts = artifacts), file.path(root, "manifest.json"),
                         auto_unbox = TRUE)
  }
  write_manifest(artifacts)
  withr::with_envvar(c(ISLET_DATA_SCHEMA = "5", ISLET_VERIFY_RELEASE_HASHES = "core"), {
    manifest <- env$islet_validate_release(root)
    expect_identical(manifest$schema_version, "5")
    meta <- env$islet_read_embedding_metadata(manifest, root)
    expect_identical(meta$composition_embedding$objects, 1215L)
    expect_identical(meta$composition_embedding$features, c("Beta", "Alpha"))
    expect_null(meta$protein_embedding)
    # A manifest without the metadata artifact gives no metadata, and fails validation under schema 5.
    without <- artifacts
    without$embedding_metadata <- NULL
    expect_null(env$islet_read_embedding_metadata(list(artifacts = without), root))
    write_manifest(without)
    expect_error(env$islet_validate_release(root), "embedding_metadata")
  })
  write_manifest(artifacts)
  withr::with_envvar(c(ISLET_DATA_SCHEMA = "1", ISLET_VERIFY_RELEASE_HASHES = "core"), {
    expect_error(env$islet_validate_release(root), "schema mismatch")
  })
})

test_that("the production loader attaches a schema-5 release's embedding metadata to prepared()", {
  skip_if_not_installed("jsonlite")
  root <- tempfile("islet-loader-")
  dir.create(file.path(root, "summary"), recursive = TRUE)
  env <- make_backend_env(root)
  sys.source(.find_source("R/data_loading.R"), envir = env)
  jsonlite::write_json(list(composition_embedding = list(objects = 7L),
                            protein_embedding = list(donor_correction = "harmony")),
                       file.path(root, env$ISLET_EMBEDDING_METADATA_ARTIFACT), auto_unbox = TRUE)
  env$islet_read_parquet <- function(path) data.frame(x = 1)
  env$islet_log_event <- function(...) invisible(NULL)
  env$islet_validate_release <- function(root) list(artifacts = list(
    embedding_metadata = list(relative_path = env$ISLET_EMBEDDING_METADATA_ARTIFACT)))
  out <- env$load_prepared_release(root)
  expect_identical(out$composition_embedding$objects, 7L)
  expect_identical(out$protein_embedding$donor_correction, "harmony")
  expect_setequal(names(out), c(names(env$ISLET_SUMMARY_ARTIFACTS), "composition_embedding", "protein_embedding"))
  # A schema-1 manifest adds nothing.
  env$islet_validate_release <- function(root) list(artifacts = list())
  expect_setequal(names(env$load_prepared_release(root)), names(env$ISLET_SUMMARY_ARTIFACTS))
})

test_that("tree digests use byte order, whatever the collation locale", {
  skip_if_not_installed("digest")
  root <- tempfile("islet-tree-order-")
  dir.create(root)
  writeLines("a", file.path(root, "a.csv"))
  writeLines("b", file.path(root, "a_1.csv"))
  env <- make_backend_env(root)
  collate <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", collate), add = TRUE)
  Sys.setlocale("LC_COLLATE", "C")
  in_c <- env$.islet_sha256_tree(root)
  # en_US collation ignores punctuation and would put a_1.csv first; the digest must not change.
  utf <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if(!nzchar(utf), "en_US.UTF-8 collation is not installed")
  expect_identical(order(c("a_1.csv", "a.csv")), c(1L, 2L))
  expect_identical(env$.islet_sha256_tree(root), in_c)
})
