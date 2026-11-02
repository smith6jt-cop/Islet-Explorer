#!/usr/bin/env Rscript
# Real-data acceptance, without printing scientific rows, filenames or exceptions.
# Rscript scripts/verify_release_parity.R --source-root data --image-root PATH \
#   --release-root PATH --output-prefix /restricted/evidence/release-parity
# Omit --release-root for a source preparation report (never a parity pass).
# Set RETICULATE_PYTHON to the approved interpreter before running this script.

`%||%` <- function(x, y) if (is.null(x)) y else x

PARITY_SUMMARY_KEYS <- list(
  core_area = c("Case ID", "Donor Status", "islet_key"),
  targets_all = c("Case ID", "Donor Status", "islet_key", "type", "class"),
  markers_all = c("Case ID", "Donor Status", "islet_key", "region_type", "marker"),
  comp = c("Case ID", "Donor Status", "islet_key"),
  combined_pairs = "stem"
)

parity_fail <- function(code) stop(structure(
  list(message = code, call = NULL), class = c("parity_failure", "error", "condition")
))

parity_column <- function(value) {
  if (is.factor(value)) return(as.character(value))
  # The H5AD bridge supplies one-dimensional R arrays. Parquet supplies plain
  # vectors; remove that storage-only attribute without changing any values.
  if (is.array(value)) {
    if (length(dim(value)) > 1L && any(dim(value)[-1L] != 1L)) {
      parity_fail("unsupported_array_column")
    }
    return(as.vector(value))
  }
  unname(value)
}

parity_compare_vectors <- function(expected, actual, exact = FALSE, tolerance = 1e-12) {
  expected <- parity_column(expected)
  actual <- parity_column(actual)
  if (length(expected) != length(actual)) parity_fail("length_mismatch")
  if (!identical(is.na(expected), is.na(actual))) parity_fail("missing_value_mismatch")
  numeric_pair <- is.numeric(expected) && is.numeric(actual)
  if (numeric_pair) {
    if (!identical(is.nan(expected), is.nan(actual))) parity_fail("nan_mismatch")
    keep <- !is.na(expected)
    expected <- as.numeric(expected[keep]); actual <- as.numeric(actual[keep])
    if (!identical(is.finite(expected), is.finite(actual))) parity_fail("nonfinite_mismatch")
    finite <- is.finite(expected)
    if (!identical(expected[!finite], actual[!finite])) parity_fail("infinity_mismatch")
    difference <- abs(expected[finite] - actual[finite])
    if (any(difference > if (exact) 0 else tolerance)) parity_fail("numeric_mismatch")
    return(if (length(difference)) max(difference) else 0)
  }
  if (!identical(unname(expected), unname(actual))) parity_fail("categorical_mismatch")
  0
}

parity_compare_frames <- function(expected, actual, keys, tolerance = 1e-12) {
  if (!is.data.frame(expected) || !is.data.frame(actual)) parity_fail("table_missing")
  if (!setequal(names(expected), names(actual))) parity_fail("column_mismatch")
  if (nrow(expected) != nrow(actual)) parity_fail("row_count_mismatch")
  if (!all(keys %in% names(expected))) parity_fail("stable_key_missing")
  actual <- actual[, names(expected), drop = FALSE]
  # All remaining columns break ties, preserving duplicate-row multiplicities.
  # Factor storage is an implementation detail; labels must match exactly.
  order_columns <- unique(c(keys, names(expected)))
  canonical <- function(frame) {
    frame[] <- lapply(frame, parity_column)
    if (!all(vapply(frame, is.atomic, logical(1)))) parity_fail("unsupported_column_type")
    frame[do.call(order, c(unname(frame[order_columns]), list(na.last = TRUE, method = "radix"))),
          , drop = FALSE]
  }
  expected <- canonical(expected); actual <- canonical(actual)
  rownames(expected) <- NULL; rownames(actual) <- NULL
  count_columns <- grep(paste0(
    "(^n$|^n_|^count$|count_|_count$|^cells_total$|^total_cells_|_pos$|_single$|",
    "^(Multi_Pos|Triple_Neg|Ins_any|Glu_any|Stt_any|binx)$)"), names(expected), value = TRUE)
  exact_columns <- unique(c(keys, count_columns))
  differences <- vapply(names(expected), function(column) {
    parity_compare_vectors(expected[[column]], actual[[column]],
                            exact = column %in% exact_columns, tolerance = tolerance)
  }, numeric(1))
  list(rows = nrow(expected), columns = ncol(expected), keys = as.list(keys),
       max_absolute_difference = if (length(differences)) max(differences) else 0)
}

# The builder's own file lists and digests (passthrough files by schema, the image bundle with its
# resolved sidecars, the JSON form of the embedding metadata), so the check cannot drift from the build.
parity_builder <- function(repo) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(repo, "scripts", "build_release_data.R"), envir = env)
  env
}

parity_runtime <- function(repo, root) {
  env <- new.env(parent = globalenv())
  env$data_root <- root
  env$data_path <- function(...) file.path(root, ...)
  env$master_path <- env$data_path("master_results.xlsx")
  env$h5ad_path <- env$data_path("islet_explorer.h5ad")
  env$project_root <- repo
  for (file in c("00_audit.R", "utils_safe_join.R", "01_runtime_data_backend.R",
                  "data_loading.R", "drilldown_helpers.R")) {
    sys.source(file.path(repo, "app", "shiny_app", "R", file), envir = env)
  }
  # Read only the canonical phenotype policy, the pixel constant and the definitions prep_data() needs for the
  # phenotype group aggregates (app commit 82ce23e), without running the UI/credential/package startup code in
  # 00_globals.R.
  policy_names <- c("PHENOTYPE_RETIRED", "phenotype_keep", "phenotype_is_kept",
                    "phenotype_filter_rows", "PIXEL_SIZE_UM",
                    "MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN", "is_macrophage_label",
                    "root_phenotype_labels", "root_has_marker_profile_labels", "TOP_MARKER_BASE",
                    "CD45_RULE_LABEL", "CD45_TOP_MARKER_LABEL", "root_has_top_marker_labels",
                    "root_cd45_is_rule_label", "top_marker_labels", "PHENOTYPE_GROUP_AGGREGATES",
                    "phenotype_group", "is_phenotype_aggregate", "phenotype_aggregate_members",
                    "phenotype_aggregates_for", "add_phenotype_aggregate_columns")
  expressions <- parse(file.path(repo, "app", "shiny_app", "R", "00_globals.R"))
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% policy_names) eval(expr, env)
  }
  if (!all(vapply(policy_names, exists, logical(1), envir = env, inherits = FALSE))) {
    parity_fail("runtime_policy_missing")
  }
  env$SF_AVAILABLE <- TRUE
  sys.source(file.path(repo, "app", "shiny_app", "R", "spatial_helpers.R"), envir = env)
  env
}

parity_hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

parity_read_stored_parquet <- function(env, path) {
  # Compare stored CSV-derived columns, excluding case_id synthesized by
  # DuckDB from the Hive-style directory name. Runtime queries are tested
  # separately through the unchanged application readers.
  con <- env$islet_duckdb_connection()
  DBI::dbGetQuery(con, paste0("SELECT * FROM read_parquet(",
    env$islet_sql_string(normalizePath(path, mustWork = TRUE), con),
    ", hive_partitioning=false)"))
}

parity_tree <- function(root, env) env$.islet_sha256_tree(normalizePath(root, mustWork = TRUE))

parity_query_result <- function(expected, actual, kind) {
  if (is.null(expected) || is.null(actual)) {
    if (!is.null(expected) || !is.null(actual)) parity_fail("query_availability_mismatch")
    return(list(rows = 0L, empty = TRUE))
  }
  keys <- if (kind == "hist") c("donor_status", "binx") else c("donor_status", "binx", "item")
  result <- parity_compare_frames(expected[[kind]], actual[[kind]], keys)
  for (field in intersect(c("lo", "hi", "bw"), names(expected))) {
    parity_compare_vectors(expected[[field]], actual[[field]])
  }
  if (!identical(expected$norm, actual$norm)) parity_fail("query_normalization_mismatch")
  result
}

parity_run <- function(options, repo) {
  source_root <- normalizePath(options$source_root %||% file.path(repo, "data"), mustWork = TRUE)
  release_root <- if (!is.null(options$release_root)) {
    normalizePath(options$release_root, mustWork = TRUE)
  } else NULL
  image_root <- normalizePath(options$image_root %||%
    file.path(repo, "app", "shiny_app", "www", "local_images"), mustWork = TRUE)
  report <- list(schema_version = 1L, created_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                 mode = if (is.null(release_root)) "source_only" else "release_comparison",
                 numeric_absolute_tolerance = 1e-12, categorical_count_na_parity = "exact",
                 identity = list(), summaries = list(), checks = list())
  record <- function(name, expr) {
    result <- tryCatch(suppressWarnings(suppressMessages(force(expr))), error = identity)
    failed <- inherits(result, "error")
    report$checks[[length(report$checks) + 1L]] <<- list(
      check = name, status = if (failed) "failed" else "passed",
      details = if (failed) list(error_code = if (inherits(result, "parity_failure")) {
        conditionMessage(result)
      } else "check_execution_failed") else result)
    if (failed) NULL else result
  }
  report$identity$git_commit <- unname(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"),
                                             stdout = TRUE, stderr = FALSE)[[1]])
  report$identity$git_dirty <- length(system2("git", c("-C", shQuote(repo), "status", "--porcelain"),
                                              stdout = TRUE, stderr = FALSE)) > 0L
  report$identity$canonical_h5ad_sha256 <- parity_hash(file.path(source_root, "islet_explorer.h5ad"))
  report$identity$donor_metadata_sha256 <- parity_hash(file.path(source_root, "donor_metadata.csv"))
  source_env <- parity_runtime(repo, source_root)
  on.exit(source_env$islet_disconnect_backend(), add = TRUE)
  prepared <- record("canonical_h5ad_preparation", {
    master <- source_env$load_master_h5ad(source_env$h5ad_path)
    if (is.null(master)) parity_fail("h5ad_preparation_failed")
    source_env$prep_data(master)
  })
  # Prepared scientific frames are never serialized to evidence.
  if (!is.null(prepared)) {
    report$checks[[length(report$checks)]]$details <- list(tables = length(prepared))
    report$summaries <- lapply(prepared, function(frame) list(
      rows = nrow(frame), fields = as.list(names(frame)),
      column_types = as.list(vapply(frame, function(x) paste(class(x), collapse = "/"), character(1)))))
  }
  if (is.null(release_root) || is.null(prepared)) {
    report$status <- if (is.null(prepared)) "failed" else "source_prepared_release_pending"
    return(report)
  }

  release_env <- parity_runtime(repo, release_root)
  on.exit(release_env$islet_disconnect_backend(), add = TRUE)
  report$identity$release_manifest_sha256 <- parity_hash(file.path(release_root, "manifest.json"))
  manifest <- record("release_manifest_and_full_hashes", {
    value <- release_env$islet_read_manifest(release_root)
    # The runtime validates a release against ISLET_DATA_SCHEMA; here the release's declared schema is the
    # one under test (the builder has already refused the wrong schema for the root, and the runtime
    # refuses an id it does not know).
    Sys.setenv(ISLET_DATA_SCHEMA = as.character(value$schema_version %||% ""))
    release_env$islet_validate_release(release_root)
    release_env$islet_verify_release_hashes(value, release_root, mode = "full")
    value
  })
  if (!is.null(manifest)) {
    report$checks[[length(report$checks)]]$details <- list(artifacts = length(manifest$artifacts))
    report$identity$release_git_commit <- manifest$git_commit
    report$identity$release_schema <- manifest$schema_version
    report$identity$canonical_sources_manifest_sha256 <- digest::digest(
      jsonlite::toJSON(manifest$canonical_sources, auto_unbox = TRUE, digits = 16),
      algo = "sha256", serialize = FALSE)
    record("release_source_identity", {
      if (!identical(manifest$git_commit, report$identity$git_commit)) parity_fail("source_commit_mismatch")
      if (!identical(manifest$git_dirty, FALSE) || isTRUE(report$identity$git_dirty)) {
        parity_fail("uncommitted_source_changes")
      }
      list(source_commit_matches = TRUE, source_worktree_clean = TRUE)
    })
    record("canonical_source_hashes", {
      for (item in manifest$canonical_sources) {
        path <- release_env$.islet_artifact_path(source_root, item$path)
        if (!identical(parity_hash(path), item$sha256)) parity_fail("source_hash_mismatch")
      }
      list(files = length(manifest$canonical_sources))
    })
  }
  for (name in names(PARITY_SUMMARY_KEYS)) {
    record(paste0("summary_", name), {
      actual <- release_env$islet_read_parquet(file.path(release_root,
        release_env$ISLET_SUMMARY_ARTIFACTS[[name]]))
      parity_compare_frames(prepared[[name]], actual, PARITY_SUMMARY_KEYS[[name]])
    })
  }
  builder <- parity_builder(repo)
  schema <- as.character(manifest$schema_version %||% "")
  record("passthrough_files_and_trees", {
    if (is.null(manifest)) parity_fail("verified_manifest_required")
    spec <- builder$release_passthrough_files(schema, source_root)
    files <- c(spec$required, spec$optional[file.exists(file.path(source_root, spec$optional))])
    for (name in files) {
      if (!identical(parity_hash(file.path(source_root, name)),
                     parity_hash(file.path(release_root, name)))) parity_fail("passthrough_hash_mismatch")
    }
    # The release trees were read in full above. Compare their verified digests
    # with the current source trees without re-reading the 100+ GiB image copy.
    # The image digest is taken over the bundle the builder materializes: the
    # source tree plus the sidecars resolved beside each linked TIFF.
    if (!identical(parity_tree(file.path(source_root, "cells"), source_env),
                   manifest$artifacts$cells$sha256)) parity_fail("passthrough_tree_mismatch")
    if (!identical(builder$sha256_image_bundle(image_root), manifest$artifacts$images$sha256)) {
      parity_fail("passthrough_tree_mismatch")
    }
    list(files = length(files), trees = 2L)
  })
  record("embedding_metadata", {
    if (is.null(manifest)) parity_fail("verified_manifest_required")
    # A schema-5 release carries the H5AD's embedding records; the frozen schema has none to carry.
    if (!identical(schema, "5")) return(list(applicable = FALSE))
    actual <- release_env$islet_read_embedding_metadata(manifest, release_root)
    if (is.null(actual)) parity_fail("embedding_metadata_missing")
    expected <- list(composition_embedding = prepared$composition_embedding,
                     protein_embedding = prepared$protein_embedding)
    canonical <- function(x) jsonlite::fromJSON(jsonlite::toJSON(
      builder$release_plain_list(x), auto_unbox = TRUE, null = "null", na = "null", digits = I(17)),
      simplifyVector = TRUE)
    if (!identical(canonical(expected), canonical(actual))) parity_fail("embedding_metadata_mismatch")
    list(applicable = TRUE, records = sum(!vapply(expected, is.null, logical(1))))
  })
  record("signed_distance_kde_table", {
    expected <- data.table::fread(file.path(source_root, "signed_distance_kde.csv"),
                                  colClasses = c(case_id = "character"), data.table = FALSE)
    actual <- release_env$islet_read_parquet(file.path(release_root, "spatial", "signed_distance_kde.parquet"))
    parity_compare_frames(expected, actual, c("case_id", "donor_status", "phenotype", "cell_region"))
  })
  record("donor_tissue_tables", {
    files <- sort(list.files(file.path(source_root, "donors"), pattern = "[.]csv$", full.names = TRUE))
    if (!length(files)) parity_fail("donor_tissue_missing")
    rows <- 0L
    for (file in files) {
      donor <- sub("[.]csv$", "", basename(file))
      expected <- data.table::fread(file, data.table = FALSE)
      actual <- parity_read_stored_parquet(release_env, file.path(release_root, "spatial", "tissue",
        paste0("case_id=", donor), paste0(donor, ".parquet")))
      parity_compare_frames(expected, actual, c("phenotype", "cell_region", "islet_name"))
      rows <- rows + nrow(expected)
    }
    list(files = length(files), rows = rows)
  })

  cells <- record("spatial_source_load", {
    frame <- source_env$load_cell_distance()
    if (is.null(frame)) parity_fail("spatial_source_missing")
    frame
  })
  if (!is.null(cells)) {
    report$checks[[length(report$checks)]]$details <- list(rows = nrow(cells))
    statuses <- sort(unique(as.character(cells$donor_status)))
    modes <- list(phenotype = sort(unique(as.character(cells$phenotype))),
                  marker = source_env$cell_distance_markers())
    # The distance to the nearest Islet or Cluster, and (roots from v9) to the nearest Islet only.
    distance_columns <- "signed_dist"
    if (source_env$cell_distance_has_islet_distance()) {
      if (!release_env$cell_distance_has_islet_distance()) parity_fail("islet_distance_availability_mismatch")
      distance_columns <- c(distance_columns, "signed_dist_islet")
    }
    for (mode in names(modes)) {
      items <- modes[[mode]]
      items <- items[!is.na(items) & nzchar(items)]
      for (norm in c("prop", "count")) {
        record(paste0("spatial_histogram_", mode, "_", norm), {
          rows <- 0L
          for (dist_col in distance_columns) for (item in c("__all__", items)) {
            result <- parity_query_result(
              source_env$.query_cell_distance_hist_memory(mode, item, statuses, norm, 120L, dist_col),
              release_env$.query_cell_distance_hist_duckdb(mode, item, statuses, norm, 120L, dist_col), "hist")
            rows <- rows + result$rows
          }
          list(selections = length(items) + 1L, distances = length(distance_columns), rows = rows)
        })
        record(paste0("spatial_composition_", mode, "_", norm), {
          rows <- 0L
          for (dist_col in distance_columns) {
            result <- parity_query_result(
              source_env$.query_cell_distance_composition_memory(mode, items, statuses, norm, dist_col),
              release_env$.query_cell_distance_composition_duckdb(mode, items, statuses, norm, dist_col), "comp")
            rows <- rows + result$rows
          }
          list(distances = length(distance_columns), rows = rows)
        })
      }
    }
  }

  record("geometry_representative_bounding_boxes", {
    sys.source(file.path(repo, "app", "shiny_app", "R", "segmentation_helpers.R"), envir = release_env)
    lookup <- read.csv(file.path(source_root, "islet_spatial_lookup.csv"), stringsAsFactors = FALSE)
    groups <- split(lookup, as.character(lookup$case_id))
    regions <- 0L; features <- 0L
    for (group in groups) {
      group <- group[order(group$islet_key), , drop = FALSE]
      case_id <- as.character(group$case_id[[1]])
      ids <- unique(c(case_id, sprintf("%04d", as.integer(case_id))))
      candidates <- unlist(lapply(ids, function(id) c(file.path(source_root, "json", paste0(id, ".geojson")),
        file.path(source_root, "gson", paste0(id, ".geojson.gz")))))
      candidates <- candidates[file.exists(candidates)]
      if (!length(candidates)) parity_fail("geometry_source_missing")
      geometry <- sf::st_read(builder$release_geometry_dsn(candidates[[1]]), quiet = TRUE,
                              options = "ARRAY_AS_STRING=YES", promote_to_multi = FALSE)
      sf::st_crs(geometry) <- NA
      for (index in unique(c(1L, ceiling(nrow(group) / 2), nrow(group)))) {
        row <- group[index, , drop = FALSE]
        x <- row$centroid_x_um / source_env$PIXEL_SIZE_UM
        y <- row$centroid_y_um / source_env$PIXEL_SIZE_UM
        b <- 200 / source_env$PIXEL_SIZE_UM
        bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = x - b, ymin = y - b, xmax = x + b, ymax = y + b)))
        expected <- geometry[lengths(sf::st_intersects(geometry, bbox)) > 0L, , drop = FALSE]
        actual <- release_env$load_case_geojson(case_id, row$centroid_x_um, row$centroid_y_um, 200)
        if (is.null(actual) || nrow(expected) != nrow(actual)) parity_fail("geometry_selection_mismatch")
        # Hash geometry bytes internally; evidence contains only aggregate counts.
        fingerprints <- function(frame) vapply(sf::st_as_binary(sf::st_geometry(frame), EWKB = FALSE),
          digest::digest, character(1), algo = "sha256", serialize = FALSE)
        a <- fingerprints(expected); z <- fingerprints(actual)
        if (!identical(sort(a), sort(z))) parity_fail("geometry_coordinate_mismatch")
        expected_fields <- sf::st_drop_geometry(expected)[order(a), , drop = FALSE]
        actual_fields <- sf::st_drop_geometry(actual)[order(z), , drop = FALSE]
        # GPKG may add fid; all source properties must be retained unchanged.
        actual_fields <- actual_fields[, names(expected_fields), drop = FALSE]
        expected_fields$.geometry_fingerprint <- a[order(a)]
        actual_fields$.geometry_fingerprint <- z[order(z)]
        parity_compare_frames(expected_fields, actual_fields, ".geometry_fingerprint")
        regions <- regions + 1L; features <- features + nrow(actual)
      }
    }
    list(donors = length(groups), regions = regions, selected_features = features, buffer_um = 200)
  })
  report$status <- if (any(vapply(report$checks, function(check) check$status != "passed", logical(1)))) {
    "failed"
  } else "passed"
  report
}

parity_write_report <- function(report, prefix) {
  paths <- paste0(prefix, c(".json", ".md"))
  if (any(file.exists(paths))) parity_fail("evidence_output_exists")
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(report, paths[[1]], auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 16)
  lines <- c("# Release data parity", "", paste("Status:", report$status), "",
             paste("Source commit:", report$identity$git_commit),
             paste("Source H5AD SHA-256:", report$identity$canonical_h5ad_sha256),
             paste("Release manifest SHA-256:", report$identity$release_manifest_sha256 %||% "not supplied"), "",
             "Numeric tolerance: absolute 1e-12. Categories, counts, missingness and infinities: exact.", "",
             "| Check | Result |", "|---|---|")
  lines <- c(lines, vapply(report$checks, function(check) {
    paste0("| ", check$check, " | ", check$status,
           if (check$status == "failed") paste0(" (", check$details$error_code, ")") else "", " |")
  }, character(1)), "", "No scientific rows, donor/image identifiers, coordinates, raw exceptions or absolute source paths are included.")
  writeLines(lines, paths[[2]], useBytes = TRUE)
}

parity_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- list()
  if (!length(args) || length(args) %% 2L != 0L) stop("Arguments require --name value pairs.", call. = FALSE)
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[i]])
    if (!(key %in% c("source-root", "image-root", "release-root", "output-prefix"))) {
      stop("Unsupported parity argument.", call. = FALSE)
    }
    options[[gsub("-", "_", key, fixed = TRUE)]] <- args[[i + 1L]]
  }
  if (is.null(options$output_prefix)) stop("--output-prefix is required.", call. = FALSE)
  if (any(file.exists(paste0(options$output_prefix, c(".json", ".md"))))) {
    stop("Evidence outputs already exist; choose a new prefix.", call. = FALSE)
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
  repo <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
  suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readxl) })
  required <- c("DBI", "duckdb", "data.table", "sf", "digest", "jsonlite", "anndata")
  if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required parity-check packages are unavailable.", call. = FALSE)
  }
  Sys.setenv(ISLET_RUNTIME_MODE = "development", ISLET_VERIFY_RELEASE_HASHES = "core")
  report <- parity_run(options, repo)
  parity_write_report(report, options$output_prefix)
  message("Release parity status: ", report$status)
  if (identical(report$status, "failed")) quit(status = 1L)
  invisible(report)
}

if (sys.nframe() == 0L) parity_main()
