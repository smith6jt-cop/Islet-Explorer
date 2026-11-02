#!/usr/bin/env Rscript
# Compatibility check for the DuckDB 1.5.2 Brotli security backport.
# Use an existing, read-only release Parquet; all output goes to a temporary
# directory. This complements deploy/verify-brotli.py's vulnerability probe.
# Usage: Rscript --vanilla check-duckdb-brotli.R /input/islets.parquet

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop("Usage: check-duckdb-brotli.R REAL_RELEASE_PARQUET", call. = FALSE)
  }
  for (package in c("DBI", "duckdb", "digest", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing required package: ", package, call. = FALSE)
    }
  }
  input <- normalizePath(args[[1L]], mustWork = TRUE)
  if (dir.exists(input) || !grepl("\\.parquet$", input, ignore.case = TRUE)) {
    stop("Input must be an existing release Parquet file", call. = FALSE)
  }
  assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  before_sha256 <- sha256(input)
  package_version <- as.character(utils::packageVersion("duckdb"))
  assert(identical(package_version, "1.5.2"),
         paste("Expected production DuckDB 1.5.2, found", package_version))

  work <- tempfile("duckdb-brotli-compat-", tmpdir = tempdir())
  assert(dir.create(work, mode = "0700"), "Cannot create temporary directory")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  output <- file.path(work, "roundtrip.parquet")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads = 1")
  DBI::dbExecute(con, "SET preserve_insertion_order = true")
  literal <- function(x) as.character(DBI::dbQuoteString(con, x))
  identifier <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
  query <- function(sql) DBI::dbGetQuery(con, sql)
  source <- paste0("read_parquet(", literal(input), ")")
  restored <- paste0("read_parquet(", literal(output), ")")
  schema <- query(paste("DESCRIBE SELECT * FROM", source))
  original <- query(paste("SELECT * FROM", source))
  assert(nrow(original) > 0L && ncol(original) > 1L,
         "Use a nonempty, representative real release table")
  numeric_columns <- schema$column_name[grepl(
    "^(U?(TINYINT|SMALLINT|INTEGER|BIGINT|HUGEINT)|FLOAT|DOUBLE|DECIMAL)",
    schema$column_type
  )]
  text_columns <- schema$column_name[schema$column_type == "VARCHAR"]
  assert(length(numeric_columns) > 0L && length(text_columns) > 0L,
         "Use a real release table containing both numeric and text columns")

  DBI::dbExecute(con, paste0(
    "COPY (SELECT * FROM ", source, ") TO ", literal(output),
    " (FORMAT PARQUET, COMPRESSION BROTLI)"
  ))
  compression <- query(paste0(
    "SELECT compression, COUNT(*) AS column_chunks, ",
    "SUM(num_values) AS values_encoded FROM parquet_metadata(",
    literal(output), ") GROUP BY compression ORDER BY compression"
  ))
  assert(nrow(compression) == 1L &&
           identical(compression$compression, "BROTLI") &&
           compression$column_chunks[[1L]] >= ncol(original) &&
           compression$values_encoded[[1L]] > 0,
         "Round-trip Parquet metadata does not prove BROTLI for every column chunk")
  restored_schema <- query(paste("DESCRIBE SELECT * FROM", restored))
  assert(identical(schema, restored_schema), "Logical schema changed after BROTLI round trip")
  roundtrip <- query(paste("SELECT * FROM", restored))
  assert(identical(original, roundtrip),
         "Ordered rows, R column types, missingness, or values changed after BROTLI round trip")
  # EXCEPT ALL independently compares SQL values while preserving duplicate counts.
  delta_rows <- function(left, right) {
    query(paste0("SELECT COUNT(*) AS n FROM ((SELECT * FROM ", left,
                 ") EXCEPT ALL (SELECT * FROM ", right, ")) AS delta"))$n[[1L]]
  }
  assert(delta_rows(source, restored) == 0 && delta_rows(restored, source) == 0,
         "SQL values or duplicate multiplicities changed after BROTLI round trip")

  # All-column missingness and all-numeric range/sum profiles exercise real data.
  profile_terms <- vapply(schema$column_name, function(column) {
    paste0("COUNT(", identifier(column), ") AS ", identifier(paste0("nonnull_", column)))
  }, character(1))
  for (column in numeric_columns) {
    for (operation in c("MIN", "MAX", "SUM")) {
      profile_terms <- c(profile_terms, paste0(
        operation, "(", identifier(column), ") AS ",
        identifier(paste0(tolower(operation), "_", column))
      ))
    }
  }
  profile_sql <- paste(c("COUNT(*) AS total_rows", profile_terms), collapse = ", ")
  source_profile <- query(paste("SELECT", profile_sql, "FROM", source))
  restored_profile <- query(paste("SELECT", profile_sql, "FROM", restored))
  assert(identical(source_profile, restored_profile),
         "Column missingness or numeric aggregates changed after BROTLI round trip")

  groups <- intersect(c("Case ID", "Donor Status", "object_class",
                        "imageid", "donor_status", "islet_class"), schema$column_name)
  if (!length(groups)) groups <- text_columns[[1L]]
  groups <- head(groups, 3L)
  metrics <- intersect(c("cells_total", "islet_diam_um", "core_region_um2",
                         "peri_region_um2", "Ins_any", "Glu_any", "Stt_any"), numeric_columns)
  if (!length(metrics)) metrics <- head(numeric_columns, 3L)
  group_sql <- paste(identifier(groups), collapse = ", ")
  aggregate_terms <- c(group_sql, "COUNT(*) AS row_count")
  for (column in metrics) {
    for (operation in c("COUNT", "SUM", "AVG")) {
      aggregate_terms <- c(aggregate_terms, paste0(
        operation, "(", identifier(column), ") AS ",
        identifier(paste0(tolower(operation), "_", column))
      ))
    }
  }
  aggregate_query <- function(relation) paste0(
    "SELECT ", paste(aggregate_terms, collapse = ", "), " FROM ", relation,
    " GROUP BY ", group_sql, " ORDER BY ", group_sql
  )
  source_aggregate <- query(aggregate_query(source))
  restored_aggregate <- query(aggregate_query(restored))
  assert(identical(source_aggregate, restored_aggregate),
         "Representative grouped aggregates changed after BROTLI round trip")
  assert(identical(before_sha256, sha256(input)), "Input file changed during validation")

  report <- list(
    result = "pass",
    scope = "Real release Parquet compatibility; vulnerability probe is separate",
    checked_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    duckdb_r_package_version = package_version,
    duckdb_sql_version = query("SELECT version() AS version")$version[[1L]],
    input = list(path = input, sha256 = before_sha256,
                 bytes = unname(file.info(input)$size), rows = nrow(original),
                 columns = ncol(original), schema = schema[, c("column_name", "column_type")]),
    output = list(bytes = unname(file.info(output)$size), sha256 = sha256(output),
                  temporary_file_removed_on_exit = TRUE, compression_metadata = compression),
    checks = list(logical_schema_identical = TRUE, ordered_values_identical = TRUE,
                  bidirectional_except_all_empty = TRUE, all_column_profiles_identical = TRUE,
                  grouped_aggregates_identical = TRUE, source_sha256_unchanged = TRUE),
    aggregates = list(group_columns = groups, numeric_columns = metrics,
                      group_count = nrow(source_aggregate), values = source_aggregate)
  )
  cat(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE,
                      dataframe = "rows", na = "null", digits = NA), "\n")
}

tryCatch(main(), error = function(error) {
  message("FAIL: ", conditionMessage(error))
  quit(save = "no", status = 1L)
})
