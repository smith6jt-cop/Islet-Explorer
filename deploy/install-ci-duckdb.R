#!/usr/bin/env Rscript

# Test-only Noble/R 4.3 binary. DuckDB recommends Posit binaries for Linux:
# https://github.com/duckdb/duckdb-r#installation-from-the-posit-public-package-manager
# Production retains its separately verified CRAN source archive and R 4.5.1.
install_ci_duckdb <- function() {
  os_release <- readLines("/etc/os-release", warn = FALSE)
  r_series <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  if (!all(c("ID=ubuntu", 'VERSION_ID="24.04"') %in% os_release) ||
      !identical(unname(Sys.info()[["machine"]]), "x86_64") || r_series != "4.3") {
    stop("The pinned CI DuckDB binary requires Ubuntu 24.04, x86_64, and R 4.3.")
  }
  if (!requireNamespace("digest", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE)) {
    stop("Install the required Ubuntu digest and DBI packages before DuckDB.")
  }

  runner_temp <- Sys.getenv("RUNNER_TEMP")
  library_dir <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(runner_temp) || !dir.exists(runner_temp) ||
      !nzchar(library_dir) || !dir.exists(library_dir)) {
    stop("RUNNER_TEMP and the isolated R_LIBS_USER directory must exist.")
  }
  if (!startsWith(normalizePath(library_dir), paste0(normalizePath(runner_temp), "/"))) {
    stop("The CI library must be inside RUNNER_TEMP.")
  }

  version <- "1.5.2"
  url <- paste0("https://packagemanager.posit.co/cran/2026-04-20/",
                "bin/linux/noble-x86_64/4.3/src/contrib/duckdb_1.5.2.tar.gz")
  expected_sha256 <- "2ab90d12564b8187db703344418c841e3e4a486ba75bcdec3613a820a032b3ae"
  expected_bytes <- 17506629
  archive <- tempfile("duckdb-ci-", tmpdir = runner_temp, fileext = ".tar.gz")
  metadata_dir <- tempfile("duckdb-ci-metadata-", tmpdir = runner_temp)
  dir.create(metadata_dir)
  on.exit(unlink(c(archive, metadata_dir), recursive = TRUE), add = TRUE)
  options(timeout = max(300, getOption("timeout")))
  utils::download.file(url, archive, method = "libcurl", mode = "wb", quiet = TRUE)
  if (file.info(archive)$size != expected_bytes ||
      digest::digest(file = archive, algo = "sha256") != expected_sha256) {
    stop("Pinned CI DuckDB binary size or SHA-256 mismatch; no source fallback is allowed.")
  }
  members <- utils::untar(archive, list = TRUE)
  if (!all(c("duckdb/DESCRIPTION", "duckdb/libs/duckdb.so") %in% members) ||
      any(startsWith(members, "duckdb/src/"))) {
    stop("Expected a compiled DuckDB package, not a source archive.")
  }
  utils::untar(archive, files = "duckdb/DESCRIPTION", exdir = metadata_dir)
  description <- read.dcf(file.path(metadata_dir, "duckdb", "DESCRIPTION"))
  if (description[1, "Package"] != "duckdb" || description[1, "Version"] != version ||
      !grepl("^R 4\\.3\\.[0-9]+; x86_64-pc-linux-gnu;", description[1, "Built"])) {
    stop("Pinned DuckDB binary package identity or R/platform build metadata is invalid.")
  }

  evidence_dir <- file.path(runner_temp, "islet-ci-r")
  dir.create(evidence_dir, showWarnings = FALSE)
  provenance <- data.frame(package = "duckdb", version = version, url = url,
                           sha256 = expected_sha256, bytes = expected_bytes,
                           built = description[1, "Built"])
  utils::write.table(provenance, file.path(evidence_dir, "duckdb-binary-provenance.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::install.packages(archive, lib = library_dir, repos = NULL, type = "source", quiet = TRUE)
  if (!requireNamespace("duckdb", quietly = TRUE, lib.loc = library_dir) ||
      as.character(utils::packageVersion("duckdb", lib.loc = library_dir)) != version) {
    stop("DuckDB 1.5.2 must load from the verified isolated CI library.")
  }
  connection <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  if (!identical(DBI::dbGetQuery(connection, "SELECT 42::INTEGER AS answer")$answer, 42L)) {
    stop("The verified DuckDB binary failed its connection/query smoke check.")
  }
  message("Verified DuckDB ", version, " binary installed; provenance recorded.")
}

install_ci_duckdb()
