# Production-only dependency checks layered on top of the application's package
# manifest. Keep version changes deliberate: update the constants, rebuild, run
# scientific-parity tests, and commit the resulting release metadata together.

source("/opt/islet-explorer-build/install_shiny_deps.R", local = new.env())

expected <- c(
  destiny = "3.22.0",
  duckdb = "1.5.2"
)

repos <- getOption("repos")
if (is.null(repos) || !length(repos) || any(is.na(repos)) || any(repos == "@CRAN@")) {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}

# R 4.5 uses Bioconductor 3.21. Setting it explicitly prevents a future default
# release from silently changing the trajectory implementation.
if (as.character(BiocManager::version()) != "3.21") {
  BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
}
if (!requireNamespace("destiny", quietly = TRUE) ||
    as.character(utils::packageVersion("destiny")) != expected[["destiny"]]) {
  BiocManager::install("destiny", ask = FALSE, update = FALSE)
}

local({
  # rocker/shiny:4.5.1 retains a 2025 CRAN snapshot, before DuckDB 1.5.2.
  # Install the reviewed CRAN source archive, checking its bytes before R can
  # execute any package configure/install code. Never substitute current CRAN.
  duckdb_source <- "https://cran.r-project.org/src/contrib/Archive/duckdb/duckdb_1.5.2.tar.gz"
  duckdb_sha256 <- "55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6"
  archive <- tempfile(pattern = "duckdb-1.5.2-", fileext = ".tar.gz")
  patched_archive <- tempfile(pattern = "duckdb-1.5.2-brotli-", fileext = ".tar.gz")
  receipt <- "/opt/islet-explorer-build/duckdb-brotli-source.json"
  tryCatch({
    utils::download.file(duckdb_source, archive, mode = "wb", method = "libcurl")
    actual_sha256 <- digest::digest(file = archive, algo = "sha256", serialize = FALSE)
    if (!identical(actual_sha256, duckdb_sha256)) {
      stop("DuckDB source SHA-256 mismatch; refusing to install.")
    }
    # The package version alone cannot distinguish the vulnerable decoder from
    # this backport. Always rebuild the verified nested source, even when the
    # base already has DuckDB 1.5.2. Do not reuse precompiled object files.
    patch_status <- system2("python3", c(
      "/opt/islet-explorer-build/patch-duckdb-brotli.py",
      "--source", shQuote(archive), "--output", shQuote(patched_archive),
      "--receipt", shQuote(receipt)
    ))
    if (patch_status != 0L) stop("DuckDB Brotli source patch verification failed.")
    # DuckDB otherwise chooses build parallelism from host resources. Bound it
    # for a predictable workstation/container build without changing the binary.
    old_makeflags <- Sys.getenv("MAKEFLAGS", unset = NA_character_)
    old_prebuilt <- Sys.getenv("DUCKDB_R_PREBUILT_ARCHIVE", unset = NA_character_)
    Sys.setenv(MAKEFLAGS = "-j2")
    Sys.unsetenv("DUCKDB_R_PREBUILT_ARCHIVE")
    tryCatch({
      install_status <- system2(file.path(R.home("bin"), "R"),
        c("CMD", "INSTALL", "--preclean", shQuote(patched_archive)))
      if (install_status != 0L) stop("Patched DuckDB source installation failed.")
      Sys.chmod(receipt, mode = "0444")
    }, finally = {
      if (is.na(old_makeflags)) Sys.unsetenv("MAKEFLAGS") else
        Sys.setenv(MAKEFLAGS = old_makeflags)
      if (is.na(old_prebuilt)) Sys.unsetenv("DUCKDB_R_PREBUILT_ARCHIVE") else
        Sys.setenv(DUCKDB_R_PREBUILT_ARCHIVE = old_prebuilt)
    })
  }, finally = unlink(c(archive, patched_archive)))
})

actual <- vapply(names(expected), function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}, character(1))

bad <- names(expected)[is.na(actual) | actual != expected]
if (length(bad)) {
  details <- paste(sprintf("%s expected=%s actual=%s", bad, expected[bad], actual[bad]), collapse = "; ")
  stop("Production R dependency pin mismatch: ", details)
}

message("Production R dependency pins verified: ",
        paste(sprintf("%s=%s", names(actual), actual), collapse = ", "))
