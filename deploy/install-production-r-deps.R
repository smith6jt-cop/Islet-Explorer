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

if (!requireNamespace("duckdb", quietly = TRUE) ||
    as.character(utils::packageVersion("duckdb")) != expected[["duckdb"]]) {
  install.packages("duckdb", repos = repos)
}

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
