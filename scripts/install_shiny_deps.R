#!/usr/bin/env Rscript

message("Installing Shiny app dependencies…")

install_if_missing <- function(pkgs, repos = getOption("repos")) {
  inst <- rownames(installed.packages())
  to_get <- setdiff(pkgs, inst)
  if (length(to_get)) {
    install.packages(to_get, repos = repos, dependencies = TRUE)
  }
}

# Base CRAN packages
cran_pkgs <- c(
  "shiny", "shinyjs", "readxl", "dplyr", "stringr", "tidyr",
  "ggplot2", "plotly", "broom", "jsonlite", "base64enc", "RColorBrewer"
)
install_if_missing(cran_pkgs)

# Optional: httr2 for AI assistant
install_if_missing(c("httr2"))

# vitessceR from GitHub (if not installed or too old)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
need_vitessceR <- TRUE
if (requireNamespace("vitessceR", quietly = TRUE)) {
  v <- tryCatch(packageVersion("vitessceR"), error = function(e) "0.0.0")
  need_vitessceR <- isTRUE(v < "0.1.0")
}
if (need_vitessceR) {
  remotes::install_github("vitessce/vitessceR", upgrade = "never", dependencies = TRUE)
}

# Optional: anndata (Bioconductor) for trajectory analysis
if (!requireNamespace("anndata", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  try(BiocManager::install("anndata", ask = FALSE, update = FALSE))
}

message("Done. Restart R before running the Shiny app if new packages were installed.")
