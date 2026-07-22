# Full runtime package manifest for the Islet Explorer Shiny app. Every name here
# is referenced by app/shiny_app/R/*.R via library() (00_globals.R) or pkg:: — verify with:
#   grep -rhoE "\b[A-Za-z][A-Za-z0-9.]+::" app/shiny_app/R/ | sed 's/:://' | sort -u
# Several packages are requireNamespace()-guarded in the app (it degrades gracefully)
# but are effectively REQUIRED in production: without `anndata` the app silently falls
# back to master_results.xlsx and the Plot "Population Measure" selector loses every
# cell type; without `data.table` Spatial Card C fails; without `sf` the segmentation
# and Viewer plots fail. Keep them installed.
pkgs <- c(
  # Core UI / data — library()'d unconditionally in 00_globals.R
  "shiny","shinyjs","readxl","dplyr","stringr","tidyr","ggplot2","plotly","broom","jsonlite",
  # Runtime-critical, but requireNamespace()-guarded in the app
  "sf","anndata","reticulate","data.table","magick","xml2","base64enc","httr2","curl",
  # Rendering / helpers
  "RColorBrewer","scales","htmltools","htmlwidgets","cluster",
  # Phase 16 statistics: mixed-effects model + post-hoc contrasts. Required by
  # utils_stats.R (lmerTest::lmer) and mod_statistics_server.R (emmeans::emmeans).
  # If missing in /usr/local/lib/R/site-library the shiny-server worker fails
  # to initialize with HTTP 500 since user `shiny` can't see personal user libs.
  "lmerTest","emmeans","lme4","car"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
message("Installed/verified: ", paste(pkgs, collapse=", "))

# --- OPTIONAL: Bioconductor + Lamian trajectory stack --------------------------
# NONE of slingshot/SingleCellExperiment/DelayedMatrixStats/destiny/tradeSeq/Lamian
# are referenced by app/shiny_app/R/*.R — the app reads a PRECOMPUTED trajectory
# (dpt_pseudotime, UMAP) from the H5AD, so it runs fine without these. They only
# support the dormant on-the-fly Lamian recompute path in data_loading.R and the
# offline scripts/ pipeline. All installs below are best-effort (try(); the script
# still succeeds if they fail), so a container build won't break on a Bioconductor
# compile error. Safe to skip entirely for a minimal app image.

# Install Bioconductor packages (slingshot and friends)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("slingshot", "SingleCellExperiment", "DelayedMatrixStats")
try({
  BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
  message("Installed/verified (Bioconductor): ", paste(bioc_pkgs, collapse=", "))
}, silent = TRUE)

# destiny (DiffusionMap) and tradeSeq from Bioconductor; lmtest from CRAN
try({
  BiocManager::install(c("destiny","tradeSeq"), ask = FALSE, update = FALSE)
  message("Installed/verified (Bioconductor): destiny, tradeSeq")
}, silent = TRUE)

if (!requireNamespace("lmtest", quietly = TRUE)) {
  install.packages("lmtest", repos = "https://cloud.r-project.org")
  message("Installed/verified: lmtest")
}

# Install vitessceR (from CRAN if available; fallback to GitHub)
install_vitesscer <- function() {
  if (requireNamespace("vitessceR", quietly = TRUE)) {
    message("vitessceR already installed")
    return(invisible(TRUE))
  }
  ok <- FALSE
  # Try CRAN first
  try({
    install.packages("vitessceR", repos = "https://cloud.r-project.org")
    ok <- requireNamespace("vitessceR", quietly = TRUE)
  }, silent = TRUE)
  if (!ok) {
    # Fallback to GitHub
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    try({
      remotes::install_github("vitessce/vitessceR")
      ok <- requireNamespace("vitessceR", quietly = TRUE)
    }, silent = TRUE)
  }
  if (ok) message("Installed/verified: vitessceR") else warning("vitessceR install did not complete; see logs")
  invisible(ok)
}

# Only install vitessceR if requested explicitly (not required for the embedded viewer)
if (identical(tolower(Sys.getenv("INSTALL_VITESSCER", "0")), "1")) {
  install_vitesscer()
}

# Install Lamian (GitHub only) — optional, gated behind INSTALL_LAMIAN=1. Only used by
# the on-the-fly Lamian pseudotime path in data_loading.R, which the app does not
# exercise (trajectory ships precomputed in the H5AD).
if (identical(tolower(Sys.getenv("INSTALL_LAMIAN", "0")), "1")) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  try({
    remotes::install_github("Winnie09/Lamian")
    if (requireNamespace("Lamian", quietly = TRUE)) message("Installed/verified: Lamian")
  }, silent = TRUE)
}
