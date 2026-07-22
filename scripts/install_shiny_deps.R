# Full runtime package manifest for the Islet Explorer Shiny app. Every name here
# is referenced by app/shiny_app/R/*.R via library() (00_globals.R) or pkg:: — verify with:
#   grep -rhoE "\b[A-Za-z][A-Za-z0-9.]+::" app/shiny_app/R/ | sed 's/:://' | sort -u
# Several packages are requireNamespace()-guarded in the app (it degrades gracefully)
# but are effectively REQUIRED in production: without `anndata` the app silently falls
# back to master_results.xlsx and the Plot "Population Measure" selector loses every
# cell type; without `data.table` Spatial Card C fails; without `sf` the segmentation
# and Viewer plots fail. Keep them installed.
#
# Repository resolution. Do NOT hardcode cloud.r-project.org: the rocker/shiny base image
# points `repos` at a Posit Package Manager snapshot serving precompiled Linux BINARIES.
# Overriding it forces every package to build from source, which both slows the container
# build enormously and drags in toolchain dependencies (cmake for nloptr -> lme4) that a
# minimal image does not carry. Honour whatever the image already configured; fall back to
# CRAN only when nothing usable is set (e.g. a bare `Rscript` on a dev box).
cran_repos <- function() {
  r <- getOption("repos")
  if (!is.null(r) && length(r) && !any(is.na(r)) && !any(r == "@CRAN@")) return(r)
  c(CRAN = "https://cloud.r-project.org")
}

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
if (length(to_install)) install.packages(to_install, repos = cran_repos())
message("Installed/verified: ", paste(pkgs, collapse=", "))

# Fail loudly if a REQUIRED package did not install. Previously this script always exited 0,
# so a container build could "succeed" and only surface the breakage as an HTTP 500 from the
# shiny-server worker (or, for anndata, as a silent fallback to master_results.xlsx).
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  stop("Required packages failed to install: ", paste(missing, collapse = ", "))
}

# --- OPTIONAL: Bioconductor + Lamian trajectory stack --------------------------
# NONE of slingshot/SingleCellExperiment/DelayedMatrixStats/destiny/tradeSeq/Lamian
# are referenced by app/shiny_app/R/*.R — the app reads a PRECOMPUTED trajectory
# (dpt_pseudotime, UMAP) from the H5AD, so it runs fine without these. They only
# support the offline scripts/ pipeline. (The dormant on-the-fly Lamian recompute path in
# data_loading.R was deleted 2026-07-22 — it had zero call sites.) All installs below are best-effort (try(); the script
# still succeeds if they fail), so a container build won't break on a Bioconductor
# compile error. Safe to skip entirely for a minimal app image.
#
# GATED behind INSTALL_BIOC=1 so the block matches the contract deploy/islet-explorer.def
# already states ("Bioconductor and Lamian are optional"). Previously it ran unconditionally
# and dominated container build time compiling a stack the app never loads.
#   INSTALL_BIOC=1 Rscript scripts/install_shiny_deps.R
if (identical(tolower(Sys.getenv("INSTALL_BIOC", "0")), "1")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repos())
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
    install.packages("lmtest", repos = cran_repos())
    message("Installed/verified: lmtest")
  }
} else {
  message("Skipping Bioconductor stack (set INSTALL_BIOC=1 to install slingshot/destiny/tradeSeq/lmtest)")
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
    install.packages("vitessceR", repos = cran_repos())
    ok <- requireNamespace("vitessceR", quietly = TRUE)
  }, silent = TRUE)
  if (!ok) {
    # Fallback to GitHub
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = cran_repos())
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

# Install Lamian (GitHub only) — optional, gated behind INSTALL_LAMIAN=1. The app has NO
# Lamian code path at all (the dormant one in data_loading.R was deleted 2026-07-22);
# trajectory ships precomputed in the H5AD. Kept only for offline experimentation.
if (identical(tolower(Sys.getenv("INSTALL_LAMIAN", "0")), "1")) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repos())
  }
  try({
    remotes::install_github("Winnie09/Lamian")
    if (requireNamespace("Lamian", quietly = TRUE)) message("Installed/verified: Lamian")
  }, silent = TRUE)
}
