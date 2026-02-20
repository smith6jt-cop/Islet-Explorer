# ===== Spatial Tab Helper Functions =====
# Provides per-donor tissue-wide cell loading for the tissue scatter plot.
# Loads per-donor CSV files from data/donors/ with env-based caching.

# ── Cache environment for per-donor tissue data ──────────────────────────
.donor_tissue_cache <- new.env(parent = emptyenv())

# Path to per-donor tissue CSVs
DONORS_DIR <- file.path("..", "..", "data", "donors")

#' Check if per-donor tissue data is available
donor_tissue_available <- function() {
  dir.exists(DONORS_DIR) && length(list.files(DONORS_DIR, pattern = "\\.csv$")) > 0
}

#' Get list of available donor image IDs
get_available_donors <- function() {
  files <- list.files(DONORS_DIR, pattern = "\\.csv$")
  gsub("\\.csv$", "", files)
}

#' Load per-donor tissue cells with caching
#' @param imageid Donor image ID (e.g., "6505")
#' @return data.frame with X_centroid, Y_centroid, phenotype, cell_region, islet_name
load_donor_tissue <- function(imageid) {
  cache_key <- as.character(imageid)

  if (exists(cache_key, envir = .donor_tissue_cache)) {
    return(get(cache_key, envir = .donor_tissue_cache))
  }

  fpath <- file.path(DONORS_DIR, paste0(imageid, ".csv"))
  if (!file.exists(fpath)) return(NULL)

  df <- tryCatch({
    d <- read.csv(fpath, stringsAsFactors = FALSE)
    if (nrow(d) == 0) return(NULL)
    d
  }, error = function(e) {
    message("[SPATIAL] Error loading donor tissue for ", imageid, ": ", e$message)
    NULL
  })

  if (!is.null(df)) {
    assign(cache_key, df, envir = .donor_tissue_cache)
    message("[SPATIAL] Cached donor tissue for ", imageid, " (", nrow(df), " cells)")
  }

  df
}
