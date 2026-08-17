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

#' Canonical phenotype vocabulary, in `PHENOTYPE_COLORS` order.
#'
#' Selectors must NOT derive their choices from whichever slice of data happens
#' to be loaded — that is how "CD4 Tcell" disappeared from the Spatial tab, since
#' the default donor (112) has zero CD4 cells and two other donors (6533, 6534)
#' do too. `PHENOTYPE_COLORS` (drilldown_helpers.R) is the documented single
#' source of truth for phenotype labels; `extra` lets an unexpected value from
#' the data still appear rather than being silently dropped.
#'
#' @param extra Optional character vector of phenotypes observed in the data.
#' @return Canonical names first (PHENOTYPE_COLORS order), then any extras.
canonical_phenotype_order <- function(extra = character(0)) {
  canon <- names(PHENOTYPE_COLORS)
  extra <- unique(as.character(extra[!is.na(extra) & nzchar(extra)]))
  c(canon, sort(setdiff(extra, canon)))
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

# ==========================================================================
# Card C-Right: per-cell signed-distance density + composition
# (ported from Islet-Explorer-Senior; adapted to CODEX's own data — all-cells
# `data/cell_distance.csv`, 31-marker CODEX phenotypes, NO RESTORE/DuckDB).
# ==========================================================================

# ── Symmetric-log distance transform ──────────────────────────────────────
# Signed distance-to-islet is negative inside the islet and ranges to ~+1.1 mm
# in bulk tissue. A linear axis would crush the biologically interesting
# islet-border transition into a sliver. T(x) = sign(x)*log10(1 + |x|/C) is
# near-linear within +/- C um and log-compressed beyond, so bins are ~1-2 um
# across the border and coarse far out.
CD_SYMLOG_C   <- 25
cd_symlog     <- function(x) sign(x) * log10(1 + abs(x) / CD_SYMLOG_C)
cd_symlog_inv <- function(t) sign(t) * CD_SYMLOG_C * (10^abs(t) - 1)
cd_symlog_ticks <- function(lo, hi) {
  cand <- c(-2000, -1000, -500, -200, -100, -50, -20, -10, 0, 10, 20, 50, 100, 200, 500, 1000, 2000)
  cand <- cand[cand >= lo & cand <= hi]
  if (!(0 %in% cand)) cand <- sort(unique(c(cand, 0)))
  list(vals = cd_symlog(cand), text = format(cand, trim = TRUE, scientific = FALSE))
}

# ── Cached all-cells cell-distance loader ─────────────────────────────────
.cell_distance_cache <- new.env(parent = emptyenv())
CELL_DISTANCE_PATH <- file.path("..", "..", "data", "cell_distance.csv")

cell_distance_available <- function() file.exists(CELL_DISTANCE_PATH)

#' Load the all-cells (~2.65M) tissue-wide signed-distance table with caching.
#' Columns: case_id, donor_status, cell_region, phenotype, signed_dist.
load_cell_distance <- function() {
  if (!is.null(.cell_distance_cache$data)) return(.cell_distance_cache$data)
  if (!file.exists(CELL_DISTANCE_PATH)) return(NULL)
  df <- tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      as.data.frame(data.table::fread(
        CELL_DISTANCE_PATH,
        colClasses = list(character = c("case_id", "donor_status", "cell_region", "phenotype"),
                          numeric = "signed_dist")))
    } else {
      read.csv(CELL_DISTANCE_PATH, stringsAsFactors = FALSE,
               colClasses = c(case_id = "character", donor_status = "character",
                              cell_region = "character", phenotype = "character",
                              signed_dist = "numeric"))
    }
  }, error = function(e) {
    message("[SPATIAL] cell_distance load failed: ", e$message); NULL
  })
  if (!is.null(df)) {
    df <- df[is.finite(df$signed_dist), , drop = FALSE]
    .cell_distance_cache$data <- df
    message("[SPATIAL] Cached cell_distance (", nrow(df), " cells)")
  }
  df
}

#' Distinct phenotypes present in the cell-distance data, in canonical order.
#' "Unknown" IS included — every phenotype must be reachable from every
#' selector; it is simply left unselected by default (see `dist_comp_items`).
cell_distance_phenotypes <- function() {
  df <- load_cell_distance()
  if (is.null(df)) return(character(0))
  canonical_phenotype_order(unique(df$phenotype))
}

#' Lineage markers with per-cell positivity (`{m}_pos` columns, scimap X>0.5).
cell_distance_markers <- function() {
  df <- load_cell_distance()
  if (is.null(df)) return(character(0))
  sub("_pos$", "", grep("_pos$", colnames(df), value = TRUE))
}

# Distinct colours for the marker-mode composition stacks.
CELL_DIST_MARKER_COLORS <- c(
  INS = "#E6194B", GCG = "#3CB44B", SST = "#B8A600", CD3e = "#4363D8",
  CD8a = "#F58231", CD4 = "#911EB4", CD20 = "#42D4F4", CD68 = "#F032E6",
  CD163 = "#9A9A00", HLADR = "#469990", CD31 = "#800000", CD34 = "#9A6324",
  SMA = "#808000", PDPN = "#000075", ECAD = "#E67E22", panCK = "#666666"
)

# Split `base` into the plotted subset. `mode` is "phenotype" (d$phenotype ==
# item) or "marker" (item's `{item}_pos` == 1). "__all__" returns everything.
.cd_subset <- function(base, mode, item) {
  if (is.null(item) || identical(item, "__all__")) return(base)
  if (identical(mode, "marker")) {
    pc <- paste0(item, "_pos")
    if (!(pc %in% colnames(base))) return(base[0, , drop = FALSE])
    base[base[[pc]] == 1, , drop = FALSE]
  } else {
    base[base$phenotype == item, , drop = FALSE]
  }
}

#' Density of signed distance-to-islet, pooled per cell, split by donor status,
#' over a symlog axis.
#'   mode  "phenotype" | "marker"; group "__all__" or a phenotype/marker.
#'   norm  "prop" = fraction of the disease class's cells (÷ that group's total,
#'         correcting for different total cell numbers per group) | "count".
#' Returns list(hist[donor_status,binx,n,x,real_x,dens], lo, hi, n_by_status,
#' norm). `dens` is per REAL um of distance either way.
query_cell_distance_hist <- function(mode, group, statuses, norm = "prop", n_bins = 120L) {
  df <- load_cell_distance()
  if (is.null(df) || length(statuses) == 0) return(NULL)
  base <- df[df$donor_status %in% statuses, , drop = FALSE]
  if (nrow(base) < 2) return(NULL)
  status_totals <- table(base$donor_status)          # class totals (denominator)
  d <- .cd_subset(base, mode, group)
  if (nrow(d) < 2) return(NULL)
  lo <- min(as.numeric(stats::quantile(d$signed_dist, 0.001, na.rm = TRUE)), -50)
  hi <- as.numeric(stats::quantile(d$signed_dist, 0.999, na.rm = TRUE))
  if (!isTRUE(hi > lo)) return(NULL)
  d <- d[d$signed_dist >= lo & d$signed_dist <= hi, , drop = FALSE]
  tlo <- cd_symlog(lo); thi <- cd_symlog(hi); dt <- (thi - tlo) / n_bins
  d$binx <- pmin(pmax(floor((cd_symlog(d$signed_dist) - tlo) / dt), 0L), n_bins - 1L)
  agg <- dplyr::count(d, donor_status, binx, name = "n")
  if (nrow(agg) == 0) return(NULL)
  tc <- tlo + (agg$binx + 0.5) * dt
  agg$x <- tc
  agg$real_x <- cd_symlog_inv(tc)
  rlo <- cd_symlog_inv(tlo + agg$binx * dt)
  rhi <- cd_symlog_inv(tlo + (agg$binx + 1) * dt)
  width <- pmax(rhi - rlo, 1e-6)
  if (identical(norm, "prop")) {
    denom <- as.numeric(status_totals[as.character(agg$donor_status)])
    agg$dens <- agg$n / denom / width
  } else {
    agg$dens <- agg$n / width
  }
  n_by_status <- tapply(agg$n, agg$donor_status, sum)
  list(hist = as.data.frame(agg), lo = lo, hi = hi, n_by_status = n_by_status, norm = norm)
}

#' Composition vs distance: per (status, distance bin, item) for stacked bars,
#' one facet per status. `mode` "phenotype" (mutually exclusive) or "marker" (a
#' cell may be + for several markers → stacks overlap). Bins ZERO-ANCHORED
#' (floor(dist/bw)) so a bar edge sits on x=0; window forced to straddle 0.
#' `norm` "prop" = fraction of the disease class's cells (÷ group total) |
#' "count". Returns list(comp[donor_status,binx,item,n,x,val], lo, hi, bw, norm).
query_cell_distance_composition <- function(mode, items, statuses, norm = "prop", n_bins = 16L) {
  df <- load_cell_distance()
  if (is.null(df) || length(statuses) == 0 || length(items) == 0) return(NULL)
  base <- df[df$donor_status %in% statuses, , drop = FALSE]
  if (nrow(base) == 0) return(NULL)
  status_totals <- table(base$donor_status)
  if (identical(mode, "marker")) {
    pcs <- paste0(items, "_pos"); pcs <- pcs[pcs %in% colnames(base)]
    if (length(pcs) == 0) return(NULL)
    sel <- base[rowSums(base[, pcs, drop = FALSE]) > 0, , drop = FALSE]
  } else {
    sel <- base[base$phenotype %in% items, , drop = FALSE]
  }
  if (nrow(sel) == 0) return(NULL)
  lo <- max(as.numeric(stats::quantile(sel$signed_dist, 0.005, na.rm = TRUE)), -100)
  hi <- min(as.numeric(stats::quantile(sel$signed_dist, 0.995, na.rm = TRUE)), 300)
  lo <- min(lo, -30); hi <- max(hi, 30)              # always straddle the islet border
  bw <- (hi - lo) / n_bins
  base <- base[base$signed_dist >= lo & base$signed_dist <= hi, , drop = FALSE]
  base$binx <- floor(base$signed_dist / bw)          # zero-anchored bin edges
  rows <- list()
  for (it in items) {
    sub <- .cd_subset(base, mode, it)
    if (nrow(sub) == 0) next
    cc <- dplyr::count(sub, donor_status, binx, name = "n")
    cc$item <- it
    rows[[length(rows) + 1L]] <- cc
  }
  if (length(rows) == 0) return(NULL)
  comp <- do.call(rbind, rows)
  comp$x <- (comp$binx + 0.5) * bw
  if (identical(norm, "prop")) {
    denom <- as.numeric(status_totals[as.character(comp$donor_status)])
    comp$val <- comp$n / denom
  } else {
    comp$val <- comp$n
  }
  list(comp = as.data.frame(comp), lo = lo, hi = hi, bw = bw, norm = norm)
}
