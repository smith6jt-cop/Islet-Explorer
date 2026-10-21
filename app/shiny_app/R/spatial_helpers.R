# ===== Spatial Tab Helper Functions =====
# Provides per-donor tissue-wide cell loading for the tissue scatter plot.
# Loads per-donor CSV files from data/donors/ with env-based caching.

# ── Cache environment for per-donor tissue data ──────────────────────────
.donor_tissue_cache <- new.env(parent = emptyenv())

# Repository-development CSVs and publication-release Parquet partitions.
DONORS_DIR <- data_path("donors")
DONOR_TISSUE_PARQUET_DIR <- data_path("spatial", "tissue")

.donor_tissue_parquet <- function(imageid) {
  imageid <- as.character(imageid)
  if (!islet_safe_path_segment(imageid)) return(NULL)
  id_num <- suppressWarnings(as.integer(imageid))
  ids <- unique(c(as.character(imageid), if (is.finite(id_num)) sprintf("%04d", id_num)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  for (id in ids) {
    d <- file.path(DONOR_TISSUE_PARQUET_DIR, paste0("case_id=", id))
    hits <- if (dir.exists(d)) list.files(d, pattern = "[.]parquet$", full.names = TRUE) else character(0)
    if (length(hits)) return(hits[[1]])
  }
  NULL
}

#' Check if per-donor tissue data is available
donor_tissue_available <- function() {
  has_parquet <- dir.exists(DONOR_TISSUE_PARQUET_DIR) &&
    length(list.files(DONOR_TISSUE_PARQUET_DIR, pattern = "[.]parquet$",
                      recursive = TRUE)) > 0
  has_parquet || (!islet_is_production() && dir.exists(DONORS_DIR) &&
                    length(list.files(DONORS_DIR, pattern = "[.]csv$")) > 0)
}

#' Get list of available donor image IDs
get_available_donors <- function() {
  pq <- if (dir.exists(DONOR_TISSUE_PARQUET_DIR)) {
    basename(list.dirs(DONOR_TISSUE_PARQUET_DIR, recursive = FALSE, full.names = TRUE))
  } else character(0)
  pq <- sub("^case_id=", "", pq[grepl("^case_id=", pq)])
  csv <- if (!islet_is_production() && dir.exists(DONORS_DIR)) {
    sub("[.]csv$", "", list.files(DONORS_DIR, pattern = "[.]csv$"))
  } else character(0)
  sort(unique(c(pq, csv)))
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
  canon <- phenotype_keep(names(PHENOTYPE_COLORS))
  # A uniform root names its own label set (uniform_release.json); the palette also carries the frozen
  # typing's labels (Ductal, SMA+, B cell, ...), which that root never assigns.
  root_labels <- if (exists("root_phenotype_labels", mode = "function")) root_phenotype_labels() else NULL
  if (length(root_labels)) canon <- c(intersect(canon, root_labels), setdiff(phenotype_keep(root_labels), canon))
  # Extras pass the same retired-label gate; otherwise a retired label observed
  # in the data (absent from `canon` because it was gated out) came back here.
  extra <- phenotype_keep(unique(as.character(extra[!is.na(extra) & nzchar(extra)])))
  c(canon, sort(setdiff(extra, canon)))
}

#' Load per-donor tissue cells with caching
#' @param imageid Donor image ID (e.g., "6505")
#' @return data.frame with X_centroid, Y_centroid, phenotype, cell_region, islet_name
load_donor_tissue <- function(imageid) {
  imageid <- as.character(imageid)
  if (!islet_safe_path_segment(imageid)) {
    islet_log_validation_failure("spatial_donor_invalid")
    return(NULL)
  }
  cache_key <- imageid

  if (exists(cache_key, envir = .donor_tissue_cache)) {
    return(get(cache_key, envir = .donor_tissue_cache))
  }

  parquet_path <- .donor_tissue_parquet(imageid)
  csv_path <- file.path(DONORS_DIR, paste0(imageid, ".csv"))
  if (is.null(parquet_path) && (islet_is_production() || !file.exists(csv_path))) return(NULL)

  df <- tryCatch({
    d <- if (!is.null(parquet_path)) {
      islet_read_parquet(parquet_path)
    } else {
      read.csv(csv_path, stringsAsFactors = FALSE)
    }
    if (nrow(d) == 0) return(NULL)
    # Retired phenotypes dropped at the LOADER. Background ("tissue")
    # cells keep their own cell_region and are unaffected.
    phenotype_filter_rows(d)
  }, error = function(e) {
    islet_log_event("spatial_tissue_load_failed")
    NULL
  })

  if (!is.null(df)) {
    assign(cache_key, df, envir = .donor_tissue_cache)
  }

  df
}

# ==========================================================================
# Card C-Right: per-cell signed-distance density + composition
# (ported from Islet-Explorer-Senior; adapted to CODEX's own data — all-cells
# Parquet/DuckDB in production or CSV in development, 31-marker CODEX
# phenotypes, no RESTORE dependency).
# ==========================================================================

# ── Symmetric-log distance transform ──────────────────────────────────────
# Signed distance-to-islet is negative inside the islet and ranges to ~+1.1 mm
# in bulk tissue. A linear axis would crush the biologically interesting
# islet-border transition into a sliver. T(x) = sign(x)*log10(1 + |x|/C) is
# near-linear within +/- C um and log-compressed beyond, so bins are ~1-2 um
# across the border and coarse far out.
CD_SYMLOG_C   <- 25
# Owner, 2026-10-20: the profile runs to 200 um outside the border, in finer bins, with a line at the 20 um band edge.
CD_MAX_OUTSIDE_UM <- 200
CD_BAND_UM        <- 20
CD_COMP_BIN_UM    <- 5      # composition bar width (zero-anchored)
CD_HIST_BINS      <- 160L   # density bins over the symlog range
cd_symlog     <- function(x) sign(x) * log10(1 + abs(x) / CD_SYMLOG_C)
cd_symlog_inv <- function(t) sign(t) * CD_SYMLOG_C * (10^abs(t) - 1)
cd_symlog_ticks <- function(lo, hi) {
  cand <- c(-2000, -1000, -500, -200, -100, -50, -20, -10, 0, 10, 20, 50, 100, 200, 500, 1000, 2000)
  cand <- cand[cand >= lo & cand <= hi]
  if (!(0 %in% cand)) cand <- sort(unique(c(cand, 0)))
  list(vals = cd_symlog(cand), text = format(cand, trim = TRUE, scientific = FALSE))
}

# ── Embedded DuckDB / Parquet cell-distance backend ───────────────────────
.cell_distance_cache <- new.env(parent = emptyenv())
CELL_DISTANCE_PATH <- data_path("cell_distance.csv")
CELL_DISTANCE_PARQUET_PATH <- data_path("spatial", "cell_distance.parquet")
SIGNED_DISTANCE_KDE_PATH <- data_path("signed_distance_kde.csv")
SIGNED_DISTANCE_KDE_PARQUET_PATH <- data_path("spatial", "signed_distance_kde.parquet")

cell_distance_available <- function() {
  file.exists(CELL_DISTANCE_PARQUET_PATH) ||
    (!islet_is_production() && file.exists(CELL_DISTANCE_PATH))
}

.cd_using_parquet <- function() file.exists(CELL_DISTANCE_PARQUET_PATH)

.cd_kept_clause <- function(con = islet_duckdb_connection()) {
  retired <- PHENOTYPE_RETIRED %||% character(0)
  if (!length(retired)) return("TRUE")
  vals <- paste(vapply(retired, islet_sql_string, character(1), con = con), collapse = ", ")
  paste0("(phenotype IS NULL OR phenotype NOT IN (", vals, "))")
}

.cd_in_clause <- function(values, con = islet_duckdb_connection()) {
  paste(vapply(as.character(values), islet_sql_string, character(1), con = con), collapse = ", ")
}

.cd_backend_failure <- function(context, error) {
  islet_log_event("spatial_backend_failed")
  if (islet_is_production()) stop("Spatial data operation failed.", call. = FALSE)
  NULL
}

# The distance a profile is measured with (owner, 2026-10-20): signed_dist is to the nearest Islet or Cluster,
# signed_dist_islet (data roots from v9) to the nearest Islet only. The queries read either through these two helpers.
CD_DISTANCE_COLUMNS <- c("signed_dist", "signed_dist_islet")
cell_distance_has_islet_distance <- function() {
  if (.cd_using_parquet()) return("signed_dist_islet" %in% islet_parquet_columns(CELL_DISTANCE_PARQUET_PATH))
  df <- load_cell_distance()
  !is.null(df) && "signed_dist_islet" %in% names(df)
}
# In memory: the table with `signed_dist` holding the chosen distance (rows without one dropped).
.cd_with_distance <- function(df, dist_col) {
  if (is.null(df) || identical(dist_col, "signed_dist")) return(df)
  if (!(dist_col %in% CD_DISTANCE_COLUMNS) || !(dist_col %in% names(df))) return(NULL)
  df$signed_dist <- df[[dist_col]]
  keep <- is.finite(df$signed_dist)
  if (all(keep)) df else df[keep, , drop = FALSE]
}
# In DuckDB: the query with the chosen distance column in place of signed_dist.
.cd_sql <- function(sql, dist_col) {
  if (!(dist_col %in% CD_DISTANCE_COLUMNS)) stop("unknown distance column", call. = FALSE)
  gsub("\\bsigned_dist\\b", dist_col, sql, perl = TRUE)
}

#' Load the all-cells (~2.65M) tissue-wide signed-distance table with caching.
#' Columns: case_id, donor_status, cell_region, phenotype, signed_dist.
load_cell_distance <- function() {
  if (!is.null(.cell_distance_cache$data)) return(.cell_distance_cache$data)
  if (!cell_distance_available()) return(NULL)
  df <- tryCatch({
    if (.cd_using_parquet()) {
      islet_read_parquet(CELL_DISTANCE_PARQUET_PATH)
    } else if (requireNamespace("data.table", quietly = TRUE)) {
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
    islet_log_event("spatial_cell_distance_load_failed"); NULL
  })
  if (!is.null(df)) {
    df <- df[is.finite(df$signed_dist), , drop = FALSE]
    # Every analysed cell is kept (2026-10-20): Marker+ counts and the composition share's denominator must include
    # the cells under a retired label (unsubtyped "T cell" and "Immune", "Endocrine"), which dropping them at the
    # loader left out (22% of CD3e+ cells). The phenotype choices are gated in cell_distance_phenotypes().
    .cell_distance_cache$data <- df
  }
  df
}

#' Load the smaller signed-distance KDE table once per worker.
#' The publication release reads Parquet through DuckDB; CSV is development-only.
load_signed_distance_kde <- function() {
  if (!is.null(.cell_distance_cache$kde)) return(.cell_distance_cache$kde)
  df <- tryCatch({
    if (file.exists(SIGNED_DISTANCE_KDE_PARQUET_PATH)) {
      islet_read_parquet(SIGNED_DISTANCE_KDE_PARQUET_PATH)
    } else if (!islet_is_production() && file.exists(SIGNED_DISTANCE_KDE_PATH)) {
      read.csv(SIGNED_DISTANCE_KDE_PATH, stringsAsFactors = FALSE)
    } else NULL
  }, error = function(e) .cd_backend_failure("signed-distance KDE load", e))
  if (!is.null(df)) .cell_distance_cache$kde <- df
  df
}

#' Distinct phenotypes present in the cell-distance data, in canonical order.
#' Retired labels (including "Unknown") are removed by the shared
#' phenotype_keep()/phenotype_filter_rows() gate before reaching a selector.
cell_distance_phenotypes <- function() {
  if (.cd_using_parquet()) {
    return(tryCatch({
      con <- islet_duckdb_connection()
      rel <- islet_parquet_relation(CELL_DISTANCE_PARQUET_PATH, con)
      sql <- paste0("SELECT DISTINCT phenotype FROM ", rel,
                    " WHERE phenotype IS NOT NULL AND ", .cd_kept_clause(con),
                    " ORDER BY phenotype")
      canonical_phenotype_order(DBI::dbGetQuery(con, sql)$phenotype)
    }, error = function(e) .cd_backend_failure("phenotype query", e)))
  }
  df <- load_cell_distance()
  if (is.null(df)) return(character(0))
  canonical_phenotype_order(phenotype_keep(unique(df$phenotype)))
}

#' Lineage markers with per-cell positivity (`{m}_pos` columns, scimap X>0.5).
cell_distance_markers <- function() {
  if (.cd_using_parquet()) {
    return(sub("_pos$", "", grep("_pos$", islet_parquet_columns(CELL_DISTANCE_PARQUET_PATH),
                                  value = TRUE)))
  }
  df <- load_cell_distance()
  if (is.null(df)) return(character(0))
  sub("_pos$", "", grep("_pos$", colnames(df), value = TRUE))
}

# Distinct colours for the marker-mode composition stacks.
CELL_DIST_MARKER_COLORS <- c(
  INS = "#E6194B", GCG = "#3CB44B", SST = "#B8A600", CD3e = "#4363D8",
  CD8a = "#F58231", CD4 = "#911EB4", CD20 = "#42D4F4", CD68 = "#F032E6",
  CD163 = "#9A9A00", HLADR = "#469990", CD31 = "#800000", CD34 = "#9A6324",
  SMA = "#808000", PDPN = "#000075", ECAD = "#E67E22", panCK = "#666666",
  # The rest of the 31-marker panel, flagged from v8 (owner, 2026-10-20: "all markers that are common to all donors").
  CD45 = "#1F77B4", CHGA = "#FF7F0E", ColIV = "#2CA02C", VIM = "#D62728", LGALS3 = "#9467BD",
  B3TUBB = "#8C564B", KRT14 = "#E377C2", GAP43 = "#BCBD22", CD35 = "#17BECF", `PGP9.5` = "#AEC7E8",
  CD44 = "#FFBB78", NaKATPase = "#98DF8A", BCatenin = "#FF9896", BActin = "#C5B0D5", CK19 = "#C49C94"
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
.query_cell_distance_hist_memory <- function(mode, group, statuses, norm = "prop", n_bins = CD_HIST_BINS,
                                             dist_col = "signed_dist") {
  df <- .cd_with_distance(load_cell_distance(), dist_col)
  if (is.null(df) || length(statuses) == 0) return(NULL)
  base <- df[df$donor_status %in% statuses, , drop = FALSE]
  if (nrow(base) < 2) return(NULL)
  status_totals <- table(base$donor_status)          # class totals (denominator)
  d <- .cd_subset(base, mode, group)
  if (nrow(d) < 2) return(NULL)
  lo <- min(as.numeric(stats::quantile(d$signed_dist, 0.001, na.rm = TRUE)), -50)
  hi <- CD_MAX_OUTSIDE_UM
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
.query_cell_distance_composition_memory <- function(mode, items, statuses, norm = "prop",
                                                    dist_col = "signed_dist") {
  df <- .cd_with_distance(load_cell_distance(), dist_col)
  if (is.null(df) || length(statuses) == 0 || length(items) == 0) return(NULL)
  base <- df[df$donor_status %in% statuses, , drop = FALSE]
  if (nrow(base) == 0) return(NULL)
  if (identical(mode, "marker")) {
    pcs <- paste0(items, "_pos"); pcs <- pcs[pcs %in% colnames(base)]
    if (length(pcs) == 0) return(NULL)
    sel <- base[rowSums(base[, pcs, drop = FALSE]) > 0, , drop = FALSE]
  } else {
    sel <- base[base$phenotype %in% items, , drop = FALSE]
  }
  if (nrow(sel) == 0) return(NULL)
  lo <- max(as.numeric(stats::quantile(sel$signed_dist, 0.005, na.rm = TRUE)), -100)
  lo <- min(lo, -30); hi <- CD_MAX_OUTSIDE_UM         # always straddle the islet border
  bw <- CD_COMP_BIN_UM
  base <- base[base$signed_dist >= lo & base$signed_dist < hi, , drop = FALSE]
  base$binx <- floor(base$signed_dist / bw)          # zero-anchored bin edges
  # Every cell in each bin: the denominator of the share (owner, 2026-10-20). Dividing by the status's total cells
  # made each bar track how many cells lie at that distance, which is nearly constant from 20 to 120 um.
  bin_totals <- dplyr::count(base, donor_status, binx, name = "total")
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
  comp$total <- as.numeric(bin_totals$total[match(paste(comp$donor_status, comp$binx),
                                                  paste(bin_totals$donor_status, bin_totals$binx))])
  comp$val <- if (identical(norm, "prop")) comp$n / comp$total else comp$n
  list(comp = as.data.frame(comp), lo = lo, hi = hi, bw = bw, norm = norm)
}

.cd_item_predicate <- function(item, mode, columns, con) {
  if (is.null(item) || identical(item, "__all__")) return("TRUE")
  if (identical(mode, "marker")) {
    col <- paste0(item, "_pos")
    if (!(col %in% columns)) return("FALSE")
    paste0(islet_sql_identifier(col, con), " = 1")
  } else {
    paste0("phenotype = ", islet_sql_string(item, con))
  }
}

.cd_base_predicate <- function(statuses, con) {
  paste0(
    "isfinite(signed_dist) AND donor_status IN (", .cd_in_clause(statuses, con), ")"
  )
}

.query_cell_distance_hist_duckdb <- function(mode, group, statuses,
                                              norm = "prop", n_bins = CD_HIST_BINS, dist_col = "signed_dist") {
  if (!length(statuses)) return(NULL)
  n_bins <- suppressWarnings(as.integer(n_bins))
  if (!is.finite(n_bins) || n_bins < 1L) n_bins <- CD_HIST_BINS
  con <- islet_duckdb_connection()
  rel <- islet_parquet_relation(CELL_DISTANCE_PARQUET_PATH, con)
  columns <- islet_parquet_columns(CELL_DISTANCE_PARQUET_PATH)
  base_where <- .cd_base_predicate(statuses, con)
  item_where <- .cd_item_predicate(group, mode, columns, con)

  totals <- DBI::dbGetQuery(con, .cd_sql(paste0(
    "SELECT donor_status, count(*) AS n FROM ", rel,
    " WHERE ", base_where, " GROUP BY donor_status"
  ), dist_col))
  if (!nrow(totals) || sum(as.numeric(totals$n)) < 2) return(NULL)

  q <- DBI::dbGetQuery(con, .cd_sql(paste0(
    "SELECT quantile_cont(signed_dist, 0.001) AS q_lo, ",
    "quantile_cont(signed_dist, 0.999) AS q_hi, count(*) AS n FROM ", rel,
    " WHERE ", base_where, " AND ", item_where
  ), dist_col))
  if (!nrow(q) || as.numeric(q$n[[1]]) < 2) return(NULL)
  lo <- min(as.numeric(q$q_lo[[1]]), -50)
  hi <- CD_MAX_OUTSIDE_UM
  if (!isTRUE(is.finite(lo) && is.finite(hi) && hi > lo)) return(NULL)

  tlo <- cd_symlog(lo)
  thi <- cd_symlog(hi)
  dt <- (thi - tlo) / n_bins
  bin_expr <- sprintf(
    paste0("LEAST(GREATEST(CAST(FLOOR(((sign(signed_dist) * log10(1 + ",
           "abs(signed_dist) / %.17g)) - %.17g) / %.17g) AS BIGINT), 0), %d)"),
    CD_SYMLOG_C, tlo, dt, n_bins - 1L
  )
  agg <- DBI::dbGetQuery(con, .cd_sql(paste0(
    "SELECT donor_status, ", bin_expr, " AS binx, count(*) AS n FROM ", rel,
    " WHERE ", base_where, " AND ", item_where,
    " AND signed_dist >= ", format(lo, digits = 17, scientific = FALSE),
    " AND signed_dist <= ", format(hi, digits = 17, scientific = FALSE),
    " GROUP BY donor_status, binx ORDER BY donor_status, binx"
  ), dist_col))
  if (!nrow(agg)) return(NULL)
  agg$binx <- as.integer(agg$binx)
  agg$n <- as.numeric(agg$n)
  tc <- tlo + (agg$binx + 0.5) * dt
  agg$x <- tc
  agg$real_x <- cd_symlog_inv(tc)
  rlo <- cd_symlog_inv(tlo + agg$binx * dt)
  rhi <- cd_symlog_inv(tlo + (agg$binx + 1) * dt)
  width <- pmax(rhi - rlo, 1e-6)
  if (identical(norm, "prop")) {
    denom <- as.numeric(totals$n[match(agg$donor_status, totals$donor_status)])
    agg$dens <- agg$n / denom / width
  } else {
    agg$dens <- agg$n / width
  }
  n_by_status <- tapply(agg$n, agg$donor_status, sum)
  list(hist = as.data.frame(agg), lo = lo, hi = hi,
       n_by_status = n_by_status, norm = norm)
}

.query_cell_distance_composition_duckdb <- function(mode, items, statuses, norm = "prop",
                                                    dist_col = "signed_dist") {
  if (!length(statuses) || !length(items)) return(NULL)
  con <- islet_duckdb_connection()
  rel <- islet_parquet_relation(CELL_DISTANCE_PARQUET_PATH, con)
  columns <- islet_parquet_columns(CELL_DISTANCE_PARQUET_PATH)
  base_where <- .cd_base_predicate(statuses, con)
  item_predicates <- vapply(items, .cd_item_predicate, character(1),
                            mode = mode, columns = columns, con = con)
  valid <- item_predicates != "FALSE"
  items <- items[valid]
  item_predicates <- item_predicates[valid]
  if (!length(items)) return(NULL)

  select_where <- paste0("(", paste(item_predicates, collapse = " OR "), ")")
  q <- DBI::dbGetQuery(con, .cd_sql(paste0(
    "SELECT quantile_cont(signed_dist, 0.005) AS q_lo, ",
    "quantile_cont(signed_dist, 0.995) AS q_hi, count(*) AS n FROM ", rel,
    " WHERE ", base_where, " AND ", select_where
  ), dist_col))
  if (!nrow(q) || as.numeric(q$n[[1]]) == 0) return(NULL)
  lo <- max(as.numeric(q$q_lo[[1]]), -100)
  lo <- min(lo, -30)
  hi <- CD_MAX_OUTSIDE_UM
  if (!isTRUE(is.finite(lo) && is.finite(hi) && hi > lo)) return(NULL)
  bw <- CD_COMP_BIN_UM
  bin_expr <- sprintf("CAST(FLOOR(signed_dist / %.17g) AS BIGINT)", bw)
  range_where <- paste0(" AND signed_dist >= ", format(lo, digits = 17, scientific = FALSE),
                        " AND signed_dist < ", format(hi, digits = 17, scientific = FALSE))
  bin_totals <- DBI::dbGetQuery(con, .cd_sql(paste0(
    "SELECT donor_status, ", bin_expr, " AS binx, count(*) AS total FROM ", rel,
    " WHERE ", base_where, range_where, " GROUP BY donor_status, binx"), dist_col))

  rows <- vector("list", length(items))
  for (i in seq_along(items)) {
    rows[[i]] <- DBI::dbGetQuery(con, .cd_sql(paste0(
      "SELECT donor_status, ", bin_expr, " AS binx, count(*) AS n FROM ", rel,
      " WHERE ", base_where, " AND ", item_predicates[[i]], range_where,
      " GROUP BY donor_status, binx ORDER BY donor_status, binx"
    ), dist_col))
    if (nrow(rows[[i]])) rows[[i]]$item <- items[[i]]
  }
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (!length(rows)) return(NULL)
  comp <- do.call(rbind, rows)
  rownames(comp) <- NULL
  comp$binx <- as.integer(comp$binx)
  comp$n <- as.numeric(comp$n)
  comp$x <- (comp$binx + 0.5) * bw
  comp$total <- as.numeric(bin_totals$total[match(paste(comp$donor_status, comp$binx),
                                                  paste(bin_totals$donor_status, as.integer(bin_totals$binx)))])
  comp$val <- if (identical(norm, "prop")) comp$n / comp$total else comp$n
  list(comp = as.data.frame(comp), lo = lo, hi = hi, bw = bw, norm = norm)
}

query_cell_distance_hist <- function(mode, group, statuses, norm = "prop", n_bins = CD_HIST_BINS,
                                     dist_col = "signed_dist") {
  if (.cd_using_parquet()) {
    return(tryCatch(
      .query_cell_distance_hist_duckdb(mode, group, statuses, norm, n_bins, dist_col),
      error = function(e) .cd_backend_failure("cell-distance density query", e)
    ))
  }
  .query_cell_distance_hist_memory(mode, group, statuses, norm, n_bins, dist_col)
}

query_cell_distance_composition <- function(mode, items, statuses, norm = "prop", dist_col = "signed_dist") {
  if (.cd_using_parquet()) {
    return(tryCatch(
      .query_cell_distance_composition_duckdb(mode, items, statuses, norm, dist_col),
      error = function(e) .cd_backend_failure("cell-distance composition query", e)
    ))
  }
  .query_cell_distance_composition_memory(mode, items, statuses, norm, dist_col)
}


#' Leiden cluster labels in numeric order (0, 1, 2, ..., 10), not text order (0, 1, 10, 11, ..., 2).
leiden_levels <- function(x) {
  lv <- unique(as.character(x[!is.na(x)]))
  lv[order(suppressWarnings(as.numeric(lv)), lv)]
}
