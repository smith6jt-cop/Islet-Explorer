# ---------------------------------------------------------------------------
# scanpy DPT + PAGA force-directed layout (Barlow et al. 2024)
# ---------------------------------------------------------------------------
#
# R side of app/shiny_app/python/paga_pseudotime.py. DPT and PAGA share a
# cached neighbor graph, while the expensive Leiden/ForceAtlas2 layout is
# deferred until the user asks to show the PAGA card.
#
# WHY A LAYOUT AND NOT JUST AN AXIS
# Barlow's Figure 2B/2E is a PAGA force-directed layout -- a 2-D graph embedding
# with one point per islet -- not a feature-vs-pseudotime scatter. That is what
# makes their figures readable. Our scatter cannot be fixed by an axis transform
# because the pile-up is at EXACT zero: 95.6% of islets have zero CD8a T cells,
# 98.2% zero B cells, and log(0) is still 0 (Barlow's own y-axis is
# log(#cells/#endocrine), so their Fig 2F has the same pile-up). Measured on the
# same occupancy grid: the layout leaves 8.1% of islets overplotted with at most
# 3 per pixel, against 84-96% and 125-129 per pixel for the scatter.
#
# The optional PAGA layout is cached per cell-count range; its full fit is ~15 s
# at n = 2,310. The main scatter computes only the DPT axis.

.paga_cache <- new.env(parent = emptyenv())
.scanpy_pt_cache <- new.env(parent = emptyenv())

paga_module_path <- function() file.path("python", "paga_pseudotime.py")

#' TRUE when reticulate's interpreter can actually run this (scanpy present).
#' Result is memoised -- the probe itself crosses the python bridge.
paga_available <- function() {
  if (!is.null(.paga_cache$available)) return(.paga_cache$available)
  ok <- tryCatch({
    file.exists(paga_module_path()) &&
      requireNamespace("reticulate", quietly = TRUE) &&
      isTRUE(reticulate::py_module_available("scanpy"))
  }, error = function(e) FALSE)
  .paga_cache$available <- ok
  if (!ok) islet_log_event("paga_unavailable")
  ok
}

.paga_module <- function() {
  if (!is.null(.paga_cache$mod)) return(.paga_cache$mod)
  m <- tryCatch(reticulate::import_from_path(
         "paga_pseudotime", path = "python", convert = TRUE),
         error = function(e) {
           islet_log_event("paga_module_import_failed"); NULL
         })
  .paga_cache$mod <- m
  m
}

#' Compute only the Scanpy diffusion-pseudotime axis for a cell-count range.
#' Keep Leiden/PAGA/ForceAtlas2 out of the initial scatter render; the graph
#' layout is computed separately only when the user requests that card.
compute_scanpy_pseudotime <- function(min_cells = 1, max_cells = NA) {
  mc <- suppressWarnings(as.numeric(min_cells))
  if (!is.finite(mc) || mc < 1) mc <- 1
  xc <- suppressWarnings(as.numeric(max_cells))
  key <- paste0("dpt_", mc, "_", if (is.finite(xc)) xc else "inf")
  cached <- .scanpy_pt_cache[[key]]
  if (!is.null(cached)) return(if (isTRUE(cached$failed)) NULL else cached$data)
  if (!paga_available()) return(NULL)
  m <- .paga_module()
  if (is.null(m)) return(NULL)

  res <- tryCatch(m$compute_pseudotime(mc, if (is.finite(xc)) xc else NULL),
                  error = function(e) {
                    islet_log_event("scanpy_dpt_failed")
                    NULL
                  })
  if (is.null(res) || !is.null(res$error)) {
    islet_log_event("scanpy_dpt_failed")
    .scanpy_pt_cache[[key]] <- list(failed = TRUE)
    return(NULL)
  }
  out <- data.frame(islet_id = as.character(res$islet_id),
                    pseudotime = as.numeric(res$pseudotime),
                    stringsAsFactors = FALSE)
  .scanpy_pt_cache[[key]] <- list(failed = FALSE, data = out)
  out
}

#' Compute (or fetch cached) PAGA pseudotime + force-directed layout.
#'
#' Returns data.frame(islet_id, case_id, donor_status, leiden, fa_1, fa_2,
#' pseudotime) with attributes n_clusters / layout_key / root_donor, or NULL.
compute_paga_layout <- function(min_cells = 1, max_cells = NA) {
  mc <- suppressWarnings(as.numeric(min_cells))
  if (!is.finite(mc) || mc < 1) mc <- 1
  xc <- suppressWarnings(as.numeric(max_cells))
  key <- paste0("paga_", mc, "_", if (is.finite(xc)) xc else "inf")
  if (!is.null(.paga_cache[[key]])) return(.paga_cache[[key]])
  if (!paga_available()) return(NULL)
  m <- .paga_module()
  if (is.null(m)) return(NULL)

  res <- tryCatch(m$compute(mc, if (is.finite(xc)) xc else NULL), error = function(e) {
    islet_log_event("paga_compute_failed")
    NULL
  })
  if (is.null(res)) return(NULL)
  if (!is.null(res$error)) {
    islet_log_event("paga_compute_failed")
    # Cache the failure so a bad range does not re-run a 15 s fit on every
    # reactive invalidation. NULL would not be cached (the is.null check above).
    .paga_cache[[key]] <- NA
    return(NULL)
  }
  df <- data.frame(
    islet_id     = as.character(res$islet_id),
    case_id      = as.character(res$case_id),
    donor_status = as.character(res$donor_status),
    leiden       = as.character(res$leiden),
    fa_1         = as.numeric(res$fa_1),
    fa_2         = as.numeric(res$fa_2),
    pseudotime   = as.numeric(res$pseudotime),
    stringsAsFactors = FALSE
  )
  attr(df, "n_clusters") <- as.integer(res$n_clusters)
  attr(df, "layout_key") <- as.character(res$layout_key)
  attr(df, "root_donor") <- as.character(res$root_donor)
  attr(df, "min_cells")  <- mc
  attr(df, "max_cells")  <- xc
  .paga_cache[[key]] <- df
  df
}

# `NA` is the cached-failure sentinel; callers want NULL.
paga_result <- function(min_cells = 1, max_cells = NA) {
  r <- compute_paga_layout(min_cells, max_cells)
  if (is.null(r) || (length(r) == 1 && is.na(r))) NULL else r
}
