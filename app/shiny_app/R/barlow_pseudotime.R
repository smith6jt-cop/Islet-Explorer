# ---------------------------------------------------------------------------
# Barlow-style islet composition pseudotime
# ---------------------------------------------------------------------------
#
# Ports the pseudotime construction from Barlow et al., eLife 2024;13:RP100535,
# "The extra-islet pancreas supports autoimmunity in human type 1 diabetes"
# (Nolan / Bollyky labs) -- Methods, "Islet segmentation and pseudotime".
#
# It is the closest published match to this project: same modality (CODEX), same
# unit (the ISLET, not the single cell), same 20 um peri-islet band, and a
# cohort smaller than ours (13 donors / ~2,000 islets vs our 15 / 5,214).
#
# Their recipe, and where each piece lives here:
#
#   "the number of each cell type inside the islet and between the islet edge
#    and 20 microns beyond"            -> data/barlow_islet_features.csv,
#                                         built by scripts/build_barlow_features.py
#   "divided by the number of endocrine cells inside the islet ... then
#    log-transformed"                  -> same script (per-islet, so it is safe
#                                         to precompute)
#   "cosine distance and 15 nearest neighbors"  -> DiffusionMap(k = 15, distance = "cosine")
#   "filtering out islets that had fewer than ten total cells"
#                                      -> the app's cell-count range, which is
#                                         why the trajectory is recomputed rather than filtered
#   "pseudotemporal distance from the centroid of non-T1D islets"
#                                      -> see `root` below
#
# NOT ported: their step of "temporarily omitting 25 islets positioned in the
# middle of the map" to isolate the inflamed path. That is a hand-picked,
# dataset-specific intervention; reproducing it here would be an unauditable
# researcher degree of freedom, and it cannot be defined on our data at all.
#
# The production path uses scanpy for the neighbor graph, DPT axis, Leiden/PAGA
# abstraction, and ForceAtlas2 layout. This file is the deterministic destiny
# fallback when scanpy is unavailable to reticulate; it produces only the DPT
# axis, not the PAGA layout card.

# Cache keyed by the lower and upper cell-count bounds: the fit is 0.3-10 s
# depending on how many islets survive, and users sweep the control repeatedly.
.barlow_cache <- new.env(parent = emptyenv())

barlow_features_path <- function() data_path("barlow_islet_features.csv")

barlow_available <- function() {
  file.exists(barlow_features_path()) &&
    requireNamespace("destiny", quietly = TRUE) &&
    requireNamespace("data.table", quietly = TRUE)
}

barlow_features <- function() {
  if (!is.null(.barlow_cache$features)) return(.barlow_cache$features)
  if (!barlow_available()) return(NULL)
  f <- tryCatch(data.table::fread(barlow_features_path()), error = function(e) NULL)
  .barlow_cache$features <- f
  f
}

#' Recompute Barlow-style pseudotime for a cell-count range.
#'
#' Returns a data.frame(islet_id, case_id, donor_status, pseudotime) covering
#' ONLY the islets in range, or NULL when unavailable. Cached per range.
compute_barlow_pseudotime <- function(min_cells = 1, max_cells = NA) {
  mc <- suppressWarnings(as.numeric(min_cells))
  if (!is.finite(mc) || mc < 1) mc <- 1
  xc <- suppressWarnings(as.numeric(max_cells))
  key <- paste0("pt_", mc, "_", if (is.finite(xc)) xc else "inf")
  if (!is.null(.barlow_cache[[key]])) return(.barlow_cache[[key]])

  d <- barlow_features()
  if (is.null(d) || !nrow(d)) return(NULL)
  keep <- d$total_cells_core >= mc
  if (is.finite(xc)) keep <- keep & d$total_cells_core <= xc
  sub <- d[keep, , drop = FALSE]
  # A diffusion map on a handful of points is meaningless; fail closed.
  if (nrow(sub) < 50) return(NULL)

  feat <- grep("^(core|peri)_", names(sub), value = TRUE)
  M <- as.matrix(sub[, feat, with = FALSE])
  keep_col <- apply(M, 2, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0)
  M <- M[, keep_col, drop = FALSE]
  if (ncol(M) < 3) return(NULL)

  # destiny drops duplicate rows internally, so eigenvectors() comes back SHORTER
  # than the input and every downstream index silently misaligns. Fit on the
  # unique rows and map back. Duplicates are common at low thresholds -- a 2-cell
  # islet's composition vector collides readily with another's.
  dup_key <- do.call(paste, c(as.data.frame(M), sep = "\r"))
  uniq <- !duplicated(dup_key)
  Mu <- M[uniq, , drop = FALSE]
  back <- match(dup_key, dup_key[uniq])

  res <- tryCatch({
    set.seed(1)  # destiny's eigendecomposition is randomised; pin it so the
                 # axis does not move between two runs at the same threshold.
    dm <- destiny::DiffusionMap(Mu, k = 15, n_eigs = 10,
                                distance = "cosine", verbose = FALSE)
    ev <- as.matrix(destiny::eigenvectors(dm))
    if (nrow(ev) != nrow(Mu)) stop("destiny returned ", nrow(ev),
                                   " eigenvector rows for ", nrow(Mu), " inputs")
    su <- sub[uniq, , drop = FALSE]

    # ROOT: "pseudotemporal distance from the centroid of non-T1D islets".
    # The centroid itself is not an islet, and DPT needs an actual vertex, so the
    # root is the ND islet NEAREST that centroid -- the medoid of the healthy
    # cloud. Candidates are restricted to ND: taking the globally nearest islet
    # let an Aab+ islet win at several thresholds, which anchors the axis in the
    # wrong group. (`non-T1D` is read as their "no T1D" control group, i.e. our
    # ND, not "every islet that is not from a T1D donor". Their groups are listed
    # as "no T1D, AA+, T1D", so AA+ is its own stage.)
    nd <- which(su$donor_status == "ND")
    if (!length(nd)) stop("no ND islets above this threshold; no root available")
    cen <- colMeans(ev[nd, 1:3, drop = FALSE])
    root <- nd[which.min(sqrt(rowSums(sweep(ev[nd, 1:3, drop = FALSE], 2, cen)^2)))]

    ptu <- destiny::DPT(dm, tips = root)[root, ]
    pt <- ptu[back]

    # Orient so that beta-cell content DECREASES with pseudotime. The diffusion
    # map's sign is arbitrary; without this the axis can flip between thresholds
    # and the same islet appears at both ends on consecutive settings.
    bc <- suppressWarnings(as.numeric(sub[["core_Beta cell"]]))
    fin <- is.finite(pt) & is.finite(bc)
    if (sum(fin) > 10 &&
        stats::cor(pt[fin], bc[fin], method = "spearman") > 0) {
      pt <- max(pt[is.finite(pt)], na.rm = TRUE) - pt
    }

    data.frame(
      islet_id     = as.character(sub$islet_id),
      case_id      = as.character(sub$case_id),
      donor_status = as.character(sub$donor_status),
      pseudotime   = as.numeric(pt),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("event=barlow_pseudotime_failed")
    NULL
  })

  .barlow_cache[[key]] <- res
  res
}

# ---------------------------------------------------------------------------
# Cell-type ratio values for the Y axis (Barlow et al. Fig 2F)
# ---------------------------------------------------------------------------
#
# Their y-axis for cell-type abundance is log(# cells / # endocrine cells), NOT
# a raw proportion and NOT a marker mean. Dividing by the islet's own endocrine
# count is what makes islets of different sizes comparable -- a raw count says
# more about how big the islet is than about how infiltrated it is.
#
# WHAT THIS DOES AND DOES NOT FIX. It removes the size confound and it spreads
# the non-zero values out. It does NOT remove the zero pile-up: log1p(0) = 0, so
# an islet with no CD8a T cell sits at exactly 0 whatever transform is applied.
# That is a property of the data, not of the scale -- Barlow's own Fig 2F has it
# too. What genuinely reduces it is sampling more cells per islet, which is why
# the region default is core+peri and why the min-cells filter matters:
#
#   CD8a T cell, fraction of islets at exactly zero
#     core only            0.956
#     core + peri          0.864
#     core + peri, >=10    0.814
#     core + peri, >=25    0.718
#
# `traj_ratio_zero_fraction()` surfaces this per selection so a flat line of
# points at zero reads as "these cells are absent", not as a broken axis.

BARLOW_RATIO_REGIONS <- c("Core + Peri" = "both", "Core only" = "core", "Peri only" = "peri")

#' Cell types available as ratio features, in a stable display order.
barlow_ratio_types <- function() {
  d <- barlow_features()
  if (is.null(d)) return(character(0))
  all <- sub("^n_core_", "", grep("^n_core_", names(d), value = TRUE))
  # The nine retired phenotypes stay out of every selector -- they are
  # gate-failure residuals and `-rest` catch-alls whose abundance tracks gate
  # quality, not biology. `phenotype_keep()` is the single gate for that
  # (00_globals.R); the CSV still carries them because `Endocrine` is the RATIO
  # DENOMINATOR even though it is not selectable as a numerator.
  sort(phenotype_keep(all))
}

#' log1p(cells of `type` / endocrine cells in the islet core), per islet.
#'
#' Returns a data.frame(islet_id, value) or NULL. log1p rather than log so an
#' absent cell type maps to 0 instead of -Inf; with ~86% of islets at zero for
#' the rarer types, plain log would drop most of the data off the plot.
barlow_ratio_values <- function(type, region = "both") {
  d <- barlow_features()
  if (is.null(d) || !nzchar(type %||% "")) return(NULL)
  cc <- paste0("n_core_", type); pc <- paste0("n_peri_", type)
  if (!all(c(cc, pc) %in% names(d))) return(NULL)
  n_core <- suppressWarnings(as.numeric(d[[cc]]))
  n_peri <- suppressWarnings(as.numeric(d[[pc]]))
  n_core[!is.finite(n_core)] <- 0
  n_peri[!is.finite(n_peri)] <- 0
  n <- switch(region, core = n_core, peri = n_peri, n_core + n_peri)
  # Denominator is the islet's own endocrine count, floored at 1 -- the same
  # denominator the PAGA features use, so the axis and the embedding agree.
  denom <- pmax(suppressWarnings(as.numeric(d$n_endocrine_core)), 1)
  data.frame(islet_id = as.character(d$islet_id),
             value = log1p(n / denom),
             stringsAsFactors = FALSE)
}

#' Fraction of islets (optionally above a min-cells threshold) at exactly zero.
traj_ratio_zero_fraction <- function(type, region = "both", min_cells = 1, max_cells = NA) {
  v <- barlow_ratio_values(type, region)
  if (is.null(v)) return(NA_real_)
  d <- barlow_features()
  n <- suppressWarnings(as.numeric(d$total_cells_core))
  keep <- n >= max(1, min_cells)
  if (is.finite(max_cells)) keep <- keep & n <= max_cells
  x <- v$value[keep]
  if (!length(x)) return(NA_real_)
  mean(x == 0)
}
