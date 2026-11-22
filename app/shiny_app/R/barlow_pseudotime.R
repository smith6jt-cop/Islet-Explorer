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
#                                      -> the app fits a new graph for its selected
#                                         core-cell range (default four or more)
#   "pseudotemporal distance from the centroid of non-T1D islets"
#                                      -> app adaptation uses a donor-supported ND peripheral root
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

# ---------------------------------------------------------------------------
# Graph inputs: the current typing format (owner, 2026-11-18)
# ---------------------------------------------------------------------------
# "The trajectory/pseudotime needs to be updated to the current typing format.
#  There should not be peri-islet features for the endocrine classes, core
#  exocrine, or 12 macrophage subtypes. Also the qc classes should not be
#  entered into the pseudotime pipeline."
#
#   endocrine classes   core only
#   exocrine classes    peri-islet band only
#   macrophages         the two Gal3 subclasses (macrophage_profile()); the twelve
#                       marker-profile labels of an older root are summed by Gal3 sign
#   QC classes          never: Conflict, Doublet suspected, Tissue, Unclassified and
#                       the unresolved lineage residuals (names(UNCLASSIFIED_RAW_TO_KEY))
#   every other named class   core and peri-islet band
#
# python/paga_pseudotime.py applies the same rules for the scanpy path; a pipeline
# test compares its label sets with the ones used here.
TRAJECTORY_INPUT_POLICY <- "phenotype-passthrough/2026-11-19"
TRAJECTORY_ENDOCRINE_CLASSES <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
# The uniform roots' Exocrine, its Acinar / Ductal split in the November revision, the frozen typing's Ductal.
TRAJECTORY_EXOCRINE_CLASSES <- c("Exocrine", "Acinar cell", "Ductal cell", "Ductal")

#' The current-format class of a data label: a macrophage label becomes its Gal3 subclass.
trajectory_class <- function(labels) macrophage_profile(as.character(labels))

#' Whether a class counted in a region ("core" or "peri") is a graph input.
trajectory_feature_allowed <- function(region, labels) {
  cls <- trajectory_class(labels)
  qc <- phenotype_is_conflict(cls) | cls %in% c(names(UNCLASSIFIED_RAW_TO_KEY),
    "Panel undected", "Conflicting markers - excluded")
  !qc & !(cls %in% PHENOTYPE_RETIRED) &
    (!(cls %in% TRAJECTORY_ENDOCRINE_CLASSES) | region == "core") &
    (!(cls %in% TRAJECTORY_EXOCRINE_CLASSES) | region == "peri")
}

#' The stored core_/peri_ ratio columns that feed the graph, in file order. Twelve-label macrophage columns
#' are returned as stored; barlow_trajectory_inputs() sums them by Gal3 sign.
barlow_trajectory_features <- function(columns) {
  features <- grep("^(core|peri)_", columns, value = TRUE)
  features[trajectory_feature_allowed(sub("_.*$", "", features), sub("^(core|peri)_", "", features))]
}

#' The graph-input matrix of a feature table, one named column per class and region.
#'
#' A class with one stored column is used as stored, so its values and the stored endocrine denominator are
#' unchanged. A class with several source labels (the pooled macrophages) is log1p(summed count / endocrine
#' cells) from the n_core_ / n_peri_ counts: the log of a sum is not the sum of the logs.
barlow_phenotype_inputs <- function(d) {
  columns <- barlow_trajectory_features(names(d))
  cls <- trajectory_class(sub("^(core|peri)_", "", columns))
  region <- sub("_.*$", "", columns)
  region[!cls %in% c(TRAJECTORY_ENDOCRINE_CLASSES, TRAJECTORY_EXOCRINE_CLASSES)] <- "combined"
  name <- paste0(region, "_", cls)
  features <- unique(name)
  num <- function(k) suppressWarnings(as.numeric(d[[k]]))
  values <- lapply(features, function(f) {
    src <- columns[name == f]
    if (identical(src, f)) return(num(f))
    if (startsWith(f, "combined_")) {
      labels <- unique(sub("^(core|peri)_", "", src))
      src <- c(paste0("core_", labels), paste0("peri_", labels))
    }
    counts <- paste0("n_", src)
    missing <- setdiff(c(counts, "n_endocrine_core"), names(d))
    if (length(missing)) {
      stop(f, " pools ", length(src), " source labels and needs their counts; missing ",
           paste(missing, collapse = ", "))
    }
    log1p(Reduce(`+`, lapply(counts, num)) / pmax(num("n_endocrine_core"), 1))
  })
  M <- if (length(values)) do.call(cbind, values) else matrix(numeric(0), nrow = nrow(d), ncol = 0)
  colnames(M) <- features
  M
}

barlow_trajectory_inputs <- function(d) {
  counts <- barlow_phenotype_inputs(d)
  keys <- unlist(lapply(colnames(counts), function(name)
    paste0(name, "__", c("nerve", "lymphatic", "vessel"))), use.names = FALSE)
  columns <- paste0("overlap__", keys)
  if (!all(columns %in% names(d))) stop("Missing phenotype-specific passthrough measurements")
  raw <- do.call(cbind, lapply(columns, function(k) as.numeric(d[[k]])))
  if (any(!is.finite(raw) | raw < 0)) stop("Invalid phenotype-specific passthrough totals")
  colnames(raw) <- keys
  cbind(counts, log1p(raw / pmax(as.numeric(d[["n_endocrine_core"]]), 1)))
}

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

BARLOW_ROOT_METHOD <- "Donor-supported ND root at a multiscale diffusion endpoint"

barlow_root_distances <- function(ev, eigenvalues, root) {
  # destiny's DPT[root, ] uses these multiscale distances. Direct evaluation
  # avoids fitting the branch tree, which the trajectory tab never requests.
  delta <- sweep(ev, 2, ev[root, ], "-")
  sqrt(rowSums(sweep(delta, 2, eigenvalues / (1 - eigenvalues), "*")^2))
}

barlow_peripheral_root <- function(dm, eigenvalues, back, sub) {
  keep <- is.finite(eigenvalues) & eigenvalues < 0.9994
  coords <- sweep(dm[, keep, drop = FALSE], 2, eigenvalues[keep] / (1 - eigenvalues[keep]), "*")
  if (!ncol(coords) || !all(is.finite(coords))) stop("No finite nonstationary diffusion coordinates")
  nd <- which(sub$donor_status == "ND")
  if (!length(nd)) stop("No ND objects for root selection")
  # Search the diameter in blocks to bound memory when small objects are included.
  n <- nrow(coords)
  norm2 <- rowSums(coords^2)
  best <- 0
  ends <- NULL
  for (start in seq.int(1L, n, by = 128L)) {
    rows <- seq.int(start, min(n, start + 127L))
    d2 <- outer(norm2, norm2[rows], "+") - 2 * tcrossprod(coords, coords[rows, , drop = FALSE])
    pos <- arrayInd(which.max(d2), dim(d2))
    if (d2[pos] > best) {
      best <- d2[pos]
      ends <- sort(c(pos[1L], rows[pos[2L]]))
    }
  }
  if (is.null(ends) || ends[1L] == ends[2L]) stop("No distinct diffusion endpoints")
  distances <- lapply(ends, function(e) sqrt(rowSums(sweep(coords, 2, coords[e, ], "-")^2)))
  k <- max(1L, ceiling(0.05 * nrow(coords)))
  nearby <- lapply(distances, function(d) head(order(d, seq_len(n)), k))
  local_nd <- lapply(nearby, function(ids) nd[back[nd] %in% ids])
  support <- vapply(local_nd, function(i) length(unique(sub$case_id[i])), integer(1))
  proximity <- vapply(distances, function(d) {
    stats::median(tapply(d[back[nd]], sub$case_id[nd], min))
  }, numeric(1))
  side <- order(-support, proximity, ends)[1]
  idx <- local_nd[[side]]
  if (!length(idx)) {
    idx <- nd
    center <- coords[ends[side], ]
  } else {
    local <- unique(data.frame(profile = back, donor = sub$case_id)[back %in% nearby[[side]], ])
    donor_centers <- lapply(split(local$profile, local$donor), function(i) colMeans(coords[i, , drop = FALSE]))
    center <- colMeans(do.call(rbind, donor_centers))
  }
  distance <- rowSums(sweep(coords[back[idx], , drop = FALSE], 2, center)^2)
  row <- idx[order(distance, as.character(sub$islet_id[idx]), method = "radix")[1]]
  list(profile = back[row], islet_id = as.character(sub$islet_id[row]), donor = as.character(sub$case_id[row]))
}

#' Refit the destiny fallback for the selected range and cache the object-level DPT.
compute_barlow_pseudotime <- function(min_cells = TRAJ_MIN_CELLS, max_cells = NA, object_classes = NULL) {
  mc <- suppressWarnings(as.numeric(min_cells))
  if (!is.finite(mc) || mc < 1) mc <- 1
  xc <- suppressWarnings(as.numeric(max_cells))
  key <- paste0("pt_", mc, "_", if (is.finite(xc)) xc else "inf", object_class_key(object_classes))
  if (!is.null(.barlow_cache[[key]])) return(.barlow_cache[[key]])

  d <- barlow_features()
  if (is.null(d) || !nrow(d)) return(NULL)
  keep <- d$total_cells_core >= mc
  if (is.finite(xc)) keep <- keep & d$total_cells_core <= xc
  if (!is.null(object_classes) && "object_class" %in% names(d)) {
    keep <- keep & d$object_class %in% object_classes
  }
  sub <- d[keep, , drop = FALSE]
  # A diffusion map on a handful of points is meaningless; fail closed.
  if (nrow(sub) < 50) return(NULL)

  # The same graph inputs as the scanpy path (current typing format, above). Stored ratios keep their
  # endocrine denominator; only the choice and pooling of input columns is made here.
  M <- tryCatch(barlow_trajectory_inputs(sub), error = function(e) NULL)
  if (is.null(M)) {
    islet_log_event("barlow_pseudotime_failed")
    return(NULL)
  }
  # An object with no cell of any input class has no make-up to place; it gets no pseudotime.
  placed <- rowSums(M != 0, na.rm = TRUE) > 0
  n_without_inputs <- sum(!placed)
  if (n_without_inputs) {
    sub <- sub[placed, , drop = FALSE]
    M <- M[placed, , drop = FALSE]
    if (nrow(sub) < 50) return(NULL)
  }
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
    set.seed(1)  # destiny's eigendecomposition is randomized; pin it so the
                 # axis does not move between two runs at the same threshold.
    dm <- destiny::DiffusionMap(Mu, k = 15, n_eigs = 10,
                                distance = "cosine", verbose = FALSE)
    ev <- as.matrix(destiny::eigenvectors(dm))
    if (nrow(ev) != nrow(Mu)) stop("destiny returned ", nrow(ev),
                                   " eigenvector rows for ", nrow(Mu), " inputs")
    root <- barlow_peripheral_root(ev, destiny::eigenvalues(dm), back, sub)
    ptu <- barlow_root_distances(ev, destiny::eigenvalues(dm), root$profile)
    pt <- ptu[back]
    root_id <- root$islet_id
    root_donor <- root$donor

    data.frame(
      islet_id     = as.character(sub$islet_id),
      case_id      = as.character(sub$case_id),
      donor_status = as.character(sub$donor_status),
      pseudotime   = as.numeric(pt),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    islet_log_event("barlow_pseudotime_failed")
    NULL
  })
  if (!is.null(res)) {
    attr(res, "root_method") <- BARLOW_ROOT_METHOD
    attr(res, "root_donor") <- root_donor
    attr(res, "root_islet_id") <- root_id
    attr(res, "features") <- colnames(M)
    attr(res, "n_without_inputs") <- n_without_inputs
  }

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

# Endocrine subtypes describe core composition. The all-endocrine core
# numerator duplicates the denominator. Keep these display rules local to
# trajectory ratios; the underlying compartment counts remain available.
barlow_ratio_choices <- function(types, region) {
  if (length(region) != 1L || !region %in% BARLOW_RATIO_REGIONS) return(character())
  excluded <- "All endocrine cells"
  if (region != "core") excluded <- c(excluded, TRAJECTORY_ENDOCRINE_CLASSES, "Endocrine")
  setdiff(types, excluded)
}

#' Cell types available as ratio features, in a stable display order.
barlow_ratio_types <- function(markers = character(0), region = "core") {
  d <- barlow_features()
  if (is.null(d)) return(character(0))
  all <- sub("^n_core_", "", grep("^n_core_", names(d), value = TRUE))
  # The nine retired phenotypes stay out of every selector -- they are
  # gate-failure residuals and `-rest` catch-alls whose abundance tracks gate
  # quality, not biology. `phenotype_keep()` is the single gate for that
  # (00_globals.R); the CSV still carries them because `Endocrine` is the RATIO
  # DENOMINATOR even though it is not selectable as a numerator.
  # The macrophages merged for the checked markers (macrophage_collapse(); owner, 2026-10-27).
  kept <- sort(unique(phenotype_collapse(phenotype_keep(all), markers)))
  # The group aggregates ("All macrophages"; owner, 2026-10-27), summed from all the members' counts, general
  # classes ("T cell", "Endocrine") included.
  types <- c(phenotype_aggregates_for(unique(phenotype_collapse(all, markers))), kept)
  barlow_ratio_choices(types, region)
}

#' log1p(cells of `type` / endocrine cells in the islet core), per islet.
#'
#' Returns a data.frame(islet_id, value) or NULL. log1p rather than log so an
#' absent cell type maps to 0 instead of -Inf; with ~86% of islets at zero for
#' the rarer types, plain log would drop most of the data off the plot.
barlow_ratio_values <- function(type, region = "both") {
  d <- barlow_features()
  if (is.null(d) || length(type) != 1L || is.na(type) || !nzchar(type) ||
      !length(barlow_ratio_choices(type, region))) return(NULL)
  # an aggregate or a merged macrophage category sums its data labels (phenotype_members())
  types <- intersect(phenotype_members(type), sub("^n_core_", "", grep("^n_core_", names(d), value = TRUE)))
  cc <- paste0("n_core_", types); pc <- paste0("n_peri_", types)
  if (!length(types) || !all(c(cc, pc) %in% names(d))) return(NULL)
  count <- function(cols) {
    m <- vapply(cols, function(k) { v <- suppressWarnings(as.numeric(d[[k]])); v[!is.finite(v)] <- 0; v },
                numeric(nrow(d)))
    if (is.matrix(m)) rowSums(m) else sum(m)
  }
  n_core <- count(cc)
  n_peri <- count(pc)
  n <- switch(region, core = n_core, peri = n_peri, n_core + n_peri)
  # Denominator is the islet's own endocrine count, floored at 1 -- the same
  # denominator the PAGA features use, so the axis and the embedding agree.
  denom <- pmax(suppressWarnings(as.numeric(d$n_endocrine_core)), 1)
  data.frame(islet_id = as.character(d$islet_id),
             value = log1p(n / denom),
             stringsAsFactors = FALSE)
}

#' Fraction of islets (optionally above a min-cells threshold) at exactly zero.
traj_ratio_zero_fraction <- function(type, region = "both", min_cells = TRAJ_MIN_CELLS, max_cells = NA) {
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
