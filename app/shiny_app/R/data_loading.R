# ---------- Data loading and wrangling ----------
# Extracted from app.R
# Dependencies: safe_left_join, add_islet_key, compute_diameter_um (from R/utils_safe_join.R)
#               master_path (from global.R)

safe_read_sheet <- function(path, sheet) {
  # Increase guess_max to reduce type misguesses on sparse columns
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
}

load_master <- function(path = master_path) {
  list(
    markers = safe_read_sheet(path, "Islet_Markers"),
    targets = safe_read_sheet(path, "Islet_Targets"),
    comp    = safe_read_sheet(path, "Islet_Composition"),
    lgals3  = safe_read_sheet(path, "LGALS3")
  )
}

# OPTIMIZATION: Consolidate donor metadata extraction
get_donor_metadata <- function(master) {
  # Priority order: composition > targets > markers
  for (sheet in c("comp", "targets", "markers")) {
    df <- master[[sheet]]
    if (!is.null(df) && nrow(df) > 0) {
      required_cols <- c("Case ID", "Donor Status")

      # Determine which AAb columns exist
      aab_candidates <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
      aab_cols <- intersect(aab_candidates, names(df))

      if (all(required_cols %in% names(df))) {
        return(df %>%
          dplyr::select(dplyr::all_of(c(required_cols, aab_cols))) %>%
          dplyr::distinct())
      }
    }
  }
  NULL
}

prep_data <- function(master) {
  # Determine which autoantibody columns are available
  aab_cols_targets <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$targets))
  aab_cols_markers <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$markers))
  aab_cols_comp    <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$comp))
  aab_cols_all     <- unique(c(aab_cols_targets, aab_cols_markers, aab_cols_comp))
  # Use new consolidated function
  donors_meta <- get_donor_metadata(master)
  # Islet size proxy per islet corefor diameter
  targets <- master$targets %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  core_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_core") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Diameter is computed ONLY from core area; if no core exists, diameter is NA
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)

  # Targets: keep all region types and later filter by user selection
  targets_all <- targets %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_targets), islet_key, type, class, area_um2, region_um2, area_density, count) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        tolower(type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(type)
      )
      # Density is already in um2, no conversion needed
    ) %>%
    { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_all:size_area") }

  # Synthesize missing union rows for counts by summing core + band (do NOT fabricate union area)
  if (nrow(targets_all) > 0) {
    keys <- c("Case ID", "Donor Status", "islet_key", "class")
    # Identify which key combos already have a union row
    have_union <- targets_all %>%
      dplyr::filter(type == "islet_union") %>%
      dplyr::select(dplyr::all_of(keys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_core = count)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_band = count)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
  { safe_left_join(., have_union, by = keys, context = "targets_union_missing:have_union") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = NA_real_, region_um2 = NA_real_,
                          area_density = NA_real_,
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band)))
      # Attach diameter via size_area
  synth <- synth %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_synth:size_area") }
      # Bind synthetic union rows to targets_all
      targets_all <- dplyr::bind_rows(targets_all, synth)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones
  if (!is.null(donors_meta)) {
    targets_all <- targets_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "targets_all:donors_meta") }
  }

  # Markers with fraction positive / mean intensity (include LGALS3 sheet)
  markers <- master$markers %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  markers_all <- markers %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_markers), islet_key, region_type, marker, n_cells, pos_count, pos_frac) %>%
    dplyr::mutate(
      region_type = dplyr::case_when(
        tolower(region_type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(region_type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(region_type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(region_type)
      )
    )
  # Add LGALS3 rows if available
  if (!is.null(master$lgals3) && nrow(master$lgals3) > 0) {
    g3 <- master$lgals3 %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key)) %>%
      mutate(marker = as.character(marker)) %>%
      dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(intersect(colnames(master$lgals3), c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"))), islet_key, region_type, marker, n_cells, pos_count, pos_frac)
    if (!is.null(g3) && nrow(g3) > 0) {
      markers_all <- bind_rows(markers_all, g3)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones (markers)
  if (!is.null(donors_meta)) {
    markers_all <- markers_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "markers_all:donors_meta") }
  }
  markers_all <- markers_all %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_all:size_area") }

  # Synthesize missing union rows for markers by summing core + band counts (do NOT fabricate area)
  if (nrow(markers_all) > 0) {
    mkeys <- c("Case ID", "Donor Status", "islet_key", "marker")
    have_union_m <- markers_all %>%
      dplyr::filter(region_type == "islet_union") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_m <- markers_all %>% dplyr::filter(region_type == "islet_core") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_core = n_cells, pos_core = pos_count)
    band_m <- markers_all %>% dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>%
      dplyr::rename(n_band = n_cells, pos_band = pos_count)
    union_m <- markers_all %>% dplyr::filter(region_type == "islet_union") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_union = n_cells, pos_union = pos_count)
    union_missing_m <- core_m %>%
      dplyr::inner_join(band_m, by = mkeys) %>%
  { safe_left_join(., have_union_m, by = mkeys, context = "markers_union_missing:have_union_m") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing_m) > 0) {
      synth_m <- union_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_union", marker,
                          n_cells = suppressWarnings(as.numeric(n_core)) + suppressWarnings(as.numeric(n_band)),
                          pos_count = suppressWarnings(as.numeric(pos_core)) + suppressWarnings(as.numeric(pos_band))) %>%
        dplyr::mutate(pos_frac = ifelse(is.finite(n_cells) & n_cells > 0,
                                        suppressWarnings(as.numeric(pos_count)) / suppressWarnings(as.numeric(n_cells)),
                                        NA_real_)) %>%
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_m:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_m)
    }

    # Backfill missing band rows when core and union exist but band is missing: band = union - core (counts only)
    have_band_m <- markers_all %>%
      dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_band = TRUE)
    band_missing_m <- core_m %>%
      dplyr::inner_join(union_m, by = mkeys) %>%
  { safe_left_join(., have_band_m, by = mkeys, context = "markers_band_missing:have_band_m") } %>%
      dplyr::filter(is.na(.has_band))
    if (nrow(band_missing_m) > 0) {
      synth_band <- band_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_band", marker,
                          n_cells = suppressWarnings(as.numeric(n_union)) - suppressWarnings(as.numeric(n_core)),
                          pos_count = suppressWarnings(as.numeric(pos_union)) - suppressWarnings(as.numeric(pos_core))) %>%
        dplyr::mutate(
          n_cells = ifelse(is.finite(n_cells), n_cells, 0),
          pos_count = ifelse(is.finite(pos_count), pos_count, 0),
          n_cells = pmax(0, n_cells),
          pos_count = pmax(0, pos_count),
          pos_frac = ifelse(n_cells > 0, pos_count / n_cells, NA_real_)
        ) %>%
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_band:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_band)
    }
  }

  # Fill pos_frac when n_cells and pos_count are present but pos_frac is NA; also compute pos_pct (0..100)
  if (nrow(markers_all) > 0) {
    markers_all <- markers_all %>%
      dplyr::mutate(
        .n = suppressWarnings(as.numeric(n_cells)),
        .p = suppressWarnings(as.numeric(pos_count)),
        pos_frac = dplyr::coalesce(pos_frac, ifelse(is.finite(.n) & .n > 0 & is.finite(.p), .p/.n, NA_real_)),
        pos_pct = ifelse(is.finite(pos_frac), 100.0 * pos_frac, NA_real_)
      ) %>%
      dplyr::select(-.n, -.p)
  }

  # Composition by islet
  comp <- master$comp %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  comp <- comp %>% dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_comp), islet_key, cells_total, Ins_single, Glu_single, Stt_single,
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any)
  comp <- comp %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "comp:size_area") }

  message("[prep-final] size_area=", paste(class(size_area), collapse='/'),
    " targets_all=", paste(class(targets_all), collapse='/'),
    " markers_all=", paste(class(markers_all), collapse='/'),
    " comp=", paste(class(comp), collapse='/'))
  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all, comp = comp)
}

# Simple NA audit
audit_na <- function(df, label) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  na_cnt <- vapply(df, function(x) sum(is.na(x)), integer(1))
  total <- nrow(df)
  pct <- ifelse(total > 0, round(100 * na_cnt / total, 2), 0)
  msg <- paste0("[NA audit] ", label, ": ", paste(names(na_cnt), paste0(na_cnt, " (", pct, "%)"), sep = "=", collapse = "; "))
  message(msg)
}

bin_islet_sizes <- function(df, diam_col, width) {
  x <- suppressWarnings(as.numeric(df[[diam_col]]))
  max_x <- max(x, na.rm = TRUE)
  max_x <- ifelse(is.finite(max_x), max_x, 0)
  # Define numeric bins: [lo, hi)
  bin_lo <- floor(x / width) * width
  bin_hi <- bin_lo + width
  diam_mid <- bin_lo + width/2
  # Build a human-readable label
  diam_bin <- paste0("[", bin_lo, ", ", bin_hi, ")")
  df$diam_bin <- factor(diam_bin, levels = unique(diam_bin[order(bin_lo)]), ordered = TRUE)
  df$diam_mid <- diam_mid
  df
}

#' Compute Lamian-based pseudotime from AnnData
#' Uses Lamian's infer_tree_structure which accounts for multi-sample design
#' @param ad AnnData object (from anndata package)
#' @param features Character vector of feature names to use (NULL = use all)
#' @return Numeric vector of pseudotime values (same length as cells)
compute_lamian_pseudotime <- function(ad, features = NULL) {
  tryCatch({
    require(Lamian, quietly = TRUE)

    cat("Computing Lamian pseudotime with infer_tree_structure...\n")

    # Extract expression matrix (cells x genes)
    expr_mat <- as.matrix(ad$X)
    rownames(expr_mat) <- ad$obs_names
    colnames(expr_mat) <- ad$var_names

    # Define curated feature set for trajectory inference
    # Focus on hormones, immune markers, and spatial features that drive disease progression
    curated_features <- c(
      # Hormone markers (islet cell types)
      "INS", "GCG",
      # Immune infiltration markers
      "CD8a", "CD4", "HLADR", "CD163", "CD68",
      # Disease-associated markers
      "LGALS3", "BCatenin",
      # Spatial features (microenvironment)
      "Dist to Closest Lymphatic", "Dist to Closest Capillary", "Dist to Closest Nerve"
    )

    # Subset to requested features if provided, otherwise use curated set
    if (!is.null(features) && length(features) > 0) {
      available_features <- intersect(features, ad$var_names)
      if (length(available_features) == 0) {
        warning("None of the requested features found in AnnData, using curated features")
        available_features <- intersect(curated_features, ad$var_names)
      } else {
        cat(sprintf("Using %d of %d requested features\n", length(available_features), length(features)))
      }
    } else {
      # Use curated feature set by default
      available_features <- intersect(curated_features, ad$var_names)
      cat(sprintf("Using curated feature set: %d features\n", length(available_features)))
      cat(sprintf("  Features: %s\n", paste(available_features, collapse=", ")))
    }

    if (length(available_features) > 0) {
      feature_idx <- match(available_features, ad$var_names)
      expr_mat <- expr_mat[, feature_idx, drop = FALSE]
      colnames(expr_mat) <- ad$var_names[feature_idx]
    } else {
      stop("No valid features found for trajectory inference")
    }

    cat("Running PCA for dimensionality reduction...\n")
    # Compute PCA (Lamian expects cells x PCs)
    # Data is already z-scored, so no need for redundant centering/scaling
    n_pcs <- min(length(available_features) - 1, nrow(expr_mat) - 1)
    pca_result <- prcomp(expr_mat, rank. = n_pcs, center = FALSE, scale. = FALSE)
    pca_coords <- pca_result$x

    # Report variance explained
    var_explained <- summary(pca_result)$importance[2, ]
    cumvar <- cumsum(var_explained)
    n_pcs_80 <- if (any(cumvar >= 0.80)) which(cumvar >= 0.80)[1] else n_pcs
    cat(sprintf("  Using %d PCs (explains %.1f%% variance, %d PCs for 80%%)\n",
                n_pcs, cumvar[n_pcs] * 100, n_pcs_80))

    cat("Building cell annotation with sample IDs...\n")
    # Build cell annotation - CRITICAL: column 2 must be sample/donor ID
    # This allows Lamian to account for multi-sample structure
    cellanno <- data.frame(
      cell = ad$obs_names,
      sample = as.character(ad$obs$imageid),  # Donor/sample ID
      stringsAsFactors = FALSE
    )

    # Add cell type if available (optional)
    if ("donor_status" %in% colnames(ad$obs)) {
      cellanno$celltype <- as.character(ad$obs$donor_status)
    }

    cat("Running Lamian trajectory inference...\n")
    cat("  This accounts for donor/sample structure to avoid artificial grouping\n")
    cat("  Origin marker: INS (clusters with highest mean INS expression)\n")

    # Use Lamian's infer_tree_structure (designed for multi-sample data)
    # Transpose expression for Lamian (genes x cells expected)
    expr_mat_t <- t(expr_mat)

    # Determine reasonable max cluster number based on data size
    n_obs <- nrow(ad$obs)
    n_samples <- length(unique(cellanno$sample))
    # Use more clusters for better resolution: ~sqrt(n) or n/50, capped at 50
    max_clusters <- min(50, max(10, ceiling(sqrt(n_obs)), ceiling(n_obs / 50)))

    cat(sprintf("  Max clusters: %d (based on %d observations, %d samples)\n",
                max_clusters, n_obs, n_samples))

    res <- Lamian::infer_tree_structure(
      pca = pca_coords,
      cellanno = cellanno,
      expression = expr_mat_t,
      origin.marker = "INS",           # Use insulin as trajectory root - Lamian finds cluster with highest mean INS
      number.cluster = NA,              # Auto-determine cluster number
      max.clunum = max_clusters,        # Maximum clusters to consider (adaptive)
      kmeans.seed = 12345              # Reproducible clustering
    )

    # Extract pseudotime from Lamian result
    pseudotime <- res$pseudotime

    # Report what Lamian found
    if (!is.null(res$clusterRes)) {
      n_clusters <- length(unique(res$clusterRes))
      cat(sprintf("  Lamian identified %d clusters\n", n_clusters))
    }

    # Ensure it's a vector aligned with cells
    if (is.matrix(pseudotime)) {
      pseudotime <- as.numeric(pseudotime[, 1])
    }

    # Handle names if present
    if (!is.null(names(pseudotime))) {
      # Reorder to match ad$obs_names
      pseudotime <- pseudotime[ad$obs_names]
    }

    # Check pseudotime distribution before normalization
    cat(sprintf("  Raw pseudotime range: %.3f - %.3f (mean: %.3f, sd: %.3f)\n",
                min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE),
                mean(pseudotime, na.rm = TRUE), sd(pseudotime, na.rm = TRUE)))

    # Normalize to 0-1 range for consistency with PAGA
    pseudotime <- (pseudotime - min(pseudotime, na.rm = TRUE)) /
                  (max(pseudotime, na.rm = TRUE) - min(pseudotime, na.rm = TRUE))

    # REVERSE pseudotime direction (biological progression ND->Aab+->T1D)
    # Lamian origin.marker='INS' finds high-INS clusters, but trajectory may run backward
    # Reversal ensures: 0=healthy (high INS), 1=diseased (low INS)
    pseudotime <- 1 - pseudotime
    cat("  Pseudotime reversed for biological progression (0=healthy, 1=diseased)\n")

    cat(sprintf("Lamian pseudotime computed (range: %.3f - %.3f)\n",
                min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE)))
    cat(sprintf("  %d cells across %d samples\n",
                length(pseudotime), length(unique(cellanno$sample))))

    return(pseudotime)

  }, error = function(e) {
    warning(sprintf("Lamian pseudotime calculation failed: %s", e$message))
    cat("Error details:", e$message, "\n")
    return(rep(NA_real_, nrow(ad$obs)))
  })
}
