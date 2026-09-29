# Pseudotime x-axis extent, derived in ONE place so the scatter and the strip heatmap
# stacked directly beneath it always share an x-extent. The scatter jitters points by
# ±1.2% of the domain and widens its axis by that much rather than clamping (clamping
# piled the outermost ~1.2% of points onto the endpoints), so the heatmap must be widened
# by the same margin or the two panels no longer line up.
#
# Returns:
#   dom    — the data domain, used for BINNING (bins must cover data, not padding)
#   jitter — the jitter half-width
#   lim    — dom expanded by `jitter`, used for AXIS LIMITS in both panels
traj_pt_extent <- function(pt, jitter_frac = 0.012) {
  # suppressWarnings: an all-NA pt column makes range() emit "no non-missing arguments"
  # and return c(Inf, -Inf); the is.finite guard below is the intended handling.
  dom <- suppressWarnings(range(pt, na.rm = TRUE))
  if (!all(is.finite(dom)) || diff(dom) <= 0) dom <- c(0, 1)
  j <- diff(dom) * jitter_frac
  list(dom = dom, jitter = j, lim = c(dom[1] - j, dom[2] + j))
}

# Resolve the plotted pseudotime column from the selected scale.
#
# "rank" (default) spreads islets uniformly, counteracting the raw-DPT compression (95% of
# islets fall between 0.25-0.70 due to diffusion-manifold geometry). "raw" shows true DPT
# spacing, so where islets actually concentrate is visible — note that on the rank axis LOESS
# curvature and bin density are artifacts of the transform, since ranks are uniform by
# construction.
#
# `pt_rank` is precomputed in traj_with_pseudotime() over the complete cell-range cohort, so
# the axis does not shift with the selected Y feature. Never re-rank the feature-complete
# subset here — that was the original bug (mean |delta rank| 0.179, max 0.357 at min_cells=25).
#
# Shared by traj_data_clean() and traj_multi_feature_data(), which are two independent
# extraction paths; a divergence here silently desynchronizes the scatter from the
# multi-feature heatmap's x-axis.
# The ONE pseudotime: an islet-composition diffusion map computed on whatever
# islets the current cell-count range selects. Islets outside the range get NA,
# which is what removes them downstream -- the app's ordinary is.finite() filters
# do it, so the axis and the plotted point set cannot disagree.
#
# The range is part of the DEFINITION, not a display filter: moving it refits the
# embedding. That is deliberate and it is announced in the UI. It does not
# reintroduce the axis-instability defect (docs/trajectory.md), which was the
# axis moving as a side effect of changing the Y-axis FEATURE -- that still
# cannot happen, and is regression-tested.
traj_pseudotime_column <- function(obs, cell_range = c(1, Inf)) {
  lo <- cell_range[1]; hi <- cell_range[2]
  # scanpy DPT is primary. The destiny implementation is the fallback if scanpy
  # is unreachable from reticulate. PAGA/ForceAtlas2 is the two-dimensional graph
  # card; it is not the pseudotime algorithm itself.
  bp <- paga_result(lo, hi)
  engine <- "scanpy diffusion pseudotime"
  if (is.null(bp)) {
    bp <- compute_barlow_pseudotime(lo, hi)
    engine <- "destiny diffusion-pseudotime fallback"
  }
  if (is.null(bp)) {
    out <- rep(NA_real_, nrow(obs))
    attr(out, "engine") <- "unavailable"
    attr(out, "n_islets") <- 0L
    return(out)
  }
  ids <- as.character(obs$combined_islet_id)
  if (all(is.na(ids)) || !any(ids %in% bp$islet_id)) {
    ids <- paste0(as.character(obs$imageid), "_",
                  sub("^Islet_Islet_", "Islet_", as.character(obs$base_islet_id)))
  }
  out <- bp$pseudotime[match(ids, bp$islet_id)]
  attr(out, "engine") <- engine
  attr(out, "n_islets") <- sum(is.finite(out))
  out
}

apply_pt_scale <- function(df, pt_scale) {
  df$pt_raw <- df$pt
  if (identical(pt_scale, "raw")) {
    df$pt <- df$pt_raw
  } else if (is.null(df$pt_rank) || all(is.na(df$pt_rank))) {
    # Defensive: older H5AD/pseudotime paths that never populated pt_rank.
    df$pt <- (rank(df$pt_raw, ties.method = "average") - 0.5) / nrow(df)
  } else {
    df$pt <- df$pt_rank
  }
  df
}

trajectory_server <- function(id, prepared, selected_islet, forced_image, active_tab = reactive("Trajectory"),
                              donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                              remove_outliers = reactive(FALSE),
                              outlier_threshold = reactive(3.0)) {
  # `remove_outliers` / `outlier_threshold` are accepted but UNUSED: app.R still
  # passes them to every analysis module, and changing that signature is a wider
  # edit than this tab. The Trajectory tab no longer exposes or applies them.
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- H5AD Path Resolution ----------
    resolve_traj_path <- function() {
      cand <- c(
        data_path("adata_ins_root.h5ad"),
        file.path("..", "..", "scripts", "adata_ins_root.h5ad"),
        if (!is.null(project_root)) file.path(project_root, "data", "adata_ins_root.h5ad") else NULL,
        if (!is.null(project_root)) file.path(project_root, "scripts", "adata_ins_root.h5ad") else NULL,
        Sys.getenv("ADATA_PATH", unset = "")
      )
      cand <- cand[nzchar(cand)]
      for (p in cand) if (file.exists(p)) return(p)
      NULL
    }

    traj_path <- resolve_traj_path()

    # ---------- Reactive state ----------
    traj <- reactiveVal(NULL)
    traj_load_attempted <- reactiveVal(FALSE)
    traj_coord_lookup <- reactiveVal(NULL)
    traj_plot_registered <- reactiveVal(FALSE)
    applied_cell_range <- reactiveVal(c(TRAJ_DEFAULT_MIN_CELLS, traj_max_islet_cells()))

    observeEvent(input$recalculate_trajectory, {
      applied_cell_range(traj_cell_range(input$cell_range, traj_max_islet_cells()))
    })

    output$traj_range_pending <- renderUI({
      pending <- traj_cell_range(input$cell_range, traj_max_islet_cells())
      if (all(pending == applied_cell_range())) return(NULL)
      tags$small("Range not applied. Press Recalculate trajectory.")
    })

    # ---------- Lazy H5AD loader ----------
    observe({
      if (traj_load_attempted()) return()

      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) return()

      traj_load_attempted(TRUE)

      if (!requireNamespace("anndata", quietly = TRUE)) {
        traj(list(error = "R package 'anndata' not installed.", error_detail = NULL))
        return()
      }

      ad <- NULL
      norm_path <- tryCatch(normalizePath(cur, mustWork = TRUE), error = function(e) cur)

      ad <- tryCatch(
        withProgress(message = "Loading data", value = 0.2, {
          result <- islet_read_h5ad_shared(norm_path)
          incProgress(0.7)
          result
        }),
        error = function(e) {
          islet_log_event("trajectory_data_load_failed")
          NULL
        }
      )

      if (is.null(ad)) {
        traj(list(error = "Trajectory data could not be loaded.", error_detail = NULL))
        return()
      }

      # Extract obs, obsm, uns
      obs  <- tryCatch(ad$obs,  error = function(e) NULL)
      obsm <- tryCatch(ad$obsm, error = function(e) NULL)
      uns  <- tryCatch(ad$uns,  error = function(e) NULL)
      var_names <- tryCatch(rownames(ad$var), error = function(e) NULL)

      if (is.null(obs)) {
        traj(list(error = "Could not extract obs data from AnnData", error_detail = NULL))
        return()
      }

      # Ensure obs is a data.frame
      if (!inherits(obs, "data.frame")) {
        obs <- tryCatch(as.data.frame(obs, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(obs)) {
          traj(list(error = "Obs coercion failed", error_detail = NULL))
          return()
        }
      }

      # Extract only the identifiers and metadata used by the reviewer view.
      # Legacy H5AD DPT columns remain in the artifact for audit history, but
      # the live axis is recomputed from the composition feature table.
      cols <- intersect(
        c("combined_islet_id", "base_islet_id", "donor_status", "Case ID",
          "case_id", "imageid", "age", "gender", "total_cells"),
        colnames(obs)
      )
      obs2 <- as.data.frame(obs[, cols, drop = FALSE])

      # Parse combined_islet_id for Case ID and islet_key
      if (!"Case ID" %in% names(obs2)) obs2$`Case ID` <- NA_character_

      if ("imageid" %in% names(obs2)) {
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, as.character(obs2$imageid))
      }

      if ("combined_islet_id" %in% names(obs2)) {
        cid <- stringr::str_extract(obs2$combined_islet_id, "^[0-9]{3,4}")
        islet_num <- stringr::str_extract(obs2$combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
        obs2$islet_key <- ifelse(!is.na(islet_num), paste0("Islet_", islet_num), NA_character_)
      }

      # Normalize donor_status
      if ("donor_status" %in% names(obs2)) {
        obs2$donor_status <- stringr::str_trim(obs2$donor_status)
        obs2$donor_status <- dplyr::case_when(
          tolower(obs2$donor_status) %in% c("nd", "non-diabetic", "nondiabetic") ~ "ND",
          tolower(obs2$donor_status) %in% c("aab+", "aab", "autoantibody") ~ "Aab+",
          tolower(obs2$donor_status) %in% c("t1d", "type1", "diabetic") ~ "T1D",
          TRUE ~ obs2$donor_status
        )
      }

      # Populated reactively from the selected cell-count cohort below.
      obs2$pseudotime <- NA_real_

      traj(list(obs = obs2, obsm = obsm, uns = uns, var_names = var_names, adata = ad))
    })

    # ---------- Pseudotime reactive ----------
    traj_with_pseudotime <- reactive({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) return(tr)

      obs <- tr$obs
      ad  <- tr$adata

      # Only the applied range defines the fitted cohort. Slider drags alone
      # must not invalidate the diffusion map or its downstream plots.
      rng <- applied_cell_range()
      pt <- withProgress(
        # Seconds to tens of seconds depending on how many islets are in range.
        # Without a progress bar the tab simply freezes and reads as a hang.
        message = "Computing pseudotime",
        detail = sprintf("islet-composition diffusion map, %g-%g cells/islet",
                         rng[1], rng[2]),
        value = 0.4,
        traj_pseudotime_column(obs, rng))
      pseudotime_engine <- attr(pt, "engine") %||% "unavailable"
      pseudotime_n <- as.integer(attr(pt, "n_islets") %||% 0L)
      obs$pseudotime <- as.numeric(pt)

      # Stable rank axis, computed ONCE over the full selected cell-range cohort.
      #
      # This used to be done downstream in traj_data_clean() AFTER the
      # !is.na(value) and min_cells filters, which made the x-axis depend on
      # both the min_cells setting and the selected Y-axis feature: the same
      # islet moved along the axis when either changed (measured mean |delta
      # rank| 0.179, max 0.357 at min_cells = 25). Ranking here, before any
      # filtering, makes the axis a fixed property of the fitted trajectory.
      #
      # (rank - 0.5)/n rather than rank/n so the axis is symmetric on (0, 1)
      # instead of (1/n, 1].
      obs$pt_rank <- NA_real_
      if (!is.null(obs$pseudotime)) {
        finite_pt <- is.finite(obs$pseudotime)
        n_finite <- sum(finite_pt)
        if (n_finite > 0) {
          obs$pt_rank[finite_pt] <-
            (rank(obs$pseudotime[finite_pt], ties.method = "average") - 0.5) / n_finite
        }
      }

      list(obs = obs, obsm = tr$obsm, uns = tr$uns, var_names = tr$var_names,
           adata = ad, pseudotime_engine = pseudotime_engine,
           pseudotime_n = pseudotime_n, cell_range = rng)
    })

    selected_measure_label <- reactive({
      if (identical(input$traj_y_mode %||% "ratio", "ratio")) {
        region <- switch(input$traj_ratio_region %||% "core",
                         core = "core", peri = "peri", "core + 20 µm peri")
        sprintf("%s — log1p(cells / endocrine cells, %s)",
                input$traj_ratio_type %||% "(none)", region)
      } else {
        sprintf("%s — normalized mean marker expression",
                input$traj_feature %||% "(none)")
      }
    })

    trajectory_provenance <- reactive({
      tr <- tryCatch(traj_with_pseudotime(), error = function(e) NULL)
      rng <- applied_cell_range()
      engine <- if (is.null(tr) || !is.null(tr$error)) "unavailable" else
        tr$pseudotime_engine %||% "unavailable"
      n <- if (is.null(tr) || !is.null(tr$error)) 0L else tr$pseudotime_n %||% 0L
      list(
        method = sprintf("Pseudotime method: %s", engine),
        scale = sprintf("Pseudotime scale: %s",
                        if (identical(input$pt_scale %||% "rank", "raw"))
                          "raw diffusion pseudotime" else "rank-normalized within fitted cohort"),
        range = sprintf(paste0("Cells per islet range: %g-%g ",
                              "(defines the %s-islet cohort used to fit pseudotime)"),
                        rng[1], rng[2], format(n, big.mark = ",")),
        measure = paste("Selected measure:", selected_measure_label())
      )
    })

    output$traj_feature_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) {
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
      }

      # Cell-type ratio mode lists the retained phenotypes, not the 31 markers.
      # Kept in a SEPARATE input id from the marker list so switching Y mode
      # cannot leave a marker name selected against the ratio data (or vice
      # versa) -- the two namespaces do not overlap and a stale value would
      # silently plot nothing.
      if (identical(input$traj_y_mode %||% "ratio", "ratio")) {
        types <- barlow_ratio_types()
        if (!length(types)) {
          return(selectInput(ns("traj_ratio_type"), "Select cell type:",
                             choices = c("Composition features unavailable" = ""),
                             selected = ""))
        }
        sel <- if (!is.null(input$traj_ratio_type) && input$traj_ratio_type %in% types)
          input$traj_ratio_type else if ("Beta cell" %in% types) "Beta cell" else types[1]
        return(selectInput(ns("traj_ratio_type"), "Select cell type:",
                           choices = types, selected = sel))
      }

      # Get available features from AnnData
      available_features <- NULL
      if (!is.null(tr$var_names)) {
        available_features <- tr$var_names
      } else {
        # Fallback - try to get from the loaded trajectory
        tryCatch({
          adata_path <- resolve_traj_path()
          if (!is.null(adata_path) && file.exists(adata_path)) {
            adata <- islet_read_h5ad_shared(adata_path)
            available_features <- rownames(adata$var)
          }
        }, error = function(e) NULL)
      }

      if (is.null(available_features) || length(available_features) == 0) {
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
      }

      # Organize features into categories
      hormone_markers  <- intersect(c("INS", "GCG", "SST"), available_features)
      immune_markers   <- intersect(c("CD3e", "CD4", "CD8a", "CD68", "CD163", "CD20", "CD45", "HLADR"), available_features)
      vascular_markers <- intersect(c("CD31", "CD34", "SMA", "ColIV"), available_features)
      neural_markers   <- intersect(c("B3TUBB", "GAP43", "PGP9.5"), available_features)
      other_markers    <- setdiff(available_features, c(hormone_markers, immune_markers, vascular_markers, neural_markers))
      # These are MARKERS, not phenotypes, so phenotype_keep() does not apply --
      # the reduced-set rule governs the cell-type selector above, not this one.

      # Build choices list with categories.
      # Note: peri-islet immune/density columns are NOT exposed here as
      # Y-axis features. Users select pseudotime mode (Core / Core+Peri)
      # to incorporate the peri-zone and structural-density signals into
      # the X-axis (pseudotime) instead, then visualise the associated
      # markers (CD3e, CD8a, CD45, CD31, CD34, PDPN, etc.) on Y.
      choices <- list()
      if (length(hormone_markers)  > 0) choices[["Hormone Markers"]]  <- setNames(hormone_markers, hormone_markers)
      if (length(immune_markers)   > 0) choices[["Immune Markers"]]   <- setNames(immune_markers, immune_markers)
      if (length(vascular_markers) > 0) choices[["Vascular Markers"]] <- setNames(vascular_markers, vascular_markers)
      if (length(neural_markers)   > 0) choices[["Neural Markers"]]   <- setNames(neural_markers, neural_markers)
      if (length(other_markers)    > 0) choices[["Other Features"]]   <- setNames(other_markers, other_markers)

      # Default selection
      default_feature <- if ("INS" %in% available_features) "INS" else available_features[1]

      selectInput(ns("traj_feature"), "Select Feature:",
                  choices = choices,
                  selected = default_feature)
    })

    # ---------- Validate selected feature ----------
    traj_validate <- reactive({
      y_mode <- input$traj_y_mode %||% "ratio"
      if (!islet_valid_scalar_choice(y_mode, c("marker", "ratio"),
                                     "trajectory_mode_invalid")) return(NULL)
      # Ratio mode selects a cell type, which is deliberately NOT in var_names --
      # validating it against the marker list would reject every valid choice.
      if (identical(y_mode, "ratio")) {
        ct <- input$traj_ratio_type
        if (is.null(ct)) return(NULL)
        if (!islet_valid_scalar_choice(ct, barlow_ratio_types(),
                                       "trajectory_ratio_invalid")) return(NULL)
        region <- input$traj_ratio_region %||% "core"
        if (!islet_valid_scalar_choice(region, c("both", "core", "peri"),
                                       "trajectory_region_invalid")) return(NULL)
        return(list(feature = ct, status = "ready", y_mode = "ratio",
                    ratio_region = region))
      }
      selected_feature <- input$traj_feature
      if (is.null(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      adata <- tryCatch(islet_read_h5ad_shared(adata_path), error = function(e) NULL)
      if (is.null(adata)) return(NULL)

      var_names <- rownames(adata$var)
      if (!islet_valid_scalar_choice(selected_feature, var_names,
                                     "trajectory_feature_invalid")) return(NULL)

      list(feature = selected_feature, status = "ready")
    })

    # ---------- Main data pipeline ----------
    traj_data_clean <- reactive({
      ms <- traj_validate()
      if (is.null(ms)) return(NULL)

      selected_feature <- ms$feature
      if (is.null(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error)) return(NULL)

      # Use cached adata object instead of re-reading
      adata <- tr$adata
      if (is.null(adata)) {
        adata <- tryCatch(islet_read_h5ad_shared(adata_path), error = function(e) NULL)
        if (is.null(adata)) return(NULL)
      }

      var_names <- rownames(adata$var)
      y_mode <- ms$y_mode %||% "marker"
      if (identical(y_mode, "marker") && !selected_feature %in% var_names) return(NULL)

      obs_df <- tr$obs
      if (is.null(obs_df)) return(NULL)

      expression_vals <- if (identical(y_mode, "ratio")) {
        # Barlow Fig 2F: log(# cells / # endocrine cells). Joined by islet_id,
        # which is the same key as `combined_islet_id` here.
        rv <- barlow_ratio_values(selected_feature, ms$ratio_region)
        if (is.null(rv)) rep(NA_real_, nrow(obs_df))
        else rv$value[match(as.character(obs_df$combined_islet_id), rv$islet_id)]
      } else {
        # Extract expression values from .X (markers only).
        tryCatch({
          idx <- which(var_names == selected_feature)
          as.numeric(adata$X[, idx])
        }, error = function(e) {
          islet_log_event("trajectory_feature_extract_failed")
          rep(NA_real_, nrow(obs_df))
        })
      }

      # Get donor_status with proper normalization
      donor_clean <- obs_df$donor_status
      if (is.null(donor_clean)) donor_clean <- rep(NA_character_, nrow(obs_df))

      # Extract UMAP coordinates
      umap_coords <- NULL
      if (!is.null(tr$obsm) && "X_umap" %in% names(tr$obsm)) {
        umap_coords <- tr$obsm$X_umap
      }

      # Get donor ID for coloring
      donor_ids <- as.character(obs_df$imageid)

      # Build result data.frame
      islet_keys <- as.character(obs_df$base_islet_id)
      islet_keys <- gsub("^Islet_Islet_", "Islet_", islet_keys)

      result_df <- data.frame(
        case_id = as.character(obs_df$imageid),
        islet_key = islet_keys,
        combined_islet_id = as.character(obs_df$combined_islet_id),
        pt = as.numeric(obs_df$pseudotime),
        pt_rank = if (!is.null(obs_df$pt_rank)) as.numeric(obs_df$pt_rank) else NA_real_,
        donor_status = donor_clean,
        donor_id = donor_ids,
        value = as.numeric(expression_vals),
        feature_name = selected_feature,
        total_cells = as.integer(obs_df$total_cells),
        stringsAsFactors = FALSE
      )

      # Guard: fill missing total_cells with 1 (safe default weight)
      if (any(is.na(result_df$total_cells))) {
        result_df$total_cells[is.na(result_df$total_cells)] <- 1L
      }

      # Add UMAP coordinates
      if (!is.null(umap_coords) && nrow(umap_coords) == nrow(result_df)) {
        result_df$umap_1 <- as.numeric(umap_coords[, 1])
        result_df$umap_2 <- as.numeric(umap_coords[, 2])
      } else {
        result_df$umap_1 <- NA_real_
        result_df$umap_2 <- NA_real_
      }

      # Merge diameter and Leiden cluster cols from prepared() comp table.
      prep_data <- prepared()
      if (!is.null(prep_data) && !is.null(prep_data$comp)) {
        leiden_cols <- intersect(c("leiden_0.3", "leiden_0.5", "leiden_0.8", "leiden_1.0"),
                                 colnames(prep_data$comp))
        lookup_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                         leiden_cols)
        size_lookup <- prep_data$comp[intersect(lookup_cols, colnames(prep_data$comp))]
        # comp$`Case ID` is integer (data_loading.R:346); match the type so merge succeeds.
        result_df$`Case ID` <- suppressWarnings(as.integer(result_df$case_id))
        result_df$`Donor Status` <- result_df$donor_status

        # Base R merge (faster for this use case)
        merged <- merge(result_df, size_lookup,
                        by = c("Case ID", "Donor Status", "islet_key"),
                        all.x = TRUE, sort = FALSE)

        match_idx <- match(
          paste(result_df$`Case ID`, result_df$`Donor Status`, result_df$islet_key),
          paste(merged$`Case ID`, merged$`Donor Status`, merged$islet_key)
        )

        if ("islet_diam_um" %in% colnames(merged)) {
          result_df$islet_diam_um <- merged$islet_diam_um[match_idx]
        }
        for (lc in leiden_cols) {
          if (lc %in% colnames(merged)) {
            result_df[[lc]] <- as.character(merged[[lc]][match_idx])
          }
        }
      }

      if (!"islet_diam_um" %in% colnames(result_df)) {
        result_df$islet_diam_um <- NA_real_
      }

      # Add spatial coordinates from segmentation data
      if (!is.null(islet_spatial_lookup)) {
        result_df$centroid_x_um <- NA_real_
        result_df$centroid_y_um <- NA_real_
        result_df$has_spatial <- FALSE

        for (i in seq_len(nrow(result_df))) {
          annot <- get_islet_annotations(result_df$case_id[i], result_df$islet_key[i])
          if (!is.null(annot)) {
            result_df$centroid_x_um[i] <- annot$centroid_x
            result_df$centroid_y_um[i] <- annot$centroid_y
            result_df$has_spatial[i] <- TRUE
          }
        }

      }

      # Filter valid rows
      valid_idx <- which(
        is.finite(result_df$pt) &
        !is.na(result_df$value) &
        !is.na(result_df$case_id) &
        !is.na(result_df$islet_key)
      )

      if (length(valid_idx) == 0) {
        return(NULL)
      }

      result_df <- result_df[valid_idx, , drop = FALSE]

      # Enforce the same cell-count cohort used to fit the trajectory. This is
      # intentionally redundant with the NA pseudotime mask and prevents a key
      # mismatch from leaking an out-of-range islet into the plot.
      cr <- applied_cell_range()
      if ("total_cells" %in% colnames(result_df)) {
        tc <- suppressWarnings(as.integer(result_df$total_cells))
        result_df <- result_df[is.finite(tc) & tc >= cr[1] & tc <= cr[2], , drop = FALSE]
        if (nrow(result_df) == 0) return(NULL)
      }

      # Pseudotime axis scale.
      #
      # "rank" (default) spreads points uniformly, which counteracts the raw-DPT
      # compression (95% of islets fall between 0.25-0.70 due to diffusion
      # manifold geometry). "raw" shows true DPT spacing, so where islets
      # actually concentrate is visible -- note that on the rank axis, LOESS
      # curvature and bin density are artifacts of the transform, since ranks
      # are uniform by construction.
      #
      # The rank itself is computed in traj_with_pseudotime() over the complete
      # selected cell-range cohort, so it does not shift with the selected
      # feature. Do not re-rank here.
      result_df <- apply_pt_scale(result_df, input$pt_scale %||% "rank")

      return(result_df)
    })

    # ---------- Status bar ----------
    output$traj_status <- renderUI({
      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) {
        return(tags$div(style = "color:#b00;",
          "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
      }
      if (!requireNamespace("anndata", quietly = TRUE)) {
        return(tags$div(style = "color:#b00;",
          "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
      }
      tr <- traj_with_pseudotime()
      if (is.null(tr)) {
        return(tags$div(style = "color:#b00;",
          "AnnData object not available or H5AD failed to load."))
      }
      if (!is.null(tr$error)) {
        install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
        return(tags$div(style = "color:#b00;",
          sprintf("Trajectory load error: %s%s", tr$error, install_hint)))
      }
      nfit <- tr$pseudotime_n %||% 0L
      src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur

      tags$div(style = "color: #555; font-size: 12px; margin-bottom: 10px;",
        title = paste("Source:", basename(src)),
        sprintf("Islets: %s total | %s fitted",
          format(nrow(tr$obs), big.mark = ","), format(nfit, big.mark = ","))
      )
    })

    # ---------- Trajectory scatter plot ----------
    output$traj_scatter <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) {
        tr <- traj()
        if (is.null(tr) || (is.null(tr$error) &&
            identical(input$traj_y_mode %||% "ratio", "ratio") &&
            is.null(input$traj_ratio_type))) {
          return(plotly_islet_empty("Loading data"))
        }
        return(plotly_islet_empty("No trajectory data available"))
      }


      # Determine color mapping and aesthetics
      color_by   <- input$traj_color_by %||% "donor_status"
      size_by    <- input$traj_point_size %||% "uniform"
      trend_type <- input$traj_show_trend %||% "by_donor"
      alpha_val  <- 1 - (input$traj_alpha %||% 0.4)
      point_size <- input$traj_point_size_slider %||% 3.0
      # 95% CI bands on the LOESS trends (Bodenmiller 2025 Fig 2D/5E convention:
      # every smoothed trajectory carries a confidence band).
      show_ci    <- isTRUE(input$traj_trend_ci %||% TRUE)
      df$loess_weight <- log1p(df$total_cells)
      trend_df <- df

      # Pseudotime extent, taken from the data rather than assumed to be [0, 1] — on the
      # "raw" scale the axis is true DPT and does not span [0, 1]. Shared with the strip
      # heatmap below via traj_pt_extent() so the two panels stay aligned.
      pt_ext <- traj_pt_extent(df$pt)
      pt_dom <- pt_ext$dom

      # Always apply jitter — use the kept-rows range as the jitter scale
      # so ghost points land in the same visual frame.
      #
      # Jitter is NOT clamped back into the domain: the previous pmax/pmin
      # clamp piled roughly the outermost 1.2% of points onto the exact
      # endpoints (~60 islets per side), producing artificial spikes at both
      # ends. The axis is widened by the jitter amount instead (see pt_lim).
      #
      # set.seed() inside a render function would reset the session-global RNG
      # stream, so the previous seed is saved and restored around the draw.
      .seed_saved <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else NULL
      set.seed(42)
      jitter_amount_x <- pt_ext$jitter
      jitter_amount_y <- diff(range(df$value, na.rm = TRUE)) * 0.08
      df$pt    <- df$pt + runif(nrow(df), -jitter_amount_x, jitter_amount_x)
      df$value <- df$value + runif(nrow(df), -jitter_amount_y, jitter_amount_y)
      if (is.null(.seed_saved)) {
        suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
      } else {
        assign(".Random.seed", .seed_saved, envir = .GlobalEnv)
      }

      # Plot limits: domain plus the jitter margin, so jittered points are not
      # clipped and no pile-up appears at the endpoints. The strip heatmap uses the
      # same padded limits (traj_pt_extent()$lim) so the panels stay aligned.
      pt_lim <- pt_ext$lim

      # Create base aesthetic mapping with cell count for hover
      df$hover_cells <- paste0("Pseudotime: ", signif(trend_df$pt, 3),
              "<br>Value: ", signif(trend_df$value, 3),
              "<br>Cells: ", df$total_cells)
      df$click_key <- paste(df$case_id, df$islet_key, sep = "|")
      aes_mapping <- aes(x = pt, y = value, text = hover_cells, key = click_key)

      # Add color aesthetic
      if (color_by == "donor_status") {
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        aes_mapping$colour <- as.name("donor_status")
        color_title <- "Donor Status"
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        aes_mapping$colour <- as.name("donor_id")
        color_title <- "Donor ID"
      } else if (startsWith(color_by, "leiden_") && color_by %in% colnames(df)) {
        # Sort factor levels numerically so legend reads 0,1,2,... instead of 0,1,10,11
        lvl_vals <- unique(stats::na.omit(df[[color_by]]))
        num_vals <- suppressWarnings(as.numeric(lvl_vals))
        sorted_lvls <- if (!any(is.na(num_vals))) lvl_vals[order(num_vals)] else sort(lvl_vals)
        df[[color_by]] <- factor(df[[color_by]], levels = sorted_lvls)
        aes_mapping$colour <- as.name(color_by)
        color_title <- paste("Leiden", sub("^leiden_", "res ", color_by))
      } else if (startsWith(color_by, "leiden_")) {
        # Column not present in df (e.g., Excel fallback) — fall back to donor_status
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        aes_mapping$colour <- as.name("donor_status")
        color_title <- "Donor Status"
        color_by <- "donor_status"
      }

      # Add size aesthetic if not uniform
      if (size_by != "uniform") {
        if (size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
          if (!all(is.na(df$islet_diam_um))) {
            aes_mapping$size <- as.name("islet_diam_um")
          }
        } else if (size_by == "total_cells" && "total_cells" %in% colnames(df)) {
          df$size_cells <- sqrt(df$total_cells)
          aes_mapping$size <- as.name("size_cells")
        }
      }

      # Create plot
      g <- ggplot(df, aes_mapping)
      if (size_by == "uniform") {
        g <- g + geom_point(alpha = alpha_val, size = point_size, stroke = 0)
      } else {
        g <- g + geom_point(alpha = alpha_val, stroke = 0)
      }
      g <- g +
        theme_islet() +
        theme(
          legend.position = "none",
          plot.margin = unit(c(5, 5, 5, 5), "pt"),
          plot.title = element_text(margin = margin(b = 15))
        )

      # Apply color scales
      if (color_by == "donor_status") {
        g <- g + scale_color_manual(
          values = donor_colors_reactive(),
          name = color_title, drop = FALSE)
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        n_donors <- length(levels(df$donor_id))
        if (n_donors <= 12) {
          colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          names(colors) <- levels(df$donor_id)[1:length(colors)]
        } else {
          colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          names(colors) <- levels(df$donor_id)
        }
        g <- g + scale_color_manual(values = colors, name = color_title, drop = FALSE)
      } else if (startsWith(color_by, "leiden_") && color_by %in% colnames(df)) {
        clust_levels <- levels(df[[color_by]])
        n_clust <- length(clust_levels)
        clust_colors <- if (n_clust <= 12) {
          RColorBrewer::brewer.pal(max(3, min(12, n_clust)), "Set3")[seq_len(n_clust)]
        } else {
          scales::hue_pal()(n_clust)
        }
        names(clust_colors) <- clust_levels
        g <- g + scale_color_manual(values = clust_colors, name = color_title, drop = FALSE)
      }

      # Size scale
      if (size_by != "uniform" && size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
        if (!all(is.na(df$islet_diam_um))) {
          size_range <- c(0.5 * point_size, 2.5 * point_size)
          g <- g + scale_size_continuous(name = "Islet Diameter (um)", range = size_range, guide = "none")
        }
      } else if (size_by == "total_cells" && "size_cells" %in% colnames(df)) {
        size_range <- c(0.3 * point_size, 3.0 * point_size)
        g <- g + scale_size_continuous(name = "Cell Count", range = size_range, guide = "none")
      }

      # Trend lines
      pt_range <- range(df$pt, na.rm = TRUE)

      if (trend_type == "overall") {
        g <- g + geom_smooth(data = trend_df, method = "loess", se = show_ci, color = "black",
                             linewidth = 1, fullrange = FALSE, span = 0.75,
                             inherit.aes = FALSE,
                             aes(x = pt, y = value, weight = loess_weight))
      } else if (trend_type == "by_donor") {
        if ("donor_status" %in% names(df) && !all(is.na(df$donor_status))) {
          df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
          trend_colors <- donor_colors_reactive()

          if (color_by == "donor_id" || startsWith(color_by, "leiden_")) {
            # Use explicit color per-status loop — adding a second scale_color_manual
            # after a Leiden/donor_id cluster scale would override it and blank out the points.
            for (status in c("ND", "Aab+", "T1D")) {
              df_subset <- trend_df[trend_df$donor_status == status & !is.na(trend_df$donor_status), ]
              if (nrow(df_subset) > 0) {
                g <- g + geom_smooth(data = df_subset,
                                     aes(x = pt, y = value, weight = loess_weight),
                                     inherit.aes = FALSE,
                                     method = "loess", se = show_ci,
                                     linewidth = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
                                     color = trend_colors[status],
                                     fill = trend_colors[status], alpha = 0.15)
              }
            }
          } else {
            g <- g + geom_smooth(data = trend_df,
              method = "loess", se = show_ci,
              linewidth = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
              inherit.aes = FALSE,
              mapping = aes(x = pt, y = value, group = donor_status,
                            color = donor_status, weight = loess_weight))
          }
        }
      }

      # Labels
      # The y-axis label MUST name the quantity, not just the feature: a ratio
      # and a denoised marker mean are different units on the same plot slot.
      if (identical(input$traj_y_mode %||% "ratio", "ratio")) {
        selected_feature <- input$traj_ratio_type %||% "Selected cell type"
        reg <- input$traj_ratio_region %||% "core"
        metric_label <- sprintf("\u2014 log(cells per endocrine cell, %s)",
                                switch(reg, core = "core", peri = "peri",
                                       "core + 20 \u00b5m peri"))
      } else {
        selected_feature <- input$traj_feature %||% "Selected feature"
        metric_label <- "\u2014 normalized expression (mean per islet, 0\u20131)"
      }
      x_label <- if ((input$pt_scale %||% "rank") == "raw") {
        "Pseudotime (raw DPT)"
      } else {
        "Pseudotime (rank-normalized)"
      }

      g <- g + labs(
        x = x_label,
        y = paste(selected_feature, metric_label),
        title = NULL
      ) +
      coord_cartesian(xlim = pt_lim, ylim = range(df$value, na.rm = TRUE) * c(1, 1) + c(-0.05, 0.05) * diff(range(df$value, na.rm = TRUE)))

      # Convert to plotly with custom data for reliable click handling
      p <- ggplotly(g, tooltip = c("colour", "text"), source = ns("traj_scatter"))
      # Strip "colour: black" leak (from overall LOESS line's color="black"
      # parameter) and collapse duplicate donor_status entries in hover text.
      p <- strip_plotly_color_leak(p)

      # Store coordinate lookup as fallback for click handling
      if ("combined_islet_id" %in% colnames(df) && "case_id" %in% colnames(df)) {
        traj_coord_lookup(data.frame(
          pt                = df$pt,
          value             = df$value,
          case_id           = df$case_id,
          islet_key         = df$islet_key,
          combined_islet_id = df$combined_islet_id,
          donor_status      = df$donor_status,
          stringsAsFactors  = FALSE
        ))
      }

      # Layout and register click events
      fonts <- plotly_islet_fonts()
      p <- p %>%
        layout(
          showlegend = FALSE,
          title = NULL,
          font = fonts$global,
          margin = list(l = 80, r = 20, t = 30, b = 70, pad = 5),
          xaxis = list(fixedrange = FALSE, automargin = TRUE, range = pt_lim,
                       title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(fixedrange = FALSE, automargin = TRUE,
                       title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(font = fonts$legend_text,
                        title = list(font = fonts$legend_title))
        ) %>%
        plotly_islet_config("trajectory_scatter") %>%
        event_register("plotly_click") %>%
        islet_click_bridge(ns("traj_native_click"))
      # Invalidating this flag schedules the click observer for a later reactive
      # pass, after renderPlotly() has recorded the registered source in the
      # session. Keeping event_data() dormant until then avoids plotly's deferred
      # false-positive registration warning.
      traj_plot_registered(TRUE)
      p
    })

    # ---------- Click handler for trajectory scatter ----------
    # Do not ask plotly for this event while the Trajectory tab (and therefore
    # its suspended plot output) is still hidden. event_data() schedules its
    # registration check for the end of the flush, outside suppressWarnings(),
    # so probing it at startup produces a false warning even though the rendered
    # plot correctly calls event_register().
    traj_click_event <- reactive({
      req(identical(active_tab(), "Trajectory"))
      req(isTRUE(traj_plot_registered()))
      event_data("plotly_click", source = ns("traj_scatter"))
    })
    observeEvent(traj_click_event(), {
      click_data <- traj_click_event()
      if (is.null(click_data)) return()

      # --- Priority 1: Parse key (format "case_id|islet_key") ---
      case_id   <- NULL
      islet_key <- NULL
      click_key <- click_data$key
      if (!is.null(click_key) && length(click_key) > 0 && !is.na(click_key[1]) && nzchar(click_key[1])) {
        key_parts <- strsplit(as.character(click_key[1]), "\\|")[[1]]
        if (length(key_parts) == 2L &&
            islet_safe_path_segment(key_parts[1]) &&
            islet_safe_path_segment(key_parts[2])) {
          case_id   <- key_parts[1]
          islet_key <- key_parts[2]
        } else {
          islet_log_validation_failure("islet_key_invalid")
        }
      }

      # --- Priority 2: Fallback to customdata ---
      if (is.null(case_id) || is.null(islet_key)) {
        custom <- click_data$customdata
        if (!is.null(custom) && length(custom) >= 2) {
          case_id   <- custom[[1]]
          islet_key <- custom[[2]]
        }
      }

      # --- Priority 3: Fallback to coordinate lookup ---
      if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
        lookup <- traj_coord_lookup()
        if (!is.null(lookup) && nrow(lookup) > 0) {
          # Tolerances are scaled to each axis's actual span rather than a flat
          # 0.01, which only made sense while pt was always rank-normalized to
          # [0, 1]. On the "raw" DPT scale the span is much narrower, and a
          # fixed 0.01 would match nearly the whole axis.
          pt_span  <- diff(range(lookup$pt, na.rm = TRUE))
          val_span <- diff(range(lookup$value, na.rm = TRUE))
          if (!is.finite(pt_span)  || pt_span  <= 0) pt_span  <- 1
          if (!is.finite(val_span) || val_span <= 0) val_span <- 1
          matches <- which(
            abs(lookup$pt - click_data$x) < 0.01 * pt_span &
            abs(lookup$value - click_data$y) < 0.01 * val_span
          )
          if (length(matches) > 0) {
            idx <- matches[1]
            case_id   <- lookup$case_id[idx]
            islet_key <- lookup$islet_key[idx]
          }
        }
      }

      # Centroid lookup + zero-pad fallback + donor join + selected_islet() —
      # shared with the Plot handler and the native-click bridge.
      select_islet_from_key(case_id, islet_key, prepared, selected_islet)
    }, ignoreInit = TRUE)

    # macOS-proof native-click bridge (see islet_click_bridge in 00_globals.R).
    # Fires on a native click even when plotly_click is suppressed by plotly's
    # click/drag discrimination on macOS.
    observeEvent(input$traj_native_click, {
      nc <- input$traj_native_click
      if (is.null(nc) || is.null(nc$key) || !nzchar(as.character(nc$key)[1])) return()
      parts <- strsplit(as.character(nc$key)[1], "\\|")[[1]]
      if (length(parts) != 2L ||
          !islet_safe_path_segment(parts[1]) ||
          !islet_safe_path_segment(parts[2])) {
        islet_log_validation_failure("islet_key_invalid")
        return()
      }
      select_islet_from_key(parts[1], parts[2], prepared, selected_islet)
    })

    # ---------- Segmentation viewer panel ----------
    output$segmentation_viewer_panel <- renderUI({
      # Only render when Trajectory tab is active to avoid duplicate non-namespaced output IDs
      if (active_tab() != "Trajectory") return(NULL)
      info <- selected_islet()
      if (is.null(info)) {
        return(div(style = "min-height: 275px; display: flex; align-items: center; justify-content: center; color: #666;",
                   "No islet selected"))
      }

      # Check if single-cell drill-down data exists for this islet
      has_cells <- drilldown_available() &&
        !is.null(load_islet_cells(info$case_id, info$islet_key))

      # Marker choices for color-by dropdown
      marker_choices <- c("phenotype" = "phenotype")
      if (has_cells) {
        cells <- load_islet_cells(info$case_id, info$islet_key)
        marker_cols <- setdiff(colnames(cells),
                               c("X_centroid", "Y_centroid", "phenotype", "cell_region",
                                 "Cell.Area", "Cell Area", "Nucleus.Area", "Nucleus Area"))
        if (length(marker_cols) > 0) {
          marker_choices <- c(marker_choices, setNames(marker_cols, marker_cols))
        }
      }

      # Build title: "Islet_N (case_id, status, age, sex)"
      title_detail <- info$case_id
      if (!is.null(info$donor_status)) title_detail <- paste0(title_detail, ", ", info$donor_status)
      if (!is.null(info$donor_age) && !is.na(info$donor_age)) title_detail <- paste0(title_detail, ", age ", round(info$donor_age))
      if (!is.null(info$donor_gender) && !is.na(info$donor_gender) && nzchar(info$donor_gender)) title_detail <- paste0(title_detail, ", ", info$donor_gender)

      div(
        # Header: title + close button
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px;",
          h4(paste0(info$islet_key, " (", title_detail, ")"),
             style = "margin: 0; color: #0066CC; font-size: 20px;"),
          actionButton(ns("clear_segmentation"), "Close", class = "btn btn-sm btn-outline-secondary")
        ),
        # Map uses the full card width; controls and composition sit beneath it.
        fluidRow(
          column(12,
            plotOutput("islet_segmentation_view", height = "430px"),
            div(style = "text-align: right; margin-top: 6px;",
              downloadButton("dl_islet_segmentation_view", "Download PNG",
                             class = "btn-sm btn-outline-primary"))
          ),
          column(5,
            div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
              if (has_cells) {
                tagList(
                  h5("Layers", style = "margin-top: 0; font-size: 15px;"),
                  checkboxInput("drilldown_show_cells", "Single Cells", value = TRUE),
                  checkboxInput("drilldown_show_peri_boundary", "Peri Boundary", value = TRUE),
                  checkboxInput("drilldown_show_structures", "Structures", value = TRUE),
                  hr(style = "margin: 8px 0;"),
                  conditionalPanel(
                    condition = "input.drilldown_show_cells",
                    selectInput("drilldown_color_by", "Color by",
                                choices = marker_choices, selected = "phenotype"),
                    selectInput("drilldown_palette", "Phenotype Palette",
                                choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                                selected = "High Contrast"),
                    checkboxInput("drilldown_show_peri", "Show peri-islet cells", value = TRUE)
                  )
                )
              } else {
                tagList(
                  h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-size: 15px;"),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #0066CC; margin-right: 8px;"),
                    span("Islet", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #00CCCC; margin-right: 8px;"),
                    span("+20\u00b5m", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #CC00CC; margin-right: 8px;"),
                    span("Nerve", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #CC0000; margin-right: 8px;"),
                    span("Capillary", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #00AA00; margin-right: 8px;"),
                    span("Lymphatic", style = "font-size: 14px;")
                  )
                )
              },
              hr(style = "margin: 8px 0;"),
              span(paste("Centroid:", round(info$centroid_x, 0), "/", round(info$centroid_y, 0), "\u00b5m"),
                   style = "color: #666; font-size: 13px;")
            )
          ),
          column(7,
            if (has_cells) {
              conditionalPanel(
                condition = "input.drilldown_show_cells",
                div(style = "padding: 12px; background-color: #f8f9fa; border-radius: 5px;",
                  h5("Cell Composition", style = "margin-top: 0; font-size: 16px;"),
                  plotOutput("islet_drilldown_summary", height = "300px"),
                  div(style = "text-align: right; margin-top: 6px;",
                    downloadButton("dl_islet_drilldown_summary", "Download PNG",
                                   class = "btn-sm btn-outline-primary")),
                  tableOutput("islet_drilldown_table"),
                  uiOutput("islet_drilldown_rules")
                )
              )
            }
          )
        )
      )
    })

    # ---------- Clear segmentation viewer ----------
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # Note: islet_segmentation_view renderPlot is at root level in app.R
    # (shared between Plot modal and Trajectory embedded panel)

    # ---------- Outlier badge (above scatter) ----------

    # ---------- Outlier table (below scatter) ----------

    # ---------- Legend for trajectory scatter ----------
    output$traj_legend <- renderUI({
      color_by <- input$traj_color_by %||% "donor_status"

      if (color_by == "donor_status") {
        dc <- donor_colors_reactive()
        tagList(
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["ND"]])),
            span("ND", style = "font-size: 15px;")
          ),
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["Aab+"]])),
            span("Aab+", style = "font-size: 15px;")
          ),
          div(style = "display: flex; align-items: center;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["T1D"]])),
            span("T1D", style = "font-size: 15px;")
          )
        )
      } else if (color_by == "donor_id") {
        df <- traj_data_clean()
        if (!is.null(df) && nrow(df) > 0) {
          donor_ids <- sort(unique(df$donor_id))
          n_donors <- length(donor_ids)

          # Generate colors matching the plot
          if (n_donors <= 12) {
            colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          } else {
            colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          }

          # Split into two columns
          mid_point <- ceiling(n_donors / 2)
          col1_ids <- donor_ids[1:mid_point]
          col2_ids <- if (n_donors > mid_point) donor_ids[(mid_point + 1):n_donors] else NULL

          col1_items <- lapply(seq_along(col1_ids), function(i) {
            div(style = "display: flex; align-items: center; margin-bottom: 6px;",
              div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[i])),
              span(col1_ids[i], style = "font-size: 14px;")
            )
          })

          col2_items <- if (!is.null(col2_ids)) {
            lapply(seq_along(col2_ids), function(i) {
              idx <- mid_point + i
              div(style = "display: flex; align-items: center; margin-bottom: 6px;",
                div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[idx])),
                span(col2_ids[i], style = "font-size: 14px;")
              )
            })
          } else {
            NULL
          }

          tagList(
            div(style = "display: flex; gap: 10px;",
              div(style = "flex: 1;", col1_items),
              if (!is.null(col2_items)) div(style = "flex: 1;", col2_items)
            )
          )
        } else {
          p("No data available", style = "font-size: 12px; color: #999;")
        }
      } else if (startsWith(color_by, "leiden_")) {
        df <- traj_data_clean()
        if (is.null(df) || !(color_by %in% colnames(df))) {
          return(p("Leiden clusters unavailable (rebuild islet_explorer.h5ad).",
                   style = "font-size: 12px; color: #999;"))
        }
        lvl_vals <- unique(stats::na.omit(df[[color_by]]))
        num_vals <- suppressWarnings(as.numeric(lvl_vals))
        clust_levels <- if (!any(is.na(num_vals))) lvl_vals[order(num_vals)] else sort(lvl_vals)
        n_clust <- length(clust_levels)
        clust_colors <- if (n_clust <= 12) {
          RColorBrewer::brewer.pal(max(3, min(12, n_clust)), "Set3")[seq_len(n_clust)]
        } else {
          scales::hue_pal()(n_clust)
        }
        mid_point <- ceiling(n_clust / 2)
        col1_idx <- seq_len(mid_point)
        col2_idx <- if (n_clust > mid_point) (mid_point + 1):n_clust else integer(0)
        make_row <- function(i) {
          div(style = "display: flex; align-items: center; margin-bottom: 6px;",
            div(style = sprintf(
              "width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;",
              clust_colors[i])),
            span(paste("Cluster", clust_levels[i]), style = "font-size: 14px;")
          )
        }
        tagList(
          div(style = "display: flex; gap: 10px;",
            div(style = "flex: 1;", lapply(col1_idx, make_row)),
            if (length(col2_idx) > 0) div(style = "flex: 1;", lapply(col2_idx, make_row))
          )
        )
      }
    })

    # ---------- UMAP: Donor Status ----------
    traj_umap_donor_plot <- reactive({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = ISLET_ANNOTATE_SIZE) +
          theme_void())
      }

      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))

      ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_color_manual(values = donor_colors_reactive()) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP: Donor Status") +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())
    })
    output$traj_umap_donor <- renderPlot({ traj_umap_donor_plot() })
    output$dl_umap_donor <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_donor"),
      content = function(file) {
        p <- traj_umap_donor_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- UMAP: Selected Feature ----------
    traj_umap_feature_plot <- reactive({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = ISLET_ANNOTATE_SIZE) +
          theme_void())
      }

      selected_feature <- selected_measure_label()
      color_label <- if (identical(input$traj_y_mode %||% "ratio", "ratio"))
        "Log cell ratio" else "Expression"

      g <- ggplot(df, aes(x = umap_1, y = umap_2, color = value)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2",
             title = paste("UMAP:", selected_feature),
             color = color_label) +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())

      # Continuous colormap scaled to data min/max
      val_range <- range(df$value, na.rm = TRUE)
      g <- g + scale_color_viridis_c(option = "inferno", limits = val_range,
                                      na.value = "#bbbbbb")

      g
    })
    output$traj_umap_feature <- renderPlot({ traj_umap_feature_plot() })
    output$dl_umap_feature <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_feature"),
      content = function(file) {
        p <- traj_umap_feature_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- Heatmap: Dominant Status × Feature Intensity ----------
    # Hue per pseudotime bin = dominant donor_status, computed as the
    # argmax of per-status NORMALIZED density (count_in_bin / total_for_status).
    # Normalization corrects for global cohort imbalance
    # (ND ~2300, Aab+ ~1668, T1D ~1246 islets) so dominance reflects
    # enrichment relative to each status's pool, not raw count.
    # Intensity per bin = bin-mean of selected feature, normalized via
    # 5–95th percentile clip, blended in [INTENSITY_FLOOR, 1] toward the
    # donor color so even the dimmest bin is a clear tint (never white).
    # Rendered with ggplot + geom_tile to guarantee bars span [0, 1].
    INTENSITY_FLOOR <- 0.50  # minimum blend factor — bars never go pale-white
    output$traj_heatmap <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      feature_name <- df$feature_name[1]
      donor_cols   <- donor_colors_reactive()
      n_bins       <- 25
      # Bin over the data's actual pseudotime domain rather than an assumed [0, 1] — on
      # the "raw" DPT scale the axis does not span [0, 1]. Bins cover `dom` (the data);
      # the AXIS uses `lim` (dom + the scatter's jitter margin) so this strip lines up
      # with the scatter above it. Same helper, so the two cannot drift apart.
      hm_ext       <- traj_pt_extent(df$pt)
      hm_dom       <- hm_ext$dom
      hm_lim       <- hm_ext$lim
      brks         <- seq(hm_dom[1], hm_dom[2], length.out = n_bins + 1)
      bin_w        <- diff(brks)[1]
      mids         <- head(brks, -1) + bin_w / 2

      df$pt_bin <- cut(df$pt, breaks = brks,
                       include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(df$pt_bin)

      status_levels <- c("ND", "Aab+", "T1D")
      status_totals <- sapply(status_levels, function(s) sum(df$donor_status == s))

      bin_rows <- lapply(seq_along(bin_levels), function(i) {
        bd <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], , drop = FALSE]
        if (nrow(bd) == 0) {
          return(data.frame(idx = i, x = mids[i], n = 0L,
                            n_ND = 0L, n_Aab = 0L, n_T1D = 0L,
                            dens_ND = 0, dens_Aab = 0, dens_T1D = 0,
                            dominant = NA_character_, feat_mean = NA_real_,
                            stringsAsFactors = FALSE))
        }
        counts <- c(ND = sum(bd$donor_status == "ND"),
                    `Aab+` = sum(bd$donor_status == "Aab+"),
                    T1D = sum(bd$donor_status == "T1D"))
        dens <- ifelse(status_totals > 0, counts / status_totals, 0)
        names(dens) <- status_levels
        dominant <- names(dens)[which.max(dens)]
        data.frame(idx = i, x = mids[i], n = nrow(bd),
                   n_ND = counts[["ND"]], n_Aab = counts[["Aab+"]], n_T1D = counts[["T1D"]],
                   dens_ND = dens[["ND"]], dens_Aab = dens[["Aab+"]], dens_T1D = dens[["T1D"]],
                   dominant = dominant,
                   feat_mean = mean(bd$value, na.rm = TRUE),
                   stringsAsFactors = FALSE)
      })
      hm <- do.call(rbind, bin_rows)
      if (is.null(hm) || nrow(hm) == 0) return(NULL)

      # Intensity from feature value, robust 5–95th percentile clip
      q   <- stats::quantile(hm$feat_mean, c(0.05, 0.95), na.rm = TRUE)
      rng <- diff(q)
      raw_intensity <- if (is.finite(rng) && rng > 0) {
        pmax(0, pmin(1, (hm$feat_mean - q[1]) / rng))
      } else {
        rep(1, nrow(hm))
      }
      raw_intensity[!is.finite(raw_intensity)] <- INTENSITY_FLOOR
      # Rescale so the lowest bin is INTENSITY_FLOOR-saturated (clearly tinted),
      # the highest bin is fully saturated.
      hm$intensity <- INTENSITY_FLOOR + (1 - INTENSITY_FLOOR) * raw_intensity

      blend_to_white <- function(hex, t) {
        rgb_full <- grDevices::col2rgb(hex) / 255
        grDevices::rgb(red   = t * rgb_full[1] + (1 - t) * 1,
                       green = t * rgb_full[2] + (1 - t) * 1,
                       blue  = t * rgb_full[3] + (1 - t) * 1)
      }
      hm$fill <- mapply(function(dom, ity) {
        if (is.na(dom) || is.null(donor_cols[[dom]])) return("#BBBBBB")
        blend_to_white(donor_cols[[dom]], ity)
      }, hm$dominant, hm$intensity, USE.NAMES = FALSE)

      hm$hover <- sprintf(
        paste0(
          "PT: %.2f<br>Dominant: %s<br>",
          "n: %d (ND %d, Aab+ %d, T1D %d)<br>",
          "Normalized density: ND %.1f%%, Aab+ %.1f%%, T1D %.1f%%<br>",
          "%s mean: %.3f"
        ),
        hm$x, ifelse(is.na(hm$dominant), "n/a", hm$dominant),
        hm$n, hm$n_ND, hm$n_Aab, hm$n_T1D,
        100 * hm$dens_ND, 100 * hm$dens_Aab, 100 * hm$dens_T1D,
        feature_name, hm$feat_mean
      )

      g <- ggplot(hm, aes(x = x, y = 0.5, fill = fill, text = hover)) +
        geom_tile(width = bin_w, height = 1, na.rm = TRUE) +
        scale_fill_identity(guide = "none") +
        scale_x_continuous(limits = hm_lim, expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        guides(fill = "none") +
        theme_void() +
        theme(legend.position = "none",
              plot.margin = margin(5, 20, 10, 60, "pt"))

      ggplotly(g, tooltip = "text") %>%
        layout(
          showlegend = FALSE,
          xaxis = list(fixedrange = TRUE, showticklabels = FALSE,
                       zeroline = FALSE, range = hm_lim),
          yaxis = list(fixedrange = TRUE, showticklabels = FALSE,
                       zeroline = FALSE, range = c(0, 1))
        ) %>%
        config(displayModeBar = FALSE)
    })

    # ---------- Multi-Feature Heatmap ----------

    # Marker selector UI
    output$multi_heatmap_marker_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$var_names)) return(NULL)
      checkboxGroupInput(ns("multi_heatmap_markers"), "Select markers:",
                         choices = tr$var_names,
                         selected = intersect(tr$var_names,
                           c("INS", "GCG", "SST", "CHGA", "CD45", "CD31", "panCK", "B3TUBB")),
                         inline = TRUE)
    })

    # Multi-feature data reactive
    traj_multi_feature_data <- reactive({
      markers <- input$multi_heatmap_markers
      if (is.null(markers) || length(markers) == 0) return(NULL)
      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$adata)) return(NULL)
      valid <- intersect(markers, tr$var_names)
      if (length(valid) == 0) return(NULL)
      idx <- match(valid, tr$var_names)
      expr <- tryCatch(as.matrix(tr$adata$X[, idx, drop = FALSE]), error = function(e) NULL)
      if (is.null(expr)) return(NULL)
      colnames(expr) <- valid
      result <- data.frame(
        pt_raw = as.numeric(tr$obs$pseudotime),
        pt_rank = if (!is.null(tr$obs$pt_rank)) as.numeric(tr$obs$pt_rank) else NA_real_,
        total_cells = as.integer(tr$obs$total_cells),
        stringsAsFactors = FALSE)
      result <- cbind(result, as.data.frame(expr))
      result <- result[is.finite(result$pt_raw), , drop = FALSE]

      # Enforce the fitted cell-count cohort (matches traj_data_clean).
      cr <- applied_cell_range()
      {
        result <- result[is.finite(result$total_cells) &
                         result$total_cells >= cr[1] &
                         result$total_cells <= cr[2], , drop = FALSE]
        if (nrow(result) == 0) return(NULL)
      }

      # Use the SAME stable axis as the scatter above, via the SAME helper. This panel used
      # to re-rank over its own row set, which silently produced a different x-axis from the
      # scatter it sits under (its denominator excluded the per-feature !is.na(value) filter
      # that traj_data_clean applies).
      result$pt <- result$pt_raw
      result <- apply_pt_scale(result, input$pt_scale %||% "rank")
      result$total_cells <- NULL
      result$pt_raw <- NULL
      result$pt_rank <- NULL

      list(data = result, markers = valid)
    })

    # The z-scored bin x marker matrix behind the multi-feature heatmap.
    # Extracted from the plot reactive so `dl_multi_heatmap_csv` exports the
    # exact numbers rendered (same binning, same clamp) instead of recomputing.
    traj_multi_heatmap_matrix <- reactive({
      mfd <- traj_multi_feature_data()
      if (is.null(mfd)) return(NULL)

      df <- mfd$data
      markers <- mfd$markers
      nb <- input$multi_heatmap_nbins %||% 20

      # Bin pseudotime
      pt_range <- range(df$pt, na.rm = TRUE)
      df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])
      brks <- seq(0, 1, length.out = nb + 1)
      df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(df$pt_bin)

      # Compute per-bin mean for each marker, then z-score
      all_mids <- head(brks, -1) + diff(brks)[1] / 2
      hm_rows <- list()

      for (mk in markers) {
        bin_means <- numeric(length(bin_levels))
        bin_counts <- integer(length(bin_levels))
        for (i in seq_along(bin_levels)) {
          bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], mk]
          bin_counts[i] <- length(bin_data)
          bin_means[i] <- if (length(bin_data) >= 1) mean(bin_data, na.rm = TRUE) else NA_real_
        }
        # Z-score across bins
        valid_means <- bin_means[is.finite(bin_means)]
        if (length(valid_means) < 2) next
        mu <- mean(valid_means)
        sd_val <- sd(valid_means)
        if (!is.finite(sd_val) || sd_val == 0) sd_val <- 1
        z <- (bin_means - mu) / sd_val
        z <- pmax(-2.5, pmin(2.5, z))  # clamp

        for (i in seq_along(bin_levels)) {
          hm_rows[[length(hm_rows) + 1]] <- data.frame(
            marker = mk,
            x = all_mids[i] * (pt_range[2] - pt_range[1]) + pt_range[1],
            z = z[i],
            bin_mean = bin_means[i],
            n_islets = bin_counts[i],
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(hm_rows) == 0) return(NULL)
      hm <- do.call(rbind, hm_rows)

      # Order markers: hormones first, then immune, then others
      hormone_order <- intersect(c("INS", "GCG", "SST"), markers)
      immune_order <- intersect(c("CD3e", "CD4", "CD8a", "CD20", "CD45", "CD68", "CD163", "HLADR"), markers)
      other_order <- setdiff(markers, c(hormone_order, immune_order))
      marker_order <- c(hormone_order, immune_order, other_order)
      hm$marker <- factor(hm$marker, levels = rev(marker_order))

      list(hm = hm, marker_order = marker_order, pt_range = pt_range, n_bins = nb)
    })

    # Multi-row heatmap renderPlot
    traj_multi_heatmap_plot <- reactive({
      mm <- traj_multi_heatmap_matrix()
      if (is.null(mm)) return(NULL)
      hm <- mm$hm; pt_range <- mm$pt_range

      ggplot(hm, aes(x = x, y = marker, fill = z)) +
        # Missing marker/bin combinations are intentionally blank. ggplot's
        # default warning reports each blank as a dropped row on every render.
        geom_tile(na.rm = TRUE) +
        scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                             midpoint = 0, limits = c(-2.5, 2.5),
                             name = "Z-score", na.value = "#f7f7f7",
                             oob = scales::squish) +
        scale_x_continuous(limits = c(pt_range[1], pt_range[2]), expand = c(0, 0)) +
        labs(x = "Pseudotime", y = "", title = "Multi-Feature Expression Along Pseudotime") +
        theme_islet() +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = ISLET_FONT_SIZES$dense_min)
        )
    })
    output$traj_multi_heatmap <- renderPlot({ traj_multi_heatmap_plot() },
      height = function() {
        n <- length(input$multi_heatmap_markers %||% character(0))
        # 44 px/row (was 36) so the dense-floor row labels have vertical room.
        max(260, 70 + n * 44)
      })
    output$dl_multi_heatmap <- downloadHandler(
      filename = function() islet_png_filename("trajectory_multi_heatmap"),
      content = function(file) {
        p <- traj_multi_heatmap_plot(); if (is.null(p)) return()
        n <- length(input$multi_heatmap_markers %||% character(0))
        h <- max(6, 1 + n * 0.4)
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_HEATMAP_PNG_DIMS$width, height = h,
                        dpi = ISLET_HEATMAP_PNG_DIMS$dpi, units = ISLET_HEATMAP_PNG_DIMS$units)
      }
    )

    # The heatmap's underlying matrix. `z` is the plotted (clamped) value;
    # `bin_mean` is the raw per-bin mean it was standardised from, so the
    # z-scoring is reproducible from the file.
    output$dl_multi_heatmap_csv <- downloadHandler(
      filename = function() paste0("trajectory_multi_heatmap_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        mm <- tryCatch(traj_multi_heatmap_matrix(), error = function(e) NULL)
        out <- NULL
        if (!is.null(mm) && !is.null(mm$hm)) {
          out <- data.frame(
            marker           = as.character(mm$hm$marker),
            pseudotime_mid   = mm$hm$x,
            bin_mean         = mm$hm$bin_mean,
            z_clamped        = mm$hm$z,
            n_islets_in_bin  = mm$hm$n_islets,
            stringsAsFactors = FALSE)
          out <- out[order(match(out$marker, mm$marker_order), out$pseudotime_mid), , drop = FALSE]
        }
        pr <- trajectory_provenance()
        write_islet_csv(file, out, islet_csv_header(
          "Trajectory — multi-feature expression along pseudotime",
          pr$method,
          pr$scale,
          pr$range,
          paste("Markers:", paste(input$multi_heatmap_markers %||% "(none)", collapse = ", ")),
          paste("Pseudotime bins:", if (is.null(mm)) NA else mm$n_bins),
          "z_clamped = (bin_mean - mean of that marker's bin means) / sd, clamped to [-2.5, 2.5]. Z-scoring is ACROSS BINS within a marker, not across islets."
        ))
      }
    )

    # UMAP coordinates + donor status + the selected feature, i.e. everything
    # the two UMAP panels encode.
    output$dl_umap_csv <- downloadHandler(
      filename = function() paste0("trajectory_umap_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        td <- tryCatch(traj_data_clean(), error = function(e) NULL)
        out <- NULL
        if (!is.null(td) && is.data.frame(td) && nrow(td) > 0) {
          want <- c("case_id", "islet_key", "combined_islet_id", "donor_status",
                    "umap_1", "umap_2", "pt", "pt_rank",
                    "feature_name", "value", "total_cells", "islet_diam_um")
          keep <- intersect(want, colnames(td))
          if (length(keep) > 0) {
            out <- td[, keep, drop = FALSE]
            # `pt` / `value` are opaque outside the app.
            names(out)[names(out) == "pt"]    <- "pseudotime"
            names(out)[names(out) == "value"] <- "feature_value"
          }
        }
        pr <- trajectory_provenance()
        write_islet_csv(file, out, islet_csv_header(
          "Trajectory — UMAP coordinates and per-islet trajectory values",
          pr$method,
          pr$scale,
          pr$measure,
          pr$range,
          paste("Rows:", if (is.null(out)) 0 else nrow(out), "islets")
        ))
      }
    )


    # =====================================================================
    # PAGA force-directed layout (Barlow et al. 2024, Fig 2B/2E)
    # =====================================================================
    paga_data <- reactive({
      cr <- applied_cell_range()
      if (!paga_available()) return(NULL)
      withProgress(
        message = "Computing PAGA layout",
        detail = sprintf("Leiden + ForceAtlas2, %g-%g cells/islet", cr[1], cr[2]),
        value = 0.4,
        paga_result(cr[1], cr[2]))
    })

    output$paga_note <- renderUI({
      if (!paga_available()) {
        return(div(style = paste0("padding: 7px 9px; border-radius: 4px; font-size: 11px;",
                                  " line-height: 1.45; color: #8a1c12; background: #fdecea;",
                                  " border-left: 3px solid #d9534f;"),
          tags$strong("PAGA unavailable. "),
          "reticulate's interpreter has no scanpy. Install it into the app venv ",
          "(", tags$code("pip install scanpy igraph leidenalg fa2-modified"), ")."))
      }
      d <- paga_data()
      if (is.null(d)) {
        return(div(style = "font-size: 11px; color: #999; padding-top: 6px;",
          "No layout for this cell-count range."))
      }
      div(style = paste0("padding: 7px 9px; border-radius: 4px; font-size: 11px;",
                         " line-height: 1.45; color: #31556e; background: #eaf3f8;",
                         " border-left: 3px solid #4477AA;"),
        sprintf("%s islets, %d Leiden clusters (resolution 1), cosine 15-NN, %s layout. ",
                format(nrow(d), big.mark = ","), attr(d, "n_clusters") %||% NA_integer_,
                if (identical(attr(d, "layout_key"), "X_draw_graph_fa"))
                  "ForceAtlas2" else sub("^X_draw_graph_", "", attr(d, "layout_key") %||% "")),
        "Axes are layout coordinates with no units — only relative position means ",
        "anything. Rebuilt when the cell-count range changes.")
    })

    output$paga_layout <- renderPlotly({
      d <- paga_data()
      if (is.null(d) || !nrow(d)) {
        return(plotly_islet_empty(
          if (paga_available()) "No PAGA layout for this cell-count range"
          else "PAGA unavailable — scanpy not installed for reticulate"))
      }
      by <- input$paga_color_by %||% "pseudotime"
      psz <- input$paga_point_size %||% 1.6

      # "Selected feature" joins the Y-axis marker from the main scatter, so the
      # layout can be read the way Barlow's Fig 2E is.
      feat_label <- NULL
      if (identical(by, "feature")) {
        td <- tryCatch(traj_data_clean(), error = function(e) NULL)
        if (!is.null(td) && "combined_islet_id" %in% names(td)) {
          d$feature <- td$value[match(d$islet_id, td$combined_islet_id)]
          feat_label <- td$feature_name[1]
        }
        if (is.null(d$feature) || all(is.na(d$feature))) {
          by <- "pseudotime"   # nothing to join; fail visible-but-usable
        }
      }

      d$hover <- paste0("Islet: ", d$islet_id,
                        "<br>Donor: ", d$case_id, " (", d$donor_status, ")",
                        "<br>Pseudotime: ", round(d$pseudotime, 4),
                        "<br>Leiden: ", d$leiden,
                        if (!is.null(d$feature)) paste0("<br>", feat_label, ": ",
                                                        round(d$feature, 4)) else "")

      g <- ggplot(d, aes(x = fa_1, y = fa_2, text = hover))
      if (by %in% c("pseudotime", "feature")) {
        g <- g + geom_point(aes(colour = .data[[by]]), size = psz, alpha = 0.85, stroke = 0) +
          scale_colour_viridis_c(option = "inferno",
                                 name = if (identical(by, "pseudotime")) "Pseudotime"
                                        else feat_label %||% "Feature")
      } else if (identical(by, "donor_status")) {
        cols <- donor_colors_reactive()
        d$donor_status <- factor(d$donor_status, levels = names(cols))
        g <- ggplot(d, aes(x = fa_1, y = fa_2, text = hover)) +
          geom_point(aes(colour = donor_status), size = psz, alpha = 0.85, stroke = 0) +
          scale_colour_manual(values = cols, name = "Donor status", drop = FALSE)
      } else {
        g <- g + geom_point(aes(colour = .data[[by]]), size = psz, alpha = 0.85, stroke = 0) +
          scale_colour_discrete(name = if (identical(by, "leiden")) "Leiden" else "Donor")
      }
      g <- g +
        labs(x = "Layout dim 1", y = "Layout dim 2",
             title = {cr <- applied_cell_range()
                      sprintf("PAGA force-directed layout — %s islets (%g–%g cells)",
                              format(nrow(d), big.mark = ","), cr[1], cr[2])}) +
        theme_islet() +
        # Layout coordinates are unitless and isotropic: 1 unit across must equal
        # 1 unit up, or the graph structure is visually sheared.
        coord_fixed() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.minor = element_blank())

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = "text") %>%
        layout(font = fonts$global,
               title = list(font = fonts$title),
               xaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick),
               yaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick),
               legend = list(font = fonts$legend_text,
                             title = list(font = fonts$legend_title)),
               margin = list(t = 70)) %>%
        plotly_islet_config("trajectory_paga_layout")
    })

    output$dl_paga_csv <- downloadHandler(
      filename = function() sprintf(
        "paga_layout_%g-%g_%s.csv",
        applied_cell_range()[1],
        applied_cell_range()[2],
        format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        d <- paga_data()
        out <- if (is.null(d)) data.frame() else
          d[, c("islet_id", "case_id", "donor_status", "leiden",
                "fa_1", "fa_2", "pseudotime")]
        write_islet_csv(file, out, islet_csv_header(
          "Trajectory — PAGA force-directed layout (Barlow et al. 2024)",
          "Features: core + 20um peri cell-type counts / endocrine cells, log1p",
          sprintf("Neighbors: cosine, %d NN | Leiden resolution %g | layout %s",
                  15L, 1.0, attr(d, "layout_key") %||% "(n/a)"),
          "Root: ND islet nearest the ND centroid in diffusion space",
          {cr <- applied_cell_range()
           sprintf("Cells per islet range: %g-%g (rebuilds the embedding, not a display filter)",
                   cr[1], cr[2])},
          "fa_1 / fa_2 are unitless layout coordinates; only relative position is meaningful",
          paste("Rows:", nrow(out), "islets")
        ))
      }
    )

  })
}
