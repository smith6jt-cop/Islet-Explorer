# Pseudotime x-axis extent, derived in ONE place so the scatter and the strip heatmap
# stacked directly beneath it always share an x-extent. Tied x positions can spread
# within the domain by up to 1.2% of its width. Both panels reserve the same edge margin.
#
# Returns:
#   dom    — the data domain, used for BINNING (bins must cover data, not padding)
#   jitter — the jitter half-width
#   lim    — dom expanded by `jitter`, used for AXIS LIMITS in both panels
traj_pt_extent <- function(pt, jitter_frac = 0.012) {
  # suppressWarnings: an all-NA pt column makes range() emit "no non-missing arguments"
  # and return c(Inf, -Inf); the is.finite guard below is the intended handling.
  dom <- suppressWarnings(range(pt[is.finite(pt)]))
  if (!all(is.finite(dom))) dom <- c(0, 1)
  if (diff(dom) == 0) {
    half <- max(0.05, abs(dom[1]) * 0.05)
    dom <- c(max(0, dom[1] - half), dom[2] + half)
  }
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
traj_pseudotime_column <- function(obs, cell_range = c(TRAJ_MIN_CELLS, Inf)) {
  lo <- cell_range[1]; hi <- cell_range[2]
  # Every object with at least TRAJ_MIN_CELLS core cells is fitted, whatever its class; the rest get no
  # pseudotime. scanpy DPT is primary. The destiny implementation is the fallback if scanpy
  # is unreachable from reticulate. PAGA/ForceAtlas2 is the two-dimensional graph
  # card; it is not the pseudotime algorithm itself.
  bp <- compute_scanpy_pseudotime(lo, hi)
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
  attr(out, "features") <- attr(bp, "features")
  attr(out, "root_islet_id") <- attr(bp, "root_islet_id")
  attr(out, "root_method") <- attr(bp, "root_method")
  attr(out, "n_without_inputs") <- attr(bp, "n_without_inputs")
  out
}

#' One provenance line naming the graph inputs a pseudotime was fitted on, by region, in display names.
#' `features` are the fitted input names ("core_Beta cell", "peri_Exocrine", ...) the engine reports.
trajectory_inputs_text <- function(features) {
  features <- as.character(features %||% character(0))
  features <- features[!is.na(features) & grepl("^(core|peri|combined)_", features)]
  if (!length(features)) return("Trajectory inputs: not reported by the pseudotime engine")
  region <- sub("_.*$", "", features)
  label <- phenotype_text(sub("^(core|peri|combined)_", "", features))
  part <- function(key, title) {
    if (!any(region == key)) return(NULL)
    sprintf("%s (%d): %s", title, sum(region == key), paste(label[region == key], collapse = ", "))
  }
  paste0("Trajectory inputs: log1p(count or summed overlap / core endocrine cells) per class | ",
         paste(c(part("core", "core"), part("peri", "20 µm peri-islet band"), part("combined", "core + peri")), collapse = " | "),
         " | endocrine classes in the core only, exocrine in the peri-islet band only,",
         " macrophages by Gal3; Conflict, Doublet suspected, Tissue, Unclassified and",
         " unresolved-lineage cells are not inputs")
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

# Both PAGA panels use one fitted layout, with independent display colorings.
PAGA_COLOR_CHOICES <- c("Pseudotime" = "pseudotime", "Donor Status" = "donor_status",
                       "Leiden Cluster" = "leiden", "Donor" = "case_id",
                       "Selected Feature" = "feature")

# Join a selected feature by immutable object identity; coloring must never
# refit, filter or reorder the PAGA cohort. Missing feature values stay missing.
paga_display_data <- function(d, colour_by, feature_data = NULL, clinical_fields = character()) {
  if (is.null(d) || !nrow(d)) return(NULL)
  clinical_fields <- clinical_fields[clinical_fields %in% names(d)]
  choices <- c(PAGA_COLOR_CHOICES, clinical_fields)
  if (!is.character(colour_by) || length(colour_by) != 1L ||
      is.na(colour_by) || !(colour_by %in% unname(choices))) colour_by <- "pseudotime"
  d$.key <- trajectory_link_key(d$case_id, sub("^[0-9]+_", "", d$islet_id))
  title <- names(choices)[match(colour_by, choices)]
  if (identical(colour_by, "feature")) {
    d$.colour <- rep(NA_real_, nrow(d))
    if (!is.null(feature_data) && all(c("combined_islet_id", "value") %in% names(feature_data))) {
      d$.colour <- feature_data$value[match(d$islet_id, feature_data$combined_islet_id)]
      if (nrow(feature_data) && "feature_name" %in% names(feature_data) &&
          !is.na(feature_data$feature_name[1]) && nzchar(feature_data$feature_name[1]))
        title <- feature_data$feature_name[1]
    }
  } else {
    d$.colour <- d[[colour_by]]
  }
  guided <- isTRUE(attr(d, "layout_age_weight") > 0)
  if (guided && identical(colour_by, "pseudotime")) {
    d$.colour <- d$age_guided_pseudotime
    title <- "Age-guided pseudotime"
  }
  if (colour_by %in% c("donor_status", "leiden", "case_id")) d$.colour <- as.character(d$.colour)
  d$.hover <- paste0("Object: ", d$islet_id,
                     "<br>Donor: ", d$case_id, " (", d$donor_status, ")",
                     if (guided) "<br>Baseline DPT: " else "<br>Pseudotime: ", round(d$pseudotime, 4),
                     if (guided) paste0("<br>Age-guided pseudotime: ", round(d$age_guided_pseudotime, 4)) else "",
                     "<br>Leiden: ", d$leiden,
                     if (identical(colour_by, "feature")) paste0("<br>", title, ": ",
                       ifelse(is.finite(d$.colour), round(d$.colour, 4), "Unavailable")) else "")
  root <- which(d$is_dpt_root %in% TRUE)
  if (length(root)) d$.hover[root] <- paste0(d$.hover[root], "<br>DPT root")
  if (colour_by %in% clinical_fields) {
    shown <- trajectory_donor_display(d, colour_by)
    d$.hover <- paste0(d$.hover, "<br>", htmltools::htmlEscape(title), ": ", htmltools::htmlEscape(shown))
    if (is.numeric(d$.colour)) {
      if (colour_by == "C_peptide_ng_mL" && "C_peptide_censored" %in% names(d))
        d$.colour[d$C_peptide_censored %in% TRUE] <- NA_real_
      if (length(unique(d$.colour[is.finite(d$.colour)])) <= 1L) d$.colour <- shown
    } else {
      d$.colour <- shown
    }
  }
  list(data = d, colour_by = colour_by, title = title)
}

# Map the selected measure on the displayed islets. The asinh scale keeps zero and
# the full upper tail while resolving differences among small positive measurements.
trajectory_umap_colour <- function(df) {
  df <- df[is.finite(df$umap_1) & is.finite(df$umap_2) & is.finite(df$value), , drop = FALSE]
  values <- sort(unique(df$value))
  if (!length(values)) return(list(data = df))
  if (length(values) == 1L) {
    label <- format(values, digits = 5, trim = TRUE)
    df$.colour <- factor(rep(label, nrow(df)))
    return(list(data = df, palette = setNames("#8c7188", label)))
  }
  unit <- stats::median(abs(df$value[df$value != 0]))
  df$.colour <- asinh(df$value / unit)
  limits <- range(df$.colour)
  ticks <- if (length(values) <= 5) asinh(values / unit) else seq(limits[1], limits[2], length.out = 5)
  list(data = df, limits = limits, ticks = ticks,
       labels = format(signif(unit * sinh(ticks), 4), trim = TRUE))
}

trajectory_donor_summary <- function(obs, pt_scale) {
  axis <- if (identical(pt_scale, "raw")) obs$pseudotime else obs$pt_rank
  keep <- is.finite(axis) & !is.na(obs$imageid)
  if (!any(keep)) return(data.frame())
  indices <- split(which(keep), as.character(obs$imageid[keep]))
  do.call(rbind, lapply(indices, function(i) {
    q <- unname(stats::quantile(axis[i], c(0.25, 0.5, 0.75)))
    data.frame(case_id = as.character(obs$imageid[i[1]]), donor_status = as.character(obs$donor_status[i[1]]),
               n_islets = length(i), first_quartile = q[1], median_pt = q[2], third_quartile = q[3])
  }))
}

trajectory_server <- function(id, prepared, selected_islet, forced_image, active_tab = reactive("Trajectory"),
                              donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                              remove_outliers = reactive(FALSE),
                              outlier_threshold = reactive(3.0), macrophage_markers = reactiveVal(character(0)),
                              phenotype_palette_name = reactive("High Contrast")) {
  # `remove_outliers` / `outlier_threshold` are accepted but UNUSED: app.R still
  # passes them to every analysis module, and changing that signature is a wider
  # edit than this tab. The Trajectory tab no longer exposes or applies them.
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    bind_macrophage_marker_inputs(input, session, "mac_markers", macrophage_markers)

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
    donor_metadata <- reactive(trajectory_donor_metadata())
    traj_load_attempted <- reactiveVal(FALSE)
    # The fixed cohort includes eligible objects from both physical classes.
    applied_cell_range <- function() c(TRAJ_MIN_CELLS, Inf)
    cohort_label <- function() sprintf("objects with %d or more total core cells", TRAJ_MIN_CELLS)

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
          "case_id", "imageid", "age", "gender", "total_cells", "total_cells_core"),
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
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
        obs2$islet_key <- stringr::str_extract(obs2$combined_islet_id, "(?<=_)(Islet|Cluster|Single)_[0-9]+$")
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

      # Display choices do not change the cohort used to build the diffusion map.
      rng <- applied_cell_range()
      pt <- withProgress(
        # Seconds to tens of seconds depending on how many islets are in range.
        # Without a progress bar the tab simply freezes and reads as a hang.
        message = "Computing pseudotime",
        detail = sprintf("object-composition diffusion map, %s", cohort_label()),
        value = 0.4,
        traj_pseudotime_column(obs, rng))
      pseudotime_engine <- attr(pt, "engine") %||% "unavailable"
      pseudotime_n <- as.integer(attr(pt, "n_islets") %||% 0L)
      pseudotime_features <- attr(pt, "features")
      pseudotime_without_inputs <- as.integer(attr(pt, "n_without_inputs") %||% 0L)
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
           root_islet_id = attr(pt, "root_islet_id"), root_method = attr(pt, "root_method"),
           pseudotime_n = pseudotime_n, pseudotime_features = pseudotime_features,
           pseudotime_without_inputs = pseudotime_without_inputs, cell_range = rng)
    })

    output$donor_annotation_field <- renderUI({
      choices <- TRAJ_DONOR_FIELDS[TRAJ_DONOR_FIELDS %in% names(donor_metadata())]
      if (!length(choices)) return(NULL)
      selectInput(ns("donor_annotation"), "Donor annotation", choices = choices,
                  selected = if ("age" %in% choices) "age" else unname(choices[1]))
    })
    output$donor_annotation_plot <- renderPlotly({
      field <- input$donor_annotation
      if (!islet_valid_scalar_choice(field, TRAJ_DONOR_FIELDS, "trajectory_annotation_invalid"))
        return(plotly_islet_empty("Donor annotations unavailable"))
      scale <- input$pt_scale %||% "rank"
      rows <- trajectory_donor_summary(traj_with_pseudotime()$obs, scale)
      if (!nrow(rows)) return(plotly_islet_empty("No fitted islets"))
      rows <- trajectory_join_donor_metadata(rows, donor_metadata())
      if (!field %in% names(rows)) return(plotly_islet_empty("Donor annotation unavailable"))
      label <- names(TRAJ_DONOR_FIELDS)[match(field, TRAJ_DONOR_FIELDS)]
      display <- trajectory_donor_display(rows, field)
      censored <- rep(FALSE, nrow(rows))
      if (identical(field, "C_peptide_ng_mL") && "C_peptide_censored" %in% names(rows)) {
        censored <- rows$C_peptide_censored %in% TRUE & is.finite(rows[[field]])
        display[censored] <- paste("Below detection limit:", display[censored])
      }
      continuous <- is.numeric(rows[[field]]) && !startsWith(field, "AAb_")
      rows$annotation <- if (continuous) rows[[field]] else display
      rows$symbol <- ifelse(censored, "triangle-down-open", "circle")
      rows$hover <- paste0("Donor: ", htmltools::htmlEscape(rows$case_id), " (", htmltools::htmlEscape(rows$donor_status), ")",
                           "<br>", htmltools::htmlEscape(label), ": ", htmltools::htmlEscape(display),
                           "<br>Median pseudotime: ", signif(rows$median_pt, 3),
                           "<br>Islets: ", rows$n_islets)
      if (continuous) rows <- rows[is.finite(rows$annotation), , drop = FALSE]
      if (!nrow(rows)) return(plotly_islet_empty("No donors have this measurement"))
      colours <- donor_colors_reactive()
      p <- plotly::plot_ly()
      for (status in intersect(names(colours), unique(rows$donor_status))) {
        d <- rows[rows$donor_status == status, , drop = FALSE]
        p <- plotly::add_markers(p, x = d$median_pt, y = d$annotation, name = status,
                                 marker = list(color = unname(colours[status]), size = 10, symbol = d$symbol),
                                 error_x = list(type = "data", symmetric = FALSE,
                                                array = d$third_quartile - d$median_pt,
                                                arrayminus = d$median_pt - d$first_quartile),
                                 text = d$hover, hoverinfo = "text")
      }
      p <- plotly::layout(p, xaxis = list(title = if (scale == "raw") "Median raw DPT" else "Median pseudotime rank"),
                         yaxis = list(title = label, type = if (continuous) "linear" else "category",
                                      automargin = TRUE))
      plotly_islet_config(p, "trajectory_donor_annotations")
    })

    selected_measure_label <- reactive({
      if (identical(input$traj_y_mode %||% "ratio", "ratio")) {
        region <- switch(input$traj_ratio_region %||% "core",
                         core = "core", peri = "peri", "core + 20 µm peri")
        sprintf("%s — log1p(cells / endocrine cells, %s)",
                phenotype_text(input$traj_ratio_type %||% "(none)"), region)
      } else {
        sprintf("%s — normalized mean marker expression",
                input$traj_feature %||% "(none)")
      }
    })

    trajectory_provenance <- reactive({
      tr <- tryCatch(traj_with_pseudotime(), error = function(e) NULL)
      engine <- if (is.null(tr) || !is.null(tr$error)) "unavailable" else
        tr$pseudotime_engine %||% "unavailable"
      n <- if (is.null(tr) || !is.null(tr$error)) 0L else tr$pseudotime_n %||% 0L
      list(
        method = paste0("Pseudotime method: ", engine,
                        if (length(tr$root_islet_id)) paste0(" | Root: ", tr$root_islet_id, " | ", tr$root_method)),
        features = trajectory_inputs_text(
          if (is.null(tr) || !is.null(tr$error)) NULL else tr$pseudotime_features),
        scale = sprintf("Pseudotime scale: %s",
                        if (identical(input$pt_scale %||% "rank", "raw"))
                          "raw diffusion pseudotime" else "rank-normalized within fitted cohort"),
        range = sprintf("Fitted cohort: %s (%s objects)",
                        cohort_label(), format(n, big.mark = ",")),
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
        region <- input$traj_ratio_region %||% "core"
        types <- phenotype_selectable(barlow_ratio_types(macrophage_markers(), region))
        if (!length(types)) {
          return(selectInput(ns("traj_ratio_type"), "Select cell type:",
                             choices = c("Composition features unavailable" = ""),
                             selected = ""))
        }
        prev <- isolate(input$traj_ratio_type)
        sel <- if (!is.null(prev) && prev %in% types) prev
               else if (!is.null(prev) && is_macrophage_label(prev) && any(is_macrophage_label(types)))
                 types[is_macrophage_label(types)][1]
               else if ("Beta cell" %in% types) "Beta cell" else types[1]
        # The macrophage markers (owner, 2026-10-27) sit under the cell type; they set the app-wide setting.
        return(tagList(
          selectInput(ns("traj_ratio_type"), "Select cell type:",
                      choices = phenotype_grouped_choices(types), selected = sel),
          if (root_has_marker_profile_labels())
            macrophage_marker_checkboxes(ns("mac_markers"), isolate(macrophage_markers()))))
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
      # the X-axis (pseudotime) instead, then visualize the associated
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
        region <- input$traj_ratio_region %||% "core"
        if (!islet_valid_scalar_choice(region, BARLOW_RATIO_REGIONS,
                                       "trajectory_region_invalid")) return(NULL)
        if (!islet_valid_scalar_choice(ct, barlow_ratio_types(macrophage_markers(), region),
                                       "trajectory_ratio_invalid")) return(NULL)
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
        # total_cells_core is the engine cutoff. Legacy core H5ADs store the count in total_cells.
        cells_core = as.integer(obs_df$total_cells_core %||% obs_df$total_cells),
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
      if ("cells_core" %in% colnames(result_df)) {
        tc <- suppressWarnings(as.integer(result_df$cells_core))
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
      result_df$link_key <- trajectory_link_key(result_df$case_id, result_df$islet_key)
      if (anyNA(result_df$link_key) || anyDuplicated(result_df$link_key)) {
        stop("Trajectory islets must have unique donor and islet identities")
      }

      return(result_df)
    })

    # ---------- Status bar ----------
    output$traj_status <- renderUI({
      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) {
        return(tags$div(style = "color:var(--islet-danger-fg);",
          "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
      }
      if (!requireNamespace("anndata", quietly = TRUE)) {
        return(tags$div(style = "color:var(--islet-danger-fg);",
          "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
      }
      tr <- traj_with_pseudotime()
      if (is.null(tr)) {
        return(tags$div(style = "color:var(--islet-danger-fg);",
          "AnnData object not available or H5AD failed to load."))
      }
      if (!is.null(tr$error)) {
        install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
        return(tags$div(style = "color:var(--islet-danger-fg);",
          sprintf("Trajectory load error: %s%s", tr$error, install_hint)))
      }
      nfit <- tr$pseudotime_n %||% 0L
      nfeat <- length(tr$pseudotime_features %||% character(0))
      nwithout <- tr$pseudotime_without_inputs %||% 0L
      src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur

      tags$div(style = "color: var(--islet-fg-muted); font-size: 12px; margin-bottom: 10px;",
        title = paste("Source:", basename(src)),
        sprintf("Endocrine objects: %s total | %s fitted (%d or more total core cells, any class)%s%s",
          format(nrow(tr$obs), big.mark = ","), format(nfit, big.mark = ","), TRAJ_MIN_CELLS,
          if (nfeat > 0) sprintf(" | %d composition inputs", nfeat) else "",
          if (nwithout > 0) sprintf(" | %s with no cell of an input class, not fitted",
                                    format(nwithout, big.mark = ",")) else "")
      )
    })

    # ---------- Trajectory scatter plot ----------
    # Built with plot_ly() directly (R/plot_figures.R, Plan A Level 2): the ggplot + ggplotly() path cost most of
    # the step's server time. Optional tie spreading affects only marker positions.
    # LOESS uses the original coordinates and log1p(cells) weights.
    output$traj_scatter <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) {
        tr <- traj()
        if (is.null(tr) || (is.null(tr$error) &&
            identical(input$traj_y_mode %||% "ratio", "ratio") &&
            is.null(input$traj_ratio_type))) {
          return(trajectory_link_plot(plotly_islet_empty("Loading data"), ns))
        }
        return(trajectory_link_plot(plotly_islet_empty("No trajectory data available"), ns))
      }

      color_by <- input$traj_color_by %||% "donor_status"
      # A Leiden column absent from df (e.g. the Excel fallback) colors by donor status instead.
      if (startsWith(color_by, "leiden_") && !color_by %in% colnames(df)) color_by <- "donor_status"
      if (color_by == "donor_id") {
        lv <- levels(factor(df$donor_id))
        colours <- if (length(lv) <= 12) {
          setNames(RColorBrewer::brewer.pal(min(12, max(3, length(lv))), "Paired")[seq_along(lv)], lv)
        } else setNames(grDevices::rainbow(length(lv), s = 0.8, v = 0.8), lv)
        colour_label <- "Donor ID"
      } else if (startsWith(color_by, "leiden_")) {
        # Levels in numeric order, so the hover reads 0, 1, 2, ... rather than 0, 1, 10, 11.
        lvl_vals <- unique(stats::na.omit(df[[color_by]]))
        num_vals <- suppressWarnings(as.numeric(lvl_vals))
        lv <- as.character(if (!any(is.na(num_vals))) lvl_vals[order(num_vals)] else sort(lvl_vals))
        colours <- setNames(if (length(lv) <= 12) {
          RColorBrewer::brewer.pal(max(3, min(12, length(lv))), "Set3")[seq_len(length(lv))]
        } else scales::hue_pal()(length(lv)), lv)
        colour_label <- paste("Leiden", sub("^leiden_", "res ", color_by))
      } else {
        color_by <- "donor_status"
        colours <- donor_colors_reactive()
        colour_label <- "Donor Status"
      }

      # Labels: the y-axis names the quantity, not just the feature (a ratio and a gate-rescaled marker mean are
      # different units on the same plot slot).
      if (identical(input$traj_y_mode %||% "ratio", "ratio")) {
        selected_feature <- phenotype_text(input$traj_ratio_type %||% "Selected cell type")
        reg <- input$traj_ratio_region %||% "core"
        metric_label <- sprintf("\u2014 log1p(cells per core endocrine cell, %s)",
                                switch(reg, core = "core", peri = "peri", "core + 20 \u00b5m peri"))
      } else {
        selected_feature <- input$traj_feature %||% "Selected feature"
        metric_label <- "\u2014 normalized expression (mean per islet, 0\u20131)"
      }
      x_label <- if ((input$pt_scale %||% "rank") == "raw") "Pseudotime (raw DPT)" else "Pseudotime (rank-normalized)"
      spread_ties <- isTRUE(input$traj_spread_ties %||% TRUE)
      y_label <- paste(selected_feature, metric_label)
      if (spread_ties) y_label <- paste0(y_label, "<br>(ties spread for display)")

      fig <- trajectory_scatter_figure(df, list(
        color_by = color_by, colours = colours, colour_label = colour_label,
        spread_ties = spread_ties,
        size_by = input$traj_point_size %||% "uniform", point_size = input$traj_point_size_slider %||% 3.0,
        alpha = 1 - (input$traj_alpha %||% 0.4), trend = input$traj_show_trend %||% "by_donor",
        # 95% CI bands on the LOESS trends (Bodenmiller 2025 Fig 2D/5E convention).
        show_ci = isTRUE(input$traj_trend_ci %||% TRUE), status_cols = donor_colors_reactive(), ink = islet_ink(),
        x_label = x_label, y_label = y_label,
        # The pseudotime extent is shared with the strip heatmap below (traj_pt_extent()), so the two stay aligned.
        pt_ext = traj_pt_extent(df$pt)))
      trajectory_link_plot(fig, ns)
    })

    # One validated selection route for the scatter, UMAPs and both PAGA panels.
    # PAGA keeps its complete fitted cohort when the selected Y feature is missing.
    linked_cohort_keys <- reactive({
      df <- traj_data_clean()
      keys <- if (is.null(df)) character(0) else df$link_key
      if (isTRUE(input$paga_show_layout)) {
        d <- paga_data()
        if (!is.null(d)) keys <- union(keys, trajectory_link_key(d$case_id, sub("^[0-9]+_", "", d$islet_id)))
      }
      keys[!is.na(keys)]
    })
    linked_selection_key <- function() {
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      key <- trajectory_link_key(info$case_id, info$islet_key)
      if (length(key) != 1L || is.na(key)) NULL else key
    }
    send_linked_state <- function() {
      session$sendCustomMessage(ns("linked_state"), list(
        key = linked_selection_key(), active = identical(active_tab(), "Trajectory")))
    }
    observeEvent(input$linked_click, {
      req(identical(active_tab(), "Trajectory"))
      event <- input$linked_click
      key <- event$key
      sources <- ns(c("traj_scatter", "traj_umap_donor", "traj_umap_feature", "paga_layout", "paga_layout_right"))
      if (!is.character(key) || length(key) != 1L || is.na(key) ||
          !is.character(event$source) || length(event$source) != 1L ||
          !(event$source %in% sources) || !(key %in% linked_cohort_keys())) {
        islet_log_validation_failure("trajectory_link_invalid")
        send_linked_state()
        return()
      }
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      select_islet_from_key(parts[1], parts[2], prepared, selected_islet)
      send_linked_state()
    }, ignoreInit = TRUE)
    observeEvent(input$linked_ready, { send_linked_state() }, ignoreInit = TRUE)
    observe({
      active <- identical(active_tab(), "Trajectory")
      key <- linked_selection_key()
      if (active && !is.null(key)) {
        df <- traj_data_clean()
        # A transient empty feature selector must not discard a pinned islet.
        # Once a feature is valid, retain any pin still present in a visible
        # fitted layout, including PAGA objects missing that feature.
        if ((!is.null(df) && !(key %in% linked_cohort_keys())) ||
            (is.null(df) && !isTRUE(input$paga_show_layout) && !is.null(traj_validate()))) {
          selected_islet(NULL)
        }
      }
      send_linked_state()
    })
    observeEvent(input$clear_link_selection, {
      selected_islet(NULL)
      send_linked_state()
    })

    # ---------- Segmentation viewer panel ----------
    output$segmentation_viewer_panel <- renderUI({
      # Only render when Trajectory tab is active to avoid duplicate non-namespaced output IDs
      if (active_tab() != "Trajectory") return(NULL)
      info <- selected_islet()
      if (is.null(info)) {
        return(div(style = "min-height: 275px; display: flex; align-items: center; justify-content: center; color: var(--islet-fg-muted);",
                   "No islet selected"))
      }
      segmentation_viewer_panel_ui(info, ns("clear_segmentation"), layout = "stacked",
                                    selected_palette = isolate(phenotype_palette_name()))
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

          # A compact responsive grid uses available width for several short
          # columns instead of making the donor legend one tall card section.
          items <- lapply(seq_along(donor_ids), function(i) {
            div(style = "display: flex; align-items: center; min-width: 0; gap: 5px;",
              div(style = sprintf("flex: 0 0 14px; width: 14px; height: 14px; background-color: %s; border-radius: 3px;", colors[i])),
              span(donor_ids[i], style = "font-size: 13px; white-space: nowrap;")
            )
          })
          div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(58px, 1fr)); gap: 5px 8px;", items)
        } else {
          p("No data available", style = "font-size: 12px; color: var(--islet-fg-subtle);")
        }
      } else if (startsWith(color_by, "leiden_")) {
        df <- traj_data_clean()
        if (is.null(df) || !(color_by %in% colnames(df))) {
          return(p("Leiden clusters unavailable (rebuild islet_explorer.h5ad).",
                   style = "font-size: 12px; color: var(--islet-fg-subtle);"))
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

    output$traj_umap_coverage <- renderText({
      df <- traj_data_clean()
      if (is.null(df)) return("")
      n <- sum(is.finite(df$umap_1) & is.finite(df$umap_2))
      cutoff <- prepared()$composition_embedding$min_core_cells
      paste0(sprintf("Saved UMAP coordinates for %s of %s displayed trajectory objects.",
                     format(n, big.mark = ","), format(nrow(df), big.mark = ",")),
             if (length(cutoff) == 1L && is.finite(as.numeric(cutoff)))
               sprintf(" The composition UMAP requires at least %d total core cells and a defined endocrine composition.",
                       as.integer(cutoff)))
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
      df <- df[islet_umap_order(nrow(df)), , drop = FALSE]  # mixed draw order, as on screen

      ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
        geom_point(alpha = 0.9, size = 1.0) +
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
    # The screen UMAPs are drawn by islet_umap_plotly() (one interleaved trace, outlined points) rather than
    # ggplotly(), which split the donor-status colors into traces drawn in turn so T1D covered the others. The
    # ggplot reactives above remain the PNG downloads.
    interactive_umap <- function(g, source, colour, pal = NULL, legend_title,
                                 df = traj_data_clean(), colorbar_ticks = list()) {
      if (is.null(g) || is.null(df) || !nrow(df)) {
        return(trajectory_link_plot(plotly_islet_empty("No trajectory data available"), ns))
      }
      df <- df[is.finite(df$umap_1) & is.finite(df$umap_2), , drop = FALSE]
      if (!nrow(df)) {
        return(trajectory_link_plot(plotly_islet_empty("UMAP coordinates not available"), ns))
      }
      df$link_hover <- paste0("Donor: ", df$case_id, "<br>Islet: ", df$islet_key,
        "<br>Status: ", df$donor_status, "<br>Pseudotime: ", signif(df$pt, 3),
        "<br>", selected_measure_label(), ": ", signif(df$value, 3))
      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
      fonts <- plotly_islet_fonts()
      axis <- function(title) list(title = list(text = title, font = fonts$axis_title), zeroline = FALSE,
                                   showgrid = FALSE, showticklabels = FALSE, ticks = "", showline = FALSE)
      p <- islet_umap_plotly(df, "umap_1", "umap_2", colour, pal, "link_hover", key = "link_key",
                             source = ns(source), colorbar_title = legend_title,
                             size = 3.5, colorbar_ticks = colorbar_ticks) %>%  # these plots are smaller than the Spatial UMAPs
        layout(margin = list(l = 40, r = 15, t = 10, b = 40),
               font = fonts$global,
               xaxis = axis("UMAP 1"),
               yaxis = c(axis("UMAP 2"), list(scaleanchor = "x", scaleratio = 1)),
               legend = list(font = fonts$legend_text,
                             title = list(text = legend_title, font = fonts$legend_title))) %>%
        plotly_islet_config(paste0("trajectory_", source))
      trajectory_link_plot(p, ns)
    }
    output$traj_umap_donor <- renderPlotly({
      interactive_umap(traj_umap_donor_plot(), "traj_umap_donor", "donor_status", donor_colors_reactive(),
                       "Donor status")
    })
    output$dl_umap_donor <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_donor"),
      content = function(file) {
        p <- traj_umap_donor_plot(); if (is.null(p)) return()
        islet_ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- UMAP: Selected Feature ----------
    output$traj_umap_feature_title <- renderText({
      card_heading_text(paste("UMAP:", selected_measure_label()))
    })
    feature_colours <- reactive({
      df <- traj_data_clean()
      if (!is.null(df)) trajectory_umap_colour(df)
    })
    traj_umap_feature_plot <- reactive({
      scale <- feature_colours()
      df <- scale$data
      if (is.null(df) || nrow(df) == 0) return(NULL)

      selected_feature <- selected_measure_label()
      color_label <- if (identical(input$traj_y_mode %||% "ratio", "ratio"))
        "Log cell ratio" else "Expression"

      df <- df[islet_umap_order(nrow(df)), , drop = FALSE]
      g <- ggplot(df, aes(x = umap_1, y = umap_2, color = .colour)) +
        geom_point(alpha = 0.9, size = 1.0) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2",
             title = paste("UMAP:", selected_feature),
             color = color_label) +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())

      if (!is.null(scale$palette)) return(g + scale_color_manual(values = scale$palette))
      g + scale_color_gradientn(colours = substr(viridisLite::inferno(11), 1, 7),
                                limits = scale$limits, breaks = scale$ticks, labels = scale$labels)
    })
    output$traj_umap_feature <- renderPlotly({
      scale <- feature_colours()
      ticks <- if (length(scale$ticks)) list(tickmode = "array", tickvals = scale$ticks,
                                             ticktext = scale$labels) else list()
      interactive_umap(traj_umap_feature_plot(), "traj_umap_feature", ".colour", scale$palette,
                       legend_title = if (identical(input$traj_y_mode %||% "ratio", "ratio"))
                         "Log cell ratio" else "Expression",
                       df = scale$data, colorbar_ticks = ticks)
    })
    output$dl_umap_feature <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_feature"),
      content = function(file) {
        p <- traj_umap_feature_plot(); if (is.null(p)) return()
        islet_ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- Heatmap: islet distribution by disease status ----------
    # Each row is one disease status; fill strength is the share of that
    # status's fitted islets in each pseudotime bin. Hover reports the raw
    # count and the status contribution within that bin.
    output$traj_heatmap <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

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
      status_totals <- vapply(status_levels, function(s) sum(df$donor_status == s, na.rm = TRUE), numeric(1))

      bin_rows <- lapply(seq_along(bin_levels), function(i) {
        bd <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], , drop = FALSE]
        counts <- vapply(status_levels, function(s) sum(bd$donor_status == s, na.rm = TRUE), integer(1))
        bin_total <- sum(counts)
        data.frame(idx = i, x = mids[i], status = status_levels,
                   n = as.integer(counts),
                   status_total = as.integer(status_totals[status_levels]),
                   status_share = ifelse(status_totals[status_levels] > 0,
                                         counts / status_totals[status_levels], 0),
                   bin_share = if (bin_total > 0) counts / bin_total else rep(0, length(status_levels)),
                   stringsAsFactors = FALSE)
      })
      hm <- do.call(rbind, bin_rows)
      if (is.null(hm) || nrow(hm) == 0) return(NULL)

      blend_to_white <- function(hex, t) {
        rgb_full <- grDevices::col2rgb(hex) / 255
        grDevices::rgb(red   = t * rgb_full[1] + (1 - t) * 1,
                       green = t * rgb_full[2] + (1 - t) * 1,
                       blue  = t * rgb_full[3] + (1 - t) * 1)
      }
      max_share <- ave(hm$status_share, hm$status, FUN = max)
      hm$intensity <- ifelse(max_share > 0, hm$status_share / max_share, 0)
      hm$fill <- mapply(function(status, intensity) {
        base <- donor_cols[[status]] %||% "#777777"
        blend_to_white(base, pmax(0, pmin(1, intensity)))
      }, hm$status, hm$intensity, USE.NAMES = FALSE)
      hm$hover <- sprintf(
        paste0("Pseudotime bin: %.3f<br>%s: %d of %d islets (%.1f%% of status)<br>",
               "Status contribution within bin: %.1f%%"),
        hm$x, hm$status, hm$n, hm$status_total, 100 * hm$status_share, 100 * hm$bin_share
      )
      hm$status <- factor(hm$status, levels = rev(status_levels))

      g <- ggplot(hm, aes(x = x, y = status, fill = fill, text = hover)) +
        geom_tile(width = bin_w * 0.98, height = 0.76, na.rm = TRUE) +
        scale_fill_identity(guide = "none") +
        scale_x_continuous(limits = hm_lim, expand = c(0, 0)) +
        scale_y_discrete(expand = expansion(add = 0.12)) +
        guides(fill = "none") +
        labs(x = NULL, y = NULL) + theme_void() +
        theme(legend.position = "none",
              axis.text.y = element_text(size = 9, color = islet_ink(), margin = margin(r = 5)),
              axis.text.x = element_blank(), axis.ticks = element_blank(),
              plot.margin = margin(2, 20, 2, 20, "pt"))

      ggplotly(g, tooltip = "text") %>%
        layout(
          showlegend = FALSE,
          xaxis = list(fixedrange = TRUE, showticklabels = FALSE,
                       zeroline = FALSE, range = hm_lim),
          yaxis = list(fixedrange = TRUE, zeroline = FALSE)
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
        cells_core = as.integer(tr$obs$total_cells_core %||% tr$obs$total_cells),
        stringsAsFactors = FALSE)
      result <- cbind(result, as.data.frame(expr))
      result <- result[is.finite(result$pt_raw), , drop = FALSE]

      # Enforce the fitted cell-count cohort (matches traj_data_clean).
      cr <- applied_cell_range()
      {
        result <- result[is.finite(result$cells_core) &
                         result$cells_core >= cr[1] &
                         result$cells_core <= cr[2], , drop = FALSE]
        if (nrow(result) == 0) return(NULL)
      }

      # Use the SAME stable axis as the scatter above, via the SAME helper. This panel used
      # to re-rank over its own row set, which silently produced a different x-axis from the
      # scatter it sits under (its denominator excluded the per-feature !is.na(value) filter
      # that traj_data_clean applies).
      result$pt <- result$pt_raw
      result <- apply_pt_scale(result, input$pt_scale %||% "rank")
      result$cells_core <- NULL
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
    output$traj_multi_heatmap <- renderPlot({ traj_multi_heatmap_plot() }, bg = "transparent",
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
        islet_ggsave(file, plot = p,
                        width = ISLET_HEATMAP_PNG_DIMS$width, height = h,
                        dpi = ISLET_HEATMAP_PNG_DIMS$dpi, units = ISLET_HEATMAP_PNG_DIMS$units)
      }
    )

    # The heatmap's underlying matrix. `z` is the plotted (clamped) value;
    # `bin_mean` is the raw per-bin mean it was standardized from, so the
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
          pr$features,
          pr$scale,
          pr$range,
          paste("Markers:", paste(mm$marker_order %||% "(none)", collapse = ", ")),
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
                    "feature_name", "value", "cells_core", "total_cells", "islet_diam_um")
          keep <- intersect(want, colnames(td))
          if (length(keep) > 0) {
            out <- td[, keep, drop = FALSE]
            out <- trajectory_join_donor_metadata(out, donor_metadata())
            # `pt` / `value` are opaque outside the app.
            names(out)[names(out) == "pt"]    <- "pseudotime"
            names(out)[names(out) == "value"] <- "feature_value"
            names(out)[names(out) == "cells_core"] <- "total_cells_core"
          }
        }
        pr <- trajectory_provenance()
        write_islet_csv(file, out, islet_csv_header(
          "Trajectory — UMAP coordinates and per-islet trajectory values",
          "Rows include the displayed pseudotime objects. Blank UMAP coordinates indicate no saved position.",
          pr$method,
          pr$features,
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
    paga_age_weight <- debounce(reactive({
      w <- input$paga_age_weight %||% 0
      req(paga_valid_age_weight(w))
      round(w, 2)
    }), 250)
    paga_data <- reactive({
      req(isTRUE(input$paga_show_layout))
      cr <- applied_cell_range()
      if (!paga_available()) return(NULL)
      withProgress(
        message = "Computing PAGA layout",
        detail = sprintf("Leiden + ForceAtlas2, %s", cohort_label()),
        value = 0.4,
        paga_result(cr[1], cr[2], age_weight = paga_age_weight()))
    })

    paga_annotated <- reactive({
      d <- paga_data()
      if (!is.null(d)) trajectory_join_donor_metadata(d, donor_metadata())
    })

    observe({
      fields <- TRAJ_DONOR_FIELDS[TRAJ_DONOR_FIELDS %in% names(donor_metadata())]
      choices <- list("Trajectory" = PAGA_COLOR_CHOICES)
      if (length(fields)) choices[["Donor demographics and clinical history"]] <- fields
      for (id in c("paga_color_by", "paga_color_right")) {
        current <- isolate(input[[id]])
        default <- if (id == "paga_color_by") "pseudotime" else "donor_status"
        updateSelectInput(session, id, choices = choices,
                          selected = if (!is.null(current) && current %in% unlist(choices)) current else default)
      }
    })

    output$paga_note <- renderUI({
      req(isTRUE(input$paga_show_layout))
      if (!paga_available()) {
        return(div(style = paste0("padding: 7px 9px; border-radius: 4px; font-size: 11px;",
                                  " line-height: 1.45; color: var(--islet-danger-fg); background: var(--islet-danger-bg);",
                                  " border-left: 3px solid #d9534f;"),
          tags$strong("PAGA unavailable. "),
          "reticulate's interpreter has no scanpy. Install it into the app venv ",
          "(", tags$code("pip install scanpy igraph leidenalg fa2-modified"), ")."))
      }
      d <- paga_data()
      if (is.null(d)) {
        return(div(style = "font-size: 11px; color: var(--islet-fg-subtle); padding-top: 6px;",
          if (paga_age_weight() > 0) "Age guidance unavailable for this cohort." else "No layout for this cell-count range."))
      }
      div(style = paste0("padding: 7px 9px; border-radius: 4px; font-size: 11px;",
                         " line-height: 1.45; color: var(--islet-info-fg); background: var(--islet-info-bg);",
                         " border-left: 3px solid #4477AA;"),
        sprintf("%s islets, %d Leiden clusters (resolution 1), cosine 15-NN, %s layout. ",
                format(nrow(d), big.mark = ","), attr(d, "n_clusters") %||% NA_integer_,
                if (identical(attr(d, "layout_key"), "X_draw_graph_fa"))
                  "ForceAtlas2" else sub("^X_draw_graph_", "", attr(d, "layout_key") %||% "")),
        sprintf("DPT root %s (ND). ", attr(d, "root_islet_id")),
        if (isTRUE(attr(d, "layout_age_weight") > 0))
          sprintf("Age guidance %.2f. The Pseudotime color shows the guided score.", attr(d, "layout_age_weight"))
        else "Age guidance 0. The Pseudotime color shows baseline DPT.",
        if (isTRUE(attr(d, "layout_age_weight") > 0) && isTRUE(attr(d, "age_direction_fallback")))
          " Donor centers give no age direction; shifts use layout dimension 1.",
        if (isTRUE(attr(d, "layout_age_weight") > 0) && isTRUE(attr(d, "n_missing_age_donors") > 0))
          sprintf(" %d donors lack a recorded age and retain their original positions.", attr(d, "n_missing_age_donors")))
    })

    paga_figure <- function(colour_by, output_id) {
      d <- paga_annotated()
      if (is.null(d) || !nrow(d)) {
        return(plotly_islet_empty(
          if (paga_available()) "No PAGA layout for this cell-count range"
          else "PAGA unavailable — scanpy not installed for reticulate"))
      }
      feature_data <- if (identical(colour_by, "feature"))
        tryCatch(traj_data_clean(), error = function(e) NULL) else NULL
      display <- paga_display_data(d, colour_by, feature_data, TRAJ_DONOR_FIELDS)
      d <- display$data
      if (identical(display$colour_by, "feature") && !any(is.finite(d$.colour)))
        return(plotly_islet_empty("No values for the selected feature in the PAGA cohort"))
      numeric_colour <- is.numeric(d$.colour)
      pal <- if (numeric_colour) NULL else if (identical(display$colour_by, "donor_status"))
        donor_colors_reactive() else {
          groups <- sort(unique(d$.colour))
          stats::setNames(scales::hue_pal()(length(groups)), groups)
        }
      # The shared UMAP renderer gives both panels the same interleaved draw
      # order and stable point keys. Only the color vector differs.
      missing_colour <- if (numeric_colour && display$colour_by %in% TRAJ_DONOR_FIELDS)
        !is.finite(d$.colour) else rep(FALSE, nrow(d))
      point_size <- plotly_point_px(input$paga_point_size %||% 1.6)
      p <- islet_umap_plotly(d[!missing_colour, , drop = FALSE], "fa_1", "fa_2", ".colour",
                            pal, ".hover", key = ".key", source = ns(output_id),
                            colorbar_title = display$title, size = point_size)
      if (any(missing_colour)) {
        # Keep every fitted islet selectable while reserving the gradient for measured values.
        m <- d[missing_colour, , drop = FALSE]
        p$x$data <- c(list(list(type = "scatter", mode = "markers", x = signif(m$fa_1, 5),
                               y = signif(m$fa_2, 5), key = m$.key, text = m$.hover, hoverinfo = "text",
                               showlegend = FALSE, marker = list(color = "#9d9d9d", size = point_size))),
                      p$x$data)
      }
      if (numeric_colour) {
        p <- plotly::plotly_build(p)
        p$x$data <- lapply(p$x$data, function(tr) {
          if (!is.null(tr$marker$colorbar))
            tr$marker$colorbar <- utils::modifyList(tr$marker$colorbar, list(
              orientation = "h", x = 0.5, xanchor = "center", y = -0.3, yanchor = "bottom",
              len = 0.8, thickness = 14, title = list(side = "top"),
              tickfont = list(size = ISLET_FONT_SIZES$dense_min)))
          tr
        })
      }
      fonts <- plotly_islet_fonts()
      limits <- attr(d, "paga_axis_limits")
      root <- which(d$is_dpt_root %in% TRUE)
      root_label <- if (length(root) == 1L) list(list(
        text = "DPT root", x = signif(d$fa_1[root], 5), y = signif(d$fa_2[root], 5),
        xref = "x", yref = "y", showarrow = TRUE, arrowhead = 2, arrowcolor = "#757575",
        ax = if (d$fa_1[root] < median(d$fa_1)) 37 else -37,
        ay = if (d$fa_2[root] < median(d$fa_2)) -24 else 24,
        font = list(size = 13, color = "#212121"), bgcolor = "#ffffff",
        bordercolor = "#757575", borderpad = 2, captureevents = FALSE))
      p <- p %>% layout(
        uirevision = paste0("paga-cohort-", cohort_label()),
        font = fonts$global,
        title = list(text = display$title, font = fonts$title),
        margin = list(l = 30, r = 10, t = 45, b = 110, autoexpand = FALSE),
        xaxis = list(range = limits$fa_1, title = list(text = "Layout Dimension 1", font = fonts$axis_title), domain = c(0, 1),
                     zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                     showline = FALSE, constrain = "domain"),
        yaxis = list(range = limits$fa_2, title = list(text = "Layout Dimension 2", font = fonts$axis_title),
                     zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, ticks = "",
                     showline = FALSE, scaleanchor = "x", scaleratio = 1, constrain = "domain"),
        legend = list(title = list(text = display$title, font = fonts$legend_title, side = "top"),
                      orientation = "h", font = fonts$legend_text,
                      x = 0.5, xanchor = "center", y = -0.3, yanchor = "bottom",
                      itemsizing = "constant", tracegroupgap = 2)
      ) %>% plotly_islet_config(paste0("trajectory_", output_id))
      # Keep the annotation's dark ink on its opaque white label in either plot theme.
      p <- plotly::layout(p, annotations = root_label)
      p <- trajectory_link_plot(p, ns)
      paired_plotly_guides(p, c(ns("paga_layout"), ns("paga_layout_right")))
    }
    output$paga_layout <- renderPlotly({
      paga_figure(input$paga_color_by %||% "pseudotime", "paga_layout")
    })
    output$paga_layout_right <- renderPlotly({
      paga_figure(input$paga_color_right %||% "donor_status", "paga_layout_right")
    })

    output$dl_paga_csv <- downloadHandler(
      filename = function() sprintf(
        "paga_layout_%s_%s.csv",
        sprintf("min%dcells", TRAJ_MIN_CELLS),
        format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        d <- paga_annotated()
        columns <- c("islet_id", "case_id", "donor_status", "total_cells_core", "is_dpt_root", "leiden", "fa_1", "fa_2",
                     "base_fa_1", "base_fa_2", "pseudotime",
                     "age_guided_pseudotime", "layout_age_weight", "age_guidance_method",
                     unname(TRAJ_DONOR_FIELDS), "C_peptide_censored")
        out <- if (is.null(d)) data.frame() else d[, intersect(columns, names(d)), drop = FALSE]
        write_islet_csv(file, out, islet_csv_header(
          "Trajectory — PAGA force-directed layout",
          trajectory_inputs_text(attr(d, "features")),
          sprintf("Age guidance weight: %.2f | method: %s", attr(d, "layout_age_weight") %||% 0,
                  attr(d, "age_guidance_method") %||% "none"),
          sprintf("Donors with age: %s | missing age: %s | age range (years): %s to %s | FA1 direction fallback: %s",
                  attr(d, "n_age_donors") %||% NA, attr(d, "n_missing_age_donors") %||% NA,
                  attr(d, "age_min_years") %||% NA, attr(d, "age_max_years") %||% NA,
                  attr(d, "age_direction_fallback") %||% NA),
          "pseudotime = baseline 56-feature DPT; age_guided_pseudotime = score shifted by donor, which may extend outside 0-1",
          sprintf("Neighbors: cosine, %d NN | Leiden resolution %g | layout %s",
                  15L, 1.0, attr(d, "layout_key") %||% "(n/a)"),
          paste("DPT root:", attr(d, "root_islet_id"), "|", attr(d, "root_method")),
          sprintf("Fitted cohort: %s, any object class", cohort_label()),
          paste("Color by — left:", input$paga_color_by %||% "pseudotime",
                "| right:", input$paga_color_right %||% "donor_status"),
          "Both panels share these fitted coordinates; coloring does not refit or subset the cohort",
          "fa_1 / fa_2 are unitless layout coordinates; only relative position is meaningful",
          paste("Rows:", nrow(out), "islets")
        ))
      }
    )

  })
}
