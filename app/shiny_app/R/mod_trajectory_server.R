# All outputs share immutable coordinates and a single donor display filter.
traj_pt_extent <- function(pt, jitter_frac = 0.012) {
  dom <- suppressWarnings(range(pt, na.rm = TRUE))
  if (!all(is.finite(dom)) || diff(dom) <= 0) dom <- c(0, 1)
  j <- diff(dom) * jitter_frac
  list(dom = dom, jitter = j, lim = c(dom[1] - j, dom[2] + j))
}
apply_pt_scale <- function(df, pt_scale) {
  df$pt_raw <- df$pt
  if (!identical(pt_scale, "raw")) df$pt <- df$pt_rank
  df
}
trajectory_read <- function(name) {
  csv_path <- data_path("trajectory", name)
  parquet_path <- sub("[.]csv$", ".parquet", csv_path)
  use_parquet <- file.exists(parquet_path) && exists("islet_read_parquet", mode = "function")
  path <- if (use_parquet) parquet_path else csv_path
  if (!file.exists(path)) return(NULL)
  # One immutable copy per worker, shared across sessions. Include release path
  # and file metadata so a different candidate cannot borrow a cached table.
  key <- normalizePath(path, mustWork = TRUE)
  info <- file.info(path)
  signature <- paste(info$size, as.numeric(info$mtime), as.numeric(info$ctime))
  cached <- .trajectory_tables[[key]]
  if (!is.null(cached) && identical(cached$signature, signature)) return(cached$value)
  value <- if (use_parquet) islet_read_parquet(path) else {
    columns <- names(utils::read.csv(path, nrows = 0, check.names = FALSE))
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = if ("case_id" %in% columns) c(case_id = "character") else NA)
  }
  if ("case_id" %in% names(value)) value$case_id <- as.character(value$case_id)
  for (column in intersect(c("eligible", "supported", "landmark", "validated", "endocrine_validated", "available", "positive_available", "geometry_valid"), names(value))) {
    flag <- tolower(as.character(value[[column]]))
    value[[column]] <- ifelse(is.na(flag), NA, flag %in% c("true", "t", "1"))
  }
  if (!"validated" %in% names(value) && "available" %in% names(value)) value$validated <- value$available
  if (!"validated" %in% names(value) && "endocrine_validated" %in% names(value)) value$validated <- value$endocrine_validated
  .trajectory_tables[[key]] <- list(signature = signature, value = value)
  value
}
.trajectory_tables <- new.env(parent = emptyenv())
trajectory_filter <- function(df, donor = "all", cell_range = c(1, Inf)) {
  if (is.null(df)) return(NULL)
  keep <- if (identical(donor, "all")) rep(TRUE, nrow(df)) else as.character(df$case_id) == donor
  if ("measured_cells" %in% names(df)) keep <- keep & df$measured_cells >= cell_range[1] & df$measured_cells <= cell_range[2]
  df[which(keep), , drop = FALSE]
}
trajectory_marker_panel <- function(features, panel = NULL) {
  if (is.null(features) && !is.null(panel)) {
    markers <- unlist(panel$markers, use.names = FALSE)
    if (length(markers) != 31L || anyDuplicated(markers)) stop("The trajectory marker panel must contain all 31 unique markers")
    return(markers)
  }
  available <- sub("^core__", "", grep("^core__", names(features), value = TRUE))
  if (is.null(panel)) return(available)
  markers <- unlist(panel$markers, use.names = FALSE)
  if (length(markers) != 31L || anyDuplicated(markers) || !setequal(markers, available))
    stop("The trajectory marker panel does not match its saved measurements")
  markers
}
trajectory_marker_columns <- function(markers, selected) {
  keep <- markers[markers %in% selected]
  if (length(keep)) paste0("core__", keep) else character(0)
}
trajectory_islet_choices <- function(d) {
  if (is.null(d) || !nrow(d)) return(character(0))
  setNames(paste(d$case_id, d$islet_key, sep = "|"),
           paste0(d$islet_id, " (", d$measured_cells, " measured core cells)"))
}
trajectory_graph_coordinates <- function(coords, layout = NULL) {
  if (is.null(coords)) return(NULL)
  if (is.null(layout)) {
    coords$graph_x <- coords$DC1; coords$graph_y <- coords$DC2
    coords$layout_kind <- "primary diffusion coordinates"
    return(coords)
  }
  if (anyDuplicated(layout$islet_id) || !setequal(coords$islet_id, layout$islet_id)) stop("Graph layout differs from trajectory islet revision")
  ix <- match(coords$islet_id, layout$islet_id)
  for (column in intersect(c("graph_x", "graph_y", "graph_landmark", "layout_kind"), names(layout))) coords[[column]] <- layout[[column]][ix]
  coords
}
trajectory_ratios <- function(coords, counts, type, region = "both", denominator = NULL) {
  if (is.null(coords)) return(NULL)
  out <- coords
  out$numerator <- out$denominator <- out$raw_ratio <- out$plotted_value <- rep(NA_real_, nrow(out))
  out$phenotype <- rep(type, nrow(out)); out$numerator_compartment <- rep(region, nrow(out))
  out$numerator_validated <- out$denominator_validated <- rep(FALSE, nrow(out))
  out$availability_reason <- rep("count_artifacts_unavailable", nrow(out))
  if (is.null(counts)) return(out)
  catalog <- if (exists("islet_typing_catalog", mode = "function")) islet_typing_catalog() else NULL
  if (is.null(catalog) && type %in% PHENOTYPE_RETIRED) stop("Retired phenotype cannot be selected")
  if (!is.null(catalog) && !type %in% islet_catalog_names(c(catalog$base_types, catalog$states))) stop("Cell type absent from release catalog")
  # The builder records unique canonical cells inside the endocrine
  # pixel-classifier core, independently of downstream phenotype labels.
  if (!is.null(denominator)) {
    ix <- match(out$islet_id, denominator$islet_id)
    out$denominator <- denominator$endocrine_count[ix]
    out$denominator_validated <- if ("validated" %in% names(denominator)) denominator$validated[ix] %in% TRUE else is.finite(out$denominator)
    out$denominator[!out$denominator_validated] <- NA_real_
  }
  cts <- counts[counts$phenotype == type, , drop = FALSE]
  core <- cts[cts$compartment == "core", , drop = FALSE]
  peri <- cts[cts$compartment == "peri", , drop = FALSE]
  a <- core$count[match(out$islet_id, core$islet_id)]
  b <- peri$count[match(out$islet_id, peri$islet_id)]
  a_ok <- if ("validated" %in% names(core)) core$validated[match(out$islet_id, core$islet_id)] %in% TRUE else is.finite(a)
  b_ok <- if ("validated" %in% names(peri)) peri$validated[match(out$islet_id, peri$islet_id)] %in% TRUE else is.finite(b)
  a[!a_ok] <- NA_real_; b[!b_ok] <- NA_real_
  out$numerator <- switch(region, core = a, peri = b, both = a + b, stop("Invalid numerator compartment"))
  out$numerator_validated <- switch(region, core = a_ok, peri = b_ok, both = a_ok & b_ok)
  valid <- is.finite(out$denominator) & out$denominator > 0 & is.finite(out$numerator)
  out$raw_ratio[valid] <- out$numerator[valid] / out$denominator[valid]
  out$plotted_value[valid] <- log1p(out$raw_ratio[valid])
  out$availability_reason <- ifelse(!is.finite(out$denominator), "classifier_core_cell_census_unavailable",
    ifelse(out$denominator == 0, "no_core_cells_in_classifier_area",
      ifelse(!is.finite(out$numerator), "selected_cell_type_count_unavailable", "available")))
  if (!is.null(denominator) && "reason" %in% names(denominator)) {
    reason <- denominator$reason[match(out$islet_id, denominator$islet_id)]
    use <- !is.finite(out$denominator) & !is.na(reason) & nzchar(reason)
    out$availability_reason[use] <- paste0("classifier_core_census: ", reason[use])
  }
  if ("reason" %in% names(counts)) {
    reason <- switch(region, core = core$reason[match(out$islet_id, core$islet_id)],
      peri = peri$reason[match(out$islet_id, peri$islet_id)],
      both = ifelse(!is.finite(a), core$reason[match(out$islet_id, core$islet_id)], peri$reason[match(out$islet_id, peri$islet_id)]))
    use <- is.finite(out$denominator) & out$denominator > 0 & !is.finite(out$numerator) & !is.na(reason) & nzchar(reason)
    out$availability_reason[use] <- paste0("selected_cell_type: ", reason[use])
  }
  out
}
trajectory_ratio_notice <- function(d, xcol = "pseudotime_percentile") {
  if (is.null(d)) return(list(title = "Trajectory unavailable", message = "Offline trajectory artifacts are not yet available for this release."))
  if (!nrow(d)) return(list(title = "No islets in this selection", message = "Choose another donor or widen the measured core-cell range."))
  if (!any(is.finite(d[[xcol]]))) return(list(title = "No supported coordinates in this selection", message = "These islets have no comparable primary pseudotime coordinates. Choose another donor or widen the measured core-cell range."))
  if (!any(is.finite(d[[xcol]]) & is.finite(d$plotted_value))) return(list(
    title = "Cell-type ratios unavailable",
    message = if (all(is.finite(d$denominator)) && all(d$denominator == 0))
      "The selected classifier regions contain no core cells, so their ratios are undefined and they are not eligible islets."
    else if (any(!is.finite(d$denominator)))
      "The classifier-core cell census is unavailable for this selection. The ratio download records the affected islets."
    else "Validated counts for the selected cell type are unavailable for this selection. Choose another cell type or donor; the ratio download records the reason."))
  NULL
}
trajectory_heatmap_notice <- function(d, xcol = "pseudotime_percentile") {
  if (is.null(d)) return(list(title = "Marker heatmap unavailable", message = "Offline marker measurements are not yet available for this release."))
  if (!length(grep("^core__", names(d))) && !("marker" %in% names(d) && length(unique(d$marker)))) return(list(title = "Choose heatmap markers", message = "Select at least one marker to display its measurements."))
  if (!nrow(d) || !any(is.finite(d[[xcol]]))) return(list(title = "No supported islets for the heatmap", message = "Choose another donor or widen the measured core-cell range."))
  NULL
}
trajectory_notice_ui <- function(notice) {
  div(class = "alert alert-info", role = "status", strong(notice$title), p(notice$message))
}
trajectory_axis_label <- function(scale = "rank") {
  if (identical(scale, "raw")) "Raw pseudotime (arbitrary direction)" else
    "Fixed reference-cohort rank (exact ties hash-spread)"
}
trajectory_zero_ratio_message <- function(d, xcol = "pseudotime_percentile", type = "cell type", region = "core") {
  if (is.null(d) || !nrow(d) || !all(c(xcol, "plotted_value") %in% names(d))) return(NULL)
  plotted <- is.finite(d[[xcol]]) & is.finite(d$plotted_value)
  if (!any(plotted) || !all(d$plotted_value[plotted] == 0)) return(NULL)
  sprintf("All %s plotted ratios are valid measured zeros for %s in %s. The points therefore lie on y = 0; choose another cell type or compartment to inspect variation.",
    format(sum(plotted), big.mark = ","), type, switch(region, core = "the core", peri = "peri-islet tissue", both = "the combined core and peri-islet compartments", region))
}
trajectory_heatmap_axis <- function(coords, status, scale = "rank") {
  if (identical(status$coordinates_available, FALSE)) {
    limits <- traj_pt_extent(log10(coords$area_um2[is.finite(coords$area_um2) & coords$area_um2 > 0]), 0)$dom
    return(list(column = "area_um2", transform = log10, limits = limits,
                label = "log10 section area (µm²); pseudotime unavailable", kind = "section_area"))
  }
  column <- if (identical(scale, "raw")) "pseudotime_raw" else "pseudotime_percentile"
  list(column = column, transform = identity,
       limits = if (identical(scale, "raw")) traj_pt_extent(coords$pseudotime_raw, 0)$dom else c(0, 1),
       label = trajectory_axis_label(scale), kind = "pseudotime")
}
trajectory_server <- function(id, prepared, selected_islet, forced_image, active_tab = reactive("Trajectory"),
                              donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                              remove_outliers = reactive(FALSE), outlier_threshold = reactive(3.0)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    typed_release <- identical(islet_release_schema(), "4")
    # Read once per session. No clinical tuning, fit functions or legacy fallback.
    coords <- trajectory_read("coordinates.csv"); counts <- trajectory_read("counts.csv")
    coords <- trajectory_graph_coordinates(coords, if (typed_release) trajectory_read("graph_layout.csv") else NULL)
    features <- if (typed_release) NULL else trajectory_read("features.csv")
    edges <- trajectory_read("edges.csv")
    denominator <- trajectory_read("endocrine_counts.csv")
    status_path <- data_path("trajectory", "validation.json")
    status <- if (file.exists(status_path)) jsonlite::read_json(status_path) else list(status = "unavailable")
    axis_available <- !identical(status$coordinates_available, FALSE)
    if (is.null(denominator) && !islet_release_schema() %in% c("3", "4")) denominator <- features
    panel_path <- data_path("trajectory", "marker_panel.json")
    panel <- if (file.exists(panel_path)) jsonlite::read_json(panel_path) else NULL
    markers <- trajectory_marker_panel(features, panel)
    meta_path <- data_path("donor_metadata.csv")
    meta <- if (file.exists(meta_path)) read.csv(meta_path, colClasses = c(case_id = "character")) else NULL
    observe({
      if (typed_release) updateSelectInput(session, "traj_ratio_type", choices = islet_phenotype_choices(), selected = "Beta cell")
      if (!is.null(meta) && "donor_status" %in% names(meta)) updateSelectInput(session, "traj_donor",
        choices = c("All donors" = "all", setNames(meta$case_id, paste0(meta$case_id, " (", meta$donor_status, ")"))))
      if (!is.null(coords)) updateSliderInput(session, "cell_range", max = max(coords$measured_cells, na.rm = TRUE), value = c(0, max(coords$measured_cells, na.rm = TRUE)))
    })
    output$traj_marker_controls <- renderUI({
      tagList(
        if (typed_release) fluidRow(
          column(4, selectInput(ns("heatmap_measure"), "Heatmap measurement", choices = c("Marker expression within measured cells" = "typed_expression", "Marker-positive cells per classifier-core cell" = "positive_ratio"), selected = "typed_expression")),
          column(4, selectInput(ns("heatmap_region"), "Heatmap cells in", choices = c("Core" = "core", "Core + peri-islet" = "both", "Peri-islet" = "peri"), selected = "core")),
          column(4, conditionalPanel(sprintf("input['%s'] === 'typed_expression'", ns("heatmap_measure")),
            selectInput(ns("heatmap_population"), "Typed population", choices = c(list("Endocrine union" = islet_endocrine_population()), islet_phenotype_choices()), selected = islet_endocrine_population())))
        ),
        tags$style(".traj-marker-panel .shiny-options-group { display: grid; grid-template-columns: repeat(auto-fit, minmax(110px, 1fr)); gap: 2px 12px; } .traj-marker-panel .checkbox { margin-top: 0; margin-bottom: 4px; }"),
        div(class = "traj-marker-panel",
          checkboxGroupInput(ns("multi_heatmap_markers"), sprintf("Heatmap markers (%d)", length(markers)), choices = markers, selected = markers, width = "100%")),
        sliderInput(ns("multi_heatmap_nbins"), "Number of pseudotime bins", min = 10, max = 40, value = 20, step = 5),
        actionButton(ns("traj_markers_all"), "Select all markers"),
        actionButton(ns("traj_markers_clear"), "Clear markers")
      )
    })
    observeEvent(input$traj_markers_all, updateCheckboxGroupInput(session, "multi_heatmap_markers", selected = markers))
    observeEvent(input$traj_markers_clear, updateCheckboxGroupInput(session, "multi_heatmap_markers", selected = character(0)))
    selected <- reactive(trajectory_filter(coords, input$traj_donor %||% "all", input$cell_range %||% c(0, Inf)))
    output$boundary_review <- renderUI(islet_boundary_notice(input$traj_donor %||% "all"))
    observe({
      d <- selected(); req(d)
      choices <- trajectory_islet_choices(d)
      updateSelectizeInput(session, "inspect_islet_id", choices = choices, selected = character(0), server = TRUE)
    })
    observeEvent(input$inspect_islet, {
      key <- input$inspect_islet_id; req(key)
      d <- selected(); req(key %in% paste(d$case_id, d$islet_key, sep = "|"))
      parts <- strsplit(key, "|", fixed = TRUE)[[1]]
      if (length(parts) == 2 && all(vapply(parts, islet_safe_path_segment, logical(1))))
        select_islet_from_key(parts[1], parts[2], prepared, selected_islet)
    })
    selected_info <- reactive({
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      donor <- input$traj_donor %||% "all"
      if (!identical(donor, "all") && !identical(as.character(info$case_id), donor)) return(NULL)
      if (!paste(info$case_id, info$islet_key, sep = "_") %in% selected()$islet_id) return(NULL)
      info
    })
    observeEvent(selected_islet(), {
      info <- selected_islet(); req(info)
      donor <- input$traj_donor %||% "all"
      if (!identical(donor, "all") && !identical(as.character(info$case_id), donor))
        updateSelectInput(session, "traj_donor", selected = as.character(info$case_id))
    }, ignoreInit = TRUE)
    ratio_data <- reactive(trajectory_ratios(selected(), counts, input$traj_ratio_type %||% "Beta cell", input$traj_ratio_region %||% "core", denominator))
    xcol <- reactive(if (identical(input$pt_scale, "raw")) "pseudotime_raw" else "pseudotime_percentile")
    xlimits <- reactive({
      if (is.null(coords)) return(c(0,1))
      if (identical(input$pt_scale, "raw")) traj_pt_extent(coords$pseudotime_raw, 0)$dom else c(0,1)
    })
    output$traj_status <- renderUI({
      path <- data_path("trajectory", "validation.json")
      status <- if (file.exists(path)) jsonlite::read_json(path) else list(status = "unavailable", interpretation = "No outcome-blind trajectory artifacts are available for this release.")
      review_path <- data_path("trajectory", "quality_review.json")
      review <- if (file.exists(review_path)) jsonlite::read_json(review_path) else NULL
      if (!axis_available) return(trajectory_notice_ui(list(title = if (typed_release) "Pseudotime unavailable for this release" else "Pseudotime awaiting a reviewed refit",
        message = paste(review$message %||% status$interpretation,
          if (typed_release) "Browse current outlines and all 31 markers below. Cell-type ratios use the classifier-core cell census as denominator and retain missing phenotype numerators." else
            "Browse current outlines and all 31 marker measurements below. Cell-type ratios separately require validated phenotype gates."))))
      div(class = "alert alert-info", role = "status", strong(paste("Trajectory validation:", status$status)),
          p(status$interpretation %||% "Validation in progress."),
          p("An islet is eligible when its endocrine pixel-classifier core does not contact the image edge and contains at least one canonical core cell. Downstream phenotype support can limit a trajectory coordinate, but cannot revoke eligibility."),
          if (!is.null(status$eligible)) p(sprintf("%s classifier-eligible islets; %s have complete supported typed-composition features; %s have supported primary coordinates.",
            format(status$eligible, big.mark = ","), format(status$feature_supported %||% NA_integer_, big.mark = ","), format(status$supported %||% NA_integer_, big.mark = ","))),
          if (!is.null(status$median_lodo_absolute_spearman)) p(sprintf("Median leave-one-donor-out |Spearman|: %.3f. Supported primary coverage: %.1f%%.", status$median_lodo_absolute_spearman, 100 * status$supported_primary_fraction)))
    })
    output$traj_y_note <- renderUI({
      d <- ratio_data()
      if (is.null(d)) return(p("Ratios unavailable while offline artifacts are being built."))
      eligible <- if ("eligible" %in% names(d)) d$eligible %in% TRUE else rep(TRUE, nrow(d))
      primary <- is.finite(d[[xcol()]])
      ratio_available <- is.finite(d$plotted_value)
      summary <- sprintf(paste0("log1p(selected cell-type count / ", if (typed_release) "canonical cell count inside the endocrine pixel-classifier core" else "core alpha–beta–delta union count", "). ",
        "%d selected islets; %d plotted points (one per islet with a primary coordinate and defined ratio); ",
        "%d excluded by classifier eligibility; %d eligible islets without a supported primary coordinate; ",
        "%d primary-coordinate ratios unavailable; %d classifier cores with no cells."),
        nrow(d), sum(primary & ratio_available), sum(!eligible), sum(eligible & !primary),
        sum(primary & !ratio_available), sum(d$denominator == 0, na.rm = TRUE))
      zero_message <- trajectory_zero_ratio_message(d, xcol(), input$traj_ratio_type %||% "Beta cell", input$traj_ratio_region %||% "core")
      tagList(p(summary), if (!is.null(zero_message)) p(role = "status", strong(zero_message)))
    })
    output$traj_scatter_panel <- renderUI({
      if (!axis_available) return(trajectory_notice_ui(list(title = "Trajectory plot unavailable",
        message = if (typed_release) "This release has no supported saved ordering. The validation status above records the reason; marker measurements remain available below." else
          "The previous fit has been withdrawn following the boundary correction. A new ordering will be shown after the remaining boundary reviews and offline validation.")))
      notice <- trajectory_ratio_notice(ratio_data(), xcol())
      if (!is.null(notice)) return(trajectory_notice_ui(notice))
      plotlyOutput(ns("traj_scatter"), height = 500)
    })
    output$traj_scatter <- renderPlotly({
      d <- ratio_data()
      req(is.null(trajectory_ratio_notice(d, xcol())))
      d$pt <- d[[xcol()]]
      d <- d[is.finite(d$pt) & is.finite(d$plotted_value), , drop = FALSE]
      d$key <- paste(d$case_id, d$islet_key, sep = "|")
      d$hover <- sprintf("%s<br>Donor: %s (%s)<br>Reference rank: %.4f<br>Raw pseudotime: %.6f<br>Numerator: %g<br>Core endocrine denominator: %g<br>log1p ratio: %.4f",
        d$islet_id, d$case_id, d$donor_status, d$pseudotime_percentile, d$pseudotime_raw,
        d$numerator, d$denominator, d$plotted_value)
      color <- input$traj_color_by %||% "donor_status"
      palette <- get_donor_color_palette(input$traj_donor_palette %||% "Bright")
      size_by <- input$traj_point_size %||% "uniform"
      point_size <- input$traj_point_size_slider %||% 3
      d$.point_size <- point_size
      if (identical(size_by, "measured_cells") && "measured_cells" %in% names(d)) d$.point_size <- scales::rescale(d$measured_cells, to = c(point_size * .55, point_size * 1.8))
      if (identical(size_by, "area_um2") && "area_um2" %in% names(d)) d$.point_size <- scales::rescale(log1p(pmax(d$area_um2, 0)), to = c(point_size * .55, point_size * 1.8))
      ymax <- max(d$plotted_value, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) ymax <- .01
      g <- ggplot(d, aes(x = pt, y = plotted_value, key = key, text = hover, color = .data[[color]], size = .point_size)) +
        geom_point(alpha = input$traj_alpha %||% .6, show.legend = color != "case_id", stroke = 0) + scale_size_identity() +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
        coord_cartesian(xlim = xlimits(), ylim = c(0, ymax * 1.05)) + theme_islet() +
        labs(x = trajectory_axis_label(input$pt_scale %||% "rank"), y = "log1p(cells / classifier-core cell)", color = "")
      if (identical(color, "donor_status")) g <- g + scale_color_manual(values = palette, drop = FALSE)
      if (identical(color, "case_id")) {
        ids <- sort(unique(as.character(d$case_id))); idpal <- setNames(grDevices::hcl.colors(length(ids), palette = "Dark 3"), ids)
        g <- g + scale_color_manual(values = idpal, drop = FALSE)
      }
      if (!identical(input$traj_show_trend, "none") && nrow(d) > 10) {
        trend_formula <- y ~ stats::poly(x, degree = 3, raw = FALSE)
        if (identical(input$traj_show_trend, "by_donor")) {
          g <- g + geom_smooth(data = d, aes(x = pt, y = plotted_value, color = donor_status, group = donor_status),
            inherit.aes = FALSE, method = "lm", se = isTRUE(input$traj_trend_ci), formula = trend_formula)
        } else {
          g <- g + geom_smooth(data = d, aes(x = pt, y = plotted_value), inherit.aes = FALSE,
            method = "lm", se = isTRUE(input$traj_trend_ci), formula = trend_formula, color = "#333333")
        }
      }
      plot <- ggplotly(g, tooltip = "text", source = ns("traj_scatter")) %>%
        plotly_islet_config("trajectory_ratio_scatter")
      plot <- event_register(plot, "plotly_click")
      if (exists("islet_click_bridge", mode = "function")) plot <- islet_click_bridge(plot, ns("traj_native_click"))
      plot
    })
    output$traj_strip_panel <- renderUI({
      if (!axis_available) return(NULL)
      notice <- trajectory_ratio_notice(ratio_data(), xcol())
      if (!is.null(notice)) return(NULL)
      tagList(
        tags$small("Pseudotime strip: hue shows the locally enriched donor status; saturation shows the selected per-endocrine-cell ratio.",
          style = "display:block; color:#666; margin-top:-8px;"),
        plotlyOutput(ns("traj_strip_heatmap"), height = 105)
      )
    })
    output$traj_strip_heatmap <- renderPlotly({
      d <- ratio_data(); req(is.null(trajectory_ratio_notice(d, xcol())))
      d$pt <- d[[xcol()]]
      d <- d[is.finite(d$pt) & is.finite(d$plotted_value), , drop = FALSE]
      limits <- xlimits(); breaks <- seq(limits[1], limits[2], length.out = 26L)
      d$bin <- cut(d$pt, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      d <- d[!is.na(d$bin), , drop = FALSE]
      statuses <- c("ND", "Aab+", "T1D")
      totals <- table(factor(d$donor_status, levels = statuses))
      rows <- lapply(seq_len(25L), function(bin) {
        current <- d[d$bin == bin, , drop = FALSE]
        counts <- table(factor(current$donor_status, levels = statuses))
        density <- ifelse(totals > 0, counts / totals, 0)
        dominant <- if (nrow(current)) statuses[which.max(density)] else NA_character_
        data.frame(bin = bin, x = mean(breaks[c(bin, bin + 1L)]), n = nrow(current),
          dominant = dominant, value = if (nrow(current)) mean(current$plotted_value) else NA_real_)
      })
      hm <- do.call(rbind, rows)
      finite <- hm$value[is.finite(hm$value)]
      span <- if (length(finite)) quantile(finite, c(.05, .95), na.rm = TRUE) else c(0, 1)
      scaled <- if (diff(span) > 0) pmax(0, pmin(1, (hm$value - span[1]) / diff(span))) else rep(1, nrow(hm))
      scaled[!is.finite(scaled)] <- 0
      donor_palette <- get_donor_color_palette(input$traj_donor_palette %||% "Bright")
      blend <- function(hex, amount) {
        rgb <- grDevices::col2rgb(hex) / 255
        grDevices::rgb((.35 + .65 * amount) * rgb[1] + (.65 - .65 * amount),
          (.35 + .65 * amount) * rgb[2] + (.65 - .65 * amount),
          (.35 + .65 * amount) * rgb[3] + (.65 - .65 * amount))
      }
      hm$fill <- mapply(function(status, amount) if (is.na(status)) "#D9D9D9" else blend(donor_palette[[status]], amount), hm$dominant, scaled)
      hm$hover <- sprintf("Trajectory coordinate: %.3f<br>Dominant status: %s<br>Islets: %d<br>Mean log1p ratio: %.4f",
        hm$x, ifelse(is.na(hm$dominant), "none", hm$dominant), hm$n, hm$value)
      g <- ggplot(hm, aes(x, 1, fill = fill, text = hover)) +
        geom_tile(width = diff(breaks)[1], height = .9) + scale_fill_identity() +
        coord_cartesian(xlim = limits, ylim = c(.5, 1.5)) +
        scale_x_continuous(expand = c(0, 0)) + theme_void() +
        theme(plot.margin = margin(0, 5, 0, 5))
      plotly::ggplotly(g, tooltip = "text", source = ns("traj_strip_heatmap")) %>%
        plotly::layout(xaxis = list(range = limits, title = "", showticklabels = FALSE),
          yaxis = list(visible = FALSE), margin = list(l = 45, r = 10, t = 2, b = 2), showlegend = FALSE) %>%
        plotly::config(displayModeBar = FALSE, displaylogo = FALSE)
    })
    output$traj_graph_note <- renderUI({
      if (!axis_available) return(NULL)
      d <- selected()
      if (is.null(d)) return(NULL)
      if (!nrow(d) || !any(is.finite(d$graph_x) & is.finite(d$graph_y))) return(
        p(role = "status", "No islets in this selection have graph coordinates. The full reference graph is shown for orientation."))
      NULL
    })
    output$traj_graph_panel <- renderUI({
      if (!axis_available || is.null(coords) || !any(is.finite(coords$graph_x) & is.finite(coords$graph_y)))
        return(trajectory_notice_ui(list(title = "Diffusion graph unavailable", message = "No previous graph is reused with the corrected boundaries. Islets can be inspected directly below.")))
      plotlyOutput(ns("paga_layout"), height = 500)
    })
    output$paga_layout <- renderPlotly({
      shiny::validate(shiny::need(!is.null(coords), "The fixed diffusion graph is not yet available."))
      d <- selected(); all <- coords[is.finite(coords$graph_x) & is.finite(coords$graph_y), , drop = FALSE]
      shiny::validate(shiny::need(nrow(all) > 0, "No saved graph coordinates are available."))
      d <- d[is.finite(d$graph_x) & is.finite(d$graph_y), , drop = FALSE]
      d$key <- paste(d$case_id, d$islet_key, sep = "|")
      g <- ggplot() + geom_point(data = all, aes(graph_x, graph_y), color = "grey80", alpha = .3, size = 1)
      if (!is.null(edges)) {
        a <- match(edges$from, coords$islet_id); b <- match(edges$to, coords$islet_id)
        ed <- data.frame(x = coords$graph_x[a], y = coords$graph_y[a], xend = coords$graph_x[b], yend = coords$graph_y[b])
        g <- g + geom_segment(data = ed, aes(x, y, xend = xend, yend = yend), color = "grey80", linewidth = .15, alpha = .2, na.rm = TRUE)
      }
      if (nrow(d)) g <- g + geom_point(data = d, aes(graph_x, graph_y, color = donor_status, key = key, text = islet_id), size = input$traj_point_size_slider %||% 3) +
        scale_color_manual(values = get_donor_color_palette(input$traj_donor_palette %||% "Bright"), drop = FALSE)
      g <- g + theme_islet() + labs(color = "Disease status", x = if (typed_release) "Graph layout 1" else "DC1", y = if (typed_release) "Graph layout 2" else "DC2")
      plot <- ggplotly(g, tooltip = "text", source = ns("paga_layout")) %>%
        plotly_islet_config("trajectory_graph")
      plot <- event_register(plot, "plotly_click")
      if (exists("islet_click_bridge", mode = "function")) plot <- islet_click_bridge(plot, ns("graph_native_click"))
      plot
    })
    heatmap_data <- reactive({
      if (typed_release) {
        measure <- input$heatmap_measure %||% "typed_expression"
        population <- if (identical(measure, "positive_ratio")) "All cells" else input$heatmap_population %||% islet_endocrine_population()
        marks <- markers[markers %in% input$multi_heatmap_markers]
        region <- input$heatmap_region %||% "core"
        summaries <- trajectory_read_marker_summaries(population, marks, input$traj_donor %||% "all", region)
        return(trajectory_typed_marker_values(selected(), summaries, denominator, marks, measure, population, region))
      }
      d <- selected(); if (is.null(d) || is.null(features)) return(NULL)
      marks <- trajectory_marker_columns(markers, input$multi_heatmap_markers)
      merge(d, features[, c("islet_id", marks), drop = FALSE], by = "islet_id", sort = FALSE)
    })
    heatmap_axis <- reactive(trajectory_heatmap_axis(coords, status, input$pt_scale %||% "rank"))
    output$heatmap_axis_note <- renderUI({
      if (typed_release) return(p(role = "status", paste(
        if (identical(heatmap_axis()$kind, "section_area")) "Pseudotime is unavailable; fixed cohort-wide section-area bins are shown." else "Marker measurements use the saved ordering and fixed cohort-wide bins. The reference rank spreads exact raw-coordinate ties deterministically by islet ID only; order inside a tie has no biological meaning.",
        "Bins average islets within each donor, then weight donors equally. Gray cells denote unavailable measurements. Marker-positive ratios use the classifier-defined core-cell denominator; typed expression is the mean log1p intensity within the selected cell population.")))
      if (identical(heatmap_axis()$kind, "section_area")) p(role = "status",
        "Pseudotime is unavailable. This descriptive heatmap groups islets by section area using fixed cohort-wide bins. Marker values are mean islet medians of log1p core intensity; unmeasured values remain missing.")
      else p("Marker measurements use the saved ordering and fixed cohort-wide bins.")
    })
    output$traj_multi_heatmap_panel <- renderUI({
      notice <- trajectory_heatmap_notice(heatmap_data(), heatmap_axis()$column)
      if (!is.null(notice)) return(trajectory_notice_ui(notice))
      n_markers <- if (typed_release) length(unique(heatmap_data()$marker)) else length(grep("^core__", names(heatmap_data())))
      plotOutput(ns("traj_multi_heatmap"), height = max(360, 140 + 44 * n_markers))
    })
    output$traj_multi_heatmap <- renderPlot({
      axis <- heatmap_axis()
      d <- heatmap_data(); req(is.null(trajectory_heatmap_notice(d, axis$column)))
      if (typed_release) {
        marks <- markers[markers %in% input$multi_heatmap_markers]
        agg <- trajectory_heatmap_bin_values(d, axis, marks, bins = input$multi_heatmap_nbins %||% 20L)
        agg$marker <- factor(agg$marker, levels = rev(marks))
        agg$scaled_value <- ave(agg$value, agg$marker, FUN = function(values) {
          finite <- is.finite(values)
          out <- rep(NA_real_, length(values))
          if (sum(finite) < 2L) return(out)
          spread <- stats::sd(values[finite])
          if (!is.finite(spread) || spread == 0) spread <- 1
          out[finite] <- pmax(-2.5, pmin(2.5, (values[finite] - mean(values[finite])) / spread))
          out
        })
        title <- if (identical(input$heatmap_measure %||% "typed_expression", "positive_ratio"))
          "Marker-positive cells per classifier-core cell along pseudotime" else
          "Typed-population marker expression along pseudotime"
        return(ggplot(agg, aes(x, marker, fill = scaled_value)) + geom_tile(width = agg$width[1]) +
          scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
            limits = c(-2.5, 2.5), name = "Within-marker\nZ-score", na.value = "grey85", oob = scales::squish) +
          coord_cartesian(xlim = axis$limits) + theme_islet() +
          theme(axis.text.y = element_text(size = ISLET_FONT_SIZES$dense_min)) +
          labs(title = title, x = axis$label, y = NULL))
      }
      marks <- grep("^core__", names(d), value = TRUE)
      d$pt <- axis$transform(d[[axis$column]]); d <- d[is.finite(d$pt), , drop = FALSE]
      breaks <- seq(axis$limits[1], axis$limits[2], length.out = 21)
      d$bin <- cut(d$pt, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      long <- tidyr::pivot_longer(d, cols = all_of(marks), names_to = "marker", values_to = "value")
      agg <- long %>% group_by(bin, marker) %>% summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
      agg$x <- (breaks[agg$bin] + breaks[agg$bin + 1])/2
      agg$marker <- factor(agg$marker, levels = rev(marks), labels = sub("^core__", "", rev(marks)))
      ggplot(agg, aes(x, marker, fill = value)) + geom_tile(width = diff(breaks)[1]) +
        scale_fill_viridis_c(name = "Mean islet median\nlog1p core intensity", na.value = "grey85") + coord_cartesian(xlim = axis$limits) + theme_minimal() + labs(x = axis$label, y = NULL)
    })
    output$dl_trajectory_csv <- downloadHandler(filename = function() "trajectory_ratios.csv", content = function(file) write.csv(ratio_data(), file, row.names = FALSE, na = ""))
    output$dl_paga_csv <- downloadHandler(filename = function() "trajectory_graph.csv", content = function(file) write.csv(selected(), file, row.names = FALSE, na = ""))
    output$dl_multi_heatmap_csv <- downloadHandler(filename = function() "trajectory_markers.csv", content = function(file) write.csv(heatmap_data(), file, row.names = FALSE, na = ""))
    for (source in c("traj_scatter", "paga_layout")) local({
      name <- source
      observeEvent(event_data("plotly_click", source = ns(name)), {
        ev <- event_data("plotly_click", source = ns(name)); req(ev$key)
        parts <- strsplit(as.character(ev$key[1]), "|", fixed = TRUE)[[1]]
        if (length(parts) == 2 && all(vapply(parts, islet_safe_path_segment, logical(1)))) select_islet_from_key(parts[1], parts[2], prepared, selected_islet)
      }, ignoreInit = TRUE)
    })
    for (source in c("traj_native_click", "graph_native_click")) local({
      name <- source
      observeEvent(input[[name]], {
        key <- input[[name]]$key; req(key)
        # Browser events may only select an islet included by the current donor filter.
        req(key %in% paste(selected()$case_id, selected()$islet_key, sep = "|"))
        parts <- strsplit(as.character(key), "|", fixed = TRUE)[[1]]
        if (length(parts) == 2 && all(vapply(parts, islet_safe_path_segment, logical(1)))) select_islet_from_key(parts[1], parts[2], prepared, selected_islet)
      }, ignoreInit = TRUE)
    })
    output$segmentation_viewer_panel <- renderUI({
      if (!identical(active_tab(), "Trajectory")) return(NULL)
      info <- selected_info()
      if (is.null(info)) return(p("Choose an islet above and click Show outline and cells, or select a point in an available graph."))
      tagList(h4(paste(info$case_id, info$islet_key)),
        selectInput("drilldown_color_by", if (typed_release) "Cell identity or marker" else "Cell marker (raw intensity)", choices = if (typed_release) islet_cell_color_choices(markers) else markers, selected = if (typed_release) "phenotype" else "INS"),
        p(if (typed_release) "Base identity, finer phenotype states and assignment support are preserved in the cell download. Unresolved cells remain explicit." else "Cell-level marker measurements are available. Subtype labels remain unavailable pending gate validation."),
        plotOutput("islet_segmentation_view", height = 500),
        downloadButton(ns("dl_islet_cells"), "Download selected islet cell measurements"),
        if (dir.exists(data_path("islet_boundaries"))) downloadButton(ns("dl_islet_boundary"), "Download selected islet outline"),
        actionButton(ns("open_image"), "Open donor image in Viewer"))
    })
    outputOptions(output, "segmentation_viewer_panel", suspendWhenHidden = FALSE)
    output$dl_islet_cells <- downloadHandler(filename = function() "islet_cell_measurements.csv", content = function(file) {
      info <- selected_info(); req(info)
      cells <- load_islet_cells(info$case_id, info$islet_key); req(cells)
      write.csv(cells, file, row.names = FALSE, na = "")
    })
    output$dl_islet_boundary <- downloadHandler(filename = function() "islet_boundary.geojson", content = function(file) {
      info <- selected_info(); req(info)
      req(islet_safe_path_segment(as.character(info$case_id)), islet_safe_path_segment(info$islet_key))
      path <- data_path("islet_boundaries", paste0(info$case_id, "_", info$islet_key, ".geojson"))
      req(file.exists(path))
      # Prebuilt single-islet files retain exact source doubles, with no full-donor read.
      file.copy(path, file, overwrite = TRUE)
    })
    observeEvent(input$open_image, {
      info <- selected_info(); req(info)
      path <- find_image_for_case(sprintf("%04d", as.integer(info$case_id)), www_local_images_dir, local_images_root)
      if (!is.null(path)) forced_image(basename(path))
      updateTabsetPanel(session$rootScope(), "tabs", selected = "Viewer")
    })
  })
}
