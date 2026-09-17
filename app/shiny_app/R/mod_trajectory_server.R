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
  path <- data_path("trajectory", name)
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(case_id = "character"))
}
trajectory_filter <- function(df, donor = "all", cell_range = c(1, Inf)) {
  if (is.null(df)) return(NULL)
  keep <- if (identical(donor, "all")) rep(TRUE, nrow(df)) else as.character(df$case_id) == donor
  if ("measured_cells" %in% names(df)) keep <- keep & df$measured_cells >= cell_range[1] & df$measured_cells <= cell_range[2]
  df[which(keep), , drop = FALSE]
}
trajectory_marker_panel <- function(features, panel = NULL) {
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
trajectory_ratios <- function(coords, counts, type, region = "both", denominator = NULL) {
  if (is.null(coords)) return(NULL)
  out <- coords
  out$numerator <- out$denominator <- out$raw_ratio <- out$plotted_value <- rep(NA_real_, nrow(out))
  out$phenotype <- rep(type, nrow(out)); out$numerator_compartment <- rep(region, nrow(out))
  out$numerator_validated <- out$denominator_validated <- rep(FALSE, nrow(out))
  out$availability_reason <- rep("count_artifacts_unavailable", nrow(out))
  if (is.null(counts)) return(out)
  if (type %in% PHENOTYPE_RETIRED) stop("Retired phenotype cannot be selected")
  # The builder records a union from unique cell IDs, not a sum of overlapping labels.
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
  out$availability_reason <- ifelse(!is.finite(out$denominator), "endocrine_denominator_unavailable",
    ifelse(out$denominator == 0, "zero_core_endocrine_denominator",
      ifelse(!is.finite(out$numerator), "selected_cell_type_count_unavailable", "available")))
  if (!is.null(denominator) && "reason" %in% names(denominator)) {
    reason <- denominator$reason[match(out$islet_id, denominator$islet_id)]
    use <- !is.finite(out$denominator) & !is.na(reason) & nzchar(reason)
    out$availability_reason[use] <- paste0("endocrine_denominator: ", reason[use])
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
      "The selected islets have zero core endocrine cells, so their ratios are undefined. Their graph coordinates and marker measurements remain available below."
    else if (any(!is.finite(d$denominator)))
      "Validated core endocrine counts are unavailable for this selection. The ratio download records the missing gate evidence for each islet. The graph and marker heatmap remain available below."
    else "Validated counts for the selected cell type are unavailable for this selection. Choose another cell type or donor; the ratio download records the reason."))
  NULL
}
trajectory_heatmap_notice <- function(d, xcol = "pseudotime_percentile") {
  if (is.null(d)) return(list(title = "Marker heatmap unavailable", message = "Offline marker measurements are not yet available for this release."))
  if (!length(grep("^core__", names(d)))) return(list(title = "Choose heatmap markers", message = "Select at least one marker to display its measurements."))
  if (!nrow(d) || !any(is.finite(d[[xcol]]))) return(list(title = "No supported islets for the heatmap", message = "Choose another donor or widen the measured core-cell range."))
  NULL
}
trajectory_notice_ui <- function(notice) {
  div(class = "alert alert-info", role = "status", strong(notice$title), p(notice$message))
}
trajectory_server <- function(id, prepared, selected_islet, forced_image, active_tab = reactive("Trajectory"),
                              donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                              remove_outliers = reactive(FALSE), outlier_threshold = reactive(3.0)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Read once per session. No clinical tuning, fit functions or legacy fallback.
    coords <- trajectory_read("coordinates.csv"); counts <- trajectory_read("counts.csv")
    features <- trajectory_read("features.csv"); edges <- trajectory_read("edges.csv")
    denominator <- trajectory_read("endocrine_counts.csv")
    if (is.null(denominator) && !identical(islet_release_schema(), "3")) denominator <- features
    panel_path <- data_path("trajectory", "marker_panel.json")
    panel <- if (file.exists(panel_path)) jsonlite::read_json(panel_path) else NULL
    markers <- trajectory_marker_panel(features, panel)
    meta_path <- data_path("donor_metadata.csv")
    meta <- if (file.exists(meta_path)) read.csv(meta_path, colClasses = c(case_id = "character")) else NULL
    observe({
      if (!is.null(meta) && "donor_status" %in% names(meta)) updateSelectInput(session, "traj_donor",
        choices = c("All donors" = "all", setNames(meta$case_id, paste0(meta$case_id, " (", meta$donor_status, ")"))))
      if (!is.null(coords)) updateSliderInput(session, "cell_range", max = max(coords$measured_cells, na.rm = TRUE), value = c(10, max(coords$measured_cells, na.rm = TRUE)))
    })
    output$traj_marker_controls <- renderUI({
      tagList(
        tags$style(".traj-marker-panel .shiny-options-group { display: grid; grid-template-columns: repeat(auto-fit, minmax(110px, 1fr)); gap: 2px 12px; } .traj-marker-panel .checkbox { margin-top: 0; margin-bottom: 4px; }"),
        div(class = "traj-marker-panel",
          checkboxGroupInput(ns("multi_heatmap_markers"), sprintf("Heatmap markers (%d)", length(markers)), choices = markers, selected = markers, width = "100%")),
        actionButton(ns("traj_markers_all"), "Select all markers"),
        actionButton(ns("traj_markers_clear"), "Clear markers")
      )
    })
    observeEvent(input$traj_markers_all, updateCheckboxGroupInput(session, "multi_heatmap_markers", selected = markers))
    observeEvent(input$traj_markers_clear, updateCheckboxGroupInput(session, "multi_heatmap_markers", selected = character(0)))
    selected <- reactive(trajectory_filter(coords, input$traj_donor %||% "all", input$cell_range %||% c(1, Inf)))
    selected_info <- reactive({
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      donor <- input$traj_donor %||% "all"
      if (!identical(donor, "all") && !identical(as.character(info$case_id), donor)) return(NULL)
      info
    })
    observeEvent(selected_islet(), {
      info <- selected_islet(); req(info)
      donor <- input$traj_donor %||% "all"
      if (!identical(donor, "all") && !identical(as.character(info$case_id), donor))
        updateSelectInput(session, "traj_donor", selected = as.character(info$case_id))
    }, ignoreInit = TRUE)
    ratio_data <- reactive(trajectory_ratios(selected(), counts, input$traj_ratio_type %||% "Beta cell", input$traj_ratio_region %||% "both", denominator))
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
      if (!is.null(review) && identical(review$status, "islet_identity_review_required")) return(
        div(class = "alert alert-warning", role = "status",
          strong("Trajectory status: islet identity review required"), p(review$message),
          p("The saved ordering is shown for inspection. Confirmed islet corrections require a new offline fit and validation."),
          if (!is.null(status$median_lodo_absolute_spearman)) p(sprintf("Previous numerical checks: median leave-one-donor-out |Spearman| %.3f; supported projection coverage %.1f%%. These checks do not establish islet identity.",
            status$median_lodo_absolute_spearman, 100 * status$supported_primary_fraction)))
      )
      div(class = "alert alert-info", role = "status", strong(paste("Trajectory validation:", status$status)),
          p(status$interpretation %||% "Validation in progress."),
          if (!is.null(status$median_lodo_absolute_spearman)) p(sprintf("Median leave-one-donor-out |Spearman|: %.3f. Supported primary coverage: %.1f%%.", status$median_lodo_absolute_spearman, 100 * status$supported_primary_fraction)))
    })
    output$traj_y_note <- renderUI({
      d <- ratio_data()
      if (is.null(d)) return(p("Ratios unavailable while offline artifacts are being built."))
      p(sprintf("log1p(selected cell-type count / core alpha–beta–delta union count). %d selected islets; %d zero endocrine denominators; %d unavailable counts or unvalidated denominators; %d without supported primary coordinates.",
        nrow(d), sum(d$denominator == 0, na.rm = TRUE), sum(!is.finite(d$numerator) | !is.finite(d$denominator)), sum(!is.finite(d$pseudotime_raw))))
    })
    output$traj_scatter_panel <- renderUI({
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
      color <- input$traj_color_by %||% "donor_status"
      g <- ggplot(d, aes(x = pt, y = plotted_value, key = key, text = islet_id, color = .data[[color]])) +
        geom_point(alpha = .5, size = 1.5) + coord_cartesian(xlim = xlimits()) + theme_minimal() +
        labs(x = if (xcol() == "pseudotime_raw") "Raw pseudotime (arbitrary direction)" else "Fixed reference-cohort percentile", y = "log1p(cells / core endocrine cell)", color = "")
      if (!identical(input$traj_show_trend, "none") && nrow(d) > 10) {
        g <- g + geom_smooth(aes(group = if (identical(input$traj_show_trend, "by_donor")) donor_status else 1), method = "loess", se = FALSE, formula = y ~ x)
      }
      ggplotly(g, tooltip = "text", source = ns("traj_scatter"))
    })
    output$traj_graph_note <- renderUI({
      d <- selected()
      if (is.null(d)) return(NULL)
      if (!nrow(d) || !any(is.finite(d$DC1) & is.finite(d$DC2))) return(
        p(role = "status", "No islets in this selection have graph coordinates. The full reference graph is shown for orientation."))
      NULL
    })
    output$paga_layout <- renderPlotly({
      shiny::validate(shiny::need(!is.null(coords), "The fixed diffusion graph is not yet available."))
      d <- selected(); all <- coords[is.finite(coords$DC1) & is.finite(coords$DC2), , drop = FALSE]
      shiny::validate(shiny::need(nrow(all) > 0, "No supported primary-component graph coordinates."))
      d <- d[is.finite(d$DC1) & is.finite(d$DC2), , drop = FALSE]
      d$key <- paste(d$case_id, d$islet_key, sep = "|")
      g <- ggplot() + geom_point(data = all, aes(DC1, DC2), color = "grey80", alpha = .3, size = 1)
      if (!is.null(edges)) {
        a <- match(edges$from, coords$islet_id); b <- match(edges$to, coords$islet_id)
        ed <- data.frame(x = coords$DC1[a], y = coords$DC2[a], xend = coords$DC1[b], yend = coords$DC2[b])
        g <- g + geom_segment(data = ed, aes(x, y, xend = xend, yend = yend), color = "grey80", linewidth = .15, alpha = .2, na.rm = TRUE)
      }
      if (nrow(d)) g <- g + geom_point(data = d, aes(DC1, DC2, color = donor_status, key = key, text = islet_id), size = 1.5) +
        scale_color_manual(values = donor_colors_reactive())
      g <- g + theme_minimal() + labs(color = "Disease status")
      ggplotly(g, tooltip = "text", source = ns("paga_layout"))
    })
    heatmap_data <- reactive({
      d <- selected(); if (is.null(d) || is.null(features)) return(NULL)
      marks <- trajectory_marker_columns(markers, input$multi_heatmap_markers)
      merge(d, features[, c("islet_id", marks), drop = FALSE], by = "islet_id", sort = FALSE)
    })
    output$traj_multi_heatmap_panel <- renderUI({
      notice <- trajectory_heatmap_notice(heatmap_data(), xcol())
      if (!is.null(notice)) return(trajectory_notice_ui(notice))
      plotOutput(ns("traj_multi_heatmap"), height = max(300, 140 + 24 * length(grep("^core__", names(heatmap_data())))))
    })
    output$traj_multi_heatmap <- renderPlot({
      d <- heatmap_data(); req(is.null(trajectory_heatmap_notice(d, xcol())))
      marks <- grep("^core__", names(d), value = TRUE)
      d$pt <- d[[xcol()]]; d <- d[is.finite(d$pt), , drop = FALSE]
      breaks <- seq(xlimits()[1], xlimits()[2], length.out = 21)
      d$bin <- cut(d$pt, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      long <- tidyr::pivot_longer(d, cols = all_of(marks), names_to = "marker", values_to = "value")
      agg <- long %>% group_by(bin, marker) %>% summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
      agg$x <- (breaks[agg$bin] + breaks[agg$bin + 1])/2
      agg$marker <- factor(agg$marker, levels = rev(marks), labels = sub("^core__", "", rev(marks)))
      ggplot(agg, aes(x, marker, fill = value)) + geom_tile(width = diff(breaks)[1]) +
        scale_fill_viridis_c(name = "Mean islet median\nlog1p core intensity") + coord_cartesian(xlim = xlimits()) + theme_minimal() + labs(x = "Fixed pseudotime", y = NULL)
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
    output$segmentation_viewer_panel <- renderUI({
      if (!identical(active_tab(), "Trajectory")) return(NULL)
      info <- selected_info()
      if (is.null(info)) return(p("Select an islet in the graph to inspect its cells or image."))
      markers <- if (is.null(features)) c("INS", "GCG", "SST") else sub("^core__", "", grep("^core__", names(features), value = TRUE))
      tagList(h4(paste(info$case_id, info$islet_key)),
        selectInput("drilldown_color_by", "Cell marker (raw intensity)", choices = markers, selected = "INS"),
        p("Cell-level marker measurements are available. Subtype labels remain unavailable pending gate validation."),
        plotOutput("islet_segmentation_view", height = 500),
        downloadButton(ns("dl_islet_cells"), "Download selected islet cell measurements"),
        actionButton(ns("open_image"), "Open donor image in Viewer"))
    })
    outputOptions(output, "segmentation_viewer_panel", suspendWhenHidden = FALSE)
    output$dl_islet_cells <- downloadHandler(filename = function() "islet_cell_measurements.csv", content = function(file) {
      info <- selected_info(); req(info)
      cells <- load_islet_cells(info$case_id, info$islet_key); req(cells)
      write.csv(cells, file, row.names = FALSE, na = "")
    })
    observeEvent(input$open_image, {
      info <- selected_info(); req(info)
      path <- find_image_for_case(sprintf("%04d", as.integer(info$case_id)), www_local_images_dir, local_images_root)
      if (!is.null(path)) forced_image(basename(path))
      updateTabsetPanel(session$rootScope(), "tabs", selected = "Viewer")
    })
  })
}
