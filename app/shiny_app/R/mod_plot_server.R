# ---------- Plot module server ----------
# Exports: plot_server(id, prepared, selected_islet)
#
# Dependencies (from utility files already sourced):
#   safe_left_join, add_islet_key, compute_diameter_um  -- utils_safe_join.R
#   bin_islet_sizes, summary_stats                      -- utils_stats.R / data_loading.R
#   islet_spatial_lookup, render_islet_segmentation_plot -- segmentation_helpers.R
#   SF_AVAILABLE, PIXEL_SIZE_UM                         -- global.R
#
# Parameters:
#   id              -- module namespace id
#   prepared        -- reactive returning list(targets_all, markers_all, comp)
#   selected_islet  -- reactiveVal for cross-module click-to-segmentation

plot_server <- function(id, prepared, selected_islet) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Donor palette helper (consistent colors across scatter & distribution) ----
    get_donor_palette <- function(ids) {
      if (is.null(ids)) return(character(0))
      ids_chr <- sort(unique(as.character(ids)))
      n <- length(ids_chr)
      if (n == 0) return(character(0))
      if (n <= 12) {
        cols <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Paired")
      } else {
        cols <- rainbow(n, s = 0.8, v = 0.8)
      }
      cols <- cols[seq_len(n)]
      names(cols) <- ids_chr
      cols
    }

    # ---- Pseudo-log break points (includes 0 + powers of 10) ----
    pseudo_log_breaks <- function(base = 10) {
      function(limits) {
        max_val <- max(limits)
        if (max_val <= 0) return(0)
        max_pow <- ceiling(log(max_val, base))
        min_pow <- if (max_val < 1) floor(log(max_val, base)) else 0
        brks <- c(0, base^seq(min_pow, max_pow))
        sort(unique(brks[brks <= max_val * 1.1]))
      }
    }

    # ---- Reactive storage for outlier table (replaces global <<- assignment) ----
    plot_outliers <- reactiveVal(NULL)

    # ===========================================================================
    # Dynamic UI outputs
    # ===========================================================================

    output$region_selector <- renderUI({
      if (identical(input$mode, "Composition")) return(NULL)
      selectInput(
        ns("region"), "Region",
        choices = c("Islet" = "core", "Peri-Islet" = "band", "Islet+20\u00b5m" = "union"),
        selected = "band"
      )
    })

    output$dynamic_selector <- renderUI({
      if (identical(input$mode, "Targets")) {
        classes <- prepared()$targets_all %>%
          dplyr::pull(class) %>% unique() %>% na.omit() %>% sort()
        if (length(classes) == 0) classes <- character(0)
        sel <- if (!is.null(input$which) && input$which %in% classes) {
          input$which
        } else if (length(classes) > 0) {
          classes[1]
        } else {
          NULL
        }
        selectInput(ns("which"), "Target class", choices = classes, selected = sel)
      } else if (identical(input$mode, "Markers")) {
        markers <- prepared()$markers_all %>%
          dplyr::pull(marker) %>% unique() %>% na.omit() %>% sort()
        if (length(markers) == 0) markers <- character(0)
        sel <- if (!is.null(input$which) && input$which %in% markers) {
          input$which
        } else if (length(markers) > 0) {
          markers[1]
        } else {
          NULL
        }
        selectInput(ns("which"), "Marker", choices = markers, selected = sel)
      } else {
        base_choices <- c("Ins_frac" = "Ins_any", "Glu_frac" = "Glu_any", "Stt_frac" = "Stt_any")
        prop_cols <- grep("^prop_", colnames(prepared()$comp), value = TRUE)
        if (length(prop_cols) > 0) {
          prop_labels <- gsub("^prop_", "", prop_cols)
          choices <- list(
            "Hormone Fractions" = base_choices,
            "Cell Type Proportions" = setNames(prop_cols, prop_labels)
          )
        } else {
          choices <- base_choices
        }
        selectInput(ns("which"), "Composition measure", choices = choices,
                    selected = if (!is.null(input$which) && input$which %in% unlist(choices)) input$which else "Ins_any")
      }
    })

    output$metric_selector <- renderUI({
      if (identical(input$mode, "Targets")) {
        radioButtons(ns("target_metric"), "Targets metric",
                     choices = c("Counts", "Density"),
                     selected = "Density", inline = TRUE)
      } else if (identical(input$mode, "Markers")) {
        radioButtons(ns("marker_metric"), "Markers metric",
                     choices = c("Counts", "% positive"),
                     selected = "% positive", inline = TRUE)
      } else {
        NULL
      }
    })

    output$aab_warning <- renderUI({
      df <- raw_df_base()
      if (!is.null(df) && nrow(df) > 0) {
        aabp_donors <- df %>% dplyr::filter(donor_status == "Aab+")
        if (nrow(aabp_donors) == 0 && "Aab+" %in% input$groups) {
          div(style = "color: #d9534f; font-weight: bold; margin-bottom: 10px;",
              "\u26a0 No matching donors.")
        } else {
          NULL
        }
      } else {
        NULL
      }
    })

    # ---- Age & Gender filter UIs (visible only when demographics available) ----
    output$age_filter_ui <- renderUI({
      pd <- prepared()
      if (is.null(pd$comp) || !("age" %in% colnames(pd$comp))) return(NULL)
      age_vals <- as.numeric(pd$comp$age)
      age_vals <- age_vals[is.finite(age_vals)]
      if (length(age_vals) == 0) return(NULL)
      sliderInput(ns("age_range"), "Donor Age (years)",
                  min = floor(min(age_vals)), max = ceiling(max(age_vals)),
                  value = c(floor(min(age_vals)), ceiling(max(age_vals))), step = 1)
    })

    output$gender_filter_ui <- renderUI({
      pd <- prepared()
      if (is.null(pd$comp) || !("gender" %in% colnames(pd$comp))) return(NULL)
      genders <- sort(unique(na.omit(as.character(pd$comp$gender))))
      if (length(genders) == 0) return(NULL)
      checkboxGroupInput(ns("gender_filter"), "Gender", choices = genders, selected = genders, inline = TRUE)
    })

    # ===========================================================================
    # Reactive data chain:  raw_df_base -> raw_df -> plot_df -> summary_df
    # ===========================================================================

    # ---- raw_df_base: per-islet dataset aligned with current selections ----
    raw_df_base <- reactive({
      pd <- prepared()
      mode <- input$mode %||% "Composition"
      groups <- input$groups
      if (is.null(groups) || !length(groups)) {
        groups <- c("ND", "Aab+", "T1D")
      }
      w <- input$which
      if (is.null(w) || length(w) == 0) w <- NA_character_

      # Resolve defaults before branching so initial renders don't error
      if (identical(mode, "Composition")) {
        valid_comp <- c("Ins_any", "Glu_any", "Stt_any", grep("^prop_", colnames(pd$comp), value = TRUE))
        if (!is.character(w) || length(w) != 1 || !nzchar(w) ||
            !(w %in% valid_comp)) {
          w <- "Ins_any"
        }
      } else if (identical(mode, "Targets")) {
        if (!is.character(w) || length(w) != 1 || !nzchar(w)) {
          candidate_classes <- pd$targets_all %>%
            dplyr::pull(class) %>% unique() %>% na.omit() %>% sort()
          if (length(candidate_classes) > 0) w <- candidate_classes[1]
        }
      } else if (identical(mode, "Markers")) {
        if (!is.character(w) || length(w) != 1 || !nzchar(w)) {
          candidate_markers <- pd$markers_all %>%
            dplyr::pull(marker) %>% unique() %>% na.omit() %>% sort()
          if (length(candidate_markers) > 0) w <- candidate_markers[1]
        }
      }

      region_val <- input$region %||% "band"

      if (identical(mode, "Targets")) {
        req(!is.null(w) && nzchar(w))
        region_tag <- paste0("islet_", tolower(region_val))
        df <- pd$targets_all %>%
          dplyr::filter(`Donor Status` %in% groups,
                        class == w,
                        tolower(type) == region_tag) %>%
          dplyr::filter(is.finite(islet_diam_um))
        df <- df %>%
          dplyr::mutate(value = as.numeric(
            if (!is.null(input$target_metric) && input$target_metric == "Counts") count
            else area_density
          ))
      } else if (identical(mode, "Markers")) {
        req(!is.null(w) && nzchar(w))
        region_tag <- paste0("islet_", tolower(region_val))
        df <- pd$markers_all %>%
          dplyr::filter(`Donor Status` %in% groups,
                        marker == w,
                        tolower(region_type) == region_tag) %>%
          dplyr::filter(is.finite(islet_diam_um))
        if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
          df <- df %>% dplyr::mutate(value = suppressWarnings(as.numeric(pos_count)))
        } else {
          num <- suppressWarnings(as.numeric(df$pos_count))
          den <- suppressWarnings(as.numeric(df$n_cells))
          df <- df %>%
            dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0,
                                         100.0 * num / den, NA_real_))
        }
      } else {
        # Composition
        df <- pd$comp %>%
          dplyr::filter(`Donor Status` %in% groups) %>%
          dplyr::filter(is.finite(islet_diam_um))
        if (startsWith(w, "prop_")) {
          # Phenotype proportions are already 0-1, scale to percentage
          df <- df %>% dplyr::mutate(value = suppressWarnings(as.numeric(.data[[w]])) * 100.0)
        } else {
          # Original hormone fractions: count / cells_total * 100
          num <- suppressWarnings(as.numeric(df[[w]]))
          den <- suppressWarnings(as.numeric(df$cells_total))
          df <- df %>%
            dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0,
                                         100.0 * num / den, NA_real_))
        }
      }

      out <- df %>%
        dplyr::mutate(donor_status = `Donor Status`) %>%
        dplyr::filter(!is.na(value))

      # Apply diameter range filter
      if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
        diam_min <- as.numeric(input$diam_range[1])
        diam_max <- as.numeric(input$diam_range[2])
        if (is.finite(diam_min) && is.finite(diam_max)) {
          out <- out %>%
            dplyr::filter(is.finite(islet_diam_um) &
                          islet_diam_um >= diam_min &
                          islet_diam_um <= diam_max)
        }
      }

      # Age filter
      if (!is.null(input$age_range) && length(input$age_range) == 2 && "age" %in% colnames(out)) {
        out <- out[is.finite(as.numeric(out$age)) &
                   as.numeric(out$age) >= input$age_range[1] &
                   as.numeric(out$age) <= input$age_range[2], , drop = FALSE]
      }
      # Gender filter
      if (!is.null(input$gender_filter) && "gender" %in% colnames(out)) {
        out <- out[as.character(out$gender) %in% input$gender_filter, , drop = FALSE]
      }

      # Apply AAb filters (only within Aab+ donor group)
      flags <- input$aab_flags
      all_aab_cols <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")

      others <- out %>% dplyr::filter(donor_status != "Aab+")
      aabp   <- out %>% dplyr::filter(donor_status == "Aab+")

      if (nrow(aabp) > 0) {
        if (is.null(flags) || length(flags) == 0) {
          aabp <- aabp[0, , drop = FALSE]
        } else {
          unchecked <- setdiff(all_aab_cols, flags)
          unchecked_avail <- intersect(unchecked, colnames(aabp))
          if (length(unchecked_avail) > 0) {
            exclude_mat <- as.data.frame(aabp[, unchecked_avail, drop = FALSE])
            for (cc in colnames(exclude_mat)) exclude_mat[[cc]] <- as.logical(exclude_mat[[cc]])
            has_unchecked <- rowSums(exclude_mat, na.rm = TRUE) > 0
            aabp <- aabp[!has_unchecked, , drop = FALSE]
          }
        }
      }

      out <- dplyr::bind_rows(others, aabp)
      attr(out, "selection_used") <- w
      attr(out, "mode_used") <- mode
      out
    })

    # ---- raw_df: applies exclude_zero_top filter ----
    raw_df <- reactive({
      out <- raw_df_base()
      sel_used <- attr(out, "selection_used")
      mode_used <- attr(out, "mode_used")
      if (isTRUE(input$exclude_zero_top)) {
        out <- out %>% dplyr::filter(value != 0)
      }
      attr(out, "selection_used") <- sel_used
      attr(out, "mode_used") <- mode_used
      out
    })

    # ---- plot_df: bins data and marks outliers ----
    plot_df <- reactive({
      df <- raw_df()
      if (is.null(df) || !nrow(df)) {
        return(df %||% data.frame())
      }
      mode <- attr(df, "mode_used") %||% (input$mode %||% "Composition")
      selection_label <- attr(df, "selection_used") %||% input$which
      if (is.null(selection_label) || length(selection_label) == 0 ||
          is.na(selection_label[1]) || !nzchar(selection_label[1])) {
        selection_label <- if (identical(mode, "Composition")) "Ins_any" else ""
      }
      selection_label <- selection_label[1]

      bw <- suppressWarnings(as.numeric(input$binwidth))
      if (length(bw) != 1 || !is.finite(bw) || bw <= 0) bw <- 50

      # Identify outliers (>3 SD) -- mark but do NOT remove
      value_mean <- mean(df$value, na.rm = TRUE)
      value_sd   <- sd(df$value, na.rm = TRUE)
      outlier_threshold <- 3

      df$is_outlier <- abs(df$value - value_mean) > (outlier_threshold * value_sd)
      outliers <- df[df$is_outlier & !is.na(df$is_outlier), ]

      # Store outliers for the sidebar table
      plot_outliers(
        if (nrow(outliers) > 0) {
          data.frame(
            Mode          = mode,
            Selection     = selection_label,
            Case_ID       = outliers$`Case ID`,
            Islet         = outliers$islet_key,
            Donor_Status  = outliers$donor_status,
            Diameter_um   = round(outliers$islet_diam_um, 1),
            Value         = round(outliers$value, 3),
            Z_Score       = round((outliers$value - value_mean) / value_sd, 2),
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      )

      cat(sprintf("[PLOT] Identified %d outliers (>3 SD) - will be colored red\n", nrow(outliers)))

      # Attach diameter bins
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      df
    })

    # ---- summary_df: mean/SE (or SD/IQR) per bin per group ----
    summary_df <- reactive({
      df <- plot_df()
      if (is.null(df) || !nrow(df)) return(data.frame())
      stat_choice <- input$stat %||% "mean_se"
      if (!stat_choice %in% c("mean_se", "mean_sd", "median_iqr")) {
        stat_choice <- "mean_se"
      }
      summary_stats(df,
                    group_cols = c("donor_status", "diam_bin", "diam_mid"),
                    value_col  = "value",
                    stat       = stat_choice)
    })

    # ===========================================================================
    # Outputs: scatter plot (plt), distribution (dist), dist_ui, outlier info
    # ===========================================================================

    # ---- Main scatter plot ----
    output$plt <- renderPlotly({
      sm <- summary_df()
      if (is.null(sm) || !nrow(sm)) {
        return(plotly_empty() %>% layout(title = "Loading data..."))
      }

      grp_levels <- c("ND", "Aab+", "T1D")
      sm$donor_status <- factor(sm$donor_status, levels = grp_levels)

      color_by <- input$plot_color_by %||% "donor_status"

      # Summary lines always use donor_status colors
      color_map <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")

      # Y-axis label
      ylab <- if (identical(input$mode, "Targets")) {
        if (!is.null(input$target_metric) && input$target_metric == "Counts") {
          "Target count"
        } else {
          "Target density (per \u00b5m\u00b2)"
        }
      } else if (identical(input$mode, "Markers")) {
        if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
          "Positive cell count"
        } else {
          "n positive / total (%)"
        }
      } else {
        "% composition"
      }

      # Title
      title_text <- if (identical(input$mode, "Targets")) {
        paste0(input$which, " vs Islet Size")
      } else if (identical(input$mode, "Markers")) {
        paste0(input$which, " vs Islet Size")
      } else {
        wsel <- input$which
        if (is.null(wsel) || length(wsel) != 1 || !nzchar(wsel)) wsel <- "Ins_any"
        if (startsWith(wsel, "prop_")) {
          nm <- paste0(gsub("^prop_", "", wsel), " proportion")
        } else {
          nm <- switch(wsel,
                       Ins_any = "Insulin+ fraction",
                       Glu_any = "Glucagon+ fraction",
                       Stt_any = "Somatostatin+ fraction",
                       wsel)
        }
        paste0(nm, " vs Islet Size")
      }

      # Base plot: summary lines + error bars
      p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status)) +
        geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1)) +
        geom_point() +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
        labs(x = "Islet diameter (\u00b5m)", y = ylab, color = NULL, title = title_text) +
        scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE) +
        theme_minimal(base_size = 14) +
        theme(legend.position = c(1, 1),
              legend.justification = c(1, 1),
              legend.direction = "vertical")

      # X-axis scale
      xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
      if (!is.finite(xmax)) xmax <- 300

      if (!is.null(input$log_scale_x) && input$log_scale_x) {
        max_pow <- ceiling(log2(xmax))
        major_breaks <- 2^(0:max_pow)
        p <- p + scale_x_continuous(trans = "log2", breaks = major_breaks)
      } else {
        major_breaks <- seq(0, ceiling(xmax / 50) * 50, by = 50)
        minor_breaks <- seq(0, ceiling(xmax / 10) * 10, by = 10)
        p <- p + scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks)
      }

      # Y-axis log scale (pseudo-log keeps zeros visible)
      if (!is.null(input$log_scale) && input$log_scale) {
        p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                                    breaks = pseudo_log_breaks(10))
      }

      # Individual points layer
      if (isTRUE(input$show_points)) {
        raw <- plot_df()
        raw$donor_status <- factor(raw$donor_status, levels = grp_levels)

        # Click key for segmentation viewer
        if ("Case ID" %in% names(raw) && "islet_key" %in% names(raw)) {
          raw$islet_click_key <- paste(raw$`Case ID`, raw$islet_key, sep = "|")
        } else {
          raw$islet_click_key <- NA_character_
        }

        # Separate normal vs outliers
        raw_normal   <- raw[!raw$is_outlier | is.na(raw$is_outlier), ]
        raw_outliers <- raw[raw$is_outlier & !is.na(raw$is_outlier), ]

        # Jitter width
        if (!is.null(input$log_scale_x) && input$log_scale_x) {
          jw <- 0.15
        } else {
          jw <- max(1, as.numeric(input$binwidth) * 0.35)
        }

        # Vertical jitter
        yr <- suppressWarnings(range(raw$value, na.rm = TRUE))
        ydiff <- if (all(is.finite(yr))) diff(yr) else 0
        is_counts <- (identical(input$mode, "Targets") &&
                        !is.null(input$target_metric) && input$target_metric == "Counts") ||
                     (identical(input$mode, "Markers") &&
                        !is.null(input$marker_metric) && input$marker_metric == "Counts")
        jh <- if (is_counts) {
          mx <- max(0.02 * ydiff, 0.5)
          if (!is.finite(mx) || mx <= 0) 0.5 else mx
        } else {
          mx <- 0.01 * ydiff
          if (!is.finite(mx) || mx < 0) 0 else mx
        }

        # Normal points layer
        if (color_by == "donor_id") {
          if ("Case ID" %in% colnames(raw_normal)) {
            raw_normal$donor_id <- factor(raw_normal$`Case ID`)
            donor_colors <- get_donor_palette(levels(raw_normal$donor_id))
            combined_colors <- c(
              "ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd",
              donor_colors
            )
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_id, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh),
                         size = input$pt_size, alpha = input$pt_alpha,
                         inherit.aes = FALSE) +
              scale_color_manual(values = combined_colors,
                                 breaks = levels(raw_normal$donor_id),
                                 name = "Donor ID")
          } else {
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh),
                         size = input$pt_size, alpha = input$pt_alpha,
                         inherit.aes = FALSE)
          }
        } else {
          p <- p +
            geom_point(data = raw_normal,
                       aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh),
                       size = input$pt_size, alpha = input$pt_alpha,
                       inherit.aes = FALSE)
        }

        # Outliers on top (slightly larger, more visible)
        if (nrow(raw_outliers) > 0) {
          p <- p +
            geom_point(data = raw_outliers,
                       aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh),
                       size = input$pt_size * 1.2,
                       alpha = min(1.0, input$pt_alpha * 1.5),
                       inherit.aes = FALSE)
        }
      }

      # Minor grid
      p <- p + theme(panel.grid.minor = element_line(size = 0.2, colour = "#eeeeee"))

      # LOESS smoothing overlay
      if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
        p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 0.9)
      }

      # Convert to plotly with click source (namespaced for module compatibility)
      gg <- ggplotly(p, source = ns("plot_scatter"))
      gg <- gg %>%
        layout(legend = list(orientation = "h", x = 0, y = -0.15)) %>%
        event_register("plotly_click")
      gg
    })

    # ---- Distribution controls ----
    output$dist_ui <- renderUI({
      tagList(
        fluidRow(
          column(4,
            sliderInput(ns("dist_pt_size"), "Point size",
                        min = 0.3, max = 4.0, value = 0.7, step = 0.1),
            sliderInput(ns("dist_pt_alpha"), "Point transparency",
                        min = 0.05, max = 1.0, value = 0.25, step = 0.05)
          ),
          column(4,
            selectInput(ns("dist_color_by"), "Color points by:",
                        choices = c("Donor Status" = "donor_status",
                                    "Donor ID" = "donor_id"),
                        selected = "donor_status"),
            radioButtons(ns("dist_type"), "Plot type",
                         choices = c("Violin", "Box"),
                         selected = "Violin", inline = TRUE)
          ),
          column(4,
            checkboxInput(ns("exclude_zero_dist"), "Exclude zero values", value = FALSE),
            checkboxInput(ns("dist_show_points"), "Show individual points", value = TRUE),
            checkboxInput(ns("dist_log_scale"), "Log scale y-axis", value = FALSE)
          )
        )
      )
    })

    # ---- Distribution violin/box plot ----
    output$dist <- renderPlotly({
      rdf <- raw_df()
      if (is.null(rdf) || nrow(rdf) == 0) return(NULL)

      req(input$mode, input$which)

      # Force reactivity on distribution-specific controls
      dist_color_by_val <- input$dist_color_by
      dist_show_points_val <- input$dist_show_points

      cat(sprintf("[DIST] Rendering with mode=%s, which=%s, region=%s, color_by=%s\n",
                  input$mode %||% "NULL",
                  input$which %||% "NULL",
                  input$region %||% "NULL",
                  dist_color_by_val %||% "NULL"))

      # Apply exclude_zero_dist separately (raw_df already applies exclude_zero_top)
      if (isTRUE(input$exclude_zero_dist) && !isTRUE(input$exclude_zero_top)) {
        rdf <- rdf %>% dplyr::filter(value != 0)
      }

      color_by <- input$dist_color_by %||% "donor_status"

      # Factor ordering
      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND", "Aab+", "T1D"))

      # Create donor_id column before it's needed
      if (color_by == "donor_id" && "Case ID" %in% colnames(rdf)) {
        rdf$donor_id <- factor(rdf$`Case ID`)
        cat(sprintf("[DIST] Created donor_id with %d levels: %s\n",
                    length(levels(rdf$donor_id)),
                    paste(levels(rdf$donor_id), collapse = ", ")))
      }

      # Y-axis label
      ylab <- if (identical(input$mode, "Targets")) {
        if (!is.null(input$target_metric) && input$target_metric == "Counts") {
          "Target count"
        } else {
          "Target density (per \u00b5m\u00b2)"
        }
      } else if (identical(input$mode, "Markers")) {
        if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
          "Positive cell count"
        } else {
          "n positive / total (%)"
        }
      } else {
        "% composition"
      }

      # Base violin or box
      g <- ggplot(rdf, aes(x = donor_status, y = value, fill = donor_status))

      if (!is.null(input$dist_type) && input$dist_type == "Box") {
        g <- g + geom_boxplot(
          outlier.shape  = NA,
          outlier.size   = 0,
          outlier.colour = NA,
          outlier.fill   = NA,
          outlier.alpha  = 0,
          outlier.stroke = 0,
          alpha = 0.8
        )
      } else {
        g <- g + geom_violin(trim = FALSE, alpha = 0.8)
      }

      g <- g +
        scale_fill_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"),
                          guide = "none")

      # Individual points
      if (isTRUE(input$dist_show_points)) {
        if (color_by == "donor_id" && "donor_id" %in% colnames(rdf)) {
          donor_colors <- get_donor_palette(levels(rdf$donor_id))
          cat(sprintf("[DIST] Applying donor_id colors for %d donors\n", length(donor_colors)))
          g <- g +
            geom_jitter(aes(x = donor_status, y = value, color = donor_id),
                        width = 0.18,
                        alpha = ifelse(is.null(input$dist_pt_alpha), 0.5, input$dist_pt_alpha),
                        size  = ifelse(is.null(input$dist_pt_size), 1.5, input$dist_pt_size),
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = donor_colors,
                               breaks = levels(rdf$donor_id),
                               name = "Donor ID",
                               guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          g <- g +
            geom_jitter(aes(x = donor_status, y = value, color = donor_status),
                        width = 0.18,
                        alpha = ifelse(is.null(input$dist_pt_alpha), 0.25, input$dist_pt_alpha),
                        size  = ifelse(is.null(input$dist_pt_size), 0.7, input$dist_pt_size),
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"),
                               guide = "none")
        }
      }

      g <- g +
        labs(x = "Donor Status", y = ylab, title = "Distribution across donor groups") +
        theme_minimal(base_size = 14) +
        theme(legend.position = if (color_by == "donor_id") "right" else "none")

      if (!is.null(input$dist_log_scale) && input$dist_log_scale) {
        g <- g + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                                    breaks = pseudo_log_breaks(10))
      }

      p <- ggplotly(g, tooltip = c("x", "y", "colour"))

      # Suppress plotly's built-in box outlier markers
      if (!is.null(input$dist_type) && input$dist_type == "Box") {
        for (i in seq_along(p$x$data)) {
          if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "box") {
            p$x$data[[i]]$marker$opacity      <- 0
            p$x$data[[i]]$marker$size          <- 0
            p$x$data[[i]]$marker$color         <- "rgba(0,0,0,0)"
            p$x$data[[i]]$marker$outliercolor  <- "rgba(0,0,0,0)"
            p$x$data[[i]]$boxpoints            <- FALSE
          }
        }
      }

      p <- p %>% layout(showlegend = isTRUE(color_by == "donor_id"))
      p
    })

    # ---- Outlier info table in sidebar ----
    output$plot_outlier_info <- renderUI({
      show_table <- input$show_plot_outlier_table %||% FALSE
      outlier_data <- plot_outliers()

      if (show_table && !is.null(outlier_data) && nrow(outlier_data) > 0) {
        tags$div(
          style = "margin-top: 15px; padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
          tags$h6(
            style = "margin-top: 0; color: #856404;",
            sprintf("\u26a0\ufe0f %d Outlier%s Removed (>3 SD)",
                    nrow(outlier_data),
                    ifelse(nrow(outlier_data) > 1, "s", ""))),
          tags$div(
            style = "max-height: 200px; overflow-y: auto;",
            renderTable(outlier_data)
          )
        )
      } else {
        NULL
      }
    })

    # ===========================================================================
    # Click handler: scatter plot -> embedded segmentation panel
    # ===========================================================================

    observeEvent(event_data("plotly_click", source = ns("plot_scatter")), {
      click_data <- event_data("plotly_click", source = ns("plot_scatter"))
      if (is.null(click_data)) return()

      # Guard: sf package required
      if (!SF_AVAILABLE) {
        showNotification("Segmentation requires 'sf' package", type = "warning", duration = 5)
        return()
      }

      if (is.null(islet_spatial_lookup)) {
        showNotification("Islet spatial lookup data not available", type = "warning", duration = 5)
        return()
      }

      # Parse key (format: "case_id|islet_key")
      click_key <- click_data$key
      if (is.null(click_key) || is.na(click_key) || !nzchar(click_key)) {
        showNotification("Click on individual points to view segmentation (enable 'Show individual points')",
                         type = "message", duration = 3)
        return()
      }

      key_parts <- strsplit(click_key, "\\|")[[1]]
      if (length(key_parts) < 2) {
        showNotification("Could not parse islet identifier", type = "warning", duration = 3)
        return()
      }

      case_id   <- key_parts[1]
      islet_key <- key_parts[2]

      # Look up centroid from spatial lookup
      spatial_match <- islet_spatial_lookup[
        islet_spatial_lookup$case_id == case_id &
          islet_spatial_lookup$islet_key == islet_key, ]

      if (nrow(spatial_match) == 0) {
        # Try zero-padded variant (e.g. "112" -> "0112")
        case_id_num <- suppressWarnings(as.integer(case_id))
        if (!is.na(case_id_num)) {
          case_id_padded <- sprintf("%04d", case_id_num)
          spatial_match <- islet_spatial_lookup[
            islet_spatial_lookup$case_id == case_id_padded &
              islet_spatial_lookup$islet_key == islet_key, ]
          if (nrow(spatial_match) > 0) case_id <- case_id_padded
        }
      }

      if (nrow(spatial_match) == 0) {
        showNotification(paste("Spatial data not found for", islet_key, "in case", case_id),
                         type = "warning", duration = 3)
        return()
      }

      centroid_x <- spatial_match$centroid_x_um[1]
      centroid_y <- spatial_match$centroid_y_um[1]

      # Update shared reactiveVal — triggers embedded segmentation panel below
      selected_islet(list(
        case_id    = case_id,
        islet_key  = islet_key,
        centroid_x = centroid_x,
        centroid_y = centroid_y
      ))
    })

    # ---- Embedded segmentation viewer panel (below Distribution card) ----
    output$segmentation_viewer_panel <- renderUI({
      info <- selected_islet()
      if (is.null(info)) return(NULL)

      div(class = "card", style = "padding: 15px; margin-bottom: 20px; border: 2px solid #0066CC;",
        fluidRow(
          column(8,
            h4(paste("Islet Segmentation:", info$islet_key, "(Case", info$case_id, ")"),
               style = "margin-top: 0; color: #0066CC;")
          ),
          column(4, style = "text-align: right;",
            actionButton(ns("clear_segmentation"), "Close", class = "btn btn-sm btn-outline-secondary")
          )
        ),
        fluidRow(
          column(8,
            plotOutput("islet_segmentation_view", height = "450px")
          ),
          column(4,
            div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; height: 100%;",
              h5("Legend", style = "margin-top: 0; margin-bottom: 15px;"),
              div(style = "margin-bottom: 10px;",
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #0066CC; margin-right: 8px;"),
                  span("Islet boundary")
                ),
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #00CCCC; margin-right: 8px;"),
                  span("Expanded (+20\u00b5m)")
                ),
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #CC00CC; margin-right: 8px;"),
                  span("Nerve")
                ),
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #CC0000; margin-right: 8px;"),
                  span("Capillary")
                ),
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #00AA00; margin-right: 8px;"),
                  span("Lymphatic")
                ),
                div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                  div(style = "width: 25px; height: 4px; background-color: #FFD700; margin-right: 8px;"),
                  span("Selected islet", style = "font-weight: bold;")
                )
              ),
              hr(),
              h6("Islet Info", style = "margin-bottom: 10px;"),
              tags$dl(style = "font-size: 13px;",
                tags$dt("Case ID:"), tags$dd(info$case_id),
                tags$dt("Islet:"), tags$dd(info$islet_key),
                tags$dt("Centroid X:"), tags$dd(paste0(round(info$centroid_x, 1), " \u00b5m")),
                tags$dt("Centroid Y:"), tags$dd(paste0(round(info$centroid_y, 1), " \u00b5m"))
              ),
              hr(),
              p(style = "font-size: 11px; color: #666; margin-bottom: 0;",
                "Click another point to view a different islet, or click Close to hide this panel.")
            )
          )
        )
      )
    })

    # ---- Clear segmentation viewer ----
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # ===========================================================================
    # Download handlers
    # ===========================================================================

    output$dl_summary <- downloadHandler(
      filename = function() {
        paste0("summary_", gsub("[^0-9A-Za-z]+", "_", Sys.time()), ".csv")
      },
      content = function(file) {
        df <- summary_df()
        if (is.null(df) || nrow(df) == 0) df <- data.frame()

        descriptions <- get_selection_description()

        # Write descriptions as comments, then the data
        writeLines(paste("#", descriptions), file)
        writeLines("", file, sep = "\n")
        write.table(df, file, append = TRUE, sep = ",", row.names = FALSE)
      }
    )

    # ===========================================================================
    # Selection description helper (also returned for Statistics download)
    # ===========================================================================

    get_selection_description <- function() {
      desc <- c()
      desc <- c(desc, paste("Generated:", Sys.time()))
      desc <- c(desc, paste("Focus:", input$mode %||% "Unknown"))
      desc <- c(desc, paste("Region:", input$region %||% "Unknown"))

      if (!is.null(input$dynamic_metric)) {
        desc <- c(desc, paste("Metric:", input$dynamic_metric))
      }

      desc <- c(desc, paste("Donor Groups:", paste(input$groups %||% c(), collapse = ", ")))

      if (length(input$aab_flags) > 0) {
        all_aab <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
        excluded <- setdiff(all_aab, input$aab_flags)
        if (length(excluded) > 0) {
          aab_desc <- paste("AAb Filter: Excluding donors with", paste(excluded, collapse = ", "))
          desc <- c(desc, aab_desc)
        } else {
          desc <- c(desc, "AAb Filter: All Aab+ donors included")
        }
      }

      desc <- c(desc, paste("Statistic:", input$stat %||% "mean_se"))
      desc <- c(desc, paste("Bin Width:", input$binwidth %||% "50", "\u00b5m"))
      if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
        desc <- c(desc, paste("Diameter Range:",
                              input$diam_range[1], "-", input$diam_range[2], "\u00b5m"))
      }

      return(desc)
    }

    # ===========================================================================
    # Return values for other modules (Statistics, etc.)
    # ===========================================================================

    list(
      raw_df                  = raw_df,
      summary_df              = summary_df,
      get_selection_description = get_selection_description
    )
  })
}
