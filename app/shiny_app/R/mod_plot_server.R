# ---------- Plot module server ----------
# Exports: plot_server(id, prepared, selected_islet, active_tab)
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
#   active_tab      -- reactive returning current tab name (e.g. "Plot")

plot_server <- function(id, prepared, selected_islet, active_tab = reactive("Plot"),
                        donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                        remove_outliers = reactive(FALSE),
                        outlier_threshold = reactive(3.0),
                        palette = reactive(PHENOTYPE_COLORS)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Build each selector vocabulary once from server-owned data, then use the
    # same values for both UI construction and browser-input validation. This
    # prevents a forged value from selecting a hidden data column or reaching a
    # label/export path after the data reactive has already fallen back.
    microenvironment_selection_choices <- function(pd) {
      classes <- pd$targets_all %>%
        dplyr::pull(class) %>% unique() %>% na.omit() %>% sort()
      list(
        "Core" = setNames(paste0(classes, "|islet_core"), classes),
        "Peri-Islet" = setNames(paste0(classes, "|islet_band"), classes),
        "Core + Peri" = setNames(paste0(classes, "|islet_union"), classes)
      )
    }

    composition_selection_choices <- function(pd) {
      plot_composition_selection_choices(
        colnames(pd$comp), pd$combined_pairs, comp = pd$comp
      )
    }

    selection_values <- function(choices) {
      unname(unlist(choices, use.names = FALSE))
    }

    selection_vocabularies <- reactive({
      pd <- prepared()
      list(
        composition = composition_selection_choices(pd),
        microenvironment = microenvironment_selection_choices(pd)
      )
    })

    validated_mode <- reactive({
      value <- input$mode %||% "Cell Populations"
      if (!islet_valid_scalar_choice(
        value, c("Cell Populations", "Microenvironment"), "plot_mode_invalid"
      )) return("Cell Populations")
      value
    })

    validated_selection <- reactive({
      vocabularies <- selection_vocabularies()
      composition_values <- selection_values(vocabularies$composition)
      microenvironment_values <- selection_values(vocabularies$microenvironment)
      mode <- validated_mode()
      value <- input$which
      allowed <- if (identical(mode, "Microenvironment")) {
        microenvironment_values
      } else {
        composition_values
      }
      stale_allowed <- if (identical(mode, "Microenvironment")) {
        composition_values
      } else {
        microenvironment_values
      }

      valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value) && value %in% allowed
      if (valid) return(value)

      # Dynamic selectInput replacement briefly leaves the preceding mode's
      # valid value in input$which. Fall back without treating that normal UI
      # transition as a validation failure.
      stale <- is.character(value) && length(value) == 1L && !is.na(value) &&
        value %in% stale_allowed
      if (!is.null(value) && !stale) {
        islet_log_validation_failure("plot_selection_invalid")
      }

      if (length(allowed)) allowed[[1]] else NULL
    })

    validated_metric <- reactive({
      # Microenvironment has a single implicit area-percentage metric; a stale
      # Composition radio value is irrelevant while that mode is active.
      if (identical(validated_mode(), "Microenvironment")) return("Percentage")
      value <- input$metric %||% "Percentage"
      if (!islet_valid_scalar_choice(
        value, c("Percentage", "Counts", "Density"), "plot_metric_invalid"
      )) return("Percentage")
      value
    })

    validated_groups <- reactive({
      value <- input$groups
      if (is.null(value) || !length(value)) return(c("ND", "Aab+", "T1D"))
      if (!islet_valid_choice_vector(
        value, c("ND", "Aab+", "T1D"), "plot_groups_invalid"
      )) return(c("ND", "Aab+", "T1D"))
      value
    })

    # The Population Measure, Metric, age and sex controls are renderUI outputs, so their inputs reach the
    # server one flush after the first render. The cards used to render on the fallbacks and again on the
    # real values, which doubled the first screen's work and payload (3.5 MB; performance baseline,
    # 2026-10-24). The data chain waits once for the controls that will exist; the flag never goes back.
    plot_controls_mounted <- reactiveVal(FALSE)
    observe({
      if (plot_controls_mounted()) return()
      if (is.null(input$which)) return()
      if (!identical(validated_mode(), "Microenvironment") && is.null(input$metric)) return()
      comp <- prepared()$comp
      if (!is.null(comp)) {
        has_age <- "age" %in% colnames(comp) && any(is.finite(suppressWarnings(as.numeric(comp$age))))
        if (has_age && is.null(input$age_range)) return()
        has_sex <- "gender" %in% colnames(comp) && length(stats::na.omit(unique(as.character(comp$gender)))) > 0
        if (has_sex && is.null(input$gender_filter)) return()
      }
      plot_controls_mounted(TRUE)
    })

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
    # Filters out breaks that are too close in pseudo-log transformed space
    # to prevent label overlap (e.g., 0 and 1 collide on large-range axes).
    pseudo_log_breaks <- function(base = 10) {
      function(limits) {
        max_val <- max(limits)
        if (max_val <= 0) return(0)
        max_pow <- ceiling(log(max_val, base))
        min_pow <- if (max_val < 1) floor(log(max_val, base)) else 0
        candidates <- sort(unique(c(0, base^seq(min_pow, max_pow))))
        candidates <- candidates[candidates <= max_val * 1.1]

        # Greedy filter: drop breaks whose transformed position is too close
        trans_pos <- asinh(candidates / 2) / log(base)
        total_range <- diff(range(trans_pos))
        if (total_range == 0) return(candidates)
        min_gap <- max(0.3, total_range * 0.15)

        keep <- logical(length(candidates))
        keep[1] <- TRUE
        last_kept <- trans_pos[1]
        for (i in seq_along(candidates)[-1]) {
          if (trans_pos[i] - last_kept >= min_gap) {
            keep[i] <- TRUE
            last_kept <- trans_pos[i]
          }
        }
        candidates[keep]
      }
    }

    # ---- Reactive storage for outlier table (replaces global <<- assignment) ----
    plot_outliers <- reactiveVal(NULL)

    # ===========================================================================
    # Dynamic UI outputs
    # ===========================================================================

    output$dynamic_selector <- renderUI({
      vocabularies <- selection_vocabularies()
      selected <- validated_selection()
      if (identical(validated_mode(), "Microenvironment")) {
        choices <- vocabularies$microenvironment
        selectInput(ns("which"), "Structure & Region", choices = choices,
                    selected = selected)
      } else {
        choices <- vocabularies$composition
        tagList(
          selectInput(ns("which"), "Population Measure", choices = choices,
                      selected = selected)
        )
      }
    })

    output$metric_selector <- renderUI({
      if (identical(validated_mode(), "Microenvironment")) {
        # Microenvironment exposes a single metric (% of region area
        # occupied by the structure). No radio shown — Nerve / Capillary /
        # Lymphatic are area-valued, so the area fraction is the only
        # measurement that makes biological sense.
        return(NULL)
      } else {
        # Composition — Density uses QuPath-measured region areas (core_region_um2,
        # peri_region_um2) for all three region variants. Peri-islet area is the
        # QuPath `_exp20um` annotation area, merged in via prep_data().
        w <- validated_selection() %||% ""
        if (startsWith(w, "combined_prop_")) {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of all (core+peri) cells"      = "Percentage",
                                     "Cells per islet (core+peri)"     = "Counts",
                                     "Cells per mm² of core+peri area" = "Density"),
                         selected = "Percentage", inline = TRUE)
          )
        } else if (startsWith(w, "peri_prop_")) {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of peri cells"          = "Percentage",
                                     "Cells per islet (peri)"   = "Counts",
                                     "Cells per mm² of peri area" = "Density"),
                         selected = "Percentage", inline = TRUE)
          )
        } else {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of core cells"            = "Percentage",
                                     "Cells per islet (core)"     = "Counts",
                                     "Cells per mm² of core area" = "Density"),
                         selected = "Percentage", inline = TRUE)
          )
        }
      }
    })

    output$aab_warning <- renderUI({
      df <- raw_df_base()
      if (!is.null(df) && nrow(df) > 0) {
        aabp_donors <- df %>% dplyr::filter(donor_status == "Aab+")
        if (nrow(aabp_donors) == 0 && "Aab+" %in% validated_groups()) {
          div(style = "color: var(--islet-danger-fg); font-weight: bold; margin-bottom: 10px;",
              "\u26a0 No matching donors.")
        } else {
          NULL
        }
      } else {
        NULL
      }
    })

    # ---- Age & Sex filter UIs (visible only when demographics available) ----
    # (data column is `gender`; user-facing label is "Sex")
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
      # User-facing label is "Sex" (project convention); the underlying data
      # column stays `gender` for backwards compatibility. This sidebar is shared
      # by the Plot AND Statistics tabs, so it was the one place in the live UI
      # still showing the word "Gender". Values stay M/F — only the labels change.
      sex_labels <- c(F = "Female", M = "Male")
      checkboxGroupInput(ns("gender_filter"), "Sex",
                         choiceNames  = unname(ifelse(genders %in% names(sex_labels),
                                                      sex_labels[genders], genders)),
                         choiceValues = as.list(genders),
                         selected = genders, inline = TRUE)
    })

    # ===========================================================================
    # Reactive data chain:  raw_df_base -> raw_df -> plot_df -> summary_df
    # ===========================================================================

    # ---- The Data controls' object filters: diameter range, age, sex, autoantibodies (Aab+ donors only) ----
    # Shared by the plot cards (raw_df_base) and the Composition by donor card. `out` needs a donor_status column.
    apply_plot_filters <- function(out) {
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
      # Sex filter (data column is `gender`)
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
      out
    }

    # ---- raw_df_base: per-islet dataset aligned with current selections ----
    raw_df_base <- reactive({
      req(plot_controls_mounted())
      pd <- prepared()
      mode <- validated_mode()
      groups <- validated_groups()
      w <- validated_selection()

      if (identical(mode, "Microenvironment")) {
        req(!is.null(w), nzchar(w))
        # Parse "className|islet_region" encoding
        parts <- strsplit(w, "\\|")[[1]]
        target_class <- parts[1]
        region_tag <- if (length(parts) >= 2) parts[2] else "islet_core"

        df <- pd$targets_all %>%
          dplyr::filter(`Donor Status` %in% groups,
                        class == target_class,
                        tolower(type) == region_tag) %>%
          dplyr::filter(is.finite(islet_diam_um))
        # Microenvironment value = % of region area occupied by the
        # structure (annotation area within the bounds of core / peri /
        # core+peri). Counts are not exposed — Nerve / Capillary / Lymphatic
        # are area-valued (a single nerve may segment into several pieces),
        # so object counts depend on segmentation morphology, not biology.
        df <- df %>%
          dplyr::mutate(value = {
            reg <- suppressWarnings(as.numeric(region_um2))
            a   <- suppressWarnings(as.numeric(area_um2))
            ifelse(is.finite(reg) & reg > 0 & is.finite(a), 100 * a / reg, NA_real_)
          })
      } else {
        # Composition
        metric <- validated_metric()
        df <- pd$comp %>%
          dplyr::filter(`Donor Status` %in% groups) %>%
          dplyr::filter(is.finite(islet_diam_um))
        # QuPath-measured annotation areas (µm²) joined by prep_data().
        # Density uses these directly rather than the π·(d/2)² estimate so
        # core / peri / core+peri are all measured on the same scale.
        core_area_um2 <- if ("core_region_um2" %in% colnames(df)) {
          suppressWarnings(as.numeric(df$core_region_um2))
        } else NA_real_
        peri_area_um2 <- if ("peri_region_um2" %in% colnames(df)) {
          suppressWarnings(as.numeric(df$peri_region_um2))
        } else NA_real_
        if (startsWith(w, "combined_prop_")) {
          stem <- sub("^combined_prop_", "", w)
          # Look up the original prop_* column via the pre-built pair table —
          # column naming diverges (prop_* keeps spaces/+, peri_prop_* sanitizes).
          pair <- if (!is.null(pd$combined_pairs)) {
            pd$combined_pairs[pd$combined_pairs$stem == stem, , drop = FALSE]
          } else NULL
          core_col <- if (!is.null(pair) && nrow(pair) >= 1) pair$core_col[1] else paste0("prop_", stem)
          peri_col <- if (!is.null(pair) && nrow(pair) >= 1) pair$peri_col[1] else paste0("peri_prop_", stem)
          core_prop <- if (core_col %in% colnames(df)) suppressWarnings(as.numeric(df[[core_col]])) else NA_real_
          peri_prop <- if (peri_col %in% colnames(df)) suppressWarnings(as.numeric(df[[peri_col]])) else NA_real_
          core_total <- suppressWarnings(as.numeric(df$total_cells_core))
          peri_total <- if ("total_cells_peri" %in% colnames(df)) suppressWarnings(as.numeric(df$total_cells_peri)) else NA_real_
          # Treat NA cell totals / proportions as 0 so an islet with no peri zone
          # still contributes its core, and vice versa
          core_n <- ifelse(is.finite(core_prop) & is.finite(core_total), core_prop * core_total, 0)
          peri_n <- ifelse(is.finite(peri_prop) & is.finite(peri_total), peri_prop * peri_total, 0)
          ctot <- ifelse(is.finite(core_total), core_total, 0)
          ptot <- ifelse(is.finite(peri_total), peri_total, 0)
          combined_n <- core_n + peri_n
          combined_total <- ctot + ptot
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(combined_n))
          } else if (metric == "Density") {
            ca <- ifelse(is.finite(core_area_um2), core_area_um2, 0)
            pa <- ifelse(is.finite(peri_area_um2), peri_area_um2, 0)
            combined_area <- ca + pa
            df <- df %>% dplyr::mutate(value = ifelse(combined_area > 0,
              combined_n / combined_area * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = ifelse(combined_total > 0,
              100.0 * combined_n / combined_total, NA_real_))
          }
        } else if (startsWith(w, "peri_prop_")) {
          raw_val <- suppressWarnings(as.numeric(df[[w]]))
          # peri_prop_ are relative to peri-islet cells, use total_cells_peri
          peri_total <- if ("total_cells_peri" %in% colnames(df)) {
            suppressWarnings(as.numeric(df$total_cells_peri))
          } else {
            suppressWarnings(as.numeric(df$cells_total))
          }
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(raw_val * peri_total))
          } else if (metric == "Density") {
            df <- df %>% dplyr::mutate(value = ifelse(is.finite(peri_area_um2) & peri_area_um2 > 0,
              raw_val * peri_total / peri_area_um2 * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = raw_val * 100.0)
          }
        } else if (startsWith(w, "prop_")) {
          raw_val <- suppressWarnings(as.numeric(df[[w]]))
          core_total <- suppressWarnings(as.numeric(df$total_cells_core))
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(raw_val * core_total))
          } else if (metric == "Density") {
            df <- df %>% dplyr::mutate(value = ifelse(is.finite(core_area_um2) & core_area_um2 > 0,
              raw_val * core_total / core_area_um2 * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = raw_val * 100.0)
          }
        } else {
          # The server validator only admits retained core/peri phenotype
          # proportion columns here. Refuse unexpected values instead of
          # interpreting a retired QuPath threshold count as a fraction.
          df <- df %>% dplyr::mutate(value = NA_real_)
        }
      }

      out <- df %>%
        dplyr::mutate(donor_status = `Donor Status`) %>%
        dplyr::filter(!is.na(value))

      # No min-cells filter since 2026-10-17: the endocrine-object classes chosen above the tabs (Islet /
      # Cluster) select the objects, through prepared().

      out <- apply_plot_filters(out)
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
      mode <- attr(df, "mode_used") %||% validated_mode()
      selection_label <- attr(df, "selection_used") %||% validated_selection()
      if (is.null(selection_label) || length(selection_label) == 0 ||
          is.na(selection_label[1]) || !nzchar(selection_label[1])) {
        selection_label <- ""
      }
      selection_label <- selection_label[1]

      bw <- suppressWarnings(as.numeric(input$binwidth))
      if (length(bw) != 1 || !is.finite(bw) || bw <= 0) bw <- 50

      # Flag outliers using the shared helper. z-score per donor status,
      # threshold from the live slider. The exclusion flag controls whether
      # outliers are dropped from summary_stats() downstream; ghost points
      # always render.
      df <- flag_outliers(df, value_col = "value", group_col = "donor_status",
                          method = "zscore", threshold = outlier_threshold(),
                          enabled = isTRUE(remove_outliers()))
      summary <- summarize_outlier_filter(df)
      outliers <- df[df$is_outlier, , drop = FALSE]

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
            Z_Score       = round(outliers$outlier_score, 2),
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      )

      attr(df, "outlier_summary") <- summary

      # Attach diameter bins
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      df
    })

    # ---- The chosen measure's axis label and short name (titles), shared by every plot type ----
    measure_axis_label <- function() {
      mode <- validated_mode()
      wsel <- validated_selection() %||% ""
      metric <- validated_metric()
      if (identical(mode, "Microenvironment")) {
        parts <- strsplit(wsel, "\\|")[[1]]
        class_name <- if (length(parts) >= 1 && nzchar(parts[1])) parts[1] else "Structure"
        region_tag <- if (length(parts) >= 2) parts[2] else "islet_core"
        region_lbl <- switch(region_tag,
                             islet_core  = "core",
                             islet_band  = "peri-islet",
                             islet_union = "core+peri",
                             "region")
        paste0(class_name, " area (% of ", region_lbl, " area)")
      } else {
        region_tag <- if (startsWith(wsel, "combined_prop_")) "core+peri"
                       else if (startsWith(wsel, "peri_prop_")) "peri"
                       else "core"
        switch(metric,
          "Counts"     = paste0("Cells per islet (", region_tag, ")"),
          "Density"    = paste0("Cells per mm\u00b2 (", region_tag, " area)"),
          "Percentage" = if (region_tag == "core+peri") "% of all (core+peri) cells"
                         else if (region_tag == "peri")    "% of peri cells"
                         else                              "% of core cells"
        )
      }
    }
    measure_name <- function() {
      wsel <- validated_selection() %||% ""
      if (identical(validated_mode(), "Microenvironment")) return(strsplit(wsel, "\\|")[[1]][1])
      # extract_phenotype_from_column() turns the sanitised peri / combined stems ("Beta_cell",
      # "CHGAplus_hormone-negative") back into the phenotype's label.
      if (startsWith(wsel, "combined_prop_")) paste0(extract_phenotype_from_column(wsel), " core+peri proportion")
      else if (startsWith(wsel, "peri_prop_")) paste0(extract_phenotype_from_column(wsel), " peri-islet proportion")
      else if (startsWith(wsel, "prop_")) paste0(extract_phenotype_from_column(wsel), " proportion")
      else wsel
    }

    # ---- Binned summary: mean/SE (or SD/IQR) per diameter bin and donor status ----
    # Excludes flagged outliers when the global filter is ON; otherwise all rows contribute.
    binned_summary <- function(df) {
      if (is.null(df) || !nrow(df)) return(data.frame())
      stat_choice <- input$stat %||% "mean_se"
      if (!stat_choice %in% c("mean_se", "mean_sd", "median_iqr")) {
        stat_choice <- "mean_se"
      }
      if ("is_outlier" %in% colnames(df) && isTRUE(remove_outliers())) {
        df <- df[!df$is_outlier | is.na(df$is_outlier), , drop = FALSE]
      }
      summary_stats(df,
                    group_cols = c("donor_status", "diam_bin", "diam_mid"),
                    value_col  = "value",
                    stat       = stat_choice)
    }
    # The left card's summary (its CSV when it plots Against diameter; returned for the Statistics tab).
    summary_df <- reactive(binned_summary(plot_df()))

    # ===========================================================================
    # Outputs: scatter plot (plt), distribution (dist), outlier info
    # ===========================================================================

    # ---- Against diameter: every object against its diameter, with the binned summary per donor status ----
    # `raw_rows` are a card's flagged, binned rows; `sm` their binned_summary().
    diameter_figure <- function(raw_rows, sm) {
      if (is.null(sm) || !nrow(sm)) {
        return(plotly_islet_empty("Loading data..."))
      }
      mode <- validated_mode()
      selection <- validated_selection() %||% ""
      metric <- validated_metric()
      # Compact hover per object, from the untransformed rows (the drawn values may be pseudo-log or pinned).
      hover_by_key <- plot_point_hover_map(raw_rows)

      grp_levels <- c("ND", "Aab+", "T1D")
      sm$donor_status <- factor(sm$donor_status, levels = grp_levels)

      color_by <- input$plot_color_by %||% "donor_status"

      # Summary lines always use donor_status colors
      color_map <- donor_colors_reactive()

      # Y-axis label
      ylab <- measure_axis_label()

      # Title (with outlier badge appended)
      base_title <- paste0(measure_name(), " vs Islet Size")
      outlier_summary <- attr(raw_rows, "outlier_summary")
      title_text <- paste0(base_title, outlier_title_suffix(outlier_summary))

      # Pseudo-log transform: apply manually to avoid ggplotly axis corruption
      log_scale_on <- !is.null(input$log_scale) && input$log_scale
      plog <- function(x) asinh(x / 2) / log(10)
      orig_y_vals <- numeric(0)
      if (log_scale_on) {
        orig_y_vals <- c(sm$y, sm$ymin, sm$ymax)
        sm$y    <- plog(sm$y)
        sm$ymin <- plog(sm$ymin)
        sm$ymax <- plog(sm$ymax)
      }
      # Axis range (View controls): fit 99% of objects; the objects above are open triangles at the top edge.
      y_upper <- if (axis_clip() && !log_scale_on) plot_axis_upper(raw_rows$value, keep = sm$y) else NULL
      y_top <- if (is.null(y_upper)) NULL else y_upper * 1.02
      beyond_note <- NULL

      # Initialize base plot
      p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status))

      # Individual points layer FIRST (drawn underneath summary lines)
      show_ind_points <- isTRUE(input$show_points)
      donor_colors <- NULL
      donor_id_breaks <- NULL

      if (show_ind_points) {
        raw <- raw_rows
        raw$donor_status <- factor(raw$donor_status, levels = grp_levels)

        # Object key shared with the distribution plot: hover links the two, a click selects (plot_link_*).
        raw$islet_click_key <- plot_link_keys(raw)

        # Transform y-values for pseudo-log if needed
        if (log_scale_on) {
          orig_y_vals <- c(orig_y_vals, raw$value)
          raw$value <- plog(raw$value)
        }

        # Objects above the fitted axis are pinned to its top (true value in the hover).
        raw$beyond <- FALSE
        raw$value_text <- paste0("Value: ", signif(raw$value, 4))
        if (!is.null(y_top)) {
          raw$beyond <- raw$value > y_top
          raw$value[raw$beyond] <- y_top
        }

        # Separate normal vs outliers
        raw_beyond   <- raw[raw$beyond & !(raw$is_outlier %in% TRUE), ]
        raw_normal   <- raw[(!raw$is_outlier | is.na(raw$is_outlier)) & !raw$beyond, ]
        raw_outliers <- raw[raw$is_outlier & !is.na(raw$is_outlier), ]
        beyond_note <- plot_beyond_note(nrow(raw_beyond))

        # Jitter width
        if (!is.null(input$log_scale_x) && input$log_scale_x) {
          jw <- 0.15
        } else {
          jw <- max(1, as.numeric(input$binwidth) * 0.35)
        }

        # Vertical jitter
        yr <- suppressWarnings(range(raw$value, na.rm = TRUE))
        ydiff <- if (all(is.finite(yr))) diff(yr) else 0
        is_counts <- identical(metric, "Counts")
        jh <- if (is_counts) {
          mx <- max(0.02 * ydiff, 0.5)
          if (!is.finite(mx) || mx <= 0) 0.5 else mx
        } else {
          mx <- 0.01 * ydiff
          if (!is.finite(mx) || mx < 0) 0 else mx
        }

        # Normal points layer. Jitter is seeded so visual positions are
        # stable across reactive re-renders (toggling zero-exclusion or the
        # outlier filter doesn't shuffle the cloud).
        if (color_by == "donor_id") {
          if ("Case ID" %in% colnames(raw_normal)) {
            raw_normal$donor_id <- factor(raw_normal$`Case ID`)
            donor_colors <- get_donor_palette(levels(raw_normal$donor_id))
            donor_id_breaks <- levels(raw_normal$donor_id)
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_id, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh, seed = 17),
                         size = input$pt_size, alpha = 1 - input$pt_alpha,
                         inherit.aes = FALSE)
          } else {
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh, seed = 17),
                         size = input$pt_size, alpha = 1 - input$pt_alpha,
                         inherit.aes = FALSE)
          }
        } else {
          p <- p +
            geom_point(data = raw_normal,
                       aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh, seed = 17),
                       size = input$pt_size, alpha = 1 - input$pt_alpha,
                       inherit.aes = FALSE)
        }

        if (nrow(raw_beyond) > 0) {
          p <- p +
            geom_point(data = raw_beyond,
                       aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key, text = value_text),
                       position = position_jitter(width = jw, height = 0, seed = 23),
                       shape = 2, size = input$pt_size * 1.3, stroke = 1, inherit.aes = FALSE)
        }

        # Ghost outliers: faint grey × markers so excluded points are still
        # visible. Drawn last (on top) so reviewers can see WHERE in the
        # distribution the exclusions land without those points pulling the
        # summary trend lines. Distinct seed so ghost markers don't sit
        # exactly on top of normal points of the same coordinate.
        if (nrow(raw_outliers) > 0) {
          p <- p +
            geom_point(data = raw_outliers,
                       aes(x = diam_mid, y = value, key = islet_click_key),
                       position = position_jitter(width = jw, height = if (is.null(y_top)) jh else 0, seed = 31),
                       shape = 4,
                       size = input$pt_size * 1.1,
                       stroke = 0.9,
                       colour = "grey45",
                       alpha = 0.55,
                       inherit.aes = FALSE,
                       show.legend = FALSE)
        }
      }

      # Summary lines + error bars ON TOP of individual points
      p <- p +
        geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1),
                  linewidth = 1.2) +
        geom_point(shape = 21, aes(fill = donor_status), color = "black", stroke = 0.4, size = 2.5) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0)

      # Color scale (must include donor_id colors when applicable)
      if (!is.null(donor_colors)) {
        combined_colors <- c(color_map, donor_colors)
        p <- p + scale_color_manual(values = combined_colors,
                                     breaks = donor_id_breaks,
                                     name = "Donor ID")
      } else {
        p <- p + scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE)
      }
      # Fill scale for summary points (shape=21 outline)
      p <- p + scale_fill_manual(values = color_map, guide = "none")

      p <- p +
        labs(x = "Islet diameter (\u00b5m)", y = ylab, color = NULL, title = title_text) +
        theme_islet() +
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

      # Y-axis log scale (manual pseudo-log — ggplotly can't handle trans objects)
      if (log_scale_on) {
        orig_range <- range(orig_y_vals[is.finite(orig_y_vals)])
        brks <- pseudo_log_breaks(10)(orig_range)
        p <- p + scale_y_continuous(
          breaks = plog(brks),
          labels = scales::label_comma()(brks)
        )
      }

      # Minor grid
      p <- p + theme(panel.grid.minor = element_line(linewidth = 0.2, colour = "#eeeeee"))

      # LOESS smoothing overlay
      if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
        p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, linewidth = 0.9)
      }

      gg <- ggplotly(p, source = ns("plot_scatter"))
      # Scrub the "colour: black" leak (from the shape-21 summary point's
      # color = "black" parameter) and any duplicate donor_status lines.
      gg <- strip_plotly_color_leak(gg)
      # Payload: compact hover on the object points, five significant figures on the coordinates.
      gg <- plotly_round_coordinates(plotly_compact_point_hover(gg, hover_by_key))
      # Add thin black outlines to all marker traces
      for (i in seq_along(gg$x$data)) {
        tr <- gg$x$data[[i]]
        if (!is.null(tr$mode) && grepl("markers", tr$mode)) {
          gg$x$data[[i]]$marker$line <- list(color = "black", width = 0.5)
        }
      }
      # Clean plotly legend: strip "(ID, N)" -> "ID" from ggplotly trace naming
      for (i in seq_along(gg$x$data)) {
        nm <- gg$x$data[[i]]$name
        if (!is.null(nm)) {
          cleaned <- gsub("^\\((.+),\\s*\\d+\\)$", "\\1", nm)
          gg$x$data[[i]]$name <- cleaned
          gg$x$data[[i]]$legendgroup <- cleaned
        }
      }
      fonts <- plotly_islet_fonts()
      gg <- gg %>%
        layout(
          font = fonts$global,
          title = list(font = fonts$title),
          xaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(orientation = "h", x = 0, y = -0.25,
                        xanchor = "left", yanchor = "top",
                        font = fonts$legend_text,
                        title = list(font = fonts$legend_title)),
          margin = list(b = 100)) %>%
        plotly_islet_config("plot_scatter")
      if (!is.null(y_top)) {
        gg <- gg %>% layout(
          title = list(text = plotly_title_with_sub(title_text, beyond_note), font = fonts$title),
          yaxis = list(range = plot_axis_range(raw_rows$value, y_top), autorange = FALSE),
          margin = list(t = 80, b = 100))
      }

      trajectory_link_plot(gg, ns)
    }

    # ---- Distribution controls ----
    # ---- Distribution violin/box plot ----
    # Note: the right panel reads `raw_df_base()` (not `raw_df()`) so its
    # "Exclude zero values" checkbox is fully independent of the left panel's
    # equivalent \u2014 toggling one no longer affects the other. Outliers are
    # flagged on the dist-specific subset so the helper sees the same
    # distribution it is summarising, and ghost \u00d7 markers render at positions
    # consistent with the current zero-exclusion state.
    # Distribution-panel data, extracted from the renderer so `dl_distribution_csv`
    # exports exactly what this panel shows. NOTE it reads `raw_df_base()` and its
    # OWN `exclude_zero_dist` — it is deliberately independent of the main
    # scatter's `exclude_zero_top`, so exporting `raw_df()` here would be wrong.
    dist_data <- reactive({
      rdf <- raw_df_base()
      if (is.null(rdf) || nrow(rdf) == 0) return(NULL)
      req(validated_mode(), validated_selection())

      # Independent zero-exclusion for the Distribution panel.
      if (isTRUE(input$exclude_zero_dist)) {
        rdf <- rdf %>% dplyr::filter(value != 0)
      }

      color_by <- input$plot_color_by %||% "donor_status"
      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND", "Aab+", "T1D"))

      if (color_by == "donor_id" && "Case ID" %in% colnames(rdf)) {
        rdf$donor_id <- factor(rdf$`Case ID`)
      }

      # Flag outliers on the dist-specific subset so badges and ghost markers
      # reflect the zero-exclusion choice. The shared helper appends
      # `is_outlier` + `outlier_score`; ghost rows are split off and drawn as
      # a separate layer.
      rdf <- flag_outliers(rdf, value_col = "value", group_col = "donor_status",
                           method = "zscore", threshold = outlier_threshold(),
                           enabled = isTRUE(remove_outliers()))
      list(all = rdf, summary = summarize_outlier_filter(rdf))
    })

    # ---- By donor status: a violin or box per donor status, every object as a point ----
    # `dd` is list(all = a card's flagged rows, summary = their outlier summary).
    status_figure <- function(dd) {
      if (is.null(dd) || is.null(dd$all) || !nrow(dd$all)) return(plotly_islet_empty("No data for the current selection"))
      mode <- validated_mode()
      selection <- validated_selection() %||% ""
      metric <- validated_metric()

      # The point controls are shared with the scatter (View controls, "Both plots").
      color_by <- input$plot_color_by %||% "donor_status"
      rdf <- dd$all
      # Compact hover per object, from the rows as flagged (before pinning or ghosting changes their values).
      hover_by_key <- plot_point_hover_map(rdf)
      rdf$donor_status <- factor(as.character(rdf$donor_status), levels = c("ND", "Aab+", "T1D"))
      if (color_by == "donor_id" && "Case ID" %in% colnames(rdf)) rdf$donor_id <- factor(rdf$`Case ID`)
      rdf$link_key <- plot_link_keys(rdf)
      dist_summary <- dd$summary
      rdf_ghost <- rdf[rdf$is_outlier, , drop = FALSE]
      rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), , drop = FALSE]
      if (nrow(rdf_ghost) > 0 && !is.null(rdf_ghost$donor_status)) {
        rdf_ghost$donor_status <- factor(rdf_ghost$donor_status,
                                         levels = c("ND", "Aab+", "T1D"))
      }

      # Y-axis label
      ylab <- measure_axis_label()

      # Axis range (View controls): fit 99% of objects. The violins / boxes are drawn from every value; the points
      # above the axis are pinned to its top as open triangles (true value in the hover).
      y_upper <- if (axis_clip() && !isTRUE(input$log_scale)) plot_axis_upper(dd$all$value) else NULL
      y_top <- if (is.null(y_upper)) Inf else y_upper * 1.02
      rdf$value_text <- paste0("Value: ", signif(rdf$value, 4))
      rdf_pts <- rdf[!(rdf$value > y_top), , drop = FALSE]
      rdf_beyond <- rdf[rdf$value > y_top, , drop = FALSE]
      if (nrow(rdf_beyond)) rdf_beyond$value <- y_top
      if (nrow(rdf_ghost)) rdf_ghost$value <- pmin(rdf_ghost$value, y_top)

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
        scale_fill_manual(values = donor_colors_reactive(),
                          guide = "none")

      # Individual points (kept = non-outliers). Jitter is seeded so the
      # visual positions are stable across reactive re-renders \u2014 toggling
      # other controls won't shuffle the cloud.
      if (isTRUE(input$show_points)) {
        if (color_by == "donor_id" && "donor_id" %in% colnames(rdf)) {
          donor_colors <- get_donor_palette(levels(rdf$donor_id))
          g <- g +
            geom_jitter(data = rdf_pts, aes(x = donor_status, y = value, color = donor_id, key = link_key),
                        position = position_jitter(width = 0.18, seed = 17),
                        alpha = 1 - (input$pt_alpha %||% 0.4),
                        size  = input$pt_size %||% 3,
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = donor_colors,
                               breaks = levels(rdf$donor_id),
                               name = "Donor ID",
                               guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          g <- g +
            geom_jitter(data = rdf_pts, aes(x = donor_status, y = value, color = donor_status, key = link_key),
                        position = position_jitter(width = 0.18, seed = 17),
                        alpha = 1 - (input$pt_alpha %||% 0.4),
                        size  = input$pt_size %||% 3,
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = donor_colors_reactive(),
                               guide = "none")
        }
      }

      if (isTRUE(input$show_points) && nrow(rdf_beyond) > 0) {
        g <- g + geom_point(data = rdf_beyond,
                            aes(x = donor_status, y = value,
                                color = .data[[if (color_by == "donor_id" && "donor_id" %in% colnames(rdf_beyond))
                                                 "donor_id" else "donor_status"]],
                                key = link_key, text = value_text),
                            position = position_jitter(width = 0.18, height = 0, seed = 23), shape = 2,
                            size = (input$pt_size %||% 3) * 1.3, stroke = 1, inherit.aes = FALSE)
      }

      # Ghost \u00d7 layer: excluded outliers always rendered on the Distribution
      # panel so users can see WHERE in the violin/box they live. Drawn last
      # so they sit on top of the violin/box body.
      if (nrow(rdf_ghost) > 0) {
        ghost_size <- (input$pt_size %||% 3) * 1.4
        g <- g + geom_jitter(
          data = rdf_ghost,
          aes(x = donor_status, y = value, key = link_key),
          inherit.aes = FALSE,
          position = position_jitter(width = 0.18, seed = 31),
          shape = 4,
          size = ghost_size,
          stroke = 0.9,
          colour = "grey45",
          alpha = 0.6,
          show.legend = FALSE
        )
      }

      title_text <- paste0("Distribution across donor groups",
                           outlier_title_suffix(dist_summary))
      g <- g +
        labs(x = "Donor Status", y = ylab, title = title_text) +
        theme_islet() +
        theme(legend.position = if (color_by == "donor_id") "right" else "none")

      if (isTRUE(input$log_scale)) {
        g <- g + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                                    breaks = pseudo_log_breaks(10),
                                    labels = scales::label_comma())
      }

      p <- ggplotly(g, tooltip = c("x", "y", "colour", "text"))
      p <- strip_plotly_color_leak(p)
      # Payload: compact hover on the object points, five significant figures on the coordinates. The violin
      # outlines arrive as line traces with one generated hover text per polygon vertex (about 70 KB each,
      # naming the vertex's density value); nothing a reviewer reads, so they carry no hover at all.
      p <- plotly_round_coordinates(plotly_compact_point_hover(p, hover_by_key))
      for (i in seq_along(p$x$data)) {
        tr <- p$x$data[[i]]
        if (identical(tr$mode, "lines") && is.null(tr$key) && length(tr$text) > 1) {
          p$x$data[[i]]$text <- NULL
          p$x$data[[i]]$hoverinfo <- "none"
        }
      }

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

      # Add thin black outlines to jitter point traces
      for (i in seq_along(p$x$data)) {
        tr <- p$x$data[[i]]
        if (!is.null(tr$mode) && grepl("markers", tr$mode)) {
          p$x$data[[i]]$marker$line <- list(color = "black", width = 0.5)
        }
      }

      # Clean plotly legend: strip "(ID, N)" -> "ID"
      for (i in seq_along(p$x$data)) {
        nm <- p$x$data[[i]]$name
        if (!is.null(nm)) {
          cleaned <- gsub("^\\((.+),\\s*\\d+\\)$", "\\1", nm)
          p$x$data[[i]]$name <- cleaned
          p$x$data[[i]]$legendgroup <- cleaned
        }
      }
      fonts <- plotly_islet_fonts()
      p <- p %>%
        layout(
          showlegend = isTRUE(color_by == "donor_id"),
          font = fonts$global,
          title = list(font = fonts$title),
          xaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(font = fonts$legend_text,
                        title = list(font = fonts$legend_title))
        ) %>%
        plotly_islet_config("plot_distribution")
      if (is.finite(y_top)) {
        p <- p %>% layout(
          title = list(text = plotly_title_with_sub(title_text, plot_beyond_note(nrow(rdf_beyond))),
                       font = fonts$title),
          yaxis = list(range = plot_axis_range(dd$all$value, y_top), autorange = FALSE),
          margin = list(t = 80))
      }
      trajectory_link_plot(p, ns)
    }

    # ===========================================================================
    # Plot types (owner, 2026-10-21): either card plots any of PLOT_TYPE_CHOICES (R/plot_types.R). The left card
    # reads plot_df() (its own "Exclude zero values" box) and opens on Against diameter; the right one reads
    # dist_data() and opens on By donor status. Output ids stay plt / dist and the CSVs dl_summary /
    # dl_distribution_csv, each exporting what its card shows.
    # ===========================================================================
    PLOT_CARD_DEFAULT <- c(left = "diameter", right = "status")
    card_type <- function(card) {
      value <- input[[paste0(card, "_type")]]
      if (is.null(value) || !islet_valid_scalar_choice(value, PLOT_TYPE_CHOICES, "plot_type_invalid"))
        PLOT_CARD_DEFAULT[[card]] else value
    }
    card_type_label <- function(card) names(PLOT_TYPE_CHOICES)[match(card_type(card), PLOT_TYPE_CHOICES)]
    output$left_title <- renderText(card_type_label("left"))
    output$right_title <- renderText(card_type_label("right"))

    # A card's rows: flagged for outliers, binned by diameter, with link keys.
    card_rows <- function(card) {
      if (card == "left") {
        df <- plot_df()
      } else {
        dd <- dist_data()
        if (is.null(dd) || is.null(dd$all) || !nrow(dd$all)) return(NULL)
        bw <- suppressWarnings(as.numeric(input$binwidth))
        if (length(bw) != 1 || !is.finite(bw) || bw <= 0) bw <- 50
        df <- bin_islet_sizes(dd$all, "islet_diam_um", bw)
        attr(df, "outlier_summary") <- dd$summary
      }
      if (is.null(df) || !nrow(df)) return(NULL)
      df$donor_status <- as.character(df$donor_status)
      df$link_key <- plot_link_keys(df)
      df
    }
    # The rows the added types draw: flagged outliers left out while the filter is on.
    card_kept_rows <- function(rows) {
      if (isTRUE(remove_outliers()) && "is_outlier" %in% names(rows))
        rows[!(rows$is_outlier %in% TRUE), , drop = FALSE] else rows
    }
    donor_id_palette <- function(rows) get_donor_palette(rows$`Case ID`)

    # The card's x measure (Against another measure) and clinical value (Against a donor clinical value).
    x_measure_choices <- reactive(plot_x_measure_choices(colnames(prepared()$comp), prepared()$combined_pairs,
                                                         comp = prepared()$comp))
    clinical_choices <- reactive({
      cols <- colnames(prepared()$comp)
      PLOT_CLINICAL_CHOICES[PLOT_CLINICAL_CHOICES %in% cols]
    })
    observe({
      xs <- x_measure_choices()
      cl <- clinical_choices()
      for (card in c("left", "right")) {
        x_id <- paste0(card, "_x_measure")
        current <- isolate(input[[x_id]])
        updateSelectInput(session, x_id, choices = xs,
                          selected = if (!is.null(current) && current %in% unlist(xs)) current else
                            intersect(c("immune_frac_peri", "islet_diam_um"), unlist(xs))[1])
        c_id <- paste0(card, "_clinical")
        current <- isolate(input[[c_id]])
        updateSelectInput(session, c_id, choices = cl,
                          selected = if (!is.null(current) && current %in% cl) current else unname(cl[1]))
      }
    })
    card_choice <- function(card, what, allowed) {
      value <- input[[paste0(card, "_", what)]]
      if (is.null(value) || !islet_valid_scalar_choice(value, allowed, "plot_card_choice_invalid")) return(NULL)
      value
    }
    flat_choices <- function(choices) { v <- unlist(unname(choices)); v }

    # View controls > Axis range: "bulk" fits 99% of objects (default), "all" every object.
    axis_clip <- function() !identical(input$axis_range, "all")

    card_figure <- function(card) {
      type <- card_type(card)
      rows <- card_rows(card)
      if (is.null(rows)) {
        return(trajectory_link_plot(plotly_islet_empty(if (card == "left") "Loading data..." else
          "No data for the current selection"), ns))
      }
      status_cols <- donor_colors_reactive()
      color_by <- input$plot_color_by %||% "donor_status"
      log_y <- isTRUE(input$log_scale)
      ylab <- measure_axis_label()
      name <- measure_name()
      kept <- card_kept_rows(rows)
      fig <- switch(type,
        diameter = diameter_figure(rows, binned_summary(rows)),
        status = status_figure(list(all = rows, summary = attr(rows, "outlier_summary"))),
        donor_medians = plot_donor_medians_figure(kept, ylab, name, status_cols, color_by, donor_id_palette(kept),
                                                  log_y, input$pt_size %||% 3),
        ecdf = plot_ecdf_figure(kept, ylab, name, status_cols, log_y, clip = axis_clip()),
        donor_curves = plot_donor_curves_figure(kept, ylab, name, status_cols, log_y, clip = axis_clip()),
        xy = {
          xs <- flat_choices(x_measure_choices())
          x_col <- card_choice(card, "x_measure", xs) %||% unname(xs[1])
          plot_xy_figure(kept, plot_x_values(kept, prepared()$comp, x_col), names(xs)[match(x_col, xs)], ylab, name,
                         status_cols, color_by, donor_id_palette(kept), log_y, input$pt_size %||% 3,
                         input$pt_alpha %||% 0.4, isTRUE(input$show_points), clip = axis_clip())
        },
        clinical = {
          cl <- clinical_choices()
          c_col <- card_choice(card, "clinical", cl) %||% unname(cl[1])
          if (is.null(c_col) || is.na(c_col)) plotly_islet_empty("No donor clinical values in this data") else
            plot_clinical_figure(kept, plot_clinical_values(kept, c_col), names(cl)[match(c_col, cl)], ylab, name,
                                 status_cols, color_by, donor_id_palette(kept), log_y, input$pt_size %||% 3,
                                 input$pt_alpha %||% 0.4, isTRUE(input$show_points), clip = axis_clip(),
                                 # ND and Aab+ donors have no diabetes, so no duration: not a missing value.
                                 missing_label = if (identical(c_col, "t1d_duration_years")) "without T1D" else
                                   "without a recorded value")
        },
        density = plot_density_figure(kept, ylab, name, status_cols, log_y, isTRUE(input$log_scale_x),
                                      clip = axis_clip())
      )
      if (type %in% c("diameter", "status")) fig else trajectory_link_plot(fig, ns)
    }
    output$plt <- renderPlotly(card_figure("left"))
    output$dist <- renderPlotly(card_figure("right"))

    # What a card's CSV holds: list(df, lines) for its current type.
    card_export <- function(card) {
      type <- card_type(card)
      rows <- tryCatch(card_rows(card), error = function(e) NULL)
      per_object <- function(df, extra = NULL) {
        if (is.null(df) || !nrow(df)) return(NULL)
        want <- c("Case ID", "islet_key", "donor_status", "value", extra, "islet_diam_um", "total_cells_core",
                  "total_cells_peri", "is_outlier", "outlier_score")
        df[, intersect(want, colnames(df)), drop = FALSE]
      }
      if (is.null(rows)) return(list(df = NULL, lines = "No rows for the current selection"))
      kept <- card_kept_rows(rows)
      switch(type,
        diameter = list(df = binned_summary(rows), lines = c(
          "Rows are per (donor status x diameter bin). y/ymin/ymax are the plotted values for the chosen Statistic.",
          "Individual objects are not in this file: choose By donor status for per-object values.")),
        status = list(df = per_object(rows), lines =
          "is_outlier marks the ghost x markers. Rows are NOT dropped, so the full distribution is exportable."),
        donor_medians = list(df = donor_median_table(kept), lines =
          "One row per donor: median and quartiles (q1, q3) of the measure over its objects."),
        ecdf = list(df = ecdf_table(kept), lines =
          "Per donor status, each object's value and the proportion of the status's objects at or below it."),
        donor_curves = list(df = per_object(kept), lines =
          "Per-object values behind each donor's curve (a density over the donor's objects)."),
        xy = {
          xs <- flat_choices(x_measure_choices())
          x_col <- card_choice(card, "x_measure", xs) %||% unname(xs[1])
          kept$x_value <- plot_x_values(kept, prepared()$comp, x_col)
          list(df = per_object(kept, "x_value"), lines = paste0("x_value: ", names(xs)[match(x_col, xs)],
                                                                 "; value: the y measure"))
        },
        clinical = {
          cl <- clinical_choices()
          c_col <- card_choice(card, "clinical", cl) %||% unname(cl[1])
          kept$clinical <- plot_clinical_values(kept, c_col)
          list(df = donor_median_table(kept, extra = "clinical"), lines = paste0(
            "One row per donor: median and quartiles of the measure over its objects; clinical: ",
            names(cl)[match(c_col, cl)], " (NA: not recorded; HbA1c 0 is treated as not recorded)"))
        },
        density = list(df = per_object(kept), lines = "Per-object values behind the density contours.")
      )
    }
    write_card_csv <- function(file, card) {
      ex <- card_export(card)
      desc <- tryCatch(get_selection_description(), error = function(e) character(0))
      desc <- desc[!grepl("^Zero-value exclusion", desc)]
      zero <- if (card == "left") input$exclude_zero_top else input$exclude_zero_dist
      write_islet_csv(file, ex$df, c(
        paste0("Plot \u2014 ", card, " card: ", card_type_label(card)),
        paste("Measure:", measure_axis_label()),
        desc,
        paste0("Zero-value exclusion (", card, " plot): ", if (isTRUE(zero)) "ON" else "OFF"),
        ex$lines))
    }

    # ===========================================================================
    # Composition by donor (R/donor_composition.R): every donor's whole cell make-up, pooled over its objects
    # ===========================================================================
    composition_table <- reactive({
      pd <- prepared()
      req(pd$comp, pd$combined_pairs)
      comp <- pd$comp %>%
        dplyr::filter(`Donor Status` %in% validated_groups(), is.finite(islet_diam_um)) %>%
        dplyr::mutate(donor_status = `Donor Status`)
      comp <- apply_plot_filters(comp)
      region <- input$composition_region %||% "core"
      if (!islet_valid_scalar_choice(region, DONOR_COMPOSITION_REGIONS, "composition_region_invalid")) region <- "core"
      donor_composition_table(comp, pd$combined_pairs, region)
    })
    # Leave out: the root's phenotypes as choices; each region's default selection when the region changes.
    composition_leave_out_choices <- reactive({
      labels <- prepared()$combined_pairs$label
      c(labels[phenotype_is_kept(labels)], DONOR_COMPOSITION_OTHER)
    })
    # Only the region resets the choice (toggling an object class must not).
    observeEvent(input$composition_region, {
      region <- input$composition_region %||% "core"
      choices <- isolate(composition_leave_out_choices())
      updateCheckboxGroupInput(session, "composition_leave_out", choices = choices, inline = TRUE,
                               selected = intersect(donor_composition_default_leave_out(region), choices))
    }, ignoreNULL = FALSE)
    composition_leave_out <- function() {
      v <- input$composition_leave_out
      if (is.null(v)) return(character(0))
      if (!islet_valid_choice_vector(v, composition_leave_out_choices(), "composition_leave_out_invalid"))
        return(character(0))
      v
    }
    output$composition_by_donor <- renderPlotly({
      req(plot_controls_mounted())
      tab <- composition_table()
      order <- input$composition_order %||% "beta"
      if (!islet_valid_scalar_choice(order, DONOR_COMPOSITION_ORDERS, "composition_order_invalid")) order <- "beta"
      region <- input$composition_region %||% "core"
      region_label <- names(DONOR_COMPOSITION_REGIONS)[match(region, DONOR_COMPOSITION_REGIONS)] %||% "Core"
      if (is.na(region_label)) region_label <- "Core"
      donor_composition_figure(tab, if (nrow(tab)) donor_composition_order(tab, order) else character(0),
                               palette(), donor_colors_reactive(), tolower(region_label), composition_leave_out())
    })
    output$dl_composition_csv <- downloadHandler(
      filename = function() paste0("composition_by_donor_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        tab <- composition_table()
        desc <- tryCatch(get_selection_description(), error = function(e) character(0))
        desc <- desc[!grepl("^(Data Layer|Selection|Metric|Statistic|Bin Width|Zero-value|Outlier)", desc)]
        write_islet_csv(file, tab, c(
          "Plot \u2014 Composition by donor",
          paste("Region:", input$composition_region %||% "core"),
          desc,
          "cells: the donor's cells of the phenotype, pooled over its selected objects; total_cells: all its cells there",
          "proportion = 100 * cells / total_cells; 'No subtype' holds the cells typed only to a general class"))
      }
    )

    # ---- Outlier info table in sidebar ----
    output$plot_outlier_info <- renderUI({
      show_table <- isTRUE(input$show_outlier_table)
      outlier_data <- plot_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (!show_table) return(NULL)
      if (is.null(outlier_data) || nrow(outlier_data) == 0) {
        return(tags$div(
          style = "margin-top: 12px; padding: 8px 10px; background-color: var(--islet-info-bg); border: 1px solid var(--islet-info-border); border-radius: 5px; color: var(--islet-info-fg); font-size: 12px;",
          "No outliers detected for the current selection."
        ))
      }

      thr <- suppressWarnings(as.numeric(outlier_threshold()))
      if (!is.finite(thr)) thr <- 3
      header <- if (filter_on) {
        sprintf("\u26a0\ufe0f %d Outlier%s Excluded (|z| > %g within donor status)",
                nrow(outlier_data),
                ifelse(nrow(outlier_data) > 1, "s", ""), thr)
      } else {
        sprintf("\u24d8 %d Outlier%s Flagged (filter OFF \u2014 still in summary)",
                nrow(outlier_data),
                ifelse(nrow(outlier_data) > 1, "s", ""))
      }

      tags$div(
        style = "margin-top: 15px; padding: 10px; background-color: var(--islet-warn-bg); border: 1px solid var(--islet-warn-border); border-radius: 5px; color: var(--islet-warn-fg);",
        tags$h6(style = "margin-top: 0; color: var(--islet-warn-fg);", header),
        tags$div(style = "max-height: 200px; overflow-y: auto;",
          renderTable(outlier_data))
      )
    })

    # ===========================================================================
    # Linked points: the size scatter and the distribution plot (Trajectory's linking, www/trajectory-links.js).
    # Hovering an object rings it in both plots in the browser; a click reaches R here, selects the object and
    # opens its segmentation panel. Keys are trajectory_link_key() ("donor|Islet_N"), as on the Trajectory tab, so
    # an object selected there stays ringed here.
    # ===========================================================================
    plot_link_keys <- function(df) {
      if (!all(c("Case ID", "islet_key") %in% names(df))) return(rep(NA_character_, nrow(df)))
      trajectory_link_key(df$`Case ID`, df$islet_key)
    }
    # One compact hover text per object key for both cards: donor and object, donor status, diameter, and the
    # true value. ggplotly() otherwise lists every mapped aesthetic with its variable name, about 110 bytes a
    # point, which was 58 % of each 1.7 MB Plot message (performance baseline, 2026-10-24).
    plot_point_hover_map <- function(df) {
      keys <- plot_link_keys(df)
      ok <- !is.na(keys)
      if (!any(ok)) return(character(0))
      diam <- if ("islet_diam_um" %in% names(df)) suppressWarnings(as.numeric(df$islet_diam_um)) else rep(NA_real_, nrow(df))
      value <- if ("value" %in% names(df)) suppressWarnings(as.numeric(df$value)) else rep(NA_real_, nrow(df))
      text <- paste0(df$`Case ID`, " · ", df$islet_key,
                     "<br>", df$donor_status,
                     ifelse(is.finite(diam), paste0(" · ", round(diam), " µm"), ""),
                     "<br>Value: ", ifelse(is.finite(value), as.character(signif(value, 4)), "NA"))
      stats::setNames(text[ok], keys[ok])
    }
    plot_link_selection_key <- function() {
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      key <- trajectory_link_key(info$case_id, info$islet_key)
      if (length(key) != 1L || is.na(key)) NULL else key
    }
    send_plot_link_state <- function() {
      session$sendCustomMessage(ns("linked_state"), list(
        key = plot_link_selection_key(), active = identical(active_tab(), "Plot")))
    }
    observeEvent(input$linked_click, {
      req(identical(active_tab(), "Plot"))
      event <- input$linked_click
      key <- event$key
      shown <- list(plot_df(), tryCatch(dist_data()$all, error = function(e) NULL))
      rows <- do.call(rbind, lapply(shown, function(d) {
        if (is.null(d) || !nrow(d) || !all(c("Case ID", "islet_key") %in% names(d))) return(NULL)
        data.frame(case_id = as.character(d$`Case ID`), islet_key = as.character(d$islet_key),
                   key = plot_link_keys(d), stringsAsFactors = FALSE)
      }))
      if (!is.character(key) || length(key) != 1L || is.na(key) ||
          !is.character(event$source) || length(event$source) != 1L ||
          !(event$source %in% ns(c("plt", "dist"))) || is.null(rows) || !(key %in% rows$key)) {
        islet_log_validation_failure("plot_link_invalid")
        send_plot_link_state()
        return()
      }
      row <- rows[match(key, rows$key), , drop = FALSE]
      select_islet_from_key(row$case_id, row$islet_key, prepared, selected_islet)
      send_plot_link_state()
    }, ignoreInit = TRUE)
    observeEvent(input$linked_ready, { send_plot_link_state() }, ignoreInit = TRUE)
    observe({ send_plot_link_state() })
    observeEvent(input$clear_link_selection, {
      selected_islet(NULL)
      send_plot_link_state()
    })

    # ---- Embedded segmentation viewer panel (below Distribution card) ----
    output$segmentation_viewer_panel <- renderUI({
      # Only render when Plot tab is active to avoid duplicate non-namespaced output IDs
      if (active_tab() != "Plot") return(NULL)
      info <- selected_islet()
      if (is.null(info)) return(NULL)
      segmentation_viewer_panel_ui(info, ns("clear_segmentation"), layout = "wide")
    })

    # ---- Clear segmentation viewer ----
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # ===========================================================================
    # Download handlers
    # ===========================================================================

    output$dl_summary <- downloadHandler(
      filename = function() paste0("plot_left_", card_type("left"), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write_card_csv(file, "left")
    )

    # The right card's rows come from `dist_data()`, which has its OWN zero-exclusion checkbox
    # (`exclude_zero_dist`), so each card's CSV is built from that card's rows (write_card_csv()).
    # ---------- Endocrine objects per donor (object_counts.R) ----------
    object_counts_data <- reactive({
      d <- load_object_counts()
      if (is.null(d)) return(NULL)
      d[d$donor_status %in% validated_groups(), , drop = FALSE]
    })

    output$object_counts_plot <- renderPlotly({
      d <- object_counts_data()
      if (is.null(d) || nrow(d) == 0) return(plotly_islet_empty("No donors selected"))
      object_counts_plot(d, input$object_count_measure %||% "singles_per_mm2", donor_colors_reactive())
    })

    output$dl_object_counts <- downloadHandler(
      filename = function() paste0("endocrine_objects_per_donor_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        d <- object_counts_data()
        write_islet_csv(file, d, islet_csv_header(
          "Plot — endocrine objects per donor",
          paste("Donor status:", paste(validated_groups(), collapse = ", ")),
          "Islet: >= 2 whole cells and >= 35 um; Cluster: >= 2 cells, smaller; Single: one endocrine cell",
          "Singles are counted only; they are not analysis objects",
          "singles_excluded: the Single's strongest hormone is under twice the next; singles_not_analysed: it has no nucleus",
          "*_per_mm2: objects per mm2 of the section's tissue outline",
          paste("Rows:", if (is.null(d)) 0 else nrow(d), "donors")))
      }
    )

    output$dl_distribution_csv <- downloadHandler(
      filename = function() paste0("plot_right_", card_type("right"), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write_card_csv(file, "right")
    )

    # ===========================================================================
    # Selection description helper (also returned for Statistics download)
    # ===========================================================================

    get_selection_description <- function() {
      desc <- c()
      desc <- c(desc, paste("Generated:", Sys.time()))
      desc <- c(desc, paste("Data Layer:", validated_mode()))
      desc <- c(desc, paste("Selection:", validated_selection() %||% "Unknown"))
      desc <- c(desc, paste("Metric:", validated_metric()))

      desc <- c(desc, paste("Donor Groups:", paste(validated_groups(), collapse = ", ")))

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

      # Demographic + size filters. These CHANGE WHICH ISLETS ARE IN THE FILE,
      # so omitting them made two CSVs with different row sets look identical in
      # provenance. `min_cells`, age, and sex were all silently absent.
      oc <- prepared()$comp$object_class
      desc <- c(desc, paste("Endocrine objects:",
                            if (is.null(oc)) "Islet" else paste(sort(unique(as.character(oc))), collapse = ", ")))
      if (!is.null(input$age_range) && length(input$age_range) == 2) {
        desc <- c(desc, paste("Donor Age Range:",
                              input$age_range[1], "-", input$age_range[2], "years"))
      }
      if (!is.null(input$gender_filter)) {
        sex_labels <- c(F = "Female", M = "Male")
        sel <- input$gender_filter
        desc <- c(desc, paste("Sex:", paste(ifelse(sel %in% names(sex_labels),
                                                   sex_labels[sel], sel), collapse = ", ")))
      }

      desc <- c(desc, paste("Statistic:", input$stat %||% "mean_se"))
      desc <- c(desc, paste("Bin Width:", input$binwidth %||% "50", "\u00b5m"))
      if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
        desc <- c(desc, paste("Diameter Range:",
                              input$diam_range[1], "-", input$diam_range[2], "\u00b5m"))
      }

      # Analysis-rule state. Without these a file with 17 outliers removed is
      # indistinguishable from one with none.
      desc <- c(desc, paste("Zero-value exclusion (left plot):",
                            if (isTRUE(input$exclude_zero_top)) "ON" else "OFF"))
      desc <- c(desc, paste("Outlier filter:",
                            if (isTRUE(remove_outliers()))
                              paste0("ON (|z| > ", outlier_threshold(), " per donor status)")
                            else "OFF"))

      return(desc)
    }

    # ===========================================================================
    # Return values for other modules (Statistics, etc.)
    # ===========================================================================

    # `raw_df_prezero` is `raw_df_base` — the per-islet frame with every sidebar
    # filter applied but BEFORE `input$exclude_zero_top`. The Statistics tab
    # consumes this one so the Plot sidebar's zero checkbox cannot silently
    # alter statistical results; Statistics applies its own `stats_exclude_zero`
    # and reports it truthfully in the CSV header.
    list(
      raw_df                  = raw_df,
      raw_df_prezero          = raw_df_base,
      summary_df              = summary_df,
      get_selection_description = get_selection_description
    )
  })
}
