# Islet Explorer - Modular App Entrypoint
# All utility functions and module code live in R/ (auto-sourced by Shiny)
# Package loading and constants live in global.R (auto-sourced by Shiny)

# ---------- UI ----------

ui <- bslib::page_fillable(
  padding = 0, gap = 0,

  # bs_theme(version = 5) UNMODIFIED, on purpose. bslib serves an unmodified
  # BS5 theme from a PRECOMPILED bundle (~0.15 s, no Sass run); any
  # customisation -- primary=, base_font=, enable-rounded= -- flips it to a
  # genuine compile, cold once per shiny-server worker. The entire design system
  # here is CSS custom properties, so bs_theme customisation would buy nothing
  # and cost the app's hard-won ~3.5 s startup.
  #
  # Bootstrap 5 is safe for this app's markup: shiny emits IDENTICAL classes
  # under 3 and 5 for selectInput / radioButtons / checkboxGroupInput /
  # actionButton / helpText / tabsetPanel (measured), and bslib's BS5 bundle
  # ships a BS3 compatibility layer that still styles .btn-default,
  # .checkbox-inline, .well, .form-group, .help-block, .control-label and
  # .irs--shiny. What is NOT safe is the old layout jQuery -- see below.
  theme = bslib::bs_theme(version = 5),

  useShinyjs(),
  tags$head(
    # Favicons for browser tab, bookmarks, and mobile pinning. Browsers pick
    # the most appropriate sized PNG based on context; the ICO is the legacy
    # default location for tools that hard-fetch /favicon.ico. Apple touch
    # icon covers iOS Safari Add-to-Home-Screen.
    tags$link(rel = "icon", type = "image/png", sizes = "16x16",  href = "favicon-16.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "32x32",  href = "favicon-32.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "48x48",  href = "favicon-48.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "192x192", href = "favicon-192.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "512x512", href = "favicon-512.png"),
    tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "apple-touch-icon.png"),
    tags$link(rel = "shortcut icon", type = "image/x-icon", href = "favicon.ico"),

    # Order is load-bearing: self-hosted face, then design-system tokens, then
    # the Shiny adapter. htmltools::renderDocument emits `result$head` AFTER
    # html-dependencies, so these already beat Bootstrap on source order --
    # which is why almost nothing in the adapter needs !important.
    tags$link(rel = "preload", as = "font", type = "font/woff2",
              href = "fonts/Archivo-Variable.woff2", crossorigin = "anonymous"),
    tags$link(rel = "stylesheet", href = "fonts/archivo.css"),
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$link(rel = "stylesheet", href = "modernist-shiny.css"),

    tags$style(HTML("
      /* Plotly modebar overlap fix, carried over verbatim from the pre-port
         sheet. The modebar sits in the same vertical band as the centred chart
         title and overlaps it on hover; reserve a 28px strip above every
         figure. The !important stays -- plotly sets `.modebar { top: 2px }` at
         equal specificity and recomputes it on resize. */
      .js-plotly-plot, .html-widget.plotly { padding-top: 28px; }
      .js-plotly-plot .plotly .modebar { top: -28px !important; }
    "))
  ),

  div(
    class = "ie-shell",

    railUI(
      "rail", IE_PIPELINE, IE_REFERENCE,
      brand = "Panc Floyd", logo = "logo.png",
      meta = "5,214 islets · 15 donors<br>nPOD CODEX · 31 markers",
      help = tagList(
        rail_help_link("help_rules", "Phenotype rules", "Marker definitions"),
        rail_help_link("help_tips",  "Tips for this tab", "How to read it")
      ),
      footer = tagList(
        # Real controls, not the delivered static text. These are the ONLY
        # writers of either palette; the sync-observer triads they replace are
        # removed in the cohort-bar phase.
        rail_palette_select("rail_phenotype_palette", "Phenotype palette",
                            names(PHENOTYPE_PALETTES), "High Contrast"),
        rail_palette_select("rail_donor_palette", "Donor palette",
                            names(DONOR_PALETTES), "Bright")
      )
    ),

    div(
      class = "ie-main",

      # The sidebar and the panels are SIBLINGS with independent scrollers, not
      # nested. That matches what .equal-height-panel gave both of them before
      # the port: the Plot sidebar is ~1720px of controls against a ~1000px
      # viewport, so if it shares the panels' scroller the whole page grows to
      # the sidebar's height and the shell stops containing anything.
      div(
        class = "ie-body",

        div(
          class = "ie-panels",
          # navset_hidden mounts EVERY panel and toggles CSS visibility, just
          # like tabsetPanel did (measured). That is why the drill-down
          # duplicate-output guards in mod_plot_server.R and
          # mod_trajectory_server.R are still required.
          bslib::navset_hidden(
            id = "panels",
            # `extra_panel` used to carry column(2, ai_assistant_ui("ai")). The
            # AI chat card was removed Jul 2026 at the user's request; the Plot
            # tab's two cards now take the full width.
            bslib::nav_panel_hidden("plot", fluidRow(plot_main_ui("plot"))),
            bslib::nav_panel_hidden("trajectory",   trajectory_ui("traj")),
            bslib::nav_panel_hidden("statistics",   statistics_ui("stats")),
            bslib::nav_panel_hidden("spatial",      spatial_ui("spatial")),
            bslib::nav_panel_hidden("inflammation", inflammation_ui("inflam")),
            bslib::nav_panel_hidden("viewer",       viewer_ui("viewer"))
          )
        ),

        # ---- Frozen control panel ----------------------------------------
        # One control surface for the whole app, in the same place on every
        # tab. Two blocks, and the split is the point: SELECTION is what each
        # analysis plots and is per-tab; FILTER is which islets are in scope and
        # is identical everywhere. Before this, donor status lived under three
        # different ids, min-cells under four with no syncing at all, and
        # Trajectory and Inflammation had no cohort filter whatsoever.
        #
        # Hidden on the Viewer, which is a read-only image surface that reads
        # none of this.
        conditionalPanel(
          condition = "input.panels != 'viewer'",
          tags$aside(
            class = "ie-aside",

            div(class = "ie-aside-block",
              div(class = "ie-aside-head", "Selection"),
              # Each tab mounts its own selectors here. They stay namespaced to
              # their module, so nothing about the modules' wiring changes --
              # only where the controls are rendered.
              conditionalPanel("input.panels == 'plot' || input.panels == 'statistics'",
                               plot_selection_ui("plot")),
              conditionalPanel("input.panels == 'trajectory'",
                               trajectory_selection_ui("traj")),
              conditionalPanel("input.panels == 'spatial'",
                               spatial_selection_ui("spatial")),
              conditionalPanel("input.panels == 'inflammation'",
                               inflammation_selection_ui("inflam"))
            ),

            div(class = "ie-aside-block",
              div(class = "ie-aside-head", "Cohort filter"),
              cohortFilterUI("cohort")
            ),

            div(class = "ie-aside-block ie-aside-foot",
              uiOutput("aside_export")
            )
          )
        )
      )
    )
  )
)

# Fail loudly at startup if a nav key is renamed without updating the guards in
# mod_plot_server.R / mod_trajectory_server.R / mod_viewer_server.R. A drifted
# key does not error at runtime -- it silently greys out a panel, which is
# exactly the failure the active-tab guard exists to prevent.
stopifnot(setequal(
  IE_TAB_KEYS,
  c("plot", "trajectory", "statistics", "spatial", "inflammation", "viewer")
))

# ---------- Server ----------

server <- function(input, output, session) {

  # ---- Shared reactive state ----
  forced_image   <- reactiveVal(NULL)
  selected_islet <- reactiveVal(NULL)

  # ---- Data loading chain ----
  # Priority: H5AD (islet_explorer.h5ad) > Excel (master_results.xlsx)
  validate_file <- reactive({
    has_h5ad <- !is.null(h5ad_path) && file.exists(h5ad_path)
    has_excel <- file.exists(master_path)
    shiny::validate(shiny::need(has_h5ad || has_excel,
                                paste("No data found. Checked:", h5ad_path, "and", master_path)))
    TRUE
  })

  master <- reactive({
    req(validate_file())
    load_master_auto(h5ad_path = h5ad_path, excel_path = master_path)
  })

  prepared <- reactive({
    pd <- prep_data(master())
    try({
      audit_na(pd$targets_all, "targets_all")
      audit_na(pd$markers_all, "markers_all")
      audit_na(pd$comp, "comp")
    }, silent = TRUE)
    pd
  })

  # NOTE: `output$theme_css` used to live here, mounted at the top of the UI. It
  # branched on `input$theme_bg`, an input that exists NOWHERE in this app -- so
  # the else branch always fired and injected
  #     body { background-color: #ffffff; color: #111111; }
  # from a renderUI in the BODY. htmltools puts head links before the body, so
  # at equal specificity that declaration won against every stylesheet: it has
  # been silently overriding the `background-color: #f8f9fa` declared a few
  # dozen lines above for as long as both have existed. Deleted rather than
  # repointed -- a real dark mode belongs in the CSS layer as a token swap on a
  # `data-theme` attribute, never as a body-level renderUI, which by
  # construction outranks the stylesheet that owns the ground colour.

  # ---- Shared cohort filter (frozen right-hand panel) ----
  # One filter for every tab. Before this, donor status lived under three ids,
  # min-cells under four with no syncing at all (and Trajectory's bound to a
  # different denominator), and Trajectory and Inflammation had no cohort filter
  # whatsoever.
  cohort <- cohortFilterServer("cohort", prepared)

  # ---- Module servers ----

  # Active tab reactive (used to prevent duplicate non-namespaced output IDs)
  # navset_hidden reports the nav_panel_hidden VALUE -- the lowercase key --
  # not a display label. See IE_PIPELINE in R/00_globals.R.
  active_tab <- reactive(input$panels)

  # ---- Rail drives the hidden tabset ----
  railServer("rail", "panels")

  # ---- Palette selectors ----
  # The rail footer is now the primary control for both palettes. The per-tab
  # selectors still exist and still work; this keeps every member of each group
  # in step until the cohort-bar phase deletes the per-tab copies, at which
  # point both of these collapse to a single `reactive(input$rail_*)`.
  #
  # Replaces five near-identical observer blocks that differed only in which ids
  # they mirrored to. Same behaviour, one definition.
  sync_selectors <- function(rv, ids) {
    lapply(ids, function(this) {
      observeEvent(input[[this]], {
        val <- input[[this]]
        if (!is.null(val) && !identical(val, rv())) {
          rv(val)
          for (other in setdiff(ids, this)) {
            updateSelectInput(session, other, selected = val)
          }
        }
      }, ignoreInit = TRUE)
    })
    invisible(NULL)
  }

  # Phenotype palette: rail footer + Spatial sidebar + the drill-down panel.
  palette_name <- reactiveVal("High Contrast")
  sync_selectors(palette_name,
                 c("rail_phenotype_palette", "spatial_palette", "drilldown_palette"))
  active_palette <- reactive(get_phenotype_palette(palette_name()))

  # Donor palette: rail footer + Plot sidebar + Trajectory + Spatial.
  donor_palette_name <- reactiveVal("Bright")
  sync_selectors(donor_palette_name,
                 c("rail_donor_palette", "sidebar_donor_palette",
                   "traj_donor_palette", "spatial_donor_palette"))
  active_donor_colors <- reactive(get_donor_color_palette(donor_palette_name()))

  # Global outlier filter (synced across Plot sidebar, Trajectory, Spatial).
  # Defaults to FALSE so the raw distribution renders out of the box; users
  # opt in via the checkbox when they want to exclude tails from trends and
  # tests. The three non-namespaced inputs `sidebar_remove_outliers`,
  # `traj_remove_outliers`, `spatial_remove_outliers` mirror this reactive;
  # observers keep them in sync. Modules read `remove_outliers_global()` and
  # the threshold reactive `outlier_z_global()` directly.
  remove_outliers_rv <- reactiveVal(FALSE)

  observeEvent(input$sidebar_remove_outliers, {
    if (!identical(input$sidebar_remove_outliers, remove_outliers_rv())) {
      remove_outliers_rv(isTRUE(input$sidebar_remove_outliers))
      updateCheckboxInput(session, "traj_remove_outliers", value = isTRUE(input$sidebar_remove_outliers))
      updateCheckboxInput(session, "spatial_remove_outliers", value = isTRUE(input$sidebar_remove_outliers))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$traj_remove_outliers, {
    if (!identical(input$traj_remove_outliers, remove_outliers_rv())) {
      remove_outliers_rv(isTRUE(input$traj_remove_outliers))
      updateCheckboxInput(session, "sidebar_remove_outliers", value = isTRUE(input$traj_remove_outliers))
      updateCheckboxInput(session, "spatial_remove_outliers", value = isTRUE(input$traj_remove_outliers))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$spatial_remove_outliers, {
    if (!identical(input$spatial_remove_outliers, remove_outliers_rv())) {
      remove_outliers_rv(isTRUE(input$spatial_remove_outliers))
      updateCheckboxInput(session, "sidebar_remove_outliers", value = isTRUE(input$spatial_remove_outliers))
      updateCheckboxInput(session, "traj_remove_outliers", value = isTRUE(input$spatial_remove_outliers))
    }
  }, ignoreInit = TRUE)

  remove_outliers_global <- reactive({
    isTRUE(remove_outliers_rv())
  })

  # Global outlier z-threshold slider (synced across Plot / Trajectory /
  # Spatial). Default 3.0 (legacy behavior). Range 0.5–10, step 0.5. Feeds
  # `flag_outliers(threshold = …)` for every z-score site. Spatial Card A's
  # IQR rule (1.5×) is intentionally NOT routed through this — enrichment z is
  # already a standardised quantity, so it uses a fixed Tukey fence on a
  # different scale.
  outlier_z_rv <- reactiveVal(3.0)

  observeEvent(input$sidebar_outlier_z, {
    val <- suppressWarnings(as.numeric(input$sidebar_outlier_z))
    if (is.finite(val) && !isTRUE(all.equal(val, outlier_z_rv()))) {
      outlier_z_rv(val)
      updateSliderInput(session, "traj_outlier_z", value = val)
      updateSliderInput(session, "spatial_outlier_z", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$traj_outlier_z, {
    val <- suppressWarnings(as.numeric(input$traj_outlier_z))
    if (is.finite(val) && !isTRUE(all.equal(val, outlier_z_rv()))) {
      outlier_z_rv(val)
      updateSliderInput(session, "sidebar_outlier_z", value = val)
      updateSliderInput(session, "spatial_outlier_z", value = val)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$spatial_outlier_z, {
    val <- suppressWarnings(as.numeric(input$spatial_outlier_z))
    if (is.finite(val) && !isTRUE(all.equal(val, outlier_z_rv()))) {
      outlier_z_rv(val)
      updateSliderInput(session, "sidebar_outlier_z", value = val)
      updateSliderInput(session, "traj_outlier_z", value = val)
    }
  }, ignoreInit = TRUE)

  outlier_z_global <- reactive({
    val <- suppressWarnings(as.numeric(outlier_z_rv()))
    if (is.finite(val) && val >= 0.5) val else 3.0
  })

  # Plot module: returns list(raw_df, raw_df_prezero, summary_df, get_selection_description)
  plot_returns <- plot_server("plot", prepared, selected_islet, active_tab, active_donor_colors,
                              remove_outliers = remove_outliers_global,
                              outlier_threshold = outlier_z_global,
                              cohort = cohort)

  # Trajectory module. `forced_image` used to be passed here as the 4th
  # positional argument and was never read in that module -- dropped, which
  # shifts `active_tab` and `active_donor_colors` up one position.
  trajectory_server("traj", prepared, selected_islet, active_tab, active_donor_colors,
                    remove_outliers = remove_outliers_global,
                    outlier_threshold = outlier_z_global)

  # Viewer module.
  # NOTE on `forced_image`: nothing in this app WRITES it. app.R creates the
  # reactiveVal, the Viewer reads it (mod_viewer_server.R:28), and no caller
  # ever sets it -- so the Viewer's forced-image branch is currently
  # unreachable. Kept because it is the correct hook for the Trajectory -> Viewer
  # hand-off, but it is a hook, not a live path. (Root CLAUDE.md described it as
  # "written by Trajectory, read by Viewer"; the first half was never true.)
  viewer_server("viewer", forced_image, active_tab)

  # Statistics module: wired to Plot module's reactive outputs.
  # NOTE: `raw_df_prezero` (not `raw_df`) — the Statistics tab owns its own
  # zero-value exclusion via `stats_exclude_zero`. Passing the Plot tab's
  # post-zero-filter frame would let a control on another tab silently change
  # every test while the CSV header still reported "Zero-value filter: OFF".
  # `summary_df` is deliberately NOT passed: the Statistics AUC used to integrate
  # it, which tied that card to the Plot sidebar's bin-width and Mean/Median
  # controls (hidden while this tab is open) and made it ignore the Statistics
  # tab's own filters. AUC now derives from `raw_df_prezero` inside the module.
  statistics_server("stats",
                    raw_df                  = plot_returns$raw_df_prezero,
                    get_selection_description = plot_returns$get_selection_description,
                    donor_colors             = active_donor_colors,
                    remove_outliers          = remove_outliers_global)

  # Spatial Neighborhood module
  spatial_server("spatial", prepared, active_palette, active_donor_colors,
                 remove_outliers = remove_outliers_global,
                 outlier_threshold = outlier_z_global)

  # Inflammation module — cohort-level edgeR abundance + infiltration score.
  # Reads data/cell_distance.csv directly (via the cached load_cell_distance()
  # in spatial_helpers.R), not the Plot module's islet frame, so no sidebar
  # filters apply. The tab is intentionally absent from the sidebar-visibility
  # condition in the layout JS above, so the sidebar stays hidden here.
  inflammation_server("inflam", prepared, active_donor_colors)

  # AI Assistant: removed from the UI Jul 2026 at the user's request. The module
  # files remain under R/ (unreferenced) and ai_helpers.R still loads the UF
  # Navigator credentials; nothing mounts them. Delete R/mod_ai_assistant_*.R
  # and R/ai_helpers.R to remove it entirely.

  # Footer of the frozen panel: the cohort readout is rendered by the module;
  # this slot carries the tab-appropriate export.
  output$aside_export <- renderUI({
    tab <- active_tab()
    tags$div(style = "margin-top: 10px; font-size: 11px; color: #777; line-height: 1.5;",
      "Filters above apply to every analysis tab. Cell-level tabs (Spatial, ",
      "Inflammation) receive donor scope only \u2014 min-cells and diameter are ",
      "islet-level and are reported as not applied rather than ignored.")
  })

  # ---- Unified segmentation + drilldown plot (root-level, used by both Plot and Trajectory panels) ----
  islet_segmentation_view_plot <- reactive({
    req(selected_islet())
    info <- selected_islet()
    show_peri_bd <- if (!is.null(input$drilldown_show_peri_boundary)) input$drilldown_show_peri_boundary else TRUE
    show_structs <- if (!is.null(input$drilldown_show_structures)) input$drilldown_show_structures else TRUE
    show_cells <- if (!is.null(input$drilldown_show_cells)) input$drilldown_show_cells else TRUE

    if (show_cells) {
      cells <- load_islet_cells(info$case_id, info$islet_key)
      if (!is.null(cells)) {
        color_by <- input$drilldown_color_by %||% "phenotype"
        show_peri_cells <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
        return(tryCatch(
          render_islet_drilldown_plot(info, cells, color_by = color_by, show_peri = show_peri_cells,
                                      show_peri_boundary = show_peri_bd, show_structures = show_structs,
                                      palette = active_palette()),
          error = function(e) {
            cat("[DRILLDOWN ERROR]", conditionMessage(e), "\n")
            ggplot2::ggplot() +
              ggplot2::annotate("text", x = 0.5, y = 0.5,
                                label = paste("Render error:", conditionMessage(e)),
                                size = ISLET_ANNOTATE_SIZE, color = "red") +
              ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
          }
        ))
      }
    }
    # Fallback: boundaries only (no cells or cells unavailable)
    tryCatch(
      render_islet_segmentation_plot(info, show_peri_boundary = show_peri_bd,
                                      show_structures = show_structs),
      error = function(e) {
        cat("[SEGMENTATION ERROR]", conditionMessage(e), "\n")
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("Render error:", conditionMessage(e)),
                            size = ISLET_ANNOTATE_SIZE, color = "red") +
          ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
      }
    )
  })
  output$islet_segmentation_view <- renderPlot({ islet_segmentation_view_plot() },
                                               bg = "transparent")
  output$dl_islet_segmentation_view <- downloadHandler(
    filename = function() {
      info <- selected_islet()
      stem <- if (!is.null(info)) paste0("islet_", info$case_id, "_", info$islet_key) else "islet_segmentation"
      islet_png_filename(stem)
    },
    content = function(file) {
      p <- islet_segmentation_view_plot(); if (is.null(p)) return()
      ggplot2::ggsave(file, plot = p,
                      width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                      dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
    }
  )

  islet_drilldown_summary_plot <- reactive({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells)
    show_peri <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
    if (!show_peri && "cell_region" %in% colnames(cells)) {
      cells <- cells[cells$cell_region == "core", , drop = FALSE]
    }
    tryCatch(render_drilldown_summary(cells, palette = active_palette()), error = function(e) NULL)
  })
  output$islet_drilldown_summary <- renderPlot({ islet_drilldown_summary_plot() },
                                               bg = "transparent")
  output$dl_islet_drilldown_summary <- downloadHandler(
    filename = function() {
      info <- selected_islet()
      stem <- if (!is.null(info)) paste0("islet_", info$case_id, "_", info$islet_key, "_composition") else "islet_composition"
      islet_png_filename(stem)
    },
    content = function(file) {
      p <- islet_drilldown_summary_plot(); if (is.null(p)) return()
      ggplot2::ggsave(file, plot = p,
                      width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                      dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
    }
  )

  # ---- Marker-rule legend for the phenotypes present in the selected islet ----
  output$islet_drilldown_rules <- renderUI({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells, "phenotype" %in% colnames(cells))
    present <- sort(unique(as.character(cells$phenotype)))
    if (length(present) == 0) return(NULL)
    pal <- active_palette()
    rows <- lapply(present, function(nm) {
      colour <- if (!is.null(pal) && nm %in% names(pal)) pal[[nm]] else "#CCCCCC"
      rule_html <- format_phenotype_rule(nm, "html")
      tags$div(
        style = "display: flex; align-items: flex-start; gap: 8px; padding: 4px 0; border-bottom: 1px solid #eef0f3;",
        tags$span(style = sprintf("display:inline-block; width:12px; height:12px; min-width:12px; margin-top:4px; background:%s; border:1px solid #888; border-radius:2px;", colour)),
        tags$div(
          tags$div(style = "font-weight:600; color:#333; font-size:13px;", nm),
          tags$div(style = "font-size:12px; color:#444; line-height:1.35;", rule_html)
        )
      )
    })
    tags$details(
      style = "margin-top: 10px; padding: 8px 10px; background:#fafbfd; border:1px solid #e2e8f0; border-radius: 4px; font-size: 13px;",
      tags$summary(style = "cursor: pointer; font-weight: 600; color: #4477AA;",
                   "Marker rules for phenotypes in this islet"),
      tags$div(style = "margin-top: 8px;", rows)
    )
  })

  output$islet_drilldown_table <- renderTable({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells)
    show_peri <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
    if (!show_peri && "cell_region" %in% colnames(cells)) {
      cells <- cells[cells$cell_region == "core", , drop = FALSE]
    }
    if ("cell_region" %in% colnames(cells)) {
      counts <- as.data.frame(table(Region = cells$cell_region), stringsAsFactors = FALSE)
      colnames(counts) <- c("Region", "Cells")
      counts$Region <- ifelse(counts$Region == "core", "Core", "Peri-islet")
      rbind(counts, data.frame(Region = "Total", Cells = sum(counts$Cells)))
    } else {
      data.frame(Region = "All", Cells = nrow(cells))
    }
  }, striped = TRUE, spacing = "s")

  # ---- Help (rail) ----
  # Was six `phenotype_rules_button("show_phenotype_rules_<tab>")` pills, one per
  # tab, plus a row of blue hint strips at the top of Plot and Trajectory. All of
  # it is in the rail now: one place, always reachable, and the first ~120px of
  # every panel goes back to the content.
  observeEvent(input$help_rules, showModal(phenotype_rules_modal_ui()), ignoreInit = TRUE)

  # Per-tab reading notes, previously a row of blue strips eating the first
  # ~120px of the Plot and Trajectory panels.
  observeEvent(input$help_tips, {
    tips <- switch(active_tab() %||% "plot",
      plot = tagList(
        tags$p(tags$strong("Click an islet."), " Turn on \u201cShow individual points\u201d, then click any point on the Islet Size Distribution scatter to open its segmentation map below."),
        tags$p(tags$strong("Isolate a donor."), " Set \u201cColour points by\u201d to Donor ID, then click a donor in the legend to show or hide it."),
        tags$p(tags$strong("Denominators are explicit."), " \u201c% of core cells\u201d, \u201ccells per islet\u201d and \u201ccells per mm\u00b2\u201d are different questions; the note under the metric selector names the denominator in force.")),
      trajectory = tagList(
        tags$p(tags$strong("Click an islet."), " Click any point on the pseudotime scatter to open its segmentation map."),
        tags$p(tags$strong("Rank vs raw DPT."), " Rank spreads islets evenly; raw shows true spacing, so where islets actually concentrate is visible. Trend curvature on the rank axis is partly an artifact of the transform.")),
      statistics = tagList(
        tags$p(tags$strong("One principled path."), " Test type is auto-selected from normality and \u03b1 is fixed at 0.05 \u2014 both were classic \u201cmove until significant\u201d levers. The sensitivity knobs sit behind a collapsed panel on purpose."),
        tags$p(tags$strong("Section 3 is exploratory."), " The per-bin analysis is hypothesis-generating, not confirmatory.")),
      spatial = tagList(
        tags$p(tags$strong("Drag to zoom."), " Double-click to reset. The tissue scatter draws every cell in the donor, so it is a ggplot rather than a plotly figure."),
        tags$p(tags$strong("Card B filters the summary, not the test."), " Unticking a phenotype changes what is drawn; the correction still runs across everything present.")),
      inflammation = tagList(
        tags$p(tags$strong("Cohort-level, not per-islet."), " Both cards read the tissue-wide single-cell table and test every cell type at once."),
        tags$p(tags$strong("Granularity matters."), " Four of the 21 labels are parents of other labels; the selector controls whether they are dropped, rolled up, or tested as given.")),
      viewer = tags$p("Click inside the viewer to enable pan and zoom. This tab is read-only \u2014 nothing here changes the data."),
      tags$p("No notes for this tab.")
    )
    showModal(modalDialog(title = "Tips for this tab", tips,
                          easyClose = TRUE, footer = modalButton("Close")))
  }, ignoreInit = TRUE)
}

shinyApp(ui, server)
