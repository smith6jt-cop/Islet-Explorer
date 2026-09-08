# Islet Explorer - Modular App Entrypoint
# All utility functions and module code live in R/ (auto-sourced by Shiny)
# Package loading and constants live in global.R (auto-sourced by Shiny)

# ---------- UI ----------

ui <- fluidPage(
  title = "Islet Explorer",
  lang = "en",
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
    tags$style(HTML("
    /* Viewer and trajectory mode styles - fix tab positioning */
    body.viewer-mode .col-sm-2 {
      display: none !important;
    }
    body.viewer-mode .col-sm-10 {
      width: 100% !important;
      max-width: 100% !important;
      flex: 0 0 100%;
    }
    body.trajectory-mode .container-fluid > .row > .col-sm-3 {
      display: none !important;
    }
    body.trajectory-mode .container-fluid > .row > .col-sm-9 {
      width: 100% !important;
      max-width: 100% !important;
      flex: 0 0 100%;
    }

    /* Global biomedical theme styling */
    body {
      background-color: #f8f9fa;
      font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
      padding-top: 0;
      min-width: 1200px;
      overflow-x: auto;
    }

    .skip-link {
      position: absolute;
      left: 12px;
      top: -80px;
      z-index: 10000;
      padding: 10px 14px;
      background: #ffffff;
      color: #1e3a72;
      border: 3px solid #1e3a72;
      border-radius: 4px;
      font-weight: 700;
    }

    .skip-link:focus {
      top: 12px;
    }

    a:focus-visible,
    button:focus-visible,
    input:focus-visible,
    select:focus-visible,
    [tabindex]:focus-visible {
      outline: 3px solid #ffbf47 !important;
      outline-offset: 2px !important;
    }

    .container-fluid {
      background-color: #ffffff;
      border-radius: 12px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.08);
      margin: 10px;
      padding: 20px;
      min-width: 1180px;
    }

    /* Logo header styling */
    .logo-header {
      display: flex;
      justify-content: flex-end;
      align-items: center;
      padding: 10px 0;
      margin-bottom: 10px;
      background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%);
      border-bottom: 2px solid #e3f2fd;
    }

    /* Enhanced card styling with biomedical color scheme */
    .card {
      background: linear-gradient(145deg, #ffffff 0%, #f8f9fa 100%);
      border: 1px solid #e3f2fd;
      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.08);
      transition: all 0.3s ease;
      border-radius: 12px;
    }
    .card:hover {
      box-shadow: 0 8px 20px rgba(44, 90, 160, 0.15);
      transform: translateY(-2px);
    }

    .card h5 {
      color: #2c5aa0;
      font-weight: 600;
      border-bottom: 2px solid #e3f2fd;
      padding-bottom: 8px;
      margin-bottom: 15px;
    }

    /* Sidebar styling with scientific theme */
    .sidebar {
      background: linear-gradient(180deg, #2c5aa0 0%, #1e3a72 100%);
      color: white;
      border-radius: 12px;
      padding: 20px;
      font-size: 14px;
      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.2);
    }

    .sidebar h4, .sidebar h5 {
      color: #ffffff;
      font-weight: 600;
      border-bottom: 1px solid rgba(255,255,255,0.2);
      padding-bottom: 8px;
    }

    .sidebar .form-group {
      margin-bottom: 18px;
    }

    .sidebar label {
      color: #e3f2fd;
      font-size: 13px;
      font-weight: 600;
    }

    .sidebar .form-control {
      background-color: rgba(255,255,255,0.9);
      border: 1px solid #b3d9ff;
      border-radius: 6px;
      color: #2c5aa0;
    }

    .sidebar .form-control:focus {
      background-color: #ffffff;
      border-color: #66b3ff;
      box-shadow: 0 0 0 0.2rem rgba(44, 90, 160, 0.25);
    }

    .sidebar .btn {
      background: linear-gradient(145deg, #1d70b8 0%, #1261a0 100%);
      border: none;
      border-radius: 6px;
      color: white;
      font-weight: 500;
    }

    .sidebar .btn:hover {
      background: linear-gradient(145deg, #1261a0 0%, #0b4778 100%);
      transform: translateY(-1px);
    }

    /* Tab styling */
    .nav-tabs {
      border-bottom: 2px solid #e3f2fd;
    }

    .nav-tabs .nav-link {
      color: #2c5aa0;
      font-weight: 500;
      border: none;
      border-radius: 8px 8px 0 0;
      margin-right: 4px;
      background-color: #f8f9fa;
    }

    .nav-tabs .nav-link.active {
      background: linear-gradient(145deg, #2c5aa0 0%, #1e3a72 100%);
      color: white;
      border-bottom: 3px solid #66b3ff;
    }

    .nav-tabs .nav-link:hover {
      background-color: #e3f2fd;
      color: #1e3a72;
    }

    /* Form controls styling */
    .form-control {
      border: 2px solid #e3f2fd;
      border-radius: 6px;
      transition: all 0.2s ease;
    }

    .form-control:focus {
      border-color: #66b3ff;
      box-shadow: 0 0 0 0.2rem rgba(44, 90, 160, 0.15);
    }

    /* Button styling */
    .btn-primary {
      background: linear-gradient(145deg, #2c5aa0 0%, #1e3a72 100%);
      border: none;
      border-radius: 8px;
      font-weight: 500;
      padding: 8px 16px;
      transition: all 0.2s ease;
    }

    .btn-primary:hover {
      background: linear-gradient(145deg, #1e3a72 0%, #0f1f3d 100%);
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.3);
    }

    /* Checkbox and radio button styling */
    input[type='checkbox'], input[type='radio'] {
      accent-color: #2c5aa0;
    }

    /* Slider styling */
    .irs--shiny {
      color: #2c5aa0;
    }

    .irs--shiny .irs-bar {
      background: linear-gradient(90deg, #66b3ff 0%, #2c5aa0 100%);
    }

    .irs--shiny .irs-handle {
      background: #2c5aa0;
      border: 3px solid #ffffff;
    }

    /* Keep the slider's small white value labels readable. */
    .irs--shiny .irs-from,
    .irs--shiny .irs-to,
    .irs--shiny .irs-single {
      background: #2c5aa0;
    }
    .irs--shiny .irs-from:before,
    .irs--shiny .irs-to:before,
    .irs--shiny .irs-single:before {
      border-top-color: #2c5aa0;
    }
    [id^='spatial-'][id$='_pheno_hint'] strong,
    [id^='spatial-show_phenotype_rules_inline'] .action-label {
      color: #315f8b !important;
    }

    /* Help text styling */
    .help-block {
      color: #b3d9ff;
      font-size: 12px;
      font-style: italic;
    }

    /* Well and panel styling */
    .well {
      background: linear-gradient(145deg, #f8f9fa 0%, #e3f2fd 100%);
      border: 1px solid #b3d9ff;
      border-radius: 8px;
    }
    .ai-chat-panel {
      background: linear-gradient(145deg, #ffffff 0%, #f8f9fa 100%);
      border: 1px solid #e3f2fd;
      border-radius: 12px;
      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.08);
      padding: 15px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .ai-chat-logo {
      display: flex;
      justify-content: center;
      align-items: center;
      padding-bottom: 8px;
      border-bottom: 1px solid rgba(44, 90, 160, 0.15);
    }
    .ai-chat-logo img {
      max-width: 100%;
      height: auto;
      max-height: 90px;
    }
    .ai-chat-header {
      display: flex;
      justify-content: center;
      align-items: center;
      margin-top: 8px;
      margin-bottom: 4px;
    }
    .ai-chat-history {
      flex: 1;
      overflow-y: auto;
      padding: 8px;
      display: flex;
      flex-direction: column;
      gap: 10px;
    }
    .ai-chat-message {
      padding: 10px 12px;
      border-radius: 8px;
      margin: 4px 0;
      max-width: 90%;
      word-wrap: break-word;
    }
    .ai-chat-message.user {
      background: linear-gradient(145deg, #e3f2fd 0%, #bbdefb 100%);
      border-left: 4px solid #2c5aa0;
      align-self: flex-end;
      margin-left: auto;
    }
    .ai-chat-message.assistant {
      background: linear-gradient(145deg, #f5f5f5 0%, #eeeeee 100%);
      border-left: 4px solid #66b3ff;
      align-self: flex-start;
      margin-right: auto;
    }
    .ai-chat-meta {
      font-weight: 600;
      font-size: 12px;
      display: block;
      margin-bottom: 4px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }
    .ai-chat-message.user .ai-chat-meta {
      color: #1e3a72;
    }
    .ai-chat-message.assistant .ai-chat-meta {
      color: #2c5aa0;
    }

    /* --- Plotly modebar positioning -------------------------------------
       The modebar (top-right hover controls) sits in the same vertical band
       as the centered chart title, overlapping it on hover. Reserve a 28px
       strip above each plotly figure so the modebar lives above the title.
       Applies to every plot_ly + ggplotly fig in the app. */
    .js-plotly-plot,
    .html-widget.plotly {
      padding-top: 28px;
    }
    .js-plotly-plot .plotly .modebar {
      top: -28px !important;
    }
  "))),
  tags$head(tags$style(HTML("
    .ai-chat-logo img {
      max-width: 100%;
      height: auto;
      max-height: 120px;
    }
    /* Equal height panels */
    .equal-height-row {
      display: flex;
      flex-wrap: nowrap;
    }
    .equal-height-panel {
      height: calc(100vh - 40px);
      overflow-y: auto;
    }
    /* AI chat panel - match card heights */
    .ai-chat-panel-container {
      height: calc(100vh - 100px);
      display: flex;
      flex-direction: column;
      overflow-y: auto;
      margin-top: 20px;
    }
  "))),
  uiOutput("theme_css"),
  tags$script(HTML("
    $(document).on('shiny:connected', function() {
      function markDecorativeIcons() {
        // These buttons already have text labels. Some Shiny/fontawesome
        // versions also label their presentation-only icons, which is invalid
        // ARIA and duplicates the button name in assistive technology.
        $('.shiny-download-link i[aria-label], .action-button i[aria-label]')
          .removeAttr('aria-label').attr('aria-hidden', 'true');
      }
      markDecorativeIcons();
      $(document).on('shiny:bound shiny:value', function() {
        window.requestAnimationFrame(markDecorativeIcons);
      });

      function updateTabAccessibility() {
        $('#tabs > li > a').each(function() {
          var active = $(this).parent().hasClass('active');
          $(this).attr({ 'aria-selected': String(active), tabindex: active ? '0' : '-1' });
        });
      }

      // Shiny supplies keyboard tab navigation. Keep its selected-state
      // semantics and a single Tab stop synchronized after activation.
      updateTabAccessibility();

      function adjustLayout() {
        var activeTab = $('#tabs li.active a').text().trim();
        var sidebar = $('.equal-height-panel').first();
        var mainPanel = $('.main-content-panel');

        // Show sidebar on Plot and Statistics tabs, hide on others
        if (activeTab === 'Plot' || activeTab === 'Statistics') {
          sidebar.show();
          mainPanel.css({
            'width': '88.33333333%',
            'max-width': '88.33333333%',
            'flex': '0 0 88.33333333%'
          });
        } else {
          sidebar.hide();
          mainPanel.css({
            'width': '100%',
            'max-width': '100%',
            'flex': '0 0 100%'
          });
        }
      }

      // Adjust on initial load
      setTimeout(adjustLayout, 100);

      // Adjust when tabs change
      $('a[data-toggle=\"tab\"]').on('shown.bs.tab', function() {
        adjustLayout();
        updateTabAccessibility();
      });
    });
  ")),
  tags$a(class = "skip-link", href = "#main-content", "Skip to main content"),
  tags$h1(class = "sr-only", "Islet Explorer"),
  fluidRow(class = "equal-height-row",
    # Left Sidebar Panel (Plot controls, hidden on other tabs via JS)
    column(width = 1.4, class = "equal-height-panel",
      role = "complementary", `aria-label` = "Plot controls",
      # 'Statistics' stays in this condition while that tab is hidden: with no
      # such tab mounted `input.tabs` can never take that value, so the branch is
      # inert, and keeping it means SHOW_STATISTICS_TAB=1 needs no second edit.
      conditionalPanel(
        condition = "input.tabs == 'Plot' || input.tabs == 'Statistics'",
        plot_sidebar_ui("plot")
      )
    ),
    # Main Panel
    column(width = 10.6, class = "equal-height-panel main-content-panel",
      id = "main-content", role = "main", tabindex = "-1",
      tabsetPanel(id = "tabs",
        tabPanel("Plot",
          fluidRow(
            plot_main_ui("plot", extra_panel = if (AI_ENABLED) column(2, ai_assistant_ui("ai")) else NULL)
          )
        ),
        tabPanel("Trajectory",
          trajectory_ui("traj")
        ),
        tabPanel("Viewer",
          viewer_ui("viewer")
        ),
        # Hidden unless SHOW_STATISTICS_TAB=1 (see 00_globals.R).
        # tabsetPanel drops NULL entries, so this mounts nothing when off.
        if (SHOW_STATISTICS_TAB) {
          tabPanel("Statistics",
            statistics_ui("stats")
          )
        },
        tabPanel("Spatial",
          spatial_ui("spatial")
        ),
        # Hidden unless SHOW_INFLAMMATION_TAB=1 (see 00_globals.R).
        # tabsetPanel drops NULL entries, so this mounts nothing when off.
        if (SHOW_INFLAMMATION_TAB) {
          tabPanel("Inflammation",
            inflammation_ui("inflam")
          )
        }
      )
    )
  ),
  tags$footer(
    class = "app-footer", role = "contentinfo",
    style = "padding: 10px 0 0; text-align: right;",
    tags$a(href = "https://accessibility.ufl.edu/", "UF Accessibility")
  )
)

# ---------- Server ----------

# Load and prepare the immutable summaries once per R worker, not once per
# browser session. In production this call also enforces the release manifest.
.shared_prepared_data <- tryCatch(
  load_shared_prepared_data(),
  error = function(e) {
    islet_log_event("release_data_validation_failed")
    stop("Application data validation failed.", call. = FALSE)
  }
)

server <- function(input, output, session) {

  # ---- Shared reactive state ----
  forced_image   <- reactiveVal(NULL)
  selected_islet <- reactiveVal(NULL)

  # Immutable process-wide object; reactive() preserves the existing module
  # interface without repeating I/O or prep_data() for each reviewer.
  prepared <- reactive(.shared_prepared_data)

  # ---- Theme CSS (light/dark toggle, if theme_bg input exists) ----
  output$theme_css <- renderUI({
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      tags$style(HTML("body { background-color: #000000; color: #e6e6e6; }
                       .well { background-color: #111111; }"))
    } else {
      tags$style(HTML("body { background-color: #ffffff; color: #111111; }"))
    }
  })

  # ---- Module servers ----

  # Active tab reactive (used to prevent duplicate non-namespaced output IDs)
  active_tab <- reactive(input$tabs)

  # Active phenotype palette (global, consumed by root renders + spatial module)
  # Two non-namespaced inputs update this: spatial_palette (Spatial tab) and drilldown_palette (drill-down panel)
  palette_name <- reactiveVal("High Contrast")

  observeEvent(input$spatial_palette, {
    if (!identical(input$spatial_palette, palette_name())) {
      palette_name(input$spatial_palette)
      updateSelectInput(session, "drilldown_palette", selected = input$spatial_palette)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$drilldown_palette, {
    if (!identical(input$drilldown_palette, palette_name())) {
      palette_name(input$drilldown_palette)
      updateSelectInput(session, "spatial_palette", selected = input$drilldown_palette)
    }
  }, ignoreInit = TRUE)

  active_palette <- reactive({
    get_phenotype_palette(palette_name())
  })

  # Active donor color palette (global, consumed by Plot, Trajectory, Statistics, Spatial)
  donor_palette_name <- reactiveVal("Bright")
  DONOR_PALETTE_CHOICES <- names(DONOR_PALETTES)

  # Sync all donor palette selectors (Plot sidebar, Trajectory, Spatial)
  observeEvent(input$sidebar_donor_palette, {
    if (!identical(input$sidebar_donor_palette, donor_palette_name())) {
      donor_palette_name(input$sidebar_donor_palette)
      updateSelectInput(session, "traj_donor_palette", selected = input$sidebar_donor_palette)
      updateSelectInput(session, "spatial_donor_palette", selected = input$sidebar_donor_palette)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$traj_donor_palette, {
    if (!identical(input$traj_donor_palette, donor_palette_name())) {
      donor_palette_name(input$traj_donor_palette)
      updateSelectInput(session, "sidebar_donor_palette", selected = input$traj_donor_palette)
      updateSelectInput(session, "spatial_donor_palette", selected = input$traj_donor_palette)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$spatial_donor_palette, {
    if (!identical(input$spatial_donor_palette, donor_palette_name())) {
      donor_palette_name(input$spatial_donor_palette)
      updateSelectInput(session, "sidebar_donor_palette", selected = input$spatial_donor_palette)
      updateSelectInput(session, "traj_donor_palette", selected = input$spatial_donor_palette)
    }
  }, ignoreInit = TRUE)

  active_donor_colors <- reactive({
    get_donor_color_palette(donor_palette_name())
  })

  # Global outlier filter (synced across Plot sidebar, Trajectory, Spatial).
  # Defaults to FALSE so the raw distribution renders out of the box; users
  # opt in via the checkbox when they want to exclude tails from trends and
  # tests. The three non-namespaced inputs `sidebar_remove_outliers`,
  # `spatial_remove_outliers` mirrors this reactive (the Trajectory tab no longer
  # exposes an outlier control, so it is not in the sync set);
  # observers keep them in sync. Modules read `remove_outliers_global()` and
  # the threshold reactive `outlier_z_global()` directly.
  remove_outliers_rv <- reactiveVal(FALSE)

  observeEvent(input$sidebar_remove_outliers, {
    if (!identical(input$sidebar_remove_outliers, remove_outliers_rv())) {
      remove_outliers_rv(isTRUE(input$sidebar_remove_outliers))
      updateCheckboxInput(session, "spatial_remove_outliers", value = isTRUE(input$sidebar_remove_outliers))
    }
  }, ignoreInit = TRUE)


  observeEvent(input$spatial_remove_outliers, {
    if (!identical(input$spatial_remove_outliers, remove_outliers_rv())) {
      remove_outliers_rv(isTRUE(input$spatial_remove_outliers))
      updateCheckboxInput(session, "sidebar_remove_outliers", value = isTRUE(input$spatial_remove_outliers))
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
      updateSliderInput(session, "spatial_outlier_z", value = val)
    }
  }, ignoreInit = TRUE)


  observeEvent(input$spatial_outlier_z, {
    val <- suppressWarnings(as.numeric(input$spatial_outlier_z))
    if (is.finite(val) && !isTRUE(all.equal(val, outlier_z_rv()))) {
      outlier_z_rv(val)
      updateSliderInput(session, "sidebar_outlier_z", value = val)
    }
  }, ignoreInit = TRUE)

  outlier_z_global <- reactive({
    val <- suppressWarnings(as.numeric(outlier_z_rv()))
    if (is.finite(val) && val >= 0.5) val else 3.0
  })

  # Plot module: returns list(raw_df, raw_df_prezero, summary_df, get_selection_description)
  plot_returns <- plot_server("plot", prepared, selected_islet, active_tab, active_donor_colors,
                              remove_outliers = remove_outliers_global,
                              outlier_threshold = outlier_z_global)

  # Trajectory module
  trajectory_server("traj", prepared, selected_islet, forced_image, active_tab, active_donor_colors,
                    remove_outliers = remove_outliers_global,
                    outlier_threshold = outlier_z_global)

  # Viewer module
  viewer_server("viewer", forced_image, reactive(input$tabs))

  # Statistics module: wired to Plot module's reactive outputs.
  # NOTE: `raw_df_prezero` (not `raw_df`) — the Statistics tab owns its own
  # zero-value exclusion via `stats_exclude_zero`. Passing the Plot tab's
  # post-zero-filter frame would let a control on another tab silently change
  # every test while the CSV header still reported "Zero-value filter: OFF".
  # `summary_df` is deliberately NOT passed: the Statistics AUC used to integrate
  # it, which tied that card to the Plot sidebar's bin-width and Mean/Median
  # controls (hidden while this tab is open) and made it ignore the Statistics
  # tab's own filters. AUC now derives from `raw_df_prezero` inside the module.
  # Gated on the same flag as its tabPanel (see 00_globals.R) so the module is
  # not wired up when the tab is hidden.
  if (SHOW_STATISTICS_TAB) {
    statistics_server("stats",
                      raw_df                  = plot_returns$raw_df_prezero,
                      get_selection_description = plot_returns$get_selection_description,
                      donor_colors             = active_donor_colors,
                      remove_outliers          = remove_outliers_global)
  }

  # Spatial Neighborhood module
  spatial_server("spatial", prepared, active_palette, active_donor_colors,
                 remove_outliers = remove_outliers_global,
                 outlier_threshold = outlier_z_global)

  # Inflammation module — cohort-level edgeR abundance + infiltration score.
  # Reads the cell-distance artifact through spatial_helpers.R (Parquet/DuckDB
  # in production; CSV only in development), not the Plot module's islet frame,
  # so no sidebar filters apply. The tab is intentionally absent from the sidebar-visibility
  # condition in the layout JS above, so the sidebar stays hidden here.
  # Gated on the same flag as its tabPanel (see 00_globals.R) so the module is
  # not wired up when the tab is hidden.
  if (SHOW_INFLAMMATION_TAB) {
    inflammation_server("inflam", prepared, active_donor_colors)
  }

  # AI Assistant module (self-contained; gated so the card is disabled AND
  # hidden by default -- see AI_ENABLED in 00_globals.R)
  if (AI_ENABLED) ai_assistant_server("ai")

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
        if (!islet_valid_scalar_choice(color_by, colnames(cells),
                                       "drilldown_color_invalid")) {
          color_by <- "phenotype"
        }
        show_peri_cells <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
        return(tryCatch(
          render_islet_drilldown_plot(info, cells, color_by = color_by, show_peri = show_peri_cells,
                                      show_peri_boundary = show_peri_bd, show_structures = show_structs,
                                      palette = active_palette()),
          error = function(e) {
            islet_log_event("drilldown_render_failed")
            ggplot2::ggplot() +
              ggplot2::annotate("text", x = 0.5, y = 0.5,
                                label = "Unable to render this view.",
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
        islet_log_event("segmentation_render_failed")
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "Unable to render this view.",
                            size = ISLET_ANNOTATE_SIZE, color = "red") +
          ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
      }
    )
  })
  output$islet_segmentation_view <- renderPlot({ islet_segmentation_view_plot() })
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
  output$islet_drilldown_summary <- renderPlot({ islet_drilldown_summary_plot() })
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

  # ---- Phenotype Rules modal (shared across all tabs) ----
  # Each tab's `phenotype_rules_button("show_phenotype_rules_<tab>")` fires one
  # of these observers; they all open the same modal so users see a single
  # source of truth for marker -> celltype gating.
  # `_inflam` stays listed even while the Inflammation tab is hidden: an
  # observeEvent on an input that never exists simply never fires, and keeping
  # it means flipping SHOW_INFLAMMATION_TAB back on needs no second edit here.
  for (id in c("show_phenotype_rules_plot",
               "show_phenotype_rules_traj",
               "show_phenotype_rules_spatial",
               "show_phenotype_rules_stats",
               "show_phenotype_rules_inflam",
               "show_phenotype_rules_viewer")) {
    local({
      input_id <- id
      observeEvent(input[[input_id]], { showModal(phenotype_rules_modal_ui()) },
                   ignoreInit = TRUE)
    })
  }
}

shinyApp(ui, server)
