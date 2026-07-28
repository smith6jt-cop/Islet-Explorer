# ---------- Plot module UI ----------
# Exports: plot_selection_ui(id), plot_main_ui(id)
# All Plot-related inputs are namespaced by this module.

#' Selection block for the frozen right-hand panel.
#'
#' Replaces plot_sidebar_ui(). What it kept is what SELECTS the data: the data
#' layer, the feature within it, the metric, and the summary statistic drawn on
#' the scatter. What it lost went to the shared cohort filter, which now owns
#' the same controls for every tab:
#'
#'   groups / age_range / gender_filter / aab_flags / min_cells / diam_range
#'   sidebar_remove_outliers / sidebar_outlier_z / sidebar_donor_palette
#'
#' `sidebar_donor_palette` moved to the rail footer with the phenotype palette.
#' The four conditionalPanels that used to hide half of these on the Statistics
#' tab are gone with them -- Statistics shares this selection block, because it
#' has always consumed the Plot module's frame.
plot_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    ie_field("Data layer",
      selectInput(ns("mode"), NULL,
                  choices = c("Microenvironment", "Cell Populations"),
                  selected = "Cell Populations", width = "100%")),
    uiOutput(ns("dynamic_selector")),
    uiOutput(ns("metric_selector")),
    # radioButtons carries its OWN label rather than being wrapped in ie_field.
    # Inside ie_field's flex column the Bootstrap-3-compat hanging indent
    # (label padding-left:20px with the input at margin-left:-20px) pushed the
    # radio circles outside the 236px panel, so the options rendered as bare
    # indented text with no visible control at all. The aside styles
    # .control-label to match .ie-field-label, so this looks identical.
    radioButtons(ns("stat"), "Summary statistic",
                 choices = c("Mean\u00b1SE" = "mean_se",
                             "Mean\u00b1SD" = "mean_sd",
                             "Median + IQR" = "median_iqr"),
                 selected = "mean_se"),
    uiOutput(ns("aab_warning")),
    ie_field("Diameter bin width",
      sliderInput(ns("binwidth"), NULL, min = 10, max = 100, value = 50, step = 5,
                  post = " \u00b5m", ticks = FALSE, width = "100%")),
    checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE),
    uiOutput(ns("plot_outlier_info")),
    downloadButton(ns("dl_summary"), "Download summary CSV", class = "btn-block")
  )
}

#' Main content UI for the Plot tab
#' Left card: scatter plot with controls. Right card: distribution comparison.
plot_main_ui <- function(id, extra_panel = NULL) {
  ns <- NS(id)
  tagList(
    # The three blue hint strips that sat here moved to the rail's
    # "Tips for this tab". They took a full row above the cards on every
    # visit, and the cards now start at the top of the panel.
    # Left card: Main scatter plot. Was col-6 with col-4 beside it and
    # col-2 held for the AI panel; with that gone the two cards split the
    # full width instead of leaving a sixth of the panel empty.
    column(7,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px;",
        h5("Islet Size Distribution", style = "margin-top: 0; color: #333;"),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("plt"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; border-top: 1px solid #eee; padding-top: 10px;",
          div(style = "display: flex; gap: 18px; align-items: stretch;",
            # Slider column 1
            div(style = "flex: 1; min-width: 0; display: flex; flex-direction: column;",
              # binwidth moved to the selection block; islet diameter range is
              # now the cohort filter's, shared with every other tab.
              selectInput(ns("plot_color_by"), "Color points by:",
                          choices = c("Donor Status" = "donor_status",
                                      "Donor ID" = "donor_id"),
                          selected = "donor_status")
            ),
            # Slider column 2
            div(style = "flex: 1; min-width: 0; display: flex; flex-direction: column;",
              sliderInput(ns("pt_size"), "Point size",
                          min = 0.3, max = 4.0, value = 3.0, step = 0.1),
              sliderInput(ns("pt_alpha"), "Point transparency",
                          min = 0, max = 0.95, value = 0.4, step = 0.05),
              radioButtons(ns("add_smooth"), "Trend line",
                           choices = c("None", "LOESS"),
                           selected = "None", inline = TRUE)
            ),
            # Checkbox column 3 \u2014 vertical stack spanning the full height of the sliders
            div(style = "flex: 0 0 auto; min-width: 165px; display: flex; flex-direction: column; justify-content: space-between; border-left: 1px solid #eee; padding-left: 16px;",
              checkboxInput(ns("exclude_zero_top"), "Exclude zero values", value = FALSE),
              checkboxInput(ns("show_points"), "Show individual points", value = FALSE),
              checkboxInput(ns("log_scale"), "Log scale y-axis", value = FALSE),
              checkboxInput(ns("log_scale_x"), "Log2 scale x-axis", value = FALSE)
            )
          )
        )
      )
    ),
    # Right card: Distribution comparison
    column(5,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px;",
        div(style = "display: flex; align-items: center; justify-content: space-between; margin-top: 0;",
          h5("Distribution Comparison", style = "margin: 0; color: #333;"),
          downloadButton(ns("dl_distribution_csv"), "CSV",
                         class = "btn-sm btn-outline-secondary",
                         style = "font-size: 11px; padding: 2px 8px;")
        ),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("dist"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; border-top: 1px solid #eee; padding-top: 10px;",
          uiOutput(ns("dist_ui"))
        )
      )
    ),
    # Extra panel (e.g., AI chat) injected between cards and full-width panels
    extra_panel,
    # Full-width card below both plot cards: segmentation/drill-down viewer
    column(12,
      uiOutput(ns("segmentation_viewer_panel"))
    )
  )
}
