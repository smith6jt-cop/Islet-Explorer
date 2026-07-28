#' Selection block for the frozen right-hand panel.
#'
#' The feature, the pseudotime definition and the colour encoding. Per-plot
#' aesthetics (point size, transparency, trend CI) stay on the card they affect.
#'
#' Gone with the cohort filter, which now owns them for every tab:
#' `traj_remove_outliers`, `traj_outlier_z` and `min_cells`. Trajectory's
#' min-cells copy bound to `total_cells`, which is core-only -- verified
#' 5214/5214 islets match core and 0/5214 match core+peri -- so the same slider
#' value gave 130 islets here and 309 everywhere else. One control, one
#' denominator.
trajectory_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("traj_feature_selector")),
    ie_field("Pseudotime mode",
      selectInput(ns("pt_mode"), NULL,
                  choices = c("Core only" = "core", "Core + Peri" = "combined"),
                  selected = "core", width = "100%")),
    ie_field("Pseudotime scale",
      selectInput(ns("pt_scale"), NULL,
                  choices = c("Rank-normalized" = "rank", "Raw DPT" = "raw"),
                  selected = "rank", width = "100%")),
    tags$small(style = "display: block; margin: -4px 0 10px; line-height: 1.35",
      "Rank spreads islets evenly; raw shows true DPT spacing, so where islets ",
      "actually concentrate is visible. Trend curvature on the rank axis is ",
      "partly an artifact of the transform."),
    ie_field("Colour points by",
      selectInput(ns("traj_color_by"), NULL,
                  choices = c("Donor Status" = "donor_status",
                              "Donor ID" = "donor_id",
                              "Leiden (res 0.3)" = "leiden_0.3",
                              "Leiden (res 0.5)" = "leiden_0.5",
                              "Leiden (res 0.8)" = "leiden_0.8",
                              "Leiden (res 1.0)" = "leiden_1.0"),
                  selected = "donor_status", width = "100%")),
    ie_field("Trend lines",
      selectInput(ns("traj_show_trend"), NULL,
                  choices = c("None" = "none", "Overall" = "overall",
                              "By Donor Status" = "by_donor"),
                  selected = "by_donor", width = "100%")),
    checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE)
  )
}

trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Hint strip moved to the rail's Help section.
    uiOutput(ns("traj_status")),
    fluidRow(
      # Left card: Scatterplot and heatmaps (narrowed from col-8)
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px",
          fluidRow(
            # The feature / pseudotime / colour selectors that used to fill the
            # first column here moved to trajectory_selection_ui(), in the frozen
            # right-hand panel. What remains is per-plot aesthetics, which stay
            # with the plot they affect.
            # 4+3+3 = 10 of 12 after the feature/pseudotime/colour selectors
            # moved to the frozen panel. Rebalanced to fill the row.
            column(5,
              selectInput(ns("traj_point_size"), "Point size by:",
                         choices = c("Cell Count" = "total_cells",
                                     "Islet Diameter" = "islet_diam_um",
                                     "Uniform" = "uniform"),
                         selected = "total_cells")
            ),
            column(3,
              sliderInput(ns("traj_alpha"), "Point transparency:",
                         min = 0, max = 0.95, value = 0.4, step = 0.05),
              sliderInput(ns("traj_point_size_slider"), "Point size:",
                         min = 0.5, max = 5.0, value = 3.0, step = 0.1),
              # min_cells moved to the shared cohort filter. It bound to
              # `total_cells` here (core-only) while every other tab used
              # core+peri, so one slider value meant two different things.
              checkboxInput(ns("traj_trend_ci"), "95% CI band on trends",
                           value = TRUE)
            ),
            column(4,
              div(class = "ie-legend-box",
                h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; font-size: 15px"),
                uiOutput(ns("traj_legend"))
              ),
              selectInput("traj_donor_palette", "Donor Status Colors",
                          choices = names(DONOR_PALETTES),
                          selected = "Bright")
            )
          ),
          uiOutput(ns("traj_outlier_badge")),
          plotlyOutput(ns("traj_scatter"), height = 650),
          br(),
          plotlyOutput(ns("traj_heatmap"), height = 100),
          hr(),
          h5("Multi-Feature Heatmap"),
          fluidRow(
            column(8, uiOutput(ns("multi_heatmap_marker_selector"))),
            column(4, sliderInput(ns("multi_heatmap_nbins"), "Number of bins:", min = 10, max = 40, value = 20, step = 5))
          ),
          plotOutput(ns("traj_multi_heatmap"), height = "auto"),
          div(style = "text-align: right; margin-top: 8px",
            downloadButton(ns("dl_multi_heatmap_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_multi_heatmap"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),
      # Right column: UMAPs stacked + segmentation viewer below
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px",
          fluidRow(
            column(6,
              h5("UMAP: Donor Status", style = "font-size: 13px; margin-top: 0"),
              plotOutput(ns("traj_umap_donor"), height = 500),
              div(style = "text-align: right; margin-top: 6px",
                downloadButton(ns("dl_umap_donor"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            ),
            column(6,
              h5("UMAP: Selected Feature", style = "font-size: 13px; margin-top: 0"),
              plotOutput(ns("traj_umap_feature"), height = 500),
              div(style = "text-align: right; margin-top: 6px",
                downloadButton(ns("dl_umap_csv"), "Download CSV",
                               class = "btn-sm btn-outline-secondary"),
                downloadButton(ns("dl_umap_feature"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          )
        ),
        uiOutput(ns("segmentation_viewer_panel"))
      )
    ),
    fluidRow(column(12, uiOutput(ns("traj_outlier_info"))))
  )
}
