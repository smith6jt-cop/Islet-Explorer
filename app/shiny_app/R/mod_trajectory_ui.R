trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    checkboxInput(ns("show_outlier_table"), "Show outlier table (if outliers detected)", value = FALSE),
    uiOutput(ns("traj_status")),
    fluidRow(
      # Left card: Scatterplot and heatmap
      column(8,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          fluidRow(
            column(3,
              uiOutput(ns("traj_feature_selector")),
              selectInput(ns("traj_show_trend"), "Trend lines:",
                         choices = c("None" = "none",
                                     "Overall" = "overall",
                                     "By Donor Status" = "by_donor"),
                         selected = "by_donor")
            ),
            column(3,
              selectInput(ns("traj_color_by"), "Color points by:",
                         choices = c("Donor Status" = "donor_status",
                                     "Donor ID" = "donor_id"),
                         selected = "donor_status"),
              selectInput(ns("traj_point_size"), "Point size by:",
                         choices = c("Uniform" = "uniform",
                                     "Islet Diameter" = "islet_diam_um"),
                         selected = "islet_diam_um")
            ),
            column(3,
              sliderInput(ns("traj_alpha"), "Point transparency:",
                         min = 0.1, max = 1.0, value = 0.3, step = 0.05),
              sliderInput(ns("traj_point_size_slider"), "Point size:",
                         min = 0.5, max = 5.0, value = 2.3, step = 0.1)
            ),
            column(3,
              div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #f9f9f9; min-height: 100px;",
                h6("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; color: #333;"),
                uiOutput(ns("traj_legend"))
              )
            )
          ),
          plotlyOutput(ns("traj_scatter"), height = 600),
          br(),
          plotOutput(ns("traj_heatmap"), height = 120)
        )
      ),
      # Right card: UMAP plots
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px; height: 950px; overflow-y: auto;",
          fluidRow(column(12, h5("UMAP: Donor Status"), plotOutput(ns("traj_umap_donor"), height = 400))),
          br(),
          fluidRow(column(12, h5("UMAP: Selected Feature"), plotOutput(ns("traj_umap_feature"), height = 400)))
        )
      )
    ),
    fluidRow(column(12, uiOutput(ns("traj_outlier_info")))),
    uiOutput(ns("segmentation_viewer_panel"))
  )
}
