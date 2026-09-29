trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style(HTML(".trajectory-controls .form-group { margin-bottom: 6px; } .trajectory-controls h5 { margin-bottom: 6px; }")),
    tab_help("Trajectory help", "show_phenotype_rules_traj",
      "Click any point on the pseudotime scatter plot to view the segmentation map for that islet."),
    fluidRow(
      # Left card: Scatterplot and heatmaps (narrowed from col-8)
      column(7,
        div(class = "card trajectory-controls", style = "padding: 15px; margin-bottom: 20px;",
          h5("Data controls", style = "margin-top: 0;"),
          fluidRow(
            column(8,
              sliderInput(ns("cell_range"), "Cells per islet:",
                          min = 1, max = traj_max_islet_cells(),
                          value = c(TRAJ_DEFAULT_MIN_CELLS, traj_max_islet_cells()),
                          step = 5),
              uiOutput(ns("traj_status"))
            ),
            column(4,
              actionButton(ns("recalculate_trajectory"), "Recalculate trajectory",
                           class = "btn-primary", style = "margin-top: 25px;"),
              uiOutput(ns("traj_range_pending"))
            )
          ),
          fluidRow(
            column(4,
              # Y-axis measure. Marker expression is a scVI-denoised mean from
              # .X; the cell-type ratio is Barlow et al.'s Fig 2F y-axis,
              # log(# cells / # endocrine cells), which removes the islet-size
              # confound that a raw count or proportion carries.
              selectInput(ns("traj_y_mode"), "Y axis:",
                          choices = c("Marker expression" = "marker",
                                      "Cell-type ratio (log cells / endocrine cell)" = "ratio"),
                          selected = "ratio")
            ),
            column(4,
              conditionalPanel(
                condition = "input.traj_y_mode == 'ratio'", ns = ns,
                selectInput(ns("traj_ratio_region"), "Count cells in:",
                            choices = BARLOW_RATIO_REGIONS, selected = "core")
              )
            ),
            column(4,
              uiOutput(ns("traj_feature_selector"))
            )
          ),
          fluidRow(
            column(4,
              selectInput(ns("traj_show_trend"), "Trend lines:",
                         choices = c("None" = "none",
                                     "Overall" = "overall",
                                     "By Donor Status" = "by_donor"),
                         selected = "by_donor")
            ),
            column(4,
              selectInput(ns("pt_scale"), "Pseudotime scale:",
                         choices = c("Rank-normalized" = "rank",
                                     "Raw DPT" = "raw"),
                         selected = "rank")
            ),
            column(4,
              checkboxInput(ns("traj_trend_ci"), "95% CI band on trends",
                            value = TRUE)
            )
          ),
          tags$details(style = "border-top: 1px solid #ddd; padding-top: 8px; margin: 8px 0 12px;",
            tags$summary("View controls", style = "cursor: pointer; font-weight: 600; color: #333;"),
            fluidRow(
              column(4,
                selectInput(ns("traj_color_by"), "Color points by:",
                           choices = c("Donor Status" = "donor_status",
                                       "Donor ID" = "donor_id",
                                       "Leiden (res 0.3)" = "leiden_0.3",
                                       "Leiden (res 0.5)" = "leiden_0.5",
                                       "Leiden (res 0.8)" = "leiden_0.8",
                                       "Leiden (res 1.0)" = "leiden_1.0"),
                           selected = "donor_status"),
                selectInput("traj_donor_palette", "Donor Status Colors",
                            choices = names(DONOR_PALETTES),
                            selected = "Bright")
              ),
              column(4,
                selectInput(ns("traj_point_size"), "Point size by:",
                           choices = c("Cell Count" = "total_cells",
                                       "Islet Diameter" = "islet_diam_um",
                                       "Uniform" = "uniform"),
                           selected = "total_cells"),
                sliderInput(ns("traj_point_size_slider"), "Point size:",
                            min = 0.5, max = 5.0, value = 1.5, step = 0.1)
              ),
              column(4,
                sliderInput(ns("traj_alpha"), "Point transparency:",
                            min = 0, max = 0.95, value = 0.4, step = 0.05),
                div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #f9f9f9; min-height: 100px;",
                  h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; color: #333; font-size: 15px;"),
                  uiOutput(ns("traj_legend"))
                )
              )
            )
          ),
          plotlyOutput(ns("traj_scatter"), height = 850),
          plotlyOutput(ns("traj_heatmap"), height = 100)
        )
      ),
      # Right column: stacked UMAPs and a persistent segmentation card.
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          h5("UMAP: Donor Status", style = "font-size: 13px; margin-top: 0;"),
          plotOutput(ns("traj_umap_donor"), height = 320),
          div(style = "text-align: right; margin-top: 6px;",
            downloadButton(ns("dl_umap_donor"), "Download PNG",
                           class = "btn-sm btn-outline-primary")),
          hr(style = "margin: 12px 0;"),
          h5("UMAP: Selected Feature", style = "font-size: 13px; margin-top: 0;"),
          plotOutput(ns("traj_umap_feature"), height = 320),
          div(style = "text-align: right; margin-top: 6px;",
            downloadButton(ns("dl_umap_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_umap_feature"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        ),
        div(class = "card", style = "padding: 15px; margin-bottom: 20px; min-height: 650px;",
          h5("Segmentation map", style = "margin-top: 0;"),
          uiOutput(ns("segmentation_viewer_panel"))
        )
      )
    ),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          h5("Multi-Feature Heatmap", style = "color: #2c5aa0; margin-top: 0;"),
          fluidRow(
            column(8, uiOutput(ns("multi_heatmap_marker_selector"))),
            column(4, sliderInput(ns("multi_heatmap_nbins"), "Number of bins:", min = 10, max = 40, value = 20, step = 5))
          ),
          plotOutput(ns("traj_multi_heatmap"), height = "auto"),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("dl_multi_heatmap_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_multi_heatmap"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      )
    ),
    # ---- PAGA force-directed layout (Barlow et al. Fig 2B/2E) ----------
    # The feature-vs-pseudotime scatter above cannot be de-overplotted by any
    # axis transform: the pile-up is at EXACT zero (95.6% of islets have no CD8a
    # T cell, 98.2% no B cell) and log(0) is still 0. This 2-D graph embedding
    # is the plot that actually separates islets -- 8.1% overplotted vs 84-96%.
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          h5("PAGA force-directed layout (Barlow et al. 2024)",
             style = "color: #2c5aa0; margin-top: 0;"),
          fluidRow(
            column(3,
              selectInput(ns("paga_color_by"), "Color islets by:",
                          choices = c("Pseudotime" = "pseudotime",
                                      "Donor status" = "donor_status",
                                      "Leiden cluster" = "leiden",
                                      "Donor" = "case_id",
                                      "Selected feature" = "feature"),
                          selected = "pseudotime")
            ),
            column(3,
              sliderInput(ns("paga_point_size"), "Point size:",
                          min = 0.3, max = 4.0, value = 1.6, step = 0.1)
            ),
            column(6, uiOutput(ns("paga_note")))
          ),
          plotlyOutput(ns("paga_layout"), height = 620),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("dl_paga_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"))
        )
      )
    )
  )
}
