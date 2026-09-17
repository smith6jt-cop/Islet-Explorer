trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("traj_status")),
    p("Fixed islet ordering with arbitrary orientation. Pseudotime is not elapsed disease time. Display filters never refit or rescale the coordinates."),
    fluidRow(
      column(3, selectInput(ns("traj_donor"), "Donor", choices = c("All donors" = "all")),
             selectInput(ns("traj_ratio_type"), "Selected cell type", choices = phenotype_keep(names(PHENOTYPE_COLORS)), selected = "Beta cell")),
      column(3, selectInput(ns("traj_ratio_region"), "Count cells in", choices = c("Core + peri-islet" = "both", "Core" = "core", "Peri-islet" = "peri"), selected = "both"),
             selectInput(ns("pt_scale"), "Pseudotime scale", choices = c("Reference-cohort percentile" = "rank", "Raw pseudotime" = "raw"))),
      column(3, sliderInput(ns("cell_range"), "Measured core cells (display filter)", min = 1, max = 3000, value = c(10,3000)),
             selectInput(ns("traj_show_trend"), "Descriptive trend", choices = c("Overall" = "overall", "By disease status" = "by_donor", "None" = "none"))),
      column(3, selectInput(ns("traj_color_by"), "Color points by", choices = c("Disease status" = "donor_status", "Donor" = "case_id")),
             downloadButton(ns("dl_trajectory_csv"), "Download ratios and coordinates"))
    ),
    uiOutput(ns("traj_y_note")), uiOutput(ns("traj_scatter_panel")),
    p("Curves are descriptive. No donor-population confidence intervals or significance tests are shown here."),
    h4("Fixed diffusion graph"),
    p("The full primary-component graph remains faintly visible; selected islets are highlighted. Other components have no comparable primary coordinates."),
    uiOutput(ns("traj_graph_note")),
    plotlyOutput(ns("paga_layout"), height = 500),
    downloadButton(ns("dl_paga_csv"), "Download selected graph coordinates"),
    h4("Marker heatmap"),
    uiOutput(ns("traj_marker_controls")),
    uiOutput(ns("traj_multi_heatmap_panel")),
    downloadButton(ns("dl_multi_heatmap_csv"), "Download selected marker measurements"),
    uiOutput(ns("segmentation_viewer_panel"))
  )
}
