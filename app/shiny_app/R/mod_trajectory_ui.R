trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("traj_status")),
    uiOutput(ns("boundary_review")),
    p("Fixed islet ordering with arbitrary orientation. Pseudotime is not elapsed disease time. Display filters never refit or rescale the coordinates. The reference rank spaces supported islets evenly from 0 to 1; exact raw-coordinate ties are ordered by an islet-ID hash and have no biological order within the tie."),
    div(class = "card trajectory-controls-card", style = "padding: 15px; margin-bottom: 18px;",
    fluidRow(
      column(3, selectInput(ns("traj_donor"), "Donor", choices = c("All donors" = "all")),
             selectInput(ns("traj_ratio_type"), "Selected cell type", choices = islet_phenotype_choices(), selected = "Beta cell")),
      column(3, selectInput(ns("traj_ratio_region"), "Count cells in", choices = c("Core" = "core", "Core + peri-islet" = "both", "Peri-islet" = "peri"), selected = "core"),
             selectInput(ns("pt_scale"), "Pseudotime scale", choices = c("Reference-cohort rank (evenly spaced)" = "rank", "Raw pseudotime" = "raw"))),
      column(3, sliderInput(ns("cell_range"), "Measured core cells (display filter)", min = 0, max = 3000, value = c(0,3000)),
             selectInput(ns("traj_show_trend"), "Descriptive trend", choices = c("Overall" = "overall", "By disease status" = "by_donor", "None" = "none"))),
      column(3, selectInput(ns("traj_color_by"), "Color points by", choices = c("Disease status" = "donor_status", "Donor" = "case_id")),
             selectInput(ns("traj_donor_palette"), "Disease-status palette", choices = names(DONOR_PALETTES), selected = "Bright"),
             downloadButton(ns("dl_trajectory_csv"), "Download ratios and coordinates"))),
    fluidRow(
      column(3, selectInput(ns("traj_point_size"), "Point size by", choices = c("Uniform" = "uniform", "Measured core cells" = "measured_cells", "Section area" = "area_um2"), selected = "uniform")),
      column(3, sliderInput(ns("traj_point_size_slider"), "Point size", min = .5, max = 6, value = 3, step = .1)),
      column(3, sliderInput(ns("traj_alpha"), "Point transparency", min = 0.1, max = 1, value = .6, step = .05)),
      column(3, checkboxInput(ns("traj_trend_ci"), "Show trend confidence band", value = FALSE))
    )),
    uiOutput(ns("traj_y_note")), uiOutput(ns("traj_scatter_panel")),
    uiOutput(ns("traj_strip_panel")),
    p("Curves are descriptive. No donor-population confidence intervals or significance tests are shown here."),
    h4("Fixed trajectory graph"),
    p("The full saved graph remains faintly visible; selected islets are highlighted. Disconnected components occupy separate graph layouts and have no comparable primary pseudotime coordinates."),
    uiOutput(ns("traj_graph_note")),
    uiOutput(ns("traj_graph_panel")),
    downloadButton(ns("dl_paga_csv"), "Download selected graph coordinates"),
    h4("Marker heatmap"),
    uiOutput(ns("heatmap_axis_note")),
    uiOutput(ns("traj_marker_controls")),
    uiOutput(ns("traj_multi_heatmap_panel")),
    downloadButton(ns("dl_multi_heatmap_csv"), "Download selected marker measurements"),
    h4("Inspect an islet"),
    selectizeInput(ns("inspect_islet_id"), "Islet (current donor and cell-range selection)", choices = NULL,
                   options = list(placeholder = "Type an islet name to search")),
    actionButton(ns("inspect_islet"), "Show outline and cells"),
    uiOutput(ns("segmentation_viewer_panel"))
  )
}
