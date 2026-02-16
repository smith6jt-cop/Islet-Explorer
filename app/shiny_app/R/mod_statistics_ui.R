statistics_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Statistical Analysis Controls"),
          fluidRow(
            column(3, selectInput(ns("alpha"), "Significance level (\u03b1):", choices = c("0.05","0.01","0.001"), selected = "0.05")),
            column(3, checkboxInput(ns("stats_remove_outliers"), "Remove outliers (>3 SD)", value = TRUE)),
            column(3, br(), actionButton(ns("run_tests"), "Run Statistics", class = "btn-primary")),
            column(3, br(), downloadButton(ns("download_stats"), "Download CSV"))
          )
        )
      )
    ),
    fluidRow(
      column(6,
        wellPanel(
          h4("Global ANOVA Results"),
          tableOutput(ns("stats_tbl"))
        )
      ),
      column(6,
        wellPanel(
          h4("Distribution by Donor Status"),
          plotlyOutput(ns("stats_plot"), height = 320)
        )
      )
    ),
    fluidRow(
      column(6,
        wellPanel(
          h4("Pairwise Comparisons"),
          plotlyOutput(ns("pairwise_plot"), height = 320)
        )
      ),
      column(6,
        wellPanel(
          h4("Area Under Curve by Donor Group"),
          plotlyOutput(ns("auc_plot"), height = 320),
          br(),
          tableOutput(ns("auc_table"))
        )
      )
    ),
    fluidRow(
      column(12,
        wellPanel(
          h4("Statistical Test Information"),
          uiOutput(ns("stats_explanation"))
        )
      )
    )
  )
}
