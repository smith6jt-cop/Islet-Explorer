# ---------- Statistics module UI ----------
# Exports: statistics_ui(id)
# 5-card layout: Overview, Hypothesis Testing, Per-Bin Heatmap,
#                Trend Analysis, Demographics (conditional), Methods

statistics_ui <- function(id) {

  ns <- NS(id)

  tagList(
    # ---- Card 1: Overview Banner with inline controls ----
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          fluidRow(
            column(8,
              h5("Statistical Analysis", style = "margin-top: 0;"),
              uiOutput(ns("overview_banner"))
            ),
            column(4, style = "text-align: right; padding-top: 5px;",
              actionButton(ns("run_tests"), "Run Statistics",
                           class = "btn-primary", style = "margin-right: 8px;"),
              downloadButton(ns("download_stats"), "Download CSV",
                             style = "margin-top: 0;")
            )
          ),
          hr(style = "margin: 10px 0;"),
          fluidRow(
            column(2,
              radioButtons(ns("test_type"), "Test Type",
                           choices = c("Parametric", "Non-parametric"),
                           selected = "Parametric", inline = FALSE)
            ),
            column(2,
              selectInput(ns("alpha"), "Significance (\u03b1)",
                          choices = c("0.05", "0.01", "0.001"), selected = "0.05")
            ),
            column(2,
              checkboxInput(ns("stats_remove_outliers"), "Remove outliers (>3 SD)", value = TRUE)
            ),
            column(3,
              sliderInput(ns("bin_width"), "Per-bin width (\u00b5m)",
                          min = 1, max = 75, value = 50, step = 1)
            ),
            column(3,
              sliderInput(ns("stats_diam_range"), "Diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10)
            )
          )
        )
      )
    ),

    # ---- Cards 2 & 3: Hypothesis Testing + Per-Bin Heatmap ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Hypothesis Testing"),
          tableOutput(ns("test_results_table")),
          hr(style = "margin: 8px 0;"),
          h6("Effect Size Forest Plot", style = "color: #555;"),
          plotlyOutput(ns("forest_plot"), height = "320px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Per-Bin Significance Heatmap"),
          plotlyOutput(ns("bin_heatmap"), height = "360px")
        )
      )
    ),

    # ---- Cards 4 & 5: Trend Analysis + Demographics ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Trend Analysis (Kendall \u03c4)"),
          plotlyOutput(ns("trend_plot"), height = "360px")
        )
      ),
      column(6,
        uiOutput(ns("demographics_card"))
      )
    ),

    # ---- Card 6: AUC Analysis ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Area Under Curve by Donor Group"),
          plotlyOutput(ns("auc_plot"), height = "320px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("AUC Summary"),
          tableOutput(ns("auc_table")),
          uiOutput(ns("auc_interpretation"))
        )
      )
    ),

    # ---- Card 7: Methods & Interpretation ----
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Methods & Interpretation"),
          uiOutput(ns("stats_explanation"))
        )
      )
    )
  )
}
