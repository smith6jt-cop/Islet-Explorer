# ---------- Statistics module UI ----------
# Exports: statistics_ui(id)
# 5-section narrative layout: Configure, Primary Results, Size-Dependent,
#                             Confounders & Deeper, Methods Reference

statistics_ui <- function(id) {

  ns <- NS(id)

  # Helper: numbered section heading with subtitle

  section_heading <- function(step, title, subtitle) {
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
      div(style = "display: flex; align-items: baseline; gap: 10px;",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 24px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 17px; color: #333;")
      ),
      p(subtitle, style = "margin: 4px 0 0 0; color: #777; font-size: 13px;")
    )
  }

  tagList(

    # ===== SECTION 1: Configure Analysis =====
    section_heading("1", "Configure Analysis",
                    "Select your test parameters, then click Run Statistics."),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          fluidRow(
            column(8,
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
                           selected = "Parametric", inline = FALSE),
              tags$small(style = "color: #888; display: block; margin-top: -8px; line-height: 1.3;",
                tags$strong("Parametric"), " (ANOVA / t-test): assumes roughly normal data within groups.",
                tags$br(),
                tags$strong("Non-parametric"), " (Kruskal-Wallis / Wilcoxon): no normality assumption; based on ranks."
              )
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

    # ===== SECTION 2: Primary Results =====
    section_heading("2", "Primary Results",
                    "Global test across all three donor groups, plus pairwise comparisons with effect sizes."),

    fluidRow(style = "display: flex; flex-wrap: wrap;",
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Hypothesis Testing"),
          tableOutput(ns("test_results_table"))
        )
      ),
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Effect Size Forest Plot"),
          plotlyOutput(ns("forest_plot"), height = "350px")
        )
      )
    ),

    # ===== SECTION 3: Size-Dependent Patterns =====
    section_heading("3", "Size-Dependent Patterns",
                    "Does the group effect depend on islet size? Each diameter bin is tested independently to detect effect modification."),

    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Stratified Tests by Islet Diameter"),
          tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
            "ANOVA and Kendall \u03c4 computed within each size bin. ",
            "P-values are BH-corrected across bins to control false discovery rate. ",
            "Identifies which size ranges drive or lack group differences."),
          plotlyOutput(ns("bin_heatmap"), height = "320px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Trend Analysis (Kendall \u03c4)"),
          tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
            "Direction and strength of the ND \u2192 Aab+ \u2192 T1D gradient within each size bin. ",
            "Positive \u03c4 = increases with disease progression."),
          plotlyOutput(ns("trend_plot"), height = "320px")
        )
      )
    ),

    # ===== SECTION 4: Confounders & Deeper Analysis =====
    section_heading("4", "Confounders & Deeper Analysis",
                    "Check whether demographics explain the effect, and compare integrated area under the curve across groups."),

    fluidRow(
      column(6,
        uiOutput(ns("demographics_card"))
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Area Under Curve by Donor Group"),
          plotlyOutput(ns("auc_plot"), height = "280px"),
          hr(style = "margin: 10px 0;"),
          tableOutput(ns("auc_table")),
          uiOutput(ns("auc_interpretation"))
        )
      )
    ),

    # ===== SECTION 5: Methods Reference =====
    section_heading("5", "Methods Reference",
                    "Technical details about tests, corrections, and guidelines."),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
          uiOutput(ns("stats_explanation"))
        )
      )
    )
  )
}
