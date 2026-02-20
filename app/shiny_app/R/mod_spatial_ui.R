# ---------- Spatial Neighborhood module UI ----------
# Exports: spatial_ui(id)
# 7-card layout for peri-islet neighborhood analysis

spatial_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ---- Card 1: Overview + Controls ----
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          fluidRow(
            column(8,
              h5("Spatial Neighborhood Analysis", style = "margin-top: 0;"),
              uiOutput(ns("overview_banner"))
            ),
            column(4, style = "text-align: right; padding-top: 5px;",
              downloadButton(ns("download_spatial"), "Download CSV")
            )
          ),
          hr(style = "margin: 10px 0;"),
          fluidRow(
            column(3,
              selectInput(ns("metric_category"), "Metric Category",
                          choices = c("Peri-Islet Composition", "Immune Infiltration",
                                      "Enrichment Z-scores", "Distance Metrics"))
            ),
            column(3,
              uiOutput(ns("feature_selector"))
            ),
            column(3,
              checkboxGroupInput(ns("groups"), "Donor Status",
                                 choices = c("ND", "Aab+", "T1D"),
                                 selected = c("ND", "Aab+", "T1D"), inline = TRUE)
            ),
            column(3,
              sliderInput(ns("diam_range"), "Diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10)
            )
          )
        )
      )
    ),

    # ---- Cards 2 & 3: Feature Distribution + Core vs Peri ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Feature Distribution by Disease Stage"),
          plotlyOutput(ns("distribution_plot"), height = "380px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Core vs Peri-Islet Comparison"),
          plotlyOutput(ns("core_peri_plot"), height = "380px")
        )
      )
    ),

    # ---- Cards 4 & 5: Peri-islet Heatmap + Immune Enrichment ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Peri-Islet Phenotype Heatmap"),
          plotlyOutput(ns("phenotype_heatmap"), height = "420px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Immune Enrichment by Disease Stage"),
          plotlyOutput(ns("enrichment_plot"), height = "420px")
        )
      )
    ),

    # ---- Cards 6 & 7: Pseudotime Correlation + Statistics ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Pseudotime Correlation"),
          plotlyOutput(ns("pseudotime_plot"), height = "380px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Statistical Tests"),
          tableOutput(ns("test_table")),
          hr(style = "margin: 8px 0;"),
          uiOutput(ns("methods_text"))
        )
      )
    )
  )
}
