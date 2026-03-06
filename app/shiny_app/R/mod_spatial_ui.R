# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# 3-card layout: Controls, Tissue Scatter + Leiden Panel

spatial_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ---- Card 1: Controls (full-width) ----
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          fluidRow(
            column(6,
              h5("Spatial Tissue Map", style = "margin-top: 0;"),
              uiOutput(ns("overview_banner"))
            ),
            column(6, style = "text-align: right; padding-top: 5px;",
              downloadButton(ns("download_spatial"), "Download Neighborhood CSV")
            )
          ),
          hr(style = "margin: 10px 0;"),
          fluidRow(
            column(3,
              uiOutput(ns("donor_selector"))
            ),
            column(2,
              selectInput(ns("color_by"), "Color by",
                          choices = c("Phenotype" = "phenotype",
                                      "Leiden Cluster" = "leiden"),
                          selected = "phenotype")
            ),
            column(2,
              uiOutput(ns("leiden_res_selector"))
            ),
            column(2,
              radioButtons(ns("region_filter"), "Show Cells",
                           choices = c("All" = "all",
                                       "Core + Peri" = "core_peri",
                                       "Core Only" = "core"),
                           selected = "all", inline = TRUE),
              checkboxInput(ns("color_background"), "Color background cells", value = FALSE)
            ),
            column(3,
              checkboxGroupInput(ns("groups"), "Donor Status",
                                 choices = c("ND", "Aab+", "T1D"),
                                 selected = c("ND", "Aab+", "T1D"), inline = TRUE)
            )
          ),
          fluidRow(style = "margin-top: 4px;",
            column(3,
              # Non-namespaced: global phenotype palette synced via app.R observers
              selectInput("spatial_palette", "Phenotype Palette",
                          choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                          selected = "High Contrast")
            ),
            column(3,
              selectInput("spatial_donor_palette", "Donor Status Colors",
                          choices = names(DONOR_PALETTES),
                          selected = "Paul Tol (default)")
            )
          )
        )
      )
    ),

    # ---- Cards 2 & 3: Tissue Scatter (col-8) + Leiden Panel (col-4) ----
    fluidRow(
      column(8,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          h5("Tissue Scatter Plot"),
          plotOutput(ns("tissue_scatter"), height = "800px")
        )
      ),
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Leiden Clustering (Islet-Level)"),
          uiOutput(ns("leiden_not_available")),
          plotlyOutput(ns("leiden_umap"), height = "380px"),
          hr(style = "margin: 8px 0;"),
          h5("Cluster Composition", style = "font-size: 15px;"),
          plotlyOutput(ns("cluster_composition"), height = "350px")
        )
      )
    )
  )
}
