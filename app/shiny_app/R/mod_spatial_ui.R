# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# 5-card layout: Controls, Tissue Scatter + Leiden Panel, Enrichment + Heatmap

spatial_ui <- function(id) {
  ns <- NS(id)

  doc_style <- "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 8px 12px; font-size: 13px; color: #856404; margin-bottom: 10px;"

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
                           selected = "all", inline = TRUE)
            ),
            column(3,
              checkboxGroupInput(ns("groups"), "Donor Status",
                                 choices = c("ND", "Aab+", "T1D"),
                                 selected = c("ND", "Aab+", "T1D"), inline = TRUE)
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
    ),

    # ---- Cards 4 & 5: Enrichment + Phenotype Heatmap ----
    fluidRow(
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Immune Enrichment by Disease Stage"),
          div(style = doc_style,
            "These z-scores compare peri-islet immune cell proportions against ",
            "tissue-wide background using a Poisson model. Higher z indicates ",
            "local enrichment of that cell type near islets relative to the ",
            "tissue as a whole."
          ),
          plotlyOutput(ns("enrichment_plot"), height = "380px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          h5("Peri-Islet Phenotype Heatmap"),
          div(style = doc_style,
            "Shows mean proportion of each cell type in the 20\u00b5m peri-islet ",
            "expansion zone only, NOT core islet composition (shown in the Plot tab). ",
            "This reveals the cellular microenvironment surrounding each islet."
          ),
          plotlyOutput(ns("phenotype_heatmap"), height = "420px")
        )
      )
    )
  )
}
