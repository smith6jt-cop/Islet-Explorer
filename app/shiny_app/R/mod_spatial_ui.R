# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4)

spatial_ui <- function(id) {
  ns <- NS(id)

  fluidRow(
    # ---- Left sidebar: Controls ----
    column(2,
      div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
        h5("Spatial Tissue Map", style = "margin-top: 0;"),
        uiOutput(ns("donor_selector")),
        hr(style = "margin: 8px 0;"),
        selectInput(ns("color_by"), "Color by",
                    choices = c("Phenotype" = "phenotype",
                                "Leiden Cluster" = "leiden"),
                    selected = "phenotype"),
        uiOutput(ns("leiden_res_selector")),
        hr(style = "margin: 8px 0;"),
        radioButtons(ns("region_filter"), "Show Cells",
                     choices = c("All" = "all",
                                 "Core + Peri" = "core_peri",
                                 "Core Only" = "core"),
                     selected = "all", inline = FALSE),
        checkboxInput(ns("color_background"), "Color background cells", value = FALSE),
        hr(style = "margin: 8px 0;"),
        checkboxGroupInput(ns("groups"), "Donor Status",
                           choices = c("ND", "Aab+", "T1D"),
                           selected = c("ND", "Aab+", "T1D"), inline = FALSE),
        hr(style = "margin: 8px 0;"),
        h5("Palettes", style = "font-size: 14px; margin-bottom: 8px;"),
        selectInput("spatial_palette", "Phenotype",
                    choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                    selected = "High Contrast"),
        selectInput("spatial_donor_palette", "Donor Status",
                    choices = names(DONOR_PALETTES),
                    selected = "Paul Tol (default)"),
        hr(style = "margin: 8px 0;"),
        downloadButton(ns("download_spatial"), "Download CSV", style = "width: 100%; font-size: 13px;")
      )
    ),

    # ---- Center: Tissue Scatter ----
    column(6,
      div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
        h5("Tissue Scatter Plot"),
        plotOutput(ns("tissue_scatter"), height = "800px")
      )
    ),

    # ---- Right: Leiden Panel ----
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
}
