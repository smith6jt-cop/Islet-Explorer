# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4)
# + 3 neighborhood analysis card rows below (A: Infiltration, B: Enrichment, C: Proximity)

#' Selection block for the frozen right-hand panel.
#'
#' These controls used to sit in a col-2 card at the left of the tab. Moving
#' them means Spatial's selectors appear in the same place as every other tab's.
#' Dropped on the way: `groups` (the cohort filter owns donor status now, for
#' all tabs) and the two palette selectInputs, which the rail footer owns.
spatial_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("donor_selector")),
    ie_field("Colour cells by",
      selectInput(ns("color_by"), NULL,
                  choices = c("Phenotype" = "phenotype", "Leiden Cluster" = "leiden"),
                  selected = "phenotype", width = "100%")),
    uiOutput(ns("leiden_res_selector")),
    ie_field("Show cells",
      radioButtons(ns("region_filter"), NULL,
                   choices = c("All" = "all", "Core + Peri" = "core_peri", "Core Only" = "core"),
                   selected = "all", inline = FALSE)),
    checkboxInput(ns("color_background"), "Colour background cells", value = FALSE),
    uiOutput(ns("phenotype_filter")),
    downloadButton(ns("download_spatial"), "Download CSV", class = "btn-block")
  )
}

spatial_ui <- function(id) {
  ns <- NS(id)

  # Helper: numbered section heading (matches statistics tab pattern)
  section_heading <- function(step, title, subtitle) {
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px",
      div(style = "display: flex; align-items: baseline; gap: 10px",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 28px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 18px")
      ),
      tags$small(subtitle, style = "display: block; margin-top: 4px")
    )
  }

  tagList(
    # ==== ROW 1 (existing): Controls + Scatter + Leiden ====
    fluidRow(style = "display: flex; flex-wrap: wrap",
      # ---- Center: Tissue Scatter ----
      # 6+4 = 10 of 12 once the col-2 control sidebar moved to the frozen
      # panel. Rebalanced so the tissue scatter and the UMAPs fill the row.
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%",
          div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 8px",
            h5("Tissue Scatter Plot", style = "margin: 0"),
            actionButton(ns("scatter_reset_zoom"), "Reset Zoom",
                         style = "font-size: 12px; padding: 2px 10px",
                         class = "btn-sm btn-outline-secondary")
          ),
          tags$small("Drag to zoom, double-click to reset.", style = "display: block; margin-bottom: 5px"),
          plotOutput(ns("tissue_scatter"), height = "700px",
                     dblclick = ns("scatter_dblclick"),
                     brush = brushOpts(id = ns("scatter_brush"), resetOnNew = TRUE)),
          div(style = "text-align: right; margin-top: 8px",
            downloadButton(ns("dl_tissue_scatter"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),

      # ---- Right: Leiden Panel ----
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%",
          h5("Islet-Level UMAPs"),
          uiOutput(ns("leiden_not_available")),
          fluidRow(
            column(6,
              h5("Leiden Clustering", style = "font-size: 13px; margin-top: 0"),
              plotlyOutput(ns("leiden_umap"), height = "400px")
            ),
            column(6,
              h5("Donor Status", style = "font-size: 13px; margin-top: 0"),
              plotOutput(ns("spatial_umap_donor"), height = 300),
              div(style = "text-align: right; margin-top: 6px",
                downloadButton(ns("dl_spatial_umap_donor"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          ),
          hr(style = "margin: 8px 0;"),
          h5("Cluster Composition", style = "font-size: 15px"),
          plotlyOutput(ns("cluster_composition"), height = "450px")
        )
      )
    ),

    # ==== ROW 2: Card A — Immune Infiltration Overview ====
    # The three neighborhood cards moved to their own tab -- see
    # spatial_neighborhood_ui(). One tab was carrying a tissue scatter, two
    # UMAPs, a composition chart AND three analysis cards, which is more than
    # one panel can hold at any width.
  )
}

#' Neighborhood analysis: the three cards that used to sit below the tissue map.
#'
#' Same module namespace as spatial_ui(), so every server output and input id is
#' unchanged -- this is a relocation, not a rewrite. `neighborhood_cards` is a
#' single renderUI holding all three cards, so it mounts here intact.
spatial_neighborhood_ui <- function(id) {
  ns <- NS(id)
  tagList(uiOutput(ns("neighborhood_cards")))
}
