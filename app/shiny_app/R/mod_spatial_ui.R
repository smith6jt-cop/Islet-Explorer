# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Data / View control panels (control_panel()), then Tissue Scatter (col-7) + Leiden Panel (col-5)
# + Neighborhood Analysis Data / View panels and card rows below (rendered server-side in
#   output$neighborhood_cards): enrichment/scatter/distance row, composition card,
#   distance-to-islet profile card. Card prose lives in the tab_help() popover.

spatial_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tab_help("Spatial help", "show_phenotype_rules_spatial",
      object_classes_help(),
      tags$p(tags$strong("Tissue scatter."),
             " Every analysed cell of the chosen donor at its position. Drag to zoom; double-click or Reset Zoom to",
             " return. Show Cells limits the plot to cells in endocrine-object cores, or cores and their 20 \u00b5m",
             " peri-islet bands."),
      cell_typing_help(),
      tags$p(tags$strong("Islet-level UMAPs."),
             " Built from each object's mean marker values. ", textOutput(ns("embedding_method"), inline = TRUE),
             " Cluster Composition",
             " shows each Leiden cluster's mean core phenotype proportions, for the checked donor statuses and",
             " phenotypes."),
      phenotype_marker_definitions_ui(),
      # Explanatory prose for the neighborhood cards lives here, not inline
      # under the plots (Oct 2026 streamlining).
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Controls"),
      tags$p("Data controls choose what is plotted; View controls change how it is drawn. Click a panel's title",
             " to open or close it. The Neighborhood Analysis section has its own pair."),
      tags$p(tags$strong("Point size / transparency"),
             " (Neighborhood View controls) apply to Peri vs core proportion and Distance to nearest cell."),
      tags$p(tags$strong("Exclude outliers / z-threshold"),
             " apply to Peri enrichment by status and Distance to nearest cell. The checkbox and slider are synced with the Plot tab."),
      tags$p(tags$strong("Objects."),
             " Every neighborhood card except the Distance-to-islet profile uses the endocrine objects selected",
             " above the tabs; the profile reads per-cell data."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Peri enrichment by status & Peri vs core proportion"),
      tags$p("Per-phenotype peri-zone enrichment across disease stages, with peri vs core proportion scatter."),
      tags$p("Immune infiltration is the share of immune cells (every T cell, B cell, macrophage and other leukocyte, subtyped or not) among all analysed cells in the 20\u00b5m peri-islet band and in the core."),
      tags$p("Outlier exclusion in Peri enrichment by status follows the Neighborhood Analysis toggle (Tukey 1.5\u00d7IQR per donor status \u2014 enrichment z is already standardised, so the IQR fence is used instead of the modified z-score)."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Phenotype composition & enrichment"),
      tags$p("Which phenotypes are present in islet core vs peri-islet zones, and which are enriched vs tissue-wide background?"),
      tags$p("Enrichment z-scores (peri-islet mode) compare each phenotype's count in the peri-islet band with the count expected from its share of all analysed cells in the cohort: z = (observed \u2212 expected) / \u221aexpected (Poisson approximation). z > 0 = enriched; z < 0 = depleted. Core and peri-islet proportion modes show raw cell type fractions (0\u20131) within each compartment."),
      tags$p("The phenotype checkboxes and the Region / Summary choices apply to both the bar chart and the heatmap. Clip extreme z > 5 only saturates the display scale; it is not an outlier filter."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Distance to nearest cell & Distance-to-islet profile"),
      tags$p(tags$strong("Distance to nearest cell:"),
             " per-object distance (\u00b5m) from the centre of the object's core cells to the nearest peri-islet cell of each phenotype. NAs indicate no cells of that type in the peri-islet zone."),
      tags$p("Outlier exclusion follows the Neighborhood Analysis toggle (robust |z| from median/MAD on log(1+distance), per donor status \u2014 two-sided; threshold from the slider). Excluded points still appear as faint \u00d7 markers."),
      tags$p(tags$strong("Distance-to-islet profile:"),
             " per-cell signed distance-to-islet on a symlog axis \u2014 a density split by disease status (per phenotype or all cells) or a phenotype composition (stacked bars, one facet per status, zero-anchored at the islet border), spanning islet interior to bulk tissue. Negative = inside islet, positive = outside; the dashed line marks the border. Pooled per cell (descriptive) \u2014 for donor-fair group comparison see Peri enrichment by status.")),
    # ==== Data / View control panels (control_panel(), the same on every tab) ====
    fluidRow(column(12,
      control_panel("Data controls", open = TRUE,
        fluidRow(
          column(3, uiOutput(ns("donor_selector"))),
          column(3,
            radioButtons(ns("region_filter"), "Show cells",
                         choices = c("All" = "all", "Core + Peri" = "core_peri", "Core Only" = "core"),
                         selected = "all", inline = TRUE),
            checkboxGroupInput(ns("groups"), "Donor status",
                               choices = c("ND", "Aab+", "T1D"),
                               selected = c("ND", "Aab+", "T1D"), inline = TRUE)
          ),
          column(4, uiOutput(ns("phenotype_filter"))),
          column(2, uiOutput(ns("leiden_res_selector")))
        )
      ),
      control_panel("View controls",
        fluidRow(
          column(3,
            selectInput(ns("color_by"), "Color tissue cells by",
                        choices = c("Phenotype" = "phenotype", "Leiden Cluster" = "leiden"),
                        selected = "phenotype"),
            checkboxInput(ns("color_background"), "Color background cells", value = FALSE)
          ),
          column(3,
            selectInput("spatial_palette", "Phenotype colors",
                        choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                        selected = "High Contrast")
          ),
          column(3,
            # Non-namespaced: global donor palette, synced with Plot and Trajectory (app.R).
            selectInput("spatial_donor_palette", "Donor status colors",
                        choices = names(DONOR_PALETTES), selected = "Bright")
          )
        )
      )
    )),
    # ==== ROW 1: Tissue Scatter + Leiden ====
    fluidRow(style = "display: flex; flex-wrap: wrap;",
      # ---- Left: Tissue Scatter ----
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 8px;",
            h5("Tissue Scatter Plot", style = "margin: 0;"),
            actionButton(ns("scatter_reset_zoom"), "Reset Zoom",
                         style = "font-size: 12px; padding: 2px 10px;",
                         class = "btn-sm btn-outline-secondary")
          ),
          plotOutput(ns("tissue_scatter"), height = "700px",
                     dblclick = ns("scatter_dblclick"),
                     brush = brushOpts(id = ns("scatter_brush"), resetOnNew = TRUE)),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("download_spatial"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_tissue_scatter"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),

      # ---- Right: Leiden Panel ----
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
          h5("Islet-Level UMAPs"),
          uiOutput(ns("leiden_not_available")),
          fluidRow(
            column(6,
              h5("Leiden Clustering", style = "font-size: 13px; margin-top: 0;"),
              plotlyOutput(ns("leiden_umap"), height = "400px")
            ),
            column(6,
              h5("Donor Status", style = "font-size: 13px; margin-top: 0;"),
              plotOutput(ns("spatial_umap_donor"), height = 300),
              div(style = "text-align: right; margin-top: 6px;",
                downloadButton(ns("dl_spatial_umap_donor"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          ),
          hr(style = "margin: 8px 0;"),
          h5("Cluster Composition", style = "font-size: 15px;"),
          plotlyOutput(ns("cluster_composition"), height = "450px")
        )
      )
    ),

    # ==== ROW 2+: Neighborhood Analysis bar + cards ====
    uiOutput(ns("neighborhood_cards"))
  )
}
