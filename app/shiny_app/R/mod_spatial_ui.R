# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4)
# + Neighborhood Analysis control bar and card rows below (rendered server-side in
#   output$neighborhood_cards): enrichment/scatter/distance row, composition card,
#   distance-to-islet profile card. Card prose lives in the tab_help() popover.

spatial_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tab_help("Spatial help", "show_phenotype_rules_spatial",
      tags$p("Drag to zoom, double-click to reset the tissue scatter plot."),
      tags$p("Phenotype pickers list the 12 retained labels. The residual and parent labels ",
             "(Endocrine, Immune, T cell, Endothelial, APCs, Acinar, ECAD+, Structural, Unknown) are excluded."),
      tags$p("Phenotype provenance: the original 15 donors retain saved labels; the six additions use reviewed reconstructed rules. Spatial phenotype comparisons span both methods."),
      phenotype_marker_definitions_ui(),
      # Explanatory prose for the neighborhood cards lives here, not inline
      # under the plots (Oct 2026 streamlining).
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Neighborhood Analysis controls"),
      tags$p(tags$strong("Point size / transparency"),
             " apply to Peri vs core proportion and Distance to nearest cell."),
      tags$p(tags$strong("Exclude outliers / z-threshold"),
             " apply to Peri enrichment by status and Distance to nearest cell. The checkbox and slider are synced with the Plot tab."),
      tags$p(tags$strong("Min cells/islet"),
             " filters islets for every neighborhood card except the Distance-to-islet profile, which reads per-cell data."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Peri enrichment by status & Peri vs core proportion"),
      tags$p("Per-phenotype peri-zone enrichment across disease stages, with peri vs core proportion scatter."),
      tags$p("Immune infiltration quantifies the proportion of immune cells among all cells in the peri-islet zone (20\u00b5m expansion) and islet core. Higher fractions indicate increased immune surveillance or active infiltration."),
      tags$p("Outlier exclusion in Peri enrichment by status follows the Neighborhood Analysis toggle (Tukey 1.5\u00d7IQR per donor status \u2014 enrichment z is already standardised, so the IQR fence is used instead of the modified z-score)."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Phenotype composition & enrichment"),
      tags$p("Which phenotypes are present in islet core vs peri-islet zones, and which are enriched vs tissue-wide background?"),
      tags$p("Enrichment z-scores (peri-islet mode) compare observed immune cell counts in the peri-islet zone vs expected counts based on tissue-wide proportions (Poisson model). z > 0 = enriched; z < 0 = depleted. Core and peri-islet proportion modes show raw cell type fractions (0\u20131) within each compartment."),
      tags$p("The phenotype checkboxes and the Region / Summary choices apply to both the bar chart and the heatmap. Clip extreme z > 5 only saturates the display scale; it is not an outlier filter."),
      tags$h6(style = "margin-top: 12px; font-weight: 700;", "Distance to nearest cell & Distance-to-islet profile"),
      tags$p(tags$strong("Distance to nearest cell:"),
             " per-islet minimum Euclidean distance (\u00b5m) from islet core centroid to nearest peri-islet cells of each phenotype. NAs indicate no cells of that type in the peri-islet zone."),
      tags$p("Outlier exclusion follows the Neighborhood Analysis toggle (robust |z| from median/MAD on log(1+distance), per donor status \u2014 two-sided; threshold from the slider). Excluded points still appear as faint \u00d7 markers."),
      tags$p(tags$strong("Distance-to-islet profile:"),
             " per-cell signed distance-to-islet on a symlog axis \u2014 a density split by disease status (per phenotype or all cells) or a phenotype composition (stacked bars, one facet per status, zero-anchored at the islet border), spanning islet interior to bulk tissue. Negative = inside islet, positive = outside; the dashed line marks the border. Pooled per cell (descriptive) \u2014 for donor-fair group comparison see Peri enrichment by status.")),
    # ==== ROW 1 (existing): Controls + Scatter + Leiden ====
    fluidRow(style = "display: flex; flex-wrap: wrap;",
      # ---- Left sidebar: Controls ----
      column(2,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
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
          uiOutput(ns("phenotype_filter")),
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
                      selected = "Bright"),
          hr(style = "margin: 8px 0;"),
          downloadButton(ns("download_spatial"), "Download CSV", style = "width: 100%; font-size: 13px;")
        )
      ),

      # ---- Center: Tissue Scatter ----
      column(6,
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
            downloadButton(ns("dl_tissue_scatter"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),

      # ---- Right: Leiden Panel ----
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
          h5("Islet-Level UMAPs"),
          textOutput(ns("embedding_method")),
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
