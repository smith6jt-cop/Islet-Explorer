# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Data / View control panels (control_panel()), then Tissue Scatter (col-7) + Leiden clusters (col-5: composition / membership),
# then the two islet-level UMAPs side by side in a full-width row
# + Neighborhood Analysis cards below (rendered server-side in output$neighborhood_cards): insulitis prevalence /
#   immune cells and beta-cell content / immune gradient row, composition card, distance-to-islet profile card
#   (owner review, 2026-10-20). Each panel and card's prose lives in its own "?" popover (card_help()).

spatial_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # ==== Data / View control panels (control_panel(), the same on every tab) ====
    fluidRow(column(12,
      control_panel("Data controls", open = TRUE,
        help = list(
          phenotype_rules_button("show_phenotype_rules_spatial"),
          tags$p("Data controls choose what the cards above Neighborhood Analysis plot. Click the panel's title to",
                 " open or close it. The Neighborhood Analysis section has its own pair."),
          object_classes_help(),
          tags$p(tags$strong("Show cells"),
                 " limits the tissue scatter to cells in endocrine-object cores, or cores and their 20 \u00b5m",
                 " peri-islet bands."),
          cell_typing_help(),
          phenotype_marker_definitions_ui()),
        fluidRow(
          column(3,
            # Non-namespaced: synced with the Plot tab's copy (app.R).
            object_class_selector("spatial_object_classes"),
            uiOutput(ns("donor_selector"))),
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
    # ==== ROW 1: Tissue Scatter + Leiden cluster composition ====
    fluidRow(style = "display: flex; flex-wrap: wrap;",
      # ---- Left: Tissue Scatter ----
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 8px;",
            card_title_help("Tissue Scatter Plot", help = list(
              tags$p("Every analysed cell of the chosen donor at its position. Drag to zoom; double-click or Reset",
                     " Zoom to return. Show cells (Data controls) limits the plot to cells in endocrine-object cores,",
                     " or cores and their 20 \u00b5m peri-islet bands."))),
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

      # ---- Right: Leiden cluster composition ----
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
          card_title_help("Leiden clusters", style = "margin: 0 0 8px;", help = list(
            tags$p(tags$strong("Composition"),
                   " shows each Leiden cluster's mean core phenotype proportions, for the checked donor statuses and",
                   " phenotypes: what each cluster is. Its Beta, Alpha, Delta and CHGA+ hormone-negative bars restate",
                   " the proportions the clusters are built from; the immune and residual bars do not."),
            tags$p(tags$strong("Donor status"), " and ", tags$strong("Donor"),
                   " show the proportion of each cluster's objects from each donor status or each donor: which",
                   " composition states are disease-associated, and whether a cluster belongs to a few donors. The",
                   " number above each bar is the cluster's donor count, and hover names its largest donor."),
            tags$p("The clusters are those of the Islet-Level UMAPs below, at the resolution chosen in the Data",
                   " controls."))),
          uiOutput(ns("leiden_not_available")),
          radioButtons(ns("cluster_view"), "Show",
                       c("Composition" = "composition", "Donor status" = "status", "Donor" = "donor"),
                       selected = "composition", inline = TRUE),
          plotlyOutput(ns("cluster_composition"), height = "650px")
        )
      )
    ),

    # ==== ROW 1b: the two islet-level UMAPs, one size (owner, 2026-10-19) ====
    # Side by side across the full width: in the old half-width card each had ~280 px, so a legend that did not
    # fit made plotly shrink that plot alone. Both use spatial_umap_plotly(); each downloads from its camera button.
    fluidRow(column(12,
      div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
        card_title_help("Islet-Level UMAPs", style = "margin: 0 0 8px;", help = list(
          tags$p(textOutput(ns("embedding_method"), inline = TRUE)),
          tags$p("The left UMAP is coloured by Leiden cluster. The right one, the same layout, is coloured by each",
                 " object's immune proportion in its 20 \u00b5m band (immune cells among all its band cells), saturating",
                 " at the 95th percentile; hover gives the value. The Trajectory tab shows the layout by donor status."))),
        fluidRow(
          column(6,
            h5("Leiden Clustering", style = "font-size: 13px; margin-top: 0;"),
            plotlyOutput(ns("leiden_umap"), height = "480px")
          ),
          column(6,
            h5("Band immune proportion", style = "font-size: 13px; margin-top: 0;"),
            plotlyOutput(ns("spatial_umap_immune"), height = "480px")
          )
        )
      )
    )),

    # ==== ROW 2+: Neighborhood Analysis bar + cards ====
    uiOutput(ns("neighborhood_cards"))
  )
}
