# plotlyOutput() returns a tag list, and tagAppendAttributes() on a tag list adds its arguments as children, so
# the class used to print as text under each plot and never reached the widget. Put it on the widget's div.
trajectory_linked_plotly <- function(output_id, height) {
  widget <- plotlyOutput(output_id, height = height)
  widget[[1]] <- tagAppendAttributes(widget[[1]], class = "trajectory-linked-plot")
  widget
}

trajectory_ui <- function(id) {
  ns <- NS(id)
  linked_help <- tags$p(tags$strong("Linked plots."),
    " Hover over an object in either UMAP or the pseudotime scatter to highlight it in all three. Click to select it",
    " and open its segmentation map; Clear selection (Data controls) releases it. An object selected on the Plot",
    " tab is ringed",
    " here too.")
  tagList(
    singleton(tags$script(src = islet_asset("trajectory-links.js"))),
    tags$style(HTML(".trajectory-linked-plot.html-widget.plotly { padding-top: 0; margin-top: 28px; }")),
    # Frozen under the page header while the tab scrolls (tab_controls()).
    tab_controls(
      control_panel("Data controls",
      help = list(
        phenotype_rules_button("show_phenotype_rules_traj"),
        tags$p("Data controls choose what is fitted and plotted. Click the panel's title to open or close it. The",
               " heatmap and PAGA cards keep their own selectors."),
        tags$p(tags$strong("Fitted objects."),
               " Every endocrine object with ten or more core cells, whatever its class, as in Barlow et al. 2024,",
               " who fitted islets of ten or more cells. The Endocrine objects selector of the other tabs does not",
               " apply here and is hidden. Objects with fewer cells carry too little make-up to order: fitted on",
               " every object, the pseudotime barely tracked the beta-cell ratio."),
        tags$p(tags$strong("Y axis."),
               " Marker expression is the mean of the object's core cells on the cohort marker scale (0-1). Cell-type",
               " ratio is log(1 + cells / endocrine cells) for the chosen phenotype in the chosen region."),
        cell_typing_help()),
      # The fitted cohort is every object with TRAJ_MIN_CELLS or more core cells (since 2026-10-18); the
      # status line reports how many objects were fitted.
      uiOutput(ns("traj_status")),
      control_grid(columns = 4,
        control_group("Measure",
          # Y-axis measure. Marker expression is the core-cell protein mean from .X
          # on the cohort marker scale; the cell-type ratio is Barlow et al.'s Fig 2F y-axis,
          # log(# cells / # endocrine cells), which removes the islet-size
          # confound that a raw count or proportion carries.
          selectInput(ns("traj_y_mode"), "Y axis",
                      choices = c("Marker expression" = "marker",
                                  "Cell-type ratio (log cells / endocrine cell)" = "ratio"),
                      selected = "ratio")
        ),
        control_group("Region",
          conditionalPanel(
            condition = "input.traj_y_mode == 'ratio'", ns = ns,
            selectInput(ns("traj_ratio_region"), "Count cells in",
                        choices = BARLOW_RATIO_REGIONS, selected = "core")
          )
        ),
        control_group("Feature",
          uiOutput(ns("traj_feature_selector"))
        ),
        # Owner, 2026-10-26: the Clear selection button lives in the Data controls.
        control_group("Selected object", clear_selection_button(ns))
      )
      ),
      control_panel("View controls",
        help = list(
          tags$p("View controls change how the pseudotime scatter and the UMAPs are drawn. Click the panel's title",
                 " to open or close it."),
          tags$p(tags$strong("Pseudotime scale."),
                 " Rank-normalized spreads objects evenly from 0 to 1; Raw DPT keeps the diffusion distances."),
          tags$p(tags$strong("Trend lines."),
                 " LOESS fits, overall or per donor status, with an optional 95% confidence band.")),
        control_grid(columns = 3,
          control_group("Pseudotime and trends",
            selectInput(ns("pt_scale"), "Pseudotime scale",
                       choices = c("Rank-normalized" = "rank",
                                   "Raw DPT" = "raw"),
                       selected = "rank"),
            selectInput(ns("traj_show_trend"), "Trend lines",
                       choices = c("None" = "none",
                                   "Overall" = "overall",
                                   "By Donor Status" = "by_donor"),
                       selected = "by_donor"),
            checkboxInput(ns("traj_trend_ci"), "95% CI band on trends",
                          value = FALSE)
          ),
          control_group("Points",
            selectInput(ns("traj_point_size"), "Point size by",
                       choices = c("Cell Count" = "total_cells",
                                   "Islet Diameter" = "islet_diam_um",
                                   "Uniform" = "uniform"),
                       selected = "total_cells"),
            sliderInput(ns("traj_point_size_slider"), "Point size",
                        min = 0.5, max = 5.0, value = 1.5, step = 0.1),
            sliderInput(ns("traj_alpha"), "Point transparency",
                        min = 0, max = 0.95, value = 0.4, step = 0.05)
          ),
          control_group("Color",
            selectInput(ns("traj_color_by"), "Color points by",
                       choices = c("Donor Status" = "donor_status",
                                   "Donor ID" = "donor_id",
                                   "Leiden (res 0.3)" = "leiden_0.3",
                                   "Leiden (res 0.5)" = "leiden_0.5",
                                   "Leiden (res 0.8)" = "leiden_0.8",
                                   "Leiden (res 1.0)" = "leiden_1.0"),
                       selected = "donor_status"),
            selectInput("traj_donor_palette", "Donor status colors",
                        choices = names(DONOR_PALETTES),
                        selected = "Bright"),
            div(style = "border: 1px solid var(--islet-border-strong); padding: 8px 10px; border-radius: 5px; background-color: var(--islet-surface-2);",
              uiOutput(ns("traj_legend"))
            )
          )
        )
      )
    ),
    fluidRow(class = "trajectory-overview-row",
      # Left card: Scatterplot and heatmaps (narrowed from col-8)
      column(7,
        div(class = "card trajectory-controls trajectory-main-card", style = "padding: 15px;",
          card_title_help("Pseudotime", help = list(
            tags$p("Each fitted object's Y-axis value (Data controls) against its pseudotime. The strip below",
                   " cuts pseudotime into 25 bins: a bin's color is the donor status holding the largest proportion of",
                   " its own objects there, and its depth is the feature's mean in the bin."),
            tags$p(tags$strong("Pseudotime."),
                   " Diffusion pseudotime (scanpy DPT) on each object's cell make-up, after Barlow et al. 2024: for",
                   " each listed phenotype, log(1 + cells / endocrine cells) in the core and in the 20 \u00b5m",
                   " peri-islet band. Features that do not vary are dropped (peri-islet endocrine values are always",
                   " zero). Objects are linked to their 15 nearest neighbors by cosine distance. The root is the ND",
                   " object nearest the center of the ND objects in diffusion components 1-3, and the direction is set",
                   " so the beta-cell ratio falls along pseudotime. Objects with fewer than ten core cells get no",
                   " pseudotime. Objects with the same make-up are fitted once and share a pseudotime."),
            tags$p(tags$strong("Reading it."),
                   " The orientation is a modeling choice, and an ordering of objects along it does not show",
                   " disease progression; donor differences can shape it."),
            linked_help)),
          trajectory_linked_plotly(ns("traj_scatter"), height = 850),
          plotlyOutput(ns("traj_heatmap"), height = 100)
        )
      ),
      # Keep the UMAP and pseudotime cards in one row; segmentation gets its own
      # full-width row so it cannot leave an empty column below pseudotime.
      column(5,
        div(class = "card trajectory-main-card", style = "padding: 15px;",
          card_title_help("UMAPs", help = list(
            tags$p(tags$strong("UMAPs and Leiden clusters."),
                   " Built from each fitted object's endocrine composition (the Endocrine Beta, Alpha, Delta and CHGA⁺",
                   " hormone-negative proportions of its core endocrine cells), with one neighbor graph for UMAP and",
                   " Leiden and no donor correction; the same layout as the Spatial tab's UMAPs. The upper UMAP is",
                   " colored by donor status, the lower by the feature chosen in the Data controls."),
            linked_help)),
          h5("UMAP: Donor Status", style = "font-size: 13px; margin-top: 8px;"),
          trajectory_linked_plotly(ns("traj_umap_donor"), height = 320),
          div(style = "text-align: right; margin-top: 6px;",
            downloadButton(ns("dl_umap_donor"), "Download PNG",
                           class = "btn-sm btn-outline-primary")),
          hr(style = "margin: 12px 0;"),
          h5(textOutput(ns("traj_umap_feature_title"), inline = TRUE),
             style = "font-size: 13px; margin-top: 0;"),
          trajectory_linked_plotly(ns("traj_umap_feature"), height = 320),
          div(style = "text-align: right; margin-top: 6px;",
            downloadButton(ns("dl_umap_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_umap_feature"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      )
    ),
    fluidRow(
      column(12,
        div(class = "card trajectory-segmentation-card", style = "padding: 15px; margin-bottom: 20px;",
          h5("Segmentation map", style = "margin-top: 0;"),
          uiOutput(ns("segmentation_viewer_panel"))
        )
      )
    ),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          card_title_help("Multi-Feature Heatmap", style = "color: var(--islet-heading); margin: 0 0 8px;",
            help = list(tags$p("Pseudotime is cut into equal bins; each row is a feature's bin means, standardized",
                               " across bins and clipped at \u00b12.5."))),
          fluidRow(
            column(8, uiOutput(ns("multi_heatmap_marker_selector"))),
            column(4, sliderInput(ns("multi_heatmap_nbins"), "Number of bins:", min = 10, max = 40, value = 20, step = 5))
          ),
          plotOutput(ns("traj_multi_heatmap"), height = "auto"),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("dl_multi_heatmap_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"),
            downloadButton(ns("dl_multi_heatmap"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      )
    ),
    # ---- PAGA force-directed layout (Barlow et al. Fig 2B/2E) ----------
    # The feature-vs-pseudotime scatter above cannot be de-overplotted by any
    # axis transform: the pile-up is at EXACT zero (95.6% of islets have no CD8a
    # T cell, 98.2% no B cell) and log(0) is still 0. This 2-D graph embedding
    # is the plot that actually separates islets -- 8.1% overplotted vs 84-96%.
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          card_title_help("PAGA force-directed layout",
            style = "color: var(--islet-heading); margin: 0 0 8px;",
            help = list(tags$p("Barlow et al.'s PAGA-initialized force-directed layout of the same graph (cosine,",
                               " 15 neighbors, Leiden resolution 1). Each point is an object; the axes are layout",
                               " coordinates with no units, so only relative position means anything. It is built on",
                               " request, on the same objects as the pseudotime."))),
          checkboxInput(ns("paga_show_layout"), "Build and show PAGA layout",
                        value = FALSE),
          conditionalPanel(
            condition = "input.paga_show_layout", ns = ns,
          fluidRow(
            column(3,
              selectInput(ns("paga_color_by"), "Color islets by:",
                          choices = c("Pseudotime" = "pseudotime",
                                      "Donor status" = "donor_status",
                                      "Leiden cluster" = "leiden",
                                      "Donor" = "case_id",
                                      "Selected feature" = "feature"),
                          selected = "pseudotime")
            ),
            column(3,
              sliderInput(ns("paga_point_size"), "Point size:",
                          min = 0.3, max = 4.0, value = 1.6, step = 0.1)
            ),
            column(6, uiOutput(ns("paga_note")))
          ),
          plotlyOutput(ns("paga_layout"), height = 620),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("dl_paga_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"))
          )
        )
      )
    )
  )
}
