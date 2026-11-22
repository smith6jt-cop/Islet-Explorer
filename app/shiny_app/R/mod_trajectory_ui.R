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
               sprintf(" Islets and clusters with at least %d total core cells are eligible.", TRAJ_MIN_CELLS),
               " All measured core cell types count toward this cutoff. An object needs a nonzero trajectory input",
               " to receive pseudotime. The class selector on the Plot and Spatial tabs applies only there."),
        tags$p(tags$strong("Y axis."),
               " Marker expression is the mean of the object's core cells on the cohort marker scale (0-1). Cell-type",
               " ratio is log(1 + cells / core endocrine cells) for the chosen phenotype and region.",
               " Select Core only to compare endocrine subtypes. All core endocrine cells form the denominator."),
        cell_typing_help(),
        unclassified_help_ui()),
      # The status line reports how many eligible objects received pseudotime.
      uiOutput(ns("traj_status")),
      control_grid(columns = 4,
        control_group("Measure",
          # Y-axis measure. Marker expression is the core-cell protein mean from .X
          # on the cohort marker scale; the cell-type ratio is Barlow et al.'s Fig 2F y-axis,
          # log(# cells / # endocrine cells), which removes the islet-size
          # confound that a raw count or proportion carries.
          selectInput(ns("traj_y_mode"), "Y axis",
                      choices = c("Marker expression" = "marker",
                                  "Cell-type ratio (log1p per core endocrine cell)" = "ratio"),
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
            checkboxInput(ns("traj_spread_ties"), "Spread tied points (display only)", value = TRUE),
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
            tags$p("Each fitted object's selected value against pseudotime. The three strips below use the same 25 bins.",
                   " Each row shows how that disease status's islets are distributed across pseudotime; hover for the",
                   " islet count and the status contribution within each bin."),
            tags$p(tags$strong("Pseudotime."),
                   " Diffusion pseudotime (scanpy DPT) on each object's cell make-up, after Barlow et al. 2024:",
                   " log(1 + cells / endocrine cells) for each cell class. The endocrine classes are counted in the",
                   " core only and exocrine cells in the 20 \u00b5m peri-islet band only; macrophages enter as Gal3\u207a",
                   " and Gal3\u207b, and the T cell, B cell, stromal, endothelial and neural classes pooled across core and peri before normalization.",
                   " Each retained phenotype additionally contributes its summed nerve, lymphatic and blood-vessel overlap fractions, normalized in the same way. Conflict, Doublet suspected, Tissue, Unclassified and unresolved-lineage cells are not inputs.",
                   " Inputs that do not vary are dropped. Objects are linked to their 15 nearest neighbors by cosine",
                   " distance. The fit finds the two most distant profiles in diffusion space and selects an ND",
                   " islet near the end supported by more ND donors. Raw diffusion pseudotime starts at zero",
                   " at this islet. Objects with the same make-up are fitted once and share a pseudotime."),
            tags$p(tags$strong("Reading it."),
                   " Tied points are spread for visibility, including equal minima and maxima.",
                   " Hover shows the measured values. Clear Spread tied points to see exact positions.",
                   " The orientation is a modeling choice, and an ordering of objects along it does not show",
                   " disease progression; donor differences can shape it."),
            linked_help,
            unclassified_help_ui())),
          trajectory_linked_plotly(ns("traj_scatter"), height = 850),
          plotlyOutput(ns("traj_heatmap"), height = 112)
        )
      ),
      # Keep the UMAP and pseudotime cards in one row; segmentation gets its own
      # full-width row so it cannot leave an empty column below pseudotime.
      column(5,
        div(class = "card trajectory-main-card", style = "padding: 15px;",
          card_title_help("UMAPs", help = list(
            tags$p(tags$strong("UMAPs and Leiden clusters."),
                   " These panels use the saved layout shown on the Spatial tab. The composition embedding uses",
                   " each object's core endocrine proportions (Beta, Alpha, Delta and CHGA⁺ hormone-negative).",
                   " Objects with a pseudotime and no saved coordinates still appear in the scatter and PAGA.",
                   " The upper UMAP shows donor status and the lower UMAP shows the selected feature."),
            linked_help,
            unclassified_help_ui())),
          tags$p(textOutput(ns("traj_umap_coverage"), inline = TRUE)),
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
          h5("Segmentation Map", style = "margin-top: 0;"),
          uiOutput(ns("segmentation_viewer_panel"))
        )
      )
    ),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          card_title_help("Multi-Feature Heatmap", style = "color: var(--islet-heading); margin: 0 0 8px;",
            help = list(tags$p("Pseudotime is cut into equal bins; each row is a feature's bin means, standardized",
                               " across bins and clipped at \u00b12.5."),
                        unclassified_help_ui())),
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
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          card_title_help("Pseudotime by donor annotation",
            help = list(tags$p("Each donor contributes one point at their median fitted-islet pseudotime.",
                               " Horizontal bars span the middle 50% of that donor's islets.",
                               " Numeric measurements with missing values are omitted.",
                               " An open downward triangle marks a C-peptide upper limit."))),
          uiOutput(ns("donor_annotation_field")),
          plotlyOutput(ns("donor_annotation_plot"), height = 390)
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
        div(class = "card paga-comparison-card", style = "padding: 15px; margin-bottom: 20px; overflow: visible;",
          card_title_help("PAGA Force-Directed Layout",
            style = "color: var(--islet-heading); margin: 0 0 8px;",
            help = list(tags$p("Barlow et al.'s PAGA-initialized force-directed layout of the same graph (cosine,",
                               " 15 neighbors, Leiden resolution 1). Each point is an object; the axes are layout",
                               " coordinates with no units, so only relative position means anything. It is built on",
                               " request, on the same objects as the pseudotime. Both panels show that same fit, each with its own color choice."),
                        tags$p("The DPT root arrow marks the ND islet where baseline diffusion pseudotime starts.",
                               " The fit selects the root in diffusion space before drawing this two-dimensional layout."),
                        tags$p("Choose a donor characteristic under Color by to see its recorded value across the layout. Each islet inherits its donor’s value. For numeric gradients, gray points have a missing, inapplicable, or censored value; hover for details."),
                        tags$p("At zero, both panels use the original PAGA positions and diffusion pseudotime.",
                               " Age guidance shifts each donor’s islets together along the direction of increasing age",
                               " estimated from donor centers. At full strength the donor medians follow age in years",
                               " while islets from the same donor retain their relative positions.",
                               " The guided score shifts each donor’s pseudotimes together and may extend outside 0–1."),
                        tags$p("Hover over an object to highlight it in both PAGA panels and the other linked plots. Click to open its segmentation map."),
                        unclassified_help_ui())),
          checkboxInput(ns("paga_show_layout"), "Build and show PAGA layout",
                        value = FALSE),
          conditionalPanel(
            condition = "input.paga_show_layout", ns = ns,
          fluidRow(
            column(3,
              sliderInput(ns("paga_point_size"), "Point size",
                          min = 0.3, max = 4.0, value = 1.6, step = 0.1)
            ),
            column(4,
              sliderInput(ns("paga_age_weight"), "Donor age guidance",
                          min = 0, max = 1, value = 0, step = 0.05)
            ),
            column(5, uiOutput(ns("paga_note")))
          ),
          fluidRow(
            column(6,
              selectInput(ns("paga_color_by"), "Color by", choices = PAGA_COLOR_CHOICES,
                          selected = "pseudotime", width = "100%"),
              plotlyOutput(ns("paga_layout"), height = 520)
            ),
            column(6,
              selectInput(ns("paga_color_right"), "Color by", choices = PAGA_COLOR_CHOICES,
                          selected = "donor_status", width = "100%"),
              plotlyOutput(ns("paga_layout_right"), height = 520)
            )
          ),
          div(style = "text-align: right; margin-top: 8px;",
            downloadButton(ns("dl_paga_csv"), "Download CSV",
                           class = "btn-sm btn-outline-secondary"))
          )
        )
      )
    )
  )
}
