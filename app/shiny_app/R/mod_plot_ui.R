# ---------- Plot module UI ----------
# Exports: plot_main_ui(id)
# All Plot-related inputs are namespaced by this module, except the app-wide synced controls
# (sidebar_remove_outliers, sidebar_outlier_z, sidebar_donor_palette), whose ids keep their old names.
# Controls sit in the tab's Data / View panels (control_panel(), 2026-10-17), as on every tab; the left
# sidebar they used to occupy is gone. The hidden Statistics tab reads the same inputs.

#' Plot tab: help, Data / View control panels, the two plot cards and the drill-down panel.
plot_main_ui <- function(id, extra_panel = NULL) {
  ns <- NS(id)
  tagList(
    singleton(tags$script(src = islet_asset("trajectory-links.js"))),
    column(12,
      tab_help("Plot help", "show_phenotype_rules_plot",
        object_classes_help(),
        tags$p(tags$strong("Data layer."),
               " Cell Populations shows the cell make-up of each endocrine object. Microenvironment shows the share",
               " of a region's area covered by Nerve, Capillary or Lymphatic structure objects."),
        tags$p(tags$strong("Regions."),
               " Core is the cells whose centre lies inside the object's outline. Peri-Islet is the cells in a 20 \u00b5m",
               " band around it, inside the tissue and outside every other endocrine object. Core+Peri pools both."),
        tags$p(tags$strong("Metrics."),
               " % is the phenotype's share of all analysed cells in the region. Counts are cells of the phenotype per",
               " object. Density divides the count by the region's area: the object's outline for Core, the 20 \u00b5m",
               " band for Peri-Islet, their sum for Core+Peri."),
        cell_typing_help(),
        tags$p(tags$strong("Islet diameter."),
               " The x axis is the diameter of a circle with the object's area. Diameter bins group objects for the",
               " summary statistic."),
        tags$p(tags$strong("Filters."),
               " Unchecking an autoantibody removes the Aab+ donors positive for it (ND and T1D donors are not",
               " affected). Exclude outliers drops objects whose robust z-score (median and MAD within",
               " each donor status) passes the threshold; the setting is shared with the Spatial tab."),
        if (!is.null(load_object_counts())) tags$p(tags$strong("Endocrine objects per donor."),
               " The card at the bottom counts Islets, Clusters and Singles in each donor's section, as a count or",
               " per mm\u00b2 of tissue, with the donor median per status; Singles by cell type splits each donor's",
               " Singles by the type of their one cell. It follows the Donor status filter, not the object selector."),
        tags$p(tags$strong("Linked points."),
               " With \"Show individual points\" on (the default), hover over an object in either plot to ring it in",
               " both. Click it in either plot to select it and open its segmentation map and single cells below;",
               " Clear selection releases it. An object selected on the Trajectory tab is ringed here too. With",
               " \"Color points by\" set to Donor ID, click a donor in the legend to show or hide it."),
        tags$p(tags$strong("Controls."),
               " Data controls choose what is plotted; View controls change how it is drawn. The point controls",
               " under \"Both plots\" apply to the size scatter and the distribution comparison alike. Click a panel's title to",
               " open or close it. Use browser zoom (Ctrl +/\u2212) to resize cards.")),
      control_panel("Data controls", open = TRUE,
        fluidRow(
          column(3,
            # Non-namespaced: synced with the Spatial tab's copy (app.R).
            object_class_selector("object_classes"),
            selectInput(ns("mode"), "Data layer",
                        choices = c("Microenvironment", "Cell Populations"),
                        selected = "Cell Populations"),
            uiOutput(ns("dynamic_selector")),
            uiOutput(ns("metric_selector"))
          ),
          column(3,
            checkboxGroupInput(ns("groups"), "Donor status",
                               choices = c("ND", "Aab+", "T1D"),
                               selected = c("ND", "Aab+", "T1D"),
                               inline = TRUE),
            tags$label("Autoantibodies (Aab+ donors)", class = "control-label"),
            uiOutput(ns("aab_warning")),
            checkboxGroupInput(ns("aab_flags"), NULL,
                               choices = c("GADA" = "AAb_GADA",
                                           "IA2A" = "AAb_IA2A",
                                           "ZnT8A" = "AAb_ZnT8A",
                                           "IAA" = "AAb_IAA",
                                           "mIAA" = "AAb_mIAA"),
                               selected = c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA"),
                               inline = TRUE)
          ),
          column(3,
            uiOutput(ns("age_filter_ui")),
            uiOutput(ns("gender_filter_ui"))
          ),
          column(3,
            sliderInput(ns("diam_range"), "Islet diameter range (\u00b5m)",
                        min = 0, max = 500, value = c(0, 350), step = 10),
            sliderInput(ns("binwidth"), "Diameter bin width (\u00b5m)",
                        min = 1, max = 75, value = 50, step = 1),
            checkboxInput(ns("exclude_zero_top"), "Exclude zero values (size distribution)", value = FALSE),
            checkboxInput(ns("exclude_zero_dist"), "Exclude zero values (distribution comparison)",
                          value = FALSE)
          )
        ),
        fluidRow(
          column(3,
            # Non-namespaced: synced with the Trajectory and Spatial outlier controls (app.R).
            tags$div(
              title = GLOBAL_OUTLIER_TOOLTIP,
              checkboxInput("sidebar_remove_outliers", "Exclude outliers", value = FALSE)
            )
          ),
          column(3,
            sliderInput("sidebar_outlier_z", "Outlier z-threshold",
                        min = 0.5, max = 10, value = 3, step = 0.5)
          ),
          column(6,
            checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE),
            uiOutput(ns("plot_outlier_info"))
          )
        )
      ),
      control_panel("View controls",
        fluidRow(
          column(3,
            radioButtons(ns("stat"), "Summary statistic (Islet Size Distribution)",
                         choices = c("Mean\u00b1SE" = "mean_se",
                                     "Mean\u00b1SD" = "mean_sd",
                                     "Median + IQR" = "median_iqr"),
                         selected = "mean_se"),
            # Non-namespaced: global donor palette, synced with Trajectory and Spatial (app.R).
            selectInput("sidebar_donor_palette", "Donor status colors",
                        choices = names(DONOR_PALETTES), selected = "Bright")
          ),
          # Point controls apply to both plots (owner, 2026-10-19: the View controls must apply to both,
          # not only the Distribution Comparison).
          column(5,
            h5("Both plots"),
            fluidRow(
              column(6,
                selectInput(ns("plot_color_by"), "Color points by:",
                            choices = c("Donor Status" = "donor_status", "Donor ID" = "donor_id"),
                            selected = "donor_status"),
                sliderInput(ns("pt_size"), "Point size", min = 0.3, max = 4.0, value = 3.0, step = 0.1),
                sliderInput(ns("pt_alpha"), "Point transparency", min = 0, max = 0.95, value = 0.4, step = 0.05)
              ),
              column(6,
                checkboxInput(ns("show_points"), "Show individual points", value = TRUE),
                checkboxInput(ns("log_scale"), "Log scale y-axis", value = FALSE)
              )
            )
          ),
          column(4,
            h5("Islet Size Distribution"),
            radioButtons(ns("add_smooth"), "Trend line", choices = c("None", "LOESS"),
                         selected = "None", inline = TRUE),
            checkboxInput(ns("log_scale_x"), "Log2 scale x-axis", value = FALSE),
            h5("Distribution Comparison"),
            radioButtons(ns("dist_type"), "Plot type", choices = c("Violin", "Box"),
                         selected = "Violin", inline = TRUE)
          )
        )
      ),
      div(style = "display: flex; align-items: center; gap: 12px; margin-bottom: 12px;",
        actionButton(ns("clear_link_selection"), "Clear selection", class = "btn-sm btn-outline-secondary"))
    ),
    # Left card: islet size distribution
    column(6,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid var(--islet-border-strong); border-radius: 8px; height: clamp(600px, 75vh, 1400px); display: flex; flex-direction: column; gap: 12px;",
        div(style = "display: flex; align-items: center; justify-content: space-between; margin-top: 0;",
          h5("Islet Size Distribution", style = "margin: 0; color: var(--islet-fg);"),
          downloadButton(ns("dl_summary"), "Summary CSV",
                         class = "btn-sm btn-outline-secondary",
                         style = "font-size: 11px; padding: 2px 8px;")
        ),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("plt"), height = "calc(100% - 28px)")
        )
      )
    ),
    # Right card: distribution comparison. Widen to fill the row when no extra panel (e.g. the AI chat
    # card) is present, otherwise its two columns are left as dead space.
    column(if (is.null(extra_panel)) 6 else 4,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid var(--islet-border-strong); border-radius: 8px; height: clamp(600px, 75vh, 1400px); display: flex; flex-direction: column; gap: 12px;",
        div(style = "display: flex; align-items: center; justify-content: space-between; margin-top: 0;",
          h5("Distribution Comparison", style = "margin: 0; color: var(--islet-fg);"),
          downloadButton(ns("dl_distribution_csv"), "CSV",
                         class = "btn-sm btn-outline-secondary",
                         style = "font-size: 11px; padding: 2px 8px;")
        ),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("dist"), height = "calc(100% - 28px)")
        )
      )
    ),
    # Extra panel (e.g., AI chat) injected between cards and full-width panels
    extra_panel,
    # Full-width card below both plot cards: segmentation/drill-down viewer
    column(12,
      uiOutput(ns("segmentation_viewer_panel"))
    ),
    # Per-donor counts of each endocrine object class (uniform root only): the one place Singles are shown.
    if (!is.null(load_object_counts())) column(12, object_counts_card_ui(ns))
  )
}
