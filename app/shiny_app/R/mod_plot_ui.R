# ---------- Plot module UI ----------
# Exports: plot_main_ui(id)
# All Plot-related inputs are namespaced by this module, except the app-wide synced controls
# (sidebar_remove_outliers, sidebar_outlier_z, sidebar_donor_palette), whose ids keep their old names.
# Controls sit in the tab's Data / View panels (control_panel(), 2026-10-17), as on every tab; the left
# sidebar they used to occupy is gone. The hidden Statistics tab reads the same inputs.

#' Plot tab: Data / View control panels, the two plot cards and the drill-down panel. Each panel and card carries
#' its own help ("?" beside its title).
plot_main_ui <- function(id, extra_panel = NULL) {
  ns <- NS(id)
  tagList(
    singleton(tags$script(src = islet_asset("trajectory-links.js"))),
    # Frozen under the page header while the tab scrolls (tab_controls()); its parent row holds the whole tab.
    column(12, class = "tab-controls",
      control_panel("Data controls", open = TRUE,
        help = list(
          phenotype_rules_button("show_phenotype_rules_plot"),
          tags$p("Data controls choose what is plotted, in both plot cards. Click the panel's title to open or close it."),
          object_classes_help(),
          tags$p(tags$strong("Data layer."),
                 " Cell Populations shows the cell make-up of each endocrine object. Microenvironment shows the proportion",
                 " of a region's area covered by Nerve, Capillary or Lymphatic structure objects."),
          tags$p(tags$strong("Regions."),
                 " Core is the cells whose centre lies inside the object's outline. Peri-Islet is the cells in a 20",
                 " \u00b5m band around it, inside the tissue and outside every other endocrine object. Core+Peri pools",
                 " both."),
          tags$p(tags$strong("Metrics."),
                 " % is the phenotype's proportion of all analysed cells in the region. Counts are cells of the phenotype per",
                 " object. Density divides the count by the region's area: the object's outline for Core, the 20 \u00b5m",
                 " band for Peri-Islet, their sum for Core+Peri."),
          cell_typing_help(),
          tags$p(tags$strong("Filters."),
                 " Unchecking an autoantibody removes the Aab+ donors positive for it (ND and T1D donors are not",
                 " affected). The diameter range keeps objects within it; the bin width sets the diameter bins of the",
                 " Against diameter summary. Exclude zero values applies to the card it names. Exclude outliers drops",
                 " objects whose robust z-score (median and MAD within each donor status) passes the threshold.")),
        # Four titled groups in equal columns (owner, 2026-10-21: a more organized control center).
        control_grid(
          control_group("Objects and measure",
            # Non-namespaced: synced with the Spatial tab's copy (app.R).
            object_class_selector("object_classes"),
            selectInput(ns("mode"), "Data layer",
                        choices = c("Microenvironment", "Cell Populations"),
                        selected = "Cell Populations"),
            uiOutput(ns("dynamic_selector")),
            uiOutput(ns("metric_selector"))
          ),
          control_group("Donors",
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
                               inline = TRUE),
            uiOutput(ns("gender_filter_ui")),
            uiOutput(ns("age_filter_ui"))
          ),
          control_group("Object size",
            sliderInput(ns("diam_range"), "Islet diameter range (\u00b5m)",
                        min = 0, max = 500, value = c(0, 350), step = 10),
            sliderInput(ns("binwidth"), "Diameter bin width (\u00b5m)",
                        min = 1, max = 75, value = 50, step = 1)
          ),
          control_group("Exclusions",
            checkboxInput(ns("exclude_zero_top"), "Exclude zero values (left plot)", value = FALSE),
            checkboxInput(ns("exclude_zero_dist"), "Exclude zero values (right plot)", value = FALSE),
            # Non-namespaced: read by app.R as the app-wide outlier setting.
            tags$div(
              title = GLOBAL_OUTLIER_TOOLTIP,
              checkboxInput("sidebar_remove_outliers", "Exclude outliers", value = FALSE)
            ),
            sliderInput("sidebar_outlier_z", "Outlier z-threshold",
                        min = 0.5, max = 10, value = 3, step = 0.5),
            checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE)
          )
        )
      ),
      control_panel("View controls",
        help = list(
          tags$p("View controls change how the plots are drawn, not what they contain. The groups marked \"both",
                 " plots\" apply to both cards, whatever they plot; the Against diameter and By donor status groups",
                 " apply to those plot types (Log2 scale x-axis also to Density against diameter). Click the panel's",
                 " title to open or close it."),
          tags$p("Use browser zoom (Ctrl +/\u2212) to resize cards.")),
        # Point controls apply to both plots (owner, 2026-10-19: the View controls must apply to both,
        # not only one card); the last two groups apply to the plot type they name.
        control_grid(
          control_group("Colour (both plots)",
            selectInput(ns("plot_color_by"), "Color points by",
                        choices = c("Donor Status" = "donor_status", "Donor ID" = "donor_id"),
                        selected = "donor_status"),
            # Non-namespaced: global donor palette, synced with Trajectory and Spatial (app.R).
            selectInput("sidebar_donor_palette", "Donor status colors",
                        choices = names(DONOR_PALETTES), selected = "Bright")
          ),
          control_group("Points and axis (both plots)",
            checkboxInput(ns("show_points"), "Show individual points", value = TRUE),
            checkboxInput(ns("log_scale"), "Log scale y-axis", value = FALSE),
            sliderInput(ns("pt_size"), "Point size", min = 0.3, max = 4.0, value = 3.0, step = 0.1),
            sliderInput(ns("pt_alpha"), "Point transparency", min = 0, max = 0.95, value = 0.4, step = 0.05)
          ),
          control_group("Against diameter",
            radioButtons(ns("stat"), "Summary statistic",
                         choices = c("Mean\u00b1SE" = "mean_se",
                                     "Mean\u00b1SD" = "mean_sd",
                                     "Median + IQR" = "median_iqr"),
                         selected = "mean_se", inline = TRUE),
            radioButtons(ns("add_smooth"), "Trend line", choices = c("None", "LOESS"),
                         selected = "None", inline = TRUE),
            checkboxInput(ns("log_scale_x"), "Log2 scale x-axis", value = FALSE)
          ),
          control_group("By donor status",
            radioButtons(ns("dist_type"), "Plot type", choices = c("Violin", "Box"),
                         selected = "Violin", inline = TRUE)
          )
        )
      ),
      div(class = "tab-controls-actions",
        actionButton(ns("clear_link_selection"), "Clear selection", class = "btn-sm btn-outline-secondary"))
    ),
    # The excluded-outlier table, when asked for, sits under the frozen controls rather than inside them.
    column(12, uiOutput(ns("plot_outlier_info"))),
    # The two plot cards: each plots any plot type (R/plot_types.R, owner 2026-10-21), opening on the two it always
    # showed. The right card widens to fill the row when no extra panel (e.g. the AI chat card) is present.
    column(6, plot_card_ui(ns, "left", "plt", "dl_summary")),
    column(if (is.null(extra_panel)) 6 else 4, plot_card_ui(ns, "right", "dist", "dl_distribution_csv")),
    # Extra panel (e.g., AI chat) injected between cards and full-width panels
    extra_panel,
    # Full-width card below both plot cards: segmentation/drill-down viewer
    column(12,
      uiOutput(ns("segmentation_viewer_panel"))
    ),
    # Every donor's whole cell make-up (R/donor_composition.R).
    column(12, donor_composition_card_ui(ns)),
    # Per-donor counts of each endocrine object class (uniform root only): the one place Singles are shown.
    if (!is.null(load_object_counts())) column(12, object_counts_card_ui(ns))
  )
}

#' One plot card: its title is the chosen plot type; the Plot type selector, the CSV button and, for the types that
#' need one, the x measure or clinical value selector sit on the card.
plot_card_ui <- function(ns, card, plot_id, dl_id) {
  type_id <- paste0(card, "_type")
  default <- if (card == "left") "diameter" else "status"
  side <- if (card == "left") "Left" else "Right"
  div(class = "card plot-type-card",
      style = paste0("margin-bottom: 20px; padding: 15px; border: 1px solid var(--islet-border-strong);",
                     " border-radius: 8px; height: clamp(600px, 75vh, 1400px); display: flex; flex-direction: column;",
                     " gap: 8px; overflow: visible;"),
    div(style = "display: flex; align-items: center; justify-content: space-between; gap: 10px; flex-wrap: wrap;",
      div(style = "display: flex; align-items: center; gap: 4px; min-width: 0;",
        h5(textOutput(ns(paste0(card, "_title")), inline = TRUE), style = "margin: 0; color: var(--islet-fg);"),
        do.call(card_help, c(list(paste(side, "plot help")), list(plot_types_help())))),
      div(style = "display: flex; align-items: center; gap: 8px;",
        div(class = "plot-type-select",
            selectInput(ns(type_id), "Plot type", choices = PLOT_TYPE_CHOICES, selected = default, width = "360px")),
        downloadButton(ns(dl_id), "CSV", class = "btn-sm btn-outline-secondary",
                       style = "font-size: 11px; padding: 2px 8px;"))
    ),
    conditionalPanel(sprintf("input['%s'] == 'xy'", ns(type_id)),
      div(class = "plot-type-select",
          selectInput(ns(paste0(card, "_x_measure")), "X axis", choices = NULL, width = "380px"))),
    conditionalPanel(sprintf("input['%s'] == 'clinical'", ns(type_id)),
      div(class = "plot-type-select",
          selectInput(ns(paste0(card, "_clinical")), "Clinical value", choices = NULL, width = "380px"))),
    div(style = "flex: 1; min-height: 0; display: flex;",
      plotlyOutput(ns(plot_id), height = "calc(100% - 28px)")
    )
  )
}
