# ---------- Inflammation tab module UI ----------
# Exports: inflammation_ui(id)
#
# Cohort-level analyses ported from Steenbuck et al. 2025 (Bodenmiller lab,
# bioRxiv 2025.03.05.641526). These do not fit the Statistics tab's flow, which
# tests ONE selected feature at a time: both cards here test every cell type at
# once and answer a cohort question rather than a per-islet one.
#
#   Card 1  Cell-type composition  donor-level tests on centred log-ratios
#   Card 2  Infiltration score     distance-weighted immune localisation
#
# The donor-cohort table (nPOD clinical metadata) and the edgeR count model were
# both removed Jul 2026 -- the first at the user's request, the second because
# its calibrations assume a transcriptome. See the Methods card and
# mod_inflammation_server.R's comp_fit().
#
# Layout note: every plot here sits in a col-6 or wider container. The Statistics
# tab's col-4 cards clipped their titles and axis labels (a 453 px title in a
# 401 px plot), so plots in this module carry no in-plot title, set `margin`
# explicitly rather than trusting `automargin`, and use `type = "category"` on
# any axis whose labels are numeric-looking strings.

#' Selection block for the frozen right-hand panel.
#'
#' This tab reads the tissue-wide single-cell table, not the islet frame, so its
#' selectors are about which CELLS enter each card. The explanatory notes that
#' sat under each control come with them.
#'
#' The cohort filter reaches this tab as a donor keep-list only: groups, age and
#' sex resolve against the islet table and land as an id membership test, while
#' min-cells and diameter are islet-level and are reported as not-applied rather
#' than silently ignored.
inflammation_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Scope is fixed to the endocrine and immune specific types; the
    # granularity selector it replaced is gone with it.
    uiOutput(ns("comp_level_note")),
    ie_field("Composition: cells included",
      selectInput(ns("edger_region"), NULL,
                  choices = c("Islet core + peri-islet" = "core_peri",
                              "Islet core only" = "core",
                              "All tissue" = "all"),
                  selected = "core_peri", width = "100%")),
    ie_field("Infiltration: denominator",
      selectInput(ns("infil_denom"), NULL,
                  choices = c("Islet + peri-islet cells" = "isletperi",
                              "Cells within 100 µm of an islet" = "near100",
                              "All cells on the section" = "all"),
                  selected = "isletperi", width = "100%")),
    tags$small(style = "display: block; margin: -4px 0 10px; line-height: 1.35",
      "The ROI analogue. Steenbuck et al. normalise per islet-centred ROI; ",
      "normalising by the whole section instead adds 4–23× donor-to-donor ",
      "variance and hid three phenotypes that are significant under the default."),
    checkboxInput(ns("infil_match_numer"), "Match numerator to denominator", value = FALSE),
    tags$small(style = "display: block; margin: -6px 0 10px; line-height: 1.35",
      "Restores the paper's exact per-cell mean. Off by default because matching ",
      "confines the 1/d gradient to the 20 µm peri band."),
    checkboxInput(ns("infil_log"), "Log10 y-axis", value = TRUE),
    uiOutput(ns("infil_pheno_ui"))
  )
}

inflammation_ui <- function(id) {
  ns <- NS(id)

  # `section_heading()` was defined locally here (and in three other
  # files) with the same hardcoded blue gradient inline. Replaced by
  # ie_section_heading() in R/ui_modernist_bits.R.
  section_heading <- ie_section_heading

  tagList(


    # Banner text moved to the rail's "Tips for this tab".

    # ===== SECTION 1: cell-type composition =====
    section_heading("1", "Cell-Type Composition",
                    "Which cell types shift between disease stages?",
                    tip = paste0("Per donor, the cell-type vector is closed to a composition ",
                                 "and centred-log-ratio transformed; each type is then tested ",
                                 "across stage on one value per donor, with Benjamini-Hochberg ",
                                 "across cell types.")),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible",
          # "Cells included" moved to inflammation_selection_ui() in the frozen
          # right-hand panel. The pairwise Contrast selector is gone entirely:
          # this is one omnibus across every stage present, so there is no
          # contrast to pick -- the table reports the log-ratio difference
          # between the extreme stages.
          div(style = "margin-bottom: 10px",
            downloadButton(ns("dl_edger"), "Download results CSV")
          ),
          uiOutput(ns("edger_note")),
          fluidRow(style = "display: flex; flex-wrap: wrap",
            column(6,
              plotlyOutput(ns("edger_volcano"), height = "460px")
            ),
            column(6,
              div(style = "max-height: 460px; overflow-y: auto",
                tableOutput(ns("edger_table"))
              )
            )
          )
        )
      )
    ),

    # ===== SECTION 2: Infiltration score =====
    section_heading("2", "Infiltration Score",
                    "Distance-weighted immune localisation relative to the islet edge.",
                    tip = paste0("Score = 1 for cells inside an islet and 1/d for cells ",
                                 "outside it (d = µm from the islet edge), normalised by ",
                                 "the donor's total cell count. Rewards peri-islet ",
                                 "accumulation, penalises distal accumulation. FDR runs ",
                                 "across all phenotypes present, not the display selection.")),

    # The card body was removed by an over-broad regex during the Jul 2026 prose
    # strip -- it took the plot, table, note and export along with the
    # explanatory box it was meant to remove, leaving a section heading with
    # nothing under it. Rebuilt to match the composition card above; the
    # denominator / numerator-matching / log-axis / phenotype controls stay in
    # the frozen right panel where they were moved to.
    fluidRow(
      column(12,
        div(class = "card",
          div(style = "margin-bottom: 10px",
            downloadButton(ns("dl_infil"), "Download scores CSV")
          ),
          uiOutput(ns("infil_note")),
          fluidRow(style = "display: flex; flex-wrap: wrap",
            column(7, plotlyOutput(ns("infil_plot"), height = "460px")),
            column(5,
              div(style = "max-height: 460px; overflow-y: auto",
                tableOutput(ns("infil_table"))
              )
            )
          )
        )
      )
    ),

    # ===== SECTION 3: Methods =====
    section_heading("3", "Methods Reference",
                    "Definitions, corrections and the adaptations made for this cohort."),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px",
          uiOutput(ns("inflam_methods"))
        )
      )
    )
  )
}
