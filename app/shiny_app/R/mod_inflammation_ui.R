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

inflammation_ui <- function(id) {
  ns <- NS(id)

  section_heading <- function(step, title, subtitle, tip = NULL) {
    tip_icon <- if (!is.null(tip) && nzchar(tip)) {
      span("?", title = tip,
           style = paste0("display: inline-block; margin-left: 6px; width: 18px; height: 18px;",
                          " line-height: 18px; text-align: center; font-size: 12px; font-weight: 700;",
                          " color: #4477AA; background: #eaf2fa; border: 1px solid #c0d4e8;",
                          " border-radius: 50%; cursor: help;"))
    } else NULL
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
      div(style = "display: flex; align-items: baseline; gap: 10px;",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 24px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 17px; color: #333;"),
        tip_icon
      ),
      p(subtitle, style = "margin: 4px 0 0 0; color: #777; font-size: 13px;")
    )
  }

  tagList(

    phenotype_rules_button("show_phenotype_rules_inflam"),

    # ===== Banner =====
    div(class = "card",
        style = paste0("padding: 14px 18px; margin-bottom: 18px; ",
                       "background: linear-gradient(135deg, #f4f8fc 0%, #eaf2fa 100%); ",
                       "border-left: 4px solid #4477AA;"),
      h4("About this tab", style = "margin: 0 0 6px 0; font-size: 16px; color: #2c4760;"),
      p(style = "margin: 0 0 6px 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("What it does."), " Two cohort-level analyses adapted from ",
        tags$em("Steenbuck et al. 2025"), " (Bodenmiller lab, bioRxiv 2025.03.05.641526), ",
        "an imaging-mass-cytometry study of 88 nPOD donors across T1D progression. ",
        "Unlike the Statistics tab, which tests one selected feature, both cards here ",
        "test ", tags$strong("every cell type at once"), " and are corrected across them."),
      p(style = "margin: 0 0 6px 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("Why it is separate."), " These are donor-level questions — ",
        "“which cell types change in abundance?” and “which immune cells move ",
        "toward islets?” — not per-islet ones. They read the tissue-wide single-cell ",
        "table (~2.65 M cells), not the islet-level frame the Statistics tab uses, so ",
        "the Plot sidebar's filters do not apply here."),
      p(style = "margin: 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("Cohort scope."), " All five T1D donors in this cohort are ",
        tags$strong("recent-onset"), " (disease duration < 1 year). There are no ",
        "long-duration donors, so nothing here speaks to established disease.")
    ),

    # ===== SECTION 1: cell-type composition =====
    section_heading("1", "Cell-Type Composition",
                    "Which cell types shift between disease stages?",
                    tip = paste0("Per donor, the cell-type vector is closed to a composition ",
                                 "and centred-log-ratio transformed; each type is then tested ",
                                 "across stage on one value per donor, with Benjamini-Hochberg ",
                                 "across cell types.")),
    div(style = paste0("background: #eef6ee; border-left: 4px solid #2e7d32; border-radius: 4px; ",
                       "padding: 10px 14px; margin-bottom: 12px; font-size: 13px; line-height: 1.5; color: #1b4620;"),
      tags$strong("Why the log-ratio, not raw proportions or a count model. "),
      "Proportions are compositionally constrained — they must sum to 1, so the loss of β cells ",
      "in T1D mechanically inflates every other cell type's share even when that type has not ",
      "moved. The centred log-ratio removes that: each value is a cell type's abundance ",
      "relative to that donor's own geometric mean, which is also invariant to how much tissue ",
      "was imaged. ",
      tags$strong("This replaced an edgeR count model in Jul 2026"), " — the counts were real, ",
      "but edgeR's dispersion moderation, TMM normalisation and library-size offset are all ",
      "calibrated for a 20,000-gene transcriptome and do not transfer to 21 cell types from 15 ",
      "donors. See Methods."),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          fluidRow(
            column(5,
              selectInput(ns("edger_region"), "Cells included",
                          choices = c("Islet core + peri-islet" = "core_peri",
                                      "Islet core only" = "core",
                                      "All tissue" = "all"),
                          selected = "core_peri")
            ),
            # The pairwise Contrast selector is gone: this is a single omnibus
            # across every stage present, so there is no contrast to pick. The
            # table reports the log-ratio difference between the extreme stages.
            column(5, style = "padding-top: 25px;",
              downloadButton(ns("dl_edger"), "Download results CSV")
            )
          ),
          uiOutput(ns("edger_note")),
          fluidRow(style = "display: flex; flex-wrap: wrap;",
            column(6,
              plotlyOutput(ns("edger_volcano"), height = "460px")
            ),
            column(6,
              div(style = "max-height: 460px; overflow-y: auto;",
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
    div(style = paste0("background: #f4f8fc; border-left: 4px solid #4477AA; border-radius: 4px; ",
                       "padding: 10px 14px; margin-bottom: 12px; font-size: 13px; line-height: 1.5; color: #2c4760;"),
      tags$strong("Definition. "),
      "Every cell scores ", tags$code("1"), " if it lies inside an islet and ",
      tags$code("1/d"), " if outside, with ", tags$em("d"),
      " the distance in µm from the islet edge; the per-donor sum is divided by that ",
      "donor's total cell count. Steenbuck et al. normalise per imaging ROI — this ",
      "cohort has whole-section CODEX images rather than ROIs, so the ROI analogue \u2014 ",
      "islet + peri-islet cells \u2014 is the default denominator here. Distances come from ",
      tags$code("data/cell_distance.csv"), " (negative = inside the islet). ",
      tags$strong("Multiple testing is corrected within a pre-specified family of 7 immune ",
                  "phenotypes"), ", fixed in source so it cannot be resized after seeing results."),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
          fluidRow(
            column(5,
              selectInput(ns("infil_denom"), "Normalise by (denominator)",
                          choices = c("Islet + peri-islet cells" = "isletperi",
                                      "Cells within 100 µm of an islet" = "near100",
                                      "All cells on the section" = "all"),
                          selected = "isletperi"),
              tags$small(style = "color: #888; display: block; margin-top: -10px; line-height: 1.35;",
                "The ROI analogue. Steenbuck et al. normalise per islet-centred ROI; ",
                "normalising by the whole section instead adds 4\u201323\u00d7 donor-to-donor ",
                "variance and hid three phenotypes that are significant under the default.")
            ),
            column(4,
              uiOutput(ns("infil_pheno_ui"))
            ),
            column(3, style = "padding-top: 25px;",
              checkboxInput(ns("infil_match_numer"),
                            "Match numerator to denominator", value = FALSE),
              tags$small(style = "color: #888; display: block; margin-top: -8px; line-height: 1.35;",
                "Restores the paper's exact per-cell mean. Off by default because matching ",
                "confines the 1/d gradient to the 20 µm peri band."),
              checkboxInput(ns("infil_log"), "Log10 y-axis", value = TRUE),
              downloadButton(ns("dl_infil"), "Download scores CSV")
            )
          ),
          fluidRow(style = "display: flex; flex-wrap: wrap;",
            column(7,
              plotlyOutput(ns("infil_plot"), height = "460px")
            ),
            column(5,
              div(style = "max-height: 460px; overflow-y: auto;",
                tableOutput(ns("infil_table"))
              )
            )
          ),
          uiOutput(ns("infil_note"))
        )
      )
    ),

    # ===== SECTION 3: Methods =====
    section_heading("3", "Methods Reference",
                    "Definitions, corrections and the adaptations made for this cohort."),
    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
          uiOutput(ns("inflam_methods"))
        )
      )
    )
  )
}
