# ---------------------------------------------------------------------------
# Shared cohort filter
#
# Mounted ONCE, above the pipeline tabs. Every analysis panel reads the
# reactives this returns — no panel owns its own donor filter any more.
# This is the change that makes the redesign feel like one app; do it before
# any of the CSS work.
#
#   cohort <- cohortFilterServer("cohort", islets, feature_cols = MARKERS)
#   output$umap <- renderPlot(plot_umap(cohort$islets(), cohort$point_size()))
# ---------------------------------------------------------------------------

GROUP_LEVELS <- c("ND", "Aab+", "T1D")

# label + control cell used by the cohort bar and the per-panel control rows
ie_field <- function(label, ..., width = NULL) {
  div(
    class = "ie-field",
    style = if (!is.null(width)) paste0("width:", width),
    tags$span(class = "ie-field-label", label),
    ...
  )
}

# Hand-rolled segmented control. Shiny binds any radio inputs sharing a name
# inside a .shiny-input-radiogroup container, so we get full markup control
# (and pure-CSS active state via :has(input:checked)) with no extra package.
ie_segmented <- function(inputId, choices, selected = choices[[1]]) {
  values <- if (is.null(names(choices))) choices else unname(choices)
  labels <- if (is.null(names(choices))) choices else names(choices)
  div(
    id = inputId, class = "shiny-input-radiogroup ie-seg",
    Map(function(v, l) {
      tags$label(
        class = "ie-seg-opt",
        tags$input(
          type = "radio", name = inputId, value = v,
          checked = if (identical(v, selected)) "checked"
        ),
        tags$span(l)
      )
    }, values, labels)
  )
}

cohortFilterUI <- function(id, age_range = c(0, 70)) {
  ns <- NS(id)
  div(
    class = "ie-cohortbar",

    ie_field(
      "Disease groups",
      # per-chip swatch colors are applied by nth-child in modernist-shiny.css,
      # so the order of GROUP_LEVELS is load-bearing
      checkboxGroupInput(
        ns("groups"), NULL,
        choices = GROUP_LEVELS, selected = GROUP_LEVELS, inline = TRUE
      )
    ),

    ie_field(
      "Donor age",
      sliderInput(
        ns("age"), NULL,
        min = age_range[1], max = age_range[2],
        value = age_range, step = 1, post = " yr",
        ticks = FALSE, width = "210px"
      )
    ),

    ie_field("Sex", ie_segmented(
      ns("sex"), c(All = "all", Female = "F", Male = "M"), "all"
    )),

    ie_field("Normalization", ie_segmented(
      ns("norm"), c(`Robust z` = "robust", `Std z` = "z", None = "none"), "robust"
    )),

    ie_field("Point size", sliderInput(
      ns("pt_size"), NULL, min = 0.4, max = 3, value = 1.2, step = 0.1,
      ticks = FALSE, width = "130px"
    )),

    ie_field("Opacity", sliderInput(
      ns("pt_alpha"), NULL, min = 0.1, max = 1, value = 0.82, step = 0.05,
      ticks = FALSE, width = "130px"
    )),

    div(class = "ie-cohort-readout", textOutput(ns("readout"), inline = TRUE))
  )
}

#' @param islets data frame or reactive returning one. Needs columns
#'   `status` (in GROUP_LEVELS), `donor_id`, `donor_age`, `donor_sex`.
#' @param feature_cols marker columns to renormalize when `norm` changes.
cohortFilterServer <- function(id, islets, feature_cols = character()) {
  moduleServer(id, function(input, output, session) {

    # debounced so dragging a slider doesn't queue a redraw per pixel
    dat <- debounce(reactive({
      d <- if (is.reactive(islets)) islets() else islets
      req(d, input$groups, input$age)

      d <- d[d$status %in% input$groups, , drop = FALSE]
      d <- d[d$donor_age >= input$age[1] & d$donor_age <= input$age[2], , drop = FALSE]
      if (!identical(input$sex, "all")) {
        d <- d[d$donor_sex == input$sex, , drop = FALSE]
      }
      normalize_features(d, intersect(feature_cols, names(d)), input$norm)
    }), 150)

    output$readout <- renderText({
      d <- dat()
      sprintf(
        "%s islets · %d donors",
        format(nrow(d), big.mark = ","), length(unique(d$donor_id))
      )
    })

    list(
      islets      = dat,
      groups      = reactive(input$groups),
      point_size  = reactive(input$pt_size),
      point_alpha = reactive(input$pt_alpha),
      norm        = reactive(input$norm)
    )
  })
}

normalize_features <- function(df, cols, method = "robust") {
  if (identical(method, "none") || !length(cols) || !nrow(df)) return(df)
  for (cl in cols) {
    x <- df[[cl]]
    df[[cl]] <- if (identical(method, "robust")) {
      m <- stats::median(x, na.rm = TRUE); s <- stats::mad(x, na.rm = TRUE)
      if (is.na(s) || s == 0) x - m else (x - m) / s
    } else {
      m <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) x - m else (x - m) / s
    }
  }
  df
}
