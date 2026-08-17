# ---------------------------------------------------------------------------
# Left rail + hidden tabset
#
# The rail is a hand-rolled Shiny radio group: one <input type=radio> per
# destination inside a .shiny-input-radiogroup container. That buys us
#   - arbitrary markup per row (gold numeral, subtitle, section headings)
#   - active state in pure CSS (.ie-nav:has(input:checked)) — no re-render,
#     no JS message handler, no actionButton counter bookkeeping
#
#   PIPELINE <- list(list(key="plot", num="01", label="Plot"), ...)
#   railUI("rail", PIPELINE, REFERENCE)
#   tab <- railServer("rail", "panels")   # keeps navset_hidden in sync
# ---------------------------------------------------------------------------

rail_row <- function(group, key, selected = FALSE, num = NULL,
                     label = NULL, sub = NULL, marker = NULL) {
  tags$label(
    class = "ie-nav",
    tags$input(
      type = "radio", name = group, value = key,
      checked = if (selected) "checked"
    ),
    tags$span(class = "ie-nav-num", num %||% marker),
    tags$span(
      class = "ie-nav-text",
      tags$span(label),
      if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub)
    )
  )
}

rail_heading <- function(text) div(class = "ie-rail-heading", text)

railUI <- function(id, pipeline, reference = list(),
                   title = c("ISLET", "EXPLORER"), meta = NULL, footer = NULL) {
  ns <- NS(id)
  first <- pipeline[[1]]$key

  tags$nav(
    class = "ie-rail",
    div(
      class = "ie-rail-brand",
      div(
        class = "ie-brandmark",
        tags$span(class = "ie-brandbar"),
        div(class = "ie-brandtext", HTML(paste(title, collapse = "<br>")))
      ),
      if (!is.null(meta)) div(class = "ie-rail-meta", HTML(meta))
    ),

    # single radio group across both sections, so exactly one row is ever lit
    div(
      id = ns("tab"), class = "shiny-input-radiogroup",

      rail_heading("Analysis pipeline"),
      lapply(pipeline, function(it) {
        rail_row(ns("tab"), it$key, identical(it$key, first),
                 num = it$num, label = it$label, sub = it$sub)
      }),

      if (length(reference)) tagList(
        div(class = "ie-rail-divider"),
        rail_heading("Reference · read-only"),
        lapply(reference, function(it) {
          rail_row(ns("tab"), it$key, FALSE,
                   marker = "\u25C7", label = it$label, sub = it$sub)
        })
      )
    ),

    if (!is.null(footer)) div(class = "ie-rail-footer", footer)
  )
}

rail_palette_note <- function(label, value) {
  div(
    div(class = "ie-rail-footer-label", label),
    div(class = "ie-rail-footer-value", value)
  )
}

#' @param tabset_id outer id of the bslib::navset_hidden holding the panels
#' @return reactive giving the active key
railServer <- function(id, tabset_id) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$tab, {
      # nav_select needs the outer namespace, hence session$rootScope()
      bslib::nav_select(tabset_id, selected = input$tab,
                        session = session$rootScope())
    })
    reactive(input$tab)
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
