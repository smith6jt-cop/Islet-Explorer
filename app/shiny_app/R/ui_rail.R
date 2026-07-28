# ---------------------------------------------------------------------------
# Left rail + hidden tabset
#
# Ported from design/delivered/shiny/R/ui_rail.R with two changes, both noted
# in place: the `%||%` redefinition is dropped, and rail_row() wraps its visible
# content in a span so the active state can be expressed with a plain sibling
# combinator instead of :has().
#
# The rail is a hand-rolled Shiny radio group: one <input type="radio"> per
# destination inside a .shiny-input-radiogroup container. Shiny's
# RadioInputBinding finds any such container and reads
# `input:radio[name=<container id>]:checked` -- a document-wide name lookup --
# so this binds with no extra package and no custom input binding. That buys:
#   - arbitrary markup per row (numeral, label, subtitle, section headings)
#   - active state in pure CSS, so switching tabs costs no re-render, no custom
#     message handler and no actionButton counter bookkeeping
#
#   railUI("rail", IE_PIPELINE, IE_REFERENCE)
#   railServer("rail", "panels")   # keeps navset_hidden in sync
# ---------------------------------------------------------------------------

rail_row <- function(group, key, selected = FALSE, num = NULL,
                     label = NULL, sub = NULL, marker = NULL) {
  tags$label(
    class = "ie-nav",
    tags$input(
      type = "radio", name = group, value = key,
      checked = if (selected) "checked"
    ),
    # WRAPPER SPAN, not in the delivered file. The original styled `.ie-nav`
    # directly and selected the active row with `.ie-nav:has(input:checked)`.
    # :has() needs Firefox 121+; on Firefox ESR 115 -- which shipped on managed
    # university desktops through Sep 2024 -- the rail would render with NO
    # active-tab indicator at all while still working, and nothing would log.
    # With the content in a sibling of the input, `input:checked ~ .ie-nav-body`
    # does the same job with no browser gate.
    tags$span(
      class = "ie-nav-body",
      tags$span(class = "ie-nav-num", num %||% marker),
      tags$span(
        class = "ie-nav-text",
        tags$span(label),
        if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub)
      )
    )
  )
}

rail_heading <- function(text) div(class = "ie-rail-heading", text)

railUI <- function(id, pipeline, reference = list(),
                   brand = "Panc Floyd", logo = "logo.png",
                   meta = NULL, help = NULL, footer = NULL) {
  ns <- NS(id)
  first <- pipeline[[1]]$key

  tags$nav(
    class = "ie-rail",
    div(
      class = "ie-rail-brand",
      div(
        class = "ie-brandmark",
        if (!is.null(logo))
          tags$img(class = "ie-brandlogo", src = logo, alt = paste(brand, "logo")),
        div(class = "ie-brandtext", brand)
      ),
      if (!is.null(meta)) div(class = "ie-rail-meta", HTML(meta))
    ),

    # A single radio group spans BOTH sections, so exactly one row is ever lit
    # -- the Viewer is a destination like any other, it is only grouped apart.
    div(
      id = ns("tab"), class = "shiny-input-radiogroup",
      role = "radiogroup", `aria-label` = "Analysis pipeline",

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
                   marker = "◇", label = it$label, sub = it$sub)
        })
      )
    ),

    # Help lives here rather than at the top of each tab. Six copies of the
    # phenotype-rules pill and a row of blue hint strips were taking the first
    # 120px of every panel; one place, always reachable, no vertical cost.
    if (!is.null(help)) tagList(
      div(class = "ie-rail-divider"),
      rail_heading("Help"),
      div(class = "ie-rail-help", help)
    ),

    if (!is.null(footer)) div(class = "ie-rail-footer", footer)
  )
}

#' A help entry in the rail. An actionLink, so the server opens a modal.
rail_help_link <- function(inputId, label, sub = NULL) {
  actionLink(inputId, class = "ie-rail-helplink",
    label = tagList(
      tags$span(class = "ie-nav-num", "?"),
      tags$span(class = "ie-nav-text",
        tags$span(label),
        if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub))
    ))
}

#' A labelled control in the rail footer.
#'
#' The delivered `rail_palette_note()` rendered palette state as static TEXT.
#' These are the real selectors -- making the rail the only writer of either
#' palette is what lets app.R's sync-observer triads go away.
rail_palette_select <- function(inputId, label, choices, selected) {
  div(
    class = "ie-rail-footer-field",
    tags$label(class = "ie-rail-footer-label", `for` = inputId, label),
    selectInput(inputId, NULL, choices = choices, selected = selected,
                width = "100%")
  )
}

#' Static label/value pair for the rail footer (kept for non-interactive notes).
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
      # rootScope() is REQUIRED and verified: without it nav_select resolves
      # `tabset_id` inside this module's namespace and silently does nothing.
      bslib::nav_select(tabset_id, selected = input$tab,
                        session = session$rootScope())
    }, ignoreInit = TRUE)
    reactive(input$tab)
  })
}

# NOTE: the delivered file ended with a `%||%` definition. Deliberately dropped.
# This app has never defined `%||%` -- it uses base R's (>= 4.4), and
# tests/testthat/helper-setup.R says so. Shiny source()s every R/*.R into the
# GLOBAL environment, which sits ahead of package:base on the search path, so a
# copy here would shadow base's for the whole app. Semantically identical today,
# but mod_statistics_server.R depends on `%||%` NOT catching NA (an emptied
# numericInput sends NA, and treating it as missing re-breaks a fixed crash),
# and that dependency is invisible from a UI file.
