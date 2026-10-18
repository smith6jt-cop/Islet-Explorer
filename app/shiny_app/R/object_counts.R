# ---------- Endocrine objects per donor ----------
# Counts of every endocrine object class per donor, from a uniform root's endocrine_objects_by_donor.csv
# (scripts/uniform/assemble_uniform.py). This card is the only place Singles enter the app: they are not analysis
# objects (owner, 2026-10-18: "Singles should not be included in any islet-level analysis except for quantifying
# or showing how many there are"). Shown in the Plot tab; absent on a root without the file.

OBJECT_COUNT_MEASURES <- c(
  "Islets per mm² of tissue" = "islets_per_mm2",
  "Clusters per mm² of tissue" = "clusters_per_mm2",
  "Singles per mm² of tissue" = "singles_per_mm2",
  "Islets (count)" = "islets",
  "Clusters (count)" = "clusters",
  "Singles (count)" = "singles",
  "Singles by cell type (count)" = "singles_by_type")

# The Single's one cell: its typed label, Excluded (strongest hormone under 2x the next), or not analysed (no nucleus;
# computed in load_object_counts()). A Single's cell is never left as Endocrine without a subtype.
SINGLE_TYPE_COLUMNS <- c("Beta cell" = "singles_Beta_cell", "Alpha cell" = "singles_Alpha_cell",
                         "Delta cell" = "singles_Delta_cell", "Excluded (hormone unresolved)" = "singles_excluded",
                         "Not analysed (no nucleus)" = "singles_not_analysed")

load_object_counts <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    path <- data_path("endocrine_objects_by_donor.csv")
    if (!file.exists(path)) return(NULL)
    d <- tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(d)) return(NULL)
    d$case_id <- as.character(d$case_id)
    typed <- c("singles_Beta_cell", "singles_Alpha_cell", "singles_Delta_cell", "singles_Endocrine",
               "singles_excluded")
    d$singles_not_analysed <- d$singles - rowSums(d[, intersect(typed, names(d)), drop = FALSE])
    d$donor_status <- factor(d$donor_status, levels = c("ND", "Aab+", "T1D"))
    cache <<- d
    d
  }
})

object_counts_card_ui <- function(ns) {
  div(class = "card",
      style = "margin-bottom: 20px; padding: 15px; border: 1px solid var(--islet-border-strong); border-radius: 8px;",
    div(style = "display: flex; align-items: center; justify-content: space-between; margin-top: 0;",
      h5("Endocrine objects per donor", style = "margin: 0; color: var(--islet-fg);"),
      downloadButton(ns("dl_object_counts"), "CSV", class = "btn-sm btn-outline-secondary",
                     style = "font-size: 11px; padding: 2px 8px;")
    ),
    selectInput(ns("object_count_measure"), "Show", choices = OBJECT_COUNT_MEASURES,
                selected = "singles_per_mm2", width = "320px"),
    plotlyOutput(ns("object_counts_plot"), height = "460px")
  )
}

object_counts_plot <- function(d, measure, colors) {
  fonts <- plotly_islet_fonts()
  label <- names(OBJECT_COUNT_MEASURES)[match(measure, OBJECT_COUNT_MEASURES)]
  if (identical(measure, "singles_by_type")) {
    d <- d[order(d$donor_status, d$case_id), , drop = FALSE]
    donors <- sprintf("%s (%s)", d$case_id, d$donor_status)
    p <- plotly::plot_ly()
    pal <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#7f7f7f")
    for (i in seq_along(SINGLE_TYPE_COLUMNS)) {
      p <- plotly::add_bars(p, x = donors, y = d[[SINGLE_TYPE_COLUMNS[i]]], name = names(SINGLE_TYPE_COLUMNS)[i],
                            marker = list(color = pal[i]),
                            hovertemplate = paste0("%{x}<br>", names(SINGLE_TYPE_COLUMNS)[i], ": %{y}<extra></extra>"))
    }
    p <- plotly::layout(p, barmode = "stack",
                        xaxis = list(title = list(text = "Donor (status)", font = fonts$axis_title),
                                     tickfont = fonts$axis_tick, categoryorder = "array", categoryarray = donors,
                                     tickangle = -45, automargin = TRUE),
                        yaxis = list(title = list(text = "Singles", font = fonts$axis_title),
                                     tickfont = fonts$axis_tick),
                        legend = list(font = fonts$legend_text))
  } else {
    set.seed(1)
    x <- as.integer(d$donor_status) + stats::runif(nrow(d), -0.12, 0.12)
    med <- as.numeric(tapply(d[[measure]], d$donor_status, stats::median))
    p <- plotly::plot_ly()
    for (s in levels(d$donor_status)) {
      k <- d$donor_status == s
      if (!any(k)) next
      p <- plotly::add_markers(p, x = x[k], y = d[[measure]][k], name = s,
                               marker = list(size = 13, color = colors[[s]] %||% "#777777",
                                             line = list(color = ISLET_THEME_PALETTE[[islet_theme_mode()]]$paper, width = 1)),
                               text = sprintf("Donor %s (%s)<br>%s: %s", d$case_id[k], s, label,
                                              format(signif(d[[measure]][k], 3), big.mark = ",")),
                               hoverinfo = "text")
    }
    lv <- seq_along(levels(d$donor_status))
    has <- is.finite(med)
    p <- plotly::add_segments(p, x = lv[has] - 0.3, xend = lv[has] + 0.3, y = med[has], yend = med[has],
                              name = "Donor median", line = list(color = islet_ink(), width = 3),
                              hoverinfo = "skip", inherit = FALSE)
    p <- plotly::layout(p,
                        xaxis = list(title = list(text = "Donor status", font = fonts$axis_title),
                                     tickmode = "array", tickvals = lv, ticktext = levels(d$donor_status),
                                     tickfont = fonts$axis_tick, range = c(0.4, length(lv) + 0.6)),
                        yaxis = list(title = list(text = label, font = fonts$axis_title),
                                     tickfont = fonts$axis_tick, rangemode = "tozero"),
                        legend = list(font = fonts$legend_text))
  }
  p %>% plotly::layout(font = fonts$global) %>% plotly_islet_config("endocrine_objects_per_donor")
}
