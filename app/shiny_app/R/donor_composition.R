# ---------- Plot tab: Composition by donor card (owner, 2026-10-21) ----------
# The plot cards show one phenotype at a time; this card shows every donor's whole make-up: the proportion of each
# retained phenotype among all the cells of the donor's selected objects, pooled (cells of the phenotype summed over
# the objects, over all their cells), in the core, the 20 um peri-islet band or both. Cells of the general classes
# (Endocrine, T cell, Immune: cells without a subtype) count in the total and are shown together as "No subtype", so
# each bar sums to 100 %. It follows the Data controls' object classes, donor status, autoantibody, age, sex and
# diameter filters; zero and outlier exclusions are per measure and do not apply.

DONOR_COMPOSITION_REGIONS <- c("Core" = "core", "Peri-islet band" = "peri", "Core + peri-islet" = "both")
DONOR_COMPOSITION_ORDERS <- c("Donor status, then beta-cell proportion" = "beta",
                              "Donor status, then donor ID" = "id")
DONOR_COMPOSITION_OTHER <- "No subtype"

# Leave out (owner, 2026-10-21: "the large proportions of endocrine obscure the other types for core, and exocrine does
# the same thing in peri"): the phenotypes left out by default in each region; the bars then show the proportions of
# the remaining cells. On v9 the four endocrine types are 99 % of core cells and Exocrine 93 % of band cells.
DONOR_COMPOSITION_ENDOCRINE <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
donor_composition_default_leave_out <- function(region) {
  switch(region, core = DONOR_COMPOSITION_ENDOCRINE, peri = "Exocrine",
         both = c(DONOR_COMPOSITION_ENDOCRINE, "Exocrine"), character(0))
}

#' Proportions among the cells shown: `tab` without the `leave_out` phenotypes, `proportion` recomputed over the
#' cells that remain (`proportion_all` keeps the proportion of all cells).
donor_composition_shown <- function(tab, leave_out = character(0)) {
  tab$proportion_all <- tab$proportion
  tab <- tab[!(tab$phenotype %in% leave_out), , drop = FALSE]
  if (!nrow(tab)) return(tab)
  shown <- stats::ave(tab$cells, tab$case_id, FUN = sum)
  tab$shown_cells <- shown
  tab$proportion <- ifelse(shown > 0, 100 * tab$cells / shown, NA_real_)
  tab
}

#' Long table, one row per donor and phenotype: cells, the donor's total cells in the region, proportion (%), objects.
#' `comp` holds one row per object (prop_* and peri_prop_* fractions, total_cells_core / total_cells_peri);
#' `pairs` is prepared()$combined_pairs (label, core_col, peri_col).
donor_composition_table <- function(comp, pairs, region = c("core", "peri", "both")) {
  region <- match.arg(region)
  empty <- data.frame(case_id = character(0), donor_status = character(0), phenotype = character(0),
                      cells = numeric(0), total_cells = numeric(0), proportion = numeric(0), n_objects = integer(0))
  if (is.null(comp) || !nrow(comp) || is.null(pairs) || !nrow(pairs)) return(empty)
  # A group aggregate ("All macrophages") sums other labels; stacked beside them it would count their cells twice.
  aggregate <- get0("is_phenotype_aggregate", mode = "function", ifnotfound = function(x) rep(FALSE, length(x)))
  pairs <- pairs[!aggregate(pairs$label), , drop = FALSE]
  num <- function(col) { v <- if (col %in% names(comp)) suppressWarnings(as.numeric(comp[[col]])) else NA_real_
                         v[!is.finite(v)] <- 0; rep_len(v, nrow(comp)) }
  core_total <- num("total_cells_core")
  peri_total <- num("total_cells_peri")
  use_core <- region %in% c("core", "both")
  use_peri <- region %in% c("peri", "both")
  total <- (if (use_core) core_total else 0) + (if (use_peri) peri_total else 0)
  cells <- sapply(seq_len(nrow(pairs)), function(k) {
    (if (use_core) num(pairs$core_col[k]) * core_total else 0) +
      (if (use_peri) num(pairs$peri_col[k]) * peri_total else 0)
  })
  cells <- matrix(cells, nrow = nrow(comp), dimnames = list(NULL, pairs$label))
  ids <- as.character(comp$`Case ID`)
  status <- as.character(comp$`Donor Status`)
  kept <- phenotype_is_kept(pairs$label)
  rows <- lapply(split(seq_len(nrow(comp)), ids), function(i) {
    tot <- sum(total[i])
    n <- colSums(cells[i, , drop = FALSE])
    listed <- n[kept]
    out <- data.frame(phenotype = c(names(listed), DONOR_COMPOSITION_OTHER),
                      cells = c(unname(listed), max(tot - sum(listed), 0)), stringsAsFactors = FALSE)
    out$case_id <- ids[i[1]]
    out$donor_status <- status[i[1]]
    out$total_cells <- tot
    out$proportion <- if (tot > 0) 100 * out$cells / tot else NA_real_
    out$n_objects <- length(i)
    out
  })
  out <- do.call(rbind, rows)
  out <- out[, names(empty)]
  rownames(out) <- NULL
  out
}

#' Donor order for the bars: by status, then beta-cell proportion (high to low) or donor ID.
donor_composition_order <- function(tab, order = c("beta", "id")) {
  order <- match.arg(order)
  donors <- unique(tab[, c("case_id", "donor_status")])
  beta <- tab$proportion[tab$phenotype == "Beta cell"][match(donors$case_id, tab$case_id[tab$phenotype == "Beta cell"])]
  if (is.null(beta) || !length(beta)) beta <- rep(0, nrow(donors))
  beta[!is.finite(beta)] <- -1
  key2 <- if (order == "beta") -beta else suppressWarnings(as.numeric(donors$case_id))
  donors$case_id[order(match(donors$donor_status, c("ND", "Aab+", "T1D")), key2)]
}

donor_composition_figure <- function(tab, order_ids, palette, status_cols, region_label, leave_out = character(0)) {
  if (!nrow(tab)) return(plotly_islet_empty("No donors for the current selection"))
  tab <- donor_composition_shown(tab, leave_out)
  if (!nrow(tab)) return(plotly_islet_empty("Every phenotype is left out"))
  status <- tab$donor_status[match(order_ids, tab$case_id)]
  labels <- paste0(order_ids, "<br>", status)
  tab$x <- factor(labels[match(tab$case_id, order_ids)], levels = labels)
  # A phenotype with no cells in any donor shown (core Exocrine, band CHGA+ hormone-negative on v9) gets no legend entry.
  cells_by <- tapply(tab$cells, tab$phenotype, sum)
  phenos <- unique(tab$phenotype)
  phenos <- phenos[cells_by[phenos] > 0]
  phenos <- c(setdiff(phenos, DONOR_COMPOSITION_OTHER), intersect(DONOR_COMPOSITION_OTHER, phenos))
  # Legend entries under their phenotype group (owner, 2026-10-27): a click on a group's title shows or hides it.
  group_of <- get0("phenotype_group", mode = "function", ifnotfound = function(x) rep("Phenotype", length(x)))
  order_groups <- get0("PHENOTYPE_GROUP_ORDER", ifnotfound = "Phenotype")
  phenos <- phenos[order(match(group_of(phenos), order_groups))]
  titled <- character(0)
  p <- plotly::plot_ly()
  for (ph in phenos) {
    t <- tab[tab$phenotype == ph, , drop = FALSE]
    g <- group_of(ph)
    col <- if (ph == DONOR_COMPOSITION_OTHER) "#8f8f8f" else unname(palette[ph]) %||% "#8f8f8f"
    if (is.na(col)) col <- "#8f8f8f"
    p <- plotly::add_trace(p, type = "bar", x = t$x, y = t$proportion, name = phenotype_text(ph),
                           legendgroup = g, legendgrouptitle = if (!(g %in% titled)) list(text = g),
                           marker = list(color = col), hoverinfo = "text",
                           text = paste0("Donor ", t$case_id, " (", t$donor_status, ")<br>", phenotype_text(ph), ": ",
                                         sprintf("%.1f", t$proportion), "% of the cells shown (",
                                         format(round(t$cells), big.mark = ","), " of ",
                                         format(round(t$shown_cells), big.mark = ","), ")<br>",
                                         sprintf("%.2f", t$proportion_all), "% of all ",
                                         format(round(t$total_cells), big.mark = ","), " cells<br>Objects: ",
                                         t$n_objects),
                           textposition = "none")
    titled <- union(titled, g)
  }
  # A strip of donor-status squares above the bars.
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = factor(labels, levels = labels),
                         y = rep(104, length(labels)), showlegend = FALSE, hoverinfo = "text",
                         text = paste0("Donor ", order_ids, ": ", status),
                         marker = list(symbol = "square", size = 14, color = unname(status_cols[status])))
  fonts <- plotly_islet_fonts()
  p %>%
    plotly::layout(
      barmode = "stack", font = fonts$global, bargap = 0.15,
      title = list(text = plotly_title_with_sub(paste0("Cell make-up of each donor's objects (", region_label, ")"),
                                                if (length(leave_out)) paste0("Leaving out ",
                                                  paste(phenotype_text(leave_out), collapse = ", "), "; bars are proportions of the",
                                                  " remaining cells") else NULL),
                   font = fonts$title, x = 0.02, xanchor = "left"),
      xaxis = list(title = list(text = "Donor", font = fonts$axis_title), tickfont = list(size = ISLET_FONT_SIZES$dense_min),
                   automargin = TRUE),
      yaxis = list(title = list(text = if (length(leave_out)) "Proportion of the cells shown (%)" else
                     "Proportion of cells (%)", font = fonts$axis_title), tickfont = fonts$axis_tick,
                   range = c(0, 108), tickvals = seq(0, 100, 25)),
      legend = list(font = fonts$legend_text, title = list(text = "Phenotype", font = fonts$legend_title),
                    groupclick = "toggleitem", tracegroupgap = 6),
      margin = list(t = 80, b = 80)) %>%
    plotly_islet_config("composition_by_donor")
}

donor_composition_card_ui <- function(ns) {
  div(class = "card", style = "margin-bottom: 20px; padding: 15px; overflow: visible;",
    div(style = "display: flex; align-items: center; justify-content: space-between; gap: 12px; flex-wrap: wrap;",
      card_title_help("Composition by donor", style = "margin: 0; color: var(--islet-fg);", help = list(
        tags$p("Each donor's cell make-up: the proportion of each phenotype among all the cells of the donor's",
               " selected endocrine objects, pooled over the objects (a large object counts for more cells), in the",
               " chosen region. The squares above the bars give each donor's status."),
        tags$p(tags$strong("No subtype"), " (gray) holds the cells typed only to a general class (Endocrine, T cell",
               " or immune cells without a subtype), so each bar sums to 100%."),
        tags$p(tags$strong("Leave out."), " One or two classes hold nearly every cell (the endocrine types in the",
               " core, Exocrine in the band) and flatten the rest, so they are left out by default: the bars then",
               " give each phenotype's proportion of the cells that remain, and hover also gives its proportion of",
               " all cells. Changing the region restores its default; uncheck everything to see all cells."),
        tags$p(tags$strong("Filters."), " The Data controls' endocrine objects, donor status, autoantibody, age, sex",
               " and diameter filters apply. The Data layer, the measure, zero and outlier exclusions do not: they",
               " concern one measure, and this card shows every phenotype."),
        tags$p("The CSV gives each donor's cells, total cells and proportion of all cells per phenotype, every",
               " phenotype included."))),
      div(style = "display: flex; align-items: flex-end; gap: 14px;",
        div(style = "width: 190px;", selectInput(ns("composition_region"), "Region", choices = DONOR_COMPOSITION_REGIONS,
                                                 selected = "core", width = "100%")),
        div(style = "width: 290px;", selectInput(ns("composition_order"), "Order donors by",
                                                 choices = DONOR_COMPOSITION_ORDERS, selected = "beta", width = "100%")),
        div(style = "padding-bottom: 6px;",
          downloadButton(ns("dl_composition_csv"), "CSV", class = "btn-sm btn-outline-secondary",
                         style = "font-size: 11px; padding: 2px 8px;")))
    ),
    # Grouped, with All | None per group and for the list (owner, 2026-10-27); drawn by the module server.
    div(class = "composition-leave-out", style = "margin: 4px 0 6px;", uiOutput(ns("composition_leave_out_ui"))),
    plotlyOutput(ns("composition_by_donor"), height = "520px")
  )
}
