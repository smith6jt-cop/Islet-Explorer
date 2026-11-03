# ---------- Plot tab: Composition by donor card (owner, 2026-10-21) ----------
# The plot cards show one phenotype at a time; this card shows every donor's whole make-up: the proportion of each
# retained phenotype among all the cells of the donor's selected objects, pooled (cells of the phenotype summed over
# the objects, over all their cells), in the core, the 20 um peri-islet band or both. Cells of the general classes
# without a resolved lineage/subtype count in the total and form an Unclassified subclass, so
# each bar sums to 100 % of the cells shown. Hidden QC-conflict cells are excluded from the displayed partition,
# but remain in the original total_cells denominator. It follows the Data controls' object classes, donor status,
# autoantibody, age, sex and
# diameter filters; zero and outlier exclusions are per measure and do not apply.

DONOR_COMPOSITION_REGIONS <- c("Core" = "core", "Peri-islet band" = "peri", "Core + peri-islet" = "both")
DONOR_COMPOSITION_ORDERS <- c("Donor status, then beta-cell proportion" = "beta",
                              "Donor status, then donor ID" = "id")
DONOR_COMPOSITION_OTHER <- unname(UNCLASSIFIED_KEYS[["subtype"]])

# The Include selector preserves the established region defaults: the abundant endocrine types in the core
# and Exocrine cells in the band start unchecked so the other populations remain visible.
DONOR_COMPOSITION_ENDOCRINE <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
donor_composition_default_include <- function(region, choices) {
  omitted <- switch(region, core = DONOR_COMPOSITION_ENDOCRINE, peri = "Exocrine",
         both = c(DONOR_COMPOSITION_ENDOCRINE, "Exocrine"), character(0))
  setdiff(choices, omitted)
}

#' Proportions among the selected phenotypes. NULL includes every category; an empty vector includes none.
#' The source count/denominator stays unchanged and `proportion_all` retains the original proportion.
donor_composition_shown <- function(tab, include = NULL) {
  tab$proportion_all <- tab$proportion
  if (!is.null(include)) tab <- tab[tab$phenotype %in% include, , drop = FALSE]
  tab$shown_cells <- numeric(nrow(tab))
  if (!nrow(tab)) return(tab)
  shown <- stats::ave(tab$cells, tab$case_id, FUN = sum)
  tab$shown_cells <- shown
  tab$proportion <- ifelse(shown > 0, 100 * tab$cells / shown, NA_real_)
  tab
}

#' Export the plotted rows while retaining the existing original-denominator `proportion` column.
donor_composition_export <- function(shown) {
  shown$proportion_selected <- shown$proportion
  shown$proportion <- shown$proportion_all
  shown$proportion_all <- NULL
  shown
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
  # The saved pair inventory can omit labels whose sanitized peri names cannot
  # be reversed (notably a slash in a Conflict label). Count hidden cells from
  # the original regional columns, independently of that display inventory.
  hidden <- numeric(nrow(comp))
  for (prefix in c(if (use_core) "prop_", if (use_peri) "peri_prop_")) {
    columns <- names(comp)[startsWith(names(comp), prefix)]
    columns <- columns[phenotype_is_conflict(substring(columns, nchar(prefix) + 1L))]
    region_total <- if (prefix == "prop_") core_total else peri_total
    for (column in columns) hidden <- hidden + num(column) * region_total
  }
  cells <- sapply(seq_len(nrow(pairs)), function(k) {
    (if (use_core) num(pairs$core_col[k]) * core_total else 0) +
      (if (use_peri) num(pairs$peri_col[k]) * peri_total else 0)
  })
  cells <- matrix(cells, nrow = nrow(comp), dimnames = list(NULL, pairs$label))
  ids <- as.character(comp$`Case ID`)
  status <- as.character(comp$`Donor Status`)
  display_labels <- unclassified_collapse(pairs$label)
  kept <- display_labels %in% phenotype_partition_labels(pairs$label)
  rows <- lapply(split(seq_len(nrow(comp)), ids), function(i) {
    tot <- sum(total[i])
    n <- colSums(cells[i, , drop = FALSE])
    leaves <- unique(display_labels[kept])
    listed <- stats::setNames(vapply(leaves, function(leaf) sum(n[kept & display_labels == leaf]), numeric(1)), leaves)
    # Residuals on older roots retain their original accounting. On the current
    # root the three broad source labels are explicit members of this leaf.
    residual <- max(tot - sum(listed) - sum(hidden[i]), 0)
    if (DONOR_COMPOSITION_OTHER %in% names(listed)) {
      listed[[DONOR_COMPOSITION_OTHER]] <- listed[[DONOR_COMPOSITION_OTHER]] + residual
    } else listed <- c(listed, stats::setNames(residual, DONOR_COMPOSITION_OTHER))
    out <- data.frame(phenotype = names(listed), cells = unname(listed), stringsAsFactors = FALSE)
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

#' Draw the same selected, normalized rows used by the CSV export.
donor_composition_figure <- function(tab, order_ids, palette, status_cols, region_label) {
  if (!nrow(tab)) return(plotly_islet_empty("No phenotypes selected or no donors match the filters"))
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
    col <- unname(palette[ph]) %||% "#8f8f8f"
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
      title = list(text = paste0("Cell make-up of each donor's objects (", region_label, ")"),
                   font = fonts$title, x = 0.02, xanchor = "left"),
      xaxis = list(title = list(text = "Donor", font = fonts$axis_title), tickfont = list(size = ISLET_FONT_SIZES$dense_min),
                   automargin = TRUE),
      yaxis = list(title = list(text = "Proportion of selected phenotypes", font = fonts$axis_title), tickfont = fonts$axis_tick,
                   range = c(0, 108), tickvals = seq(0, 100, 25)),
      legend = list(font = fonts$legend_text, title = list(text = "Phenotype", font = fonts$legend_title),
                    groupclick = "toggleitem", tracegroupgap = 6),
      margin = list(t = 48, b = 80)) %>%
    plotly_islet_config("composition_by_donor")
}

donor_composition_card_ui <- function(ns) {
  div(class = "card", style = "margin-bottom: 20px; padding: 15px; overflow: visible;",
    div(style = "display: flex; align-items: center; justify-content: space-between; gap: 12px; flex-wrap: wrap;",
      card_title_help("Composition By Donor", style = "margin: 0; color: var(--islet-fg);", help = list(
        tags$p("Each donor's cell make-up: the proportion of each phenotype among all the cells of the donor's",
               " selected endocrine objects, pooled over the objects (a large object counts for more cells), in the",
               " chosen region. The squares above the bars give each donor's status."),
        tags$p("The Unclassified family has four mutually exclusive subclasses describing why a more specific",
               " identity was not assigned. Hidden QC-conflict cells are excluded from these categories.",
               " Each bar sums to 100% of the cells shown; hover also reports the proportion",
               " of all original cells."),
        tags$p(tags$strong("Include."), " Open Include Phenotypes and check the populations to compare. Each bar",
               " gives their proportions among the selected phenotypes. Endocrine types in the core and Exocrine",
               " cells in the band start unchecked because they otherwise obscure the other populations. Changing",
               " the region restores its default. All selects every retained category; None clears the selection.",
               " Phenotype Colors controls the palette and stays synchronized with Spatial and the cell viewer."),
        tags$p(tags$strong("Filters."), " The Data controls' endocrine objects, donor status, autoantibody, age, sex",
               " and diameter filters apply. The Data layer, the measure, zero and outlier exclusions do not: they",
               " concern one measure, and this card shows every phenotype."),
        tags$p("The CSV contains the plotted, included categories. cells and total_cells retain the original",
               " counts; proportion retains the percentage of all original cells, while proportion_selected gives",
               " the plotted percentage of the selected phenotypes and shown_cells its denominator. Original",
               " proportions need not sum to 100%."),
        unclassified_help_ui())),
      downloadButton(ns("dl_composition_csv"), "CSV", class = "btn-sm btn-outline-secondary")
    ),
    div(class = "composition-toolbar", style = "display: flex; align-items: center; gap: 8px 16px; flex-wrap: wrap; margin: 8px 0;",
      div(class = "plot-type-select", style = "flex: 0 0 200px; min-width: 0;",
        selectInput(ns("composition_region"), "Region", choices = DONOR_COMPOSITION_REGIONS,
                    selected = "core", width = "100%")),
      div(class = "plot-type-select", style = "flex: 2 1 360px; min-width: 0;",
        selectInput(ns("composition_order"), "Order Donors By", choices = DONOR_COMPOSITION_ORDERS,
                    selected = "beta", width = "100%")),
      div(class = "plot-type-select", style = "flex: 1.2 1 280px; min-width: 0;",
        selectInput("composition_palette", "Phenotype Colors", choices = names(PHENOTYPE_PALETTES),
                    selected = "High Contrast", width = "100%"))
    ),
    # The same grouped selection as Spatial, kept collapsed to reserve room for the plot.
    div(class = "composition-include", style = "margin: 0 0 4px;",
        control_panel("Include Phenotypes", uiOutput(ns("composition_include_ui")))),
    plotlyOutput(ns("composition_by_donor"), height = "520px")
  )
}
