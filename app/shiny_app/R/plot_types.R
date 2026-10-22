# ---------- Plot tab: plot types for either plot card (owner, 2026-10-21) ----------
# "Add the plot types" -- each of the Plot tab's two cards has a Plot type selector offering every type below; the
# left card opens on Against diameter and the right on By donor status, the two plots the cards always showed. Those
# two are drawn in mod_plot_server.R (diameter_figure(), status_figure()); the six added here take the card's rows
# (one per endocrine object: Case ID, islet_key, donor_status, value, islet_diam_um, is_outlier, the donor columns)
# and return plotly figures. The helpers that turn rows into donor or curve tables are pure, so they are tested
# (tests/testthat/test-plot-types.R).

PLOT_TYPE_CHOICES <- c(
  "Against diameter"              = "diameter",
  "By donor status"               = "status",
  "Donor medians"                 = "donor_medians",
  "Cumulative distribution"       = "ecdf",
  "One curve per donor"           = "donor_curves",
  "Against another measure"       = "xy",
  "Against a donor clinical value" = "clinical",
  "Density against diameter"      = "density"
)

PLOT_CLINICAL_CHOICES <- c(
  "Age (years)"                      = "age",
  "BMI"                              = "BMI",
  "HbA1c (%)"                        = "HbA1c_percent",
  "C-peptide (ng/mL)"                = "C_peptide_ng_mL",
  "T1D duration (years)"             = "t1d_duration_years",
  "Autoantibodies (number positive)" = "n_aab",
  # From the nPOD export (owner, 2026-10-22; v10 donor_metadata.csv).
  "Height (cm)"                      = "height_cm",
  "Weight (kg)"                      = "weight_kg",
  "BMI percentile"                   = "BMI_percentile",
  "Admission glucose (mg/dL)"        = "admission_glucose_mg_dL",
  "Peak glucose (mg/dL)"             = "peak_glucose_mg_dL",
  "ICU stay (days)"                  = "ICU_days",
  "Organ transit time (min)"         = "transit_time_min",
  "Pancreas weight (g)"              = "pancreas_weight_g",
  "RNA integrity (RIN)"              = "RIN"
)

PLOT_STATUS_LEVELS <- c("ND", "Aab+", "T1D")

#' The x measures of "Against another measure": object measures, then each retained phenotype's core and peri-islet
#' proportion (%). Values are prepared()$comp columns; `pairs` is prepared()$combined_pairs.
plot_x_measure_choices <- function(comp_cols, pairs = NULL, comp = NULL) {
  object <- c("Diameter (\u00b5m)" = "islet_diam_um", "Core cells" = "total_cells_core",
              "Peri-islet cells" = "total_cells_peri", "Core immune proportion (%)" = "immune_frac_core",
              "Band immune proportion (%)" = "immune_frac_peri")
  object <- object[object %in% comp_cols]
  core <- peri <- character(0)
  if (!is.null(pairs) && nrow(pairs)) {
    pairs <- pairs[phenotype_is_kept(pairs$label), , drop = FALSE]
    # As the measure list, without the near-empty proportions (plot_composition_selection_choices()).
    enough <- function(cols) if (is.null(comp)) rep(TRUE, length(cols)) else
      vapply(cols, function(k) plot_nonzero_fraction(comp, k) >= PLOT_MIN_NONZERO_FRACTION, logical(1))
    has_core <- pairs$core_col %in% comp_cols & enough(pairs$core_col)
    has_peri <- pairs$peri_col %in% comp_cols & enough(pairs$peri_col)
    core <- setNames(pairs$core_col[has_core], paste0(pairs$label[has_core], " (core %)"))
    peri <- setNames(pairs$peri_col[has_peri], paste0(pairs$label[has_peri], " (peri-islet %)"))
  }
  out <- list("Object" = object, "Core proportion" = core, "Peri-islet proportion" = peri)
  out[lengths(out) > 0]
}

#' The x values of `col` for the rows of `rows`, looked up in `comp` by donor and object. Proportions are given in %.
plot_x_values <- function(rows, comp, col) {
  if (identical(col, "islet_diam_um") && "islet_diam_um" %in% names(rows)) return(as.numeric(rows$islet_diam_um))
  key <- paste(comp$`Case ID`, comp$islet_key)
  v <- suppressWarnings(as.numeric(comp[[col]]))[match(paste(rows$`Case ID`, rows$islet_key), key)]
  if (grepl("^(prop_|peri_prop_|immune_frac_)", col)) v <- 100 * v
  v
}

#' A donor's value of a clinical column (constant within a donor). HbA1c 0 is "not recorded" (donors 6090 and 6232).
plot_clinical_values <- function(rows, col) {
  v <- suppressWarnings(as.numeric(rows[[col]]))
  if (identical(col, "HbA1c_percent")) v[!is.na(v) & v <= 0] <- NA
  v
}

#' One row per donor: status, objects, and the median and quartiles of `value` over its objects.
donor_median_table <- function(rows, extra = NULL) {
  if (is.null(rows) || !nrow(rows)) return(data.frame())
  rows <- rows[is.finite(rows$value), , drop = FALSE]
  if (!nrow(rows)) return(data.frame())
  ids <- as.character(rows$`Case ID`)
  out <- do.call(rbind, lapply(split(seq_len(nrow(rows)), ids), function(i) {
    v <- rows$value[i]
    q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
    d <- data.frame(case_id = ids[i[1]], donor_status = as.character(rows$donor_status[i[1]]),
                    n_objects = length(i), median = q[2], q1 = q[1], q3 = q[3], stringsAsFactors = FALSE)
    if (!is.null(extra)) d[[extra]] <- extra_value(rows[[extra]][i])
    d
  }))
  out <- out[order(match(out$donor_status, PLOT_STATUS_LEVELS), suppressWarnings(as.numeric(out$case_id))), ,
             drop = FALSE]
  rownames(out) <- NULL
  out
}
extra_value <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }

#' Cumulative distribution per donor status: each status's sorted values with the proportion of its objects at or
#' below each.
ecdf_table <- function(rows) {
  rows <- rows[is.finite(rows$value), , drop = FALSE]
  out <- do.call(rbind, lapply(intersect(PLOT_STATUS_LEVELS, unique(as.character(rows$donor_status))), function(s) {
    v <- sort(rows$value[rows$donor_status == s])
    data.frame(donor_status = s, value = v, cum_prop = seq_along(v) / length(v), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) data.frame(donor_status = character(0), value = numeric(0), cum_prop = numeric(0)) else out
}

#' One density curve per donor, over the range of all donors' values (`n` points each). A donor whose values are all
#' equal (or fewer than three objects) has no curve; it is drawn as its points.
donor_density_curves <- function(rows, n = 200) {
  rows <- rows[is.finite(rows$value), , drop = FALSE]
  if (!nrow(rows)) return(list())
  lim <- range(rows$value)
  if (diff(lim) == 0) lim <- lim + c(-0.5, 0.5)
  med <- donor_median_table(rows)
  lapply(seq_len(nrow(med)), function(k) {
    v <- rows$value[as.character(rows$`Case ID`) == med$case_id[k]]
    curve <- NULL
    if (length(v) >= 3 && length(unique(v)) > 1) {
      # Half the default bandwidth: proportions pile up at 0 and 100, which the default smooths into a flat band.
      d <- stats::density(v, from = lim[1], to = lim[2], n = n, adjust = 0.5)
      curve <- data.frame(x = d$x, y = d$y / max(d$y))
    }
    list(case_id = med$case_id[k], donor_status = med$donor_status[k], n_objects = med$n_objects[k],
         median = med$median[k], values = v, curve = curve)
  })
}

#' Spearman correlation of two vectors, NA-safe; NULL when fewer than four complete pairs.
plot_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NULL)
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

plot_format_p <- function(p) if (!is.finite(p)) "n/a" else if (p < 1e-16) "< 1e-16" else
  if (p < 0.001) format(p, digits = 2, scientific = TRUE) else
  formatC(p, digits = 3, format = "fg")

# ---- Drawing helpers ----

# Axis range (owner, 2026-10-21: "these outliers compact the rest of the data"). With View controls > Axis range at
# "Fit 99% of objects" (the default), an axis whose few largest values stretch it ends a little above the 99th
# percentile; the objects beyond it are drawn as open triangles at its edge, with their true value in the hover.
PLOT_AXIS_QUANTILE <- 0.99

#' The upper end of an axis fitted to the bulk of `values`, or NULL when every value fits without squeezing the rest
#' (the largest value is within 25 % of the 99th percentile). `keep` are values the axis must still show (summary
#' means). When 99 % of the values are 0, the axis fits the median of the nonzero ones.
plot_axis_upper <- function(values, keep = NULL) {
  v <- values[is.finite(values)]
  if (length(v) < 20) return(NULL)
  q <- stats::quantile(v, PLOT_AXIS_QUANTILE, names = FALSE)
  if (q <= 0) {
    nz <- v[v > 0]
    if (!length(nz)) return(NULL)
    q <- stats::median(nz)
  }
  if (max(v) <= q * 1.25) return(NULL)
  lim <- q * 1.1
  keep <- keep[is.finite(keep)]
  if (length(keep)) lim <- max(lim, max(keep) * 1.05)
  if (max(v) <= lim) NULL else lim
}

#' An axis range from the lowest value to `upper`, padded so markers at either end are not cut.
plot_axis_range <- function(values, upper) {
  low <- min(c(0, values[is.finite(values)]))
  pad <- 0.03 * (upper - low)
  c(low - pad, upper + pad)
}

#' Subtitle note for the objects drawn at an axis edge.
plot_beyond_note <- function(n) if (n > 0) sprintf("%d off scale (\u25b3 at the edge)", n) else NULL

#' Pseudo-log axis (asinh, as the Against diameter plot uses): the transform and tick positions/labels for `values`.
plot_pseudo_log <- function(x) asinh(x / 2) / log(10)
plot_pseudo_log_axis <- function(values) {
  r <- range(values[is.finite(values)], 0)
  brks <- c(0, 10^(0:ceiling(log10(max(abs(r), 1)))))
  brks <- brks[brks <= max(r) * 1.05 | brks == 0]
  list(tickvals = plot_pseudo_log(brks), ticktext = format(brks, big.mark = ",", scientific = FALSE, trim = TRUE))
}

#' Colours of the object points: donor status, or each donor (the Plot tab's Donor ID palette, `donor_pal`).
plot_point_colours <- function(rows, color_by, status_cols, donor_pal) {
  if (identical(color_by, "donor_id")) unname(donor_pal[as.character(rows$`Case ID`)])
  else unname(status_cols[as.character(rows$donor_status)])
}

plot_donor_label <- function(case_id, status) paste0(case_id, " (", status, ")")

#' The common layout of the added plot types: fonts, title (with an optional second line), axes, legend.
plot_type_layout <- function(p, title, subtitle = NULL, xaxis = list(), yaxis = list(), legend_title = NULL,
                             filename_root, margin = list(t = if (is.null(subtitle)) 60 else 85, b = 70)) {
  fonts <- plotly_islet_fonts()
  p %>%
    plotly::layout(
      font = fonts$global,
      title = list(text = plotly_title_with_sub(title, subtitle), font = fonts$title, x = 0.02, xanchor = "left"),
      xaxis = utils::modifyList(list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick,
                                     automargin = TRUE, zeroline = FALSE), xaxis),
      yaxis = utils::modifyList(list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick,
                                     automargin = TRUE, zeroline = FALSE), yaxis),
      legend = list(font = fonts$legend_text, title = list(text = legend_title %||% "", font = fonts$legend_title),
                    orientation = "h", x = 0, y = -0.18, xanchor = "left", yanchor = "top"),
      margin = margin, hovermode = "closest") %>%
    plotly_islet_config(filename_root)
}

#' Status legend entries (empty traces), so a plot whose points carry per-point colours still has a legend.
plot_status_legend <- function(p, status_cols, levels) {
  for (s in levels) {
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = list(NA), y = list(NA), name = s,
                           marker = list(color = unname(status_cols[s]), size = 11), hoverinfo = "skip",
                           showlegend = TRUE, inherit = FALSE)
  }
  p
}

# ---- The six added figures. `ylab` names the measure; `name` is its short form for titles. ----

#' Donor medians: one point per donor (median of its objects, bars the quartiles), a line per status at the median of
#' its donors, and a Kruskal-Wallis test on the donor medians, so every donor counts once.
plot_donor_medians_figure <- function(rows, ylab, name, status_cols, color_by = "donor_status", donor_pal = NULL,
                                      log_scale = FALSE, point_size = 3) {
  med <- donor_median_table(rows)
  if (!nrow(med)) return(plotly_islet_empty("No data for the current selection"))
  levels <- intersect(PLOT_STATUS_LEVELS, med$donor_status)
  tf <- if (log_scale) plot_pseudo_log else identity
  # Spread the donors of a status evenly across its slot, in ID order, so the points never overlap.
  med$x <- NA_real_
  for (s in levels) {
    i <- which(med$donor_status == s)
    med$x[i] <- match(s, levels) + if (length(i) > 1) seq(-0.28, 0.28, length.out = length(i)) else 0
  }
  med$colour <- if (identical(color_by, "donor_id") && !is.null(donor_pal)) unname(donor_pal[med$case_id]) else
    unname(status_cols[med$donor_status])
  med$hover <- paste0("Donor ", plot_donor_label(med$case_id, med$donor_status), "<br>Objects: ", med$n_objects,
                      "<br>Median: ", signif(med$median, 3), "<br>Quartiles: ", signif(med$q1, 3), " \u2013 ",
                      signif(med$q3, 3))
  p <- plotly::plot_ly()
  for (s in levels) {
    m <- stats::median(med$median[med$donor_status == s])
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = match(s, levels) + c(-0.38, 0.38),
                           y = tf(c(m, m)), line = list(color = unname(status_cols[s]), width = 3),
                           hoverinfo = "text", text = paste0(s, " median of donor medians: ", signif(m, 3)),
                           showlegend = FALSE)
  }
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = med$x, y = tf(med$median), text = med$hover,
                         hoverinfo = "text", showlegend = FALSE,
                         error_y = list(type = "data", symmetric = FALSE, array = tf(med$q3) - tf(med$median),
                                        arrayminus = tf(med$median) - tf(med$q1), color = "#9d9d9d",
                                        thickness = 1.2, width = 0),
                         marker = list(color = med$colour, size = 4 + 3 * point_size, line = list(color = "#000000",
                                       width = 0.6)))
  kw <- if (length(levels) > 1) tryCatch(stats::kruskal.test(med$median, factor(med$donor_status))$p.value,
                                         error = function(e) NA_real_) else NA_real_
  sub <- paste0(nrow(med), " donors; bars are each donor's quartiles",
                if (is.finite(kw)) paste0("; Kruskal-Wallis on donor medians p = ", plot_format_p(kw)) else "")
  yaxis <- list(title = list(text = paste0(ylab, ", donor median")))
  if (log_scale) yaxis <- c(yaxis, plot_pseudo_log_axis(c(med$q1, med$q3)))
  plot_type_layout(p, paste0(name, ": donor medians"), sub,
                   xaxis = list(title = list(text = "Donor status"), tickvals = seq_along(levels), ticktext = levels,
                                range = c(0.4, length(levels) + 0.6), showgrid = FALSE),
                   yaxis = yaxis, filename_root = "plot_donor_medians")
}

#' Cumulative distribution: for each donor status, the proportion of its objects at or below each value.
plot_ecdf_figure <- function(rows, ylab, name, status_cols, log_scale = FALSE, clip = FALSE) {
  tab <- ecdf_table(rows)
  if (!nrow(tab)) return(plotly_islet_empty("No data for the current selection"))
  tf <- if (log_scale) plot_pseudo_log else identity
  p <- plotly::plot_ly()
  for (s in unique(tab$donor_status)) {
    t <- tab[tab$donor_status == s, , drop = FALSE]
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = tf(t$value), y = 100 * t$cum_prop,
                           line = list(color = unname(status_cols[s]), width = 2.5, shape = "hv"),
                           name = paste0(s, " (", nrow(t), " objects)"), hoverinfo = "text",
                           text = paste0(s, ": ", sprintf("%.1f", 100 * t$cum_prop), "% of objects at or below ",
                                         signif(t$value, 3)))
  }
  xaxis <- list(title = list(text = ylab))
  if (log_scale) xaxis <- c(xaxis, plot_pseudo_log_axis(tab$value))
  upper <- if (clip && !log_scale) plot_axis_upper(tab$value) else NULL
  if (!is.null(upper)) xaxis$range <- plot_axis_range(tab$value, upper)
  zero_txt <- paste(vapply(unique(tab$donor_status), function(s) sprintf("%s %.0f%%", s,
    100 * mean(tab$value[tab$donor_status == s] == 0)), ""), collapse = ", ")
  plot_type_layout(p, paste0(name, ": cumulative distribution"),
                   paste0(c(paste0("Objects at 0: ", zero_txt), if (!is.null(upper))
                     sprintf("axis ends at %s; %d objects lie beyond it", signif(upper, 3), sum(tab$value > upper))),
                     collapse = "; "),
                   xaxis = xaxis, yaxis = list(title = list(text = "Objects at or below (%)"), range = c(0, 102)),
                   filename_root = "plot_cumulative")
}

#' One curve per donor: each donor's density of the measure over its objects, scaled to the same height, stacked in
#' donor-status order (a ridgeline). A tick marks each donor's median.
plot_donor_curves_figure <- function(rows, ylab, name, status_cols, log_scale = FALSE, clip = FALSE) {
  tf <- if (log_scale) plot_pseudo_log else identity
  upper <- if (clip && !log_scale) plot_axis_upper(rows$value) else NULL
  rows$value <- tf(rows$value)
  curves <- donor_density_curves(rows)
  if (!length(curves)) return(plotly_islet_empty("No data for the current selection"))
  n <- length(curves)
  p <- plotly::plot_ly()
  labels <- character(n)
  for (k in seq_len(n)) {
    cv <- curves[[k]]
    off <- n - k  # first donor (ND, lowest ID) at the top
    col <- unname(status_cols[cv$donor_status])
    labels[k] <- plot_donor_label(cv$case_id, cv$donor_status)
    hover <- paste0("Donor ", labels[k], "<br>Objects: ", cv$n_objects, "<br>Median: ",
                    signif(if (log_scale) 2 * sinh(cv$median * log(10)) else cv$median, 3))
    if (!is.null(cv$curve)) {
      p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = c(cv$curve$x, rev(cv$curve$x)),
                             y = c(off + 0.9 * cv$curve$y, rep(off, nrow(cv$curve))), fill = "toself",
                             fillcolor = grDevices::adjustcolor(col, alpha.f = 0.45), line = list(color = col, width = 1.2),
                             hoverinfo = "text", text = hover, hoveron = "fills+points", showlegend = FALSE)
    } else {
      p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = cv$values, y = rep(off + 0.2, length(cv$values)),
                             marker = list(color = col, size = 7), hoverinfo = "text", text = hover, showlegend = FALSE)
    }
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = rep(cv$median, 2), y = off + c(0, 0.9),
                           line = list(color = "#000000", width = 1.5, dash = "dot"), hoverinfo = "text",
                           text = hover, showlegend = FALSE)
  }
  p <- plot_status_legend(p, status_cols, intersect(PLOT_STATUS_LEVELS, vapply(curves, `[[`, "", "donor_status")))
  xaxis <- list(title = list(text = ylab))
  if (log_scale) xaxis <- c(xaxis, plot_pseudo_log_axis(2 * sinh(rows$value * log(10))))
  if (!is.null(upper)) xaxis$range <- plot_axis_range(rows$value, upper)
  plot_type_layout(p, paste0(name, ": one curve per donor"),
                   paste0(c("Each curve is scaled to its own peak; the dotted tick is the donor's median",
                            if (!is.null(upper)) sprintf("axis ends at %s; %d objects lie beyond it",
                                                         signif(upper, 3), sum(rows$value > upper, na.rm = TRUE))),
                          collapse = "; "),
                   xaxis = xaxis,
                   yaxis = list(title = list(text = ""), tickvals = (n - seq_len(n)) + 0.35, ticktext = labels,
                                showgrid = FALSE, range = c(-0.2, n + 0.1),
                                tickfont = list(size = ISLET_FONT_SIZES$dense_min)),
                   filename_root = "plot_donor_curves")
}

#' Against another measure: every object's measure (y) against a second per-object measure (x), with Spearman's rho
#' over objects and over donor medians. Object points carry link keys.
plot_xy_figure <- function(rows, x, xlab, ylab, name, status_cols, color_by = "donor_status", donor_pal = NULL,
                           log_scale = FALSE, point_size = 3, point_alpha = 0.4, show_points = TRUE,
                           clip = FALSE) {
  rows$x <- x
  rows <- rows[is.finite(rows$x) & is.finite(rows$value), , drop = FALSE]
  if (!nrow(rows)) return(plotly_islet_empty("No objects have both values"))
  tf <- if (log_scale) plot_pseudo_log else identity
  rows$hover <- paste0("Donor ", plot_donor_label(rows$`Case ID`, rows$donor_status), "<br>", rows$islet_key,
                       "<br>", xlab, ": ", signif(rows$x, 3), "<br>", ylab, ": ", signif(rows$value, 3))
  x_upper <- if (clip) plot_axis_upper(rows$x) else NULL
  y_upper <- if (clip && !log_scale) plot_axis_upper(rows$value) else NULL
  px <- if (is.null(x_upper)) rows$x else pmin(rows$x, x_upper * 1.02)
  py <- if (is.null(y_upper)) rows$value else pmin(rows$value, y_upper * 1.02)
  beyond <- px != rows$x | py != rows$value
  p <- plotly::plot_ly()
  if (isTRUE(show_points)) {
    colours <- plot_point_colours(rows, color_by, status_cols, donor_pal)
    for (out in c(FALSE, TRUE)) {
      i <- which(beyond == out)
      if (!length(i)) next
      p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = px[i], y = tf(py[i]), key = rows$link_key[i],
                             text = rows$hover[i], hoverinfo = "text", showlegend = FALSE,
                             marker = if (out) list(symbol = "triangle-up-open", color = colours[i],
                                                    size = 6 + 2 * point_size, line = list(width = 2))
                                      else list(color = colours[i], size = 3 + 2 * point_size,
                                                opacity = 1 - point_alpha, line = list(color = "#000000", width = 0.4)))
    }
  }
  med <- donor_median_table(transform(rows, value = x))
  med_y <- donor_median_table(rows)
  med$y <- med_y$median[match(med$case_id, med_y$case_id)]
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = med$median, y = tf(med$y), showlegend = FALSE,
                         hoverinfo = "text", text = paste0("Donor ", plot_donor_label(med$case_id, med$donor_status),
                                                           " median<br>", xlab, ": ", signif(med$median, 3), "<br>",
                                                           ylab, ": ", signif(med$y, 3)),
                         marker = list(symbol = "diamond", size = 13, color = unname(status_cols[med$donor_status]),
                                       line = list(color = "#000000", width = 1.2)))
  p <- plot_status_legend(p, status_cols, intersect(PLOT_STATUS_LEVELS, rows$donor_status))
  r_obj <- plot_spearman(rows$x, rows$value)
  r_don <- plot_spearman(med$median, med$y)
  sub <- paste0(c(if (!is.null(r_obj)) sprintf("Spearman \u03c1 %.2f over %s objects", r_obj$rho,
                                               format(r_obj$n, big.mark = ",")),
                  if (!is.null(r_don)) sprintf("%.2f over %d donor medians (p = %s)", r_don$rho, r_don$n,
                                               plot_format_p(r_don$p))), collapse = "; ")
  yaxis <- list(title = list(text = ylab))
  if (log_scale) yaxis <- c(yaxis, plot_pseudo_log_axis(rows$value))
  if (!is.null(y_upper)) yaxis$range <- plot_axis_range(rows$value, y_upper * 1.02)
  xaxis <- list(title = list(text = xlab))
  if (!is.null(x_upper)) xaxis$range <- plot_axis_range(rows$x, x_upper * 1.02)
  sub <- paste0(c(if (nzchar(sub)) sub, "diamonds are donor medians", plot_beyond_note(sum(beyond))), collapse = "; ")
  plot_type_layout(p, paste0(name, " against ", xlab), sub,
                   xaxis = xaxis, yaxis = yaxis, filename_root = "plot_against_measure")
}

#' Against a donor clinical value: each donor's median (bars its quartiles) at its clinical value, its objects behind
#' it, and Spearman's rho over the donors that have a value.
plot_clinical_figure <- function(rows, clinical, clinical_label, ylab, name, status_cols, color_by = "donor_status",
                                 donor_pal = NULL, log_scale = FALSE, point_size = 3, point_alpha = 0.4,
                                 show_points = TRUE, clip = FALSE, missing_label = "without a recorded value") {
  rows$clinical <- clinical
  med <- donor_median_table(rows, extra = "clinical")
  if (!nrow(med)) return(plotly_islet_empty("No data for the current selection"))
  missing <- sum(!is.finite(med$clinical))
  med <- med[is.finite(med$clinical), , drop = FALSE]
  rows <- rows[is.finite(rows$clinical), , drop = FALSE]
  if (!nrow(med)) return(plotly_islet_empty(paste("No selected donor has a recorded", clinical_label)))
  tf <- if (log_scale) plot_pseudo_log else identity
  span <- diff(range(med$clinical))
  jitter <- if (span > 0) span * 0.012 else 0.05
  y_upper <- if (clip && !log_scale) plot_axis_upper(rows$value, keep = med$median) else NULL
  py <- if (is.null(y_upper)) rows$value else pmin(rows$value, y_upper * 1.02)
  beyond <- py != rows$value
  p <- plotly::plot_ly()
  if (isTRUE(show_points) && nrow(rows)) {
    set.seed(17)
    jx <- rows$clinical + stats::runif(nrow(rows), -jitter, jitter)
    colours <- plot_point_colours(rows, color_by, status_cols, donor_pal)
    hover <- paste0("Donor ", plot_donor_label(rows$`Case ID`, rows$donor_status), "<br>",
                    rows$islet_key, "<br>", ylab, ": ", signif(rows$value, 3))
    for (out in c(FALSE, TRUE)) {
      i <- which(beyond == out)
      if (!length(i)) next
      p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = jx[i], y = tf(py[i]),
                             key = rows$link_key[i], hoverinfo = "text", text = hover[i], showlegend = FALSE,
                             marker = if (out) list(symbol = "triangle-up-open", color = colours[i],
                                                    size = 6 + 1.5 * point_size, line = list(width = 2))
                                      else list(color = colours[i], size = 2 + 1.5 * point_size,
                                                opacity = 0.6 * (1 - point_alpha)))
    }
  }
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = med$clinical, y = tf(med$median),
                         showlegend = FALSE, hoverinfo = "text",
                         text = paste0("Donor ", plot_donor_label(med$case_id, med$donor_status), "<br>",
                                       clinical_label, ": ", signif(med$clinical, 3), "<br>Objects: ", med$n_objects,
                                       "<br>Median: ", signif(med$median, 3)),
                         error_y = list(type = "data", symmetric = FALSE, array = tf(med$q3) - tf(med$median),
                                        arrayminus = tf(med$median) - tf(med$q1), color = "#9d9d9d",
                                        thickness = 1.5, width = 4),
                         marker = list(symbol = "diamond", size = 14, color = unname(status_cols[med$donor_status]),
                                       line = list(color = "#000000", width = 1.2)))
  p <- plot_status_legend(p, status_cols, intersect(PLOT_STATUS_LEVELS, med$donor_status))
  r <- plot_spearman(med$clinical, med$median)
  sub <- paste0(c(if (!is.null(r)) sprintf("Spearman \u03c1 %.2f over %d donors (p = %s)", r$rho, r$n,
                                           plot_format_p(r$p)) else sprintf("%d donors", nrow(med)),
                  if (missing > 0) sprintf("%d donor%s %s left out", missing,
                                           if (missing == 1) "" else "s", missing_label)), collapse = "; ")
  yaxis <- list(title = list(text = ylab))
  if (log_scale) yaxis <- c(yaxis, plot_pseudo_log_axis(c(med$q1, med$q3, rows$value)))
  if (!is.null(y_upper)) yaxis$range <- plot_axis_range(rows$value, y_upper * 1.02)
  sub <- paste0(c(sub, "diamonds are donor medians", plot_beyond_note(sum(beyond & isTRUE(show_points)))),
                collapse = "; ")
  plot_type_layout(p, paste0(name, " against ", clinical_label), sub,
                   xaxis = list(title = list(text = clinical_label)), yaxis = yaxis,
                   filename_root = "plot_against_clinical")
}

#' Density against diameter: contour lines of each donor status's object density in (diameter, measure), each status
#' normalised to its own objects, where the scatter's points pile up.
plot_density_figure <- function(rows, ylab, name, status_cols, log_scale = FALSE, log_x = FALSE, clip = FALSE) {
  rows <- rows[is.finite(rows$value) & is.finite(rows$islet_diam_um), , drop = FALSE]
  if (nrow(rows) < 10) return(plotly_islet_empty("Too few objects for a density"))
  tf <- if (log_scale) plot_pseudo_log else identity
  p <- plotly::plot_ly()
  for (s in intersect(PLOT_STATUS_LEVELS, rows$donor_status)) {
    r <- rows[rows$donor_status == s, , drop = FALSE]
    col <- unname(status_cols[s])
    p <- plotly::add_trace(p, type = "histogram2dcontour", x = r$islet_diam_um, y = tf(r$value),
                           histnorm = "probability", ncontours = 8, name = paste0(s, " (", nrow(r), " objects)"),
                           contours = list(coloring = "lines", showlabels = FALSE),
                           colorscale = list(c(0, col), c(1, col)), showscale = FALSE,
                           line = list(width = 1.6), showlegend = TRUE, hoverinfo = "name")
  }
  yaxis <- list(title = list(text = ylab))
  if (log_scale) yaxis <- c(yaxis, plot_pseudo_log_axis(rows$value))
  upper <- if (clip && !log_scale) plot_axis_upper(rows$value) else NULL
  if (!is.null(upper)) yaxis$range <- plot_axis_range(rows$value, upper)
  plot_type_layout(p, paste0(name, ": density against diameter"),
                   "Contours of each donor status's objects, each status normalised to its own count",
                   xaxis = list(title = list(text = "Islet diameter (\u00b5m)"), type = if (log_x) "log" else "linear"),
                   yaxis = yaxis, filename_root = "plot_density")
}

#' Help paragraphs for the plot cards: what each plot type shows.
plot_types_help <- function() {
  shiny::tagList(
    shiny::tags$p("Choose what this card plots with its Plot type selector; both cards offer every type and plot the",
                  " measure chosen in the Data controls."),
    shiny::tags$p(shiny::tags$strong("Against diameter."),
                  " Every object against its diameter, with each donor status's summary statistic per diameter bin",
                  " (View controls) and the bin width in the Data controls. The diameter is that of a circle with",
                  " the object's area."),
    shiny::tags$p(shiny::tags$strong("By donor status."), " A violin or box per donor status (View controls)."),
    shiny::tags$p(shiny::tags$strong("Donor medians."),
                  " One point per donor, the median of its objects, with bars at its quartiles; the line is the",
                  " median of each status's donors. The Kruskal-Wallis p compares the donor medians, so every donor",
                  " counts once however many objects it has."),
    shiny::tags$p(shiny::tags$strong("Cumulative distribution."),
                  " For each donor status, the proportion of its objects at or below each value, which shows the",
                  " objects at zero that a violin smooths over."),
    shiny::tags$p(shiny::tags$strong("One curve per donor."),
                  " Each donor's distribution over its objects, scaled to the same height and stacked in donor-status",
                  " order, with its median marked; a donor with fewer than three objects or one value shows points."),
    shiny::tags$p(shiny::tags$strong("Against another measure."),
                  " Every object's measure against a second one chosen on the card (an object measure, or a",
                  " phenotype's core or peri-islet proportion), with donor medians as diamonds and Spearman's",
                  " \u03c1 over objects and over donor medians."),
    shiny::tags$p(shiny::tags$strong("Against a donor clinical value."),
                  " Each donor's median (diamond, bars its quartiles) at a donor value from nPOD (age, BMI, HbA1c,",
                  " C-peptide, T1D duration, autoantibody count, height, weight, BMI percentile, admission and peak",
                  " glucose, ICU stay, organ transit time, pancreas weight, RNA integrity), its objects behind it, and Spearman's \u03c1 over the donors. Donors",
                  " without a value are left out: T1D duration applies only to T1D donors (nPOD records five of the six",
                  " as diagnosed during their terminal admission, duration 0, and one at 7 months), and the two donors",
                  " whose HbA1c nPOD records as 0 were not measured."),
    shiny::tags$p(shiny::tags$strong("Density against diameter."),
                  " Contour lines of each donor status's objects in diameter and measure, each status normalised to",
                  " its own count, for where the scatter's points pile up."),
    shiny::tags$p(shiny::tags$strong("Axis range."),
                  " With Fit 99% of objects (View controls, the default), an axis that a few extreme objects would",
                  " stretch ends just above the 99th percentile, and those objects are drawn as open triangles at its",
                  " edge, with their true value in the hover; the count is in the subtitle. The violins, boxes, summary",
                  " statistics and correlations still use every value. Fit every object shows the whole range."),
    shiny::tags$p(shiny::tags$strong("Exclusions."), " The card's Exclude zero values box (Data controls) and the outlier",
                  " filter apply to every type; outliers are left out of the added types when the filter is on."),
    shiny::tags$p(shiny::tags$strong("Linked points."),
                  " In the types that draw objects (Against diameter, By donor status, Against another measure, Against",
                  " a donor clinical value), hover over an object to ring it in both cards; click it to open its",
                  " segmentation map and single cells below. Clear selection releases it. An object selected on",
                  " another tab is ringed here too."),
    shiny::tags$p("The CSV button exports what the card shows: binned summaries, per-object values or one row per",
                  " donor, as the type plots."))
}
