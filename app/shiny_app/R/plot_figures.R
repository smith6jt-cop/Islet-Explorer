# ---------- Direct plotly figures: the Plot cards and the Trajectory scatter (Plan A Level 2, B1) ----------
# ggplotly() conversion was 36 % of the server time of a Plot change and about 90 % of each Plot step (performance
# baseline, 2026-10-24): it walks the ggplot build with one R-level loop per layer and point and composes a hover
# text from every mapped aesthetic. The three figures are markers, summary points with error bars, LOESS lines and
# violin or box layers, so since 2026-10-26 they are built here with plot_ly() on the data frames the modules
# prepare. Each builder is a pure function of its rows and an options list (tests/testthat/test-plot-figures.R);
# the module passes the inputs and attaches the linked-plot behaviour (trajectory_link_plot()). The statistics are
# those of the ggplot layers they replace: the same jitter widths and seeds (ggplot2::position_jitter), the same
# LOESS (stats::loess with the layer's span and weights) and the same 95 % band (the t quantile on the loess
# degrees of freedom, as ggplot2's predictdf.loess).

PLOTLY_PX_PER_MM <- 96 / 25.4

#' A ggplot2 point size (mm) as a plotly marker size (px): the conversion ggplotly() applied.
plotly_point_px <- function(size_mm) size_mm * PLOTLY_PX_PER_MM

#' Evaluate `expr` with the RNG seeded at `seed`, then put the session's RNG stream back as it was (a set.seed()
#' inside a render function would otherwise reset the process-wide stream).
islet_with_seed <- function(seed, expr) {
  saved <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(if (is.null(saved)) suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)) else
    assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(expr)
}

#' ggplot2::position_jitter(width, height, seed) on one layer's points: x then y, each uniform within +/- the
#' amount, none for an amount of 0. A layer's rows go in together, as ggplot jitters them in one draw.
plot_jitter <- function(x, y, width, height, seed) {
  islet_with_seed(seed, {
    if (isTRUE(width > 0)) x <- x + stats::runif(length(x), -width, width)
    if (isTRUE(height > 0)) y <- y + stats::runif(length(y), -height, height)
    list(x = x, y = y)
  })
}

#' position_jitter()'s default height: 40 % of the resolution of the layer's y values (so integer counts are
#' spread by +/- 0.4 and continuous values barely at all).
plot_jitter_default_height <- function(y) {
  y <- y[is.finite(y)]
  if (length(y) < 2) return(0)
  0.4 * ggplot2::resolution(y, zero = FALSE, discrete = TRUE)
}

#' The pseudo-log (asinh) y-axis breaks of the Plot cards: 0 and the powers of `base` up to the data, thinned where
#' the labels would overlap in transformed space.
plot_pseudo_log_breaks <- function(limits, base = 10) {
  max_val <- suppressWarnings(max(limits[is.finite(limits)]))
  if (!is.finite(max_val) || max_val <= 0) return(0)
  max_pow <- ceiling(log(max_val, base))
  min_pow <- if (max_val < 1) floor(log(max_val, base)) else 0
  candidates <- sort(unique(c(0, base^seq(min_pow, max_pow))))
  candidates <- candidates[candidates <= max_val * 1.1]
  trans_pos <- asinh(candidates / 2) / log(base)
  total_range <- diff(range(trans_pos))
  if (total_range == 0) return(candidates)
  min_gap <- max(0.3, total_range * 0.15)
  keep <- logical(length(candidates))
  keep[1] <- TRUE
  last_kept <- trans_pos[1]
  for (i in seq_along(candidates)[-1]) {
    if (trans_pos[i] - last_kept >= min_gap) {
      keep[i] <- TRUE
      last_kept <- trans_pos[i]
    }
  }
  candidates[keep]
}

#' tickvals / ticktext of a pseudo-log y-axis over `values` (data scale).
plot_pseudo_log_yaxis <- function(values) {
  brks <- plot_pseudo_log_breaks(range(values[is.finite(values)], 0))
  list(tickvals = plot_pseudo_log(brks), ticktext = scales::label_comma()(brks))
}

#' The object keys of rows: `link_key` when the module set it, else from Case ID and islet_key.
plot_row_keys <- function(rows) {
  if ("link_key" %in% names(rows)) return(as.character(rows$link_key))
  if (all(c("Case ID", "islet_key") %in% names(rows))) return(trajectory_link_key(rows$`Case ID`, rows$islet_key))
  rep(NA_character_, nrow(rows))
}

#' The compact hover text of each key (plot_point_hover_map()); "" where a key has none.
plot_hover_text <- function(hover, key) {
  text <- if (length(hover)) unname(hover[key]) else rep(NA_character_, length(key))
  text[is.na(text)] <- ""
  text
}

#' One marker trace of linked objects: keys and compact hover, five significant figures on the coordinates.
plot_add_points <- function(p, x, y, key, text, colour, size, opacity = 1, symbol = NULL, line = NULL) {
  marker <- list(color = colour, size = size, opacity = opacity)
  if (!is.null(symbol)) marker$symbol <- symbol
  if (!is.null(line)) marker$line <- line
  plotly::add_trace(p, type = "scatter", mode = "markers", x = signif(x, 5), y = signif(y, 5), key = key,
                    text = text, hoverinfo = "text", showlegend = FALSE, marker = marker)
}

#' A LOESS curve as ggplot2::geom_smooth(method = "loess") fits it: stats::loess with `span` (and `weights`),
#' predicted at `n` points over the x range; with `se`, the 95 % band from the t quantile on the fit's degrees of
#' freedom (ggplot2's predictdf.loess). NULL when the fit fails (too few points for the span).
plot_loess_curve <- function(x, y, weights = NULL, span = 0.75, se = FALSE, n = 80) {
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(weights)) ok <- ok & is.finite(weights)
  d <- data.frame(x = x[ok], y = y[ok])
  if (!is.null(weights)) d$w <- weights[ok]
  if (nrow(d) < 3 || length(unique(d$x)) < 2) return(NULL)
  fit <- tryCatch(suppressWarnings(
    if (is.null(weights)) stats::loess(y ~ x, data = d, span = span)
    else stats::loess(y ~ x, data = d, weights = w, span = span)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  grid <- seq(min(d$x), max(d$x), length.out = n)
  pred <- tryCatch(suppressWarnings(stats::predict(fit, data.frame(x = grid), se = isTRUE(se))),
                   error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  out <- if (isTRUE(se)) {
    data.frame(x = grid, y = as.numeric(pred$fit), ci = as.numeric(pred$se.fit) * stats::qt(0.975, pred$df))
  } else data.frame(x = grid, y = as.numeric(pred), ci = NA_real_)
  out[is.finite(out$y), , drop = FALSE]
}

#' The curve of plot_loess_curve() as a line, with its band when it has one. `fill` is the band colour (rgba).
plot_add_loess <- function(p, curve, colour, width, fill = NULL, legendgroup = NULL) {
  if (is.null(curve) || !nrow(curve)) return(p)
  if (!is.null(fill) && any(is.finite(curve$ci))) {
    b <- curve[is.finite(curve$ci), , drop = FALSE]
    p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = signif(c(b$x, rev(b$x)), 5),
                           y = signif(c(b$y + b$ci, rev(b$y - b$ci)), 5), fill = "toself", fillcolor = fill,
                           line = list(width = 0), hoverinfo = "skip", showlegend = FALSE, legendgroup = legendgroup)
  }
  plotly::add_trace(p, type = "scatter", mode = "lines", x = signif(curve$x, 5), y = signif(curve$y, 5),
                    line = list(color = colour, width = width), hoverinfo = "skip", showlegend = FALSE,
                    legendgroup = legendgroup)
}

#' A colour with an alpha, as plotly wants it.
plot_rgba <- function(colour, alpha) {
  rgb <- grDevices::col2rgb(colour)[, 1]
  sprintf("rgba(%d,%d,%d,%.2f)", rgb[1], rgb[2], rgb[3], alpha)
}

# ---- Against diameter ----

#' Against diameter: every object against its diameter (jittered, coloured by donor status or donor), the objects
#' above a fitted axis as open triangles at its top, flagged outliers as grey crosses, and each donor status's
#' binned summary (points, error bars, a line or a LOESS trend) on top.
#' `rows`: the card's flagged, binned rows (value, diam_mid, donor_status, is_outlier, Case ID, islet_key,
#' link_key); `sm`: their binned summary (donor_status, diam_mid, y, ymin, ymax, n); `hover`: text by key.
#' `o`: list(color_by, status_cols, donor_pal, show_points, pt_size, pt_alpha, log_x, log_y, binwidth, counts,
#' smooth, stat_label, ylab, title, clip).
plot_diameter_figure <- function(rows, sm, hover, o) {
  if (is.null(sm) || !nrow(sm)) return(plotly_islet_empty("Loading data..."))
  levels <- PLOT_STATUS_LEVELS
  log_x <- isTRUE(o$log_x)
  log_y <- isTRUE(o$log_y)
  tf <- if (log_y) plot_pseudo_log else identity
  sm <- sm[is.finite(sm$diam_mid), , drop = FALSE]
  sm$donor_status <- as.character(sm$donor_status)
  sm <- sm[order(match(sm$donor_status, levels), sm$diam_mid), , drop = FALSE]
  sm_y <- tf(sm$y)
  sm_lo <- tf(sm$ymin)
  sm_hi <- tf(sm$ymax)
  value <- suppressWarnings(as.numeric(rows$value))
  axis_values <- c(sm$y, sm$ymin, sm$ymax)
  # Axis range (View controls): fit 99 % of objects; the objects above are open triangles at the top edge.
  y_upper <- if (isTRUE(o$clip) && !log_y) plot_axis_upper(value, keep = sm$y) else NULL
  y_top <- if (is.null(y_upper)) NULL else y_upper * 1.02
  donor_id <- identical(o$color_by, "donor_id") && !is.null(o$donor_pal)
  pt_size <- o$pt_size %||% 3
  px <- plotly_point_px(pt_size)
  jx <- function(v) if (log_x) log2(v) else v
  ujx <- function(v) if (log_x) 2^v else v
  n_beyond <- 0L
  p <- plotly::plot_ly()

  if (isTRUE(o$show_points) && nrow(rows)) {
    x <- suppressWarnings(as.numeric(rows$diam_mid))
    y <- tf(value)
    axis_values <- c(axis_values, value)
    shown <- is.finite(x) & is.finite(y) & x > 0
    outlier <- rows$is_outlier %in% TRUE
    beyond <- if (is.null(y_top)) rep(FALSE, length(y)) else y > y_top
    y[beyond] <- y_top
    key <- plot_row_keys(rows)
    text <- plot_hover_text(hover, key)
    colours <- plot_point_colours(rows, if (donor_id) "donor_id" else "donor_status", o$status_cols, o$donor_pal)
    colours[is.na(colours)] <- "#9d9d9d"
    # The jitter of the ggplot layers: width 0.15 in log2 units on a log axis, else 35 % of the bin width; height
    # 2 % of the drawn range (at least 0.5) for counts, 1 % otherwise.
    jw <- if (log_x) 0.15 else max(1, suppressWarnings(as.numeric(o$binwidth %||% 50)) * 0.35)
    if (!is.finite(jw)) jw <- 17.5
    yr <- suppressWarnings(range(y[shown]))
    ydiff <- if (all(is.finite(yr))) diff(yr) else 0
    jh <- if (isTRUE(o$counts)) max(0.02 * ydiff, 0.5) else max(0.01 * ydiff, 0)
    if (!is.finite(jh)) jh <- if (isTRUE(o$counts)) 0.5 else 0
    normal <- shown & !outlier & !beyond
    if (any(normal)) {
      j <- plot_jitter(jx(x[normal]), y[normal], jw, jh, 17)
      p <- plot_add_points(p, ujx(j$x), j$y, key[normal], text[normal], colours[normal], px,
                           opacity = 1 - (o$pt_alpha %||% 0.4), line = list(color = "#000000", width = 0.5))
    }
    edge <- shown & beyond & !outlier
    n_beyond <- sum(edge)
    if (n_beyond > 0) {
      j <- plot_jitter(jx(x[edge]), y[edge], jw, 0, 23)
      p <- plot_add_points(p, ujx(j$x), j$y, key[edge], text[edge], colours[edge], px * 1.3,
                           symbol = "triangle-up-open", line = list(color = colours[edge], width = 1.5))
    }
    ghost <- shown & outlier
    if (any(ghost)) {
      # Excluded outliers stay visible as grey crosses, so a reviewer sees where the exclusions land.
      j <- plot_jitter(jx(x[ghost]), y[ghost], jw, if (is.null(y_top)) jh else 0, 31)
      p <- plot_add_points(p, ujx(j$x), j$y, key[ghost], text[ghost], "#737373", px * 1.1, opacity = 0.55,
                           symbol = "x-thin-open", line = list(color = "#737373", width = 1.3))
    }
  }

  # Summary per donor status on top of the points: a line (or a LOESS trend over the bins), the summary points
  # with their error bars. The summary points carry the status legend unless the points are coloured by donor.
  for (s in intersect(levels, unique(sm$donor_status))) {
    i <- which(sm$donor_status == s)
    col <- unname(o$status_cols[s])
    if (is.na(col)) col <- "#9d9d9d"
    group <- paste0("status-", s)
    if (isTRUE(o$smooth)) {
      # Fitted over the bins in the axis's units (log2 on a log axis), as geom_smooth() fitted it.
      curve <- plot_loess_curve(jx(sm$diam_mid[i]), sm_y[i], span = 0.6)
      if (!is.null(curve)) curve$x <- ujx(curve$x)
      p <- plot_add_loess(p, curve, col, 2.5, legendgroup = group)
    } else {
      p <- plotly::add_trace(p, type = "scatter", mode = "lines", x = signif(sm$diam_mid[i], 5),
                             y = signif(sm_y[i], 5), line = list(color = col, width = 3), hoverinfo = "skip",
                             showlegend = FALSE, legendgroup = group)
    }
    p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = signif(sm$diam_mid[i], 5),
                           y = signif(sm_y[i], 5), name = s, legendgroup = group, showlegend = !donor_id,
                           hoverinfo = "text",
                           text = paste0(s, "<br>Bin centre: ", signif(sm$diam_mid[i], 4), " µm<br>",
                                         o$stat_label %||% "Summary", ": ", signif(sm$y[i], 4), " (",
                                         signif(sm$ymin[i], 4), " to ", signif(sm$ymax[i], 4), ")<br>Objects: ",
                                         sm$n[i]),
                           error_y = list(type = "data", symmetric = FALSE, array = signif(sm_hi[i] - sm_y[i], 5),
                                          arrayminus = signif(sm_y[i] - sm_lo[i], 5), color = col,
                                          thickness = 1.5, width = 0),
                           marker = list(color = col, size = plotly_point_px(2.5),
                                         line = list(color = "#000000", width = 1)))
  }
  if (donor_id) p <- plot_status_legend(p, o$donor_pal, names(o$donor_pal))

  xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
  if (!is.finite(xmax)) xmax <- 300
  xaxis <- list(title = list(text = "Islet diameter (µm)"))
  if (log_x) {
    pow <- 0:ceiling(log2(xmax))
    xaxis <- c(xaxis, list(type = "log", tickvals = 2^pow, ticktext = format(2^pow, scientific = FALSE, trim = TRUE)))
  } else {
    xaxis <- c(xaxis, list(tick0 = 0, dtick = 50))
  }
  yaxis <- list(title = list(text = o$ylab))
  if (log_y) yaxis <- c(yaxis, plot_pseudo_log_yaxis(axis_values))
  if (!is.null(y_top)) yaxis <- c(yaxis, list(range = plot_axis_range(value, y_top), autorange = FALSE))
  plot_type_layout(p, o$title, plot_beyond_note(n_beyond), xaxis = xaxis, yaxis = yaxis,
                   legend_title = if (donor_id) "Donor ID" else NULL, filename_root = "plot_scatter",
                   margin = list(t = if (n_beyond > 0) 85 else 60, b = 90), legend = list(y = -0.2))
}

# ---- By donor status ----

#' By donor status: a violin or box per donor status drawn from every kept value, every object as a jittered point
#' (coloured by status or donor), the objects above a fitted axis as open triangles at its top, flagged outliers as
#' grey crosses on top. `rows`: the card's flagged rows (value, donor_status, is_outlier, Case ID, islet_key,
#' link_key); `hover`: text by key. `o`: list(color_by, status_cols, donor_pal, show_points, pt_size, pt_alpha,
#' log_y, box, ylab, title, clip).
plot_status_figure <- function(rows, hover, o) {
  if (is.null(rows) || !nrow(rows)) return(plotly_islet_empty("No data for the current selection"))
  status <- as.character(rows$donor_status)
  levels <- intersect(PLOT_STATUS_LEVELS, unique(status))
  if (!length(levels)) return(plotly_islet_empty("No data for the current selection"))
  log_y <- isTRUE(o$log_y)
  tf <- if (log_y) plot_pseudo_log else identity
  value <- suppressWarnings(as.numeric(rows$value))
  finite <- is.finite(value) & status %in% levels
  pos <- match(status, levels)
  outlier <- rows$is_outlier %in% TRUE
  y_upper <- if (isTRUE(o$clip) && !log_y) plot_axis_upper(value) else NULL
  y_top <- if (is.null(y_upper)) Inf else y_upper * 1.02
  donor_id <- identical(o$color_by, "donor_id") && !is.null(o$donor_pal)
  pt_size <- o$pt_size %||% 3
  px <- plotly_point_px(pt_size)
  key <- plot_row_keys(rows)
  text <- plot_hover_text(hover, key)
  colours <- plot_point_colours(rows, if (donor_id) "donor_id" else "donor_status", o$status_cols, o$donor_pal)
  colours[is.na(colours)] <- "#9d9d9d"
  p <- plotly::plot_ly()

  # The distribution of each status over its kept values (every value, including those above the axis).
  for (k in seq_along(levels)) {
    i <- which(finite & !outlier & pos == k)
    if (!length(i)) next
    col <- unname(o$status_cols[levels[k]])
    if (is.na(col)) col <- "#9d9d9d"
    fill <- plot_rgba(col, 0.8)
    if (isTRUE(o$box)) {
      p <- plotly::add_trace(p, type = "box", x0 = k, y = signif(tf(value[i]), 5), name = levels[k],
                             fillcolor = fill, line = list(color = "#333333", width = 1), boxpoints = FALSE,
                             width = 0.75, hoverinfo = if (log_y) "skip" else "y", showlegend = FALSE)
    } else {
      p <- plotly::add_trace(p, type = "violin", x0 = k, y = signif(tf(value[i]), 5), name = levels[k],
                             fillcolor = fill, line = list(color = "#333333", width = 1), points = FALSE,
                             spanmode = "soft", box = list(visible = FALSE), meanline = list(visible = FALSE),
                             width = 0.8, hoverinfo = "skip", showlegend = FALSE)
    }
  }

  kept <- finite & !outlier
  beyond <- kept & value > y_top
  n_beyond <- sum(beyond)
  if (isTRUE(o$show_points)) {
    pts <- kept & !beyond
    if (any(pts)) {
      y <- tf(value[pts])
      j <- plot_jitter(pos[pts], y, 0.18, plot_jitter_default_height(y), 17)
      p <- plot_add_points(p, j$x, j$y, key[pts], text[pts], colours[pts], px,
                           opacity = 1 - (o$pt_alpha %||% 0.4), line = list(color = "#000000", width = 0.5))
    }
    if (n_beyond > 0) {
      j <- plot_jitter(pos[beyond], rep(y_top, n_beyond), 0.18, 0, 23)
      p <- plot_add_points(p, j$x, j$y, key[beyond], text[beyond], colours[beyond], px * 1.3,
                           symbol = "triangle-up-open", line = list(color = colours[beyond], width = 1.5))
    }
  }
  ghost <- finite & outlier
  if (any(ghost)) {
    y <- tf(pmin(value[ghost], y_top))
    j <- plot_jitter(pos[ghost], y, 0.18, plot_jitter_default_height(y), 31)
    p <- plot_add_points(p, j$x, j$y, key[ghost], text[ghost], "#737373", px * 1.4, opacity = 0.6,
                         symbol = "x-thin-open", line = list(color = "#737373", width = 1.3))
  }
  if (donor_id) p <- plot_status_legend(p, o$donor_pal, names(o$donor_pal))

  yaxis <- list(title = list(text = o$ylab))
  if (log_y) yaxis <- c(yaxis, plot_pseudo_log_yaxis(value))
  if (is.finite(y_top)) yaxis <- c(yaxis, list(range = plot_axis_range(value, y_top), autorange = FALSE))
  plot_type_layout(p, o$title, plot_beyond_note(n_beyond),
                   xaxis = list(title = list(text = "Donor Status"), tickvals = seq_along(levels), ticktext = levels,
                                range = c(0.4, length(levels) + 0.6), showgrid = FALSE),
                   yaxis = yaxis, legend_title = if (donor_id) "Donor ID" else NULL,
                   filename_root = "plot_distribution",
                   margin = list(t = if (n_beyond > 0) 85 else 60, b = 70),
                   legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left", yanchor = "top"),
                   showlegend = donor_id)
}

# ---- Trajectory scatter ----

#' A ggplot2 area size scale (scale_size_continuous(range = range_mm)) as plotly marker sizes in px.
plot_area_size_px <- function(v, range_mm) {
  v <- suppressWarnings(as.numeric(v))
  lim <- suppressWarnings(range(v[is.finite(v)]))
  r <- if (all(is.finite(lim)) && diff(lim) > 0) (v - lim[1]) / diff(lim) else rep(0.5, length(v))
  r[!is.finite(r)] <- 0
  plotly_point_px(range_mm[1] + sqrt(r) * diff(range_mm))
}

#' The Trajectory scatter: every fitted object's value against its pseudotime (jittered), coloured by donor status,
#' donor or Leiden cluster, sized uniformly or by diameter or cell count, with the LOESS trend (overall, or one per
#' donor status) and its 95 % band, weighted by log1p(cells) as the ggplot layers were. `df`: traj_data_clean()
#' rows with link_key; `o`: list(color_by, colours (named over the colour column's values), colour_label, size_by,
#' point_size, alpha, trend, show_ci, status_cols, ink, x_label, y_label, pt_ext (traj_pt_extent(df$pt))).
trajectory_scatter_figure <- function(df, o) {
  n <- nrow(df)
  colour_col <- if (!is.null(o$color_by) && o$color_by %in% names(df)) o$color_by else "donor_status"
  group <- as.character(df[[colour_col]])
  col <- unname(o$colours[group])
  col[is.na(col)] <- "#9d9d9d"
  weight <- log1p(as.numeric(df$total_cells))
  pt <- as.numeric(df$pt)
  value <- as.numeric(df$value)
  vr <- suppressWarnings(range(value[is.finite(value)]))
  if (!all(is.finite(vr))) vr <- c(0, 1)
  jy <- diff(vr) * 0.08
  j <- islet_with_seed(42, {
    x <- pt + stats::runif(n, -o$pt_ext$jitter, o$pt_ext$jitter)
    y <- value + stats::runif(n, -jy, jy)
    list(x = x, y = y)
  })
  hover <- paste0(o$colour_label %||% "Group", ": ", group,
                  "<br>Donor: ", df$case_id, "<br>Islet: ", df$islet_key,
                  "<br>Pseudotime: ", signif(pt, 3), "<br>Value: ", signif(value, 3),
                  "<br>Cells: ", df$total_cells)
  point_size <- o$point_size %||% 3
  size <- if (identical(o$size_by, "islet_diam_um") && "islet_diam_um" %in% names(df) &&
              any(is.finite(suppressWarnings(as.numeric(df$islet_diam_um))))) {
    plot_area_size_px(df$islet_diam_um, c(0.5, 2.5) * point_size)
  } else if (identical(o$size_by, "total_cells") && "total_cells" %in% names(df)) {
    plot_area_size_px(sqrt(as.numeric(df$total_cells)), c(0.3, 3.0) * point_size)
  } else plotly_point_px(point_size)
  # One trace in a fixed interleaved order, so no donor status is painted over another (islet_umap_plotly()).
  ord <- islet_umap_order(n)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p, type = "scatter", mode = "markers", x = signif(j$x[ord], 5), y = signif(j$y[ord], 5),
                         key = plot_row_keys(df)[ord], text = hover[ord], hoverinfo = "text", showlegend = FALSE,
                         marker = list(color = col[ord], size = if (length(size) > 1) size[ord] else size,
                                       opacity = o$alpha %||% 0.6, line = list(width = 0)))
  trend <- o$trend %||% "by_donor"
  if (identical(trend, "overall")) {
    p <- plot_add_loess(p, plot_loess_curve(pt, value, weight, span = 0.75, se = isTRUE(o$show_ci)),
                        o$ink %||% "#000000", 2.5, fill = "rgba(153,153,153,0.4)")
  } else if (identical(trend, "by_donor")) {
    status <- as.character(df$donor_status)
    for (s in intersect(PLOT_STATUS_LEVELS, unique(status))) {
      i <- which(status == s)
      colour <- unname(o$status_cols[s])
      if (is.na(colour)) colour <- "#9d9d9d"
      p <- plot_add_loess(p, plot_loess_curve(pt[i], value[i], weight[i], span = 0.75, se = isTRUE(o$show_ci)),
                          colour, 2, fill = plot_rgba(colour, 0.15))
    }
  }
  fonts <- plotly_islet_fonts()
  p <- plotly::layout(p,
    showlegend = FALSE, title = NULL, font = fonts$global, hovermode = "closest",
    margin = list(l = 80, r = 20, t = 30, b = 70, pad = 5),
    xaxis = list(fixedrange = FALSE, automargin = TRUE, range = o$pt_ext$lim, autorange = FALSE,
                 title = list(text = o$x_label, font = fonts$axis_title), tickfont = fonts$axis_tick),
    yaxis = list(fixedrange = FALSE, automargin = TRUE, range = vr + c(-0.05, 0.05) * diff(vr), autorange = FALSE,
                 title = list(text = o$y_label, font = fonts$axis_title), tickfont = fonts$axis_tick))
  plotly_islet_config(p, "trajectory_scatter")
}
