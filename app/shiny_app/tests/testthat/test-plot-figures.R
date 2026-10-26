# The direct plotly figures of the Plot cards and the Trajectory scatter (R/plot_figures.R, Plan A Level 2,
# 2026-10-26): traces, object keys, pinned and ghost points, axes, legends, the jitter of the ggplot layers they
# replaced and the LOESS trends, on synthetic rows.

.plot_figures_env <- function() {
  env <- new.env(parent = globalenv())
  env$PHENOTYPE_RETIRED <- c("Endocrine", "T cell", "Immune")
  env$phenotype_is_kept <- function(x) !(as.character(x) %in% env$PHENOTYPE_RETIRED)
  env$ISLET_FONT_SIZES <- list(base = 18, title = 18, subtitle = 15, axis_title = 16, axis_text = 16,
                               legend_title = 16, legend_text = 16, strip = 16, dense_min = 14, hover = 15,
                               annotation = 15)
  env$plotly_islet_fonts <- function() list(global = list(size = 18), title = list(size = 18),
                                            axis_title = list(size = 16), axis_tick = list(size = 16),
                                            legend_title = list(size = 16), legend_text = list(size = 16),
                                            colorbar_title = list(size = 16), hover = list(size = 15),
                                            annotation = list(size = 15))
  env$plotly_islet_config <- function(p, ...) p
  env$plotly_title_with_sub <- function(title, subtitle = NULL)
    if (is.null(subtitle) || !nzchar(trimws(subtitle))) title else paste0(title, "<br><sup>", subtitle, "</sup>")
  env$plotly_islet_empty <- function(msg) plotly::layout(plotly::plotly_empty(type = "scatter", mode = "markers"),
                                                        title = list(text = msg))
  env$islet_umap_order <- function(n) order((seq_len(n) * 0.6180339887498949) %% 1)
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  sys.source(.find_source("R/trajectory_links.R"), envir = env)
  sys.source(.find_source("R/plot_types.R"), envir = env)
  sys.source(.find_source("R/plot_figures.R"), envir = env)
  env
}

.status_cols <- c("ND" = "#3366CC", "Aab+" = "#FF6600", "T1D" = "#109618")

# 120 objects of three donors; row 1 an extreme flagged outlier, row 2 an unflagged object above the fitted axis.
.figure_rows <- function(n_per = 40) {
  set.seed(1)
  status <- rep(c("ND", "Aab+", "T1D"), each = n_per)
  diam <- stats::runif(3 * n_per, 20, 300)
  value <- abs(stats::rnorm(3 * n_per, 10, 4))
  value[1] <- 500
  value[2] <- 150
  rows <- data.frame(`Case ID` = rep(c(6001, 6002, 6003), each = n_per), islet_key = paste0("Islet_", seq_len(3 * n_per)),
                     donor_status = status, value = value, islet_diam_um = diam,
                     diam_mid = 12.5 + 25 * floor(diam / 25), is_outlier = c(TRUE, rep(FALSE, 3 * n_per - 1)),
                     check.names = FALSE, stringsAsFactors = FALSE)
  rows$link_key <- paste0(rows$`Case ID`, "|", rows$islet_key)
  rows
}
.figure_summary <- function(rows) {
  rows <- rows[!rows$is_outlier, ]
  out <- do.call(rbind, lapply(split(rows, list(rows$donor_status, rows$diam_mid), drop = TRUE), function(d) {
    m <- mean(d$value); se <- stats::sd(d$value) / sqrt(nrow(d)); if (!is.finite(se)) se <- 0
    data.frame(donor_status = d$donor_status[1], diam_mid = d$diam_mid[1], y = m, ymin = m - se, ymax = m + se,
               n = nrow(d), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
.hover <- function(rows) stats::setNames(paste0("hover ", rows$islet_key), rows$link_key)
.diameter_opts <- function(...) utils::modifyList(
  list(color_by = "donor_status", status_cols = .status_cols, donor_pal = NULL, show_points = TRUE, pt_size = 3,
       pt_alpha = 0.4, log_x = FALSE, log_y = FALSE, binwidth = 25, counts = FALSE, smooth = FALSE,
       stat_label = "Mean ± SE", ylab = "Value", title = "Value vs Islet Size", clip = TRUE), list(...))
.status_opts <- function(...) utils::modifyList(
  list(color_by = "donor_status", status_cols = .status_cols, donor_pal = NULL, show_points = TRUE, pt_size = 3,
       pt_alpha = 0.4, log_y = FALSE, box = FALSE, ylab = "Value", title = "Distribution", clip = TRUE), list(...))
.traces <- function(fig) plotly::plotly_build(fig)$x$data
.keyed <- function(traces) Filter(function(t) !is.null(t$key) && identical(t$mode, "markers"), traces)
.keys <- function(trace) as.character(unlist(trace$key))
# plotly_build() gives every trace an error_y color; the summary points are those with an error array.
.with_error <- function(traces) Filter(function(t) !is.null(t$error_y$array), traces)

test_that("the diameter figure draws the objects with their keys, pins the object above the axis and ghosts the outlier", {
  skip_if_not_installed("plotly")
  env <- .plot_figures_env()
  rows <- .figure_rows()
  fig <- env$plot_diameter_figure(rows, .figure_summary(rows), .hover(rows), .diameter_opts())
  b <- plotly::plotly_build(fig)
  keyed <- .keyed(b$x$data)
  expect_length(keyed, 3)
  normal <- keyed[[1]]
  expect_equal(sort(.keys(normal)), sort(rows$link_key[-(1:2)]))
  expect_equal(unlist(normal$text)[1], paste0("hover ", rows$islet_key[3]))
  triangle <- keyed[[2]]
  expect_equal(.keys(triangle), rows$link_key[2])
  expect_equal(triangle$marker$symbol, "triangle-up-open")
  ghost <- keyed[[3]]
  expect_equal(.keys(ghost), rows$link_key[1])
  expect_equal(ghost$marker$symbol, "x-thin-open")
  # Both sit at the top edge of the fitted axis, below the true values.
  expect_lt(max(unlist(triangle$y)), 150)
  expect_equal(unlist(triangle$y), unlist(ghost$y))
  expect_false(isTRUE(b$x$layout$yaxis$autorange))
  expect_match(b$x$layout$title$text, "1 off scale")
  # One summary line and one summary point trace (with error bars and the legend) per donor status.
  summary <- .with_error(b$x$data)
  expect_equal(vapply(summary, `[[`, "", "name"), c("ND", "Aab+", "T1D"))
  expect_true(all(vapply(summary, function(t) isTRUE(t$showlegend), logical(1))))
  lines <- Filter(function(t) identical(t$mode, "lines"), b$x$data)
  expect_length(lines, 3)
})

test_that("the diameter figure follows its View options", {
  skip_if_not_installed("plotly")
  env <- .plot_figures_env()
  rows <- .figure_rows()
  sm <- .figure_summary(rows)
  hover <- .hover(rows)
  # Log axes: a log2 x-axis with power-of-two ticks and a pseudo-log y-axis with its own ticks; the log y-axis is
  # never clipped, so nothing is pinned.
  b <- plotly::plotly_build(env$plot_diameter_figure(rows, sm, hover, .diameter_opts(log_x = TRUE, log_y = TRUE)))
  expect_equal(b$x$layout$xaxis$type, "log")
  expect_equal(b$x$layout$xaxis$tickvals[1:3], c(1, 2, 4))
  expect_true(length(b$x$layout$yaxis$tickvals) >= 2)
  expect_length(.keyed(b$x$data), 2)
  expect_false(grepl("off scale", b$x$layout$title$text))
  # LOESS: a smooth curve per status replaces the straight summary line.
  b <- plotly::plotly_build(env$plot_diameter_figure(rows, sm, hover, .diameter_opts(smooth = TRUE)))
  lines <- Filter(function(t) identical(t$mode, "lines"), b$x$data)
  expect_length(lines, 3)
  expect_true(all(vapply(lines, function(t) length(t$x) == 80, logical(1))))
  # Color by donor: the points carry the donor colors, the legend lists the donors, not the statuses.
  pal <- c("6001" = "#111111", "6002" = "#222222", "6003" = "#333333")
  b <- plotly::plotly_build(env$plot_diameter_figure(rows, sm, hover, .diameter_opts(color_by = "donor_id", donor_pal = pal)))
  normal <- .keyed(b$x$data)[[1]]
  expect_true(all(unlist(normal$marker$color) %in% pal))
  legend <- Filter(function(t) isTRUE(t$showlegend), b$x$data)
  expect_equal(vapply(legend, `[[`, "", "name"), names(pal))
  expect_equal(b$x$layout$legend$title$text, "Donor ID")
  # No points: only the summaries.
  b <- plotly::plotly_build(env$plot_diameter_figure(rows, sm, hover, .diameter_opts(show_points = FALSE)))
  expect_length(.keyed(b$x$data), 0)
  expect_false(grepl("off scale", b$x$layout$title$text))
})

test_that("the status figure draws a violin or box per status present, the points with keys, and the ghosts on top", {
  skip_if_not_installed("plotly")
  env <- .plot_figures_env()
  rows <- .figure_rows()
  b <- plotly::plotly_build(env$plot_status_figure(rows, .hover(rows), .status_opts()))
  violins <- Filter(function(t) identical(t$type, "violin"), b$x$data)
  expect_equal(vapply(violins, `[[`, "", "name"), c("ND", "Aab+", "T1D"))
  expect_equal(vapply(violins, `[[`, 1, "x0"), c(1, 2, 3))
  # Every kept value is in its violin, the flagged outlier in none.
  expect_equal(sum(vapply(violins, function(t) length(t$y), 1L)), nrow(rows) - 1L)
  keyed <- .keyed(b$x$data)
  expect_length(keyed, 3)
  expect_equal(length(.keys(keyed[[1]])), nrow(rows) - 2L)
  expect_equal(.keys(keyed[[2]]), rows$link_key[2])
  expect_equal(.keys(keyed[[3]]), rows$link_key[1])
  expect_equal(b$x$layout$xaxis$ticktext, c("ND", "Aab+", "T1D"))
  expect_false(isTRUE(b$x$layout$showlegend))
  expect_match(b$x$layout$title$text, "1 off scale")
  # Boxes instead of violins; a status without rows has no slot.
  b <- plotly::plotly_build(env$plot_status_figure(rows[rows$donor_status != "Aab+", ], .hover(rows), .status_opts(box = TRUE)))
  boxes <- Filter(function(t) identical(t$type, "box"), b$x$data)
  expect_equal(vapply(boxes, `[[`, "", "name"), c("ND", "T1D"))
  expect_equal(b$x$layout$xaxis$ticktext, c("ND", "T1D"))
  # Log y: ticks on the pseudo-log scale, every point shown.
  b <- plotly::plotly_build(env$plot_status_figure(rows, .hover(rows), .status_opts(log_y = TRUE)))
  expect_true(length(b$x$layout$yaxis$tickvals) >= 2)
  expect_length(.keyed(b$x$data), 2)
})

test_that("the jitter matches ggplot2::position_jitter() with the same seed, and leaves the RNG stream alone", {
  skip_if_not_installed("ggplot2")
  env <- .plot_figures_env()
  d <- data.frame(x = rep(1:3, each = 5), y = rep(c(1, 2, 3, 5, 8), 3))
  layer <- ggplot2::layer_data(ggplot2::ggplot(d, ggplot2::aes(x, y)) +
                                 ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.18, height = 0.4, seed = 17)))
  j <- env$plot_jitter(d$x, d$y, 0.18, 0.4, 17)
  expect_equal(j$x, layer$x)
  expect_equal(j$y, layer$y)
  expect_equal(env$plot_jitter_default_height(d$y), 0.4 * ggplot2::resolution(d$y, zero = FALSE))
  set.seed(99); before <- stats::runif(3)
  set.seed(99); invisible(env$plot_jitter(d$x, d$y, 1, 1, 17)); after <- stats::runif(3)
  expect_equal(before, after)
})

test_that("the LOESS curve is the weighted loess fit with ggplot2's 95 % band", {
  env <- .plot_figures_env()
  set.seed(3)
  x <- stats::runif(200); y <- sin(3 * x) + stats::rnorm(200, 0, 0.2); w <- log1p(stats::rpois(200, 20))
  curve <- env$plot_loess_curve(x, y, w, span = 0.75, se = TRUE)
  fit <- stats::loess(y ~ x, data = data.frame(x = x, y = y, w = w), weights = w, span = 0.75)
  pred <- stats::predict(fit, data.frame(x = curve$x), se = TRUE)
  expect_equal(curve$y, as.numeric(pred$fit))
  expect_equal(curve$ci, as.numeric(pred$se.fit) * stats::qt(0.975, pred$df))
  expect_equal(nrow(curve), 80)
  expect_null(env$plot_loess_curve(1:2, c(1, 2), span = 0.75))
})

.traj_df <- function(n = 300) {
  set.seed(5)
  status <- sample(c("ND", "Aab+", "T1D"), n, replace = TRUE)
  data.frame(case_id = as.character(6000 + match(status, c("ND", "Aab+", "T1D"))), islet_key = paste0("Islet_", seq_len(n)),
             pt = stats::runif(n), value = stats::rnorm(n), donor_status = status, donor_id = status,
             total_cells = stats::rpois(n, 30) + 10L, islet_diam_um = stats::runif(n, 20, 200),
             link_key = paste0(6000 + match(status, c("ND", "Aab+", "T1D")), "|Islet_", seq_len(n)),
             stringsAsFactors = FALSE)
}
.traj_opts <- function(df, ...) utils::modifyList(
  list(color_by = "donor_status", colours = .status_cols, colour_label = "Donor Status", size_by = "uniform",
       point_size = 3, alpha = 0.6, trend = "by_donor", show_ci = TRUE, status_cols = .status_cols, ink = "#000000",
       x_label = "Pseudotime", y_label = "Value", pt_ext = list(dom = c(0, 1), jitter = 0.012, lim = c(-0.012, 1.012))),
  list(...))

test_that("the trajectory scatter is one keyed trace with the trends and bands of the chosen kind", {
  skip_if_not_installed("plotly")
  env <- .plot_figures_env()
  df <- .traj_df()
  b <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df)))
  keyed <- .keyed(b$x$data)
  expect_length(keyed, 1)
  expect_equal(sort(.keys(keyed[[1]])), sort(df$link_key))
  expect_equal(length(unlist(keyed[[1]]$marker$color)), nrow(df))
  bands <- Filter(function(t) identical(t$fill, "toself"), b$x$data)
  lines <- Filter(function(t) identical(t$mode, "lines") && is.null(t$fill), b$x$data)
  expect_length(bands, 3)
  expect_length(lines, 3)
  expect_equal(b$x$layout$xaxis$range, c(-0.012, 1.012))
  expect_false(isTRUE(b$x$layout$showlegend))
  # The overall trend: one band and one line in the ink color; no trend at all when asked.
  b <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df, trend = "overall")))
  lines <- Filter(function(t) identical(t$mode, "lines") && is.null(t$fill), b$x$data)
  expect_length(lines, 1)
  expect_equal(lines[[1]]$line$color, "#000000")
  expect_length(Filter(function(t) identical(t$fill, "toself"), b$x$data), 1)
  b <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df, trend = "none", show_ci = FALSE)))
  expect_length(Filter(function(t) identical(t$mode, "lines"), b$x$data), 0)
})

test_that("the trajectory scatter sizes points by cell count on ggplot2's area scale and is deterministic", {
  skip_if_not_installed("plotly")
  env <- .plot_figures_env()
  df <- .traj_df()
  b <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df, size_by = "total_cells", point_size = 2)))
  size <- unlist(.keyed(b$x$data)[[1]]$marker$size)
  expect_length(size, nrow(df))
  expect_equal(min(size), 0.3 * 2 * 96 / 25.4)
  expect_equal(max(size), 3.0 * 2 * 96 / 25.4)
  first <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df)))$x$data[[1]]
  second <- plotly::plotly_build(env$trajectory_scatter_figure(df, .traj_opts(df)))$x$data[[1]]
  expect_equal(unlist(first$x), unlist(second$x))
  expect_equal(unlist(first$y), unlist(second$y))
  # Jittered, within the padded domain, and never piled on the endpoints.
  expect_true(all(unlist(first$x) >= -0.012 & unlist(first$x) <= 1.012))
  expect_false(any(unlist(first$x) %in% c(0, 1)))
})
