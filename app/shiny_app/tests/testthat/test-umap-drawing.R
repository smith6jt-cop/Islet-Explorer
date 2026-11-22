# islet_umap_plotly(): one interleaved, keyed marker trace, so no color group is painted over another.

.umap_env <- function() {
  env <- new.env(parent = globalenv())
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% c("ISLET_THEME_PALETTE", "islet_theme_mode", "ISLET_UMAP_MARKER",
                                       "islet_umap_order", "islet_umap_plotly")) eval(expr, env)
  }
  env
}

.umap_df <- function() {
  data.frame(x = seq_len(30), y = rev(seq_len(30)),
             status = factor(rep(c("ND", "Aab+", "T1D"), each = 10), levels = c("ND", "Aab+", "T1D")),
             value = c(NA, seq_len(29) / 29), hover = paste0("o", seq_len(30)), key = paste0("k", seq_len(30)),
             stringsAsFactors = FALSE)
}

test_that("groups are drawn as one interleaved trace with legend-only stubs", {
  env <- .umap_env()
  pal <- c(ND = "#4477AA", "Aab+" = "#CC6633", T1D = "#228833")
  traces <- plotly::plotly_build(env$islet_umap_plotly(.umap_df(), "x", "y", "status", pal, "hover",
                                                       key = "key"))$x$data
  expect_length(traces, 4L)
  points <- traces[[1]]
  expect_false(points$showlegend)
  expect_setequal(points$key, paste0("k", 1:30))
  # Rows come grouped; the first ten drawn must already mix the groups.
  expect_length(unique(points$marker$color[1:10]), 3L)
  drawn <- .umap_df()[as.integer(sub("k", "", points$key)), ]
  expect_equal(unname(pal[as.character(drawn$status)]), as.vector(points$marker$color))
  expect_equal(drawn$x, as.vector(points$x))
  expect_identical(vapply(traces[-1], `[[`, "", "name"), c("ND", "Aab+", "T1D"))
  expect_true(all(vapply(traces[-1], function(t) is.null(t$key), TRUE)))
  expect_identical(env$islet_umap_order(30), env$islet_umap_order(30))
})

test_that("a numeric color is one trace with a color bar and no legend stubs", {
  env <- .umap_env()
  traces <- plotly::plotly_build(env$islet_umap_plotly(.umap_df(), "x", "y", "value", text = "hover",
                                                       colorbar_title = "Log cell ratio"))$x$data
  expect_length(traces, 1L)
  expect_true(traces[[1]]$marker$showscale)
  expect_equal(c(traces[[1]]$marker$cmin, traces[[1]]$marker$cmax), c(1 / 29, 1))
})

test_that("selected-feature UMAP resolves sparse values while keeping measured legend units", {
  env <- .umap_env()
  df <- data.frame(umap_1 = c(1:23, NA), umap_2 = 1:24,
                   value = c(rep(0, 15), .01, .02, .03, .04, .05, 4, Inf, NA, 100),
                   hover = "islet")
  scale <- trajectory_umap_colour(df)
  d <- scale$data
  expect_equal(nrow(d), 21L)
  expect_equal(d$value[1:15], rep(0, 15))
  expect_equal(length(unique(d$.colour)), length(unique(d$value)))
  expect_gt(d$.colour[d$value == .02] / max(d$.colour), .02 / 4)
  expect_equal(as.numeric(scale$labels[c(1, 5)]), c(0, 4))
  trace <- plotly::plotly_build(env$islet_umap_plotly(
    d, "umap_1", "umap_2", ".colour", text = "hover", colorbar_title = "Log cell ratio",
    colorbar_ticks = list(tickmode = "array", tickvals = scale$ticks, ticktext = scale$labels)
  ))$x$data[[1]]
  expect_equal(as.character(trace$marker$colorbar$ticktext), scale$labels)
  expect_true(all(grepl("^#[[:xdigit:]]{6}$", vapply(trace$marker$colorscale, `[[`, "", 2))))
})

test_that("constant or unavailable feature values cannot produce a colour gradient", {
  env <- .umap_env()
  for (value in c(0, log(2))) {
    scale <- trajectory_umap_colour(data.frame(umap_1 = 1:4, umap_2 = 4:1, value = value, h = "o"))
    traces <- plotly::plotly_build(env$islet_umap_plotly(
      scale$data, "umap_1", "umap_2", ".colour", scale$palette, "h"
    ))$x$data
    expect_null(traces[[1]]$marker$colorscale)
    expect_length(unique(traces[[1]]$marker$color), 1L)
    expect_equal(traces[[2]]$name, format(value, digits = 5, trim = TRUE))
  }
  expect_equal(nrow(trajectory_umap_colour(data.frame(umap_1 = 1, umap_2 = 1, value = NA_real_))$data), 0L)
})
