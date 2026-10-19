# islet_umap_plotly(): one interleaved, keyed marker trace, so no colour group is painted over another.

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

test_that("a numeric colour is one trace with a colour bar and no legend stubs", {
  env <- .umap_env()
  traces <- plotly::plotly_build(env$islet_umap_plotly(.umap_df(), "x", "y", "value", text = "hover",
                                                       colorbar_title = "Log cell ratio"))$x$data
  expect_length(traces, 1L)
  expect_true(traces[[1]]$marker$showscale)
  expect_equal(c(traces[[1]]$marker$cmin, traces[[1]]$marker$cmax), c(1 / 29, 1))
})
