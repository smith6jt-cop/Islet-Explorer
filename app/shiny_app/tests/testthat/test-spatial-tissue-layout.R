# Measure actual raster viewports: dense legends must add height instead of
# taking width from the fixed-aspect tissue map.
.spatial_layout_env <- function() {
  env <- new.env(parent = globalenv())
  env$ISLET_FONT_SIZES <- list(legend_text = 16)
  env$phenotype_text <- identity
  sys.source(.find_source("R/mod_spatial_server.R"), envir = env)
  env
}

.spatial_layout_plot <- function(env, labels, width) {
  cells <- data.frame(X_centroid = rep(c(0, 6000), length(labels)),
                      Y_centroid = rep(c(0, 4000), length(labels)),
                      phenotype = rep(labels, each = 2))
  ggplot2::ggplot() +
    ggplot2::geom_point(data = cells, ggplot2::aes(X_centroid, Y_centroid, colour = phenotype)) +
    ggplot2::scale_colour_manual(values = stats::setNames(scales::hue_pal()(length(labels)), labels),
      breaks = labels, labels = function(x) env$spatial_tissue_legend_labels(x, width),
      guide = ggplot2::guide_legend(ncol = env$spatial_tissue_legend_columns(width), byrow = TRUE)) +
    ggplot2::coord_fixed() + ggplot2::scale_y_reverse() +
    ggplot2::labs(title = "Donor 112 (ND)", subtitle = "Foreground and background cells", colour = "Phenotype") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
      legend.title.position = "top", legend.text = ggplot2::element_text(size = 16),
      plot.margin = ggplot2::margin(r = 24, t = 5.5, b = 5.5, l = 5.5))
}

.spatial_raster_dimensions <- function(env, plot, width) {
  height <- env$spatial_tissue_plot_height(plot, width)
  ragg::agg_capture(width = width, height = height, res = 72)
  on.exit(grDevices::dev.off(), add = TRUE)
  grob <- ggplot2::ggplotGrob(plot)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  panel <- grob$layout[grob$layout$name == "panel", ]
  grid::seekViewport(paste0("panel.", panel$t, "-", panel$l, "-", panel$b, "-", panel$r))
  map_width <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE) * 72
  map_height <- grid::convertHeight(grid::unit(1, "npc"), "in", valueOnly = TRUE) * 72
  guide <- grob$grobs[[which(grob$layout$name == "guide-box-bottom")]]
  legend_width <- sum(grid::convertWidth(guide$widths, "in", valueOnly = TRUE)) * 72
  c(height = height, map_width = map_width, map_height = map_height, legend_width = legend_width)
}

test_that("responsive dense tissue legends preserve the full-width map and fit the raster", {
  env <- .spatial_layout_env()
  basic <- c("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative")
  dense <- c(basic, "Unclassified — Insufficient typing evidence",
    "Unclassified — Insufficient evidence after reassessment",
    "Unclassified — Lineage or subtype unresolved",
    "Unclassified — Competing compartment or lineage evidence",
    paste("Macrophage", c("CD68⁺CD163⁺LGALS3⁺CD4⁺", "CD68⁺CD163⁺LGALS3⁻CD4⁺",
      "CD68⁺CD163⁻LGALS3⁺CD4⁺", "CD68⁻CD163⁻LGALS3⁻CD4⁺")))
  for (width in c(600, 850, 1200)) {
    one <- .spatial_raster_dimensions(env, .spatial_layout_plot(env, basic, width), width)
    all <- .spatial_raster_dimensions(env, .spatial_layout_plot(env, dense, width), width)
    expect_equal(unname(one["map_width"]), unname(all["map_width"]), tolerance = 1)
    expect_gt(unname(all["height"]), unname(one["height"]))
    expect_gt(unname(all["map_width"]), width - 140)
    expect_lte(unname(all["legend_width"]), width)
    expect_equal(unname(all["map_height"] / all["map_width"]), 2 / 3, tolerance = 0.001)
  }
})

test_that("initial hidden-output widths and empty tissue extents remain valid", {
  env <- .spatial_layout_env()
  for (width in list(NULL, 0, NA_real_, Inf)) {
    expect_true(is.finite(env$spatial_tissue_width(width)))
    expect_gt(env$spatial_tissue_width(width), 0)
  }
  empty <- ggplot2::ggplot() + ggplot2::theme_minimal()
  expect_gt(env$spatial_tissue_plot_height(empty, 600), 0)
  expect_equal(env$spatial_tissue_map_aspect(empty), 0.75)
})


test_that("zoomed tissue panels size from the visible spatial aspect", {
  env <- .spatial_layout_env()
  labels <- c("Beta cell", "Alpha cell")
  for (width in c(600, 850)) {
    original <- .spatial_layout_plot(env, labels, width)
    zoomed <- suppressMessages(original + ggplot2::coord_fixed(
      xlim = c(1000, 5000), ylim = c(1000, 2000)))
    one <- .spatial_raster_dimensions(env, original, width)
    zoom <- .spatial_raster_dimensions(env, zoomed, width)
    expect_equal(env$spatial_tissue_map_aspect(zoomed), 1 / 4)
    expect_gt(unname(zoom["map_width"]), width - 140)
    expect_equal(unname(zoom["map_height"] / zoom["map_width"]), 1 / 4, tolerance = 0.001)
    expect_lt(unname(zoom["height"]), unname(one["height"]))
  }
})
