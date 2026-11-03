# A palette-only filter used to leave retired/QC cells on the map as gray points.
test_that("drilldown removes hidden phenotype points without changing the loaded evidence", {
  skip_if_not_installed("ggplot2")
  env <- new.env(parent = globalenv())
  defs <- c("PHENOTYPE_RETIRED", "phenotype_keep", "phenotype_retired", "phenotype_is_kept",
            "phenotype_is_conflict")
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% defs) eval(expr, env)
  }
  for (expr in parse(.find_source("R/drilldown_helpers.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("render_islet_drilldown_plot"))) eval(expr, env)
  }
  env$build_segmentation_base_plot <- function(...) ggplot2::ggplot()
  env$theme_islet <- ggplot2::theme_minimal
  env$ISLET_FONT_SIZES <- list(legend_text = 14, legend_title = 14)
  env$PIXEL_SIZE_UM <- 0.5
  env$phenotype_text <- identity
  env$PHENOTYPE_COLORS <- c("Beta cell" = "red", "CD4 Tcell" = "blue", "Tissue" = "gray")
  cells <- data.frame(phenotype = c("Beta cell", "Conflict B/T", "T cell", "CD4 Tcell", "Tissue"),
                      X_centroid = 1:5, Y_centroid = 6:10, INS = 11:15,
                      cell_region = c("core", "core", "peri", "peri", "peri"))
  original <- cells
  info <- list(islet_key = "Islet_1")
  categorical <- env$render_islet_drilldown_plot(info, cells)
  expect_identical(categorical$layers[[1]]$data$phenotype, c("Beta cell", "CD4 Tcell", "Tissue"))
  expect_equal(categorical$layers[[1]]$data$x_px, c(2, 8, 10))
  continuous <- env$render_islet_drilldown_plot(info, cells, color_by = "INS")
  expect_equal(continuous$layers[[1]]$data$marker_val, c(11, 13, 14, 15))
  core <- env$render_islet_drilldown_plot(info, cells, show_peri = FALSE)
  expect_identical(core$layers[[1]]$data$phenotype, "Beta cell")
  expect_identical(cells, original)
})
