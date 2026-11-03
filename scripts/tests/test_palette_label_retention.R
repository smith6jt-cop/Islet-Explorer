# Plot-scale metadata only: dynamic labels must retain their names and legend entries.
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
drilldown <- parse("app/shiny_app/R/drilldown_helpers.R")
spatial <- parse("app/shiny_app/R/mod_spatial_server.R")
palettes <- list()
for (expr in drilldown) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      startsWith(as.character(expr[[2]]), "PHENOTYPE_COLORS") &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("c"))) {
    palettes[[as.character(expr[[2]])]] <- eval(expr[[3]])
  }
}
assignments <- list()
walk <- function(expr) {
  if (!is.call(expr) && !is.expression(expr)) return(invisible(NULL))
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("pal", "bg_pal") &&
      grepl("\\bpalette[[(]", paste(deparse(expr[[3]]), collapse = " "))) {
    assignments[[length(assignments) + 1L]] <<- expr
  }
  for (item in as.list(expr)) {
    if (missing(item)) next
    if (is.call(item) || is.expression(item)) walk(item)
  }
}
walk(drilldown); walk(spatial)
stopifnot(length(palettes) == 4L, length(assignments) == 5L)
labels <- c("Beta cell", "Unclassified")
for (base in palettes) for (assignment in assignments) {
  e <- new.env(parent = globalenv())
  e$pheno_present <- e$legend_phenos <- e$bg_phenos <- labels
  e$palette <- if (grepl("palette()", paste(deparse(assignment), collapse = " "), fixed = TRUE))
    local({ value <- base; function() value }) else base
  eval(assignment, e)
  values <- e[[as.character(assignment[[2]])]]
  values[is.na(values)] <- "#CCCCCC"
  stopifnot(identical(names(values), labels), values[["Beta cell"]] == base[["Beta cell"]],
            values[["Unclassified"]] == "#CCCCCC")
  for (scale in list(ggplot2::scale_fill_manual(values = values),
                    ggplot2::scale_colour_manual(values = values))) {
    scale$train(labels)
    stopifnot(setequal(scale$get_limits(), labels), scale$map("Unclassified") == "#CCCCCC")
  }
}
cat("All five palette lookups retain dynamic labels across all four palettes and both plot scales.\n")
