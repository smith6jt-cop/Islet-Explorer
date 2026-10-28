# Retained-phenotype gate and marker definitions (Spatial help popover).
# Only the needed top-level definitions are evaluated: 00_globals.R loads
# packages and data paths at source time.

.load_defs <- function(file, names, env) {
  for (expr in parse(.find_source(file))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% names) eval(expr, env)
  }
  env
}

.retained_env <- function() {
  env <- new.env(parent = globalenv())
  env$PHENOTYPE_COLORS <- c("Beta cell" = 1, "Endocrine" = 1, "Ductal" = 1,
                            "Macrophage" = 1, "CD8a Tcell" = 1, "T cell" = 1,
                            "Immune" = 1, "Acinar" = 1)
  .load_defs("R/00_globals.R", c("PHENOTYPE_RETIRED", "phenotype_keep", "phenotype_retired",
                                 "phenotype_is_kept", "phenotype_marker_definitions"), env)
  .load_defs("R/spatial_helpers.R", "canonical_phenotype_order", env)
}

test_that("a top-marker root keeps ECAD+ (v17), other roots keep it retired", {
  env <- .retained_env()
  env$TOP_MARKER_BASE <- c("ECAD+", "CD3e+")
  env$root_has_top_marker_labels <- function() FALSE
  expect_false("ECAD+" %in% env$phenotype_keep(c("ECAD+", "Beta cell")))
  expect_false(env$phenotype_is_kept("ECAD+"))
  env$root_has_top_marker_labels <- function() TRUE
  expect_identical(env$phenotype_keep(c("ECAD+", "Beta cell", "Acinar")), c("ECAD+", "Beta cell"))
  expect_true(env$phenotype_is_kept("ECAD+"))
  env$root_has_top_marker_labels <- function() stop("no root")
  expect_false(env$phenotype_is_kept("ECAD+"))
})

test_that("observed retired labels do not re-enter the canonical order", {
  env <- .retained_env()
  out <- env$canonical_phenotype_order(c("Acinar", "Unknown", "Beta cell", "Novel label"))
  expect_false(any(out %in% env$PHENOTYPE_RETIRED))
  expect_equal(out, c("Beta cell", "Ductal", "Macrophage", "CD8a Tcell", "Novel label"))
})

test_that("protein rules compose parent constraints for retained labels only", {
  env <- .retained_env()
  rules <- tempfile(fileext = ".csv")
  writeLines(c(
    "parent,label,all_positive,any_positive,all_negative,note",
    "all,Endocrine,,INS|GCG|SST,,x",
    "Endocrine,Beta cell,INS,,GCG|SST,x",
    "all,Immune,,CD20|CD3e|CD68,INS|GCG|SST,x",
    "Immune,Macrophage,,CD68|CD163,CD20|CD3e,x",
    "Immune,T cell,CD3e,,CD20|CD68,x",
    "T cell,CD8a Tcell,CD8a,,CD4,x",
    "all,Ductal,,CK19|panCK,INS,x"), rules)
  defs <- env$phenotype_marker_definitions(rules)
  expect_equal(names(defs), c("Beta cell", "Ductal", "Macrophage", "CD8a Tcell"))
  expect_equal(defs[["Beta cell"]], "INS+  —  GCG−, SST−")
  expect_equal(defs[["Macrophage"]],
               "any of CD68+ / CD163+  —  CD20−, CD3e−, INS−, GCG−, SST−")
  expect_equal(defs[["CD8a Tcell"]],
               "CD8a+, CD3e+  —  CD4−, CD20−, CD68−, INS−, GCG−, SST−")
})
