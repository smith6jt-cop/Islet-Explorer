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
  .load_defs("R/00_globals.R", c("UNCLASSIFIED_KEYS", "UNCLASSIFIED_PARENT", "UNCLASSIFIED_RAW_TO_KEY", "unclassified_collapse",
                                 "PHENOTYPE_DISPLAY_NAMES", "phenotype_text", "MACROPHAGE_LABEL_PATTERN",
                                 "PHENOTYPE_RETIRED", "phenotype_keep", "phenotype_retired",
                                 "phenotype_is_kept", "phenotype_is_conflict", "phenotype_filter_conflict_columns",
                                 "phenotype_filter_rows", "phenotype_marker_definitions", "write_islet_csv"), env)
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

test_that("Conflict labels stay hidden while current pan residuals project into Unclassified", {
  env <- .retained_env()
  labels <- c("Conflict B/T", "Conflict_Endocrine/Exocrine", " conflict B/Myeloid", "Hematopoietic CD45+",
              "Hematopoietic CD45-", "T cell", "CD4 Tcell", "Tissue", "Unclassified", NA, "")
  expect_identical(env$phenotype_is_conflict(labels),
    c(TRUE, TRUE, TRUE, rep(FALSE, 8)))
  expect_identical(env$phenotype_keep(labels), labels[4:9])
  frame <- data.frame(phenotype = labels, source_row = seq_along(labels))
  out <- env$phenotype_filter_rows(frame)
  expect_identical(out$source_row, 4:9)
  expect_identical(out$phenotype, c(rep(unname(env$UNCLASSIFIED_KEYS["subtype"]), 3),
    "CD4 Tcell", unname(env$UNCLASSIFIED_KEYS[c("tissue", "reassessed")])) )
  expect_identical(frame$phenotype, labels)
  expect_false(env$phenotype_is_conflict("NonConflict marker"))
})

test_that("wide app exports omit Conflict columns without changing source totals or unrelated fields", {
  env <- .retained_env()
  fields <- c("total_cells_core", "prop_Conflict B/T", "peri_prop_Conflict_B/T", "peri_count_Conflict_B/T",
              "enrich_z_Conflict_B/T", "min_dist_Conflict_B/T", "core_Conflict B/T", "peri_Conflict B/T",
              "Conflict B/T", "prop_Beta cell", "prop_T cell", "conflict_count_audit")
  frame <- as.data.frame(stats::setNames(as.list(seq_along(fields)), fields), check.names = FALSE)
  out <- env$phenotype_filter_conflict_columns(frame)
  expect_identical(names(out), c("total_cells_core", "prop_Beta cell", "prop_T cell", "conflict_count_audit"))
  expect_identical(out$total_cells_core, frame$total_cells_core)
  expect_identical(names(frame), fields)
})

test_that("identity reference delegates the Unclassified family once to its shared explanation", {
  env <- .retained_env()
  .load_defs("R/00_globals.R", "identity_phenotype_definitions_ui", env)
  env$root_has_selected_endocrine <- function() TRUE
  env$identity_structure_text <- function() "Structure-qualified cell rules."
  env$phenotype_text <- identity
  calls <- 0L
  env$unclassified_help_ui <- function() {
    calls <<- calls + 1L
    shiny::tags$p("Shared Unclassified family explanation and cohort census")
  }
  html <- as.character(env$identity_phenotype_definitions_ui())
  expect_false(grepl("Hematopoietic CD45", html, fixed = TRUE))
  expect_false(grepl("<strong>T cell</strong>", html, fixed = TRUE))
  expect_match(html, "<strong>CD4 Tcell</strong>", fixed = TRUE)
  expect_false(grepl("<strong>Unclassified</strong>", html, fixed = TRUE))
  expect_equal(calls, 1L)
  expect_match(html, "Shared Unclassified family explanation and cohort census", fixed = TRUE)
})

test_that("the shared CSV writer projects phenotype fields and retains unrelated provenance", {
  env <- .retained_env()
  frame <- data.frame(phenotype = c("Beta cell", "Conflict B/T", "Tissue"),
                      label = c("Conflict review note", "other note", "source note"),
                      total_cells = c(3, 3, 3), check.names = FALSE)
  frame[["prop_Conflict B/T"]] <- c(0.1, 0.2, 0.3)
  original <- frame
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  env$write_islet_csv(path, frame, "Source provenance retained")
  lines <- readLines(path)
  expect_true(any(grepl("omitted 1 Conflict-labeled rows and 1 Conflict phenotype columns", lines, fixed = TRUE)))
  actual <- read.csv(path, comment.char = "#", check.names = FALSE)
  expect_identical(actual$phenotype, c("Beta cell", "Unclassified \u2014 Insufficient typing evidence"))
  expect_identical(actual$phenotype_source_label, c("Beta cell", "Tissue"))
  expect_identical(actual$label, c("Conflict review note", "source note"))
  expect_identical(names(actual), c("phenotype", "label", "total_cells", "phenotype_source_label"))
  expect_identical(frame, original)
})

test_that("observed retired labels do not re-enter the canonical order", {
  env <- .retained_env()
  out <- env$canonical_phenotype_order(c("Acinar", "Unknown", "Beta cell", "Novel label"))
  expect_false(any(out %in% env$PHENOTYPE_RETIRED))
  expect_equal(out, c("Beta cell", "Ductal", "Macrophage", "CD8a Tcell",
                      unname(env$UNCLASSIFIED_KEYS["subtype"]), "Novel label"))
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
  expect_equal(names(defs), c("Beta cell", "Ductal", "Macrophage", "CD8a Tcell", "T cell"))
  expect_equal(defs[["Beta cell"]], "INS+  —  GCG−, SST−")
  expect_equal(defs[["Macrophage"]],
               "any of CD68+ / CD163+  —  CD20−, CD3e−, INS−, GCG−, SST−")
  expect_equal(defs[["CD8a Tcell"]],
               "CD8a+, CD3e+  —  CD4−, CD20−, CD68−, INS−, GCG−, SST−")
})
