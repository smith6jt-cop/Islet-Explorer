# Presentation taxonomy is separate from immutable source labels and scientific query lists.
.unclassified_env <- function() {
  env <- new.env(parent = globalenv())
  wanted <- c("UNCLASSIFIED_KEYS", "UNCLASSIFIED_PARENT", "UNCLASSIFIED_RAW_TO_KEY", "unclassified_collapse",
    "unclassified_collapse_comp", "phenotype_collapse", "phenotype_members", "MACROPHAGE_LABELS",
    "MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN", "MACROPHAGE_MARKERS", "macrophage_signs",
    "macrophage_collapse", "is_macrophage_label", "PHENOTYPE_DISPLAY_NAMES", "phenotype_text",
    "PHENOTYPE_GROUP_AGGREGATES", "PHENOTYPE_GROUP_ORDER", "phenotype_group", "is_phenotype_aggregate",
    "phenotype_aggregate_members", "phenotype_partition_labels", "phenotype_aggregates_for",
    "add_phenotype_aggregate_columns", "phenotype_groups_split", "phenotype_grouped_choices",
    "phenotype_is_conflict", "phenotype_is_unclassified", "phenotype_selectable",
    "PHENOTYPE_RETIRED", "phenotype_retired", "phenotype_is_kept", "phenotype_keep",
    "phenotype_filter_rows", "phenotype_filter_conflict_columns", "write_islet_csv")
  for (expr in parse(.find_source("R/00_globals.R")))
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && as.character(expr[[2]]) %in% wanted) eval(expr, env)
  env$TOP_MARKER_ORDER <- character(0)
  env$root_has_top_marker_labels <- function() FALSE
  env$root_has_marker_profile_labels <- function() TRUE
  env$root_cd45_is_rule_label <- function() FALSE
  env$root_phenotype_labels <- function() c(names(env$UNCLASSIFIED_RAW_TO_KEY), "Beta cell", "CD4 Tcell", "Conflict B/T")
  env
}

test_that("the four distinct leaves retain source meaning and biological raw queries", {
  e <- .unclassified_env(); k <- unname(e$UNCLASSIFIED_KEYS)
  raw <- c("Tissue", "Unclassified", "Hematopoietic CD45+", "Hematopoietic CD45-", "T cell", "Doublet suspected")
  expect_identical(e$unclassified_collapse(raw), k[c(1, 2, 3, 3, 3, 4)])
  expect_identical(e$unclassified_collapse(c(k, "Beta cell", NA)), c(k, "Beta cell", NA))
  expect_identical(e$phenotype_members("Unclassified"), "Unclassified")
  expect_identical(e$phenotype_members("T cell"), "T cell")
  expect_setequal(e$phenotype_members(c("T cell", "Hematopoietic CD45+")), c("T cell", "Hematopoietic CD45+"))
  expect_setequal(e$phenotype_members(c(e$UNCLASSIFIED_PARENT, k, raw)), raw)
  expect_setequal(e$phenotype_members(k[3]), raw[3:5])
  expect_false("Conflict B/T" %in% e$phenotype_members(e$UNCLASSIFIED_PARENT))
  expect_equal(e$phenotype_text(c("Unclassified", e$UNCLASSIFIED_PARENT)),
    c("Unclassified \u2014 Insufficient evidence after reassessment", "Unclassified (all)"))
})

test_that("cell projections and additive selections retain each row once and exclude parent aliases", {
  e <- .unclassified_env(); k <- unname(e$UNCLASSIFIED_KEYS)
  raw <- c(names(e$UNCLASSIFIED_RAW_TO_KEY), "Conflict B/T", "Beta cell")
  df <- data.frame(phenotype = raw, cell_id = seq_along(raw)); original <- df
  out <- e$phenotype_filter_rows(df)
  expect_identical(out$cell_id, c(1:6, 8L))
  expect_identical(df, original)
  expect_identical(e$phenotype_filter_rows(out), out)
  expect_identical(e$phenotype_partition_labels(c(raw, k, "No subtype", e$UNCLASSIFIED_PARENT)), "Beta cell")
  choices <- e$phenotype_grouped_choices(c(e$UNCLASSIFIED_PARENT, raw, k))
  expect_null(choices$Unclassified)
  expect_false(any(c("Tissue", "Unclassified", "Doublet suspected", "No subtype", "T cell") %in% unlist(choices)))
  e$root_has_marker_profile_labels <- function() FALSE
  expect_identical(e$phenotype_aggregates_for(k), e$UNCLASSIFIED_PARENT)
})

test_that("object metrics conserve subtype counts, pool enrichment and use nearest distance without aliases", {
  e <- .unclassified_env(); k <- unname(e$UNCLASSIFIED_KEYS)
  raw <- names(e$UNCLASSIFIED_RAW_TO_KEY)
  safe <- function(x) gsub("+", "plus", gsub(" ", "_", x, fixed = TRUE), fixed = TRUE)
  comp <- data.frame(total_cells_core = c(100, 100, 100), source_row = 1:3)
  for (i in seq_along(raw)) {
    comp[[paste0("prop_", raw[i])]] <- c(i / 100, NA, 0)
    comp[[paste0("peri_prop_", safe(raw[i]))]] <- c(i / 100, NA, 0)
    comp[[paste0("peri_count_", safe(raw[i]))]] <- c(i, NA, 0)
    comp[[paste0("enrich_z_", safe(raw[i]))]] <- c((i - i / 2) / sqrt(i / 2), NA, -sqrt(i / 2))
    comp[[paste0("min_dist_", safe(raw[i]))]] <- c(10 - i, NA, i)
  }
  # Both stale aliases and fresh raw columns may be supplied; only raw members count.
  comp[[paste0("prop_", k[3])]] <- c(999, 999, 999)
  comp$`prop_No subtype` <- c(999, 999, 999)
  original <- comp; out <- e$unclassified_collapse_comp(comp)
  expect_identical(comp, original)
  expect_equal(out[[paste0("prop_", k[3])]], c(0.12, NA, 0))
  expect_equal(out[[paste0("peri_count_", safe(k[3]))]], c(12, NA, 0))
  expect_equal(out[[paste0("enrich_z_", safe(k[3]))]], c(6 / sqrt(6), NA, -sqrt(6)))
  expect_equal(out[[paste0("min_dist_", safe(k[3]))]], c(5, NA, 3))
  expect_equal(out[[paste0("prop_", e$UNCLASSIFIED_PARENT)]], c(0.21, NA, 0))
  expect_equal(out[[paste0("peri_count_", safe(e$UNCLASSIFIED_PARENT))]], c(21, NA, 0))
  expect_equal(out[[paste0("enrich_z_", safe(e$UNCLASSIFIED_PARENT))]], c(10.5 / sqrt(10.5), NA, -sqrt(10.5)))
  expect_equal(out[[paste0("min_dist_", safe(e$UNCLASSIFIED_PARENT))]], c(4, NA, 1))
  expect_false(any(paste0("prop_", c(raw, "No subtype")) %in% names(out)))
  expect_equal(e$unclassified_collapse_comp(out), out, tolerance = 1e-12)
  expect_equal(e$add_phenotype_aggregate_columns(out), out, tolerance = 1e-12)
  expect_identical(out$total_cells_core, original$total_cells_core)
  expect_identical(out$source_row, original$source_row)
})

test_that("exports display family names while preserving rows, measurements and source-label provenance", {
  e <- .unclassified_env(); k <- unname(e$UNCLASSIFIED_KEYS)
  input <- data.frame(phenotype = c("Tissue", "Unclassified"), cell_type = c("T cell", k[3]),
    item = c("Doublet suspected", "Beta cell"), feature_name = c("prop_Unclassified", "peri_prop_Tissue"),
    z = c(1.25, -0.4), cell_id = c("first", "second"), source_column = c("raw_a", "raw_b"),
    phenotype_source_label = c("existing-audit-a", "existing-audit-b"), check.names = FALSE)
  input[[paste0("prop_", e$UNCLASSIFIED_PARENT)]] <- c(0.15, 0.2)
  original <- input; path <- tempfile(fileext = ".csv"); on.exit(unlink(path), add = TRUE)
  e$write_islet_csv(path, input, "Existing provenance")
  output <- read.csv(path, comment.char = "#", check.names = FALSE)
  expect_identical(input, original)
  expect_identical(output$phenotype, e$phenotype_text(k[1:2]))
  expect_identical(output$cell_type, rep(e$phenotype_text(k[3]), 2))
  expect_identical(output$item, c(e$phenotype_text(k[4]), "Beta cell"))
  expect_identical(output$feature_name, paste0(c("prop_", "peri_prop_"), e$phenotype_text(k[c(2, 1)])))
  expect_identical(output$phenotype_source_label, input$phenotype_source_label)
  expect_identical(output$phenotype_source_label_original, input$phenotype)
  expect_identical(output$source_column, input$source_column)
  expect_identical(output$z, input$z)
  expect_identical(output$cell_id, input$cell_id)
  expect_equal(output[["prop_Unclassified (all)"]], c(0.15, 0.2))
  expect_true(any(grepl("Rows and measured values are not merged or recomputed", readLines(path), fixed = TRUE)))
})

test_that("wide export headers require pre-aggregated metrics instead of silently combining statistics", {
  e <- .unclassified_env(); path <- tempfile(fileext = ".csv"); on.exit(unlink(path), add = TRUE)
  df <- data.frame(`enrich_z_T_cell` = 1, `enrich_z_Hematopoietic_CD45plus` = 2, check.names = FALSE)
  expect_error(e$write_islet_csv(path, df), "prepare mutually exclusive display metrics")
  expect_false(file.exists(path))
  df <- data.frame(peri_prop_Tissue = 0.2, count_Unclassified = 5, check.names = FALSE)
  e$write_islet_csv(path, df)
  out <- read.csv(path, comment.char = "#", check.names = FALSE)
  expect_identical(names(out), c("peri_prop_Unclassified \u2014 Insufficient typing evidence",
                               "count_Unclassified \u2014 Insufficient evidence after reassessment"))
  expect_equal(unname(unlist(out)), c(0.2, 5))
  expect_true(any(grepl("Source column: peri_prop_Tissue ->", readLines(path), fixed = TRUE)))
})

test_that("All T cells retains biological membership before and after display projection", {
  e <- .unclassified_env(); subtype <- unname(e$UNCLASSIFIED_KEYS["subtype"])
  e$root_phenotype_labels <- function() c("CD4 Tcell", "CD8a Tcell", "T cell", "Hematopoietic CD45+",
                                       "Hematopoietic CD45-", "Tissue", "Conflict T/B")
  expect_setequal(e$phenotype_members("All T cells"), c("CD4 Tcell", "CD8a Tcell", "T cell"))
  expect_identical(e$phenotype_group("T cell"), "Unclassified")
  expect_identical(e$phenotype_members("T cell"), "T cell")
  expect_setequal(e$phenotype_members(c("T cell", "Hematopoietic CD45+")), c("T cell", "Hematopoietic CD45+"))
  raw <- data.frame(`prop_CD4 Tcell` = c(0.1, NA), `prop_CD8a Tcell` = c(0.2, NA),
    `prop_T cell` = c(0.3, NA), `prop_Hematopoietic CD45+` = c(0.05, NA),
    `prop_Hematopoietic CD45-` = c(0.07, NA), peri_prop_CD4_Tcell = c(0.01, NA),
    peri_prop_CD8a_Tcell = c(0.02, NA), peri_prop_T_cell = c(0.03, NA),
    peri_prop_Hematopoietic_CD45plus = c(0.04, NA), `peri_prop_Hematopoietic_CD45-` = c(0.05, NA),
    check.names = FALSE)
  original <- raw
  biological <- e$add_phenotype_aggregate_columns(raw)
  expect_equal(biological[["prop_All T cells"]], c(0.6, NA))
  expect_equal(biological$peri_prop_All_T_cells, c(0.06, NA))
  projected <- e$unclassified_collapse_comp(biological)
  expect_false("prop_T cell" %in% names(projected))
  expect_equal(projected[[paste0("prop_", subtype)]], c(0.42, NA))
  final <- e$add_phenotype_aggregate_columns(projected)
  expect_identical(final[["prop_All T cells"]], biological[["prop_All T cells"]])
  expect_identical(final$peri_prop_All_T_cells, biological$peri_prop_All_T_cells)
  expect_identical(e$add_phenotype_aggregate_columns(final), final)
  expect_identical(raw, original)
})
