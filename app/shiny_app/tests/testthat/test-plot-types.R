# The Plot tab's plot types and Composition by donor card (owner, 2026-10-21): the pure tables behind them.

.plot_types_env <- function() {
  env <- new.env(parent = globalenv())
  defs <- c("PHENOTYPE_RETIRED", "phenotype_is_kept", "phenotype_keep", "phenotype_retired",
            "phenotype_is_conflict", "phenotype_partition_labels", "PHENOTYPE_GROUP_AGGREGATES",
            "is_phenotype_aggregate", "MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN", "is_macrophage_label",
            "UNCLASSIFIED_KEYS", "UNCLASSIFIED_PARENT", "UNCLASSIFIED_RAW_TO_KEY", "unclassified_collapse")
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% defs) eval(expr, env)
  }
  sys.source(.find_source("R/plot_types.R"), envir = env)
  sys.source(.find_source("R/donor_composition.R"), envir = env)
  env
}

.plot_rows <- function() {
  data.frame(`Case ID` = c(1, 1, 1, 2, 2, 3, 3, 3, 3), islet_key = paste0("Islet_", 1:9),
             donor_status = c("ND", "ND", "ND", "T1D", "T1D", "Aab+", "Aab+", "Aab+", "Aab+"),
             value = c(0, 10, 20, 5, NA, 1, 2, 3, 4), age = c(10, 10, 10, 20, 20, 30, 30, 30, 30),
             HbA1c_percent = c(0, 0, 0, 6, 6, 5, 5, 5, 5), check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("donor medians are one row per donor in status order, over finite values", {
  env <- .plot_types_env()
  med <- env$donor_median_table(.plot_rows(), extra = "age")
  expect_equal(med$case_id, c("1", "3", "2"))
  expect_equal(med$donor_status, c("ND", "Aab+", "T1D"))
  expect_equal(med$n_objects, c(3, 4, 1))
  expect_equal(med$median, c(10, 2.5, 5))
  expect_equal(med$q1, c(5, 1.75, 5))
  expect_equal(med$age, c(10, 30, 20))
})

test_that("the cumulative table reaches 1 within each status", {
  env <- .plot_types_env()
  tab <- env$ecdf_table(.plot_rows())
  expect_equal(unique(tab$donor_status), c("ND", "Aab+", "T1D"))
  expect_equal(tab$cum_prop[tab$donor_status == "ND"], c(1, 2, 3) / 3)
  expect_equal(tab$value[tab$donor_status == "ND"], c(0, 10, 20))
  expect_equal(tapply(tab$cum_prop, tab$donor_status, max)[["T1D"]], 1)
})

test_that("a donor with fewer than three objects has no curve; curves peak at 1", {
  env <- .plot_types_env()
  curves <- env$donor_density_curves(.plot_rows(), n = 50)
  expect_equal(vapply(curves, `[[`, "", "case_id"), c("1", "3", "2"))
  expect_null(curves[[3]]$curve)
  expect_equal(max(curves[[1]]$curve$y), 1)
  expect_equal(nrow(curves[[2]]$curve), 50)
})

test_that("HbA1c 0 is not recorded; Spearman needs four pairs", {
  env <- .plot_types_env()
  expect_equal(env$plot_clinical_values(.plot_rows(), "HbA1c_percent")[1:4], c(NA, NA, NA, 6))
  expect_null(env$plot_spearman(1:3, 3:1))
  expect_equal(env$plot_spearman(1:5, c(2, 4, 6, 8, 10))$rho, 1)
})

test_that("x values are looked up per object and proportions given in %", {
  env <- .plot_types_env()
  comp <- data.frame(`Case ID` = c(1, 1), islet_key = c("Islet_2", "Islet_1"), `prop_Beta cell` = c(0.25, 0.5),
                     check.names = FALSE)
  rows <- .plot_rows()[1:3, ]
  expect_equal(env$plot_x_values(rows, comp, "prop_Beta cell"), c(50, 25, NA))
})

test_that("donor composition pools cells over objects and closes each bar with No subtype", {
  env <- .plot_types_env()
  comp <- data.frame(`Case ID` = c(1, 1, 2), `Donor Status` = c("ND", "ND", "T1D"),
                     total_cells_core = c(10, 30, 20), total_cells_peri = c(100, 0, NA),
                     `prop_Beta cell` = c(0.5, 0.5, 0), `prop_Immune` = c(0.1, 0, 0.5),
                     peri_prop_Beta_cell = c(0.01, 0, NA), peri_prop_Immune = c(0.2, 0, NA), check.names = FALSE)
  pairs <- data.frame(label = c("Beta cell", "Immune"), core_col = c("prop_Beta cell", "prop_Immune"),
                      peri_col = c("peri_prop_Beta_cell", "peri_prop_Immune"), stringsAsFactors = FALSE)
  core <- env$donor_composition_table(comp, pairs, "core")
  d1 <- core[core$case_id == "1", ]
  expect_equal(d1$phenotype, c("Beta cell", env$DONOR_COMPOSITION_OTHER))  # Immune is a general class: in No subtype
  expect_equal(d1$cells, c(20, 20))
  expect_equal(d1$proportion, c(50, 50))
  expect_equal(unique(d1$n_objects), 2)
  both <- env$donor_composition_table(comp, pairs, "both")
  b1 <- both[both$case_id == "1", ]
  expect_equal(b1$total_cells, c(140, 140))
  expect_equal(b1$cells[b1$phenotype == "Beta cell"], 21)
  expect_equal(sum(both$proportion[both$case_id == "2"]), 100)
  expect_equal(env$donor_composition_order(core, "beta"), c("1", "2"))
})

test_that("donor categories count leaves once and never turn hidden QC cells into No subtype", {
  env <- .plot_types_env()
  labels <- c("Beta cell", "T cell", "CD4 Tcell", "Hematopoietic CD45+", "Hematopoietic CD45-",
              "Macrophage CD68+CD163-", "Macrophage CD68-CD163+", "Macrophage", "Conflict B/T", "All T cells")
  core <- c(30, 10, 20, 5, 5, 10, 10, 20, 10, 30)
  peri <- c(20, 20, 20, 30, 20, 20, 20, 40, 50, 40)
  comp <- data.frame(`Case ID` = "1", `Donor Status` = "ND", total_cells_core = 100,
                     total_cells_peri = 200, check.names = FALSE)
  pairs <- data.frame(label = labels, core_col = paste0("prop_", labels),
                     peri_col = paste0("peri_prop_", labels))
  for (i in seq_along(labels)) {
    comp[[pairs$core_col[i]]] <- core[i] / 100
    comp[[pairs$peri_col[i]]] <- peri[i] / 200
  }
  original <- comp
  tab <- env$donor_composition_table(comp, pairs, "both")
  expect_identical(tab$phenotype, c(labels[1], env$DONOR_COMPOSITION_OTHER, labels[c(3, 6, 7)]))
  expect_equal(tab$cells, c(50, 90, 40, 30, 30))
  expect_true(all(tab$total_cells == 300))
  expect_equal(sum(tab$proportion), 80) # 60 hidden cells stay in the original denominator.
  shown <- env$donor_composition_shown(tab)
  expect_equal(sum(shown$proportion), 100)
  expect_true(all(shown$shown_cells == 240))
  # Frozen releases may have no pair for labels with a sanitized slash. Hidden
  # counts still come from the complete source regional columns.
  without_conflict_pair <- pairs[pairs$label != "Conflict B/T", , drop = FALSE]
  expect_equal(env$donor_composition_table(comp, without_conflict_pair, "both"), tab)
  core_tab <- env$donor_composition_table(comp, without_conflict_pair, "core")
  peri_tab <- env$donor_composition_table(comp, without_conflict_pair, "peri")
  expect_equal(core_tab$cells[core_tab$phenotype == env$DONOR_COMPOSITION_OTHER], 20)
  expect_equal(peri_tab$cells[peri_tab$phenotype == env$DONOR_COMPOSITION_OTHER], 70)
  expect_equal(sum(core_tab$cells), 90)
  expect_equal(sum(peri_tab$cells), 150)
  expect_identical(comp, original)
})

test_that("the fitted axis ends above the 99th percentile only when a few values stretch it", {
  env <- .plot_types_env()
  bulk <- c(rep(0, 60), seq(0.5, 5, length.out = 139))
  expect_null(env$plot_axis_upper(bulk))
  stretched <- c(bulk, 100)
  upper <- env$plot_axis_upper(stretched)
  expect_true(upper < 10 && upper > 4)
  expect_equal(env$plot_axis_upper(stretched, keep = 50), 52.5)
  # 99 % zeros: the axis fits the nonzero values' median.
  expect_equal(env$plot_axis_upper(c(rep(0, 300), 2, 4, 100)), 4.4)
  expect_null(env$plot_axis_upper(1:10))
})

test_that("included phenotypes share one displayed denominator without changing source proportions", {
  env <- .plot_types_env()
  tab <- data.frame(case_id = "1", donor_status = "ND", phenotype = c("Beta cell", "Macrophage", env$DONOR_COMPOSITION_OTHER),
                    cells = c(90, 6, 4), total_cells = 100, proportion = c(90, 6, 4), n_objects = 2)
  original <- tab
  included <- c("Macrophage", env$DONOR_COMPOSITION_OTHER)
  shown <- env$donor_composition_shown(tab, included)
  expect_equal(shown$phenotype, included)
  expect_equal(shown$proportion, c(60, 40))
  expect_equal(shown$proportion_all, c(6, 4))
  expect_equal(shown$shown_cells, c(10, 10))
  expect_identical(tab, original)

  # The CSV shares the plotted rows but keeps the historical original-denominator column.
  exported <- env$donor_composition_export(shown)
  expect_equal(exported$phenotype, included)
  expect_equal(exported$proportion, c(6, 4))
  expect_equal(exported$proportion_selected, c(60, 40))
  expect_equal(exported$total_cells, c(100, 100))
  expect_equal(exported$shown_cells, c(10, 10))

  # All and None are distinct, including a well-formed empty CSV.
  expect_equal(env$donor_composition_shown(tab)$proportion, tab$proportion)
  empty <- env$donor_composition_shown(tab, character(0))
  expect_equal(nrow(empty), 0L)
  expect_equal(nrow(env$donor_composition_export(empty)), 0L)
  expect_true(all(c("shown_cells", "proportion_selected") %in% names(env$donor_composition_export(empty))))

  choices <- c(env$DONOR_COMPOSITION_ENDOCRINE, "Exocrine", included)
  expect_equal(env$donor_composition_default_include("core", choices), c("Exocrine", included))
  expect_equal(env$donor_composition_default_include("peri", choices), c(env$DONOR_COMPOSITION_ENDOCRINE, included))
  expect_equal(env$donor_composition_default_include("both", choices), included)
})

test_that("composition inclusion normalizes each donor independently and handles zero selected cells", {
  env <- .plot_types_env()
  tab <- data.frame(case_id = rep(c("1", "2", "3"), each = 2), donor_status = rep(c("ND", "Aab+", "T1D"), each = 2),
    phenotype = rep(c("Beta cell", "Macrophage"), 3), cells = c(90, 10, 25, 25, 10, 0),
    total_cells = rep(c(100, 50, 10), each = 2), proportion = c(90, 10, 50, 50, 100, 0), n_objects = 1)
  shown <- env$donor_composition_shown(tab, "Macrophage")
  expect_equal(shown$proportion, c(100, 100, NA_real_))
  expect_equal(shown$shown_cells, c(10, 25, 0))
  expect_equal(shown$proportion_all, c(10, 50, 0))
  expect_equal(shown$total_cells, c(100, 50, 10))
})

test_that("composition bars render only included rows and use the selected palette", {
  skip_if_not_installed("plotly")
  env <- .plot_types_env()
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% c("ISLET_FONT_SIZES", "plotly_islet_fonts")) eval(expr, env)
  }
  env$plotly_islet_config <- function(p, ...) p
  env$phenotype_text <- identity
  tab <- data.frame(case_id = "1", donor_status = "ND", phenotype = c("Beta cell", "Macrophage", "B cell"),
    cells = c(90, 6, 4), total_cells = 100, proportion = c(90, 6, 4), n_objects = 2)
  shown <- env$donor_composition_shown(tab, c("Macrophage", "B cell"))
  render <- function(pal) plotly::plotly_build(env$donor_composition_figure(
    shown, "1", pal, c(ND = "#3366CC"), "core"))
  first <- render(c(Macrophage = "#112233", `B cell` = "#445566"))
  bars <- Filter(function(trace) identical(trace$type, "bar"), first$x$data)
  expect_equal(vapply(bars, `[[`, "", "name"), c("Macrophage", "B cell"))
  expect_equal(vapply(bars, function(trace) as.numeric(trace$y), numeric(1)), c(60, 40))
  expect_equal(vapply(bars, function(trace) trace$marker$color, ""), c("#112233", "#445566"))
  second <- render(c(Macrophage = "#778899", `B cell` = "#AABBCC"))
  next_bars <- Filter(function(trace) identical(trace$type, "bar"), second$x$data)
  expect_equal(vapply(next_bars, function(trace) trace$marker$color, ""), c("#778899", "#AABBCC"))
  expect_equal(first$x$layout$yaxis$title$text, "Proportion of selected phenotypes")
  expect_false(grepl("<br>|<sup>", first$x$layout$title$text))
})

test_that("a measure nonzero in under 1% of objects is not offered", {
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  env$phenotype_is_kept <- function(x) TRUE
  env$extract_phenotype_from_column <- function(x) sub("^(peri_prop_|prop_)", "", x)
  env$phenotype_text <- function(x) x
  sys.source(.find_source("R/plot_helpers.R"), envir = env)
  comp <- data.frame(`prop_Beta cell` = c(rep(0, 150), 0.5, 0.2), `prop_B cell` = c(rep(0, 151), 0.5),
                     peri_prop_Beta_cell = 0, check.names = FALSE)
  pairs <- data.frame(stem = "Beta_cell", label = "Beta cell", core_col = "prop_Beta cell",
                      peri_col = "peri_prop_Beta_cell")
  choices <- env$plot_composition_selection_choices(names(comp), pairs, comp = comp)
  expect_equal(unname(choices[["Core — Phenotype Proportions"]]), "prop_Beta cell")
  expect_null(choices[["Peri-Islet — Phenotype Proportions"]])
  expect_equal(unname(choices[["Core+Peri — Phenotype Proportions"]]), "combined_prop_Beta_cell")
})
