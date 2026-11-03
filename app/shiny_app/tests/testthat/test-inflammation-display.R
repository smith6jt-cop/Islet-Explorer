inflammation_display_env <- function() {
  env <- new.env(parent = globalenv())
  wanted <- c("PHENOTYPE_RETIRED", "phenotype_retired", "phenotype_keep", "phenotype_is_kept",
              "phenotype_is_conflict", "phenotype_partition_labels", "PHENOTYPE_GROUP_AGGREGATES",
              "is_phenotype_aggregate", "MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN",
              "is_macrophage_label", "inflammation_display_rows", "inflammation_fit_display")
  for (file in c("R/00_globals.R", "R/mod_inflammation_server.R")) {
    for (expr in parse(.find_source(file))) {
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
          as.character(expr[[2]]) %in% wanted) eval(expr, env)
    }
  }
  env
}

test_that("Inflammation presentation preserves fitted statistics and source totals", {
  env <- inflammation_display_env()
  profile <- "CD68+CD163+Macrophage_LGALS3+CD4+"
  fit <- list(
    res = data.frame(cell_type = c("B cell", "Conflict Myeloid", "Hematopoietic CD45+", profile),
                     PValue = c(0.003, 0.006, 0.04, 0.2), FDR = c(0.012, 0.012, 0.053333, 0.2)),
    dropped = c("Conflict Endocrine", "Macrophage", "CD4 Tcell"),
    n_tested = 4L, n_dropped = 3L, n_cells = 1000L, dispersion = 0.12,
    group = factor(c("ND", "T1D"))
  )
  original <- serialize(fit, NULL)
  shown <- env$inflammation_fit_display(fit)
  expect_identical(shown$res, fit$res[c(1, 4), , drop = FALSE])
  expect_identical(shown$dropped, "CD4 Tcell")
  retained <- setdiff(names(fit), c("res", "dropped"))
  expect_identical(shown[retained], fit[retained])
  expect_identical(serialize(fit, NULL), original)
})

test_that("Inflammation stale selections and two-block export projection hide forbidden labels", {
  env <- inflammation_display_env()
  scores <- data.frame(
    phenotype = c("B cell", " conflict_t ", "Endocrine", "Hematopoietic CD45-", "All T cells"),
    score = seq_len(5), n_cells = rep(1000L, 5),
    in_family = c(TRUE, FALSE, FALSE, FALSE, FALSE), q_kw = c(0.017, NA, NA, NA, NA)
  )
  original <- serialize(scores, NULL)
  expect_identical(env$inflammation_display_rows(scores), scores[1, , drop = FALSE])
  expect_identical(env$inflammation_display_rows(scores, c("B cell", " conflict_t ", "Endocrine")),
                   scores[1, , drop = FALSE])
  expect_equal(nrow(env$inflammation_display_rows(scores, " conflict_t ")), 0)
  expect_equal(nrow(env$inflammation_display_rows(scores, character(0))), 0)
  expect_identical(serialize(scores, NULL), original)
})
