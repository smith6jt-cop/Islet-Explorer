# Donor-level islet immune measures behind the Spatial tab's Neighborhood cards (owner review, 2026-10-20).

.immune_env <- function() {
  env <- new.env(parent = globalenv())
  env$leiden_levels <- function(x) {
    lv <- unique(as.character(x[!is.na(x)]))
    lv[order(suppressWarnings(as.numeric(lv)), lv)]
  }
  sys.source(.find_source("R/spatial_islet_immune.R"), envir = env)
  for (expr in parse(.find_source("R/plot_types.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% c("PLOT_AXIS_QUANTILE", "plot_axis_upper")) eval(expr, env)
  }
  env
}

.immune_comp <- function() {
  data.frame(`Case ID` = c(1, 1, 1, 2, 2, 3), `Donor Status` = c("ND", "ND", "ND", "T1D", "T1D", "T1D"),
             immune_count_core = c(0, 10, 0, 20, 1, 0), immune_count_peri = c(2, 5, 14, 0, 1, 15),
             immune_frac_peri = c(0.01, 0.05, 0.10, 0.30, 0.02, 0.04),
             `prop_Beta cell` = c(0.5, 0, 0.05, 0.6, 0, 0), `prop_Alpha cell` = c(0.5, 1, 0.95, 0.4, 1, 1),
             `prop_Delta cell` = 0, `prop_CHGA+ hormone-negative` = 0,
             leiden_0.5 = c("0", "1", "1", "0", NA, "10"),
             peri_count_CD8a_Tcell = c(0, 2, 4, 1, 0, 3), total_cells_peri = c(100, 100, 200, 50, 50, 100),
             check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("enrichment donor means retain finite extremes before summarizing", {
  env <- .immune_env()
  values <- c(1, 9, 100, -20, 4, NA_real_, Inf)
  donors <- c("1", "1", "1", "2", "2", "2", "2")
  result <- env$enrichment_donor_means(values, donors, c(rep("ND", 3), rep("T1D", 4)))
  expect_identical(result$case_id, c("1", "2"))
  expect_identical(result$donor_status, c("ND", "T1D"))
  expect_equal(result$value, c(mean(c(1, 9, 100)), mean(c(-20, 4))))
  expect_null(env$enrichment_donor_means(c(NA_real_, Inf), c("1", "2"), c("ND", "T1D")))
})

test_that("enrichment axis range changes display only and handles both tails", {
  env <- .immune_env()
  values <- c(seq(-2, 2, length.out = 100), -100, 100)
  summaries <- c(-3, 0, 3)
  axis <- env$enrichment_axis_limits(values, summaries, "bulk")
  expect_true(axis$lower < min(summaries) && axis$upper > max(summaries))
  expect_true(axis$lower > min(values) && axis$upper < max(values))
  expect_lt(axis$range[1], axis$lower)
  expect_gt(axis$range[2], axis$upper)
  expect_null(env$enrichment_axis_limits(values, summaries, "all"))
  expect_null(env$enrichment_axis_limits(1:10, c(2, 3), "bulk"))
  expect_identical(values[c(101, 102)], c(-100, 100))
})

test_that("v20 all-immune gradient includes unnamed CD45-positive cells but not CD45-negative or conflicts", {
  env <- .immune_env()
  labels <- c("Hematopoietic CD45+", "T cell", "Hematopoietic CD45-", "Conflict B/T", "Exocrine")
  profile <- list(
    bins = data.frame(case_id = "1", donor_status = "ND", phenotype = labels,
                      bin = "0-20", n = c(8, 2, 10, 10, 70)),
    tissue = data.frame(case_id = "1", donor_status = "ND", phenotype = labels,
                        n = c(2, 3, 20, 20, 55)))
  g <- env$immune_gradient(profile, env$IMMUNE_LABELS)
  expect_equal(g$count, 10)
  expect_equal(g$share, 0.10)
  expect_equal(g$tissue_share, 0.05)
  expect_equal(g$value, 2)
})

test_that("insulitis prevalence is each donor's share of objects at or over the threshold", {
  env <- .immune_env()
  d <- env$insulitis_prevalence(.immune_comp(), 15)
  expect_identical(d$case_id, c("1", "2", "3"))
  expect_equal(d$value, c(1 / 3, 1 / 2, 1))
  expect_identical(as.character(d$status), c("ND", "T1D", "T1D"))
  expect_equal(env$insulitis_prevalence(.immune_comp(), 100)$value, c(0, 0, 0))
})

test_that("band immune share is split by beta-cell content per donor", {
  env <- .immune_env()
  d <- env$immune_by_beta_content(.immune_comp(), 0.1)
  one <- d[d$case_id == "1", ]
  expect_equal(one$value[one$kind == "Insulin-containing Islets"], 0.01)
  # 0.05 beta share is under the cut: insulin-deficient with the beta-free object.
  expect_equal(one$value[one$kind == "Insulin-deficient Islets"], mean(c(0.05, 0.10)))
  expect_equal(one$objects[one$kind == "Insulin-deficient Islets"], 2)
  expect_false(any(d$case_id == "3" & d$kind == "Insulin-containing Islets"))
})

test_that("cluster membership counts objects by status and donor, clusters in numeric order", {
  env <- .immune_env()
  m <- env$cluster_membership(.immune_comp(), "leiden_0.5")
  expect_identical(levels(m$donors$cluster), c("0", "1", "10"))
  expect_equal(m$donors$donors, c(2, 1, 1))
  expect_equal(m$donors$top_donor_share, c(0.5, 1, 1))
  zero <- m$by_status[m$by_status$cluster == "0", ]
  expect_equal(sum(zero$share), 1)
})

test_that("the immune gradient and own-tissue z divide by each donor's own tissue share", {
  env <- .immune_env()
  profile <- list(
    bins = data.frame(case_id = c("1", "1", "1", "1"), donor_status = "ND",
                      phenotype = c("CD8a Tcell", "Exocrine", "CD8a Tcell", "Exocrine"),
                      bin = c("0-20", "0-20", ">400", ">400"), n = c(2, 8, 1, 19)),
    tissue = data.frame(case_id = c("1", "1", "2", "2", "3"), donor_status = c("ND", "ND", "T1D", "T1D", "T1D"),
                        phenotype = c("CD8a Tcell", "Exocrine", "CD8a Tcell", "Exocrine", "Exocrine"),
                        n = c(5, 95, 10, 90, 100)))
  share <- env$donor_tissue_share(profile, "CD8a Tcell")
  expect_equal(unname(share[c("1", "2", "3")]), c(0.05, 0.10, 0))
  g <- env$immune_gradient(profile, "CD8a Tcell")
  expect_equal(g$value[g$bin == "0-20"], 0.2 / 0.05)
  expect_equal(g$count[g$bin == "0-20"], 2)
  expect_equal(g$value[g$bin == ">400"], 0.05 / 0.05)
  z <- env$own_tissue_enrich_z(.immune_comp(), "CD8a_Tcell", "CD8a Tcell", profile)
  expect_equal(z[1], (0 - 5) / sqrt(5))
  expect_equal(z[4], (1 - 5) / sqrt(5))
  # A donor with none in its tissue has no expectation, so no z.
  expect_true(is.na(z[6]))
})
