# Donor-level islet immune measures behind the Spatial tab's Neighborhood cards (owner review, 2026-10-20).

.immune_env <- function() {
  env <- new.env(parent = globalenv())
  env$leiden_levels <- function(x) {
    lv <- unique(as.character(x[!is.na(x)]))
    lv[order(suppressWarnings(as.numeric(lv)), lv)]
  }
  sys.source(.find_source("R/spatial_islet_immune.R"), envir = env)
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
  expect_equal(one$value[one$kind == "Beta cells present"], 0.01)
  # 0.05 beta share is under the cut: insulin-deficient with the beta-free object.
  expect_equal(one$value[one$kind == "Insulin-deficient"], mean(c(0.05, 0.10)))
  expect_equal(one$objects[one$kind == "Insulin-deficient"], 2)
  expect_false(any(d$case_id == "3" & d$kind == "Beta cells present"))
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
  expect_equal(g$value[g$bin == ">400"], 0.05 / 0.05)
  z <- env$own_tissue_enrich_z(.immune_comp(), "CD8a_Tcell", "CD8a Tcell", profile)
  expect_equal(z[1], (0 - 5) / sqrt(5))
  expect_equal(z[4], (1 - 5) / sqrt(5))
  # A donor with none in its tissue has no expectation, so no z.
  expect_true(is.na(z[6]))
})
