source(.find_source("R/trajectory_marker_helpers.R"))

typed_heatmap_fixture <- function() {
  coords <- data.frame(case_id = c("1", "1", "2"), islet_id = c("1_a", "1_b", "2_c"),
    measured_cells = c(10, 30, 100), pseudotime_percentile = c(.2, .3, .4), pseudotime_raw = c(.7,.8,.9))
  summaries <- expand.grid(islet_id = coords$islet_id, compartment = c("core", "peri"), marker = c("INS", "CD45"), stringsAsFactors = FALSE)
  summaries$case_id <- coords$case_id[match(summaries$islet_id, coords$islet_id)]
  summaries$population <- "All cells"
  summaries$n_cells <- ifelse(summaries$compartment == "core", 10, 2)
  summaries$mean_log1p <- ifelse(summaries$compartment == "core", 1, 4)
  summaries$positive_count <- ifelse(summaries$compartment == "core", 2, 1)
  summaries$positive_available <- TRUE
  denominator <- data.frame(case_id = coords$case_id, islet_id = coords$islet_id, endocrine_count = c(3, 0, 8), validated = TRUE)
  list(coords = coords, summaries = summaries, denominator = denominator)
}

test_that("typed heatmap exports exact positive ratios and retains zero denominators missing", {
  f <- typed_heatmap_fixture()
  d <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, c("INS", "CD45"))
  expect_equal(d$positive_count, rep(3, 6))
  expect_equal(d$raw_ratio[d$islet_id == "1_a"], c(1,1))
  expect_equal(d$value[d$islet_id == "1_a"], rep(log1p(1), 2))
  expect_true(all(is.na(d$value[d$islet_id == "1_b"])))
  expect_true(all(d$availability_reason[d$islet_id == "1_b"] == "no_core_cells_in_classifier_area"))
  one <- trajectory_filter(f$coords, "2", c(10, Inf))
  selected <- trajectory_typed_marker_values(one, f$summaries, f$denominator, "INS")
  expect_identical(selected$case_id, "2")
  expect_equal(selected$pseudotime_percentile, .4)
  expect_equal(selected$pseudotime_raw, .9)
})

test_that("typed expression weights compartments by cells and preserves unavailable calls", {
  f <- typed_heatmap_fixture()
  f$summaries$positive_available[f$summaries$marker == "INS"] <- FALSE
  ratio <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, "INS")
  expect_true(all(is.na(ratio$value)))
  expression <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, "INS", "typed_expression")
  expect_equal(expression$value, rep(1.5, 3))
  # The alternate measurement is available even when endocrine denominator is zero.
  expect_equal(expression$value[2], 1.5)
  f$summaries$n_cells[f$summaries$compartment == "peri"] <- 0
  f$summaries$mean_log1p[f$summaries$compartment == "peri"] <- NA_real_
  expression <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, "INS", "typed_expression")
  expect_equal(expression$value, rep(1,3))
  f$summaries <- f$summaries[!(f$summaries$islet_id == "1_a" & f$summaries$compartment == "peri"), ]
  expression <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, "INS", "typed_expression")
  expect_true(is.na(expression$value[1]))
})

test_that("heatmap bins weight donors equally and retain unavailable panel rows", {
  d <- data.frame(case_id = c("1","1","1","2"), islet_id = letters[1:4], marker = "INS", value = c(0,0,0,10), pseudotime_percentile = .2)
  axis <- list(column = "pseudotime_percentile", transform = identity, limits = c(0,1))
  bins <- trajectory_heatmap_bin_values(d, axis, c("INS","CD45"), bins = 2)
  expect_equal(bins$value[bins$marker == "INS" & bins$bin == 1], 5)
  expect_equal(bins$n_donors[bins$marker == "INS" & bins$bin == 1], 2)
  expect_equal(bins$n_islets[bins$marker == "INS" & bins$bin == 1], 4)
  expect_true(all(is.na(bins$value[bins$marker == "CD45"])))
  expect_equal(nrow(bins), 4)
})

test_that("cleared markers and duplicate summaries cannot produce phantom values", {
  f <- typed_heatmap_fixture()
  empty <- trajectory_typed_marker_values(f$coords, f$summaries, f$denominator, character())
  expect_equal(nrow(empty), 0)
  expect_true(all(c("marker", "value", "positive_count", "denominator") %in% names(empty)))
  expect_error(trajectory_typed_marker_values(f$coords, rbind(f$summaries, f$summaries[1, ]), f$denominator, "INS"), "Duplicate")
  expect_identical(trajectory_marker_panel(NULL, list(markers = paste0("M",1:31))), paste0("M",1:31))
})
