test_that("donor filters leave fixed coordinates and percentiles intact", {
  d <- data.frame(islet_id = c("1_a", "2_b", "1_c"), case_id = c("1", "2", "1"),
                  measured_cells = c(20, 50, 100), pseudotime_raw = c(.1,.8,.9),
                  pseudotime_percentile = c(.05,.4,.95))
  one <- trajectory_filter(d, "1", c(10, 70))
  expect_equal(one$pseudotime_raw, .1)
  expect_equal(one$pseudotime_percentile, .05)
  expect_equal(trajectory_filter(d, "all", c(1, Inf)), d)
})

test_that("secondary-component graph layout does not invent pseudotime", {
  coords <- data.frame(islet_id=c("a","b"),case_id=c("1","2"),DC1=c(1,NA),DC2=c(2,NA),pseudotime_raw=c(.2,NA),pseudotime_percentile=c(.4,NA))
  layout <- data.frame(islet_id=c("b","a"),graph_x=c(10,1),graph_y=c(12,2),layout_kind="componentwise spectral graph layout")
  result <- trajectory_graph_coordinates(coords,layout)
  expect_equal(result$graph_x,c(1,10))
  expect_identical(result$pseudotime_raw,coords$pseudotime_raw)
  expect_true(is.na(result$pseudotime_raw[2]))
  expect_error(trajectory_graph_coordinates(coords,layout[1,,drop=FALSE]),"islet revision")
})

test_that("withdrawn coordinates use explicit fixed section-area heatmap bins", {
  coords <- data.frame(case_id = c("78", "1"), area_um2 = c(10, 1000),
                       pseudotime_raw = NA_real_, pseudotime_percentile = NA_real_)
  unavailable <- list(coordinates_available = FALSE)
  axis <- trajectory_heatmap_axis(coords, unavailable)
  expect_identical(axis$kind, "section_area")
  expect_identical(axis$column, "area_um2")
  expect_equal(axis$limits, c(1, 3))
  expect_match(axis$label, "pseudotime unavailable")
  expect_equal(axis$transform(coords$area_um2), c(1, 3))
  expect_identical(trajectory_heatmap_axis(coords, unavailable, "raw"), axis)
  expect_true(all(is.na(coords$pseudotime_raw)))
  expect_null(trajectory_heatmap_notice(data.frame(area_um2 = 10, core__INS = 1), axis$column))
  existing <- trajectory_heatmap_axis(coords, list(status = "stable"))
  expect_identical(existing$kind, "pseudotime")
  expect_equal(existing$limits, c(0, 1))
})
test_that("ratios use the classifier-core cell denominator and preserve missingness", {
  assign("PHENOTYPE_RETIRED", c("Endocrine", "Unknown"), envir = globalenv())
  d <- data.frame(islet_id = c("1_a", "1_b", "2_c"), case_id = c("1", "1", "2"))
  counts <- data.frame(islet_id = rep(d$islet_id, 2), phenotype = "CD8a Tcell",
                       compartment = rep(c("core", "peri"), each = 3), count = c(2,3,NA,4,5,2))
  denom <- data.frame(islet_id = d$islet_id, endocrine_count = c(3,0,NA))
  r <- trajectory_ratios(d, counts, "CD8a Tcell", "both", denom)
  expect_equal(r$numerator[1], 6)
  expect_equal(r$raw_ratio[1], 2)
  expect_equal(r$plotted_value[1], log1p(2))
  expect_true(all(is.na(r$raw_ratio[2:3])))
  expect_error(trajectory_ratios(d, counts, "Endocrine", "core", denom), "Retired")
})

test_that("filters with no matching islets preserve empty ratio download columns", {
  d <- data.frame(islet_id = "1_a", case_id = "1", measured_cells = 20,
                  pseudotime_raw = .1, pseudotime_percentile = .2)
  counts <- data.frame(islet_id = rep("1_a", 2), phenotype = "Beta cell",
                       compartment = c("core", "peri"), count = c(3, 1))
  denom <- data.frame(islet_id = "1_a", endocrine_count = 6)
  empty <- trajectory_filter(d, "1", c(100, 200))
  for (region in c("core", "peri", "both")) {
    r <- trajectory_ratios(empty, counts, "Beta cell", region, denom)
    expect_equal(nrow(r), 0L)
    expect_true(all(c("numerator", "denominator", "raw_ratio", "plotted_value", "phenotype", "numerator_compartment") %in% names(r)))
    expect_type(r$plotted_value, "double")
  }
  expect_equal(nrow(trajectory_ratios(empty, NULL, "Beta cell")), 0L)
  restored <- trajectory_ratios(trajectory_filter(d, "1", c(10, 200)), counts, "Beta cell", "both", denom)
  expect_equal(restored$raw_ratio, 4/6)
  expect_equal(restored$pseudotime_raw, d$pseudotime_raw)
})

test_that("islet choices clear safely when donor filters yield zero rows", {
  d <- data.frame(case_id = "78", islet_key = "AdaptiveIslet_78_00011", islet_id = "78_AdaptiveIslet_78_00011", measured_cells = 24)
  expect_identical(unname(trajectory_islet_choices(d)), "78|AdaptiveIslet_78_00011")
  expect_identical(trajectory_islet_choices(d[FALSE, ]), character(0))
  expect_identical(trajectory_islet_choices(NULL), character(0))
})

test_that("expected trajectory missingness has explicit informational states", {
  d <- data.frame(pseudotime_percentile = c(.1, .2), plotted_value = c(NA_real_, NA_real_))
  expect_match(trajectory_ratio_notice(d)$title, "ratios unavailable")
  expect_match(trajectory_ratio_notice(d[FALSE, ])$title, "No islets")
  expect_match(trajectory_ratio_notice(NULL)$title, "unavailable")
  d$pseudotime_percentile[] <- NA_real_
  expect_match(trajectory_ratio_notice(d)$title, "No supported coordinates")
  d$pseudotime_percentile[1] <- .1; d$plotted_value[1] <- log1p(2)
  expect_null(trajectory_ratio_notice(d))
})

test_that("cleared heatmap markers and empty donor ranges have informational states", {
  d <- data.frame(pseudotime_percentile = .1, core__INS = 2)
  expect_null(trajectory_heatmap_notice(d))
  expect_match(trajectory_heatmap_notice(d[FALSE, ])$title, "No supported islets")
  expect_match(trajectory_heatmap_notice(d[, "pseudotime_percentile", drop = FALSE])$title, "Choose heatmap markers")
  d$pseudotime_percentile <- NA_real_
  expect_match(trajectory_heatmap_notice(d)$title, "No supported islets")
})

test_that("all 31 markers retain panel order and clear without phantom columns", {
  markers <- paste0("M", 1:31)
  features <- as.data.frame(setNames(rep(list(1), 31), paste0("core__", markers)))
  expect_identical(trajectory_marker_panel(features, list(markers = markers)), markers)
  expect_identical(trajectory_marker_columns(markers, rev(markers)), paste0("core__", markers))
  expect_identical(trajectory_marker_columns(markers, character(0)), character(0))
  expect_identical(trajectory_marker_columns(markers, NULL), character(0))
  expect_error(trajectory_marker_panel(features, list(markers = markers[-1])), "does not match")
})

test_that("uncertified numeric counts never produce ratios and valid zeros do", {
  d <- data.frame(islet_id = c("1_a", "1_b", "1_c"), case_id = "1")
  counts <- data.frame(islet_id = rep(d$islet_id, 2), phenotype = "Beta cell",
    compartment = rep(c("core", "peri"), each = 3), count = c(0, 5, 2, 0, 1, 1),
    validated = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE), reason = c("", "gate_pending", "", "", "", ""))
  denominator <- data.frame(islet_id = d$islet_id, endocrine_count = c(4, 7, 0), validated = TRUE, reason = "")
  r <- trajectory_ratios(d, counts, "Beta cell", "both", denominator)
  expect_identical(r$plotted_value[1], 0)
  expect_true(is.na(r$plotted_value[2]))
  expect_true(is.na(r$plotted_value[3]))
  expect_match(r$availability_reason[2], "gate_pending")
  expect_identical(r$availability_reason[3], "no_core_cells_in_classifier_area")
  denominator$validated[1] <- FALSE
  denominator$reason[1] <- "CD20 review missing"
  r <- trajectory_ratios(d, counts, "Beta cell", "core", denominator)
  expect_true(is.na(r$denominator[1]))
  expect_match(r$availability_reason[1], "CD20 review missing")
})
