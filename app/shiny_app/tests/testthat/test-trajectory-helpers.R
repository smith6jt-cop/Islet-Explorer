source(.find_source("R/trajectory_links.R"))

test_that("destiny uses the app's retained phenotype policy without changing ratios", {
  env <- new.env(parent = globalenv())
  definitions <- parse(.find_source("R/00_globals.R"))
  for (expr in definitions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% c("PHENOTYPE_RETIRED", "phenotype_keep", "phenotype_retired",
                                      "phenotype_is_kept", "phenotype_is_conflict")) {
      eval(expr, envir = env)
    }
  }
  sys.source(.find_source("R/barlow_pseudotime.R"), envir = env)
  columns <- c("n_endocrine_core", "n_core_Beta cell", "core_Beta cell", "peri_Beta cell",
               "core_Ductal", paste0("core_", env$PHENOTYPE_RETIRED),
               paste0("peri_", env$PHENOTYPE_RETIRED), "core_Conflict B/T", "peri_Conflict_Endocrine/Exocrine")
  expect_identical(env$barlow_trajectory_features(columns),
                   c("core_Beta cell", "peri_Beta cell", "core_Ductal"))
})

test_that("trajectory identities preserve complete donor and islet names", {
  keys <- trajectory_link_key(c("0078", "112", "78", "6356"),
    c("Islet_1234", "Islet_1", "Islet_Islet_1", "Islet_1"))
  expect_identical(keys, c("78|Islet_1234", "112|Islet_1", "78|Islet_1", "6356|Islet_1"))
  expect_identical(trajectory_link_key(c(NA, "../78", "78"),
                                      c("Islet_1", "Islet_1", "Islet_1/../2")),
                   rep(NA_character_, 3))
})

test_that("linked point keys survive ggplotly grouping and row reordering", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("ggplot2")
  df <- data.frame(donor = c("78", "112", "78", "112"),
                   islet = c("Islet_1234", "Islet_1", "Islet_1", "Islet_1234"),
                   x = c(3, 1, 4, 2), y = c(11, 12, 13, 14))
  df$key <- trajectory_link_key(df$donor, df$islet)
  build <- function(rows) {
    g <- ggplot2::ggplot(rows, ggplot2::aes(x, y, color = donor, key = key)) +
      ggplot2::geom_point() + ggplot2::geom_hline(yintercept = 12)
    trajectory_link_plot(plotly::ggplotly(g), function(x) paste0("traj-", x))
  }
  for (rows in list(df, df[c(4, 1, 2, 3), ])) {
    p <- build(rows)
    points <- Filter(function(t) identical(t$meta$trajectory_role, "islet"), p$x$data)
    actual <- do.call(rbind, lapply(points, function(t) data.frame(key = t$key, x = t$x, y = t$y)))
    expect_setequal(actual$key, df$key)
    expect_equal(as.numeric(actual$x), df$x[match(actual$key, df$key)])
    expect_equal(as.numeric(actual$y), df$y[match(actual$key, df$key)])
    expect_length(p$preRenderHook(p)$x$shinyEvents, 0)
    expect_length(p$x$shinyEvents, 0)
    expect_false(p$x$layout$xaxis$autorange)
    expect_false(p$x$layout$yaxis$autorange)
  }
})

test_that("trajectory extent keeps jitter and panel limits synchronized", {
  ext <- traj_pt_extent(c(0, 10), jitter_frac = 0.02)
  expect_equal(ext$dom, c(0, 10))
  expect_equal(ext$jitter, 0.2)
  expect_equal(ext$lim, c(-0.2, 10.2))

  fallback <- traj_pt_extent(c(NA_real_, NA_real_))
  expect_equal(fallback$dom, c(0, 1))
})

test_that("pseudotime scale helper preserves raw values and stable ranks", {
  df <- data.frame(pt = c(4, 2, 8), pt_rank = c(0.5, 0.2, 0.9))
  raw <- apply_pt_scale(df, "raw")
  ranked <- apply_pt_scale(df, "rank")

  expect_equal(raw$pt, c(4, 2, 8))
  expect_equal(raw$pt_raw, c(4, 2, 8))
  expect_equal(ranked$pt, c(0.5, 0.2, 0.9))
  expect_equal(ranked$pt_raw, c(4, 2, 8))
})


test_that("PAGA color choices share a cohort and preserve the fitted values", {
  d <- data.frame(islet_id = c("0078_Islet_1", "112_Cluster_2", "6356_Islet_3"),
                  case_id = c("0078", "112", "6356"), donor_status = c("ND", "Aab+", "T1D"),
                  leiden = c("0", "1", "0"), fa_1 = c(1.125, -3.625, 7.75),
                  fa_2 = c(-2.25, 0.5, 4.875), pseudotime = c(0, 0.25, 1))
  original <- d
  features <- data.frame(combined_islet_id = d$islet_id[c(3, 1)], value = c(9, 2),
                         feature_name = "CD8a Tcell")
  for (by in unname(PAGA_COLOR_CHOICES)) {
    actual <- paga_display_data(d, by, features)
    expect_identical(actual$data[, names(d)], original)
    expect_identical(actual$data$.key, c("78|Islet_1", "112|Cluster_2", "6356|Islet_3"))
    expect_identical(actual$colour_by, by)
    if (by == "feature") {
      expect_equal(actual$data$.colour, c(2, NA, 9))
      expect_identical(actual$title, "CD8a Tcell")
      expect_match(actual$data$.hover[2], "Unavailable", fixed = TRUE)
    } else {
      expect_identical(actual$data$.colour, original[[by]])
    }
  }
  expect_identical(d, original)
})

test_that("missing PAGA features stay explicit instead of silently changing coloring", {
  d <- data.frame(islet_id = "112_Islet_1", case_id = "112", donor_status = "ND", leiden = "0",
                  fa_1 = 1, fa_2 = 2, pseudotime = 0.5)
  missing <- paga_display_data(d, "feature")
  expect_true(all(is.na(missing$data$.colour)))
  expect_identical(missing$colour_by, "feature")
  expect_identical(missing$data[, names(d)], d)
  invalid <- paga_display_data(d, "unrecognized")
  expect_identical(invalid$colour_by, "pseudotime")
  expect_identical(invalid$data$.colour, d$pseudotime)
})
