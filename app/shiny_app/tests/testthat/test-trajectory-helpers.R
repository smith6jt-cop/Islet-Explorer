source(.find_source("R/trajectory_links.R"))
source(.find_source("R/trajectory_donor_metadata.R"))

# The destiny fallback with only the app definitions its graph inputs use.
trajectory_input_env <- function() {
  env <- new.env(parent = globalenv())
  need <- c("TRAJ_MIN_CELLS", "object_class_key", "object_class_arg", "PHENOTYPE_RETIRED", "phenotype_is_conflict", "UNCLASSIFIED_KEYS", "UNCLASSIFIED_RAW_TO_KEY",
            "MACROPHAGE_LABELS", "MACROPHAGE_LABEL_PATTERN", "MACROPHAGE_KEY_PATTERN", "is_macrophage_label",
            "macrophage_profile", "unclassified_collapse", "PHENOTYPE_DISPLAY_NAMES", "phenotype_text")
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && as.character(expr[[2]])[1] %in% need) {
      eval(expr, envir = env)
    }
  }
  sys.source(.find_source("R/barlow_pseudotime.R"), envir = env)
  env
}

test_that("the default trajectory includes four-core-cell objects in DPT and PAGA", {
  env <- trajectory_input_env()
  sys.source(.find_source("R/paga_layout.R"), envir = env)
  sys.source(.find_source("R/mod_trajectory_server.R"), envir = env)
  expect_identical(env$TRAJ_MIN_CELLS, 4L)
  cohort <- data.frame(
    islet_id = c("112_Cluster_1", "112_Cluster_2", "112_Islet_3", "112_Islet_4"),
    total_cells_core = c(3L, 4L, 9L, 10L), pseudotime = c(.7, 0, .5, 1),
    stringsAsFactors = FALSE)
  # A bridge stub supplies fitted IDs and lets the assertions detect a different
  # cutoff in any R entry point. The real engine's feature/eligibility tests run in Python.
  result <- function(min_cells, max_cells, object_classes, age_weight = 0) {
    d <- cohort[cohort$total_cells_core >= min_cells, ]
    as.list(c(d, list(case_id = rep("112", nrow(d)), donor_status = rep("ND", nrow(d)),
      leiden = rep("0", nrow(d)), fa_1 = d$pseudotime, fa_2 = -d$pseudotime,
      guided_pseudotime = d$pseudotime, root_islet_id = d$islet_id[1], root_donor = "112")))
  }
  env$paga_available <- function() TRUE
  env$.paga_module <- function() list(compute_pseudotime = result, compute = result)
  shuffled <- cohort[c(4, 2, 1, 3), , drop = FALSE]
  obs <- data.frame(combined_islet_id = shuffled$islet_id)
  expect_identical(as.numeric(env$traj_pseudotime_column(obs)), c(1, 0, NA_real_, .5))
  expect_identical(env$compute_scanpy_pseudotime()$islet_id, cohort$islet_id[-1])
  paga <- env$paga_result()
  expect_identical(paga$islet_id, cohort$islet_id[-1])
  expect_identical(paga$total_cells_core, c(4L, 9L, 10L))
  expect_equal(paga$pseudotime, env$compute_scanpy_pseudotime()$pseudotime)
  expect_identical(env$paga_result(age_weight = .5)$islet_id, paga$islet_id)
  expect_identical(env$compute_scanpy_pseudotime(10L)$islet_id, "112_Islet_4")
  expect_identical(env$paga_result(10L)$islet_id, "112_Islet_4")
})

test_that("trajectory ratios offer and compute endocrine subtypes only for core", {
  env <- trajectory_input_env()
  subtypes <- env$TRAJECTORY_ENDOCRINE_CLASSES
  retained <- c("All T cells", "Exocrine", "Stromal")
  labels <- c("All endocrine cells", subtypes, retained)
  expect_identical(env$barlow_ratio_choices(labels, "core"), c(subtypes, retained))
  for (region in c("both", "peri"))
    expect_identical(env$barlow_ratio_choices(labels, region), retained)
  expect_length(env$barlow_ratio_choices(labels, "invalid"), 0)

  d <- data.frame(islet_id = c("one", "two"), n_endocrine_core = c(4, 2),
                  "n_core_Beta cell" = c(2, 1), "n_peri_Beta cell" = c(1, 0),
                  n_core_Stromal = c(1, 2), n_peri_Stromal = c(3, 1), check.names = FALSE)
  env$barlow_features <- function() d
  env$phenotype_members <- identity
  for (region in env$BARLOW_RATIO_REGIONS)
    expect_null(env$barlow_ratio_values("All endocrine cells", region))
  for (region in c("peri", "both"))
    expect_null(env$barlow_ratio_values("Beta cell", region))
  expect_equal(env$barlow_ratio_values("Beta cell", "core")$value, rep(log1p(0.5), 2))
  expect_equal(env$barlow_ratio_values("Stromal", "peri")$value, log1p(c(3/4, 1/2)))
  expect_equal(env$barlow_ratio_values("Stromal", "both")$value, log1p(c(1, 3/2)))
})

test_that("destiny graph inputs follow the current typing format", {
  env <- trajectory_input_env()
  qc <- c(names(env$UNCLASSIFIED_RAW_TO_KEY), "Conflict B/T", "Conflict_Endocrine/Exocrine")
  expect_setequal(names(env$UNCLASSIFIED_RAW_TO_KEY),
                  c("Tissue", "Unclassified", "Doublet suspected", "Hematopoietic CD45+",
                    "Hematopoietic CD45-", "T cell"))
  columns <- c("n_endocrine_core", "n_core_Beta cell", "core_Beta cell", "peri_Beta cell",
               "core_CHGA+ hormone-negative", "peri_CHGA+ hormone-negative",
               "core_Exocrine", "peri_Exocrine", "core_Acinar cell", "peri_Acinar cell",
               "core_Ductal cell", "peri_Ductal cell", "core_Stromal", "peri_Stromal",
               paste0("core_", env$PHENOTYPE_RETIRED), paste0("peri_", env$PHENOTYPE_RETIRED),
               paste0("core_", qc), paste0("peri_", qc))
  expect_identical(env$barlow_trajectory_features(columns),
                   c("core_Beta cell", "core_CHGA+ hormone-negative", "peri_Exocrine", "peri_Acinar cell",
                     "peri_Ductal cell", "core_Stromal", "peri_Stromal"))
  expect_false(any(env$trajectory_feature_allowed("core", qc)))
  expect_false(any(env$trajectory_feature_allowed("peri", qc)))
  expect_identical(env$trajectory_feature_allowed(c("core", "peri"), "Alpha cell"), c(TRUE, FALSE))
  expect_identical(env$trajectory_feature_allowed(c("core", "peri"), "Exocrine"), c(FALSE, TRUE))
})

test_that("destiny sums the older macrophage labels by Gal3 before the log", {
  env <- trajectory_input_env()
  legacy <- c("CD68+CD163+Macrophage_LGALS3+CD4+", "CD68+CD163-Macrophage_LGALS3+CD4-",
              "CD68-CD163+Macrophage_LGALS3-CD4+", "CD68+CD163+Macrophage_LGALS3-CD4-")
  d <- data.frame(n_endocrine_core = c(10, 0, 4), check.names = FALSE)
  counts <- list(c(1, 0, 2), c(3, 1, 0), c(0, 5, 1), c(2, 2, 2))
  for (region in c("core", "peri")) {
    d[[paste0("n_", region, "_Beta cell")]] <- c(9, 0, 3)
    d[[paste0(region, "_Beta cell")]] <- log1p(c(9, 0, 3) / pmax(d$n_endocrine_core, 1))
    for (i in seq_along(legacy)) {
      d[[paste0("n_", region, "_", legacy[i])]] <- counts[[i]]
      d[[paste0(region, "_", legacy[i])]] <- log1p(counts[[i]] / pmax(d$n_endocrine_core, 1))
    }
  }
  expect_error(env$barlow_trajectory_inputs(d), "Missing phenotype-specific")
  m <- env$barlow_phenotype_inputs(d)
  expect_identical(colnames(m), c("core_Beta cell", "combined_Macrophage_LGALS3+", "combined_Macrophage_LGALS3-"))
  expect_equal(m[, "combined_Macrophage_LGALS3+"], log1p(2*c(4, 1, 2) / c(10, 1, 4)), tolerance = 1e-15)
  expect_equal(m[, "combined_Macrophage_LGALS3-"], log1p(2*c(2, 7, 3) / c(10, 1, 4)), tolerance = 1e-15)
  expect_identical(unname(m[, "core_Beta cell"]), d[["core_Beta cell"]])
  expect_error(env$barlow_trajectory_inputs(d[, !grepl("^n_peri_", names(d))]), "needs their counts")
})

test_that("the provenance line names the fitted inputs by region", {
  inputs_text <- trajectory_inputs_text
  environment(inputs_text) <- trajectory_input_env()      # the app's display names
  text <- inputs_text(c("core_Beta cell", "peri_Exocrine", "core_Macrophage_LGALS3+",
                        "peri_Macrophage_LGALS3+"))
  expect_match(text, "core (2): Endocrine Beta, Macrophage Gal3⁺", fixed = TRUE)
  expect_match(text, "peri-islet band (2): Exocrine, Macrophage Gal3⁺", fixed = TRUE)
  expect_match(text, "Unclassified", fixed = TRUE)
  expect_match(inputs_text(NULL), "not reported", fixed = TRUE)
})

test_that("donor annotations summarize each donor's fitted islets on the selected scale", {
  obs <- data.frame(imageid = c("112", "112", "112", "6505", "6505"),
                    donor_status = c("ND", "ND", "ND", "Aab+", "Aab+"),
                    pseudotime = c(.03, .07, NA, .12, .48), pt_rank = c(.125, .375, NA, .625, .875))
  ranks <- trajectory_donor_summary(obs, "rank")
  expect_identical(ranks$case_id, c("112", "6505"))
  expect_identical(ranks$n_islets, c(2L, 2L))
  expect_equal(ranks$median_pt, c(.25, .75))
  expect_equal(ranks$third_quartile - ranks$first_quartile, c(.125, .125))
  expect_equal(trajectory_donor_summary(obs, "raw")$median_pt, c(.05, .30))
  expect_equal(nrow(trajectory_donor_summary(NULL, "rank")), 0)
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

test_that("PAGA clinical overlays preserve the fit and distinguish censored or missing values", {
  d <- data.frame(islet_id = paste0(1:4, "_Islet_1"), case_id = as.character(1:4),
                  donor_status = c("ND","Aab+","T1D","T1D"), leiden = "1",
                  pseudotime = c(0,.3,.7,1), fa_1 = 1:4, fa_2 = 4:1,
                  C_peptide_ng_mL = c(.02, 2, 4, NA_real_),
                  C_peptide_censored = c(TRUE, FALSE, FALSE, FALSE),
                  AAb_GADA = c(FALSE,TRUE,NA,TRUE))
  result <- paga_display_data(d, "C_peptide_ng_mL", clinical_fields = TRAJ_DONOR_FIELDS)
  expect_identical(result$data[names(d)], d)
  expect_equal(result$data$.colour, c(NA_real_, 2, 4, NA_real_))
  expect_match(result$data$.hover[1], "&lt;0.02", fixed = TRUE)
  expect_identical(paga_display_data(d, "AAb_GADA", clinical_fields = TRAJ_DONOR_FIELDS)$data$.colour,
                   c("Negative","Positive","Unavailable","Positive"))
})

test_that("the PAGA weight blends cached endpoints with shared viewport and baseline DPT", {
  env <- new.env(parent = globalenv())
  sys.source(.find_source("R/paga_layout.R"), envir = env)
  base <- data.frame(islet_id = c("112_Islet_1", "6334_Islet_2"), case_id = c("112","6334"),
                     donor_status = c("ND","T1D"), leiden = c("0","1"),
                     fa_1 = c(-1, 1), fa_2 = c(1, 0), pseudotime = c(0, 1))
  end <- list(islet_id = base$islet_id, fa_1 = c(-2, 3), fa_2 = c(0, 3), guided_pseudotime = c(-.2, 1.3))
  zero <- env$paga_interpolate_layout(base, end, 0)
  half <- env$paga_interpolate_layout(base, end, .5)
  full <- env$paga_interpolate_layout(base, end, 1)
  expect_identical(zero$fa_1, base$fa_1)
  expect_identical(zero$fa_2, base$fa_2)
  expect_identical(zero$age_guided_pseudotime, base$pseudotime)
  expect_identical(full$fa_1, end$fa_1)
  expect_equal(half$fa_2, c(.5, 1.5))
  expect_equal(half$age_guided_pseudotime, c(-.1, 1.15))
  expect_identical(half$pseudotime, base$pseudotime)
  expect_identical(attr(zero, "paga_axis_limits"), attr(full, "paga_axis_limits"))
  expect_true(all(c(base$fa_1,end$fa_1) >= attr(half, "paga_axis_limits")$fa_1[1]))
  expect_identical(paga_display_data(half, "pseudotime")$data$.colour, half$age_guided_pseudotime)
  expect_identical(paga_display_data(zero, "pseudotime")$data$.colour, base$pseudotime)
  expect_match(paga_display_data(half, "pseudotime")$title, "Age-guided")
  for (w in list(NA, Inf, -1, 1.1, 0.005, "0", c(0, 1))) expect_false(env$paga_valid_age_weight(w))
  expect_null(env$paga_interpolate_layout(base, NULL, .05))
  expect_identical(env$paga_interpolate_layout(base, NULL, 0)$fa_1, base$fa_1)
  end$islet_id <- rev(end$islet_id)
  expect_error(env$paga_interpolate_layout(base, end, 1), "baseline")
})

test_that("terminal review panel and exclusion labels never enter trajectory features", {
  env <- trajectory_input_env()
  for (region in c("core", "peri"))
    expect_identical(env$trajectory_feature_allowed(region, c("Panel undected", "Conflicting markers - excluded")), c(FALSE, FALSE))
})
