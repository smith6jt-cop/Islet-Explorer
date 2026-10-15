source(.find_source("R/mod_plot_server.R"))

test_that("schema-4 Plot defaults to a measured typed population", {
  comp <- data.frame(
    Ins_any = c(NA_real_, NA_real_),
    Glu_any = c(NA_real_, NA_real_),
    Stt_any = c(NA_real_, NA_real_),
    check.names = FALSE
  )
  comp[["prop_Alpha cell"]] <- c(0.2, 0.3)
  comp[["prop_Beta cell"]] <- c(0.6, 0.5)
  comp[["peri_prop_Alpha cell"]] <- c(0.1, 0.2)
  comp[["peri_prop_Beta cell"]] <- c(0.4, 0.3)
  pd <- list(comp = comp, combined_pairs = data.frame())

  choices <- plot_composition_selection_choices(
    pd, typed_release = TRUE,
    app_labels = c("Alpha cell", "Beta cell", "Unavailable reviewed type")
  )
  values <- unname(unlist(choices, use.names = FALSE))

  expect_equal(values[[1]], "prop_Beta cell")
  expect_true(all(c(
    "prop_Alpha cell", "prop_Beta cell", "prop_Unavailable reviewed type",
    "peri_prop_Alpha cell", "peri_prop_Beta cell",
    "peri_prop_Unavailable reviewed type"
  ) %in% values))
  expect_false(any(c("Ins_any", "Glu_any", "Stt_any") %in% values))
})

test_that("unavailable Plot identities produce missing values without errors", {
  comp <- data.frame(existing = 1:3)
  expect_equal(plot_numeric_column(comp, "prop_Unavailable"), rep(NA_real_, 3))
  expect_false(plot_column_has_values(comp, "prop_Unavailable"))
})
