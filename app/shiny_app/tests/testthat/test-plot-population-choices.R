test_that("Plot Cell Populations applies the short phenotype list to every region", {
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  env$phenotype_is_kept <- function(x) {
    !(as.character(x) %in% c("Endocrine", "Acinar", "Endothelial", "T cell",
                              "Immune", "Unknown", "ECAD+", "Structural", "APCs"))
  }
  env$phenotype_is_unclassified <- function(x) rep(FALSE, length(x))
  env$extract_phenotype_from_column <- function(x) {
    stem <- sub("^(peri_prop_|prop_)", "", x)
    gsub("_", " ", gsub("plus", "+", stem, fixed = TRUE), fixed = TRUE)
  }
  env$phenotype_text <- function(x) x
  sys.source(.find_source("R/plot_helpers.R"), envir = env)

  pairs <- data.frame(
    stem = c("Beta_cell", "Unknown", "Macrophage"),
    label = c("Beta cell", "Unknown", "Macrophage"),
    stringsAsFactors = FALSE
  )
  choices <- env$plot_composition_selection_choices(c(
    "prop_Beta cell", "prop_Unknown", "prop_APCs", "prop_Macrophage",
    "peri_prop_Beta_cell", "peri_prop_Unknown", "peri_prop_APCs",
    "peri_prop_Macrophage"
  ), pairs)

  expect_named(choices, c(
    "Core — Phenotype Proportions",
    "Peri-Islet — Phenotype Proportions",
    "Core+Peri — Phenotype Proportions"
  ))
  expect_setequal(names(choices[["Core — Phenotype Proportions"]]),
                  c("Beta cell", "Macrophage"))
  expect_setequal(names(choices[["Peri-Islet — Phenotype Proportions"]]),
                  c("Beta cell", "Macrophage"))
  expect_setequal(names(choices[["Core+Peri — Phenotype Proportions"]]),
                  c("Beta cell", "Macrophage"))
  expect_false(any(c("Unknown", "APCs", "Ins_any", "Glu_any", "Stt_any") %in%
                   c(unlist(lapply(choices, names)), unlist(choices, use.names = FALSE))))
})
