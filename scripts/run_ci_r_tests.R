#!/usr/bin/env Rscript
# Preserve machine-readable outcomes and reject every failure/error/skip in CI.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1L] %in% c("app", "parity", "extraction")) {
  stop("Usage: Rscript scripts/run_ci_r_tests.R app|parity|extraction OUTPUT_DIRECTORY")
}
suite <- args[1L]
dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[2L], mustWork = TRUE)
collector <- testthat::ListReporter$new()
reporter <- testthat::MultiReporter$new(list(
  testthat::SummaryReporter$new(), collector,
  testthat::JunitReporter$new(file = file.path(output, paste0(suite, "-tests.xml")))
))
if (suite == "app") {
  testthat::test_dir("app/shiny_app/tests/testthat", reporter = reporter, stop_on_failure = FALSE)
} else {
  path <- switch(suite, parity = "scripts/tests/test_release_parity.R", extraction = "scripts/tests/test_browser_extraction.R")
  testthat::with_reporter(reporter, {
    reporter$start_file(path)
    source(path, local = new.env(parent = globalenv()))
    reporter$end_file()
  })
}
results <- as.data.frame(collector$get_results())
if (!nrow(results)) stop("No R tests executed")
summary <- list(suite = suite, tests = nrow(results), assertions = sum(results$nb),
                passed = sum(results$passed), failed = sum(results$failed),
                errors = sum(results$error), skipped = sum(results$skipped))
jsonlite::write_json(summary, file.path(output, paste0(suite, "-tests.json")), pretty = TRUE, auto_unbox = TRUE)
if (summary$failed || summary$errors || summary$skipped) stop("CI requires zero failed, errored or skipped R tests")
