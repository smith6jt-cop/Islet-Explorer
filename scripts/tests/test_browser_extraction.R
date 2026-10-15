#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(testthat))
source("scripts/extract_browser_scripts.R")

test_that("R syntax is inspected without evaluating arbitrary expressions", {
  fixture <- tempfile(fileext = ".R")
  canary <- tempfile()
  writeLines(c(sprintf("system(%s)", deparse(paste("touch", canary))),
               'tags$script(HTML("document.body.textContent = \'safe\';"))',
               'style <- "body { color: red; }"'), fixture)
  result <- browser_extract_file(fixture, "fixture.R")
  expect_false(file.exists(canary))
  expect_length(result, 1L)
  expect_identical(result[[1L]]$origin, "fixture.R")
  expect_identical(result[[1L]]$line_start, 2L)
  expect_match(result[[1L]]$code, "document.body", fixed = TRUE)
  unlink(fixture)
})

test_that("nested literal functions and sprintf templates retain auditable transforms", {
  fixture <- tempfile(fileext = ".R")
  writeLines('f <- function() htmlwidgets::onRender(p, sprintf("function(el,x){ Shiny.setInputValue(\'%s\', %d); var p=\'100%%\'; }", id, n))', fixture)
  record <- browser_extract_file(fixture)[[1L]]
  expect_true(record$anonymous_function_wrapped)
  expect_true(record$sprintf_placeholders_replaced)
  expect_match(record$code, "Shiny.setInputValue('ISLET_DYNAMIC_VALUE', 0)", fixed = TRUE)
  expect_match(record$code, "100%", fixed = TRUE)
  js <- tempfile(fileext = ".js")
  writeLines(record$code, js)
  expect_identical(system2("node", c("--check", shQuote(js))), 0L)
  unlink(c(fixture, js))
})

test_that("malformed R and malformed extracted JavaScript fail validation", {
  fixture <- tempfile(fileext = ".R")
  writeLines('f <- function(', fixture)
  expect_error(browser_extract_file(fixture))
  writeLines('x <- "document.body = ;"', fixture)
  result <- browser_extract_file(fixture)
  js <- tempfile(fileext = ".js")
  writeLines(result[[1L]]$code, js)
  expect_gt(system2("node", c("--check", shQuote(js)), stdout = FALSE, stderr = FALSE), 0L)
  unlink(c(fixture, js))
})

test_that("long multiline literals are recovered from source without parser display truncation", {
  fixture <- tempfile(fileext = ".R")
  script <- paste0("document.body.textContent = '", paste(rep("safe", 400), collapse = ""), "';\nwindow.ready = true;")
  writeLines(paste0("x <- ", deparse(script, width.cutoff = 500L)), fixture)
  record <- browser_extract_file(fixture)[[1L]]
  expect_identical(record$code, script)
  expect_false(record$sprintf_placeholders_replaced)
  unlink(fixture)
})
