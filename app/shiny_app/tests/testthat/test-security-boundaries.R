.capture_messages <- function(expr) {
  messages <- character(0)
  withCallingHandlers(
    invisible(force(expr)),
    message = function(cnd) {
      messages <<- c(messages, sub("\n$", "", conditionMessage(cnd)))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

test_that("validation failures are one fixed-schema line without rejected values", {
  malicious <- "../donor\nforged-log-line"
  lines <- .capture_messages(
    islet_valid_scalar_choice(malicious, c("safe"), "spatial_donor_invalid")
  )

  expect_length(lines, 1L)
  event <- jsonlite::fromJSON(lines[[1]])
  expect_identical(event$event, "input_validation_failed")
  expect_identical(event$result, "rejected")
  expect_identical(event$error_code, "spatial_donor_invalid")
  expect_false(grepl("../donor", lines[[1]], fixed = TRUE))
  expect_false(grepl("forged-log-line", lines[[1]], fixed = TRUE))
  expect_false(grepl("[\r\n]", lines[[1]]))
})

test_that("unknown validation codes fail closed", {
  lines <- .capture_messages(islet_log_validation_failure("attacker-controlled"))
  expect_length(lines, 1L)
  event <- jsonlite::fromJSON(lines[[1]])
  expect_identical(event$error_code, "invalid_input")
  expect_false(grepl("attacker-controlled", lines[[1]], fixed = TRUE))
})

test_that("security flags recognize explicit true values", {
  expect_true(islet_flag_enabled("1"))
  expect_true(islet_flag_enabled("TRUE"))
  expect_true(islet_flag_enabled(" yes "))
  expect_false(islet_flag_enabled("0"))
  expect_false(islet_flag_enabled("FALSE"))
  expect_false(islet_flag_enabled("unexpected"))
})

test_that("trajectory mode failures retain their fixed allowlisted code", {
  lines <- .capture_messages(islet_log_validation_failure("trajectory_mode_invalid"))
  event <- jsonlite::fromJSON(lines[[1]])
  expect_identical(event$error_code, "trajectory_mode_invalid")
})

test_that("browser choices require exact authoritative membership", {
  expect_true(islet_valid_scalar_choice("6505", c("112", "6505"),
                                        "spatial_donor_invalid"))
  expect_true(islet_valid_choice_vector(c("ND", "T1D"),
                                        c("ND", "Aab+", "T1D"),
                                        "plot_groups_invalid"))
  expect_false(suppressMessages(islet_valid_choice_vector(
    c("ND", "bogus"), c("ND", "Aab+", "T1D"), "plot_groups_invalid"
  )))
})

test_that("numeric browser controls require a finite bounded scalar", {
  expect_true(islet_valid_bounded_number(25, 1, 200,
                                         "plot_min_cells_invalid"))
  expect_false(suppressMessages(islet_valid_bounded_number(
    Inf, 1, 200, "plot_min_cells_invalid"
  )))
  expect_false(suppressMessages(islet_valid_bounded_number(
    c(10, 20), 1, 200, "plot_min_cells_invalid"
  )))
  expect_false(suppressMessages(islet_valid_bounded_number(
    201, 1, 200, "plot_min_cells_invalid"
  )))
})

test_that("viewer selections cannot traverse outside an image root", {
  images <- c("sample_0112.ome.tif", "sample 6505.ome.tiff")
  expect_true(islet_valid_asset_basename(images[[1]], images))
  expect_true(islet_valid_asset_basename(images[[2]], images))
  expect_false(suppressMessages(islet_valid_asset_basename(
    "../sample_0112.ome.tif", images
  )))
  expect_false(suppressMessages(islet_valid_asset_basename(
    "folder\\sample_0112.ome.tif", images
  )))
  expect_false(suppressMessages(islet_valid_asset_basename(
    "not-listed.ome.tif", images
  )))
})

test_that("path segments reject separators and dot traversal", {
  expect_true(islet_safe_path_segment("Islet_284"))
  expect_true(islet_safe_path_segment("6505"))
  expect_false(islet_safe_path_segment("../6505"))
  expect_false(islet_safe_path_segment("6505/part"))
  expect_false(islet_safe_path_segment("6505\\part"))
})

.find_repo_file <- function(rel) {
  candidates <- c(
    rel,
    file.path("..", "..", rel),
    file.path("..", "..", "..", "..", rel),
    file.path("app", "shiny_app", rel)
  )
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) stop("Could not locate repository file: ", rel)
  normalizePath(hits[[1]])
}

test_that("production configuration rejects debug modes and limits worker logs", {
  entrypoint <- paste(readLines(.find_repo_file("deploy/container-entrypoint.sh"),
                                warn = FALSE), collapse = "\n")
  config <- paste(readLines(.find_repo_file("deploy/shiny-server.conf"),
                            warn = FALSE), collapse = "\n")

  expect_match(entrypoint, "DEBUG_CREDENTIALS:-0", fixed = TRUE)
  expect_match(entrypoint, "VIEWER_DEBUG:-0", fixed = TRUE)
  expect_match(config, "preserve_logs false;", fixed = TRUE)
  expect_match(config, "log_file_mode 0640;", fixed = TRUE)
  expect_match(config, "sanitize_errors true;", fixed = TRUE)
})

test_that("viewer markup has an accessible title and no identifier console log", {
  viewer <- paste(readLines(.find_source("R/mod_viewer_server.R"), warn = FALSE),
                  collapse = "\n")
  globals <- paste(readLines(.find_source("R/00_globals.R"), warn = FALSE),
                   collapse = "\n")

  expect_match(viewer, "title = \"Interactive pancreatic tissue image viewer\"",
               fixed = TRUE)
  expect_false(grepl("[VIEWER] sel:", viewer, fixed = TRUE))
  expect_false(grepl("console.log", globals, fixed = TRUE))
})

test_that("AI diagnostics never emit credential paths, endpoints, or raw errors", {
  helpers <- paste(readLines(.find_source("R/ai_helpers.R"), warn = FALSE),
                   collapse = "\n")
  server <- paste(readLines(.find_source("R/mod_ai_assistant_server.R"),
                            warn = FALSE), collapse = "\n")

  expect_match(helpers, "if (islet_is_production())", fixed = TRUE)
  expect_false(grepl("message(\"[OPENAI]", helpers, fixed = TRUE))
  expect_false(grepl("API Base:", helpers, fixed = TRUE))
  expect_false(grepl("Current working directory:", helpers, fixed = TRUE))
  expect_false(grepl("Home directory:", helpers, fixed = TRUE))
  expect_false(grepl("message(\"[OPENAI] Assistant request failed:",
                     server, fixed = TRUE))
  expect_false(grepl("sprintf(\"The assistant couldn't reach %s.",
                     server, fixed = TRUE))
})

test_that("Plot consumers use only centralized validated selection state", {
  plot_server <- paste(readLines(.find_source("R/mod_plot_server.R"),
                                 warn = FALSE), collapse = "\n")
  validator_end <- regexpr("# ---- Donor palette helper", plot_server,
                           fixed = TRUE)[[1]]
  expect_gt(validator_end, 0L)
  consumers <- substring(plot_server, validator_end)

  expect_match(plot_server, "selection_vocabularies <- reactive", fixed = TRUE)
  expect_match(plot_server, "validated_selection <- reactive", fixed = TRUE)
  expect_match(plot_server, "validated_metric <- reactive", fixed = TRUE)
  expect_match(plot_server, "validated_min_cells <- reactive", fixed = TRUE)
  expect_false(grepl("input\\$(mode|which|metric|min_cells|groups)", consumers))
})
