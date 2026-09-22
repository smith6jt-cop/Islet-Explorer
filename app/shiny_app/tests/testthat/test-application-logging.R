.capture_log_lines <- function(expr) {
  lines <- character(0)
  withCallingHandlers(
    invisible(force(expr)),
    message = function(cnd) {
      lines <<- c(lines, strsplit(conditionMessage(cnd), "\n", fixed = TRUE)[[1]])
      invokeRestart("muffleMessage")
    }
  )
  lines
}

test_that("application events use one physical JSON line and a fixed schema", {
  codes <- names(.ISLET_LOG_EVENT_LEVELS)
  lines <- .capture_log_lines(for (code in codes) islet_log_event(code))
  expect_length(lines, length(codes))
  expect_false(any(grepl("[\r\n]", lines)))

  events <- lapply(lines, jsonlite::fromJSON)
  fields <- c("ts", "level", "event", "component", "result")
  expect_true(all(vapply(events, function(event) {
    identical(names(event), c(fields, if (event$event == "input_validation_failed") {
      "error_code"
    })) &&
      grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$", event$ts) &&
      event$level %in% c("INFO", "WARN", "ERROR") &&
      identical(event$component, "shiny") &&
      event$result %in% c("ok", "failed", "rejected")
  }, logical(1))))
  expect_identical(vapply(events, `[[`, character(1), "event"), codes)
})

test_that("malformed event codes cannot inject log records or raw values", {
  canary <- "SECRET_CANARY\r\n{\"event\":\"forged\"}\t\033[31m"
  invalid <- list(canary, "", NULL, NA_character_, 42, c("app_loaded", canary))
  lines <- .capture_log_lines(for (code in invalid) islet_log_event(code, canary))
  expect_length(lines, length(invalid))
  expect_false(any(grepl("SECRET_CANARY|forged|[\r\n\t\033]", lines)))
  expect_true(all(vapply(lapply(lines, jsonlite::fromJSON), function(event) {
    identical(event$event, "application_event_invalid") &&
      identical(event$level, "WARN") && identical(event$result, "failed") &&
      is.null(event$error_code)
  }, logical(1))))

  lines <- .capture_log_lines({
    islet_log_event("release_data_validation_failed", canary)
    islet_log_event("input_validation_failed", canary)
  })
  expect_length(lines, 2L)
  expect_false(any(grepl("SECRET_CANARY|forged|[\r\n\t\033]", lines)))
  expect_null(jsonlite::fromJSON(lines[[1]])$error_code)
  expect_identical(jsonlite::fromJSON(lines[[2]])$error_code, "invalid_input")
})

test_that("production debug rejection can log before any package loads", {
  withr::local_envvar(c(ISLET_RUNTIME_MODE = "production",
                      DEBUG_CREDENTIALS = "1", VIEWER_DEBUG = "0"))
  env <- new.env(parent = baseenv())
  env$library <- function(...) stop("Unexpected package loading")
  sys.source(.find_source("R/00_audit.R"), envir = env)
  lines <- .capture_log_lines(expect_error(
    sys.source(.find_source("R/00_globals.R"), envir = env),
    "Production does not permit credential or viewer debug modes.", fixed = TRUE
  ))
  expect_length(lines, 1L)
  event <- jsonlite::fromJSON(lines[[1]])
  expect_identical(event$event, "runtime_profile_rejected")
  expect_identical(event$level, "ERROR")
})

test_that("production lookup failures keep schema and source details out of logs", {
  withr::local_envvar(c(ISLET_RUNTIME_MODE = "production"))
  root <- tempfile("SOURCE_CANARY-")
  dir.create(root)
  withr::defer(unlink(root, recursive = TRUE))
  env <- new.env(parent = globalenv())
  env$data_path <- function(...) file.path(root, ...)

  missing <- .capture_log_lines(
    sys.source(.find_source("R/segmentation_helpers.R"), envir = env)
  )
  expect_length(missing, 1L)
  expect_identical(jsonlite::fromJSON(missing[[1]])$event,
                   "islet_spatial_lookup_missing")

  write.csv(data.frame(SCIENTIFIC_CANARY = "private-value"),
            file.path(root, "islet_spatial_lookup.csv"), row.names = FALSE)
  lines <- .capture_log_lines(value <- env$load_islet_spatial_lookup())
  expect_null(value)
  expect_length(lines, 1L)
  event <- jsonlite::fromJSON(lines[[1]])
  expect_identical(event$event, "islet_spatial_lookup_schema_invalid")
  expect_identical(event$level, "WARN")
  expect_false(any(grepl("SOURCE_CANARY|SCIENTIFIC_CANARY|private-value", c(missing, lines))))
})

test_that("production spatial failures discard raw exceptions and context", {
  env <- new.env(parent = globalenv())
  env$data_path <- function(...) file.path(tempdir(), ...)
  env$islet_is_production <- function() TRUE
  sys.source(.find_source("R/spatial_helpers.R"), envir = env)
  canary <- "SOURCE_CANARY\nscientific-row and raw exception"
  lines <- .capture_log_lines(failure <- tryCatch(
    env$.cd_backend_failure(canary, simpleError(canary)), error = identity
  ))
  expect_identical(conditionMessage(failure), "Spatial data operation failed.")
  expect_length(lines, 1L)
  expect_identical(jsonlite::fromJSON(lines[[1]])$event, "spatial_backend_failed")
  expect_false(grepl("SOURCE_CANARY|scientific-row|raw exception", lines[[1]]))
})
