#!/usr/bin/env Rscript
# Verify exact producer export and independent raw/served tamper detection.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- dirname(dirname(dirname(normalizePath(script, mustWork = TRUE))))
env <- new.env(parent = globalenv())
sys.source(file.path(repo, "scripts", "verify_release_parity.R"), envir = env)
checks <- 0L
check <- function(value) { stopifnot(isTRUE(value)); checks <<- checks + 1L }
reject <- function(expr, code = NULL) {
  error <- tryCatch(force(expr), error = identity)
  check(inherits(error, "error") && (is.null(code) || identical(conditionMessage(error), code)))
}
fixture <- tempfile("producer-fixture-")
dir.create(fixture)
local({
  on.exit(unlink(fixture, recursive = TRUE, force = TRUE), add = TRUE)
  git <- function(...) {
    args <- c("-C", shQuote(fixture), ...)
    output <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) stop("fixture_git_failed")
    output
  }
  git("init", "--quiet")
  dir.create(file.path(fixture, "app"))
  dir.create(file.path(fixture, "scripts"))
  writeLines("source_value <- 7", file.path(fixture, "app", "values.R"))
  verifier <- file.path(fixture, "scripts", "verify_release_parity.R")
  writeLines(c(
    'parity_runtime <- function(repo, root) {',
    '  out <- new.env(parent = globalenv())',
    '  sys.source(file.path(repo, "app", "values.R"), envir = out)',
    '  out$h5ad_path <- file.path(root, "unused.h5ad")',
    '  out$load_master_h5ad <- function(path) list()',
    '  out$prep_data <- function(master) list(comp = data.frame(key = "synthetic", count = out$source_value))',
    '  out$islet_disconnect_backend <- function() invisible(TRUE)',
    '  out',
    '}',
    'if (sys.nframe() == 0L) stop("HISTORICAL_MAIN_MUST_NOT_RUN")'
  ), verifier)
  git("add", "app", "scripts")
  git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "--quiet", "-m", "fixture")
  revision <- git("rev-parse", "HEAD")[[1]]
  stages <- function() sort(list.files(tempdir(), pattern = "^parity-producer-", full.names = TRUE))
  before <- stages()
  # Dirty producer working files cannot enter the exported committed artifact.
  writeLines("source_value <- 999", file.path(fixture, "app", "values.R"))
  result <- env$parity_producer_prepared(fixture, revision, tempdir())
  check(identical(result$prepared$comp$count, 7))
  check(identical(result$producer_commit, revision))
  check(grepl("^[0-9a-f]{64}$", result$producer_factory_sha256))
  check(identical(before, stages()))
  actual <- result$prepared$comp
  check(identical(env$parity_compare_frames(result$prepared$comp, actual, "key")$max_absolute_difference, 0))
  actual$count <- 8
  reject(env$parity_compare_frames(result$prepared$comp, actual, "key"), "numeric_mismatch")

  reject(env$parity_producer_prepared(fixture, "HEAD", tempdir()), "invalid_producer_commit")
  reject(env$parity_producer_prepared(fixture, paste(rep("f", 40), collapse = ""), tempdir()), "producer_export_failed")
  check(identical(before, stages()))
  writeLines("not_a_factory <- TRUE", verifier)
  git("add", "scripts")
  git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "--quiet", "-m", "missing-factory")
  reject(env$parity_producer_prepared(fixture, git("rev-parse", "HEAD")[[1]], tempdir()), "producer_factory_missing")
  check(identical(before, stages()))
})

# The current shared entrypoint is invoked in each actual runtime mode. The
# adapter preserves its full output and restores process mode on success/error.
previous <- Sys.getenv("ISLET_RUNTIME_MODE", unset = NA_character_)
local({
  on.exit(if (is.na(previous)) Sys.unsetenv("ISLET_RUNTIME_MODE") else
    Sys.setenv(ISLET_RUNTIME_MODE = previous), add = TRUE)
  Sys.setenv(ISLET_RUNTIME_MODE = "previous")
  runtime <- new.env(parent = globalenv())
  modes <- character()
  runtime$load_shared_prepared_data <- function(root) {
    mode <- Sys.getenv("ISLET_RUNTIME_MODE")
    modes <<- c(modes, mode)
    list(comp = data.frame(key = "synthetic", count = if (mode == "production") 9 else 7))
  }
  source <- env$parity_served_prepared(runtime, "/synthetic/source", "development")
  actual <- env$parity_served_prepared(runtime, "/synthetic/release", "production")
  check(identical(modes, c("development", "production")))
  reject(env$parity_compare_frames(source$comp, actual$comp, "key"), "numeric_mismatch")
  check(identical(Sys.getenv("ISLET_RUNTIME_MODE"), "previous"))
  runtime$load_shared_prepared_data <- function(root) stop("intentional_failure")
  reject(env$parity_served_prepared(runtime, "/synthetic/release", "production"), "intentional_failure")
  check(identical(Sys.getenv("ISLET_RUNTIME_MODE"), "previous"))
})
cat("Split parity boundaries:", checks, "checks passed.\n")
