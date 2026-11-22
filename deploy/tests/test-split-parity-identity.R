#!/usr/bin/env Rscript
# Provenance-only synthetic fixtures; no real data or optional packages loaded.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- dirname(dirname(dirname(normalizePath(script, mustWork = TRUE))))
env <- new.env(parent = globalenv())
sys.source(file.path(repo, "scripts", "verify_release_parity.R"), envir = env)
commit <- function(char) paste(rep(char, 40L), collapse = "")
app <- commit("a")
data <- commit("b")
other <- commit("c")
options <- list(app_commit = app, data_commit = data, release_root = "/synthetic/release")
manifest <- list(git_commit = data, git_dirty = FALSE)
checks <- 0L
check <- function(value) {
  stopifnot(isTRUE(value))
  checks <<- checks + 1L
}
reject <- function(expr, code) {
  error <- tryCatch(force(expr), error = identity)
  check(inherits(error, "parity_failure") && identical(conditionMessage(error), code))
}
identify <- function(opts = options, validator = app, dirty = FALSE,
                     producer = manifest, changed = character()) {
  env$parity_source_identity(opts, validator, dirty, producer, changed)
}

# Different app/data commits are accepted only when the exact pair is supplied.
result <- identify()
check(identical(result$identity_mode, "independent_app_data"))
check(identical(result$app_git_commit, app))
check(identical(result$validator_git_commit, app))
check(identical(result$data_producer_git_commit, data))
check(identical(manifest$git_commit, data))  # Never relabel the producer.
reject(identify(validator = other), "app_commit_mismatch")
reject(identify(producer = list(git_commit = other, git_dirty = FALSE)), "data_commit_mismatch")
reject(identify(dirty = TRUE), "uncommitted_source_changes")
reject(identify(dirty = NA), "uncommitted_source_changes")
reject(identify(producer = list(git_commit = data, git_dirty = TRUE)), "uncommitted_source_changes")
reject(identify(producer = list(git_commit = data)), "uncommitted_source_changes")

# Flags must be paired, full lowercase hashes, and used with an actual release.
reject(identify(opts = list(app_commit = app)), "split_commit_pair_required")
reject(identify(opts = list(data_commit = data)), "split_commit_pair_required")
for (bad in c("HEAD", substr(app, 1L, 7L), toupper(app), paste0(app, "x"), "")) {
  bad_options <- options
  bad_options$app_commit <- bad
  reject(identify(opts = bad_options), "invalid_explicit_commit")
  bad_options <- options
  bad_options$data_commit <- bad
  reject(identify(opts = bad_options), "invalid_explicit_commit")
}
without_release <- options
without_release$release_root <- NULL
reject(identify(opts = without_release), "split_release_root_required")

# Existing same-commit and verifier-only legacy behavior is retained.
legacy <- identify(opts = list(), validator = data)
check(identical(legacy$identity_mode, "legacy"))
check(identical(legacy$data_producer_git_commit, data))
reject(identify(opts = list()), "source_commit_mismatch")
legacy_fix <- identify(opts = list(), changed = "scripts/verify_release_parity.R")
check(identical(legacy_fix$identity_mode, "legacy"))
reject(identify(opts = list(), changed = c("scripts/verify_release_parity.R", "app/shiny_app/app.R")),
       "source_commit_mismatch")
reject(identify(opts = list(), validator = data, dirty = TRUE), "uncommitted_source_changes")
reject(identify(opts = list(), validator = data, producer = list(git_commit = data, git_dirty = TRUE)),
       "uncommitted_source_changes")

# CLI rejects malformed explicit identities before loading any scientific data.
reject(env$parity_main(c("--app-commit", app, "--output-prefix", "/synthetic/evidence")),
       "split_commit_pair_required")
reject(env$parity_main(c("--app-commit", "HEAD", "--data-commit", data,
                        "--release-root", "/synthetic/release", "--output-prefix", "/synthetic/evidence")),
       "invalid_explicit_commit")
cat("Split parity identity:", checks, "checks passed. Scientific/hash checks remain unchanged.\n")
