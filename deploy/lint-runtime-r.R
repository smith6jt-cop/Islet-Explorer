#!/usr/bin/env Rscript

if (!requireNamespace("lintr", quietly = TRUE)) {
  stop("The lintr package is required for the deployed-runtime lint gate.")
}

runtime_root <- file.path("app", "shiny_app")
if (!dir.exists(runtime_root)) {
  stop("Run this script from the repository root; missing ", runtime_root, ".")
}

dangerous_functions <- c(
  eval = "Dynamic evaluation is forbidden in the deployed runtime.",
  parse = "Dynamic parsing is forbidden in the deployed runtime.",
  system = "Use a fixed-argument system2 call when a subprocess is unavoidable.",
  shell = "Shell execution is forbidden in the deployed runtime.",
  "download.file" = "Runtime network downloads are forbidden in production."
)

runtime_linters <- list(
  duplicate_argument_linter = lintr::duplicate_argument_linter(),
  missing_argument_linter = lintr::missing_argument_linter(),
  unreachable_code_linter = lintr::unreachable_code_linter(),
  T_and_F_symbol_linter = lintr::T_and_F_symbol_linter(),
  undesirable_function_linter = lintr::undesirable_function_linter(dangerous_functions)
)

lints <- lintr::lint_dir(
  runtime_root,
  linters = runtime_linters,
  exclusions = list("tests"),
  show_progress = FALSE
)

if (length(lints)) {
  print(lints)
  stop("Deployed R runtime lint failed with ", length(lints), " finding(s).")
}

message("Deployed R runtime lint passed for ", runtime_root, ".")
