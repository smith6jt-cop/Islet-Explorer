#!/usr/bin/env Rscript
# Delegate to the single outcome-blind offline implementation. Both roots are
# mandatory, with islet identity decisions bound to source geometry hashes.
args <- commandArgs(trailingOnly = TRUE)
if (!all(c("--input-root", "--output-root", "--islet-review") %in% args)) {
  stop("Supply --input-root CANDIDATE --output-root CANDIDATE/trajectory --islet-review DECISIONS.csv")
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
python <- Sys.getenv("RETICULATE_PYTHON", unset = Sys.which("python3"))
status <- system2(python, c(shQuote(file.path(script_dir, "research", "trajectory.py")), vapply(args, shQuote, character(1))))
quit(status = status)
