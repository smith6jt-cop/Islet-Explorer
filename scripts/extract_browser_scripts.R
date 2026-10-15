#!/usr/bin/env Rscript
# Parse literal browser programs from R syntax without sourcing/evaluating R.
# This is a syntax/SAST input collector, not cross-language taint analysis.
browser_literal_candidates <- function(parsed) {
  found <- list()
  visit <- function(node, formatted = FALSE) {
    if (is.character(node)) {
      for (value in node) {
        if (grepl("function\\s*[(]|function\\s+[[:alnum:]_$]+\\s*[(]|\\b(document|window|Shiny)\\s*[.$(]", value, perl = TRUE)) {
          found[[length(found) + 1L]] <<- list(value = value, formatted = formatted)
        }
      }
    } else if (is.call(node) || is.expression(node) || is.pairlist(node)) {
      children <- as.list(node)
      is_sprintf <- is.call(node) && identical(node[[1L]], as.name("sprintf"))
      for (i in seq_along(children)) {
        if (!identical(children[[i]], quote(expr = ))) visit(children[[i]], is_sprintf && i == 2L)
      }
    }
  }
  visit(parsed)
  found
}

browser_prepare_literal <- function(value, formatted) {
  if (formatted) {
    # Keep literal percent signs separate while replacing format arguments.
    value <- gsub("%%", "ISLET_LITERAL_PERCENT_SENTINEL", value, fixed = TRUE)
    value <- gsub("%([0-9]+\\$)?[-+ #0]*[0-9]*(\\.[0-9]+)?s", "ISLET_DYNAMIC_VALUE", value, perl = TRUE)
    value <- gsub("%([0-9]+\\$)?[-+ #0]*[0-9]*(\\.[0-9]+)?[diouxXeEfFgGaA]", "0", value, perl = TRUE)
    value <- gsub("ISLET_LITERAL_PERCENT_SENTINEL", "%", value, fixed = TRUE)
  }
  wrapped <- grepl("^\\s*(async\\s+)?function\\s*[(]", value, perl = TRUE)
  if (wrapped) value <- paste0("(\n", value, "\n);\n")
  list(code = value, anonymous_function_wrapped = wrapped, sprintf_placeholders_replaced = formatted)
}

browser_extract_file <- function(path, origin = path) {
  parsed <- parse(file = path, keep.source = TRUE) # A syntax error must abort.
  candidates <- browser_literal_candidates(parsed)
  if (!length(candidates)) return(list())
  parse_data <- utils::getParseData(parsed, includeText = TRUE)
  tokens <- parse_data
  tokens <- tokens[tokens$terminal & tokens$token == "STR_CONST", , drop = FALSE]
  tokens <- tokens[order(tokens$line1, tokens$col1), , drop = FALSE]
  records <- list()
  for (i in seq_len(nrow(tokens))) {
    literal <- parse(text = utils::getParseText(parse_data, tokens$id[i]), keep.source = FALSE)[[1L]]
    if (!is.character(literal) || length(literal) != 1L) next
    matched <- which(vapply(candidates, function(x) identical(x$value, literal), logical(1)))
    if (!length(matched)) next
    formatted <- any(vapply(candidates[matched], `[[`, logical(1), "formatted"))
    records[[length(records) + 1L]] <- c(
      list(origin = origin, line_start = tokens$line1[i], line_end = tokens$line2[i]),
      browser_prepare_literal(literal, formatted)
    )
  }
  records
}

browser_extract_tree <- function(source_root, output_dir) {
  if (dir.exists(output_dir)) stop("Extraction output already exists")
  files <- sort(c(file.path(source_root, "app.R"), list.files(file.path(source_root, "R"), "[.]R$", full.names = TRUE)))
  if (!length(files) || any(!file.exists(files))) stop("Application R sources missing")
  records <- list()
  for (path in files) {
    origin <- substring(path, nchar(source_root) + 2L)
    records <- c(records, browser_extract_file(path, origin))
  }
  if (!length(records)) stop("No browser literals found in application R sources")
  dir.create(output_dir, recursive = TRUE)
  for (i in seq_along(records)) {
    record <- records[[i]]
    name <- sprintf("%03d-%s-L%d.js", i, gsub("[^[:alnum:]_.-]", "_", record$origin), record$line_start)
    writeLines(c(paste0("// R literal origin: ", record$origin, ":", record$line_start), record$code), file.path(output_dir, name), useBytes = TRUE)
    records[[i]]$code <- NULL
    records[[i]]$target <- name
  }
  report <- list(schema_version = 1L, r_files_parsed = length(files), snippets = records,
                 scope = "Literal programs containing function/document/window/Shiny syntax; R is parsed, never evaluated.",
                 limitations = c("Runtime R substitutions use synthetic values; cross-language taint/data flow is not modeled.",
                                 "Fragments assembled dynamically, package-generated scripts and third-party bundles require separate review."))
  jsonlite::write_json(report, file.path(output_dir, "origins.json"), pretty = TRUE, auto_unbox = TRUE)
  invisible(report)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: Rscript scripts/extract_browser_scripts.R APP_ROOT OUTPUT_DIR")
  report <- browser_extract_tree(sub("/+$", "", args[1L]), args[2L])
  message("Parsed ", report$r_files_parsed, " R files; extracted ", length(report$snippets), " literal browser programs")
}
