# Shared explanations and whole-cohort counts for the Unclassified family.
# Counts come from raw tissue records, before any card/region/phenotype filter.
# The only database operation is a SELECT over the read-only Parquet files.
.unclassified_census_cache <- new.env(parent = emptyenv())

unclassified_census_sources <- function(root) {
  root <- normalizePath(root, mustWork = FALSE)
  tissue_dir <- file.path(root, "spatial", "tissue")
  tissue <- sort(list.files(tissue_dir, pattern = "[.]parquet$", recursive = TRUE,
                            full.names = TRUE))
  paths <- c(file.path(root, "manifest.json"), file.path(root, "donor_metadata.csv"), tissue)
  info <- file.info(paths)
  # File identities, including missing files, invalidate both successful and
  # unavailable results. The release itself is immutable in production.
  signature <- paste(root, paste(paths, info$size, as.numeric(info$mtime), sep = "|", collapse = "\n"),
                     sep = "\n")
  list(root = root, tissue = tissue, metadata = paths[[2]], signature = signature)
}

unclassified_census_read <- function(sources) {
  if (!length(sources$tissue) || !file.exists(sources$metadata)) {
    stop("Complete tissue Parquet and donor metadata are required.", call. = FALSE)
  }
  metadata <- utils::read.csv(sources$metadata, stringsAsFactors = FALSE, check.names = FALSE,
                              colClasses = "character")
  if (!all(c("case_id", "donor_type") %in% names(metadata))) {
    stop("Donor metadata lacks case_id or donor_type.", call. = FALSE)
  }
  ids <- islet_canonical_donor_id(metadata$case_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Donor metadata must have unique, nonempty canonical IDs.", call. = FALSE)
  }
  status_map <- c("ND" = "ND", "No Diabetes" = "ND", "Aab+" = "Aab+",
                  "Autoantibody Positive" = "Aab+", "T1D" = "T1D", "Type 1 Diabetes" = "T1D")
  status <- unname(status_map[trimws(metadata$donor_type)])

  con <- islet_duckdb_connection()
  files_sql <- paste(vapply(sources$tissue, function(path) {
    islet_sql_string(normalizePath(path, mustWork = TRUE), con)
  }, character(1)), collapse = ", ")
  # Materialize only donor/label aggregates, never the full cell table. Explicit
  # file names avoid selecting files outside the fingerprinted snapshot.
  counts <- DBI::dbGetQuery(con, paste0(
    "SELECT CAST(case_id AS VARCHAR) AS case_id, phenotype, COUNT(*) AS cells ",
    "FROM read_parquet([", files_sql, "], hive_partitioning = true) ",
    "GROUP BY case_id, phenotype"))
  if (!nrow(counts) || any(!is.finite(counts$cells)) || any(counts$cells < 0)) {
    stop("Tissue census is empty or has invalid counts.", call. = FALSE)
  }
  counts$case_id <- islet_canonical_donor_id(counts$case_id)
  donor_match <- match(counts$case_id, ids)
  if (anyNA(donor_match) || anyNA(status[donor_match])) {
    stop("Every tissue donor must have an ND, Aab+ or T1D metadata status.", call. = FALSE)
  }
  if (length(islet_excluded_donors_in(counts$case_id))) {
    stop("Tissue census contains an excluded donor.", call. = FALSE)
  }
  counts$donor_status <- status[donor_match]
  counts$subclass <- unclassified_collapse(counts$phenotype)
  keys <- unname(UNCLASSIFIED_KEYS)
  selected <- counts[counts$subclass %in% keys, , drop = FALSE]
  # In particular, literal raw Unclassified is the reassessed subclass, not the
  # parent aggregate. Conflict labels are not members of any of these four keys.
  totals <- matrix(0, nrow = length(keys), ncol = 3L,
                   dimnames = list(keys, c("ND", "Aab+", "T1D")))
  for (i in seq_len(nrow(selected))) {
    totals[selected$subclass[[i]], selected$donor_status[[i]]] <-
      totals[selected$subclass[[i]], selected$donor_status[[i]]] + selected$cells[[i]]
  }
  totals <- rbind(colSums(totals), totals)
  rownames(totals)[[1]] <- UNCLASSIFIED_PARENT
  totals <- cbind(totals, Total = rowSums(totals))
  table <- data.frame(key = rownames(totals), totals, row.names = NULL, check.names = FALSE)
  list(available = TRUE, table = table, raw_counts = selected,
       tissue_cells = sum(counts$cells), donors = length(unique(counts$case_id)),
       source_root = sources$root, tissue_partitions = length(sources$tissue),
       metadata_path = sources$metadata,
       scope = "All tissue cells in the loaded cohort; independent of every card selection.")
}

#' Cached, read-only cohort census used by every Unclassified help menu.
#' Missing/invalid input is unavailable, never an invented zero-count census.
unclassified_census <- function(root = data_root) {
  sources <- unclassified_census_sources(root)
  if (identical(.unclassified_census_cache$signature, sources$signature)) {
    return(.unclassified_census_cache$value)
  }
  value <- tryCatch(unclassified_census_read(sources), error = function(e) {
    if (exists("islet_log_event", mode = "function")) islet_log_event("unclassified_census_unavailable")
    list(available = FALSE, reason = conditionMessage(e), source_root = sources$root)
  })
  .unclassified_census_cache$signature <- sources$signature
  .unclassified_census_cache$value <- value
  value
}

#' Shared popover content. Safe to append to static or dynamically built card
#' help: no session/output IDs, no writes, and no per-card filtered totals.
unclassified_help_ui <- function(root = data_root) {
  census <- unclassified_census(root)
  table_ui <- if (isTRUE(census$available)) {
    headings <- c("Classification", "ND", "Aab+", "T1D", "Total")
    rows <- lapply(seq_len(nrow(census$table)), function(i) {
      row <- census$table[i, , drop = FALSE]
      shiny::tags$tr(
        shiny::tags$th(scope = "row", style = "text-align: left; font-weight: normal; padding: 4px;",
                       phenotype_text(row$key)),
        lapply(c("ND", "Aab+", "T1D", "Total"), function(status) {
          shiny::tags$td(style = "text-align: right; white-space: nowrap; padding: 4px;",
                         format(row[[status]], big.mark = ",", scientific = FALSE, trim = TRUE))
        }))
    })
    shiny::tags$div(style = "overflow-x: auto;",
      shiny::tags$table(class = "unclassified-census-table",
        style = "width: 100%; font-size: 12px; border-collapse: collapse;",
        shiny::tags$caption(style = "text-align: left; color: inherit;", "Global cell counts by donor status"),
        shiny::tags$thead(shiny::tags$tr(lapply(headings, function(x) {
          shiny::tags$th(scope = "col", style = "padding: 4px;", x)
        }))),
        shiny::tags$tbody(rows)))
  } else {
    shiny::tags$p("Global counts are unavailable: the complete tissue census and donor-status metadata could not be verified.")
  }
  shiny::tagList(
    shiny::tags$hr(),
    shiny::tags$p(shiny::tags$strong("Unclassified cells."),
      " This family preserves cells whose available evidence does not resolve a named phenotype; it is not one biological cell type. ",
      shiny::tags$strong("Unclassified (all)"),
      " is the sum of the four mutually exclusive subclasses below. Do not add that parent total to its subclasses."),
    shiny::tags$dl(
      shiny::tags$dt(phenotype_text(UNCLASSIFIED_KEYS[["tissue"]])),
      shiny::tags$dd("The available marker pattern was insufficient for a specific typing call. Raw source label: ",
        shiny::tags$code("Tissue"), ". This does not mean that no cell or no marker signal was present."),
      shiny::tags$dt(phenotype_text(UNCLASSIFIED_KEYS[["reassessed"]])),
      shiny::tags$dd("Reassessment did not support a resolved phenotype from the current marker and compartment evidence. ",
        "Some cells retain shared, endocrine-associated or structure-associated signal; others have no marker above background. ",
        "Endocrine-object membership alone is not a typing call. Raw source label: ", shiny::tags$code("Unclassified"),
        "; that literal source label is only this subclass, not the whole family."),
      shiny::tags$dt(phenotype_text(UNCLASSIFIED_KEYS[["subtype"]])),
      shiny::tags$dd("Broad T-cell or hematopoietic evidence is present, but a named lineage or subtype is unresolved. ",
        "Raw source labels: ", shiny::tags$code("T cell"), ", ", shiny::tags$code("Hematopoietic CD45+"),
        " and ", shiny::tags$code("Hematopoietic CD45-"),
        ". The CD45 sign records the source marker state; it does not supply a resolved subtype."),
      shiny::tags$dt(phenotype_text(UNCLASSIFIED_KEYS[["competing"]])),
      shiny::tags$dd("Signals meet criteria in at least three lineage or compartment groups. ",
        "This indicates competing marker support, not proof of a physical doublet. Raw source label: ",
        shiny::tags$code("Doublet suspected"), ".")),
    shiny::tags$p("The counts below cover all tissue cells across the loaded cohort, including core, peri-islet and background tissue. ",
      "They are global totals, not the current card's selected donors, regions, objects or phenotypes. ",
      "ND means no diabetes; Aab+ means autoantibody positive; T1D means type 1 diabetes. Conflict records are excluded from this family."),
    table_ui)
}
