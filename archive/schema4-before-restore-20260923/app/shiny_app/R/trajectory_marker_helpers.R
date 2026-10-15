# Schema-4 heatmaps read precomputed typed-cell summaries, never whole-islet
# median intensities or trajectory feature columns. DuckDB pushes selection into
# the immutable Parquet scan instead of loading every population into R.
trajectory_read_marker_summaries <- function(population, markers, donor = "all", region = "both") {
  path <- data_path("trajectory", "marker_summaries.parquet")
  regions <- if (identical(region, "both")) c("core", "peri") else region
  if (file.exists(path)) {
    con <- islet_duckdb_connection()
    quoted <- function(x) paste(vapply(x, islet_sql_string, character(1), con = con), collapse = ",")
    clauses <- c(paste0("population = ", quoted(population)),
      paste0("marker IN (", if (length(markers)) quoted(markers) else "NULL", ")"),
      paste0("compartment IN (", quoted(regions), ")"))
    if (!identical(donor, "all")) clauses <- c(clauses, paste0("CAST(case_id AS VARCHAR) = ", quoted(donor)))
    result <- DBI::dbGetQuery(con, paste0("SELECT case_id,islet_id,compartment,population,marker,n_cells,mean_log1p,positive_count,positive_available FROM ",
      islet_parquet_relation(path, con), " WHERE ", paste(clauses, collapse = " AND ")))
  } else {
    result <- trajectory_read("marker_summaries.csv")
    if (is.null(result)) return(NULL)
    result <- result[result$population == population & result$marker %in% markers & result$compartment %in% regions &
      (identical(donor, "all") | as.character(result$case_id) == donor), , drop = FALSE]
  }
  result$case_id <- as.character(result$case_id)
  result$positive_available <- tolower(as.character(result$positive_available)) %in% c("true", "t", "1")
  result
}

trajectory_typed_marker_values <- function(coords, summaries, denominator, markers,
    measure = "positive_ratio", population = "All cells", region = "both") {
  if (is.null(coords) || is.null(summaries)) return(NULL)
  stopifnot(measure %in% c("positive_ratio", "typed_expression"), region %in% c("core", "peri", "both"))
  key <- c("case_id", "islet_id", "marker")
  grid <- coords[rep(seq_len(nrow(coords)), each = length(markers)), , drop = FALSE]
  grid$marker <- rep(markers, times = nrow(coords))
  grid$population <- rep(population, nrow(grid))
  grid$measurement <- rep(measure, nrow(grid))
  grid$numerator_compartment <- rep(region, nrow(grid))
  summaries <- summaries[summaries$population == population & summaries$marker %in% markers, , drop = FALSE]
  if (anyDuplicated(summaries[c(key, "compartment")])) stop("Duplicate typed marker summaries")
  lookup <- function(compartment) {
    s <- summaries[summaries$compartment == compartment, c(key, "n_cells", "mean_log1p", "positive_count", "positive_available"), drop = FALSE]
    indices <- match(do.call(paste, c(grid[key], sep = "\r")), do.call(paste, c(s[key], sep = "\r")))
    s[indices, , drop = FALSE]
  }
  a <- lookup(if (identical(region, "peri")) "peri" else "core")
  n_cells <- a$n_cells
  means <- a$mean_log1p
  positives <- a$positive_count
  available <- a$positive_available %in% TRUE
  if (identical(region, "both")) {
    b <- lookup("peri")
    n_cells <- a$n_cells + b$n_cells
    # An assessed compartment with zero cells contributes zero to the sum, but
    # an unassessed compartment remains NA rather than silently disappearing.
    sum_a <- ifelse(a$n_cells == 0, 0, a$n_cells * a$mean_log1p)
    sum_b <- ifelse(b$n_cells == 0, 0, b$n_cells * b$mean_log1p)
    means <- ifelse(n_cells > 0, (sum_a + sum_b) / n_cells, NA_real_)
    positives <- a$positive_count + b$positive_count
    available <- available & b$positive_available %in% TRUE
  }
  grid$n_cells <- n_cells
  grid$mean_log1p <- means
  grid$positive_count <- ifelse(available, positives, NA_real_)
  grid$positive_available <- available
  grid$denominator <- rep(NA_real_, nrow(grid))
  if (!is.null(denominator)) {
    keys <- intersect(c("case_id", "islet_id"), names(denominator))
    ix <- match(do.call(paste, c(grid[keys], sep = "\r")), do.call(paste, c(denominator[keys], sep = "\r")))
    validated <- if ("validated" %in% names(denominator)) denominator$validated[ix] %in% TRUE else
      if ("available" %in% names(denominator)) denominator$available[ix] %in% TRUE else is.finite(denominator$endocrine_count[ix])
    grid$denominator <- ifelse(validated, denominator$endocrine_count[ix], NA_real_)
  }
  valid_ratio <- is.finite(grid$denominator) & grid$denominator > 0 & is.finite(grid$positive_count)
  grid$raw_ratio <- ifelse(valid_ratio, grid$positive_count / grid$denominator, NA_real_)
  grid$value <- if (identical(measure, "positive_ratio")) log1p(grid$raw_ratio) else ifelse(n_cells > 0, means, NA_real_)
  grid$availability_reason <- if (identical(measure, "positive_ratio"))
    ifelse(!available, "marker_positivity_unavailable", ifelse(!is.finite(grid$denominator), "classifier_core_cell_census_unavailable",
      ifelse(grid$denominator == 0, "no_core_cells_in_classifier_area", "available"))) else
    ifelse(!is.finite(n_cells), "population_measurements_unavailable", ifelse(n_cells == 0, "no_cells_in_selected_population",
      ifelse(!is.finite(means), "marker_measurement_unavailable", "available")))
  grid
}

trajectory_heatmap_bin_values <- function(d, axis, markers, bins = 20L) {
  breaks <- seq(axis$limits[1], axis$limits[2], length.out = bins + 1L)
  # Complete marker/bin grid keeps unavailable markers visibly gray.
  out <- expand.grid(bin = seq_len(bins), marker = markers, stringsAsFactors = FALSE)
  out$value <- NA_real_; out$n_donors <- 0L; out$n_islets <- 0L
  if (!is.null(d) && nrow(d)) {
    d$bin <- cut(axis$transform(d[[axis$column]]), breaks = breaks, include.lowest = TRUE, labels = FALSE)
    available <- d[!is.na(d$bin) & is.finite(d$value), , drop = FALSE]
    if (nrow(available)) {
      donor <- aggregate(value ~ case_id + bin + marker, available, mean)
      means <- aggregate(value ~ bin + marker, donor, mean)
      donors <- aggregate(case_id ~ bin + marker, donor, length)
      islets <- aggregate(islet_id ~ bin + marker, available, length)
      ix <- match(paste(out$bin, out$marker), paste(means$bin, means$marker))
      out$value <- means$value[ix]
      ix <- match(paste(out$bin, out$marker), paste(donors$bin, donors$marker))
      out$n_donors <- ifelse(is.na(ix), 0L, donors$case_id[ix])
      ix <- match(paste(out$bin, out$marker), paste(islets$bin, islets$marker))
      out$n_islets <- ifelse(is.na(ix), 0L, islets$islet_id[ix])
    }
  }
  out$x <- (breaks[out$bin] + breaks[out$bin + 1L]) / 2
  out$width <- diff(breaks)[1]
  out
}
