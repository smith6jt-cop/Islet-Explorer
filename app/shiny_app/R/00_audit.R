# Logging and security-boundary helpers are deliberately dependency-free so
# Shiny can source them before package loading in 00_globals.R.

.ISLET_LOG_EVENT_LEVELS <- c(
  app_loaded = "INFO",
  sf_ready = "INFO",
  vitessce_ready = "INFO",
  vitessce_routes_ready = "INFO",
  python_configured = "INFO",
  release_data_loaded = "INFO",
  application_event_invalid = "WARN",
  input_validation_failed = "WARN",
  sf_unavailable = "WARN",
  httr2_unavailable = "WARN",
  anndata_unavailable = "WARN",
  vitessce_unavailable = "WARN",
  vitessce_version_unsupported = "WARN",
  vitessce_check_failed = "WARN",
  vitessce_routes_failed = "WARN",
  app_load_flag_rejected = "WARN",
  python_anndata_unavailable = "WARN",
  python_unavailable = "WARN",
  phenotype_rules_unavailable = "WARN",
  viewer_channels_unavailable = "WARN",
  paga_unavailable = "WARN",
  scanpy_dpt_failed = "WARN",
  barlow_fallback_too_large = "WARN",
  segmentation_asset_missing = "WARN",
  segmentation_annotations_missing = "WARN",
  segmentation_annotations_schema_invalid = "WARN",
  islet_spatial_lookup_missing = "WARN",
  islet_spatial_lookup_schema_invalid = "WARN",
  runtime_profile_rejected = "ERROR",
  release_data_validation_failed = "ERROR",
  drilldown_render_failed = "ERROR",
  segmentation_render_failed = "ERROR",
  trajectory_data_load_failed = "ERROR",
  trajectory_feature_extract_failed = "ERROR",
  drilldown_cells_load_failed = "ERROR",
  barlow_pseudotime_failed = "ERROR",
  paga_module_import_failed = "ERROR",
  paga_compute_failed = "ERROR",
  segmentation_query_failed = "ERROR",
  segmentation_asset_read_failed = "ERROR",
  segmentation_bbox_failed = "ERROR",
  segmentation_intersection_failed = "ERROR",
  segmentation_annotations_load_failed = "ERROR",
  islet_spatial_lookup_load_failed = "ERROR",
  spatial_tissue_load_failed = "ERROR",
  spatial_backend_failed = "ERROR",
  spatial_cell_distance_load_failed = "ERROR"
)

.ISLET_VALIDATION_ERROR_CODES <- c(
  "viewer_image_invalid",
  "spatial_donor_invalid",
  "plot_mode_invalid",
  "plot_groups_invalid",
  "plot_selection_invalid",
  "plot_metric_invalid",
  "plot_min_cells_invalid",
  "trajectory_mode_invalid",
  "trajectory_feature_invalid",
  "trajectory_ratio_invalid",
  "trajectory_region_invalid",
  "islet_key_invalid",
  "drilldown_color_invalid"
)

islet_flag_enabled <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    tolower(trimws(value)) %in% c("1", "true", "yes", "on")
}

# Emit exactly one JSON line. The only arguments are allowlisted codes: raw
# messages, exceptions, paths, browser values and arbitrary context fields are
# intentionally unsupported. Every interpolated string is fixed or validated,
# so serialization is safe before any JSON/package dependency has loaded.
islet_log_event <- function(event_code, error_code = NULL) {
  if (!is.character(event_code) || length(event_code) != 1L ||
      is.na(event_code) || !(event_code %in% names(.ISLET_LOG_EVENT_LEVELS))) {
    event_code <- "application_event_invalid"
  }
  level <- .ISLET_LOG_EVENT_LEVELS[[event_code]]
  result <- if (identical(event_code, "input_validation_failed")) {
    "rejected"
  } else if (identical(level, "INFO")) {
    "ok"
  } else {
    "failed"
  }
  error_field <- ""
  if (identical(event_code, "input_validation_failed")) {
    if (!is.character(error_code) || length(error_code) != 1L ||
        is.na(error_code) || !(error_code %in% .ISLET_VALIDATION_ERROR_CODES)) {
      error_code <- "invalid_input"
    }
    error_field <- sprintf(',"error_code":"%s"', error_code)
  }
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  message(sprintf(
    paste0('{"ts":"%s","level":"%s","event":"%s",',
           '"component":"shiny","result":"%s"%s}'),
    timestamp, level, event_code, result, error_field
  ))
  invisible(NULL)
}

# Preserve the validation helper's FALSE result for callers that reject input.
islet_log_validation_failure <- function(error_code) {
  islet_log_event("input_validation_failed", error_code)
  invisible(FALSE)
}

islet_valid_scalar_choice <- function(value, allowed, error_code,
                                      log_missing = FALSE) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value) && value %in% as.character(allowed)
  if (!valid && (log_missing || !is.null(value))) {
    islet_log_validation_failure(error_code)
  }
  valid
}

islet_valid_choice_vector <- function(value, allowed, error_code,
                                      allow_empty = FALSE) {
  valid <- is.character(value) && !anyNA(value) &&
    (allow_empty || length(value) > 0L) && all(value %in% as.character(allowed))
  if (!valid && !is.null(value)) islet_log_validation_failure(error_code)
  valid
}

islet_valid_bounded_number <- function(value, lower, upper, error_code) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= lower && value <= upper
  if (!valid && !is.null(value)) islet_log_validation_failure(error_code)
  valid
}

islet_valid_asset_basename <- function(value, allowed) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value) && !grepl("[/\\\\]", value) &&
    grepl("[.]ome[.]tiff?$", value, ignore.case = TRUE) &&
    value %in% as.character(allowed)
  if (!valid && !is.null(value)) {
    islet_log_validation_failure("viewer_image_invalid")
  }
  valid
}

# Defense in depth for values that become a single filesystem path segment
# after first being checked against an authoritative application lookup.
islet_safe_path_segment <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", value) &&
    !(value %in% c(".", ".."))
}
