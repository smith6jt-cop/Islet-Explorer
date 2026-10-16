# Safe wrapper around dplyr::left_join to handle non-data.frame inputs
safe_left_join <- function(x, y, by, context = "join") {
  # Ensure inputs are data.frames
  if (!inherits(x, "data.frame")) {
    x <- tryCatch(as.data.frame(x), error = function(e) {
      stop(sprintf("[safe_left_join] %s: Cannot coerce LHS to data.frame (class: %s)",
                   context, paste(class(x), collapse = "/")))
    })
  }
  if (!inherits(y, "data.frame")) {
    y <- tryCatch(as.data.frame(y), error = function(e) {
      stop(sprintf("[safe_left_join] %s: Cannot coerce RHS to data.frame (class: %s)",
                   context, paste(class(y), collapse = "/")))
    })
  }

  # Perform join with error context
  tryCatch(
    dplyr::left_join(x, y, by = by),
    error = function(e) {
      stop(sprintf("[safe_left_join] %s: %s", context, e$message))
    }
  )
}

# Endocrine object names: Islet_N, Cluster_N, Single_N (the uniform workflow's three object classes).
ISLET_KEY_PATTERN <- "(Islet|Cluster|Single)_\\d+"

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)

  if (!"islet_key" %in% names(df)) df$islet_key <- NA_character_

  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, ISLET_KEY_PATTERN)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
  }

  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, ISLET_KEY_PATTERN)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name)
  }

  # A bare number names an islet only in the frozen release, which had no other object class.
  if ("islet_id" %in% names(df)) {
    id_str <- suppressWarnings(as.character(df$islet_id))
    key_from_id <- stringr::str_extract(id_str, ISLET_KEY_PATTERN)
    id_digits <- ifelse(grepl("^[0-9]+(\\.0+)?$", id_str), sub("\\.0+$", "", id_str), NA_character_)
    fallback_id <- ifelse(!is.na(id_digits), paste0("Islet_", id_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_id, fallback_id)
  }

  df
}

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
}
