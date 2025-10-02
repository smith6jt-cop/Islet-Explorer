#!/usr/bin/env Rscript
# R helper functions to load Python pickle files created by preprocess_geojson.py
# Uses reticulate package to read pickle format

# Load pickle file using reticulate (or fall back to system call)
load_geojson_cache <- function(pkl_file) {
  if (!file.exists(pkl_file)) {
    stop(paste("Cache file not found:", pkl_file))
  }
  
  # Try using reticulate if available
  if (requireNamespace("reticulate", quietly = TRUE)) {
    tryCatch({
      pickle <- reticulate::import("pickle")
      con <- file(pkl_file, "rb")
      on.exit(close(con))
      data <- pickle$load(con)
      return(data)
    }, error = function(e) {
      # Fall through to JSON approach
    })
  }
  
  # Fallback: use Python to convert pickle to JSON
  temp_json <- tempfile(fileext = ".json")
  on.exit(unlink(temp_json), add = TRUE)
  
  # Create a temporary Python script
  temp_py <- tempfile(fileext = ".py")
  on.exit(unlink(temp_py), add = TRUE)
  
  writeLines(sprintf("
import pickle
import json
with open('%s', 'rb') as f:
    data = pickle.load(f)
with open('%s', 'w') as f:
    json.dump(data, f)
", pkl_file, temp_json), temp_py)
  
  system2("python3", args = temp_py)
  
  if (!file.exists(temp_json) || file.size(temp_json) == 0) {
    stop("Failed to load pickle file. Ensure python3 is available.")
  }
  
  data <- jsonlite::fromJSON(temp_json, simplifyVector = FALSE)
  return(data)
}

# Query features within bounding box
query_geojson_bbox <- function(data, xmin, ymin, xmax, ymax, max_features = 5000) {
  # Find intersecting features using spatial index
  spatial_idx <- data$spatial_index
  
  intersects <- vapply(spatial_idx, function(bbox) {
    bbox[[3]] >= xmin &&  # xmax >= xmin
    bbox[[1]] <= xmax &&  # xmin <= xmax
    bbox[[4]] >= ymin &&  # ymax >= ymin
    bbox[[2]] <= ymax     # ymin <= ymax
  }, logical(1))
  
  matching_ids <- which(intersects)
  
  # Limit to max_features
  if (length(matching_ids) > max_features) {
    matching_ids <- sample(matching_ids, max_features)
  }
  
  # Extract matching features
  features <- lapply(matching_ids, function(i) {
    f <- data$features[[i]]
    list(
      type = "Feature",
      id = f$id,
      geometry = list(
        type = "Polygon",
        coordinates = f$geometry
      ),
      properties = f$properties
    )
  })
  
  list(
    type = "FeatureCollection",
    features = features,
    n_total = data$n_features,
    n_returned = length(features)
  )
}

# Test function
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0 && args[1] == "test") {
    pkl_file <- if (length(args) >= 2) args[2] else "data/gson_cache/0112_simplified.pkl"
    cat("Testing pickle load from:", pkl_file, "\n")
    data <- load_geojson_cache(pkl_file)
    cat("Loaded", data$n_features, "features\n")
    cat("Global bbox:", paste(names(data$global_bbox), data$global_bbox, collapse = ", "), "\n")
  }
}
