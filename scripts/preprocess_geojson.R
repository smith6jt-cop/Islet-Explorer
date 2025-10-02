#!/usr/bin/env Rscript
# Preprocess massive GeoJSON files for efficient spatial queries
# Creates simplified, spatially-indexed versions for the viewer

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
})

#' Extract minimal geometry from GeoJSON (strip heavy properties)
#' @param geojson_gz_path Path to compressed .geojson.gz file
#' @param output_rds Path to save simplified RDS output
simplify_geojson_for_viewer <- function(geojson_gz_path, output_rds) {
  cat("[PREPROCESS] Processing:", basename(geojson_gz_path), "\n")
  
  # Stream parse to avoid loading entire file into memory
  con <- gzfile(geojson_gz_path, "rt")
  on.exit(close(con))
  
  # Read in chunks and extract only essential data
  features <- list()
  chunk_size <- 1000
  
  tryCatch({
    # Use stream_in for large files (avoids 2GB string limit)
    cat("[PREPROCESS] Streaming GeoJSON data...\n")
    
    # For compressed GeoJSON files, we need a different approach
    # Use jq or python to extract features if available, otherwise use ndjson streaming
    
    # Try using jsonlite's stream_in (requires ndjson format)
    # Since we have standard GeoJSON, we'll use a hybrid approach
    temp_ndjson <- tempfile(fileext = ".ndjson")
    on.exit(unlink(temp_ndjson), add = TRUE)
    
    # Convert GeoJSON to NDJSON using system tools if available
    # Otherwise fall back to line-by-line processing
    cat("[PREPROCESS] Converting to streaming format...\n")
    
    # Use zcat + sed to extract features one per line (works for most GeoJSON)
    # This is a workaround for R's 2GB string limit
    system2("bash", args = c("-c", sprintf(
      "zcat '%s' | python3 -c 'import sys, json; data=json.load(sys.stdin); [print(json.dumps(f)) for f in data[\"features\"]]' > '%s'",
      geojson_gz_path, temp_ndjson
    )))
    
    if (!file.exists(temp_ndjson) || file.size(temp_ndjson) == 0) {
      stop("Failed to convert GeoJSON to streaming format. Ensure python3 is available.")
    }
    
    cat("[PREPROCESS] Reading features from stream...\n")
    con_ndjson <- file(temp_ndjson, "r")
    on.exit(close(con_ndjson), add = TRUE)
    
    # Read features in chunks
    features_list <- list()
    i <- 0
    while (length(line <- readLines(con_ndjson, n = 1, warn = FALSE)) > 0) {
      i <- i + 1
      if (i %% 5000 == 0) {
        cat(sprintf("[PREPROCESS] Read %d features...\n", i))
      }
      tryCatch({
        feature <- jsonlite::fromJSON(line)
        features_list[[i]] <- feature
      }, error = function(e) {
        cat("[PREPROCESS] Warning: skipped malformed feature at line", i, "\n")
      })
    }
    
    cat("[PREPROCESS] Read", length(features_list), "features total\n")
    
    # Create structure similar to full_data
    full_data <- list(
      type = "FeatureCollection",
      features = list(
        geometry = lapply(features_list, function(f) f$geometry),
        properties = lapply(features_list, function(f) f$properties)
      )
    )
    
    if (is.null(full_data$features) || length(full_data$features$geometry) == 0) {
      stop("No features found in GeoJSON")
    }
    
    n_features <- length(full_data$features$geometry)
    cat("[PREPROCESS] Found", n_features, "features\n")
    
    # Extract minimal data: geometry + essential properties only
    cat("[PREPROCESS] Simplifying features...\n")
    simplified <- list()
    
    for (i in seq_len(n_features)) {
      if (i %% 5000 == 0) {
        cat(sprintf("[PREPROCESS] Processed %d / %d features (%.1f%%)\n", 
                    i, n_features, 100 * i / n_features))
      }
      
      # Extract geometry
      geom_obj <- full_data$features$geometry[[i]]
      
      # Get coordinates based on structure
      if (!is.null(geom_obj$coordinates)) {
        geom <- geom_obj$coordinates
      } else {
        next  # Skip if no coordinates
      }
      
      # Get bounding box for spatial indexing
      if (is.list(geom) && length(geom) > 0) {
        coords <- geom[[1]]  # First ring for polygon
        if (is.matrix(coords) || is.data.frame(coords)) {
          coords <- as.matrix(coords)
          bbox <- c(
            xmin = min(coords[, 1], na.rm = TRUE),
            ymin = min(coords[, 2], na.rm = TRUE),
            xmax = max(coords[, 1], na.rm = TRUE),
            ymax = max(coords[, 2], na.rm = TRUE)
          )
          
          props <- full_data$features$properties[[i]]
          
          # Store minimal feature
          simplified[[i]] <- list(
            id = i,
            geometry = geom,
            bbox = bbox,
            # Keep only essential properties (cell type, etc.)
            properties = list(
              classification = props$classification %||% props$Classification %||% NA,
              name = props$name %||% props$Name %||% NA
            )
          )
        }
      }
    }
    
    cat("[PREPROCESS] Creating spatial index...\n")
    # Create spatial grid index for faster queries
    all_bboxes <- do.call(rbind, lapply(simplified, function(f) f$bbox))
    
    result <- list(
      features = simplified,
      n_features = length(simplified),
      global_bbox = c(
        xmin = min(all_bboxes[, "xmin"], na.rm = TRUE),
        ymin = min(all_bboxes[, "ymin"], na.rm = TRUE),
        xmax = max(all_bboxes[, "xmax"], na.rm = TRUE),
        ymax = max(all_bboxes[, "ymax"], na.rm = TRUE)
      ),
      spatial_index = all_bboxes,
      processed_date = Sys.time()
    )
    
    cat("[PREPROCESS] Saving to RDS...\n")
    saveRDS(result, output_rds, compress = "xz")
    
    # Report size savings
    original_size <- file.size(geojson_gz_path) / 1024^2
    new_size <- file.size(output_rds) / 1024^2
    cat(sprintf("[PREPROCESS] Complete! Size: %.1f MB -> %.1f MB (%.1f%% reduction)\n",
                original_size, new_size, 100 * (1 - new_size / original_size)))
    
    return(invisible(result))
    
  }, error = function(e) {
    cat("[PREPROCESS] ERROR:", conditionMessage(e), "\n")
    stop(e)
  })
}

#' Query features within a bounding box (for viewport queries)
#' @param rds_data Preprocessed RDS data
#' @param xmin,ymin,xmax,ymax Viewport bounding box
query_bbox <- function(rds_data, xmin, ymin, xmax, ymax, max_features = 5000) {
  # Find features that intersect viewport
  intersects <- (
    rds_data$spatial_index[, "xmax"] >= xmin &
    rds_data$spatial_index[, "xmin"] <= xmax &
    rds_data$spatial_index[, "ymax"] >= ymin &
    rds_data$spatial_index[, "ymin"] <= ymax
  )
  
  matching_ids <- which(intersects)
  
  # Limit to max_features for performance
  if (length(matching_ids) > max_features) {
    # Sample spatially distributed features
    matching_ids <- sample(matching_ids, max_features)
  }
  
  features <- rds_data$features[matching_ids]
  
  return(list(
    type = "FeatureCollection",
    features = lapply(features, function(f) {
      list(
        type = "Feature",
        id = f$id,
        geometry = list(
          type = "Polygon",
          coordinates = f$geometry
        ),
        properties = f$properties
      )
    }),
    n_total = rds_data$n_features,
    n_returned = length(features)
  ))
}

# CLI interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript preprocess_geojson.R <input.geojson.gz> [output.rds]\n")
    cat("       Rscript preprocess_geojson.R --batch <gson_dir> <output_dir>\n")
    quit(status = 1)
  }
  
  if (args[1] == "--batch") {
    # Batch process all files in directory
    gson_dir <- args[2]
    output_dir <- args[3] %||% file.path(dirname(gson_dir), "gson_cache")
    
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    geojson_files <- list.files(gson_dir, pattern = "\\.geojson\\.gz$", full.names = TRUE)
    cat("[BATCH] Found", length(geojson_files), "files to process\n")
    
    for (gf in geojson_files) {
      base_name <- sub("\\.geojson\\.gz$", "", basename(gf))
      output_rds <- file.path(output_dir, paste0(base_name, "_simplified.rds"))
      
      if (file.exists(output_rds)) {
        cat("[BATCH] Skipping (already exists):", base_name, "\n")
        next
      }
      
      simplify_geojson_for_viewer(gf, output_rds)
    }
    
    cat("[BATCH] Complete!\n")
    
  } else {
    # Single file processing
    input_file <- args[1]
    output_file <- if (length(args) >= 2) args[2] else sub("\\.geojson\\.gz$", "_simplified.rds", input_file)
    
    simplify_geojson_for_viewer(input_file, output_file)
  }
}
