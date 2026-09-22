#!/usr/bin/env Rscript
# Build only missing geometry indexes from the explicit candidate's own exports.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: index_candidate_geometry.R CANDIDATE_ROOT [CASE_ID ...]")
root <- normalizePath(args[1], mustWork = TRUE)
if (!file.exists(file.path(root, "candidate.json"))) stop("Explicit candidate identity required")
cohort <- read.csv(file.path(root, "cohort.csv"), colClasses = "character")
cases <- if (length(args) > 1L) args[-1] else cohort$case_id
if (!all(cases %in% cohort$case_id)) stop("Case outside candidate cohort")
dir.create(file.path(root, "geometry"), showWarnings = FALSE)
for (cid in cases) {
  key <- sprintf("%04d", as.integer(cid))
  target <- file.path(root, "geometry", paste0(key, ".gpkg"))
  if (file.exists(target)) stop("Refusing to overwrite geometry index: ", target)
  source <- file.path(root, "gson", paste0(key, ".geojson.gz"))
  g <- sf::st_read(paste0("/vsigzip/", source), quiet = TRUE, options = "ARRAY_AS_STRING=YES")
  sf::st_crs(g) <- NA
  sf::st_write(g, target, layer = "segmentation", layer_options = "SPATIAL_INDEX=YES", quiet = TRUE)
  message("Geometry indexed: ", key, " (", nrow(g), " objects)")
}
