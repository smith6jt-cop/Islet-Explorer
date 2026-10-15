#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: benchmark_candidate_runtime.R REPO ROOT OUTPUT_JSON")
repo <- normalizePath(args[1], mustWork = TRUE)
data_root <- normalizePath(args[2], mustWork = TRUE)
output <- args[3]
Sys.setenv(ISLET_DATA_SCHEMA = "3", ISLET_RUNTIME_MODE = "development")
data_path <- function(...) file.path(data_root, ...)
source(file.path(repo, "app/shiny_app/R/01_runtime_data_backend.R"))
source(file.path(repo, "app/shiny_app/R/mod_trajectory_server.R"))
measure <- function(expr) unname(system.time(force(expr))[["elapsed"]])
manifest_seconds <- measure(islet_validate_release(data_root))
names <- c("coordinates", "features", "counts", "endocrine_counts", "edges")
tables <- setNames(vector("list", length(names)), names)
cold <- measure(for (name in names) tables[[name]] <- trajectory_read(paste0(name, ".csv")))
warm <- replicate(20, measure(for (name in names) trajectory_read(paste0(name, ".csv"))))
d <- tables$coordinates
selected <- d[d$case_id == "78" & d$measured_cells >= 20, ][1, ]
single <- data_path("cells", paste0(selected$islet_id, ".csv"))
cells_seconds <- measure(cells <- read.csv(single))
pixel <- .3774; x <- selected$centroid_x_um / pixel; y <- selected$centroid_y_um / pixel; b <- 150 / pixel
filter <- sprintf("POLYGON ((%.17g %.17g, %.17g %.17g, %.17g %.17g, %.17g %.17g, %.17g %.17g))",
                  x-b,y-b,x+b,y-b,x+b,y+b,x-b,y+b,x-b,y-b)
gpkg <- data_path("geometry", "0078.gpkg")
geometry_seconds <- measure(crop <- sf::st_read(gpkg, layer = "segmentation", wkt_filter = filter, quiet = TRUE))
report <- list(manifest_core_verification_seconds = manifest_seconds,
  trajectory_tables_cold_seconds = cold, trajectory_tables_cached_median_seconds = median(warm),
  trajectory_tables_cached_max_seconds = max(warm), trajectory_tables_bytes_in_R = sum(vapply(tables, function(d) as.numeric(object.size(d)), numeric(1))),
  table_rows = lapply(tables, nrow), inspected_islet = selected$islet_id,
  single_islet_cells_seconds = cells_seconds, single_islet_cell_rows = nrow(cells),
  single_islet_cells_bytes_in_R = as.numeric(object.size(cells)),
  indexed_geometry_query_seconds = geometry_seconds, cropped_geometry_rows = nrow(crop),
  cropped_geometry_bytes_in_R = as.numeric(object.size(crop)), donor_geometry_file_bytes = file.info(gpkg)$size,
  R_version = R.version.string, note = "Single-worker microbenchmark; warm OS cache is possible. Browser timings and worker memory are separate measurements. No concurrent-user capacity claim.")
jsonlite::write_json(report, output, pretty = TRUE, auto_unbox = TRUE)
islet_disconnect_backend()
cat(jsonlite::toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
