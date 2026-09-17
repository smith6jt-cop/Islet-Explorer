# Runtime fitting was retired. Frozen schema-1 data roll back with their frozen
# application image; this version requires the offline trajectory artifacts.
compute_barlow_pseudotime <- function(min_cells = 1, max_cells = NA) {
  d <- trajectory_read("coordinates.csv")
  if (is.null(d)) return(NULL)
  d$pseudotime <- d$pseudotime_raw
  d
}
