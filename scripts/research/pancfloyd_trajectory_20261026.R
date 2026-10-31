#!/usr/bin/env Rscript
# Export the current app's descriptive LOESS and point-jitter calculations.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1]); out <- normalizePath(args[2])
source(file.path(repo, "app/shiny_app/R/plot_figures.R"))
d <- read.csv(file.path(out, "data/trajectory_points.csv"), check.names = FALSE)
points <- list(); curves <- list()
for (feature in unique(d$feature)) {
  z <- d[d$feature == feature, , drop = FALSE]
  width <- diff(range(z$pt_rank)) * 0.012
  beta_fraction <- identical(feature, "Beta endocrine fraction")
  height <- if (beta_fraction) 0 else diff(range(z$value)) * 0.08
  j <- plot_jitter(z$pt_rank, z$value, width, height, seed = 42)
  z$x_plot <- j$x; z$y_plot <- j$y
  points[[feature]] <- z
  if (beta_fraction) next
  for (status in c("ND", "Aab+", "T1D")) {
    s <- z[z$donor_status == status, , drop = FALSE]
    cr <- plot_loess_curve(s$pt_rank, s$value, weights = log1p(s$total_cells_core),
                           span = 0.75, se = TRUE, n = 80)
    if (is.null(cr)) next
    cr$feature <- feature; cr$donor_status <- status
    cr$lower <- cr$y - cr$ci; cr$upper <- cr$y + cr$ci
    curves[[paste(feature, status)]] <- cr
  }
}
write.csv(do.call(rbind, points), file.path(out, "data/trajectory_plotted_points.csv"), row.names = FALSE)
write.csv(do.call(rbind, curves), file.path(out, "data/trajectory_loess_curves.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_session_info.txt"))
