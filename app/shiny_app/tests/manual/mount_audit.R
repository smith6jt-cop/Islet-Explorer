#!/usr/bin/env Rscript
# Every `output$x <-` in a module must have a matching `ns("x")` SOMEWHERE it can
# actually be reached from -- its own UI file, OR another server-side renderUI
# (several cards here are built that way: spatial's `neighborhood_cards` mounts
# ~17 outputs from inside the server, and statistics builds Section 4 the same
# way). Searching only the UI file reports those as unmounted, which buries the
# real hits.
#
# Nothing else catches this class: parse() is happy, the testthat suite is happy,
# and the session smoke test is happy, because an unmounted output simply never
# renders. It is how the infiltration card's plot, table, note and download
# button vanished behind a section heading with nothing under it.
#
#   cd app/shiny_app && Rscript tests/manual/mount_audit.R
mods <- c("plot", "trajectory", "statistics", "spatial", "inflammation", "viewer")
all_src <- paste(unlist(lapply(list.files("R", "[.]R$", full.names = TRUE), readLines)),
                 collapse = "\n")
app_src <- if (file.exists("app.R")) paste(readLines("app.R"), collapse = "\n") else ""
bad <- 0
for (m in mods) {
  sv <- sprintf("R/mod_%s_server.R", m)
  if (!file.exists(sv)) next
  s <- readLines(sv)
  outs <- unique(sub(".*output\\$([A-Za-z0-9_.]+).*", "\\1",
                     grep("output\\$[A-Za-z0-9_.]+ *<-", s, value = TRUE)))
  mounted <- function(o) {
    grepl(sprintf('ns("%s")', o), all_src, fixed = TRUE) ||
    grepl(sprintf('"%s"', o), app_src, fixed = TRUE)
  }
  miss <- outs[!vapply(outs, mounted, TRUE)]
  if (length(miss)) {
    cat(sprintf("  %-14s NOT MOUNTED ANYWHERE: %s\n", m, paste(miss, collapse = ", ")))
    bad <- bad + length(miss)
  }
}
if (bad == 0) cat("  every server output is reachable\n") else quit(status = 1)
