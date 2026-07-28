#!/usr/bin/env Rscript
# Every `output$x <-` in a module server must have a matching `ns("x")` in its
# UI file. Nothing else catches this: parse() is happy, the testthat suite is
# happy, and the session smoke test is happy, because an unmounted output simply
# never renders. It is how the infiltration card's plot, table, note and
# download button vanished behind a section heading with nothing under it.
#
#   cd app/shiny_app && Rscript tests/manual/mount_audit.R
mods <- c("plot", "trajectory", "statistics", "spatial", "inflammation", "viewer")
bad <- 0
for (m in mods) {
  sv <- sprintf("R/mod_%s_server.R", m); ui <- sprintf("R/mod_%s_ui.R", m)
  if (!file.exists(sv) || !file.exists(ui)) next
  s <- readLines(sv); u <- paste(readLines(ui), collapse = "\n")
  outs <- unique(sub(".*output\\$([A-Za-z0-9_.]+).*", "\\1",
                     grep("output\\$[A-Za-z0-9_.]+ *<-", s, value = TRUE)))
  miss <- outs[!vapply(outs, function(o) grepl(sprintf('ns("%s")', o), u, fixed = TRUE), TRUE)]
  if (length(miss)) {
    cat(sprintf("  %-14s UNMOUNTED: %s\n", m, paste(miss, collapse = ", ")))
    bad <- bad + length(miss)
  }
}
if (bad == 0) cat("  all server outputs are mounted\n") else quit(status = 1)
