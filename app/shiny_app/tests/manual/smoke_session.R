# A REAL smoke test: open a session, click every tab, and fail on any server
# error. HTTP 200 only proves the static page was served -- the R session can
# die on connect and the page still returns 200, which is exactly how a broken
# spatial_server signature reached the user.
suppressMessages(library(chromote))
b <- ChromoteSession$new(width=1600, height=1000)
b$Page$navigate("http://127.0.0.1:7780/"); Sys.sleep(30)
ev <- function(x) b$Runtime$evaluate(x, returnByValue=TRUE)$result$value
bad <- character(0)
alive <- function(where) {
  d <- ev("document.querySelectorAll('#shiny-disconnected-overlay').length")
  e <- ev("document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)').length")
  if (isTRUE(d > 0)) bad <<- c(bad, paste(where, "DISCONNECTED"))
  if (isTRUE(e > 0)) bad <<- c(bad, paste(where, "-", e, "output errors"))
  cat(sprintf("  %-14s disconnected=%s output-errors=%s\n", where, d, e))
}
alive("plot")
keys <- ev("Array.from(document.querySelectorAll('.ie-nav input')).map(function(i){return i.value}).join(',')")
for (i in seq_len(length(strsplit(keys, ",")[[1]])) - 1) {
  ev(sprintf("document.querySelectorAll('.ie-nav input')[%d].click()", i)); Sys.sleep(7)
  alive(ev("document.querySelector('.ie-nav input:checked').value"))
}
Sys.sleep(1); b$close()
cat(if (length(bad)) paste("SMOKE_FAIL:", paste(bad, collapse="; ")) else "SMOKE_PASS", "\n")
