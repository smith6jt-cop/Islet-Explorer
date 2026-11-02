# app.R's process warm-up draws one small PNG to build the font cache (islet_warm_png(), 00_globals.R). With the
# default margins plot.new() stopped on its 64 px device ("figure margins too large"), so every production
# warm-up failed before it reached the trajectory caches (v20 release container, 2026-11-01).

test_that("the warm-up PNG draw runs on its 64 px device and closes it", {
  env <- new.env(parent = globalenv())
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        identical(as.character(expr[[2]]), "islet_warm_png")) {
      eval(expr, envir = env)
    }
  }
  devices <- grDevices::dev.list()
  expect_true(env$islet_warm_png())
  expect_identical(grDevices::dev.list(), devices)
})
