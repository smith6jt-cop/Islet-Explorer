# Every app source file must parse. Nothing else in the suite reads app.R, so on 2026-10-19 a double-quoted word
# inside the CSS string in app.R passed every test while shiny-server refused to start the app.
test_that("app.R and every R/ file parse", {
  files <- c(.find_source("app.R"), list.files(dirname(.find_source("R/00_globals.R")), pattern = "\\.R$",
                                                full.names = TRUE))
  for (f in files) expect_silent(parse(f, keep.source = FALSE))
})
