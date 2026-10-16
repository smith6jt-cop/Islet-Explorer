test_that("viewer waits for session origin and preserves the iframe across image changes", {
  skip_if_not_installed("shiny")
  image_root <- tempfile("viewer-readiness-")
  dir.create(image_root)
  on.exit(unlink(image_root, recursive = TRUE), add = TRUE)
  file.create(file.path(image_root, c("first.ome.tiff", "second.ome.tiff", "T1D_0078.ome.tiff")))

  env <- new.env(parent = asNamespace("shiny"))
  env$resolve_avivator_base <- function() "avivator/index.html"
  env$local_images_root <- image_root
  env$islet_valid_asset_basename <- function(value, choices) value %in% choices
  env$default_image_url <- NULL
  env$channel_names_vec <- character()
  env$build_channel_config_b64 <- function(...) NULL
  env$VIEWER_DEBUG_ENABLED <- FALSE
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        as.character(expr[[2]]) %in% c("ISLET_EXCLUDED_DONORS", "islet_canonical_donor_id")) {
      eval(expr, envir = env)
    }
  }
  env$build_app_absolute_url <- function(session, rel_path) {
    paste0(session$clientData$url_pathname, rel_path)
  }
  env$build_public_http_url <- function(session, rel_path) {
    paste0(session$clientData$url_protocol, "//", session$clientData$url_hostname,
           ":", session$clientData$url_port, env$build_app_absolute_url(session, rel_path))
  }
  sys.source(.find_source("R/mod_viewer_server.R"), envir = env)
  mock <- shiny::MockShinySession$new()
  mock$clientData <- shiny::reactiveValues(
    url_protocol = NULL, url_hostname = NULL, url_port = NULL, url_pathname = NULL
  )
  shiny::testServer(env$viewer_server,
    args = list(forced_image = shiny::reactiveVal(NULL),
                current_tab = shiny::reactiveVal("Viewer")),
    session = mock, {
      session$flushReact()
      # Donor 78 was removed from the app; its image is never offered.
      expect_setequal(available_images_r(), c("first.ome.tiff", "second.ome.tiff"))
      expect_false(has_image_stable())
      expect_null(viewer_info()$iframe_src)
      expect_false(grepl("<iframe", output$viewer_frame$html, fixed = TRUE))

      session$clientData$url_protocol <- "http:"
      session$clientData$url_hostname <- "127.0.0.1"
      session$clientData$url_port <- "7778"
      session$flushReact()
      expect_false(has_image_stable())
      session$clientData$url_pathname <- "/islet-explorer/"
      session$flushReact()
      expect_true(has_image_stable())
      expect_true(viewer_info()$ok)
      expect_match(viewer_info()$iframe_src,
                   "http://127.0.0.1:7778/islet-explorer/images/first.ome.tiff", fixed = TRUE)
      initial_frame <- output$viewer_frame$html
      expect_true(grepl("<iframe", initial_frame, fixed = TRUE))

      session$setInputs(selected_image = "second.ome.tiff")
      expect_match(viewer_info()$iframe_src, "second.ome.tiff", fixed = TRUE)
      expect_identical(output$viewer_frame$html, initial_frame)
      session$clientData$url_hostname <- "later-host.invalid"
      session$flushReact()
      expect_match(viewer_info()$iframe_src, "127.0.0.1:7778", fixed = TRUE)
      expect_identical(output$viewer_frame$html, initial_frame)
    })
})
