viewer_server <- function(id, forced_image, current_tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Avivator static build location is fixed for the session (file existence
    # check) — compute once, non-reactively.
    avivator_base <- resolve_avivator_base()
    iframe_id <- ns("avivator_viewer")
    wrap_id   <- ns("viewer_wrap")

    # Wait for the browser's actual origin and application path before making
    # an image URL. Latch this session metadata once complete so subsequent
    # clientData updates cannot reload the WebGL iframe.
    session_origin <- reactiveVal(NULL)
    observe({
      if (!is.null(isolate(session_origin()))) return()
      fields <- c("url_protocol", "url_hostname", "url_port", "url_pathname")
      origin <- setNames(lapply(fields, function(field) session$clientData[[field]]), fields)
      if (any(vapply(origin, is.null, logical(1)))) return()
      if (!origin$url_protocol %in% c("http:", "https:") ||
          !nzchar(origin$url_hostname) || !nzchar(origin$url_pathname)) return()
      session_origin(origin)
    })

    # Images available to pick from (www/local_images + LOCAL_IMAGE_ROOT).
    available_images_r <- reactive({
      imgs <- character(0)
      images_dir <- file.path("www", "local_images")
      if (dir.exists(images_dir)) {
        imgs <- c(imgs, list.files(images_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE))
      }
      local_root <- local_images_root
      if (!is.null(local_root) && nzchar(local_root) && dir.exists(local_root)) {
        imgs <- unique(c(imgs, list.files(local_root, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)))
      }
      imgs
    })

    # The basename currently being viewed (forced image wins, then the user's
    # dropdown choice, then the first available image).
    selected_basename_r <- reactive({
      ai <- available_images_r()
      fi <- forced_image()
      if (!is.null(fi) && islet_valid_asset_basename(fi, ai)) return(fi)
      si <- input$selected_image
      if (!is.null(si) && islet_valid_asset_basename(si, ai)) return(si)
      if (length(ai) > 0) return(ai[[1]])
      NULL
    })

    # Coarse "is there something to show" flag. Computed reactively, then
    # mirrored into a reactiveVal so it only invalidates dependents when the
    # boolean actually FLIPS. Switching between images keeps it TRUE, so
    # output$viewer_frame does NOT re-render — which would destroy + recreate
    # the WebGL <iframe> and cause the macOS WebRender flicker. Image-URL
    # changes are instead pushed into the existing iframe (see observeEvent).
    has_image_now <- reactive({
      sb <- selected_basename_r()
      info <- viewer_info()
      !is.null(avivator_base) && !is.null(sb) && nzchar(sb) &&
        isTRUE(info$ok) && !is.null(info$image_url) && nzchar(info$image_url)
    })
    has_image_stable <- reactiveVal(FALSE)
    observe({ has_image_stable(isTRUE(has_image_now())) })

    viewer_info <- reactive({
      info <- list(
        base = avivator_base, selection = NULL, ok = FALSE,
        iframe_src = NULL, image_url = NULL, image_public_url = NULL,
        image_rel = NULL, image_app_url = NULL, asset_diag = NULL, env = NULL
      )
      if (is.null(avivator_base)) return(info)
      origin <- session_origin()
      if (is.null(origin)) return(info)

      sel_basename <- selected_basename_r()

      rel_url <- NULL
      if (!is.null(sel_basename) && nzchar(sel_basename)) {
        www_candidate <- file.path("www", "local_images", sel_basename)
        resource_candidate <- NULL
        if (!is.null(local_images_root) && dir.exists(local_images_root)) {
          resource_candidate <- file.path(local_images_root, sel_basename)
        }
        if (file.exists(www_candidate)) {
          rel_url <- paste("local_images", sel_basename, sep = "/")
        } else if (!is.null(resource_candidate) && file.exists(resource_candidate)) {
          rel_url <- paste("images", sel_basename, sep = "/")
        } else {
          rel_url <- paste("local_images", sel_basename, sep = "/")
        }
        info$selection <- sel_basename
      }

      resolved_app_url <- NULL
      resolved_public_url <- NULL
      if (!is.null(rel_url) && nzchar(rel_url)) {
        info$image_rel <- rel_url
        origin_session <- list(clientData = origin)
        resolved_app_url    <- build_app_absolute_url(origin_session, rel_url)
        resolved_public_url <- build_public_http_url(origin_session, rel_url)
        info$image_app_url    <- resolved_app_url
        info$image_public_url <- resolved_public_url
      }

      sel_url <- NULL
      if (!is.null(resolved_public_url) && nzchar(resolved_public_url)) {
        sel_url <- resolved_public_url
      } else if (!is.null(resolved_app_url) && nzchar(resolved_app_url)) {
        sel_url <- resolved_app_url
      } else if (!is.null(default_image_url) && nzchar(default_image_url)) {
        sel_url <- default_image_url
      }

      if (is.null(sel_url) || !nzchar(sel_url)) return(info)
      params <- list(image_url = utils::URLencode(sel_url, reserved = FALSE))
      info$image_url <- resolved_app_url %||% sel_url

      ch_b64 <- tryCatch(build_channel_config_b64(channel_names_vec), error = function(e) NULL)
      if (!is.null(ch_b64) && nzchar(ch_b64)) {
        params[["channel_config"]] <- utils::URLencode(ch_b64, reserved = TRUE)
      }

      query <- NULL
      if (length(params)) {
        parts <- vapply(names(params), function(nm) sprintf("%s=%s", nm, params[[nm]]), character(1))
        query <- paste(parts, collapse = "&")
      }
      info$iframe_src <- if (!is.null(query)) paste0(avivator_base, "?", query) else avivator_base
      info$ok <- TRUE

      # Diagnostics only — detect_environment() and probe_viewer_asset() both
      # read session$clientData / do a blocking curl HEAD, so keep them off the
      # hot path unless explicitly debugging.
      if (VIEWER_DEBUG_ENABLED) {
        info$env <- isolate(detect_environment(session))
        info$asset_diag <- isolate(probe_viewer_asset(session, rel_url))
      }
      info
    })

    # ---- Image picker (mounts once; independent of the iframe) --------------
    output$image_selector <- renderUI({
      tab <- current_tab()
      if (is.null(tab) || tab != "Viewer") return(NULL)
      if (is.null(avivator_base)) return(NULL)  # missing-build message lives in viewer_frame

      available_images <- available_images_r()
      # isolate(): the selector reflects the current choice on first mount but
      # must not take a reactive dependency on it, or picking an image would
      # rebuild the dropdown.
      current_selection <- isolate(selected_basename_r())

      tagList(
        if (length(available_images) > 0) {
          fluidRow(
            column(12,
              selectInput(ns("selected_image"), "Select Image:",
                          choices = available_images, selected = current_selection, width = "100%")
            )
          )
        } else {
          tags$div(
            style = "padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; margin-bottom: 10px;",
            tags$strong("No images found."),
            " Place OME-TIFF files in ", tags$code("app/shiny_app/www/local_images/"),
            " or set ", tags$code("LOCAL_IMAGE_ROOT"), " environment variable."
          )
        },
        if (VIEWER_DEBUG_ENABLED) {
          vi <- isolate(viewer_info())
          tags$details(
            style = "margin: 10px 0;",
            tags$summary("Viewer debug info"),
            tags$pre(
              style = "max-height:200px; overflow:auto;",
              jsonlite::toJSON(list(
                env = vi$env, image_rel = vi$image_rel,
                image_app = vi$image_app_url %||% vi$image_url,
                image_public = vi$image_public_url,
                iframe_src = vi$iframe_src,
                asset_diag = vi$asset_diag
              ), auto_unbox = TRUE, pretty = TRUE)
            )
          )
        }
      )
    })

    # ---- Avivator iframe (mounts once when an image first becomes available) -
    output$viewer_frame <- renderUI({
      tab <- current_tab()
      if (is.null(tab) || tab != "Viewer") return(NULL)

      if (is.null(avivator_base)) {
        return(tagList(
          tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
          tags$div(style = "color:#666; font-size:90%;",
                   "To install: run scripts/install_avivator.sh (requires Node >= 18) or place a prebuilt bundle under shiny/www/avivator.")
        ))
      }

      # Depend on the FLIP-only boolean (mounts once); read the URL via isolate
      # so this renderUI does NOT re-fire — and recreate the iframe — when the
      # src changes. URL changes are pushed in place via the observeEvent below.
      ready <- has_image_stable()
      if (!isTRUE(ready)) {
        return(tags$div(
          style = "padding: 20px; text-align: center; color: #666;",
          "Select an image from the dropdown above to view it in AVIVATOR."
        ))
      }

      src0 <- isolate(viewer_info()$iframe_src)
      if (is.null(src0) || !nzchar(src0)) {
        return(tags$div(
          style = "padding: 20px; text-align: center; color: #666;",
          "Select an image from the dropdown above to view it in AVIVATOR."
        ))
      }

      tags$div(
        id = wrap_id,
        style = "position: relative;",
        tags$div(
          style = "padding: 6px 10px; background-color: #f6f8fa; color: #666; font-size: 90%; margin-bottom: 6px; border: 1px solid #e1e4e8; border-radius: 3px;",
          "Click inside the viewer to enable pan/zoom and wheel zoom."
        ),
        tags$iframe(
          id = iframe_id, src = src0, `data-relsrc` = src0,
          title = "Interactive pancreatic tissue image viewer",
          width = "100%", frameBorder = 0, allowfullscreen = NA,
          tabindex = "0", allow = "fullscreen",
          style = "width:100%; height:calc(100vh - 180px); border:0;"
        ),
        tags$script(HTML(sprintf("(function(){
  var wrap = document.getElementById('%s');
  if (!wrap) return;
  var iframe = wrap.querySelector('#%s');
  if (!iframe) return;
  function focusViewer(){ try { if (iframe.contentWindow) iframe.contentWindow.focus(); } catch (e) {} }
  iframe.addEventListener('load', focusViewer);
  wrap.addEventListener('mouseenter', focusViewer);
  wrap.addEventListener('wheel', function(){
    try {
      if (document.activeElement !== iframe) focusViewer();
    } catch (e) {}
  }, { passive: true });
  setTimeout(focusViewer, 200);
})();", wrap_id, iframe_id)))
      )
    })

    # ---- Push src changes into the existing iframe (no DOM recreation) -------
    observeEvent(viewer_info()$iframe_src, {
      src <- viewer_info()$iframe_src
      if (is.null(src) || !nzchar(src)) return()
      session$sendCustomMessage("islet_update_viewer_src",
                                list(id = iframe_id, src = src))
    }, ignoreNULL = TRUE)
  })
}
