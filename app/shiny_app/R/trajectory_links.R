# Stable islet identity shared by the trajectory scatter and both UMAPs.
trajectory_link_key <- function(case_id, islet_key) {
  if (length(case_id) != length(islet_key)) stop("Islet identity lengths differ")
  donor <- as.character(case_id)
  islet <- sub("^Islet_Islet_", "Islet_", as.character(islet_key))
  donor <- sub("^0+(?=[0-9])", "", donor, perl = TRUE)
  valid <- !is.na(donor) & !is.na(islet) &
    grepl("^[0-9]+$", donor) & grepl("^Islet_[0-9]+$", islet)
  out <- rep(NA_character_, length(donor))
  out[valid] <- paste(donor[valid], islet[valid], sep = "|")
  out
}

trajectory_link_plot <- function(p, ns) {
  # Build once so keys stay aligned with ggplotly's color/group trace splits.
  p <- plotly::plotly_build(p)
  # Plotly's R binding otherwise forwards every hover/afterplot event to Shiny.
  p$x$shinyEvents <- character(0)
  for (axis in c("xaxis", "yaxis")) {
    if (length(p$x$layout[[axis]]$range) == 2L) {
      p$x$layout[[axis]]$autorange <- FALSE
    }
  }
  for (i in seq_along(p$x$data)) {
    trace <- p$x$data[[i]]
    if (grepl("markers", trace$mode %||% "", fixed = TRUE) &&
        length(trace$key) == length(trace$x) && length(trace$key) > 0L) {
      trace$meta <- list(trajectory_role = "islet")
      p$x$data[[i]] <- trace
    }
  }
  # prepareWidget()/preRenderHook can rebuild the widget and restore the R
  # binding's default events. Disable them after the final build as well.
  p$preRenderHook <- function(widget) {
    widget <- plotly::plotly_build(widget)
    widget$x$shinyEvents <- character(0)
    widget
  }
  htmlwidgets::onRender(p, "function(el, x, config) {
    window.IsletTrajectoryLinks.attach(el, config);
  }", data = list(group = ns("linked_plots"), input = ns("linked_click"),
                  ready = ns("linked_ready"), message = ns("linked_state")))
}
