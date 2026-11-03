# Stable islet identity shared by linked plots: the trajectory scatter and both UMAPs, and the Plot tab's size
# scatter and distribution plot (trajectory_link_plot() serves both tabs).
trajectory_link_key <- function(case_id, islet_key) {
  if (length(case_id) != length(islet_key)) stop("Islet identity lengths differ")
  donor <- as.character(case_id)
  islet <- sub("^Islet_Islet_", "Islet_", as.character(islet_key))
  donor <- sub("^0+(?=[0-9])", "", donor, perl = TRUE)
  valid <- !is.na(donor) & !is.na(islet) &
    grepl("^[0-9]+$", donor) & grepl("^(Islet|Cluster|Single)_[0-9]+$", islet)
  out <- rep(NA_character_, length(donor))
  out[valid] <- paste(donor[valid], islet[valid], sep = "|")
  out
}

#' Build a plotly figure and drop the build inputs plotly keeps beside the built traces (`x$attrs`, `x$visdat`,
#' `x$cur_data`). The widget sent them to the browser as well, which doubled a plot_ly() figure's payload
#' (2026-10-26: 311 of 657 KB on a 5,000-point trace; the Plot cards 884 and 951 KB instead of about 450 each).
#' The browser reads only data, layout and config, and later layout(), config() and plotly_build() calls still work
#' on the stripped figure. Scalars plotly recycled to one value per point (hoverinfo "text" on every point, 65 KB a
#' card) go back to scalars.
plotly_strip_build_inputs <- function(p) {
  p <- plotly::plotly_build(p)
  p$x$attrs <- NULL
  p$x$visdat <- NULL
  p$x$cur_data <- NULL
  scalar_fields <- c("hoverinfo", "hovertemplate", "legendgroup", "name", "showlegend", "visible", "mode", "type",
                     "xaxis", "yaxis")
  for (i in seq_along(p$x$data)) {
    for (f in scalar_fields) {
      v <- p$x$data[[i]][[f]]
      if (length(v) > 1 && is.atomic(v) && length(unique(v)) == 1) p$x$data[[i]][[f]] <- v[1]
    }
  }
  p
}

trajectory_link_plot <- function(p, ns) {
  # Build once so keys stay aligned with the trace splits, without the build inputs (payload).
  p <- plotly_strip_build_inputs(p)
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

# Keep paired plots equally sized when one has a longer legend. The rendered
# guides determine one shared bottom margin, in Plotly 2.11-compatible paper
# coordinates; title clearance is measured again after color changes and resize.
paired_plotly_guides <- function(p, ids) {
  htmlwidgets::onRender(p, "function(el, x, ids) {
        var prior = el._isletGuideSync;
        if (prior) {
          el.removeListener('plotly_afterplot', prior.schedule);
          if (prior.frame) window.cancelAnimationFrame(prior.frame);
        }
        var state = {frame: null};
        state.schedule = function() {
          if (state.frame) return;
          state.frame = window.requestAnimationFrame(function() {
            state.frame = null;
            var plots = ids.map(function(id) { return document.getElementById(id); });
            if (plots.some(function(plot) { return !plot || !plot._fullLayout; })) return;
            var tallest = 0;
            var clearance = 40;
            plots.forEach(function(plot) {
              plot.querySelectorAll('g.colorbar').forEach(function(guide) {
                var title = guide.querySelector('g.cbtitle > text');
                var fill = guide.querySelector('.cbfill');
                if (!title || !fill) return;
                // Reset first: repeated draws must not accumulate an offset.
                title.setAttribute('dy', '0');
                var gap = fill.getBoundingClientRect().top - title.getBoundingClientRect().bottom;
                title.setAttribute('dy', String(Math.min(0, gap - 6)));
              });
              plot.querySelectorAll('g.legend, g.colorbar').forEach(function(guide) {
                tallest = Math.max(tallest, guide.getBoundingClientRect().height);
              });
              var axisTitle = plot.querySelector('.g-xtitle');
              if (axisTitle) {
                var layout = plot._fullLayout;
                var paperBottom = plot.getBoundingClientRect().top + layout._size.t + layout._size.h;
                // Keep the guide clear of the actual axis title too; font and
                // title offsets differ between the supported Plotly versions.
                clearance = Math.max(clearance,
                  axisTitle.getBoundingClientRect().bottom - paperBottom + 8 + 0.02 * layout.height);
              }
            });
            var bottom = Math.ceil(Math.max(110, tallest + clearance));
            plots.forEach(function(plot) {
              var layout = plot._fullLayout;
              // Both supported Plotly versions default to paper coordinates;
              // Plotly.js 2.11 (the pinned image) has no yref attribute.
              // Use paper coordinates to put each guide bottom 2% above the
              // widget bottom, using the height AFTER the shared margin.
              var height = layout.height - layout.margin.t - bottom;
              if (height <= 0) return;
              var y = (0.02 * layout.height - bottom) / height;
              var update = {};
              if (Math.abs(layout.margin.b - bottom) > 0.5) update['margin.b'] = bottom;
              if (layout.legend && Math.abs(layout.legend.y - y) > 0.000001)
                update['legend.y'] = y;
              var bars = [];
              plot.data.forEach(function(trace, index) {
                if (trace.marker && trace.marker.colorbar &&
                    Math.abs(trace.marker.colorbar.y - y) > 0.000001) bars.push(index);
              });
              if (bars.length)
                Plotly.update(plot, {'marker.colorbar.y': y}, update, bars);
              else if (Object.keys(update).length) Plotly.relayout(plot, update);
            });
          });
        };
        el._isletGuideSync = state;
        el.on('plotly_afterplot', state.schedule);
        state.schedule();
      }", data = ids)
}
