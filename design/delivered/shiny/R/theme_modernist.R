# ---------------------------------------------------------------------------
# ggplot theme + palettes
#
# Font: Archivo must be registered before the first plot draws, or ggplot
# silently falls back to sans. Call register_modernist_fonts() once in
# global.R / app.R (needs sysfonts + showtext, or a system install + ragg).
#
# The donor and phenotype palettes are the app's existing scientific colors,
# copied verbatim — the redesign does NOT touch data-encoding color.
# ---------------------------------------------------------------------------

IE_INK      <- "#201e1d"
IE_GROUND   <- "#f3f2f2"
IE_RULE     <- "#c9c5c2"
IE_MUTED    <- "#6f6a67"
IE_CRIMSON  <- "#b71533"
IE_GOLD     <- "#c8a02e"

PAL_DONOR <- c("ND" = "#3366CC", "Aab+" = "#FF6600", "T1D" = "#109618")

PAL_PHENO <- c(
  Beta = "#D32F2F", Alpha = "#FF8F00", Delta = "#FDD835",
  Acinar = "#00ACC1", Ductal = "#43A047", Structural = "#9E9E9E",
  `CD8 T` = "#7B1FA2", `CD4 T` = "#1565C0", Macrophage = "#5C6BC0",
  Endothelial = "#FF7043", SMA = "#F06292", Vessel = "#E91E63"
)

PAL_LEIDEN <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377")

register_modernist_fonts <- function() {
  if (!requireNamespace("sysfonts", quietly = TRUE)) return(invisible(FALSE))
  ok <- try(sysfonts::font_add_google("Archivo", "Archivo"), silent = TRUE)
  if (inherits(ok, "try-error")) return(invisible(FALSE))
  if (requireNamespace("showtext", quietly = TRUE)) {
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 96)
  }
  invisible(TRUE)
}

#' Flat, flush-left, 2px-rule plot theme.
#' @param grid which gridlines to keep — Modernist wants structure visible but
#'   thin; the axis lines carry the weight.
theme_modernist <- function(base_size = 12, grid = "y", family = "Archivo") {
  ggplot2::theme_minimal(base_size = base_size, base_family = family) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = IE_GROUND, colour = NA),
      panel.background  = ggplot2::element_rect(fill = IE_GROUND, colour = NA),
      panel.border      = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) {
        ggplot2::element_line(colour = IE_RULE, linewidth = 0.3)
      } else ggplot2::element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) {
        ggplot2::element_line(colour = IE_RULE, linewidth = 0.3)
      } else ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(colour = IE_INK, linewidth = 0.7),
      axis.ticks        = ggplot2::element_line(colour = IE_INK, linewidth = 0.5),
      axis.ticks.length = grid::unit(4, "pt"),
      axis.text         = ggplot2::element_text(colour = IE_MUTED, size = base_size - 2),
      axis.title        = ggplot2::element_text(
        colour = IE_INK, size = base_size - 2, hjust = 0,
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 3, hjust = 0, colour = IE_INK,
        margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1, hjust = 0, colour = IE_MUTED,
        margin = ggplot2::margin(b = 10)
      ),
      plot.caption = ggplot2::element_text(size = base_size - 3, hjust = 0,
                                           colour = IE_MUTED),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      legend.position   = "bottom",
      legend.justification = "left",
      legend.title      = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = base_size - 2, colour = IE_INK),
      strip.background  = ggplot2::element_rect(fill = "#e4e1df", colour = NA),
      strip.text        = ggplot2::element_text(
        face = "bold", size = base_size - 2, hjust = 0, colour = IE_INK,
        margin = ggplot2::margin(4, 6, 4, 6)
      ),
      plot.margin = ggplot2::margin(6, 10, 4, 4)
    )
}

scale_colour_donor <- function(...) {
  ggplot2::scale_colour_manual(values = PAL_DONOR, drop = FALSE, ...)
}
scale_fill_donor <- function(...) {
  ggplot2::scale_fill_manual(values = PAL_DONOR, drop = FALSE, ...)
}
scale_colour_pheno <- function(...) {
  ggplot2::scale_colour_manual(values = PAL_PHENO, drop = FALSE, ...)
}
scale_colour_leiden <- function(...) {
  ggplot2::scale_colour_manual(values = PAL_LEIDEN, ...)
}

# Emphasis scale for stats output: crimson = significant, grey = n.s.
scale_colour_significance <- function(...) {
  ggplot2::scale_colour_manual(
    values = c(`TRUE` = IE_CRIMSON, `FALSE` = "#9b9693"),
    labels = c(`TRUE` = "significant", `FALSE` = "n.s."), ...
  )
}
