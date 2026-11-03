#!/usr/bin/env Rscript
# Real renderer checks: marker signs must be present, raised, and distinguishable.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[1]), capabilities("cairo"))
out <- args[1]
scratch <- tempfile("raster-font-check-")
dir.create(scratch)
faces <- c(plain = 1L, bold = 2L, italic = 3L, bolditalic = 4L)
render <- function(text, face, device, file, width = 160L, height = 80L) {
  if (device == "ragg") ragg::agg_png(file, width, height, background = "white", res = 96) else
    grDevices::png(file, width, height, bg = "white", res = 96, type = "cairo")
  grid::grid.text(text, gp = grid::gpar(fontfamily = "sans", fontface = face, fontsize = 28))
  grDevices::dev.off()
  pixels <- png::readPNG(file)
  apply(pixels[, , 1:3, drop = FALSE], c(1, 2), min) < 0.8
}
checks <- list()
for (name in names(faces)) {
  face <- faces[[name]]
  glyphs <- systemfonts::glyph_info(c("\u207A", "\u207B"), family = "sans",
                                  weight = if (face %in% c(2L, 4L)) "bold" else "normal",
                                  italic = face %in% c(3L, 4L))
  stopifnot(all(glyphs$index > 0L), all(glyphs$width > 0), all(glyphs$height > 0))
  for (device in c("ragg", "cairo")) {
    maps <- lapply(c("\u207A", "\u207B", "+", "-"), function(sign)
      render(sign, face, device, file.path(scratch, "glyph.png")))
    stopifnot(all(vapply(maps, any, logical(1))), !identical(maps[[1]], maps[[2]]))
    tops <- vapply(maps, function(m) min(which(m, arr.ind = TRUE)[, 1]), integer(1))
    stopifnot(tops[1] < tops[3], tops[2] < tops[4])
    checks[[paste(device, name, sep = "_")]] <- list(
      superscript_glyph_indices = glyphs$index, superscript_top_px = tops[1:2],
      baseline_top_px = tops[3:4], nonempty = TRUE, distinct = TRUE, raised = TRUE)
  }
  grDevices::cairo_pdf(file.path(out, paste0("cairo-", name, ".pdf")), width = 5, height = 1.5)
  grid::grid.text("CD68\u207ACD163\u207BLGALS3\u207ACD4\u207B",
                  gp = grid::gpar(fontfamily = "sans", fontface = face, fontsize = 20))
  grDevices::dev.off()
}

# All twelve actual macrophage profile combinations must stay distinct as pixels.
profiles <- expand.grid(pair = c("CD68+CD163+", "CD68+CD163-", "CD68-CD163+"),
                        lgals3 = c("LGALS3+", "LGALS3-"), cd4 = c("CD4+", "CD4-"))
labels <- paste0("Macrophage ", profiles$pair, profiles$lgals3, profiles$cd4)
labels <- chartr("+-", "\u207A\u207B", labels)
profile_checks <- list()
for (device in c("ragg", "cairo")) {
  maps <- lapply(labels, function(label)
    render(label, 1L, device, file.path(scratch, "profile.png"), width = 1000L))
  hashes <- vapply(maps, digest::digest, character(1), algo = "sha256")
  stopifnot(length(unique(hashes)) == 12L)
  file <- file.path(out, paste0("profiles-", device, ".png"))
  if (device == "ragg") ragg::agg_png(file, 1000, 560, background = "white", res = 96) else
    grDevices::png(file, 1000, 560, bg = "white", res = 96, type = "cairo")
  specimen <- c(labels, "CHGA\u207A hormone-negative", "T Cell CD8\u207A", "T Cell CD4\u207A")
  grid::grid.text(specimen, x = 0.03, y = seq(0.96, 0.04, length.out = length(specimen)),
                  just = "left", gp = grid::gpar(fontfamily = "sans", fontsize = 20))
  grDevices::dev.off()
  profile_checks[[device]] <- list(distinct_profiles = length(unique(hashes)),
                                   labels = labels, bitmap_sha256 = hashes)
}
unlink(scratch, recursive = TRUE)
jsonlite::write_json(list(schema = "islet.raster_font_validation/1", result = "pass",
  default_sans = systemfonts::match_fonts("sans")$path,
  faces_and_devices = checks, macrophage_profiles = profile_checks),
  file.path(out, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Raster font checks passed: four faces, ragg/Cairo PNG, twelve distinct profiles.\n")
