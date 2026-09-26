#!/usr/bin/env Rscript

# Combine the four untitled LocusZoom regional plots into one 2 x 2 figure
# with panel letters a-d, ordered by chromosome.
#
# Input: the untitled PNGs written by awigen_locuszoom_plots.sh to
#   <plots_dir>/panels/awigen_<trait>_<gene>_*.png
# Output: <plots_dir>/AWI-Gen_regional_plots_combined.{pdf,png,tiff}
#
# Usage:
#   Rscript combine_regional_plots.R <metal_results_dir> [plots_dir]
#   (plots_dir defaults to <metal_results_dir>/awigen_regional_assoc_plots)
#
# Needs the R package "png" (install.packages("png")).

user_library <- path.expand("~/R/library")
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

if (!requireNamespace("png", quietly = TRUE)) {
  stop('R package "png" is missing. Install it with install.packages("png").')
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript combine_regional_plots.R <metal_results_dir> [plots_dir]")
}
plots_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(args[[1L]], "awigen_regional_assoc_plots")
panel_dir <- file.path(plots_dir, "panels")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
sentinels <- utils::read.delim(
  file.path(script_dir, "top4_sentinels_awigen.tsv"),
  stringsAsFactors = FALSE
)

# Panels a-d in chromosome order.
sentinels$chromosome <- as.integer(sub("^chr([0-9]+):.*$", "\\1", sentinels$chrpos))
sentinels$position <- as.integer(sub("^chr[0-9]+:", "", sentinels$chrpos))
sentinels <- sentinels[order(sentinels$chromosome, sentinels$position), ]

panel_files <- vapply(seq_len(nrow(sentinels)), function(i) {
  pattern <- paste0("^awigen_", sentinels$trait[[i]], "_", sentinels$gene[[i]], "_.*\\.png$")
  hits <- list.files(panel_dir, pattern = pattern, full.names = TRUE)
  if (length(hits) != 1L) {
    stop("Expected one untitled PNG for ", sentinels$gene[[i]], " in ", panel_dir,
         " but found ", length(hits), ". Run awigen_locuszoom_plots.sh first.")
  }
  hits[[1L]]
}, character(1))

images <- lapply(panel_files, png::readPNG)
for (i in seq_along(panel_files)) {
  message(letters[[i]], ": ", sentinels$gene[[i]], " (chr", sentinels$chromosome[[i]],
          ")  <- ", basename(panel_files[[i]]))
}

# Figure size: 180 mm wide (two 90 mm columns); height follows the panels'
# aspect ratio, plus a small strip above each row for the panel letter.
ncol <- 2L
nrow <- ceiling(length(images) / ncol)
width_mm <- 180
panel_width_mm <- width_mm / ncol
aspect <- max(vapply(images, function(img) dim(img)[1] / dim(img)[2], numeric(1)))
letter_strip_mm <- 5
row_height_mm <- panel_width_mm * aspect + letter_strip_mm
height_mm <- nrow * row_height_mm
dpi <- 600

draw_figure <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  layout <- grid::grid.layout(nrow = nrow, ncol = ncol)
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(images)) {
    r <- (i - 1L) %/% ncol + 1L
    c <- (i - 1L) %% ncol + 1L
    grid::pushViewport(grid::viewport(layout.pos.row = r, layout.pos.col = c))
    grid::grid.raster(
      images[[i]],
      y = grid::unit(0, "npc"), just = "bottom",
      height = grid::unit(1, "npc") - grid::unit(letter_strip_mm, "mm"),
      interpolate = TRUE
    )
    grid::grid.text(
      letters[[i]],
      x = grid::unit(2, "mm"), y = grid::unit(1, "npc") - grid::unit(1, "mm"),
      just = c("left", "top"),
      gp = grid::gpar(fontface = "bold", fontsize = 12, fontfamily = "sans")
    )
    grid::popViewport()
  }
  grid::popViewport()
}

stem <- file.path(plots_dir, "AWI-Gen_regional_plots_combined")
w_in <- width_mm / 25.4
h_in <- height_mm / 25.4
use_ragg <- requireNamespace("ragg", quietly = TRUE)

grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in)
draw_figure()
invisible(grDevices::dev.off())

if (use_ragg) {
  ragg::agg_png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in", res = dpi)
} else {
  grDevices::png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in", res = dpi)
}
draw_figure()
invisible(grDevices::dev.off())

if (use_ragg) {
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w_in, height = h_in, units = "in",
                 res = dpi, compression = "lzw")
} else {
  grDevices::tiff(paste0(stem, ".tiff"), width = w_in, height = h_in, units = "in",
                  res = dpi, compression = "lzw")
}
draw_figure()
invisible(grDevices::dev.off())

message("Combined figure written to ", stem, ".{pdf,png,tiff}")
