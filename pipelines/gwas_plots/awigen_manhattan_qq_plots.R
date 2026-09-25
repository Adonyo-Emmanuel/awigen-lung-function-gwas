#!/usr/bin/env Rscript

# Publication Manhattan and QQ plots for the four AWI-Gen lung function
# traits, drawn directly with ggplot2 (no topr dependency).
#
# Manhattan plots
#   * Autosomal SNPs with MAF >= 1%.
#   * Chromosomes alternate blue and grey.
#   * The genome-wide significant loci listed in `genome_wide_annotations`
#     are highlighted in red (SNPs within +/- `locus_half_window` of the
#     sentinel with P < `locus_highlight_p`) and the sentinel is labelled with
#     its nearest gene in red italics.
#   * SNPs with P > `downsample_p` are randomly thinned for plotting only;
#     every SNP below that P value is drawn.
#
# QQ plots
#   * Observed -log10(P) against expected -log10(P) under the null, using
#     every autosomal SNP tested in the meta-analysis with a valid P value.
#   * Expected quantiles are (i - 0.5) / n; the grey band is the pointwise 95%
#     confidence interval from the Beta(i, n - i + 1) order-statistic
#     distribution.
#   * lambda GC (median chi-squared / 0.4549) is printed on each panel.
#   * Overlapping points in the dense null region are removed after rounding
#     to the plotting resolution, which does not change the appearance.
#
# Usage:
#   Rscript awigen_manhattan_qq_plots.R [input_dir] [output_dir]

user_library <- path.expand("~/R/library")
if (dir.exists(user_library)) {
  .libPaths(c(user_library, .libPaths()))
}

required_packages <- c("data.table", "ggplot2", "ggrepel", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them with install.packages() before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1L) args[[1L]] else
  "/data/gen1/ADONYO_PhD/PhD_Resources/AWIGEN_METAL_RESULTS/UPDATED"

output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(input_dir, "publication_manhattan_qq_plots")

manhattan_dir <- file.path(output_dir, "manhattan")
qq_dir <- file.path(output_dir, "qq")

for (directory in c(manhattan_dir, qq_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

trait_files <- c(
  FEV1 = "updated_awigen_metal_FEV11.tbl",
  FVC  = "updated_awigen_metal_FVC1.tbl",
  FF   = "updated_awigen_metal_FF1.tbl",
  PEF  = "updated_awigen_metal_PEF1.tbl"
)

plot_order <- c("FEV1", "FVC", "FF", "PEF")

trait_labels <- list(
  FEV1 = expression(FEV[1]),
  FVC  = expression(FVC),
  FF   = expression(FEV[1] / FVC),
  PEF  = expression(PEF)
)

trait_panel_tags <- list(
  FEV1 = expression(bold("a") ~ ~FEV[1]),
  FVC  = expression(bold("b") ~ ~FVC),
  FF   = expression(bold("c") ~ ~FEV[1] / FVC),
  PEF  = expression(bold("d") ~ ~PEF)
)

# These are the nearest-gene assignments from the final AWI-Gen V2G table.
# TMEM266 and LINC01899 are the current symbols for C15orf27 and
# LOC102724913, respectively.
genome_wide_annotations <- data.table(
  trait = c("FF", "FF", "FEV1", "FF"),
  CHROM = c(8L, 11L, 15L, 18L),
  POS = c(54811191L, 98327692L, 76357354L, 69622284L),
  sentinel = c(
    "8_54811191_A_G",
    "11_98327692_A_G",
    "15_76357354_A_G",
    "18_69622284_A_G"
  ),
  gene = c("RGS20", "CNTN5", "TMEM266", "LINC01899")
)

genome_wide_threshold <- 5e-8
suggestive_threshold <- 5e-6
maf_threshold <- 0.01

# Manhattan appearance.
chromosome_colors <- c("#1F4E99", "#A6A6A6")  # odd = blue, even = grey
locus_color <- "#D7191C"
locus_half_window <- 500000L
locus_highlight_p <- 1e-3
downsample_p <- 0.05
downsample_keep <- 0.10
chromosome_gap <- 15e6

# QQ appearance.
qq_point_color <- "#1F4E99"
qq_line_color <- "#D7191C"
qq_band_fill <- "#D9D9D9"

# Output formats: PNG for review, PDF (vector axes/text) and TIFF for journals.
output_formats <- c("png", "pdf")
plot_dpi <- 600

base_font <- "sans"

set.seed(2026)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

read_trait_data <- function(trait, path) {
  if (!file.exists(path)) {
    stop("Input file not found for ", trait, ": ", path)
  }

  dat <- fread(
    path,
    select = c("chromosome", "position", "p", "maf"),
    showProgress = FALSE
  )

  setnames(
    dat,
    c("chromosome", "position", "p", "maf"),
    c("CHROM", "POS", "P", "MAF")
  )

  dat[, `:=`(
    CHROM = suppressWarnings(as.integer(sub("^chr", "", CHROM))),
    POS = as.integer(POS),
    P = as.numeric(P),
    MAF = as.numeric(MAF)
  )]

  invalid_p <- dat[is.finite(P) & (P <= 0 | P > 1), .N]
  if (invalid_p > 0L) {
    stop(
      trait, " contains ", invalid_p,
      " finite P values outside the interval (0,1]."
    )
  }

  dat <- dat[
    !is.na(CHROM) & CHROM %between% c(1L, 22L) &
      !is.na(POS) &
      is.finite(P) & P > 0 & P <= 1
  ]

  if (nrow(dat) == 0L) {
    stop("No tested autosomal SNPs with valid P values remained for ", trait, ".")
  }

  dat
}

# Chromosome offsets are computed once from all traits so that the x axes of
# the stacked Manhattan panels line up exactly.
compute_chromosome_offsets <- function(chrom_max) {
  chrom_max <- chrom_max[, .(max_pos = max(max_pos)), by = CHROM]
  setorder(chrom_max, CHROM)
  chrom_max[, offset := c(0, cumsum(as.numeric(max_pos) + chromosome_gap)[-.N])]
  chrom_max[, center := offset + max_pos / 2]
  chrom_max
}

validate_sentinels <- function(trait, dat) {
  trait_name <- trait
  expected <- genome_wide_annotations[trait == trait_name]

  if (nrow(expected) == 0L) {
    return(expected)
  }

  observed <- dat[expected, on = .(CHROM, POS), .(
    sentinel = i.sentinel, gene = i.gene, CHROM = i.CHROM, POS = i.POS, P = x.P
  )]

  missing_sentinels <- observed[is.na(P), sentinel]
  if (length(missing_sentinels) > 0L) {
    stop(
      "Expected genome-wide sentinel(s) missing after filtering for ", trait,
      ": ", paste(missing_sentinels, collapse = ", ")
    )
  }

  non_significant <- observed[P >= genome_wide_threshold, sentinel]
  if (length(non_significant) > 0L) {
    stop(
      "Expected sentinel(s) do not satisfy P<5e-8 for ", trait, ": ",
      paste(non_significant, collapse = ", ")
    )
  }

  observed
}

prepare_manhattan_data <- function(trait, dat, offsets) {
  mdat <- dat[is.finite(MAF) & MAF >= maf_threshold, .(CHROM, POS, P)]

  if (nrow(mdat) == 0L) {
    stop("No SNPs with MAF>=1% remained for the ", trait, " Manhattan plot.")
  }

  sentinels <- validate_sentinels(trait, mdat)

  # Thin the null bulk for plotting only.
  mdat <- mdat[P <= downsample_p | stats::runif(.N) < downsample_keep]

  mdat[offsets, on = "CHROM", BP_CUM := POS + i.offset]
  mdat[, `:=`(
    LOGP = -log10(P),
    COLOR_GROUP = ifelse(CHROM %% 2L == 1L, "odd", "even")
  )]

  if (nrow(sentinels) > 0L) {
    for (i in seq_len(nrow(sentinels))) {
      mdat[
        CHROM == sentinels$CHROM[[i]] &
          abs(POS - sentinels$POS[[i]]) <= locus_half_window &
          P < locus_highlight_p,
        COLOR_GROUP := "locus"
      ]
    }
    sentinels[offsets, on = "CHROM", BP_CUM := POS + i.offset]
    sentinels[, LOGP := -log10(P)]
  }

  # Draw highlighted loci last so they sit on top of the grey/blue points.
  mdat[, draw_order := fifelse(COLOR_GROUP == "locus", 2L, 1L)]
  setorder(mdat, draw_order, BP_CUM)

  list(points = mdat, sentinels = sentinels)
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

publication_theme <- function(base_size = 8) {
  theme_classic(base_family = base_font, base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 1, hjust = 0, vjust = 1,
        margin = margin(0, 0, 1.5, 0, unit = "mm")
      ),
      plot.title.position = "plot",
      legend.position = "none",
      axis.title = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.margin = margin(3, 3, 2, 2, unit = "mm")
    )
}

create_manhattan_plot <- function(trait, mplot, offsets, y_max) {
  pts <- mplot$points
  sentinels <- mplot$sentinels

  p <- ggplot(pts, aes(x = BP_CUM, y = LOGP, colour = COLOR_GROUP)) +
    geom_hline(
      yintercept = -log10(suggestive_threshold),
      colour = "#7F7F7F", linetype = "dashed", linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(genome_wide_threshold),
      colour = locus_color, linetype = "dashed", linewidth = 0.35
    ) +
    geom_point(size = 0.45, shape = 16, stroke = 0) +
    scale_colour_manual(
      values = c(
        odd = chromosome_colors[[1L]],
        even = chromosome_colors[[2L]],
        locus = locus_color
      )
    ) +
    scale_x_continuous(
      breaks = offsets$center,
      labels = offsets$CHROM,
      limits = c(-chromosome_gap / 2,
                 max(offsets$offset + offsets$max_pos) + chromosome_gap / 2),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = scales::breaks_pretty(n = 5),
      expand = c(0, 0)
    ) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(P))),
      title = trait_labels[[trait]]
    ) +
    coord_cartesian(clip = "off") +
    publication_theme() +
    theme(axis.text.x = element_text(size = 6))

  if (nrow(sentinels) > 0L) {
    p <- p +
      geom_point(
        data = sentinels,
        aes(x = BP_CUM, y = LOGP),
        inherit.aes = FALSE,
        shape = 23, size = 1.8, stroke = 0.35,
        colour = "black", fill = locus_color
      ) +
      geom_text_repel(
        data = sentinels,
        aes(x = BP_CUM, y = LOGP, label = gene),
        inherit.aes = FALSE,
        colour = locus_color,
        fontface = "italic",
        family = base_font,
        size = 2.6,
        nudge_y = 0.9,
        direction = "x",
        min.segment.length = 0,
        segment.size = 0.25,
        segment.colour = "#4D4D4D",
        box.padding = 0.3,
        max.overlaps = Inf,
        seed = 2026
      )
  }

  p
}

calculate_lambda_gc <- function(p_values) {
  stats::median(stats::qchisq(p_values, df = 1, lower.tail = FALSE)) /
    stats::qchisq(0.5, df = 1, lower.tail = FALSE)
}

prepare_qq_data <- function(p_values) {
  n <- length(p_values)
  observed <- sort(p_values)
  i <- seq_len(n)

  qq <- data.table(
    i = i,
    expected = -log10((i - 0.5) / n),
    observed = -log10(observed)
  )

  # Remove points that would overlap at plotting resolution. Rank order is
  # preserved, so the tail (the informative part) is untouched.
  qq[, key := paste(round(expected, 3), round(observed, 3))]
  qq <- qq[!duplicated(key)]
  qq[, key := NULL]

  qq[, `:=`(
    ci_lower = -log10(stats::qbeta(0.975, i, n - i + 1)),
    ci_upper = -log10(stats::qbeta(0.025, i, n - i + 1))
  )]

  qq[]
}

create_qq_plot <- function(trait, qq, lambda_gc, n_snps) {
  axis_max <- ceiling(max(qq$expected, qq$observed, qq$ci_upper) + 0.2)

  label <- sprintf(
    "lambda[GC] == %s",
    formatC(lambda_gc, format = "f", digits = 3)
  )
  n_label <- sprintf(
    "italic(n) == \"%s\"",
    format(n_snps, big.mark = ",")
  )

  ggplot(qq, aes(x = expected, y = observed)) +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper),
      fill = qq_band_fill, colour = NA
    ) +
    geom_abline(
      intercept = 0, slope = 1,
      colour = qq_line_color, linewidth = 0.4
    ) +
    geom_point(colour = qq_point_color, size = 0.6, shape = 16, stroke = 0) +
    annotate(
      "text", x = 0.97 * axis_max, y = 0.13 * axis_max,
      label = label, parse = TRUE,
      hjust = 1, vjust = 0, size = 3.2, family = base_font
    ) +
    annotate(
      "text", x = 0.97 * axis_max, y = 0.05 * axis_max,
      label = n_label, parse = TRUE,
      hjust = 1, vjust = 0, size = 2.6, family = base_font,
      colour = "#4D4D4D"
    ) +
    scale_x_continuous(
      limits = c(0, axis_max), expand = expansion(mult = c(0, 0.02)),
      breaks = scales::breaks_pretty(n = 5)
    ) +
    scale_y_continuous(
      limits = c(0, axis_max), expand = expansion(mult = c(0, 0.02)),
      breaks = scales::breaks_pretty(n = 5)
    ) +
    coord_equal(clip = "off") +
    labs(
      x = expression(Expected ~ -log[10](italic(P))),
      y = expression(Observed ~ -log[10](italic(P))),
      title = trait_labels[[trait]]
    ) +
    publication_theme(base_size = 9)
}

save_plot <- function(plot, stem, width_mm, height_mm) {
  for (fmt in output_formats) {
    filename <- paste0(stem, ".", fmt)
    device <- switch(
      fmt,
      pdf = grDevices::cairo_pdf,
      png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
      tiff = if (requireNamespace("ragg", quietly = TRUE)) {
        function(...) ragg::agg_tiff(..., compression = "lzw")
      } else {
        function(...) grDevices::tiff(..., compression = "lzw")
      },
      fmt
    )
    ggsave(
      filename = filename,
      plot = plot,
      device = device,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = plot_dpi,
      bg = "white",
      limitsize = FALSE
    )
    message("  wrote ", filename)
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

trait_data <- list()
chrom_max <- list()

for (trait in plot_order) {
  input_file <- file.path(input_dir, trait_files[[trait]])
  message("Reading ", trait, " from ", input_file)
  trait_data[[trait]] <- read_trait_data(trait, input_file)
  chrom_max[[trait]] <- trait_data[[trait]][, .(max_pos = max(POS)), by = CHROM]
}

offsets <- compute_chromosome_offsets(rbindlist(chrom_max))

manhattan_data <- list()
qq_plots <- list()

for (trait in plot_order) {
  dat <- trait_data[[trait]]

  lambda_gc <- calculate_lambda_gc(dat$P)
  message(
    trait, ": ", format(nrow(dat), big.mark = ","),
    " tested SNPs used for QQ; lambda GC = ", sprintf("%.3f", lambda_gc)
  )

  qq <- prepare_qq_data(dat$P)
  qq_plots[[trait]] <- create_qq_plot(trait, qq, lambda_gc, nrow(dat))
  save_plot(
    qq_plots[[trait]],
    file.path(qq_dir, paste0("AWI-Gen_qq_", trait)),
    width_mm = 85, height_mm = 85
  )

  manhattan_data[[trait]] <- prepare_manhattan_data(trait, dat, offsets)

  trait_data[[trait]] <- NULL
  rm(dat, qq)
  invisible(gc())
}

# A common y-axis ceiling across traits keeps the stacked panels comparable.
y_max <- max(
  vapply(manhattan_data, function(x) max(x$points$LOGP), numeric(1)),
  -log10(genome_wide_threshold) + 1
) * 1.15

manhattan_plots <- list()
for (trait in plot_order) {
  manhattan_plots[[trait]] <- create_manhattan_plot(
    trait, manhattan_data[[trait]], offsets, y_max
  )
  save_plot(
    manhattan_plots[[trait]],
    file.path(manhattan_dir, paste0("AWI-Gen_manhattan_", trait)),
    width_mm = 180, height_mm = 65
  )
}

# Combined Manhattan figure: panels a-d stacked, x-axis only on the bottom.
combined_manhattan <- lapply(plot_order, function(trait) {
  p <- manhattan_plots[[trait]] + labs(title = trait_panel_tags[[trait]])
  if (trait != utils::tail(plot_order, 1L)) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  p
})

save_plot(
  wrap_plots(combined_manhattan, ncol = 1L),
  file.path(manhattan_dir, "AWI-Gen_manhattan_all_traits"),
  width_mm = 180, height_mm = 200
)

# Combined QQ figure: 2 x 2 grid.
combined_qq <- lapply(plot_order, function(trait) {
  qq_plots[[trait]] + labs(title = trait_panel_tags[[trait]])
})

save_plot(
  wrap_plots(combined_qq, ncol = 2L),
  file.path(qq_dir, "AWI-Gen_qq_all_traits"),
  width_mm = 170, height_mm = 170
)

message("All publication plots written under: ", output_dir)
