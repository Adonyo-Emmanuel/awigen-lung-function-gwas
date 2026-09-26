#!/usr/bin/env Rscript

# Format the METAL results for LocusZoom. Only the traits and regions listed
# in top4_sentinels_awigen.tsv are written (sentinel +/- 600 kb, slightly
# wider than the 500 kb plotting flank), so the output files are small.
#
# Output: <trait>_locuszoom_formated.txt with columns rsid ("chr<CHROM>:<POS>")
# and p. Variants are filtered to MAF >= 1% so the regional plots show the
# same variants as the Manhattan plots, and duplicate positions (e.g. a SNP
# and an indel at one position) keep the variant with the smallest P value.
#
# Usage:
#   Rscript format_locuszoom_input.R <metal_results_dir>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript format_locuszoom_input.R <metal_results_dir>")
}
file_dir <- args[[1L]]

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
sentinels <- fread(file.path(script_dir, "top4_sentinels_awigen.tsv"))
sentinels[, `:=`(
  chromosome = as.integer(sub("^chr([0-9]+):.*$", "\\1", chrpos)),
  position = as.integer(sub("^chr[0-9]+:", "", chrpos))
)]

trait_files <- c(
  FEV1 = "updated_awigen_metal_FEV11.tbl",
  FVC  = "updated_awigen_metal_FVC1.tbl",
  FF   = "updated_awigen_metal_FF1.tbl",
  PEF  = "updated_awigen_metal_PEF1.tbl"
)

maf_threshold <- 0.01
region_half_width <- 600000L

for (trait_name in unique(sentinels$trait)) {
  input_file <- file.path(file_dir, trait_files[[trait_name]])
  message("Reading ", input_file)

  dt <- fread(input_file, select = c("chromosome", "position", "p", "maf"))
  dt <- dt[is.finite(p) & p > 0 & p <= 1 & is.finite(maf) & maf >= maf_threshold]

  regions <- sentinels[trait == trait_name]
  in_region <- Reduce(`|`, lapply(seq_len(nrow(regions)), function(i) {
    dt$chromosome == regions$chromosome[[i]] &
      abs(dt$position - regions$position[[i]]) <= region_half_width
  }))
  dt <- dt[in_region]

  dt[, rsid := paste0("chr", chromosome, ":", position)]
  setorder(dt, p)
  dt <- dt[!duplicated(rsid)]

  missing <- setdiff(regions$chrpos, dt$rsid)
  if (length(missing) > 0L) {
    stop("Sentinel(s) missing after filtering for ", trait_name, ": ",
         paste(missing, collapse = ", "))
  }

  output_file <- file.path(file_dir, paste0(trait_name, "_locuszoom_formated.txt"))
  fwrite(dt[, .(rsid, p)], output_file, sep = "\t", quote = FALSE)
  message(trait_name, ": ", format(nrow(dt), big.mark = ","),
          " variants in ", nrow(regions), " region(s) -> ", output_file)
}
