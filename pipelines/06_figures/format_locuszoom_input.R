#!/usr/bin/env Rscript

# Format the METAL results for LocusZoom: one file per trait with columns
# rsid ("chr<CHROM>:<POS>") and p. Variants are filtered to MAF >= 1% so the
# regional plots show the same variants as the Manhattan plots, and duplicate
# positions (e.g. a SNP and an indel at one position) are reduced to the
# variant with the smallest P value.
#
# Usage:
#   Rscript format_locuszoom_input.R <metal_results_dir>
#
# Writes <trait>_locuszoom_formated.txt into <metal_results_dir>.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript format_locuszoom_input.R <metal_results_dir>")
}
file_dir <- args[[1L]]

trait_files <- c(
  FEV1 = "updated_awigen_metal_FEV11.tbl",
  FVC  = "updated_awigen_metal_FVC1.tbl",
  FF   = "updated_awigen_metal_FF1.tbl",
  PEF  = "updated_awigen_metal_PEF1.tbl"
)

maf_threshold <- 0.01

for (trait in names(trait_files)) {
  input_file <- file.path(file_dir, trait_files[[trait]])
  message("Reading ", input_file)

  dt <- fread(input_file, select = c("chromosome", "position", "p", "maf"))
  dt <- dt[is.finite(p) & p > 0 & p <= 1 & is.finite(maf) & maf >= maf_threshold]
  dt[, rsid := paste0("chr", chromosome, ":", position)]

  setorder(dt, p)
  dt <- dt[!duplicated(rsid)]

  output_file <- file.path(file_dir, paste0(trait, "_locuszoom_formated.txt"))
  fwrite(dt[, .(rsid, p)], output_file, sep = "\t", quote = FALSE)
  message(trait, ": ", format(nrow(dt), big.mark = ","), " variants -> ", output_file)
}
