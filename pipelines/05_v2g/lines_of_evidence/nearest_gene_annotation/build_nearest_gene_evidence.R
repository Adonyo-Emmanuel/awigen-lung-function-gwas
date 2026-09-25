#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_nearest_gene_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
source_file <- config_value(cfg, "annovar_nearest_gene_file")
nearest_gene <- read_required_table(source_file)

mapping <- list(
  sentinel_snp = c("CHR_POS_AMIN_AMAX"),
  variant_id = c("CHR_POS_AMIN_AMAX"),
  rsid = c("RSID"),
  chromosome_b37 = c("CHROM"),
  position_b37 = c("POS"),
  effect_allele = c("ALLELE1"),
  other_allele = c("ALLELE2"),
  trait = c("TRAIT"),
  traits = c("TRAITS"),
  gwas_p_value = c("P"),
  gwas_z_score = c("ZSCORE"),
  gene_symbol = c("Gene"),
  gene_or_protein = c("Gene"),
  distance_from_sentinel_bp = c("dist", "distance"),
  evidence_detail = c("Location", "Genes")
)

out <- standardise_source_table(
  nearest_gene,
  mapping,
  constants = list(
    evidence_line = "nearest_gene_annotation",
    evidence_source = "ANNOVAR_nearest_gene",
    evidence_class = "positional"
  ),
  source_file = source_file
)

write_standard_evidence(out, cfg, "nearest_gene_annotation.tsv")

