#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_rare_variant_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
source_file <- config_value(cfg, "rare_variant_file")
rare_variant <- read_required_table(source_file)

mapping <- list(
  sentinel_snp = c("awigen_CHR_POS_AMIN_AMAX"),
  variant_id = c("single_variant_id", "awigen_CHR_POS_AMIN_AMAX"),
  rsid = c("awigen_RSID"),
  chromosome_b37 = c("awigen_CHROM_b37"),
  position_b37 = c("awigen_POS_b37"),
  chromosome_b38 = c("awigen_CHROM_b38"),
  position_b38 = c("awigen_POS_b38"),
  effect_allele = c("awigen_ALLELE1"),
  other_allele = c("awigen_ALLELE2"),
  trait = c("single_variant_trait", "collapsing_trait", "awigen_TRAIT"),
  traits = c("awigen_TRAITS"),
  gwas_p_value = c("awigen_P"),
  gwas_z_score = c("awigen_ZSCORE"),
  gene_id = c("gene_id"),
  gene_symbol = c("gene_symbol"),
  gene_or_protein = c("gene_symbol", "gene_id"),
  molecular_trait_id = c("gene_id", "gene_symbol"),
  p_value = c("single_variant_p", "collapsing_gene_p"),
  z_score = c("single_variant_zscore"),
  distance_from_sentinel_bp = c("distance_from_sentinel_bp"),
  consequence = c("single_variant_consequence"),
  impact = c("single_variant_impact"),
  sift = c("single_variant_sift"),
  polyphen = c("single_variant_polyphen"),
  cadd_phred = c("single_variant_cadd_phred"),
  evidence_detail = c("evidence_type", "evidence_source", "collapsing_best_test")
)

out <- standardise_source_table(
  rare_variant,
  mapping,
  constants = list(
    evidence_line = "rare_variant_associations",
    evidence_source = "UK_Biobank_WES_rare_variant_association",
    evidence_class = "rare_variant_association",
    molecular_trait_type = "gene"
  ),
  source_file = source_file
)

write_standard_evidence(out, cfg, "rare_variant_associations.tsv")

