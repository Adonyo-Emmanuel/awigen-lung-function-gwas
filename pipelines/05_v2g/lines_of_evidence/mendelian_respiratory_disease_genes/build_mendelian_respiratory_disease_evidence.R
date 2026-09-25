#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_mendelian_respiratory_disease_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
source_file <- config_value(cfg, "rare_disease_file")
rare_disease <- read_required_table(source_file)

mapping <- list(
  sentinel_snp = c("CHR_POS_AMIN_AMAX", "sentinel"),
  variant_id = c("CHR_POS_AMIN_AMAX"),
  rsid = c("RSID"),
  chromosome_b37 = c("CHROM", "chr"),
  position_b37 = c("POS", "pos"),
  effect_allele = c("ALLELE1"),
  other_allele = c("ALLELE2"),
  trait = c("TRAIT"),
  traits = c("TRAITS"),
  gwas_p_value = c("P"),
  gwas_z_score = c("ZSCORE"),
  gene_symbol = c("Symbol", "Gene"),
  gene_or_protein = c("Symbol", "Gene"),
  molecular_trait_id = c("Symbol", "Gene"),
  distance_from_sentinel_bp = c("distance"),
  evidence_detail = c("Diseases", "HPOTerms", "Evidence", "Validation")
)

out <- standardise_source_table(
  rare_disease,
  mapping,
  constants = list(
    evidence_line = "mendelian_respiratory_disease_genes",
    evidence_source = "Orphanet_respiratory_disease_gene",
    evidence_class = "Mendelian_disease",
    molecular_trait_type = "gene"
  ),
  source_file = source_file
)

write_standard_evidence(out, cfg, "mendelian_respiratory_disease_genes.tsv")

