#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_pqtl_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])

standardise_pqtl <- function(config_key, evidence_source, panel, mapping) {
  source_file <- config_value(cfg, config_key)
  dt <- read_required_table(source_file)
  standardise_source_table(
    dt,
    mapping,
    constants = list(
      evidence_line = "pqtl_associations",
      evidence_source = evidence_source,
      evidence_class = "pQTL",
      molecular_trait_type = "protein_abundance",
      tissue_or_protein_panel = panel
    ),
    source_file = source_file
  )
}

olink <- standardise_pqtl(
  "pqtl_olink_file",
  "UK_Biobank_Olink_African_ancestry_blood",
  "Olink",
  list(
    sentinel_snp = c("sentinel_snp"),
    variant_id = c("chr_pos_amin_amax"),
    rsid = c("rsid"),
    chromosome_b37 = c("chromosome"),
    position_b37 = c("position"),
    effect_allele = c("allele1"),
    other_allele = c("allele2"),
    trait = c("trait"),
    gwas_p_value = c("gwas_p"),
    gwas_z_score = c("gwas_zscore"),
    gwas_beta = c("gwas_beta"),
    gwas_se = c("gwas_se"),
    posterior_probability = c("postprob"),
    gene_symbol = c("protein"),
    gene_or_protein = c("protein"),
    molecular_trait_id = c("protein"),
    p_value = c("pqtl_p"),
    effect_estimate = c("pqtl_beta"),
    standard_error = c("pqtl_se"),
    evidence_detail = c("ancestry")
  )
)

masc <- standardise_pqtl(
  "pqtl_masc_file",
  "MASC_plasma",
  "MASC",
  list(
    sentinel_snp = c("sentinel_snp"),
    variant_id = c("masc_id", "chr_pos_amin_amax"),
    rsid = c("rsid"),
    chromosome_b37 = c("chromosome_b37", "masc_chrom"),
    position_b37 = c("position_b37", "masc_genpos"),
    chromosome_b38 = c("chromosome_b38"),
    position_b38 = c("position_b38"),
    effect_allele = c("allele1"),
    other_allele = c("allele2"),
    trait = c("trait"),
    gwas_p_value = c("p"),
    gwas_z_score = c("zscore"),
    gwas_beta = c("beta"),
    gwas_se = c("se"),
    posterior_probability = c("postprob"),
    gene_symbol = c("masc_protein"),
    gene_or_protein = c("masc_protein"),
    molecular_trait_id = c("masc_protein"),
    p_value = c("masc_p"),
    effect_estimate = c("masc_beta"),
    standard_error = c("masc_se"),
    evidence_detail = c("masc_test", "match_key_b37_unordered")
  )
)

out <- rbindlist(list(olink, masc), use.names = TRUE, fill = TRUE)
write_standard_evidence(out, cfg, "pqtl_associations.tsv")

