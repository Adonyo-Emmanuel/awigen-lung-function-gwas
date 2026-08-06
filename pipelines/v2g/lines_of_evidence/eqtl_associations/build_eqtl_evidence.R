#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_eqtl_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])

standardise_eqtl <- function(config_key, evidence_source, panel, mapping) {
  source_file <- config_value(cfg, config_key)
  dt <- read_required_table(source_file)
  standardise_source_table(
    dt,
    mapping,
    constants = list(
      evidence_line = "eqtl_associations",
      evidence_source = evidence_source,
      evidence_class = "eQTL",
      molecular_trait_type = "gene_expression",
      tissue_or_protein_panel = panel
    ),
    source_file = source_file
  )
}

eqtlgen <- standardise_eqtl(
  "eqtlgen_file",
  "eQTLGen_blood",
  "blood",
  list(
    sentinel_snp = c("sentinel_snp"),
    variant_id = c("SNP", "chr_pos_amin_amax"),
    rsid = c("rsid"),
    chromosome_b37 = c("chromosome_b37"),
    position_b37 = c("position_b37"),
    chromosome_b38 = c("chromosome_b38", "SNPChr"),
    position_b38 = c("position_b38", "SNPPos"),
    effect_allele = c("allele1"),
    other_allele = c("allele2"),
    trait = c("trait"),
    gwas_p_value = c("p"),
    gwas_z_score = c("zscore"),
    posterior_probability = c("postprob"),
    gene_id = c("Gene"),
    gene_symbol = c("GeneSymbol"),
    gene_or_protein = c("GeneSymbol", "Gene"),
    molecular_trait_id = c("Gene"),
    p_value = c("Pvalue"),
    fdr = c("FDR"),
    threshold_p_value = c("BonferroniP"),
    z_score = c("eqtlgen_zscore_raw", "Zscore"),
    evidence_detail = c("eqtlgen_allele_alignment")
  )
)

gtex <- standardise_eqtl(
  "gtex_file",
  "GTEx_v10",
  NA_character_,
  list(
    sentinel_snp = c("sentinel_snp"),
    variant_id = c("id_b38", "chr_pos_amin_amax"),
    rsid = c("rsid"),
    chromosome_b37 = c("chromosome_b37"),
    position_b37 = c("position_b37"),
    chromosome_b38 = c("chromosome_b38", "chr"),
    position_b38 = c("position_b38", "pos"),
    effect_allele = c("allele1"),
    other_allele = c("allele2"),
    trait = c("trait"),
    gwas_p_value = c("p"),
    gwas_z_score = c("zscore"),
    posterior_probability = c("postprob"),
    gene_id = c("gene_id"),
    gene_symbol = c("gene_name"),
    gene_or_protein = c("gene_name", "gene_id"),
    molecular_trait_id = c("gene_id"),
    tissue_or_protein_panel = c("tissue"),
    p_value = c("pval_nominal"),
    threshold_p_value = c("pval_nominal_threshold"),
    effect_estimate = c("gtex_slope_aligned_to_allele1", "slope"),
    standard_error = c("slope_se"),
    z_score = c("gtex_z_aligned_to_allele1", "gtex_z_from_slope"),
    distance_from_sentinel_bp = c("tss_distance"),
    evidence_detail = c("gtex_allele_alignment")
  )
)

sabr <- standardise_eqtl(
  "sabr_file",
  "SABR_blood",
  "blood",
  list(
    sentinel_snp = c("gwas_sentinel_snp"),
    variant_id = c("sabr_variant_id", "gwas_chr_pos_amin_amax"),
    rsid = c("gwas_rsid"),
    chromosome_b37 = c("gwas_chromosome_b37"),
    position_b37 = c("gwas_position_b37"),
    chromosome_b38 = c("gwas_chromosome_b38", "sabr_locus_contig"),
    position_b38 = c("gwas_position_b38", "sabr_locus_position"),
    effect_allele = c("gwas_allele1"),
    other_allele = c("gwas_allele2"),
    trait = c("gwas_trait"),
    gwas_p_value = c("gwas_p"),
    gwas_z_score = c("gwas_zscore"),
    posterior_probability = c("gwas_postprob"),
    gene_id = c("sabr_gene_id"),
    gene_symbol = c("sabr_s7_gene_symbol"),
    gene_or_protein = c("sabr_s7_gene_symbol", "sabr_gene_id"),
    molecular_trait_id = c("sabr_gene_id"),
    p_value = c("sabr_pval_nominal"),
    fdr = c("sabr_fdr", "sabr_s7_qval"),
    q_value = c("sabr_s7_qval"),
    threshold_p_value = c("sabr_s7_pval_nominal_threshold"),
    effect_estimate = c("sabr_slope_aligned_to_gwas_allele1", "sabr_slope"),
    standard_error = c("sabr_slope_se"),
    z_score = c("sabr_z_aligned_to_gwas_allele1", "sabr_z_from_slope", "sabr_sabr_z"),
    distance_from_sentinel_bp = c("sabr_tss_distance"),
    evidence_detail = c("sabr_allele_alignment")
  )
)

african_american <- standardise_eqtl(
  "amr_file",
  "African_American_blood_cis_eQTL",
  "blood",
  list(
    sentinel_snp = c("sentinel_snp"),
    variant_id = c("Variant", "chr_pos_amin_amax"),
    rsid = c("rsid"),
    chromosome_b37 = c("chromosome_b37"),
    position_b37 = c("position_b37"),
    chromosome_b38 = c("chromosome_b38"),
    position_b38 = c("position_b38"),
    effect_allele = c("allele1"),
    other_allele = c("allele2"),
    trait = c("trait"),
    gwas_p_value = c("p"),
    gwas_z_score = c("zscore"),
    posterior_probability = c("postprob"),
    gene_id = c("Gene Ensembl ID"),
    gene_symbol = c("Gene Symbol"),
    gene_or_protein = c("Gene Symbol", "Gene Ensembl ID"),
    molecular_trait_id = c("Gene Ensembl ID"),
    p_value = c("P-Value"),
    fdr = c("FDR"),
    effect_estimate = c("Beta"),
    z_score = c("Zscore"),
    evidence_detail = c("direction")
  )
)

out <- rbindlist(list(eqtlgen, gtex, sabr, african_american), use.names = TRUE, fill = TRUE)
write_standard_evidence(out, cfg, "eqtl_associations.tsv")

