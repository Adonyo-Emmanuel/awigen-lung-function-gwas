#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_credible_set_variant_annotation_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
source_file <- config_value(cfg, "vep_pp50_gene_annotation_file")
vep <- read_required_table(source_file)

pp_col <- first_existing_col(vep, c("postprob", "PP", "posterior_probability"))
if (!is.na(pp_col)) {
  vep <- vep[num(vep[[pp_col]]) > 0.5]
}

mapping <- list(
  sentinel_snp = c("sentinel_snp"),
  variant_id = c("#Uploaded_variation", "chr_pos_amin_amax"),
  rsid = c("rsid", "Existing_variation"),
  chromosome_b37 = c("chromosome"),
  position_b37 = c("position"),
  effect_allele = c("allele1"),
  other_allele = c("allele2"),
  trait = c("trait"),
  traits = c("traits", "trait"),
  gwas_p_value = c("p"),
  gwas_z_score = c("zscore"),
  gwas_beta = c("beta"),
  gwas_se = c("se"),
  posterior_probability = c("postprob", "PP"),
  gene_id = c("Gene"),
  gene_symbol = c("SYMBOL"),
  gene_or_protein = c("SYMBOL", "Gene"),
  p_value = c("p"),
  distance_from_sentinel_bp = c("DISTANCE"),
  consequence = c("Consequence"),
  impact = c("IMPACT"),
  sift = c("SIFT"),
  polyphen = c("PolyPhen"),
  cadd_phred = c("CADD_PHRED"),
  evidence_detail = c("HGVSp", "HGVSc", "Feature")
)

out <- standardise_source_table(
  vep,
  mapping,
  constants = list(
    evidence_line = "credible_set_variant_annotation",
    evidence_source = "VEP_PP50_credible_set_variant_annotation",
    evidence_class = "variant_annotation",
    molecular_trait_type = "gene"
  ),
  source_file = source_file
)

write_standard_evidence(out, cfg, "credible_set_variant_annotation.tsv")

