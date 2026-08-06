#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_mouse_knockout_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
source_file <- config_value(cfg, "mouse_knockout_file")
mouse_knockout <- read_required_table(source_file)

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
  gene_symbol = c("Symbol"),
  gene_or_protein = c("Symbol"),
  molecular_trait_id = c("mouse_symbol", "MKO"),
  distance_from_sentinel_bp = c("distance"),
  evidence_detail = c("MKO", "mouse_symbol", "overlap")
)

out <- standardise_source_table(
  mouse_knockout,
  mapping,
  constants = list(
    evidence_line = "mouse_knockout_respiratory_phenotypes",
    evidence_source = "IMPC_respiratory_system_phenotype",
    evidence_class = "mouse_knockout",
    molecular_trait_type = "gene"
  ),
  source_file = source_file
)

write_standard_evidence(out, cfg, "mouse_knockout_respiratory_phenotypes.tsv")

