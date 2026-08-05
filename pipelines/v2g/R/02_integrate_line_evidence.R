#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 02_integrate_line_evidence.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
source(file.path(script_dir, "v2g_utils.R"))

cfg <- read_config(args[[1]])
out_dir <- pipeline_output_dir(cfg)

evidence_sources <- data.table(
  evidence_source = c(
    "nearest_gene",
    "vep_pp50_annotation",
    "eqtlgen",
    "gtex",
    "sabr",
    "african_american_eqtl",
    "pqtl_olink",
    "pqtl_masc",
    "rare_variant",
    "rare_disease",
    "mouse_knockout"
  ),
  evidence_class = c(
    "nearest_gene",
    "credible_set_annotation",
    "eQTL",
    "eQTL",
    "eQTL",
    "eQTL",
    "pQTL",
    "pQTL",
    "rare_variant",
    "rare_disease",
    "mouse_knockout"
  ),
  config_key = c(
    "annovar_nearest_gene_file",
    "vep_pp50_gene_annotation_file",
    "eqtlgen_file",
    "gtex_file",
    "sabr_file",
    "amr_file",
    "pqtl_olink_file",
    "pqtl_masc_file",
    "rare_variant_file",
    "rare_disease_file",
    "mouse_knockout_file"
  )
)

read_evidence_file <- function(source_row) {
  path <- config_value(cfg, source_row$config_key)
  if (is.na(path) || !file.exists(path)) {
    stop("Missing evidence file for ", source_row$evidence_source, ": ", path)
  }

  dt <- fread(path, check.names = FALSE)
  if (nrow(dt) == 0) {
    return(data.table())
  }

  sentinel <- pull_col(dt, c(
    "sentinel_snp", "gwas_sentinel_snp", "CHR_POS_AMIN_AMAX", "chr_pos_amin_amax",
    "sentinel", "awigen_CHR_POS_AMIN_AMAX"
  ))
  snp_id <- pull_col(dt, c(
    "snp_id", "chr_pos_amin_amax", "gwas_chr_pos_amin_amax", "CHR_POS_AMIN_AMAX",
    "rare_variant_single_variant_id", "vep_uploaded_variation"
  ))
  trait <- pull_col(dt, c("trait", "gwas_trait", "TRAIT", "traits", "TRAITS", "rare_variant_single_variant_trait"))
  gene_symbol <- pull_col(dt, c(
    "gene_symbol", "GeneSymbol", "Gene Symbol", "sabr_s7_gene_symbol",
    "SYMBOL", "Symbol", "protein", "PROTEIN", "masc_protein",
    "collapsing_gene_symbol", "rare_variant_gene", "mouse_symbol",
    "gene", "Gene"
  ))
  gene_id <- pull_col(dt, c("gene_id", "Gene Ensembl ID", "sabr_gene_id", "Gene", "ENSG", "ensembl_gene_id"))
  protein <- pull_col(dt, c("protein", "PROTEIN", "masc_protein"))
  p_value <- pull_col(dt, c(
    "p_value", "p", "P", "gwas_p", "Pvalue", "P-Value",
    "eqtl_pvalue", "eqtl_pvalue_nominal", "sabr_pval_nominal", "pqtl_p",
    "masc_p", "rare_variant_single_variant_p", "rare_variant_collapsing_gene_p",
    "collapsing_gene_p"
  ), default = NA_real_)
  fdr <- pull_col(dt, c("fdr", "FDR", "eqtl_fdr", "sabr_fdr", "qval", "q_value"), default = NA_real_)
  postprob <- pull_col(dt, c("postprob", "gwas_postprob", "PP", "PPA"), default = NA_real_)
  distance <- pull_col(dt, c(
    "distance_to_gene", "distance", "dist", "distance_from_sentinel_bp",
    "rare_variant_distance_from_sentinel_bp"
  ), default = NA_real_)

  gene_symbol <- clean_text(gene_symbol)
  protein <- clean_text(protein)
  gene_symbol[is.na(gene_symbol) & !is.na(protein)] <- protein[is.na(gene_symbol) & !is.na(protein)]
  gene_id <- strip_ensembl_version(gene_id)

  out <- data.table(
    evidence_source = source_row$evidence_source,
    evidence_class = source_row$evidence_class,
    source_file = basename(path),
    sentinel_snp = clean_text(sentinel),
    snp_id = clean_text(snp_id),
    trait = clean_text(trait),
    gene_symbol = gene_symbol,
    gene_id = gene_id,
    p_value = num(p_value),
    fdr = num(fdr),
    postprob = num(postprob),
    distance_from_sentinel_bp = num(distance)
  )

  out[, gene_label := gene_symbol]
  out[is.na(gene_label) & !is.na(gene_id), gene_label := gene_id]
  out[!is.na(gene_label)]
}

message("Reading line-of-evidence files...")
master <- rbindlist(
  lapply(seq_len(nrow(evidence_sources)), function(i) read_evidence_file(evidence_sources[i])),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(master) == 0) {
  stop("No evidence rows were created. Check input files and gene/protein columns.")
}

setorder(master, evidence_class, evidence_source, sentinel_snp, gene_label)
write_tsv(master, file.path(out_dir, "awigen_v2g_master_evidence_table.tsv"))

source_summary <- master[, .(
  n_rows = .N,
  n_sentinels = uniqueN(sentinel_snp),
  n_genes_or_proteins = uniqueN(gene_label)
), by = .(evidence_source, evidence_class, source_file)]
setorder(source_summary, evidence_class, evidence_source)
write_tsv(source_summary, file.path(out_dir, "awigen_v2g_source_summary.tsv"))

gene_summary <- master[, .(
  n_evidence_sources = uniqueN(evidence_source),
  n_evidence_classes = uniqueN(evidence_class),
  evidence_sources = join_unique(evidence_source),
  evidence_classes = join_unique(evidence_class),
  sentinels = join_unique(sentinel_snp),
  traits = join_unique(trait),
  min_p_value = minimum_numeric(p_value),
  min_fdr = minimum_numeric(fdr),
  max_postprob = maximum_numeric(postprob),
  min_distance_from_sentinel_bp = minimum_numeric(distance_from_sentinel_bp)
), by = .(gene_label)]
setorder(gene_summary, -n_evidence_classes, -n_evidence_sources, gene_label)
write_tsv(gene_summary, file.path(out_dir, "awigen_v2g_gene_summary.tsv"))

signal_summary <- master[, .(
  n_evidence_sources = uniqueN(evidence_source),
  n_evidence_classes = uniqueN(evidence_class),
  evidence_sources = join_unique(evidence_source),
  evidence_classes = join_unique(evidence_class),
  genes_or_proteins = join_unique(gene_label),
  n_genes_or_proteins = uniqueN(gene_label)
), by = .(sentinel_snp)]
setorder(signal_summary, sentinel_snp)
write_tsv(signal_summary, file.path(out_dir, "awigen_v2g_signal_summary.tsv"))

threshold_summary <- data.table(
  metric = c(
    "master_rows",
    "unique_genes_or_proteins",
    "genes_or_proteins_with_at_least_2_evidence_classes",
    "genes_or_proteins_with_at_least_3_evidence_classes",
    "sentinels_with_at_least_2_evidence_classes",
    "sentinels_with_at_least_3_evidence_classes"
  ),
  value = c(
    nrow(master),
    uniqueN(master$gene_label),
    nrow(gene_summary[n_evidence_classes >= 2]),
    nrow(gene_summary[n_evidence_classes >= 3]),
    nrow(signal_summary[n_evidence_classes >= 2]),
    nrow(signal_summary[n_evidence_classes >= 3])
  )
)
write_tsv(threshold_summary, file.path(out_dir, "awigen_v2g_threshold_summary.tsv"))

message("Integrated V2G outputs written to: ", out_dir)
