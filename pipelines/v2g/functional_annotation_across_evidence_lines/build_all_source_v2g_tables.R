#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript build_all_source_v2g_tables.R <config.tsv>")
}

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
pipeline_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(pipeline_dir, "R", "v2g_utils.R"))

cfg <- read_config(args[[1]])
out_dir <- pipeline_output_dir(cfg)
line_dir <- standard_evidence_dir(cfg)

expected_line_files <- c(
  "nearest_gene_annotation.tsv",
  "credible_set_variant_annotation.tsv",
  "eqtl_associations.tsv",
  "pqtl_associations.tsv",
  "rare_variant_associations.tsv",
  "mendelian_respiratory_disease_genes.tsv",
  "mouse_knockout_respiratory_phenotypes.tsv"
)

paths <- file.path(line_dir, expected_line_files)
missing <- paths[!file.exists(paths)]
if (length(missing) > 0) {
  stop(
    "Missing standardised evidence file(s). Run the line-of-evidence step first:\n",
    paste(missing, collapse = "\n")
  )
}

message("Reading standardised line-of-evidence files...")
master <- rbindlist(
  lapply(paths, function(path) {
    dt <- fread(path, check.names = FALSE)
    add_missing_standard_cols(dt)
    dt[, source_file := clean_text(source_file)]
    dt
  }),
  use.names = TRUE,
  fill = TRUE
)

for (col in intersect(numeric_evidence_columns(), names(master))) {
  master[, (col) := num(get(col))]
}

master[, gene_or_protein := clean_text(gene_or_protein)]
master[, gene_id := strip_ensembl_version(gene_id)]
master[, gene_symbol := clean_text(gene_symbol)]
master <- master[!is.na(sentinel_snp) & !is.na(gene_or_protein)]

if (nrow(master) == 0) {
  stop("No V2G evidence rows were available after standardisation.")
}

setorder(master, evidence_line, evidence_source, sentinel_snp, gene_or_protein)
write_tsv(master, file.path(out_dir, "awigen_v2g_master_evidence_table.tsv"))

source_summary <- master[, .(
  n_rows = .N,
  n_sentinels = uniqueN(sentinel_snp),
  n_genes_or_proteins = uniqueN(gene_or_protein)
), by = .(evidence_line, evidence_source, evidence_class, source_file)]
setorder(source_summary, evidence_line, evidence_source)
write_tsv(source_summary, file.path(out_dir, "awigen_v2g_source_summary.tsv"))

gene_summary <- master[, .(
  n_evidence_lines = uniqueN(evidence_line),
  n_evidence_sources = uniqueN(evidence_source),
  evidence_lines = join_unique(evidence_line),
  evidence_sources = join_unique(evidence_source),
  sentinels = join_unique(sentinel_snp),
  traits = join_unique(trait),
  gene_ids = join_unique(gene_id),
  gene_symbols = join_unique(gene_symbol),
  min_gwas_p_value = minimum_numeric(gwas_p_value),
  max_posterior_probability = maximum_numeric(posterior_probability),
  min_p_value = minimum_numeric(p_value),
  min_fdr = minimum_numeric(fdr),
  min_q_value = minimum_numeric(q_value),
  min_distance_from_sentinel_bp = minimum_numeric(distance_from_sentinel_bp),
  max_cadd_phred = maximum_numeric(cadd_phred)
), by = .(gene_or_protein)]
setorder(gene_summary, -n_evidence_lines, -n_evidence_sources, gene_or_protein)
write_tsv(gene_summary, file.path(out_dir, "awigen_v2g_gene_summary.tsv"))

signal_summary <- master[, .(
  n_evidence_lines = uniqueN(evidence_line),
  n_evidence_sources = uniqueN(evidence_source),
  evidence_lines = join_unique(evidence_line),
  evidence_sources = join_unique(evidence_source),
  genes_or_proteins = join_unique(gene_or_protein),
  n_genes_or_proteins = uniqueN(gene_or_protein),
  traits = join_unique(trait)
), by = .(sentinel_snp)]
setorder(signal_summary, sentinel_snp)
write_tsv(signal_summary, file.path(out_dir, "awigen_v2g_signal_summary.tsv"))

presence <- unique(master[, .(gene_or_protein, evidence_line)])
presence[, present := 1L]
presence_table <- dcast(
  presence,
  gene_or_protein ~ evidence_line,
  value.var = "present",
  fill = 0L
)
presence_table <- merge(
  gene_summary[, .(gene_or_protein, n_evidence_lines, n_evidence_sources)],
  presence_table,
  by = "gene_or_protein",
  all.x = TRUE,
  sort = FALSE
)
setorder(presence_table, -n_evidence_lines, -n_evidence_sources, gene_or_protein)
write_tsv(presence_table, file.path(out_dir, "awigen_v2g_gene_evidence_presence_table.tsv"))

threshold_summary <- data.table(
  metric = c(
    "master_rows",
    "unique_sentinels",
    "unique_genes_or_proteins",
    "genes_or_proteins_with_at_least_2_evidence_lines",
    "genes_or_proteins_with_at_least_3_evidence_lines",
    "sentinels_with_at_least_2_evidence_lines",
    "sentinels_with_at_least_3_evidence_lines"
  ),
  value = c(
    nrow(master),
    uniqueN(master$sentinel_snp),
    uniqueN(master$gene_or_protein),
    nrow(gene_summary[n_evidence_lines >= 2]),
    nrow(gene_summary[n_evidence_lines >= 3]),
    nrow(signal_summary[n_evidence_lines >= 2]),
    nrow(signal_summary[n_evidence_lines >= 3])
  )
)
write_tsv(threshold_summary, file.path(out_dir, "awigen_v2g_threshold_summary.tsv"))

message("All-source V2G tables written to: ", out_dir)

