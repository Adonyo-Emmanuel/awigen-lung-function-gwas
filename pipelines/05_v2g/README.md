# Variant-to-Gene and Variant-to-Protein Pipeline

This pipeline implements the AWI-Gen lung-function variant-to-gene and variant-to-protein analysis. It is organised by evidence line, with one module per scientific approach and a separate module for functional annotation across all evidence lines.

## Evidence-Line Modules

- `lines_of_evidence/nearest_gene_annotation/build_nearest_gene_evidence.R`
- `lines_of_evidence/credible_set_variant_annotation/build_credible_set_variant_annotation_evidence.R`
- `lines_of_evidence/eqtl_associations/build_eqtl_evidence.R`
- `lines_of_evidence/pqtl_associations/build_pqtl_evidence.R`
- `lines_of_evidence/rare_variant_associations/build_rare_variant_evidence.R`
- `lines_of_evidence/mendelian_respiratory_disease_genes/build_mendelian_respiratory_disease_evidence.R`
- `lines_of_evidence/mouse_knockout_respiratory_phenotypes/build_mouse_knockout_evidence.R`

Each script reads the relevant configured input file(s), converts source-specific column names to the shared V2G schema and writes one table to `standardised_line_evidence/` under the configured output directory.

## Final Integration

`functional_annotation_across_evidence_lines/build_all_source_v2g_tables.R` combines the seven standardised evidence tables and writes:

- `awigen_v2g_master_evidence_table.tsv`
- `awigen_v2g_source_summary.tsv`
- `awigen_v2g_gene_summary.tsv`
- `awigen_v2g_signal_summary.tsv`
- `awigen_v2g_gene_evidence_presence_table.tsv`
- `awigen_v2g_threshold_summary.tsv`

## Standard Columns

The canonical column definitions are listed in `docs/standard_evidence_schema.tsv`. Source-specific names such as `FDR`, `qval`, `sabr_fdr` and `eqtl_fdr` are converted to the standard column `fdr` inside the evidence-line scripts, before final integration.

## Run

```bash
cp pipelines/05_v2g/config/v2g_config.example.tsv pipelines/05_v2g/config/v2g_config.local.tsv
bash pipelines/05_v2g/bin/run_v2g_pipeline.sh --config pipelines/05_v2g/config/v2g_config.local.tsv --step all
```

Available steps are:

- `check`
- `manifest`
- `standardise-lines`
- `integrate`
- `all`

The local config is not committed to Git.

