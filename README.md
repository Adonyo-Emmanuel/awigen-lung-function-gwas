# AWI-Gen Lung Function GWAS

This repository contains analysis code for the AWI-Gen lung-function genome-wide association study, including the variant-to-gene and variant-to-protein analyses used to prioritise candidate genes, proteins and variants at associated loci.

The repository is organised around analysis functionality. Each analysis module has a defined input contract, produces tabular outputs with standardised column names and can be run from a local configuration file.

## Repository Structure

- `pipelines/v2g/`: variant-to-gene and variant-to-protein pipeline.
- `pipelines/v2g/lines_of_evidence/`: scripts for the seven V2G evidence lines.
- `pipelines/v2g/functional_annotation_across_evidence_lines/`: all-source integration and summary tables.
- `pipelines/v2g/config/`: configuration template for local file paths.
- `pipelines/v2g/docs/`: standard output schema and current script inventory.
- `pipelines/gwas_plots/`: publication Manhattan and QQ plots (ggplot2) for the meta-analysis results.
- `data/`: placeholder documenting data availability and access restrictions.
- `results/`: placeholder for generated outputs.

## V2G Evidence Lines

The V2G pipeline currently implements seven evidence lines:

1. Nearest-gene annotation.
2. Annotation of credible-set variants with posterior probability >50%.
3. eQTL associations.
4. pQTL associations.
5. Nearby rare-variant associations from exome sequencing.
6. Nearby Mendelian respiratory-disease genes.
7. Nearby mouse-knockout orthologs with respiratory phenotypes.

Each evidence-line script writes a standardised evidence table. The final integration script combines these tables into an all-source master table, gene-level summaries, signal-level summaries and an evidence-presence matrix.

## Running The V2G Pipeline

Create a local configuration file from the template:

```bash
cp pipelines/v2g/config/v2g_config.example.tsv pipelines/v2g/config/v2g_config.local.tsv
```

Edit `pipelines/v2g/config/v2g_config.local.tsv` so that each input key points to the corresponding local analysis file, then run:

```bash
bash pipelines/v2g/bin/run_v2g_pipeline.sh --config pipelines/v2g/config/v2g_config.local.tsv --step all
```

Local configuration files and generated results are ignored by Git.

## Data Availability

This repository does not distribute controlled-access cohort data, UK Biobank files, downloaded QTL resources or generated result tables. The code is intended to be run in an approved analysis environment where the required input files are available. Publicly redistributable lookup tables may be added when appropriate.

