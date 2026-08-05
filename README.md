# AWI-Gen Lung Function GWAS

This repository contains reproducible code for the AWI-Gen lung-function genome-wide association study.

The current repository starts with the variant-to-gene and variant-to-protein pipeline used to prioritise putative candidate genes and variants from AWI-Gen lung-function signals.

## Repository Contents

- `pipelines/v2g/`: R and Bash pipeline for the AWI-Gen variant-to-gene and variant-to-protein analyses.
- `data/`: placeholder for local input data. Restricted cohort data and large external resources should not be committed.
- `results/`: placeholder for generated outputs. Pipeline results should be regenerated from code and local configuration.

## V2G Pipeline

Copy the example V2G configuration, edit paths for the local analysis environment, and run:

```bash
cp pipelines/v2g/config/v2g_config.example.tsv pipelines/v2g/config/v2g_config.local.tsv
bash pipelines/v2g/bin/run_v2g_pipeline.sh --config pipelines/v2g/config/v2g_config.local.tsv --step all
```

The local configuration file is ignored by Git so that machine-specific paths do not enter the repository.

## Data Policy

Do not commit individual-level cohort data, UK Biobank resources, downloaded QTL summary statistics, generated intermediate files or final result tables. Commit code, configuration templates, small public lookup tables when redistributable, documentation and reproducibility manifests.

