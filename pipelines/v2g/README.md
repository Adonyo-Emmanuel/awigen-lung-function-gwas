# AWI-Gen V2G Pipeline

This folder contains a reproducible pipeline scaffold for the AWI-Gen variant-to-gene and variant-to-protein analysis.

The pipeline is designed to be committed to the AWI-Gen GWAS GitHub repository while keeping private data paths in a local config file that does not need to be committed.

## What This Pipeline Does

The current version performs four reproducibility tasks:

1. Checks that all configured V2G input files are present.
2. Creates an audit manifest with file sizes, modification times, MD5 checksums, R session information and git status.
3. Integrates the final line-of-evidence files into a master V2G evidence table and gene/sentinel summaries.
4. Provides a wrapper for the current authoritative finalisation scripts already present under `QTL_analyses/functional_annotations_across_all_lines_of_evidences`.

The line-of-evidence inputs correspond to:

- nearest-gene annotation
- credible-set variant annotation with PP >50%
- eQTLGen, GTEx, SABR and African American cis-eQTL resources
- UK Biobank Olink and MASC pQTL resources
- nearby rare-variant associations from exome sequencing
- nearby Mendelian respiratory-disease genes
- nearby mouse-knockout orthologs with respiratory phenotypes

## How To Run

Copy the example config and edit the paths:

```bash
cp pipelines/v2g/config/v2g_config.example.tsv pipelines/v2g/config/v2g_config.local.tsv
```

Then run:

```bash
bash pipelines/v2g/bin/run_v2g_pipeline.sh --config pipelines/v2g/config/v2g_config.local.tsv --step all
```

Outputs are written to the `output_dir` specified in the config.

For the current local `QTL_analyses` tree, create a local config from the example file and set the input paths to the existing analysis files.

To reproduce the current authoritative final tables from the existing final scripts, set `QTL_ROOT` to the local `QTL_analyses` directory:

```bash
QTL_ROOT=/path/to/QTL_analyses bash pipelines/v2g/bin/run_current_qtl_v2g_finalisation.sh
```

## Recommended GitHub Layout

For the full AWI-Gen GWAS repository, I recommend:

```text
awigen-gwas/
  README.md
  data/
    README.md
  scripts/
    gwas/
    fine_mapping/
    v2g/
  pipelines/
    v2g/
  results/
    README.md
  docs/
    methods/
```

Large cohort data, UK Biobank files and intermediate QTL resources should not be committed. Instead, commit the pipeline, config templates, small lookup tables that are public or redistributable, and a manifest describing the exact input versions used.

The current script inventory is in:

```text
pipelines/v2g/docs/current_v2g_script_inventory.tsv
```

## Notes

This scaffold intentionally separates the reproducible integration pipeline from older exploratory scripts. The historical scripts can be moved into `scripts/v2g/` and gradually refactored so each stage accepts a config path and writes into a controlled output directory.
