# AWI-Gen Lung Function GWAS

Analysis code for the AWI-Gen genome-wide association study of lung function
(FEV1, FVC, FEV1/FVC and PEF). It covers quality control, association
testing, meta-analysis, fine-mapping, variant-to-gene (V2G) prioritisation and
the manuscript figures.

> **Manuscript:** _citation and DOI to be added on publication._

## Repository Structure

The analysis is organised as numbered steps. Each step's outputs feed the next.

| Step | Folder | Contents |
|------|--------|----------|
| 1 | [`pipelines/01_qc/`](pipelines/01_qc/) | Sample and variant QC, imputation filters, phenotype preparation |
| 2 | [`pipelines/02_association/`](pipelines/02_association/) | Per-group GWAS of the four lung-function traits |
| 3 | [`pipelines/03_meta_analysis/`](pipelines/03_meta_analysis/) | METAL meta-analysis and post-processing |
| 4 | [`pipelines/04_fine_mapping/`](pipelines/04_fine_mapping/) | Independent signals and credible sets |
| 5 | [`pipelines/05_v2g/`](pipelines/05_v2g/) | Variant-to-gene and variant-to-protein prioritisation (seven evidence lines) |
| 6 | [`pipelines/06_figures/`](pipelines/06_figures/) | Manuscript figures |

Supporting folders:

- [`docs/`](docs/): software versions and other documentation.
- [`data/`](data/): data availability and access (no data are distributed).
- [`results/`](results/): placeholder for generated outputs (not version-controlled).

## Manuscript Figures and Tables

| Manuscript item | Description | Script |
|-----------------|-------------|--------|
| Figure _X_ | Manhattan plots for FEV1, FVC, FEV1/FVC and PEF | [`pipelines/06_figures/awigen_manhattan_qq_plots.R`](pipelines/06_figures/awigen_manhattan_qq_plots.R) |
| Supplementary Figure _X_ | QQ plots with λGC | [`pipelines/06_figures/awigen_manhattan_qq_plots.R`](pipelines/06_figures/awigen_manhattan_qq_plots.R) |
| Table _X_ | V2G gene summary | [`pipelines/05_v2g/functional_annotation_across_evidence_lines/build_all_source_v2g_tables.R`](pipelines/05_v2g/functional_annotation_across_evidence_lines/build_all_source_v2g_tables.R) |
| Supplementary Table _X_ | V2G master evidence table | [`pipelines/05_v2g/functional_annotation_across_evidence_lines/build_all_source_v2g_tables.R`](pipelines/05_v2g/functional_annotation_across_evidence_lines/build_all_source_v2g_tables.R) |

_Replace X with the final figure and table numbers once the manuscript is accepted._

## Running The Analyses

Scripts take their input locations from command-line arguments or from local
configuration files (`*.local.tsv`), which Git ignores. No cluster-specific
paths are stored in the repository.

### V2G pipeline

```bash
cp pipelines/05_v2g/config/v2g_config.example.tsv pipelines/05_v2g/config/v2g_config.local.tsv
# edit v2g_config.local.tsv so each key points to your local input file
bash pipelines/05_v2g/bin/run_v2g_pipeline.sh --config pipelines/05_v2g/config/v2g_config.local.tsv --step all
```

### Manhattan and QQ plots

```bash
Rscript pipelines/06_figures/awigen_manhattan_qq_plots.R <metal_results_dir> [output_dir]
```

See each folder's README for details.

## Software

Software and package versions are listed in [`docs/software_versions.md`](docs/software_versions.md).

## Data Availability

This repository contains code only. See [`data/README.md`](data/README.md) for
how to access the AWI-Gen data and the GWAS summary statistics.

## Citation

If you use this code, please cite the manuscript (see
[`CITATION.cff`](CITATION.cff), or use GitHub's "Cite this repository" button).

## License

The code is released under the MIT License; see [`LICENSE`](LICENSE).
