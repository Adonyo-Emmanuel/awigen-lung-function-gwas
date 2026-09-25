# Software Versions

Record the version of every tool used for the manuscript analyses. Fill in the
_to be added_ entries from the job logs or by running the tool's version
command on the cluster.

## External tools

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| 01_qc | PLINK | _to be added_ | Genotype QC |
| 01_qc | Imputation software / server and reference panel | _to be added_ | Genotype imputation |
| 02_association | Association software (e.g. REGENIE/SAIGE/BOLT-LMM) | _to be added_ | Per-group GWAS |
| 03_meta_analysis | METAL | _to be added_ | Meta-analysis |
| 04_fine_mapping | Fine-mapping software (e.g. GCTA-COJO, SuSiE) | _to be added_ | Conditional analysis and credible sets |
| 05_v2g | ANNOVAR | _to be added_ | Nearest-gene annotation |
| 05_v2g | Ensembl VEP | _to be added_ | Credible-set variant annotation |
| 05_v2g | UCSC liftOver | _to be added_ | GRCh37 to GRCh38 conversion |

## R

| Package | Used in | Version tested |
|---------|---------|----------------|
| R | all R scripts | _to be added_ |
| data.table | 05_v2g, 06_figures | 1.14.10 or later |
| ggplot2 | 06_figures | 3.4.4 or later |
| ggrepel | 06_figures | 0.9.5 or later |
| ragg (optional) | 06_figures | 1.2.7 or later |

The V2G pipeline's `manifest` step writes `R_session_info.txt` to its output
directory, recording the exact R and package versions used for that run.

To record the R environment for the figure scripts on the cluster, run:

```bash
Rscript -e 'library(data.table); library(ggplot2); library(ggrepel); sessionInfo()'
```

For full reproducibility you can also snapshot the R library with
[`renv`](https://rstudio.github.io/renv/) (`renv::init()`, then
`renv::snapshot()`) and commit the resulting `renv.lock`.
