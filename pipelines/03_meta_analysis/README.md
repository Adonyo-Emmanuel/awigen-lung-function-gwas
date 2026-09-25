# 03 Meta-Analysis

Meta-analysis of the per-group association results with METAL, producing the
`*.tbl` result files for each trait that are used by `04_fine_mapping/`,
`05_v2g/` and `06_figures/`.

> **Status:** scripts to be added.

## Expected contents

- METAL parameter files for FEV1, FVC, FEV1/FVC and PEF, including the scheme
  (for example `SCHEME STDERR`), genomic control settings, allele-frequency
  tracking and heterogeneity analysis.
- Pre-processing of per-group results into METAL input format.
- Post-processing of the METAL output: MAF filters, the minimum number of
  contributing groups, and genome-wide significance calls
  (P < 5 × 10⁻⁸).

## Output used downstream

`06_figures/awigen_manhattan_qq_plots.R` reads the columns `chromosome`,
`position`, `p` and `maf` from the per-trait METAL `.tbl` files.
