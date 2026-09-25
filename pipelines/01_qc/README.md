# 01 Quality Control

Sample and variant quality control of the AWI-Gen genotype data, and the
post-imputation filters applied before association testing.

> **Status:** scripts to be added.

## Expected contents

- Sample QC: call rate, sex check, heterozygosity, relatedness, ancestry /
  principal components, phenotype exclusions.
- Variant QC: call rate, Hardy-Weinberg equilibrium, MAF, duplicate and
  strand checks.
- Imputation: reference panel, imputation server or software, and the
  post-imputation filters (for example INFO/R² and MAF thresholds).
- Lung-function phenotype preparation: spirometry quality criteria, outlier
  handling and trait transformations (for example rank-based inverse normal
  transformation of residuals).

## Conventions

- Scripts take input and output locations as command-line arguments or a
  `*.local.tsv` config file. Do not hard-code cluster paths.
- Do not commit genotype files, sample lists, participant IDs or
  individual-level phenotypes. `.gitignore` blocks the common formats.
- Record every threshold in the script or config file so that the numbers in
  the Methods section can be traced to code.
