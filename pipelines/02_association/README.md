# 02 Association Analysis

Genome-wide association testing of FEV1, FVC, FEV1/FVC and PEF within each
AWI-Gen study site or analysis group.

> **Status:** scripts to be added.

## Expected contents

- Association software commands (for example REGENIE, SAIGE or BOLT-LMM),
  including step 1 (null model) and step 2 (single-variant tests).
- Covariate definitions: age, age², sex, height, smoking status, principal
  components and any site-specific covariates.
- The job-submission wrappers (for example SLURM scripts) used on the cluster,
  with site-specific paths moved into a local config file.
- Per-group output formatting required as input to `03_meta_analysis/`.

## Conventions

- Pin the software version used, or record it in the output logs.
- Per-group summary statistics are generated outputs and are not committed.
