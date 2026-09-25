# 06 Figures

Scripts that produce the manuscript figures from the analysis outputs.

## Manhattan and QQ plots

`awigen_manhattan_qq_plots.R` draws the Manhattan and QQ plots for FEV1, FVC,
FEV1/FVC and PEF from the METAL meta-analysis results.

```bash
Rscript pipelines/06_figures/awigen_manhattan_qq_plots.R <metal_results_dir> [output_dir]
```

`<metal_results_dir>` must contain the per-trait METAL files named in
`trait_files` at the top of the script. Outputs are written to
`<output_dir>/manhattan/` and `<output_dir>/qq/` (default:
`<metal_results_dir>/publication_manhattan_qq_plots/`), as PNG and PDF.

- **Manhattan:** autosomal SNPs with MAF ≥ 1%. Chromosomes alternate blue and
  grey. The genome-wide significant loci are highlighted in red and labelled
  with their nearest gene. Dashed lines mark P = 5 × 10⁻⁸ and
  P = 5 × 10⁻⁶. SNPs with P > 0.05 are randomly thinned for plotting only.
- **QQ:** observed against expected −log10(P) for all tested autosomal SNPs,
  with the pointwise 95% confidence band and λGC.

R packages: `data.table`, `ggplot2`, `ggrepel`, and optionally `ragg` for
higher-quality PNG output.
