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

## Regional association plots

LocusZoom 1.4 plots of the four genome-wide significant loci listed in
`top4_sentinels_awigen.tsv`, using ±500 kb windows and LD from 1000 Genomes
Phase 3 AFR (GRCh37). The lead variant is labelled with its rsID
and the plot title names the trait and nearest gene. Each plot is written as
PDF and PNG.

```bash
export PATH=${PATH}:<locuszoom_install>/bin
Rscript pipelines/06_figures/format_locuszoom_input.R <metal_results_dir>
bash pipelines/06_figures/awigen_locuszoom_plots.sh <metal_results_dir> [output_dir]
```

The formatting step keeps variants with MAF ≥ 1% (matching the Manhattan
plots) and keeps one variant per position (the one with the smallest P value).
