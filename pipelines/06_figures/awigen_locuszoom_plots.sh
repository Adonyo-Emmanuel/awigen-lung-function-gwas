#!/bin/bash
# Regional association plots for the AWI-Gen genome-wide significant loci,
# drawn with LocusZoom 1.4 standalone.
#
# LD: 1000 Genomes Phase 3 (1000G_Nov2014), AFR super-population, GRCh37.
# Lines: red = genome-wide significance (P = 5e-8),
#        grey = suggestive significance (P = 5e-6).
#
# Usage:
#   bash awigen_locuszoom_plots.sh <metal_results_dir> [output_dir]
#
# <metal_results_dir> must contain <trait>_locuszoom_formated.txt, made by
# format_locuszoom_input.R. The loci are listed in top4_sentinels_awigen.tsv
# (next to this script). LocusZoom must be on PATH.
#
# LD is looked up by chr:pos; the lead variant is labelled with its rsID.
# Each plot is written as both PDF and PNG.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <metal_results_dir> [output_dir]" >&2
  exit 1
fi

metal_dir="$1"
output_dir="${2:-${metal_dir}/awigen_regional_assoc_plots}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sentinel_file="${script_dir}/top4_sentinels_awigen.tsv"

command -v locuszoom >/dev/null || { echo "locuszoom is not on PATH" >&2; exit 1; }
mkdir -p "${output_dir}"

trait_label() {
  case "$1" in
    FEV1) echo "FEV1" ;;
    FVC)  echo "FVC" ;;
    FF)   echo "FEV1/FVC" ;;
    PEF)  echo "PEF" ;;
    *)    echo "$1" ;;
  esac
}

# Columns: chrpos, rsid, trait, gene (tab- or space-separated; header skipped).
tail -n +2 "${sentinel_file}" | tr -d '\r' | while read -r chrpos rsid trait gene; do
  [[ -z "${chrpos}" ]] && continue

  metal_file="${metal_dir}/${trait}_locuszoom_formated.txt"
  [[ -f "${metal_file}" ]] || { echo "Missing ${metal_file}" >&2; exit 1; }

  out_prefix="${output_dir}/awigen_${trait}_${gene}"
  echo "Plotting ${gene} (${trait}, ${rsid}, ${chrpos})"

  locuszoom --metal "${metal_file}" --markercol rsid --pvalcol p \
    --refsnp "${chrpos}" --flank 500kb \
    --build hg19 --pop AFR --source 1000G_Nov2014 \
    --plotonly --snpset NULL --no-date --prefix "${out_prefix}" \
    title="$(trait_label "${trait}"): ${gene} locus" \
    refsnpName="${rsid}" \
    format=both \
    theme="publication" \
    showPartialGenes=TRUE \
    geneFontSize=0.7 \
    rfrows=8 \
    refsnpTextSize=0.8 \
    signifLine="7.30103,5.30103" \
    signifLineColor="red,grey50" \
    signifLineWidth=1 \
    ldColors="gray60,navy,lightskyblue,green,orange,red,purple3" \
    recombColor="blue" \
    recombAxisColor="black"
done

echo "Plots written to ${output_dir}"
