#!/bin/bash
# Regional association plots for the AWI-Gen genome-wide significant loci,
# drawn with LocusZoom 1.4 standalone.
#
# LD: 1000 Genomes Phase 3, African (AFR) samples, GRCh37. Run
# prepare_1000g_phase3_afr_ld.sh first; its VCF
# (<metal_results_dir>/1000G_phase3_AFR_ld/1000G_phase3_AFR_loci.vcf.gz) is
# used automatically. Without it, LocusZoom's built-in reference named by
# LD_SOURCE is used (default 1000G_Nov2014).
# Lines: red solid = genome-wide significance (P = 5e-8),
#        grey dashed = suggestive significance (P = 5e-6).
#
# Usage:
#   bash awigen_locuszoom_plots.sh <metal_results_dir> [output_dir]
#
# <metal_results_dir> must contain <trait>_locuszoom_formated.txt, made by
# format_locuszoom_input.R. The loci are listed in top4_sentinels_awigen.tsv
# (next to this script). LocusZoom and PLINK must be on PATH, e.g.
#   module load R/4.3.1 plink
#   export PATH=${PATH}:<locuszoom_install>/bin
#
# LD is looked up by chr:pos; the lead variant is labelled with its rsID.
# LocusZoom writes a PDF (page 1 = plot, page 2 = settings log); page 1 is
# also converted to a 300 dpi PNG if pdftoppm or Ghostscript is available.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <metal_results_dir> [output_dir]" >&2
  exit 1
fi

metal_dir="$1"
output_dir="${2:-${metal_dir}/awigen_regional_assoc_plots}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sentinel_file="${script_dir}/top4_sentinels_awigen.tsv"
ld_vcf="${metal_dir}/1000G_phase3_AFR_ld/1000G_phase3_AFR_loci.vcf.gz"
ld_source="${LD_SOURCE:-1000G_Nov2014}"

command -v locuszoom >/dev/null || { echo "locuszoom is not on PATH" >&2; exit 1; }

# LocusZoom 1.4 is written in Python 2 and otherwise runs whichever "python"
# is first on PATH (modules or conda may put Python 3 there). Find Python 2.
lz_python=""
for candidate in "${LZ_PYTHON:-}" python2.7 python2 python /usr/bin/python2.7 /usr/bin/python2 /usr/bin/python; do
  [[ -n "${candidate}" ]] || continue
  if command -v "${candidate}" >/dev/null 2>&1 &&
     "${candidate}" -c 'import sys; sys.exit(0 if sys.version_info[0] == 2 else 1)' 2>/dev/null; then
    lz_python="$(command -v "${candidate}")"
    break
  fi
done
[[ -n "${lz_python}" ]] || {
  echo "LocusZoom needs Python 2.7, but none was found. Set LZ_PYTHON=/path/to/python2.7" >&2
  exit 1
}
echo "Python for LocusZoom: ${lz_python}"
# LocusZoom calls PLINK to compute LD from the 1000 Genomes reference.
command -v plink >/dev/null || { echo "plink is not on PATH (e.g. run: module load plink)" >&2; exit 1; }
mkdir -p "${output_dir}"
if [[ -f "${ld_vcf}" ]]; then
  command -v tabix >/dev/null || { echo "tabix is not on PATH (e.g. module load tabix)" >&2; exit 1; }
  ld_args=(--ld-vcf "${ld_vcf}")
  echo "LD reference: 1000 Genomes Phase 3 AFR VCF (${ld_vcf})"
else
  ld_args=(--pop AFR --source "${ld_source}")
  echo "LD reference: LocusZoom built-in ${ld_source} AFR"
fi

trait_label() {
  case "$1" in
    FEV1) echo "FEV1" ;;
    FVC)  echo "FVC" ;;
    FF)   echo "FEV1/FVC" ;;
    PEF)  echo "PEF" ;;
    *)    echo "$1" ;;
  esac
}

mkdir -p "${output_dir}/panels"

# Run LocusZoom for one locus and convert page 1 of the PDF to PNG.
#   $1 = output prefix; remaining arguments are extra plot settings (e.g. title).
plot_locus() {
  local prefix="$1"; shift
  "${lz_python}" "$(command -v locuszoom)" --metal "${metal_file}" --markercol rsid --pvalcol p \
    --refsnp "${chrpos}" --flank 1000kb \
    --build hg19 "${ld_args[@]}" \
    --plotonly --snpset NULL --no-date --prefix "${prefix}" \
    "$@" \
    refsnpName="${rsid}" \
    theme="publication" \
    showPartialGenes=TRUE \
    geneFontSize=0.7 \
    rfrows=8 \
    refsnpTextSize=0.8 \
    signifLine="7.30103,5.30103" \
    signifLineColor="red,grey50" \
    signifLineWidth="1,1" \
    signifLineType="1,2" \
    ldColors="gray60,navy,lightskyblue,green,orange,red,purple3" \
    recombColor="blue" \
    recombAxisColor="black"

  local pdf png
  for pdf in "${prefix}"_*.pdf; do
    [[ -f "${pdf}" ]] || { echo "No PDF was produced for ${prefix}" >&2; exit 1; }
    png="${pdf%.pdf}.png"
    if command -v pdftoppm >/dev/null; then
      pdftoppm -png -r 300 -f 1 -l 1 -singlefile "${pdf}" "${png%.png}"
    elif command -v gs >/dev/null; then
      gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
        -dFirstPage=1 -dLastPage=1 -sOutputFile="${png}" "${pdf}"
    else
      echo "  (no pdftoppm or Ghostscript found; PNG not made)"
      continue
    fi
    echo "  wrote ${pdf} and ${png}"
  done
}

# Columns: chrpos, rsid, trait, gene, alleles (tab- or space-separated; header skipped).
tail -n +2 "${sentinel_file}" | tr -d '\r' | while read -r chrpos rsid trait gene _; do
  [[ -z "${chrpos}" ]] && continue

  metal_file="${metal_dir}/${trait}_locuszoom_formated.txt"
  [[ -f "${metal_file}" ]] || { echo "Missing ${metal_file}" >&2; exit 1; }

  echo "Plotting ${gene} (${trait}, ${rsid}, ${chrpos})"
  # Titled plot for individual use.
  plot_locus "${output_dir}/awigen_${trait}_${gene}" title="$(trait_label "${trait}"): ${gene} locus"
  # Untitled copy for the combined multi-panel figure (panel letters added later).
  plot_locus "${output_dir}/panels/awigen_${trait}_${gene}"
done

echo "Plots written to ${output_dir}"
echo "Untitled panels for the combined figure: ${output_dir}/panels"
