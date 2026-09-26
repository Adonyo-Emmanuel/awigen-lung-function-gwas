#!/bin/bash
# Download 1000 Genomes Phase 3 genotypes for the African (AFR) samples in the
# regions around the AWI-Gen sentinel variants, for LocusZoom's --ld-vcf
# option. Only the regions are fetched (remote tabix queries), so the output
# is a few MB.
#
# Usage:
#   bash prepare_1000g_phase3_afr_ld.sh <metal_results_dir>
#
# Writes <metal_results_dir>/1000G_phase3_AFR_ld/1000G_phase3_AFR_loci.vcf.gz
# (+ .tbi), which awigen_locuszoom_plots.sh uses automatically.
#
# Needs internet access (run on a login node) and bcftools and tabix on PATH.
# Source: 1000 Genomes Project Phase 3 release 20130502 (GRCh37).

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <metal_results_dir>" >&2
  exit 1
fi

for tool in bcftools tabix; do
  command -v "${tool}" >/dev/null || { echo "${tool} is not on PATH (e.g. module load ${tool})" >&2; exit 1; }
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sentinel_file="${script_dir}/top4_sentinels_awigen.tsv"
out_dir="$(cd "$1" && pwd)/1000G_phase3_AFR_ld"
out_vcf="${out_dir}/1000G_phase3_AFR_loci.vcf.gz"
base_url="${KG_BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502}"
flank=600000

mkdir -p "${out_dir}"
cd "${out_dir}"

echo "Downloading the 1000 Genomes Phase 3 sample list.."
curl -sSfL -o samples.panel "${base_url}/integrated_call_samples_v3.20130502.ALL.panel"
awk '$3 == "AFR" {print $1}' samples.panel > afr_samples.txt
echo "AFR samples: $(wc -l < afr_samples.txt)"

# One region per sentinel, sorted by chromosome so the parts concatenate in order.
tail -n +2 "${sentinel_file}" | tr -d '\r' | awk -v f="${flank}" 'NF {
  split($1, a, ":"); chr = a[1]; sub(/^chr/, "", chr); pos = a[2];
  start = pos - f; if (start < 1) start = 1;
  print chr "\t" start "\t" pos + f "\t" pos "\t" $2
}' | sort -k1,1n -k2,2n > regions.tsv

parts=()
while IFS=$'\t' read -r chr start end pos rsid; do
  part="region_chr${chr}_${pos}.vcf.gz"
  echo "Fetching chr${chr}:${start}-${end} (${rsid})"
  bcftools view -r "${chr}:${start}-${end}" -S afr_samples.txt -Ou \
    "${base_url}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" |
    bcftools view -m2 -M2 -c 1 -Oz -o "${part}"
  parts+=("${part}")

  if [[ -z "$(bcftools view -H "${part}" 2>/dev/null | awk -v p="${pos}" '$2 == p' | head -1)" ]]; then
    echo "  WARNING: ${rsid} (chr${chr}:${pos}) is not in the 1000G Phase 3 AFR data;" \
         "LD cannot be shown for this locus." >&2
  fi
done < regions.tsv

bcftools concat -Oz -o "${out_vcf}" "${parts[@]}"
tabix -f -p vcf "${out_vcf}"
rm -f "${parts[@]}" ./*.tbi.* 2>/dev/null || true
rm -f ALL.chr*.vcf.gz.tbi

echo "LD reference written to ${out_vcf}"
