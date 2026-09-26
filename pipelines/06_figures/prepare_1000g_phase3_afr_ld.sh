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
: > bcftools.log

# bcftools messages go to bcftools.log; show them if any step fails.
trap 'status=$?; if [[ ${status} -ne 0 ]]; then
  echo "ERROR: a step failed (exit ${status}). Last bcftools messages:" >&2
  tail -20 "${out_dir}/bcftools.log" >&2
fi' EXIT

echo "Downloading the 1000 Genomes Phase 3 sample list.."
curl -sSfL -o samples.panel "${base_url}/integrated_call_samples_v3.20130502.ALL.panel"
awk '$3 == "AFR" {print $1}' samples.panel > afr_samples.txt
echo "AFR samples: $(wc -l < afr_samples.txt)"

# One region per sentinel, sorted by chromosome so the parts concatenate in order.
tail -n +2 "${sentinel_file}" | tr -d '\r' | awk -v f="${flank}" 'NF {
  split($1, a, ":"); chr = a[1]; sub(/^chr/, "", chr); pos = a[2];
  start = pos - f; if (start < 1) start = 1;
  print chr "\t" start "\t" pos + f "\t" pos "\t" $2 "\t" $5
}' | sort -k1,1n -k2,2n > regions.tsv

# LocusZoom matches variants by chr:pos only and uses the first VCF record at a
# position, so each position must have exactly one record:
#   - multi-allelic sites are split into biallelic records;
#   - at the sentinel position, the record matching the sentinel rsID and
#     alleles is kept (else the rsID alone, else the first SNP);
#   - at other positions a SNP is preferred over an indel;
#   - indels that span the sentinel position (e.g. an upstream deletion) are
#     dropped, because LocusZoom's lookup of the sentinel would return them too.
parts=()
while IFS=$'\t' read -r chr start end pos rsid alleles; do
  a1="${alleles%/*}"; a2="${alleles#*/}"
  tag="chr${chr}_${pos}"
  echo "Fetching chr${chr}:${start}-${end} (${rsid}); this can take a few minutes.."
  bcftools view -r "${chr}:${start}-${end}" -S afr_samples.txt -Ou \
    "${base_url}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" |
    bcftools norm -m -any -Ou 2>>bcftools.log |
    bcftools view -c 1 -Oz -o "${tag}.split.vcf.gz"
  tabix -f -p vcf "${tag}.split.vcf.gz"

  # Sentinel record: rsID and alleles, else rsID, else the first SNP there.
  for filter in \
    "ID==\"${rsid}\" && ((REF==\"${a1}\" && ALT==\"${a2}\") || (REF==\"${a2}\" && ALT==\"${a1}\"))" \
    "ID==\"${rsid}\"" \
    "TYPE==\"snp\""; do
    bcftools view -r "${chr}:${pos}" -i "${filter}" -Ou "${tag}.split.vcf.gz" |
      bcftools norm -d all -Oz -o "${tag}.index.vcf.gz" 2>>bcftools.log
    [[ -n "$(bcftools view -H "${tag}.index.vcf.gz" | head -1)" ]] && break
  done
  if [[ -z "$(bcftools view -H "${tag}.index.vcf.gz" | head -1)" ]]; then
    echo "  WARNING: ${rsid} (chr${chr}:${pos}) is not a polymorphic variant in the" \
         "1000G Phase 3 AFR samples; LD cannot be shown for this locus." >&2
  else
    echo "  sentinel found: $(bcftools view -H "${tag}.index.vcf.gz" | cut -f1-5 | head -1)"
  fi

  # Other positions: SNPs first, then indels only where no SNP exists.
  bcftools view -e "POS==${pos}" -v snps "${tag}.split.vcf.gz" -Ou |
    bcftools norm -d all -Oz -o "${tag}.snps.vcf.gz" 2>>bcftools.log
  bcftools query -f '%CHROM\t%POS\n' "${tag}.snps.vcf.gz" > "${tag}.snp_positions.txt"
  echo -e "${chr}\t${pos}" >> "${tag}.snp_positions.txt"
  bcftools view -V snps -T "^${tag}.snp_positions.txt" "${tag}.split.vcf.gz" -Ou |
    bcftools view -e "POS<=${pos} && POS+strlen(REF)-1>=${pos}" -Ou |
    bcftools norm -d all -Oz -o "${tag}.other.vcf.gz" 2>>bcftools.log

  bcftools concat "${tag}.index.vcf.gz" "${tag}.snps.vcf.gz" "${tag}.other.vcf.gz" -Ou 2>>bcftools.log |
    bcftools sort -Oz -o "${tag}.vcf.gz" 2>>bcftools.log
  rm -f "${tag}".split.vcf.gz* "${tag}".index.vcf.gz "${tag}".snps.vcf.gz "${tag}".other.vcf.gz "${tag}".snp_positions.txt
  parts+=("${tag}.vcf.gz")
done < regions.tsv

bcftools concat -Oz -o "${out_vcf}" "${parts[@]}" 2>>bcftools.log
tabix -f -p vcf "${out_vcf}"

# LocusZoom looks up each sentinel with tabix and needs exactly one record back.
while IFS=$'\t' read -r chr start end pos rsid alleles; do
  n=$(tabix "${out_vcf}" "${chr}:${pos}-${pos}" | wc -l)
  if [[ "${n}" -gt 1 ]]; then
    echo "ERROR: ${n} records overlap ${rsid} (chr${chr}:${pos}); LocusZoom cannot use them:" >&2
    tabix "${out_vcf}" "${chr}:${pos}-${pos}" | cut -f1-5 >&2
    exit 1
  fi
done < regions.tsv
rm -f "${parts[@]}" ./*.tbi.* 2>/dev/null || true
rm -f ALL.chr*.vcf.gz.tbi

echo "LD reference written to ${out_vcf}"
