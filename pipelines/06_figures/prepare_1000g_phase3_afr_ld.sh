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
# Needs internet access (run on a login node) and bcftools, tabix and bgzip
# on PATH (bgzip and tabix come with htslib/tabix).
# Source: 1000 Genomes Project Phase 3 release 20130502 (GRCh37).

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 <metal_results_dir>" >&2
  exit 1
fi

for tool in bcftools tabix bgzip; do
  command -v "${tool}" >/dev/null || { echo "${tool} is not on PATH (e.g. module load ${tool})" >&2; exit 1; }
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sentinel_file="${script_dir}/top4_sentinels_awigen.tsv"
out_dir="$(cd "$1" && pwd)/1000G_phase3_AFR_ld"
out_vcf="${out_dir}/1000G_phase3_AFR_loci.vcf.gz"
base_url="${KG_BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502}"
flank=1100000

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
#   - structural variants (ALT like <CN0>) are dropped: they can span the
#     sentinel and LocusZoom cannot compute r2 for them;
#   - multi-allelic sites are split into biallelic records;
#   - at the sentinel position, the record matching the sentinel rsID and
#     alleles is kept (else the rsID alone, else the first SNP);
#   - at other positions a SNP is preferred over an indel;
#   - indels that span the sentinel position (e.g. an upstream deletion) are
#     dropped, because LocusZoom's lookup of the sentinel would return them too.
body_files=()
while IFS=$'\t' read -r chr start end pos rsid alleles; do
  a1="${alleles%/*}"; a2="${alleles#*/}"
  tag="chr${chr}_${pos}"
  echo "Fetching chr${chr}:${start}-${end} (${rsid}); this can take a few minutes.."
  bcftools view -r "${chr}:${start}-${end}" -S afr_samples.txt -Ou \
    "${base_url}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
    2>>bcftools.log |
    bcftools norm -m -any -Ou 2>>bcftools.log |
    bcftools view -c 1 2>>bcftools.log |
    awk -F'\t' '/^#/ || $5 !~ /^</' > "${tag}.split.vcf"   # drop structural variants (<CN0>, <DEL>, ...)

  [[ -f header.vcf ]] || grep '^#' "${tag}.split.vcf" > header.vcf

  # Sentinel record: rsID and alleles, else rsID, else the first SNP there.
  grep -v '^#' "${tag}.split.vcf" | awk -F'\t' -v p="${pos}" -v rs="${rsid}" -v a1="${a1}" -v a2="${a2}" '
    $2 == p {
      snp = (length($4) == 1 && length($5) == 1)
      if ($3 == rs && (($4 == a1 && $5 == a2) || ($4 == a2 && $5 == a1))) { if (!m1) m1 = $0 }
      else if ($3 == rs) { if (!m2) m2 = $0 }
      else if (snp)      { if (!m3) m3 = $0 }
    }
    END { if (m1) print m1; else if (m2) print m2; else if (m3) print m3 }' > "${tag}.index.txt"

  if [[ -s "${tag}.index.txt" ]]; then
    echo "  sentinel found: $(cut -f1-5 "${tag}.index.txt")"
  else
    echo "  WARNING: ${rsid} (chr${chr}:${pos}) is not a polymorphic variant in the" \
         "1000G Phase 3 AFR samples; LD cannot be shown for this locus." >&2
  fi

  # Other positions: one record per position, SNPs before indels; drop
  # anything at, or spanning, the sentinel position.
  grep -v '^#' "${tag}.split.vcf" | awk -F'\t' -v OFS='\t' -v p="${pos}" '
    $2 != p && !($2 < p && $2 + length($4) - 1 >= p) {
      print ((length($4) == 1 && length($5) == 1) ? 0 : 1), $0
    }' | sort -t$'\t' -k3,3n -k1,1n -s | awk -F'\t' '!seen[$3]++' | cut -f2- > "${tag}.others.txt"

  cat "${tag}.index.txt" "${tag}.others.txt" | sort -t$'\t' -k2,2n -s > "${tag}.body.txt"
  rm -f "${tag}.split.vcf" "${tag}.index.txt" "${tag}.others.txt"
  body_files+=("${tag}.body.txt")
done < regions.tsv

# Regions are already in chromosome order.
cat header.vcf "${body_files[@]}" | bgzip -c > "${out_vcf}"
rm -f header.vcf "${body_files[@]}" ./*.tbi 2>/dev/null || true
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

echo "LD reference written to ${out_vcf}"
