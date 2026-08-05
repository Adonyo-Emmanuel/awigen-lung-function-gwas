#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${QTL_ROOT:-}" ]]; then
  CANDIDATE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
  if [[ -d "${CANDIDATE_ROOT}/functional_annotations_across_all_lines_of_evidences" ]]; then
    QTL_ROOT="${CANDIDATE_ROOT}"
  else
    echo "Set QTL_ROOT to the local QTL_analyses directory before running this wrapper." >&2
    echo "Example: QTL_ROOT=/path/to/QTL_analyses bash pipelines/v2g/bin/run_current_qtl_v2g_finalisation.sh" >&2
    exit 1
  fi
fi

FINAL_DIR="${QTL_ROOT}/functional_annotations_across_all_lines_of_evidences"

UPDATE_SCRIPT="${FINAL_DIR}/10_update_all_sources_with_rare_variant_associations.R"
TABLE6_SCRIPT="${FINAL_DIR}/11_create_awigen_shrine_like_table6_sentinel_v2g_evidence.R"

if [[ ! -f "${UPDATE_SCRIPT}" ]]; then
  echo "Missing final integration script: ${UPDATE_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${TABLE6_SCRIPT}" ]]; then
  echo "Missing sentinel evidence script: ${TABLE6_SCRIPT}" >&2
  exit 1
fi

echo "[v2g-current] Running final all-source integration"
Rscript "${UPDATE_SCRIPT}"

echo "[v2g-current] Running Shrine-like sentinel evidence table"
Rscript "${TABLE6_SCRIPT}"

echo "[v2g-current] Final outputs"
wc -l \
  "${FINAL_DIR}/awigen_v2g_all_source_evidence_master_table.tsv" \
  "${FINAL_DIR}/awigen_v2g_prioritised_genes_shrine_like.tsv" \
  "${FINAL_DIR}/awigen_shrine_like_table6_sentinel_v2g_evidence.tsv"
