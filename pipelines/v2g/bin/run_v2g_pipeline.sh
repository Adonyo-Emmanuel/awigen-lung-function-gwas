#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash pipelines/v2g/bin/run_v2g_pipeline.sh --config <config.tsv> [--step all|check|manifest|standardise-lines|integrate]

Steps:
  check              Validate required config values and input files.
  manifest           Write input file checksums, config snapshot, R session info and git status.
  standardise-lines  Build one standard evidence table for each V2G line of evidence.
  integrate          Combine standard evidence tables into all-source V2G summary tables.
  all                Run check, manifest, standardise-lines and integrate.
USAGE
}

CONFIG=""
STEP="all"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --step)
      STEP="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "${CONFIG}" ]]; then
  echo "Missing --config" >&2
  usage
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
R_DIR="${PIPELINE_DIR}/R"
LINE_DIR="${PIPELINE_DIR}/lines_of_evidence"
FINAL_DIR="${PIPELINE_DIR}/functional_annotation_across_evidence_lines"

run_r_script() {
  local script_path="$1"
  echo "[v2g] ${script_path#${PIPELINE_DIR}/}"
  Rscript "${script_path}" "${CONFIG}"
}

run_line_evidence_scripts() {
  run_r_script "${LINE_DIR}/nearest_gene_annotation/build_nearest_gene_evidence.R"
  run_r_script "${LINE_DIR}/credible_set_variant_annotation/build_credible_set_variant_annotation_evidence.R"
  run_r_script "${LINE_DIR}/eqtl_associations/build_eqtl_evidence.R"
  run_r_script "${LINE_DIR}/pqtl_associations/build_pqtl_evidence.R"
  run_r_script "${LINE_DIR}/rare_variant_associations/build_rare_variant_evidence.R"
  run_r_script "${LINE_DIR}/mendelian_respiratory_disease_genes/build_mendelian_respiratory_disease_evidence.R"
  run_r_script "${LINE_DIR}/mouse_knockout_respiratory_phenotypes/build_mouse_knockout_evidence.R"
}

case "${STEP}" in
  check)
    run_r_script "${R_DIR}/check_v2g_config.R"
    ;;
  manifest)
    run_r_script "${R_DIR}/write_input_file_manifest.R"
    ;;
  standardise-lines)
    run_line_evidence_scripts
    ;;
  integrate)
    run_r_script "${FINAL_DIR}/build_all_source_v2g_tables.R"
    ;;
  all)
    run_r_script "${R_DIR}/check_v2g_config.R"
    run_r_script "${R_DIR}/write_input_file_manifest.R"
    run_line_evidence_scripts
    run_r_script "${FINAL_DIR}/build_all_source_v2g_tables.R"
    ;;
  *)
    echo "Unknown step: ${STEP}" >&2
    usage
    exit 1
    ;;
esac
