#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash pipelines/v2g/bin/run_v2g_pipeline.sh --config <config.tsv> [--step all|check|manifest|integrate]

Steps:
  check      Validate required config values and input files.
  manifest   Write config snapshot, file checksums, R session info and git status.
  integrate  Combine line-of-evidence files into V2G master and summary tables.
  all        Run check, manifest and integrate.
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

run_r() {
  local script="$1"
  echo "[v2g] ${script}"
  Rscript "${R_DIR}/${script}" "${CONFIG}"
}

case "${STEP}" in
  check)
    run_r "00_check_config.R"
    ;;
  manifest)
    run_r "01_create_reproducibility_manifest.R"
    ;;
  integrate)
    run_r "02_integrate_line_evidence.R"
    ;;
  all)
    run_r "00_check_config.R"
    run_r "01_create_reproducibility_manifest.R"
    run_r "02_integrate_line_evidence.R"
    ;;
  *)
    echo "Unknown step: ${STEP}" >&2
    usage
    exit 1
    ;;
esac
