#!/bin/bash
# Run the whole lucsib scFv pipeline: download inputs and tools, QC, preprocessing,
# annotation, count tables, specificity analyses, TIM-3 binder analysis, results.
# Settings: config/settings.env. Zenodo token, if the inputs are downloaded: .env (see .env.example).
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${REPO_DIR}/config/settings.env"

usage() {
  cat <<'EOF'
Usage:
  bash run_all.sh [--from STEP] [--to STEP]

Steps:
  00  inputs and tools
  01  read QC
  02  read orientation
  03  scFv annotation (seq2scfv)
  04  count tables (figures 5-13)
  05  sequence families (figures 1-4)
  06  TIM-3 binders
  07  results
EOF
}

FROM=00
TO=07
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    --to) TO="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# resolved to a real absolute path, so the Docker steps mount real folders
WORK_DIR="$(cd "${REPO_DIR}" && mkdir -p "${WORK_DIR}" && realpath "${WORK_DIR}")"
INPUTS="${WORK_DIR}/inputs"
TOOLS="${WORK_DIR}/tools"
QC_DIR="${WORK_DIR}/01_qc"
PREPROCESS_DIR="${WORK_DIR}/02_preprocess"
ANNOTATION_DIR="${WORK_DIR}/03_annotate"
COUNT_TABLES_DIR="${WORK_DIR}/04_count_tables"
SPECIFICITY_DIR="${WORK_DIR}/05_specificity"
TIM3_DIR="${WORK_DIR}/06_tim3_binders"
RESULTS_DIR="${WORK_DIR}/results"

run_step() {
  local step="$1"
  [[ "10#${step}" -ge "10#${FROM}" && "10#${step}" -le "10#${TO}" ]]
}

# a step always starts from an empty output folder, so no partial or older files are reused
fresh_dir() {
  rm -rf "$1"
}

if run_step 00; then
  echo "PROGRESS 0/8 step 00: inputs and tools"
  if [[ -f "${REPO_DIR}/.env" ]]; then
    set -a; source "${REPO_DIR}/.env"; set +a
  fi
  python3 "${REPO_DIR}/steps/00_download/get_inputs.py" --record-id "${ZENODO_RECORD_ID}" \
    --manifest "${REPO_DIR}/config/inputs_sha256.tsv" --out-dir "${INPUTS}"
  bash "${REPO_DIR}/steps/00_download/get_tools.sh" --tools-dir "${TOOLS}" \
    --seq2scfv-url "${SEQ2SCFV_REPO_URL}" --seq2scfv-commit "${SEQ2SCFV_COMMIT}" \
    --dog-ref-url "${DOG_REF_REPO_URL}" --dog-ref-commit "${DOG_REF_COMMIT}" \
    --qc-image "${QC_IMAGE}" --seq2scfv-image "${SEQ2SCFV_IMAGE}"
fi

if run_step 01; then
  echo "PROGRESS 1/8 step 01: read QC"
  fresh_dir "${QC_DIR}"
  bash "${REPO_DIR}/steps/01_qc/run_qc.sh" --input "${INPUTS}/reads" --output "${QC_DIR}" \
    --cpus "${QC_THREADS_PER_JOB}" --jobs "${QC_JOBS}" --image "${QC_IMAGE}"
fi

if run_step 02; then
  echo "PROGRESS 2/8 step 02: read orientation"
  fresh_dir "${PREPROCESS_DIR}"
  bash "${REPO_DIR}/steps/02_preprocess/run_preprocess.sh" --input "${QC_DIR}/2_filtered/reads" \
    --output "${PREPROCESS_DIR}" --germlines "${TOOLS}/germlines/dog" \
    --cpus "${CPUS}" --image "${SEQ2SCFV_IMAGE}"
fi

if run_step 03; then
  echo "PROGRESS 3/8 step 03: scFv annotation"
  fresh_dir "${ANNOTATION_DIR}"
  bash "${REPO_DIR}/steps/03_annotate/run_seq2scfv.sh" --input-preprocessed "${PREPROCESS_DIR}" \
    --output-dir "${ANNOTATION_DIR}" --germlines "${TOOLS}/germlines/dog" --threads "${CPUS}" \
    --image "${SEQ2SCFV_IMAGE}"
fi

if run_step 04; then
  echo "PROGRESS 4/8 step 04: count tables"
  fresh_dir "${COUNT_TABLES_DIR}"
  bash "${REPO_DIR}/steps/04_count_tables/run_count_tables.sh" --annotation-dir "${ANNOTATION_DIR}" \
    --output-dir "${COUNT_TABLES_DIR}" --r-env "${TOOLS}/envs/r"
fi

if run_step 05; then
  echo "PROGRESS 5/8 step 05: sequence families"
  fresh_dir "${SPECIFICITY_DIR}"
  bash "${REPO_DIR}/steps/05_specificity/run_specificity.sh" --tables-dir "${COUNT_TABLES_DIR}/r_output/tables" \
    --samples-csv "${INPUTS}/samples.csv" --output-dir "${SPECIFICITY_DIR}" --python-env "${TOOLS}/envs/python"
fi

if run_step 06; then
  echo "PROGRESS 6/8 step 06: TIM-3 binders"
  fresh_dir "${TIM3_DIR}"
  bash "${REPO_DIR}/steps/06_tim3_binders/run_tim3.sh" \
    --counts-tsv "${ANNOTATION_DIR}/7.Counts/in_frame_igBLAST_paired_delim_linker_scored_flags_counts.tsv" \
    --sanger-tsv "${INPUTS}/TIM-3_sanger_binding.tsv" --output-dir "${TIM3_DIR}" --python-env "${TOOLS}/envs/python"
fi

if run_step 07; then
  echo "PROGRESS 7/8 step 07: results"
  fresh_dir "${RESULTS_DIR}"
  "${TOOLS}/envs/python/bin/python" "${REPO_DIR}/steps/07_collect_results/collect_results.py" \
    --work-dir "${WORK_DIR}" --output-dir "${RESULTS_DIR}" --repo-dir "${REPO_DIR}"
fi

echo "PROGRESS 8/8 done"
