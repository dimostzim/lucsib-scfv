#!/bin/bash
# Render count_tables.Rmd on the annotation output. All outputs go to --output-dir:
#   r_output/tables/   count tables, attrition tables, per-read filtered table (filt_df.csv)
#   r_output/*/        figures
#   count_tables.html  rendered notebook
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash steps/04_count_tables/run_count_tables.sh --annotation-dir DIR --output-dir DIR --r-env DIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --annotation-dir) ANNOTATION_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --r-env) R_ENV="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

COUNTS="${ANNOTATION_DIR}/7.Counts/in_frame_igBLAST_paired_delim_linker_scored_flags_counts.tsv"
ATTRITION="${ANNOTATION_DIR}/8.Attrition/attrition_table.tsv"
FOCUS_LIBS="${ANNOTATION_DIR}/focus_libs.txt"
for f in "${COUNTS}" "${ATTRITION}" "${FOCUS_LIBS}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: missing annotation output ${f}" >&2
    exit 1
  fi
done

# the notebook reads its inputs from scfv_extraction_output/ relative to the output directory
mkdir -p "${OUTPUT_DIR}/scfv_extraction_output"
OUTPUT_DIR="$(realpath "${OUTPUT_DIR}")"
ln -sfn "$(realpath "${COUNTS}")" "${OUTPUT_DIR}/scfv_extraction_output/$(basename "${COUNTS}")"
cp "${ATTRITION}" "${OUTPUT_DIR}/scfv_extraction_output/attrition_table.tsv"
cp "${FOCUS_LIBS}" "${OUTPUT_DIR}/scfv_extraction_output/focus_libs.txt"

PATH="${R_ENV}/bin:${PATH}" "${R_ENV}/bin/Rscript" --vanilla -e "rmarkdown::render(
  input = '${SCRIPT_DIR}/count_tables.Rmd',
  output_dir = '${OUTPUT_DIR}',
  knit_root_dir = '${OUTPUT_DIR}',
  intermediates_dir = '${OUTPUT_DIR}'
)"

for f in r_output/tables/VH_VL.csv r_output/tables/CDR3s.csv r_output/tables/HCDR3.csv r_output/tables/filt_df.csv; do
  if [[ ! -s "${OUTPUT_DIR}/${f}" ]]; then
    echo "ERROR: notebook did not write ${f}" >&2
    exit 1
  fi
done
echo "Count tables written to ${OUTPUT_DIR}"
