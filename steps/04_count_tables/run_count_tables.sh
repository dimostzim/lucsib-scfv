#!/bin/bash
# Render count_tables.Rmd on the annotation output. All outputs go to --output-dir:
#   r_output/tables/   count tables, attrition table, per-read filtered table (filt_df.csv)
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

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(realpath "${OUTPUT_DIR}")"
ANNOTATION_DIR="$(realpath "${ANNOTATION_DIR}")"

PATH="${R_ENV}/bin:${PATH}" "${R_ENV}/bin/Rscript" --vanilla -e "rmarkdown::render(
  input = '${SCRIPT_DIR}/count_tables.Rmd',
  output_dir = '${OUTPUT_DIR}',
  knit_root_dir = '${OUTPUT_DIR}',
  intermediates_dir = '${OUTPUT_DIR}',
  params = list(annotation_dir = '${ANNOTATION_DIR}')
)"

echo "Count tables written to ${OUTPUT_DIR}"
