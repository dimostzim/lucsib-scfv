#!/bin/bash
# TIM-3 binder analysis on the step-03 counts table. Outputs in --output-dir:
#   matching/              Sanger matching, cutoff sensitivity, first-pass and final labels
#   report/                label report, top clones, conflict audit, labelled FASTA, full labelled table
#   diversity_sweep/       abundance-first non-redundant strict-binder sets at 95-99% identity
#   representative_sets/   the selected sets (protein 95, 96, 97%; DNA 99%)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash steps/06_tim3_binders/run_tim3.sh --counts-tsv FILE --sanger-tsv FILE --output-dir DIR --python-env DIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --counts-tsv) COUNTS_TSV="$2"; shift 2 ;;
    --sanger-tsv) SANGER_TSV="$2"; shift 2 ;;
    --output-dir) OUT="$2"; shift 2 ;;
    --python-env) PY_ENV="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

PY="${PY_ENV}/bin/python"
MATCHING="${OUT}/matching"

echo "PROGRESS 1/5 Sanger matching and first-pass labels"
"${PY}" "${SCRIPT_DIR}/sanger_matching.py" --sanger-tsv "${SANGER_TSV}" --counts-tsv "${COUNTS_TSV}" \
  --out-dir "${MATCHING}"

echo "PROGRESS 2/5 full-read context for ambiguous clones and final labels"
"${PY}" "${SCRIPT_DIR}/ambiguous_context.py" --min-reads 10 --sanger-tsv "${SANGER_TSV}" \
  --counts-tsv "${COUNTS_TSV}" --labels-tsv "${MATCHING}/tim3_conservative_labels.tsv" --out-dir "${MATCHING}"

echo "PROGRESS 3/5 label report and exports"
"${PY}" "${SCRIPT_DIR}/label_report.py" --analysis-dir "${MATCHING}" --sanger-tsv "${SANGER_TSV}" \
  --out-dir "${OUT}/report"

echo "PROGRESS 4/5 diversity sweep"
"${PY}" "${SCRIPT_DIR}/diversity_sweep.py" \
  --input-tsv "${OUT}/report/tim3_full_labeled_scfv_results_with_sequences.tsv" --out-dir "${OUT}/diversity_sweep"

echo "PROGRESS 5/5 representative sets"
"${PY}" "${SCRIPT_DIR}/representative_sets.py" --sweep-dir "${OUT}/diversity_sweep" \
  --out-dir "${OUT}/representative_sets"

echo "TIM-3 outputs written to ${OUT}"
