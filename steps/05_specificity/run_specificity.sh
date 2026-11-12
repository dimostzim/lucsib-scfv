#!/bin/bash
# Specificity analyses on the step-04 tables. Outputs in --output-dir:
#   candidates/              sequences present in H460 PCNB R4 (barcode09) and absent from huAEC R4 (barcode12)
#   r4_h460pcnb_families/    sequence families of those candidates
#   all_rounds/              families, sample centroids and sample overlap over all rounds and targets
#   figures/                 figures 1 to 4
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash steps/05_specificity/run_specificity.sh --tables-dir DIR --samples-csv FILE --output-dir DIR --python-env DIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tables-dir) TABLES_DIR="$2"; shift 2 ;;
    --samples-csv) SAMPLES_CSV="$2"; shift 2 ;;
    --output-dir) OUT="$2"; shift 2 ;;
    --python-env) PY_ENV="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

PY="${PY_ENV}/bin/python"
FILT_DF="${TABLES_DIR}/filt_df.csv"
mkdir -p "${OUT}/figures"

echo "PROGRESS 1/6 candidates present in barcode09 and absent from barcode12"
"${PY}" "${SCRIPT_DIR}/extract_candidates.py" --input-csv "${FILT_DF}" --outdir "${OUT}/candidates"

echo "PROGRESS 2/6 sequence families of the candidates"
"${PY}" "${SCRIPT_DIR}/r4_h460pcnb_families.py" --candidates-csv "${OUT}/candidates/candidates.csv" \
  --outdir "${OUT}/r4_h460pcnb_families"

echo "PROGRESS 3/6 all rounds and targets"
"${PY}" "${SCRIPT_DIR}/all_rounds_map.py" --input-csv "${FILT_DF}" --barcode-map "${SAMPLES_CSV}" \
  --outdir "${OUT}/all_rounds"

echo "PROGRESS 4/6 figures 1 and 2"
"${PY}" "${SCRIPT_DIR}/make_family_figures.py" --allseq-dir "${OUT}/r4_h460pcnb_families" \
  --allrounds-dir "${OUT}/all_rounds" --outdir "${OUT}/figures"

echo "PROGRESS 5/6 figure 3"
"${PY}" "${SCRIPT_DIR}/make_trajectory_plot.py" --centroids-csv "${OUT}/all_rounds/sample_centroids_pca.csv" \
  --out "${OUT}/figures/figure3_target_trajectories.png"

echo "PROGRESS 6/6 figure 4"
"${PY}" "${SCRIPT_DIR}/cdr3_sharing_upset.py" --cdr3-csv "${TABLES_DIR}/CDR3s.csv" --samples-csv "${SAMPLES_CSV}" \
  --out-dir "${OUT}/figures"

echo "Specificity outputs written to ${OUT}"
