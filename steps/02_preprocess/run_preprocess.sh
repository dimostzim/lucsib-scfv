#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PREPROCESS_DIR="${SCRIPT_DIR}"

INPUT_DIR="${INPUT_DIR:-}"
OUTPUT_DIR="${ROOT_DIR}/preprocess_output"
GERMLINE_DIR="${ROOT_DIR}/germlines/dog"
THREADS=8
MIN_LEN=0
MAX_LEN=0
ORIENT="igblast"
ANCHOR=""
KEEP_IGBLAST=0
FORCE=0
SKIP_PREPROCESS=0
SKIP_CLEAN=0

SEQ2SCFV_IMAGE_NAME="${SEQ2SCFV_IMAGE_NAME:-seq2scfv}"
ORGANISM="dog"
IG_SEQTYPE="${IG_SEQTYPE:-Ig}"
DOMAIN_SYSTEM="${DOMAIN_SYSTEM:-imgt}"

usage() {
  cat <<'EOF'
Usage:
  bash scfv_analysis/preprocess/run_preprocess_pipeline.sh [options]

Pipeline steps:
  1) preprocess_reads.sh (FASTQ -> FASTA with optional orientation)

Options:
  --input DIR                 Input FASTQ directory (required unless --skip-preprocess)
  --output DIR                Output FASTA directory (default: scfv_analysis/preprocess_output)
  --germlines DIR             Germline output/input directory (default: scfv_analysis/germlines/dog)
  --cpus N                    Threads for IgBLAST (default: 8)
  --threads N                 Alias for --cpus
  --orient MODE               none | igblast (default: igblast)
  --keep-igblast              Keep IgBLAST/debug files in output
  --force                     Rebuild output FASTA even if present
  --skip-preprocess           Skip step 1
  --skip-clean                Skip cleanup of intermediate files
  --help                      Show this help

Environment overrides:
  SEQ2SCFV_IMAGE_NAME,
  IG_SEQTYPE, DOMAIN_SYSTEM
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT_DIR="$2"; shift 2 ;;
    --output)
      OUTPUT_DIR="$2"; shift 2 ;;
    --germlines)
      GERMLINE_DIR="$2"; shift 2 ;;
    --cpus|--threads)
      THREADS="$2"; shift 2 ;;
    --orient)
      ORIENT="$2"; shift 2 ;;
    --anchor)
      ANCHOR="$2"; shift 2 ;;
    --min-len)
      MIN_LEN="$2"; shift 2 ;;
    --max-len)
      MAX_LEN="$2"; shift 2 ;;
    --keep-igblast)
      KEEP_IGBLAST=1; shift ;;
    --force)
      FORCE=1; shift ;;
    --skip-preprocess)
      SKIP_PREPROCESS=1; shift ;;
    --skip-clean)
      SKIP_CLEAN=1; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

if [[ "${ORIENT}" == "anchor" && -z "${ANCHOR}" ]]; then
  echo "ERROR: --orient anchor requires --anchor" >&2
  exit 1
fi

if [[ "${SKIP_PREPROCESS}" != "1" ]]; then
  if [[ -z "${INPUT_DIR}" ]]; then
    echo "ERROR: --input is required unless --skip-preprocess is set." >&2
    exit 1
  fi
  if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: Input directory not found: ${INPUT_DIR}" >&2
    exit 1
  fi
fi

mkdir -p "${OUTPUT_DIR}"

echo "Pipeline config:"
echo "  input: ${INPUT_DIR}"
echo "  output: ${OUTPUT_DIR}"
echo "  germlines: ${GERMLINE_DIR}"
echo "  orient: ${ORIENT}"
echo "  threads: ${THREADS}"
echo "  min_len: ${MIN_LEN}"
echo "  max_len: ${MAX_LEN}"
echo "  image: ${SEQ2SCFV_IMAGE_NAME}"

clean_intermediate_output() {
  local out_dir="$1"
  local before after removed
  before=$(find "${out_dir}" -maxdepth 1 -type f \
    \( -name '*.tmp.fasta' -o -name '*.igblast.tsv' -o -name '*.keep_ids.txt' -o -name '*.rc_ids.txt' -o -name '*.ambig_ids.txt' \) | wc -l)
  if [[ "${before}" -eq 0 ]]; then
    echo "No intermediate files to clean."
    return
  fi
  find "${out_dir}" -maxdepth 1 -type f \
    \( -name '*.tmp.fasta' -o -name '*.igblast.tsv' -o -name '*.keep_ids.txt' -o -name '*.rc_ids.txt' -o -name '*.ambig_ids.txt' \) -delete
  after=$(find "${out_dir}" -maxdepth 1 -type f \
    \( -name '*.tmp.fasta' -o -name '*.igblast.tsv' -o -name '*.keep_ids.txt' -o -name '*.rc_ids.txt' -o -name '*.ambig_ids.txt' \) | wc -l)
  removed=$((before - after))
  echo "Cleaned intermediate files: removed=${removed}, remaining=${after}"
}

if [[ "${SKIP_PREPROCESS}" != "1" ]]; then
  echo "[1/2] Preprocessing reads"
  READS_DIR="${INPUT_DIR}" \
  OUT_DIR="${OUTPUT_DIR}" \
  MIN_LEN="${MIN_LEN}" \
  MAX_LEN="${MAX_LEN}" \
  ORIENT="${ORIENT}" \
  ANCHOR="${ANCHOR}" \
  THREADS="${THREADS}" \
  SEQ2SCFV_IMAGE_NAME="${SEQ2SCFV_IMAGE_NAME}" \
  GERMLINE_DIR="${GERMLINE_DIR}" \
  ORGANISM="${ORGANISM}" \
  IG_SEQTYPE="${IG_SEQTYPE}" \
  DOMAIN_SYSTEM="${DOMAIN_SYSTEM}" \
  KEEP_IGBLAST="${KEEP_IGBLAST}" \
  FORCE="${FORCE}" \
    bash "${PREPROCESS_DIR}/preprocess_reads.sh"
else
  echo "[1/2] Skipped preprocessing"
fi

if [[ "${SKIP_CLEAN}" == "1" ]]; then
  echo "[2/2] Skipped cleanup (--skip-clean)"
elif [[ "${KEEP_IGBLAST}" == "1" ]]; then
  echo "[2/2] Skipped cleanup (--keep-igblast keeps debug/intermediate files)"
else
  echo "[2/2] Cleaning intermediate files"
  clean_intermediate_output "${OUTPUT_DIR}"
fi

echo "Pipeline finished."
