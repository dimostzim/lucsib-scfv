#!/bin/bash
# Convert each barcode FASTQ in --input to --output/<barcode>.fasta, reverse-complementing the reads that IgBLAST
# aligns to the dog germline genes in reverse orientation.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash steps/02_preprocess/run_preprocess.sh --input DIR --output DIR --germlines DIR --cpus N --image NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT_DIR="$2"; shift 2 ;;
    --output) OUTPUT_DIR="$2"; shift 2 ;;
    --germlines) GERMLINE_DIR="$2"; shift 2 ;;
    --cpus) THREADS="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

GERMLINE_DIR="$(realpath -e "${GERMLINE_DIR}")"
# IgBLAST query and result of one barcode at a time, kept apart from the oriented FASTA files
IGBLAST_DIR="${OUTPUT_DIR}/igblast"
mkdir -p "${IGBLAST_DIR}"
IGBLAST_DIR="$(realpath "${IGBLAST_DIR}")"

for fastq in "${INPUT_DIR}"/*.fastq; do
  sample="$(basename "${fastq}" .fastq)"
  echo "Orienting ${sample}"
  python3 "${SCRIPT_DIR}/fastq_to_fasta.py" --input "${fastq}" --output "${IGBLAST_DIR}/${sample}.fasta"
  # -organism dog reads the dog germline annotation from IgBLAST's internal_data folder
  docker run --rm --init --user "$(id -u):$(id -g)" \
    -v "${IGBLAST_DIR}:/work" \
    -v "${GERMLINE_DIR}:/germ:ro" \
    -v "${GERMLINE_DIR}/internal_data/dog:/igblast/ncbi-igblast-1.20.0/internal_data/dog:ro" \
    "${IMAGE}" igblastn \
    -germline_db_V /germ/databases/dog/dog_V \
    -germline_db_D /germ/databases/dog/dog_D \
    -germline_db_J /germ/databases/dog/dog_J \
    -auxiliary_data /germ/optional_file/dog_gl.aux \
    -organism dog -ig_seqtype Ig -domain_system imgt \
    -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
    -extend_align5end -extend_align3end \
    -outfmt 19 -num_threads "${THREADS}" \
    -query "/work/${sample}.fasta" -out "/work/${sample}.tsv"
  python3 "${SCRIPT_DIR}/fastq_to_fasta.py" --input "${fastq}" --output "${OUTPUT_DIR}/${sample}.fasta" \
    --igblast-tsv "${IGBLAST_DIR}/${sample}.tsv"
  rm "${IGBLAST_DIR}/${sample}.fasta" "${IGBLAST_DIR}/${sample}.tsv"
done
rmdir "${IGBLAST_DIR}"
