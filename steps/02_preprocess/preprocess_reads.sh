#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set by run_preprocess.sh
READS_DIR=${READS_DIR:?}
OUT_DIR=${OUT_DIR:?}
MIN_LEN=${MIN_LEN:-0}
MAX_LEN=${MAX_LEN:-0}
ORIENT=${ORIENT:-"none"}  # none | rc | anchor | igblast
ANCHOR=${ANCHOR:-""}
THREADS=${THREADS:-8}
SEQ2SCFV_IMAGE_NAME=${SEQ2SCFV_IMAGE_NAME:-"seq2scfv"}
GERMLINE_DIR=${GERMLINE_DIR:?}
ORGANISM="dog"
IG_SEQTYPE=${IG_SEQTYPE:-"Ig"}
DOMAIN_SYSTEM=${DOMAIN_SYSTEM:-"imgt"}
GERMLINE_V_DB="databases/dog/dog_V"
GERMLINE_D_DB="databases/dog/dog_D"
GERMLINE_J_DB="databases/dog/dog_J"
GERMLINE_AUX="optional_file/dog_gl.aux"
KEEP_IGBLAST=${KEEP_IGBLAST:-0}
FORCE=${FORCE:-0}

mkdir -p "${OUT_DIR}"
OUT_DIR="$(realpath "${OUT_DIR}")"

# OUT_DIR is mounted as /work in the container
to_container_path() {
  echo "/work/$(basename "$1")"
}

if [ "${ORIENT}" = "igblast" ]; then
  GERMLINE_DIR="$(realpath "${GERMLINE_DIR}")"
  for required in "${GERMLINE_V_DB}" "${GERMLINE_D_DB}" "${GERMLINE_J_DB}" "${GERMLINE_AUX}"; do
    if [ ! -f "${GERMLINE_DIR}/${required}" ]; then
      echo "ERROR: Missing ${GERMLINE_DIR}/${required}" >&2
      exit 1
    fi
  done
fi

IGBLAST_INTERNAL_DATA_DIR="${GERMLINE_DIR}/internal_data/dog"
if [ "${ORIENT}" = "igblast" ] && [ ! -d "${IGBLAST_INTERNAL_DATA_DIR}" ]; then
  echo "ERROR: Missing ${IGBLAST_INTERNAL_DATA_DIR}" >&2
  exit 1
fi

shopt -s nullglob
for fastq in "${READS_DIR}"/*.fastq; do
  base=$(basename "${fastq}" .fastq)
  out_fa="${OUT_DIR}/${base}.fasta"
  if [ "${FORCE}" != "1" ] && [ -s "${out_fa}" ]; then
    echo "${base}.fastq: output exists, skipping"
    continue
  fi
  if [ "${ORIENT}" = "igblast" ]; then
    tmp_fa="${OUT_DIR}/${base}.tmp.fasta"
    ig_tsv="${OUT_DIR}/${base}.igblast.tsv"
    keep_ids="${OUT_DIR}/${base}.keep_ids.txt"
    rc_ids="${OUT_DIR}/${base}.rc_ids.txt"
    ambig_ids="${OUT_DIR}/${base}.ambig_ids.txt"

    "${SCRIPT_DIR}/preprocess_fastq_to_fasta.py" \
      --input "${fastq}" \
      --output "${tmp_fa}" \
      --min-len "${MIN_LEN}" \
      --max-len "${MAX_LEN}" \
      --orient none

    docker run --rm --init \
      --user "$(id -u):$(id -g)" \
      -v "${OUT_DIR}:/work" \
      -v "${GERMLINE_DIR}:/germ:ro" \
      -v "${IGBLAST_INTERNAL_DATA_DIR}:/igblast/ncbi-igblast-1.20.0/internal_data/dog:ro" \
      -w /work \
      "${SEQ2SCFV_IMAGE_NAME}" \
      bash -lc "igblastn \
        -germline_db_V /germ/${GERMLINE_V_DB} \
        -germline_db_D /germ/${GERMLINE_D_DB} \
        -germline_db_J /germ/${GERMLINE_J_DB} \
        -organism ${ORGANISM} \
        -ig_seqtype ${IG_SEQTYPE} \
        -num_threads ${THREADS} \
        -domain_system ${DOMAIN_SYSTEM} \
        -num_alignments_V 1 \
        -num_alignments_D 1 \
        -num_alignments_J 1 \
        -outfmt 19 \
        -extend_align5end \
        -extend_align3end \
        -auxiliary_data /germ/${GERMLINE_AUX} \
        -query $(to_container_path "${tmp_fa}") \
        -out $(to_container_path "${ig_tsv}")"

    "${SCRIPT_DIR}/extract_igblast_orientation.py" \
      --igblast-tsv "${ig_tsv}" \
      --keep-out "${keep_ids}" \
      --rc-out "${rc_ids}" \
      --ambig-out "${ambig_ids}"

    "${SCRIPT_DIR}/preprocess_fastq_to_fasta.py" \
      --input "${fastq}" \
      --output "${out_fa}" \
      --min-len "${MIN_LEN}" \
      --max-len "${MAX_LEN}" \
      --orient none \
      --keep-list "${keep_ids}" \
      --rc-list "${rc_ids}"

    if [ "${KEEP_IGBLAST}" != "1" ]; then
      rm -f "${tmp_fa}" "${ig_tsv}" "${keep_ids}" "${rc_ids}" "${ambig_ids}"
    fi
  else
    "${SCRIPT_DIR}/preprocess_fastq_to_fasta.py" \
      --input "${fastq}" \
      --output "${out_fa}" \
      --min-len "${MIN_LEN}" \
      --max-len "${MAX_LEN}" \
      --orient "${ORIENT}" \
      --anchor "${ANCHOR}"
  fi
done
