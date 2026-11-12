#!/bin/bash
set -euo pipefail

INPUT_PREPROCESSED=""
RUN_DIR=""
LINKER_MODE="infer"            # infer | known
LINKER_SEQ=""
LINKER_MIN_ID="90"
MISMATCHES="2"
LINKER_OVERHANG="undefined"
CHAIN_LENGTH="undefined"
THREADS=8
LINKER_SAMPLE="0.1"
FORCE=0

IMAGE_NAME="lucsib-scfv-seq2scfv"
GERMLINE_DIR=""
ORGANISM="dog"
IG_SEQTYPE="${IG_SEQTYPE:-Ig}"
DOMAIN_SYSTEM="${DOMAIN_SYSTEM:-imgt}"
GERMLINE_V_DB="databases/dog/dog_V"
GERMLINE_D_DB="databases/dog/dog_D"
GERMLINE_J_DB="databases/dog/dog_J"
GERMLINE_AUX="optional_file/dog_gl.aux"
IGBLAST_INTERNAL_DATA_DIR=""

usage() {
  cat <<'EOF'
Usage:
  bash steps/03_annotate/run_seq2scfv.sh [options]

Required:
  --input-preprocessed DIR      Directory with preprocessed FASTA files (barcode*.fasta)
  --output-dir DIR              Output directory for this annotation run
  --germlines DIR               Dog germline directory (staged by steps/00_download/get_tools.sh)

Optional:
  --image NAME                  seq2scfv Docker image (default: lucsib-scfv-seq2scfv)
  --threads N                   Threads for IgBLAST (default: 8)
  --linker-mode MODE            infer | known (default: infer)
  --linker-seq AASEQ            Amino-acid linker sequence (required when --linker-mode known)
  --linker-sample P             Fraction for linker sampling before MSA, 0<P<=1 (default: 0.1)
  --linker-min-id N             scfv_flags threshold (default: 90)
  --mismatches N                scfv_flags threshold (default: 2)
  --linker-overhang N|undefined scfv_flags threshold (default: undefined)
  --chain-length N|undefined    scfv_flags threshold (default: undefined)
  --force                       Remove existing run dir if it exists
  --help                        Show this help

Environment overrides:
  IG_SEQTYPE, DOMAIN_SYSTEM
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-preprocessed)
      INPUT_PREPROCESSED="$2"; shift 2 ;;
    --output-dir)
      RUN_DIR="$2"; shift 2 ;;
    --germlines)
      GERMLINE_DIR="$2"; shift 2 ;;
    --image)
      IMAGE_NAME="$2"; shift 2 ;;
    --threads|--cpus)
      THREADS="$2"; shift 2 ;;
    --linker-mode)
      LINKER_MODE="$2"; shift 2 ;;
    --linker-seq)
      LINKER_SEQ="$2"; shift 2 ;;
    --linker-sample)
      LINKER_SAMPLE="$2"; shift 2 ;;
    --linker-min-id)
      LINKER_MIN_ID="$2"; shift 2 ;;
    --mismatches)
      MISMATCHES="$2"; shift 2 ;;
    --linker-overhang)
      LINKER_OVERHANG="$2"; shift 2 ;;
    --chain-length)
      CHAIN_LENGTH="$2"; shift 2 ;;
    --force)
      FORCE=1; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

if [[ -z "${INPUT_PREPROCESSED}" ]]; then
  echo "ERROR: --input-preprocessed is required." >&2
  exit 1
fi
if [[ -z "${RUN_DIR}" || -z "${GERMLINE_DIR}" ]]; then
  echo "ERROR: --output-dir and --germlines are required." >&2
  exit 1
fi
if [[ ! -d "${INPUT_PREPROCESSED}" ]]; then
  echo "ERROR: Input directory not found: ${INPUT_PREPROCESSED}" >&2
  exit 1
fi
if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --threads must be a positive integer." >&2
  exit 1
fi
if [[ "${LINKER_MODE}" != "infer" && "${LINKER_MODE}" != "known" ]]; then
  echo "ERROR: --linker-mode must be infer or known." >&2
  exit 1
fi
if [[ "${LINKER_MODE}" == "known" && -z "${LINKER_SEQ}" ]]; then
  echo "ERROR: --linker-seq is required when --linker-mode known." >&2
  exit 1
fi

# Normalize host paths for Docker mounts.
INPUT_PREPROCESSED="$(realpath "${INPUT_PREPROCESSED}")"
GERMLINE_DIR="$(realpath "${GERMLINE_DIR}")"
RUN_DIR="$(realpath -m "${RUN_DIR}")"

mapfile -t FASTA_FILES < <(find "${INPUT_PREPROCESSED}" -maxdepth 1 -type f -name '*.fasta' ! -name '*.tmp.fasta' | sort)
if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No FASTA files found in ${INPUT_PREPROCESSED}" >&2
  exit 1
fi

if [[ ! -f "${GERMLINE_DIR}/${GERMLINE_V_DB}" || ! -f "${GERMLINE_DIR}/${GERMLINE_D_DB}" || ! -f "${GERMLINE_DIR}/${GERMLINE_J_DB}" || ! -f "${GERMLINE_DIR}/${GERMLINE_AUX}" ]]; then
  echo "ERROR: Missing germline files in ${GERMLINE_DIR}." >&2
  echo "Required: ${GERMLINE_V_DB}, ${GERMLINE_D_DB}, ${GERMLINE_J_DB}, ${GERMLINE_AUX}" >&2
  echo "Run: bash steps/00_download/get_tools.sh" >&2
  exit 1
fi

if ! docker image inspect "${IMAGE_NAME}" >/dev/null 2>&1; then
  echo "ERROR: Docker image '${IMAGE_NAME}' not found." >&2
  echo "Run: bash steps/00_download/get_tools.sh" >&2
  exit 1
fi

IGBLAST_INTERNAL_DATA_DIR="${GERMLINE_DIR}/internal_data/dog"
if [[ ! -d "${IGBLAST_INTERNAL_DATA_DIR}" ]]; then
  echo "ERROR: Missing ${IGBLAST_INTERNAL_DATA_DIR}" >&2
  exit 1
fi

if [[ -e "${RUN_DIR}" ]]; then
  if [[ "${FORCE}" == "1" ]]; then
    rm -rf "${RUN_DIR}"
  else
    echo "ERROR: Run directory already exists: ${RUN_DIR} (use --force to replace)" >&2
    exit 1
  fi
fi

mkdir -p "${RUN_DIR}"/{1.Preprocessed,2.Catalogued,3.vscan,4.scFv_split,5.Linkers,6.Flags,7.Counts,8.Attrition,logs}

# Record inputs and per-library identifiers used downstream.
{
  echo -e "library\tfasta_path"
  for f in "${FASTA_FILES[@]}"; do
    b="$(basename "$f" .fasta)"
    echo -e "${b}\t${f}"
  done
} > "${RUN_DIR}/1.Preprocessed/input_files.tsv"

for f in "${FASTA_FILES[@]}"; do
  basename "$f" .fasta
done > "${RUN_DIR}/focus_libs.txt"

CHECK_LINKER="undefined"
if [[ "${LINKER_MODE}" == "known" ]]; then
  CHECK_LINKER="${LINKER_SEQ}"
fi

cat > "${RUN_DIR}/run_config.env" <<EOF
RUN_DIR=${RUN_DIR}
INPUT_PREPROCESSED=${INPUT_PREPROCESSED}
IMAGE_NAME=${IMAGE_NAME}
GERMLINE_DIR=${GERMLINE_DIR}
ORGANISM=${ORGANISM}
IG_SEQTYPE=${IG_SEQTYPE}
DOMAIN_SYSTEM=${DOMAIN_SYSTEM}
GERMLINE_V_DB=${GERMLINE_V_DB}
GERMLINE_D_DB=${GERMLINE_D_DB}
GERMLINE_J_DB=${GERMLINE_J_DB}
GERMLINE_AUX=${GERMLINE_AUX}
THREADS=${THREADS}
LINKER_MODE=${LINKER_MODE}
LINKER_SAMPLE=${LINKER_SAMPLE}
LINKER_MIN_ID=${LINKER_MIN_ID}
MISMATCHES=${MISMATCHES}
LINKER_OVERHANG=${LINKER_OVERHANG}
CHAIN_LENGTH=${CHAIN_LENGTH}
EOF

cat > "${RUN_DIR}/_container_pipeline.sh" <<'EOF'
#!/bin/bash
set -euo pipefail

echo "[2/7] Catalog: SEGUID + merge + dedup"
cd /work/2.Catalogued
shopt -s nullglob
for fasta in /input/*.fasta; do
  # unoriented copies left by preprocessing must not be annotated as extra libraries
  [[ "$fasta" == *.tmp.fasta ]] && continue
  base="$(basename "$fasta" .fasta)"
  uniq_id.py "$base" "$fasta"
done
cat ./*.nt_correspondence.tsv > merged.nt_correspondence.tsv
seqkit rmdup -n < <(cat ./*.nt_seguid.fa) > merged.nt_seguid.fa

echo "[3/7] V-scan: iterative IgBLAST"
cd /work/3.vscan
vscan_functions.py \
  --query /work/2.Catalogued/merged.nt_seguid.fa \
  --germline_V /germ/${GERMLINE_V_DB} \
  --germline_D /germ/${GERMLINE_D_DB} \
  --germline_J /germ/${GERMLINE_J_DB} \
  --organism "${ORGANISM}" \
  --ig_seqtype "${IG_SEQTYPE}" \
  --domain_system "${DOMAIN_SYSTEM}" \
  --auxiliary_data /germ/${GERMLINE_AUX} \
  --num_threads "${THREADS}"
rm -f ddDNA*tsv ddDNA*fasta || true

echo "[4/7] scFv delimitation"
cd /work/4.scFv_split
scFv_splitter.py \
  --igBLAST-hits /work/3.vscan/igBLAST.tsv \
  --domain-type "${DOMAIN_SYSTEM}"

if [[ ! -s /work/4.scFv_split/in_frame_igBLAST_paired_delim.tsv ]]; then
  echo "ERROR: in_frame_igBLAST_paired_delim.tsv is missing or empty." >&2
  exit 1
fi

echo "[5/7] Linker inference/scoring"
cd /work/5.Linkers
seqkit sample -p "${LINKER_SAMPLE}" /work/4.scFv_split/nt_inferred_linkers.fasta -o nt_inferred_linkers_sampled.fasta
clustalo --use-kimura -i nt_inferred_linkers_sampled.fasta -o nt_inferred_linkers.aln
seqkit sample -p "${LINKER_SAMPLE}" /work/4.scFv_split/aa_inferred_linkers.fasta -o aa_inferred_linkers_sampled.fasta
clustalo --use-kimura -i aa_inferred_linkers_sampled.fasta -o aa_inferred_linkers.aln
get_logos.py --nt-aln nt_inferred_linkers.aln --aa-aln aa_inferred_linkers.aln --check-linker "${CHECK_LINKER}"
linker_scorer.py --scFv /work/4.scFv_split/in_frame_igBLAST_paired_delim.tsv --linker aa_reference_linker.fasta

echo "[6/7] Flags"
cd /work/6.Flags
scfv_flags.py \
  --scFv /work/5.Linkers/in_frame_igBLAST_paired_delim_linker_scored.tsv \
  --linker-min-id "${LINKER_MIN_ID}" \
  --mismatches "${MISMATCHES}" \
  --linker-overhang "${LINKER_OVERHANG}" \
  --chain-length "${CHAIN_LENGTH}"

echo "[7/7] Per-library counts"
cd /work/7.Counts
ntqual_counts.py \
  --scFv /work/6.Flags/in_frame_igBLAST_paired_delim_linker_scored_flags.tsv \
  --correspondence /work/2.Catalogued/merged.nt_correspondence.tsv \
  --focus-libraries /work/focus_libs.txt
EOF

chmod +x "${RUN_DIR}/_container_pipeline.sh"

echo "Run directory: ${RUN_DIR}"
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e THREADS="${THREADS}" \
  -e ORGANISM="${ORGANISM}" \
  -e IG_SEQTYPE="${IG_SEQTYPE}" \
  -e DOMAIN_SYSTEM="${DOMAIN_SYSTEM}" \
  -e GERMLINE_V_DB="${GERMLINE_V_DB}" \
  -e GERMLINE_D_DB="${GERMLINE_D_DB}" \
  -e GERMLINE_J_DB="${GERMLINE_J_DB}" \
  -e GERMLINE_AUX="${GERMLINE_AUX}" \
  -e LINKER_SAMPLE="${LINKER_SAMPLE}" \
  -e CHECK_LINKER="${CHECK_LINKER}" \
  -e LINKER_MIN_ID="${LINKER_MIN_ID}" \
  -e MISMATCHES="${MISMATCHES}" \
  -e LINKER_OVERHANG="${LINKER_OVERHANG}" \
  -e CHAIN_LENGTH="${CHAIN_LENGTH}" \
  -v "${INPUT_PREPROCESSED}:/input:ro" \
  -v "${GERMLINE_DIR}:/germ:ro" \
  -v "${IGBLAST_INTERNAL_DATA_DIR}:/igblast/ncbi-igblast-1.20.0/internal_data/dog:ro" \
  -v "${RUN_DIR}:/work" \
  "${IMAGE_NAME}" \
  bash /work/_container_pipeline.sh | tee "${RUN_DIR}/logs/pipeline.log"

count_tsv_rows() {
  local tsv="$1"
  if [[ -f "${tsv}" ]]; then
    local n
    n=$(wc -l < "${tsv}")
    if [[ "${n}" -gt 0 ]]; then
      echo $((n - 1))
    else
      echo 0
    fi
  else
    echo 0
  fi
}

count_fasta_headers() {
  local fa="$1"
  if [[ -f "${fa}" ]]; then
    grep -c '^>' "${fa}" || true
  else
    echo 0
  fi
}

{
  echo -e "Category\tCount"
  echo -e "Merged libraries\t$(count_tsv_rows "${RUN_DIR}/2.Catalogued/merged.nt_correspondence.tsv")"
  echo -e "De-duplicated reads\t$(count_fasta_headers "${RUN_DIR}/2.Catalogued/merged.nt_seguid.fa")"
  echo -e "Ig domains\t$(count_tsv_rows "${RUN_DIR}/3.vscan/igBLAST.tsv")"
  echo -e "Non-paired / multi VH/VL\t$(count_tsv_rows "${RUN_DIR}/4.scFv_split/igBLAST_non_paired_nondelim.tsv")"
  echo -e "Out-of-frame scFvs\t$(count_tsv_rows "${RUN_DIR}/4.scFv_split/out_of_frame_igBLAST_paired_delim.tsv")"
  echo -e "Non-standard AAs\t$(count_tsv_rows "${RUN_DIR}/4.scFv_split/in_frame_igBLAST_delim_nonstandard_aas.tsv")"
  echo -e "Validated scFvs\t$(count_tsv_rows "${RUN_DIR}/4.scFv_split/in_frame_igBLAST_paired_delim.tsv")"
} > "${RUN_DIR}/8.Attrition/attrition_table.tsv"

echo "Completed. Key outputs:"
echo "  /work/3.vscan/igBLAST.tsv"
echo "  /work/4.scFv_split/in_frame_igBLAST_paired_delim.tsv"
echo "  /work/5.Linkers/in_frame_igBLAST_paired_delim_linker_scored.tsv"
echo "  /work/6.Flags/in_frame_igBLAST_paired_delim_linker_scored_flags.tsv"
echo "  /work/7.Counts/in_frame_igBLAST_paired_delim_linker_scored_flags_counts.tsv"
echo "  /work/8.Attrition/attrition_table.tsv"
echo "Host path: ${RUN_DIR}"
