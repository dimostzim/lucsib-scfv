#!/bin/bash
# Read QC of each barcode FASTQ in --input. Outputs in --output:
#   2_filtered/reads/    reads after malformed-record removal (seqkit sana), adapter trimming (Porechop) and
#                        length/quality filtering (Filtlong)
#   */nanoplot/          NanoPlot report of each barcode after each of the three stages
#   */logs/              seqkit sana, Porechop and Filtlong logs
#   3_multiqc/           MultiQC report of all NanoPlot reports
#   pipeline_stats.tsv   reads per barcode before QC and after each stage
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash steps/01_qc/run_qc.sh --input DIR --output DIR --cpus N --jobs N --image NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) READS_DIR="$2"; shift 2 ;;
    --output) OUT_DIR="$2"; shift 2 ;;
    --cpus) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# on the host: run this script again inside the QC image
if [[ ! -f /.dockerenv ]]; then
  READS_DIR="$(realpath -e "${READS_DIR}")"
  mkdir -p "${OUT_DIR}"
  OUT_DIR="$(realpath "${OUT_DIR}")"
  exec docker run --rm --user "$(id -u):$(id -g)" \
    -v "${READS_DIR}:/data/reads:ro" \
    -v "${OUT_DIR}:/data/qc_results" \
    -v "$(realpath "${BASH_SOURCE[0]}"):/opt/run_qc.sh:ro" \
    "${IMAGE}" bash /opt/run_qc.sh --input /data/reads --output /data/qc_results --cpus "${THREADS}" --jobs "${JOBS}"
fi

# The prefix names the stage (cleaned, porechop, filtlong); with the bare barcode the three stages collide in
# MultiQC and only one is shown ("_trimmed" would be stripped by MultiQC, hence "porechop").
run_nanoplot() {
  local fastq="$1" outdir="$2" prefix="$3"
  NanoPlot --fastq "${fastq}" --outdir "${outdir}" --prefix "${prefix}" --threads "${THREADS}" --plots dot \
    > "${outdir}.nanoplot.log" 2>&1
}

qc_barcode() {
  local sample="$1"
  local raw="${READS_DIR}/${sample}.fastq"
  local cleaned="${OUT_DIR}/0_cleaned/reads/${sample}.fastq"
  local trimmed="${OUT_DIR}/1_trimmed/reads/${sample}.fastq"
  local filtered="${OUT_DIR}/2_filtered/reads/${sample}.fastq"
  local raw_count stage_counts

  echo "QC ${sample}"
  seqkit sana "${raw}" -o "${cleaned}" 2> "${OUT_DIR}/0_cleaned/logs/${sample}.log"
  # reads before QC: records that pass seqkit sana plus its discarded malformed lines / 4
  raw_count=$(awk '/Pass records:/ {printf "%d", $(NF-3) + int($NF / 4)}' "${OUT_DIR}/0_cleaned/logs/${sample}.log")
  run_nanoplot "${cleaned}" "${OUT_DIR}/0_cleaned/nanoplot/${sample}" "${sample}_cleaned_"
  porechop -i "${cleaned}" -o "${trimmed}" --threads "${THREADS}" --discard_middle \
    > "${OUT_DIR}/1_trimmed/logs/${sample}.log" 2>&1
  run_nanoplot "${trimmed}" "${OUT_DIR}/1_trimmed/nanoplot/${sample}" "${sample}_porechop_"
  filtlong --min_length 1000 --keep_percent 90 "${trimmed}" \
    > "${filtered}" 2> "${OUT_DIR}/2_filtered/logs/${sample}.log"
  run_nanoplot "${filtered}" "${OUT_DIR}/2_filtered/nanoplot/${sample}" "${sample}_filtlong_"
  stage_counts=$(seqkit stats --quiet -T "${cleaned}" "${trimmed}" "${filtered}" | awk 'NR > 1 {printf "\t%s", $4}')
  echo -e "${sample}\t${raw_count}${stage_counts}" > "${OUT_DIR}/read_counts/${sample}.tsv"
}

mkdir -p "${OUT_DIR}"/{0_cleaned,1_trimmed,2_filtered}/{reads,logs,nanoplot} "${OUT_DIR}/read_counts"

# at most --jobs barcodes at a time; a failed barcode stops the run
pids=()
for fastq in "${READS_DIR}"/*.fastq; do
  if (( ${#pids[@]} >= JOBS )); then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
  qc_barcode "$(basename "${fastq}" .fastq)" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "${pid}"
done

multiqc "${OUT_DIR}/0_cleaned/nanoplot" "${OUT_DIR}/1_trimmed/nanoplot" "${OUT_DIR}/2_filtered/nanoplot" \
  --outdir "${OUT_DIR}/3_multiqc" --title "Nanopore QC Report" --filename multiqc_report --force

{
  echo -e "sample\traw\tcleaned\ttrimmed\tfiltered"
  cat "${OUT_DIR}"/read_counts/*.tsv
} > "${OUT_DIR}/pipeline_stats.tsv"
rm -r "${OUT_DIR}/read_counts" "${OUT_DIR}/0_cleaned/reads" "${OUT_DIR}/1_trimmed/reads"
cat "${OUT_DIR}/pipeline_stats.tsv"
