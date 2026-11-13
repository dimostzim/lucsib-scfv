#!/bin/bash
# Annotate the preprocessed reads with seq2scfv in its Docker image. Outputs in --output-dir:
#   2.Catalogued/ to 6.Flags/  seq2scfv intermediate files
#   7.Counts/                  read counts of each scFv per barcode
#   8.Attrition/               attrition table
#   focus_libs.txt             barcode names, one per line
#   logs/pipeline.log          seq2scfv log
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash steps/03_annotate/run_seq2scfv.sh --input-preprocessed DIR --output-dir DIR --germlines DIR --threads N
    --image NAME
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-preprocessed) INPUT_PREPROCESSED="$2"; shift 2 ;;
    --output-dir) RUN_DIR="$2"; shift 2 ;;
    --germlines) GERMLINE_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --image) IMAGE_NAME="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# Docker mounts need absolute paths
INPUT_PREPROCESSED="$(realpath "${INPUT_PREPROCESSED}")"
GERMLINE_DIR="$(realpath "${GERMLINE_DIR}")"
RUN_DIR="$(realpath -m "${RUN_DIR}")"

mapfile -t FASTA_FILES < <(find "${INPUT_PREPROCESSED}" -maxdepth 1 -type f -name '*.fasta' | sort)
if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No FASTA files found in ${INPUT_PREPROCESSED}" >&2
  exit 1
fi

mkdir -p "${RUN_DIR}"/{2.Catalogued,3.vscan,4.scFv_split,5.Linkers,6.Flags,7.Counts,8.Attrition,logs}

for f in "${FASTA_FILES[@]}"; do
  basename "$f" .fasta
done > "${RUN_DIR}/focus_libs.txt"

cat > "${RUN_DIR}/_container_pipeline.sh" <<'EOF'
#!/bin/bash
set -euo pipefail

echo "[2/7] Catalog: SEGUID + merge + dedup"
cd /work/2.Catalogued
for fasta in /input/*.fasta; do
  uniq_id.py "$(basename "$fasta" .fasta)" "$fasta"
done
cat ./*.nt_correspondence.tsv > merged.nt_correspondence.tsv
seqkit rmdup -n < <(cat ./*.nt_seguid.fa) > merged.nt_seguid.fa

echo "[3/7] V-scan: iterative IgBLAST"
cd /work/3.vscan
vscan_functions.py \
  --query /work/2.Catalogued/merged.nt_seguid.fa \
  --germline_V /germ/databases/dog/dog_V \
  --germline_D /germ/databases/dog/dog_D \
  --germline_J /germ/databases/dog/dog_J \
  --organism dog \
  --ig_seqtype Ig \
  --domain_system imgt \
  --auxiliary_data /germ/optional_file/dog_gl.aux \
  --num_threads "${THREADS}"
rm -f ddDNA*tsv ddDNA*fasta

echo "[4/7] scFv delimitation"
cd /work/4.scFv_split
scFv_splitter.py \
  --igBLAST-hits /work/3.vscan/igBLAST.tsv \
  --domain-type imgt

if [[ ! -s /work/4.scFv_split/in_frame_igBLAST_paired_delim.tsv ]]; then
  echo "ERROR: in_frame_igBLAST_paired_delim.tsv is missing or empty." >&2
  exit 1
fi

echo "[5/7] Linker inference/scoring"
cd /work/5.Linkers
seqkit sample -p 0.1 /work/4.scFv_split/nt_inferred_linkers.fasta -o nt_inferred_linkers_sampled.fasta
clustalo --use-kimura -i nt_inferred_linkers_sampled.fasta -o nt_inferred_linkers.aln
seqkit sample -p 0.1 /work/4.scFv_split/aa_inferred_linkers.fasta -o aa_inferred_linkers_sampled.fasta
clustalo --use-kimura -i aa_inferred_linkers_sampled.fasta -o aa_inferred_linkers.aln
get_logos.py --nt-aln nt_inferred_linkers.aln --aa-aln aa_inferred_linkers.aln --check-linker undefined
linker_scorer.py --scFv /work/4.scFv_split/in_frame_igBLAST_paired_delim.tsv --linker aa_reference_linker.fasta

echo "[6/7] Flags"
cd /work/6.Flags
scfv_flags.py \
  --scFv /work/5.Linkers/in_frame_igBLAST_paired_delim_linker_scored.tsv \
  --linker-min-id 90 \
  --mismatches 2 \
  --linker-overhang undefined \
  --chain-length undefined

echo "[7/7] Per-library counts"
cd /work/7.Counts
ntqual_counts.py \
  --scFv /work/6.Flags/in_frame_igBLAST_paired_delim_linker_scored_flags.tsv \
  --correspondence /work/2.Catalogued/merged.nt_correspondence.tsv \
  --focus-libraries /work/focus_libs.txt
EOF

docker run --rm --init \
  --user "$(id -u):$(id -g)" \
  -e THREADS="${THREADS}" \
  -v "${INPUT_PREPROCESSED}:/input:ro" \
  -v "${GERMLINE_DIR}:/germ:ro" \
  -v "${GERMLINE_DIR}/internal_data/dog:/igblast/ncbi-igblast-1.20.0/internal_data/dog:ro" \
  -v "${RUN_DIR}:/work" \
  "${IMAGE_NAME}" \
  bash /work/_container_pipeline.sh | tee "${RUN_DIR}/logs/pipeline.log"

# rows below the header line; fails if the table is missing
count_rows() {
  tail -n +2 "$1" | wc -l
}

# the correspondence table has no header: one line per read
reads=$(wc -l < "${RUN_DIR}/2.Catalogued/merged.nt_correspondence.tsv")
unique_seqs=$(grep -c '^>' "${RUN_DIR}/2.Catalogued/merged.nt_seguid.fa")
igblast_rows=$(count_rows "${RUN_DIR}/3.vscan/igBLAST.tsv")
non_paired=$(count_rows "${RUN_DIR}/4.scFv_split/igBLAST_non_paired_nondelim.tsv")
out_of_frame=$(count_rows "${RUN_DIR}/4.scFv_split/out_of_frame_igBLAST_paired_delim.tsv")
nonstandard=$(count_rows "${RUN_DIR}/4.scFv_split/in_frame_igBLAST_delim_nonstandard_aas.tsv")
validated=$(count_rows "${RUN_DIR}/4.scFv_split/in_frame_igBLAST_paired_delim.tsv")
{
  echo -e "Category\tCount"
  echo -e "Reads (all barcodes)\t${reads}"
  echo -e "Unique read sequences\t${unique_seqs}"
  echo -e "Ig domain hits (IgBLAST rows)\t${igblast_rows}"
  echo -e "Non-paired / multi VH/VL\t${non_paired}"
  echo -e "Out-of-frame scFvs\t${out_of_frame}"
  echo -e "Non-standard AAs\t${nonstandard}"
  echo -e "Validated scFvs\t${validated}"
} > "${RUN_DIR}/8.Attrition/attrition_table.tsv"

echo "Annotation written to ${RUN_DIR}"
