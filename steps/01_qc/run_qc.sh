#!/bin/bash
set -e
set -o pipefail

IMAGE_NAME="${QC_IMAGE:-lucsib-scfv-qc}"

THREADS="${THREADS:-64}"
JOBS="${JOBS:-1}"
MIN_LENGTH="${MIN_LENGTH:-1000}"
KEEP_PERCENT="${KEEP_PERCENT:-90}"
CLEAN_FASTQ="${CLEAN_FASTQ:-1}"

usage() {
    cat <<'EOF'
Usage: bash steps/01_qc/run_qc.sh --input <dir> --output <dir> [options]

Options:
  --input <dir>    Input directory containing *.fastq files
  --output <dir>   Output directory for QC results
  --cpus <n>       Number of CPU threads to use (maps to THREADS)
  --jobs <n>       Number of samples processed in parallel (maps to JOBS)
  --image <name>   QC Docker image (built by steps/00_download/get_tools.sh)
  --help           Show this help message

Environment overrides (optional):
  MIN_LENGTH, KEEP_PERCENT, CLEAN_FASTQ
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: --input requires a directory path" >&2
                exit 1
            fi
            READS_DIR="$2"
            shift 2
            ;;
        --output)
            QC_OUT_DIR="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        --image)
            IMAGE_NAME="$2"
            shift 2
            ;;
        --cpus|--threads)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: --cpus requires an integer value" >&2
                exit 1
            fi
            THREADS="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option '$1'" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS/--cpus must be a positive integer (got '${THREADS}')" >&2
    exit 1
fi

if [ -f /.dockerenv ]; then
    BASEDIR="/data"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    if [ -z "${READS_DIR:-}" ] || [ -z "${QC_OUT_DIR:-}" ]; then
        usage >&2
        exit 1
    fi

    if [ ! -d "${READS_DIR}" ]; then
        echo "Reads directory not found: ${READS_DIR}" >&2
        echo "Set READS_DIR to the folder containing .fastq files." >&2
        exit 1
    fi
    READS_DIR="$(realpath "${READS_DIR}")"

    OUT_HOST="${QC_OUT_DIR}"
    mkdir -p "${OUT_HOST}"

    echo "Input: ${READS_DIR}"
    echo "Output: ${OUT_HOST}"

    exec docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e THREADS \
        -e JOBS \
        -e MIN_LENGTH \
        -e KEEP_PERCENT \
        -e CLEAN_FASTQ \
        -v "${READS_DIR}:/data/reads:ro" \
        -v "${OUT_HOST}:/data/qc_results" \
        -v "${SCRIPT_DIR}/run_qc.sh:/opt/pipeline/run_pipeline.sh:ro" \
        "$IMAGE_NAME" \
        bash /opt/pipeline/run_pipeline.sh
fi

READS_RAW="${BASEDIR}/reads"
OUT_DIR="${BASEDIR}/qc_results"

STAGE_CLEANED="${OUT_DIR}/0_cleaned"
STAGE_TRIMMED="${OUT_DIR}/1_trimmed"
STAGE_FILTERED="${OUT_DIR}/2_filtered"
STAGE_MULTIQC="${OUT_DIR}/3_multiqc"

STATS_FILE="${OUT_DIR}/pipeline_stats.tsv"

DROP_CLEANED="${STAGE_CLEANED}/dropped"
DROP_TRIMMED="${STAGE_TRIMMED}/dropped"
DROP_FILTERED="${STAGE_FILTERED}/dropped"

log_stage() {
    echo "$1"
}

log_step() {
    echo "$1"
}

log_sample() {
    echo "[$1/$2] $3"
}

count_reads() {
    local count
    count=$(seqkit stats -T "$1" 2>/dev/null | awk 'NR>1{print $4; exit}')
    if [[ -z "$count" || ! "$count" =~ ^[0-9]+$ ]]; then
        echo "Warning: could not parse read count for $1 (got '${count:-<empty>}'), treating as 0" >&2
        count=0
    fi
    echo "$count"
}

write_dropped_reads() {
    local infile="$1"
    local keptfile="$2"
    local dropfile="$3"
    local tmp tmpout

    if [ ! -s "$keptfile" ]; then
        mkdir -p "$(dirname "$dropfile")"
        seqkit seq "$infile" -o "$dropfile"
        return
    fi

    tmp="$(mktemp)"
    tmpout="$(mktemp)"
    seqkit seq -n "$keptfile" > "$tmp"
    seqkit grep -v -n -f "$tmp" "$infile" -o "$tmpout"
    if [ -s "$tmpout" ]; then
        mkdir -p "$(dirname "$dropfile")"
        mv "$tmpout" "$dropfile"
    else
        rm -f "$tmpout"
    fi
    rm -f "$tmp"
}

echo "Threads: ${THREADS}"
echo "Jobs: ${JOBS}"

SAMPLES=($(ls "${READS_RAW}"/*.fastq 2>/dev/null))
NUM_SAMPLES=${#SAMPLES[@]}

echo "Samples: ${NUM_SAMPLES}"

echo -e "sample\traw\tcleaned\ttrimmed\tfiltered" > "$STATS_FILE"

log_stage "STAGE 0-2"

mkdir -p "${STAGE_CLEANED}/reads"
mkdir -p "${STAGE_CLEANED}/nanoplot"
mkdir -p "${STAGE_CLEANED}/logs"
mkdir -p "${STAGE_TRIMMED}/reads"
mkdir -p "${STAGE_TRIMMED}/logs"
mkdir -p "${STAGE_TRIMMED}/nanoplot"
mkdir -p "${STAGE_FILTERED}/reads"
mkdir -p "${STAGE_FILTERED}/logs"
mkdir -p "${STAGE_FILTERED}/nanoplot"

log_step "Processing samples"

PIDS=()

run_with_limit() {
    while (( ${#PIDS[@]} >= JOBS )); do
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[i]}" 2>/dev/null; then
                unset 'PIDS[i]'
            fi
        done
        PIDS=("${PIDS[@]}")
        sleep 1
    done
    "$@" &
    PIDS+=($!)
}

wait_all() {
    local pid
    for pid in "${PIDS[@]}"; do
        wait "$pid"
    done
}

process_sample() {
    local fastq="$1"
    local sample outfile dropfile outdir
    local cleaned trimmed filtered raw_count_file
    local raw_count cleaned_count pass_count disc_lines dropped_reads

    sample="$(basename "$fastq")"
    sample="${sample%.fastq.gz}"
    sample="${sample%.fastq}"

    cleaned="${STAGE_CLEANED}/reads/${sample}.fastq"
    trimmed="${STAGE_TRIMMED}/reads/${sample}.fastq"
    filtered="${STAGE_FILTERED}/reads/${sample}.fastq"
    raw_count_file="${STAGE_CLEANED}/logs/${sample}.raw_count"

    # Stage 0: Clean
    if [ "$CLEAN_FASTQ" -eq 1 ]; then
        if [ ! -f "$cleaned" ]; then
            read pass_count disc_lines < <(
                seqkit sana "$fastq" -o "$cleaned" 2>&1 | awk '
                    {
                        if (match($0, /Pass records: [0-9]+/)) {
                            s=substr($0, RSTART, RLENGTH);
                            gsub(/[^0-9]/, "", s);
                            pass=s;
                        }
                        if (match($0, /Discarded lines: [0-9]+/)) {
                            s=substr($0, RSTART, RLENGTH);
                            gsub(/[^0-9]/, "", s);
                            disc=s;
                        }
                    }
                    END {
                        if (pass == "") pass=0;
                        if (disc == "") disc=0;
                        printf "%s %s\n", pass, disc;
                    }
                '
            )
            pass_count="${pass_count:-0}"
            disc_lines="${disc_lines:-0}"
            dropped_reads=$((disc_lines / 4))
            cleaned_count=$(count_reads "$cleaned")
            if [ "$pass_count" -eq 0 ]; then
                pass_count="$cleaned_count"
            fi
            raw_count=$((pass_count + dropped_reads))
            echo "$raw_count" > "$raw_count_file"
        fi
    else
        if [ ! -f "$cleaned" ]; then
            ln -sf "$fastq" "$cleaned"
        fi
        cleaned_count=$(count_reads "$cleaned")
        echo "$cleaned_count" > "$raw_count_file"
    fi

    if [ ! -f "$raw_count_file" ]; then
        if [ "$CLEAN_FASTQ" -eq 1 ]; then
            read pass_count disc_lines < <(
                seqkit sana "$fastq" -o /dev/null 2>&1 | awk '
                    {
                        if (match($0, /Pass records: [0-9]+/)) {
                            s=substr($0, RSTART, RLENGTH);
                            gsub(/[^0-9]/, "", s);
                            pass=s;
                        }
                        if (match($0, /Discarded lines: [0-9]+/)) {
                            s=substr($0, RSTART, RLENGTH);
                            gsub(/[^0-9]/, "", s);
                            disc=s;
                        }
                    }
                    END {
                        if (pass == "") pass=0;
                        if (disc == "") disc=0;
                        printf "%s %s\n", pass, disc;
                    }
                '
            )
            pass_count="${pass_count:-0}"
            disc_lines="${disc_lines:-0}"
            dropped_reads=$((disc_lines / 4))
            cleaned_count=$(count_reads "$cleaned")
            if [ "$pass_count" -eq 0 ]; then
                pass_count="$cleaned_count"
            fi
            raw_count=$((pass_count + dropped_reads))
            echo "$raw_count" > "$raw_count_file"
        else
            cleaned_count=$(count_reads "$cleaned")
            echo "$cleaned_count" > "$raw_count_file"
        fi
        echo "Warning: raw count recomputed for ${sample}" >&2
    fi

    outdir="${STAGE_CLEANED}/nanoplot/${sample}"
    if [ ! -f "${outdir}/${sample}_NanoPlot-report.html" ]; then
        NanoPlot --fastq "$cleaned" \
            --outdir "$outdir" \
            --prefix "${sample}_" \
            --threads "$THREADS" \
            --plots dot 2>&1 | grep -v "WARNING:root:Plotly\|Kaleido\|pip install" || true
    fi

    # Stage 1: Trim
    dropfile="${DROP_TRIMMED}/${sample}.fastq"
    if [ -f "$trimmed" ]; then
        if [ ! -f "$dropfile" ]; then
            write_dropped_reads "$cleaned" "$trimmed" "$dropfile"
        fi
    else
        porechop -i "$cleaned" \
            -o "$trimmed" \
            --threads "$THREADS" \
            --discard_middle \
            > "${STAGE_TRIMMED}/logs/${sample}.log" 2>&1
        write_dropped_reads "$cleaned" "$trimmed" "$dropfile"
    fi

    outdir="${STAGE_TRIMMED}/nanoplot/${sample}"
    if [ ! -f "${outdir}/${sample}_NanoPlot-report.html" ]; then
        NanoPlot --fastq "$trimmed" \
            --outdir "$outdir" \
            --prefix "${sample}_" \
            --threads "$THREADS" \
            --plots dot 2>&1 | grep -v "WARNING:root:Plotly\|Kaleido\|pip install" || true
    fi

    # Stage 2: Filter
    dropfile="${DROP_FILTERED}/${sample}.fastq"
    if [ -f "$filtered" ]; then
        if [ ! -f "$dropfile" ]; then
            write_dropped_reads "$trimmed" "$filtered" "$dropfile"
        fi
    else
        filtlong \
            --min_length "$MIN_LENGTH" \
            --keep_percent "$KEEP_PERCENT" \
            "$trimmed" > "$filtered" 2> "${STAGE_FILTERED}/logs/${sample}.log"
        write_dropped_reads "$trimmed" "$filtered" "$dropfile"
    fi

    outdir="${STAGE_FILTERED}/nanoplot/${sample}"
    if [ ! -f "${outdir}/${sample}_NanoPlot-report.html" ]; then
        NanoPlot --fastq "$filtered" \
            --outdir "$outdir" \
            --prefix "${sample}_" \
            --threads "$THREADS" \
            --plots dot 2>&1 | grep -v "WARNING:root:Plotly\|Kaleido\|pip install" || true
    fi
}

if [ "$JOBS" -le 1 ]; then
    for fastq in "${SAMPLES[@]}"; do
        process_sample "$fastq"
    done
else
    for fastq in "${SAMPLES[@]}"; do
        run_with_limit process_sample "$fastq"
    done
    wait_all
fi

log_stage "STAGE 3"

mkdir -p "${STAGE_MULTIQC}"

multiqc \
    "${STAGE_CLEANED}/nanoplot" \
    "${STAGE_TRIMMED}/nanoplot" \
    "${STAGE_FILTERED}/nanoplot" \
    --outdir "${STAGE_MULTIQC}" \
    --title "Nanopore QC Report" \
    --filename multiqc_report \
    --force

log_stage "Stats & Cleanup"

for fastq in "${SAMPLES[@]}"; do
    sample="$(basename "$fastq")"
    sample="${sample%.fastq.gz}"
    sample="${sample%.fastq}"

    cleaned="${STAGE_CLEANED}/reads/${sample}.fastq"
    trimmed="${STAGE_TRIMMED}/reads/${sample}.fastq"
    filtered="${STAGE_FILTERED}/reads/${sample}.fastq"

    raw_count=0
    cleaned_count=0
    trimmed_count=0
    filtered_count=0

    if [ -f "${STAGE_CLEANED}/logs/${sample}.raw_count" ]; then
        raw_count=$(cat "${STAGE_CLEANED}/logs/${sample}.raw_count")
    fi
    if [ -s "$cleaned" ]; then
        cleaned_count=$(count_reads "$cleaned")
    fi
    if [ -s "$trimmed" ]; then
        trimmed_count=$(count_reads "$trimmed")
    fi
    if [ -s "$filtered" ]; then
        filtered_count=$(count_reads "$filtered")
    fi

    echo -e "${sample}\t${raw_count}\t${cleaned_count}\t${trimmed_count}\t${filtered_count}" >> "$STATS_FILE"
done

head -1 "$STATS_FILE" > "${STATS_FILE}.tmp"
tail -n +2 "$STATS_FILE" | sort >> "${STATS_FILE}.tmp"
mv "${STATS_FILE}.tmp" "$STATS_FILE"

log_step "Cleanup"
rm -rf "${STAGE_CLEANED}/reads"
rm -rf "${STAGE_TRIMMED}/reads"

echo "Pipeline complete"
echo "Read counts per stage:"
cat "$STATS_FILE"
echo "Output (container path): ${OUT_DIR}"
echo "Final reads (container path): ${OUT_DIR}/2_filtered/reads/*.fastq"
echo "Dropped reads (container path): ${OUT_DIR}/*/dropped/*.fastq"
echo "QC report (container path): ${OUT_DIR}/3_multiqc/multiqc_report.html"
echo "Stats (container path): ${OUT_DIR}/pipeline_stats.tsv"
