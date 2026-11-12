#!/usr/bin/env python3
"""Try to resolve ambiguous TIM-3 labels using full Nanopore read context.

The first-pass labeler compares Sanger sequences to Nanopore nt_scFv only.
For ambiguous clones, this script compares the best yes/no Sanger references
against the full Nanopore `sequence` field and aggregates read support.
"""

from __future__ import annotations

import argparse
import csv
import sys
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

try:
    import edlib
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency: edlib. Use the Python environment from envs/python.yml."
    ) from exc


# Set from the command line in main()
SANGER_TSV: Path
COUNTS_TSV: Path
LABELS_TSV: Path
OUT_DIR: Path

# minimum identity of the yes Sanger clone to the full reads for a bind_context call
MIN_YES_FULL_IDENTITY = 91.0

# full reads can be longer than the csv module's default field limit
csv.field_size_limit(sys.maxsize)

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_RE = re.compile(r"(\d+)([=XID])")


@dataclass
class LabelRecord:
    nt_scfv_id: str
    best_yes_clone: str
    best_no_clone: str
    tim3_reads: int


@dataclass
class AlignMetrics:
    identity: float
    edit_distance: int
    query_cov: float
    aligned_len: int
    strand: str


@dataclass
class CloneAggregate:
    nt_scfv_id: str
    best_yes_clone: str
    best_no_clone: str
    first_pass_reads: int
    source_rows: int = 0
    observed_reads: int = 0
    margin05_reads: Counter = field(default_factory=Counter)
    margin10_reads: Counter = field(default_factory=Counter)
    margin05_rows: Counter = field(default_factory=Counter)
    margin10_rows: Counter = field(default_factory=Counter)
    weighted_delta_sum: float = 0.0
    best_yes_identity: float = -1.0
    best_no_identity: float = -1.0
    best_yes_cov: float = 0.0
    best_no_cov: float = 0.0

    @property
    def mean_delta(self) -> float:
        if self.observed_reads == 0:
            return 0.0
        return self.weighted_delta_sum / self.observed_reads


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Full-context resolution for ambiguous TIM-3 clones.")
    parser.add_argument("--min-reads", type=int, default=10, help="Minimum ambiguous TIM-3 reads to analyze.")
    parser.add_argument(
        "--output-prefix",
        default="tim3_ambiguous_full_context_minreads",
        help="Output filename prefix under --out-dir.",
    )
    parser.add_argument("--sanger-tsv", type=Path, required=True)
    parser.add_argument("--counts-tsv", type=Path, required=True, help="step-03 counts table")
    parser.add_argument("--labels-tsv", type=Path, required=True, help="tim3_conservative_labels.tsv")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def parse_int(value: str) -> int:
    if not value:
        return 0
    return int(float(value))


def parse_cigar(cigar: str) -> Tuple[int, int, int, int, int, int]:
    matches = mismatches = insertions = deletions = q_consumed = 0
    for count_s, op in CIGAR_RE.findall(cigar):
        count = int(count_s)
        if op == "=":
            matches += count
            q_consumed += count
        elif op == "X":
            mismatches += count
            q_consumed += count
        elif op == "I":
            insertions += count
            q_consumed += count
        elif op == "D":
            deletions += count
    aligned_len = matches + mismatches + insertions + deletions
    return matches, mismatches, insertions, deletions, q_consumed, aligned_len


def align_sanger_to_full_read(sanger_seq: str, full_seq: str) -> Optional[AlignMetrics]:
    best: Optional[AlignMetrics] = None
    for strand, query in (("+", sanger_seq), ("-", revcomp(sanger_seq))):
        result = edlib.align(query, full_seq, mode="HW", task="path")
        if result["editDistance"] < 0 or not result.get("cigar"):
            continue
        matches, mismatches, ins, dels, q_consumed, aligned_len = parse_cigar(result["cigar"])
        if aligned_len == 0:
            continue
        metrics = AlignMetrics(
            identity=100.0 * matches / aligned_len,
            edit_distance=int(result["editDistance"]),
            query_cov=100.0 * q_consumed / len(query),
            aligned_len=aligned_len,
            strand=strand,
        )
        key = (metrics.identity, metrics.query_cov, -metrics.edit_distance)
        best_key = (best.identity, best.query_cov, -best.edit_distance) if best else None
        if best is None or key > best_key:
            best = metrics
    return best


def load_sanger() -> Dict[str, str]:
    refs = {}
    with SANGER_TSV.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            refs[row["clone_number"]] = row["sequence"].strip().upper()
    return refs


def load_ambiguous_labels(min_reads: int) -> Dict[str, LabelRecord]:
    selected = {}
    with LABELS_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reads = int(row["tim3_reads"])
            if row["assigned_label"] != "ambiguous" or reads < min_reads:
                continue
            if not row["best_yes_clone"] or not row["best_no_clone"]:
                continue
            selected[row["nt_scfv_id"]] = LabelRecord(
                nt_scfv_id=row["nt_scfv_id"],
                best_yes_clone=row["best_yes_clone"],
                best_no_clone=row["best_no_clone"],
                tim3_reads=reads,
            )
    return selected


def call_by_delta(delta: float, margin: float) -> str:
    if delta > margin:
        return "bind"
    if delta < -margin:
        return "no_bind"
    return "still_ambiguous"


def main() -> None:
    global SANGER_TSV, COUNTS_TSV, LABELS_TSV, OUT_DIR
    args = parse_args()
    SANGER_TSV, COUNTS_TSV, LABELS_TSV, OUT_DIR = args.sanger_tsv, args.counts_tsv, args.labels_tsv, args.out_dir
    sanger = load_sanger()
    labels = load_ambiguous_labels(args.min_reads)
    aggregates = {
        nt_id: CloneAggregate(
            nt_scfv_id=record.nt_scfv_id,
            best_yes_clone=record.best_yes_clone,
            best_no_clone=record.best_no_clone,
            first_pass_reads=record.tim3_reads,
        )
        for nt_id, record in labels.items()
    }

    print(
        f"Analyzing {len(labels)} ambiguous nt_scFv clones with >= {args.min_reads} reads",
        flush=True,
    )

    with COUNTS_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for i, row in enumerate(reader, start=1):
            nt_id = row.get("nt_scFv_id", "")
            record = labels.get(nt_id)
            if record is None:
                continue
            reads = parse_int(row.get("barcode02", "0")) + parse_int(row.get("barcode03", "0"))
            if reads == 0:
                continue
            full_seq = (row.get("sequence") or "").strip().upper()
            yes = align_sanger_to_full_read(sanger[record.best_yes_clone], full_seq)
            no = align_sanger_to_full_read(sanger[record.best_no_clone], full_seq)
            if yes is None or no is None:
                continue

            agg = aggregates[nt_id]
            agg.source_rows += 1
            agg.observed_reads += reads
            delta = yes.identity - no.identity
            agg.weighted_delta_sum += delta * reads
            agg.margin05_reads[call_by_delta(delta, 0.5)] += reads
            agg.margin10_reads[call_by_delta(delta, 1.0)] += reads
            agg.margin05_rows[call_by_delta(delta, 0.5)] += 1
            agg.margin10_rows[call_by_delta(delta, 1.0)] += 1
            if yes.identity > agg.best_yes_identity:
                agg.best_yes_identity = yes.identity
                agg.best_yes_cov = yes.query_cov
            if no.identity > agg.best_no_identity:
                agg.best_no_identity = no.identity
                agg.best_no_cov = no.query_cov

    rows = []
    for agg in sorted(aggregates.values(), key=lambda x: (-x.first_pass_reads, x.nt_scfv_id)):
        dom05, dom05_reads = agg.margin05_reads.most_common(1)[0] if agg.margin05_reads else ("not_seen", 0)
        dom10, dom10_reads = agg.margin10_reads.most_common(1)[0] if agg.margin10_reads else ("not_seen", 0)
        rows.append(
            {
                "nt_scfv_id": agg.nt_scfv_id,
                "best_yes_clone": agg.best_yes_clone,
                "best_no_clone": agg.best_no_clone,
                "first_pass_reads": agg.first_pass_reads,
                "observed_reads": agg.observed_reads,
                "source_rows": agg.source_rows,
                "mean_full_identity_delta_yes_minus_no": round(agg.mean_delta, 4),
                "best_yes_full_identity": round(agg.best_yes_identity, 4),
                "best_no_full_identity": round(agg.best_no_identity, 4),
                "best_yes_full_cov": round(agg.best_yes_cov, 4),
                "best_no_full_cov": round(agg.best_no_cov, 4),
                "margin05_bind_reads": agg.margin05_reads["bind"],
                "margin05_no_bind_reads": agg.margin05_reads["no_bind"],
                "margin05_still_ambiguous_reads": agg.margin05_reads["still_ambiguous"],
                "margin05_dominant_call": dom05,
                "margin05_dominant_fraction": round(dom05_reads / agg.observed_reads, 4)
                if agg.observed_reads
                else 0.0,
                "margin10_bind_reads": agg.margin10_reads["bind"],
                "margin10_no_bind_reads": agg.margin10_reads["no_bind"],
                "margin10_still_ambiguous_reads": agg.margin10_reads["still_ambiguous"],
                "margin10_dominant_call": dom10,
                "margin10_dominant_fraction": round(dom10_reads / agg.observed_reads, 4)
                if agg.observed_reads
                else 0.0,
            }
        )

    out_tsv = OUT_DIR / f"{args.output_prefix}{args.min_reads}.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    totals05 = Counter()
    totals10 = Counter()
    for agg in aggregates.values():
        totals05.update(agg.margin05_reads)
        totals10.update(agg.margin10_reads)

    out_summary = OUT_DIR / f"{args.output_prefix}{args.min_reads}_summary.md"
    out_summary.write_text(
        "\n".join(
            [
                "# TIM-3 Ambiguous Full-Context Resolution",
                "",
                f"- Input ambiguous clones with >= {args.min_reads} reads: {len(labels)}",
                f"- First-pass ambiguous reads represented: {sum(r.tim3_reads for r in labels.values())}",
                f"- Observed reads in counts table: {sum(a.observed_reads for a in aggregates.values())}",
                "",
                "## Full-Context Delta Calls",
                "Delta is `full-context yes-reference identity - no-reference identity`.",
                "",
                "Margin 0.5 percentage points:",
                f"- bind reads: {totals05['bind']}",
                f"- no_bind reads: {totals05['no_bind']}",
                f"- still_ambiguous reads: {totals05['still_ambiguous']}",
                "",
                "Margin 1.0 percentage points:",
                f"- bind reads: {totals10['bind']}",
                f"- no_bind reads: {totals10['no_bind']}",
                f"- still_ambiguous reads: {totals10['still_ambiguous']}",
                "",
                f"Detailed table: `{out_tsv}`",
                "",
            ]
        )
        + "\n"
    )

    # Final labels: the first-pass label, except that ambiguous clones whose reads mostly favour the
    # yes Sanger clone at the 0.5 margin, with the yes clone aligning to the full reads at
    # >= MIN_YES_FULL_IDENTITY, become bind_context.
    context_rows = {row["nt_scfv_id"]: row for row in rows}
    first_pass_confidence = {"bind": "first_pass", "no_bind": "first_pass", "ambiguous": "ambiguous", "unknown": "unknown"}
    with LABELS_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        label_fields = list(reader.fieldnames) + [
            "final_label_margin05_context",
            "full_context_mean_delta_yes_minus_no",
            "confidence",
        ]
        final_rows = []
        for row in reader:
            context = context_rows.get(row["nt_scfv_id"])
            if context is None:
                row["final_label_margin05_context"] = row["assigned_label"]
                row["full_context_mean_delta_yes_minus_no"] = ""
                row["confidence"] = first_pass_confidence[row["assigned_label"]]
            elif (
                context["margin05_dominant_call"] == "bind"
                and context["best_yes_full_identity"] >= MIN_YES_FULL_IDENTITY
            ):
                row["final_label_margin05_context"] = "bind_context"
                row["full_context_mean_delta_yes_minus_no"] = context["mean_full_identity_delta_yes_minus_no"]
                row["confidence"] = "context_tentative_margin05"
            else:
                row["final_label_margin05_context"] = row["assigned_label"]
                row["full_context_mean_delta_yes_minus_no"] = context["mean_full_identity_delta_yes_minus_no"]
                row["confidence"] = "ambiguous_after_context"
            final_rows.append(row)
    out_labels = OUT_DIR / f"tim3_labels_with_full_context_minreads{args.min_reads}.tsv"
    with out_labels.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=label_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(final_rows)

    print(f"Wrote {out_tsv}", flush=True)
    print(f"Wrote {out_summary}", flush=True)
    print(f"Wrote {out_labels}", flush=True)


if __name__ == "__main__":
    main()
