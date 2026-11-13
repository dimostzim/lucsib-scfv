#!/usr/bin/env python3
"""Try to resolve ambiguous TIM-3 labels using full Nanopore read context.

The first-pass labeler compares Sanger sequences to Nanopore nt_scFv only.
For ambiguous clones, this script compares the best yes/no Sanger references
against the full Nanopore `sequence` field and aggregates read support.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

from common import Alignment, align, read_tsv, revcomp, write_tsv


# minimum first-pass TIM-3 reads of an ambiguous clone to analyse it
MIN_READS = 10
# minimum identity of the yes Sanger clone to the full reads for a bind_context call
MIN_YES_FULL_IDENTITY = 91.0
# a read favours the yes (no) Sanger clone when its identity is higher (lower) by more than this margin
MARGIN = 0.5


@dataclass
class CloneAggregate:
    nt_scfv_id: str
    best_yes_clone: str
    best_no_clone: str
    first_pass_reads: int
    source_rows: int = 0
    observed_reads: int = 0
    margin05_reads: Counter = field(default_factory=Counter)
    weighted_delta_sum: float = 0.0
    best_yes_identity: float = -1.0
    best_no_identity: float = -1.0
    best_yes_cov: float = 0.0
    best_no_cov: float = 0.0


def align_sanger_to_full_read(sanger_seq: str, full_seq: str) -> Alignment:
    best = None
    for query in (sanger_seq, revcomp(sanger_seq)):
        aln = align(query, full_seq, mode="HW")
        if best is None or (aln.identity, aln.query_cov, -aln.edit_distance) > (
            best.identity,
            best.query_cov,
            -best.edit_distance,
        ):
            best = aln
    return best


def call_by_delta(delta: float) -> str:
    if delta > MARGIN:
        return "bind"
    if delta < -MARGIN:
        return "no_bind"
    return "still_ambiguous"


def main() -> None:
    parser = argparse.ArgumentParser(description="Full-context resolution for ambiguous TIM-3 clones.")
    parser.add_argument("--sanger-tsv", type=Path, required=True)
    parser.add_argument("--counts-tsv", type=Path, required=True, help="step-03 counts table")
    parser.add_argument("--labels-tsv", type=Path, required=True, help="tim3_conservative_labels.tsv")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    sanger = {row["clone_number"]: row["sequence"] for row in read_tsv(args.sanger_tsv)}
    labels = read_tsv(args.labels_tsv)
    aggregates: Dict[str, CloneAggregate] = {
        row["nt_scfv_id"]: CloneAggregate(
            nt_scfv_id=row["nt_scfv_id"],
            best_yes_clone=row["best_yes_clone"],
            best_no_clone=row["best_no_clone"],
            first_pass_reads=int(row["tim3_reads"]),
        )
        for row in labels
        if row["assigned_label"] == "ambiguous" and int(row["tim3_reads"]) >= MIN_READS
    }
    print(f"Analyzing {len(aggregates)} ambiguous nt_scFv clones with >= {MIN_READS} reads", flush=True)

    with args.counts_tsv.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            agg = aggregates.get(row["nt_scFv_id"])
            if agg is None:
                continue
            reads = int(row["barcode02"]) + int(row["barcode03"])
            if reads == 0:
                continue
            yes = align_sanger_to_full_read(sanger[agg.best_yes_clone], row["sequence"])
            no = align_sanger_to_full_read(sanger[agg.best_no_clone], row["sequence"])
            agg.source_rows += 1
            agg.observed_reads += reads
            delta = yes.identity - no.identity
            agg.weighted_delta_sum += delta * reads
            agg.margin05_reads[call_by_delta(delta)] += reads
            if yes.identity > agg.best_yes_identity:
                agg.best_yes_identity = yes.identity
                agg.best_yes_cov = yes.query_cov
            if no.identity > agg.best_no_identity:
                agg.best_no_identity = no.identity
                agg.best_no_cov = no.query_cov

    rows = []
    for agg in sorted(aggregates.values(), key=lambda x: (-x.first_pass_reads, x.nt_scfv_id)):
        dominant_call, dominant_reads = agg.margin05_reads.most_common(1)[0]
        rows.append(
            {
                "nt_scfv_id": agg.nt_scfv_id,
                "best_yes_clone": agg.best_yes_clone,
                "best_no_clone": agg.best_no_clone,
                "first_pass_reads": agg.first_pass_reads,
                "observed_reads": agg.observed_reads,
                "source_rows": agg.source_rows,
                "mean_full_identity_delta_yes_minus_no": round(agg.weighted_delta_sum / agg.observed_reads, 4),
                "best_yes_full_identity": round(agg.best_yes_identity, 4),
                "best_no_full_identity": round(agg.best_no_identity, 4),
                "best_yes_full_cov": round(agg.best_yes_cov, 4),
                "best_no_full_cov": round(agg.best_no_cov, 4),
                "margin05_bind_reads": agg.margin05_reads["bind"],
                "margin05_no_bind_reads": agg.margin05_reads["no_bind"],
                "margin05_still_ambiguous_reads": agg.margin05_reads["still_ambiguous"],
                "margin05_dominant_call": dominant_call,
                "margin05_dominant_fraction": round(dominant_reads / agg.observed_reads, 4),
            }
        )
    write_tsv(args.out_dir / f"tim3_ambiguous_full_context_minreads{MIN_READS}.tsv", rows)

    # Final labels: the first-pass label, except that ambiguous clones whose reads mostly favour the
    # yes Sanger clone, with the yes clone aligning to the full reads at >= MIN_YES_FULL_IDENTITY,
    # become bind_context.
    context_rows = {row["nt_scfv_id"]: row for row in rows}
    first_pass_confidence = {"bind": "first_pass", "no_bind": "first_pass", "ambiguous": "ambiguous", "unknown": "unknown"}
    for row in labels:
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
    write_tsv(args.out_dir / f"tim3_labels_with_full_context_minreads{MIN_READS}.tsv", labels)


if __name__ == "__main__":
    main()
