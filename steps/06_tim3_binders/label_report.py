#!/usr/bin/env python3
"""Create final TIM-3 label summaries, FASTA exports, and report."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence

from common import read_tsv, write_tsv


DIRECT_LABELS = ("bind", "bind_context", "no_bind", "ambiguous", "unknown")
FASTA_LABELS = ("bind", "bind_context", "no_bind")


def write_fasta(path: Path, rows: Sequence[dict], seq_field: str, label: str) -> None:
    with path.open("w") as handle:
        for row in rows:
            header = (
                f"{row['nt_scfv_id']}|label={label}|aa_id={row['aa_scfv_id']}|"
                f"barcode02={row['barcode02']}|barcode03={row['barcode03']}|"
                f"tim3_reads={row['tim3_reads']}|best_yes={row['best_yes_clone']}|"
                f"best_no={row['best_no_clone']}"
            )
            handle.write(f">{header}\n")
            seq = row[seq_field]
            for i in range(0, len(seq), 80):
                handle.write(f"{seq[i : i + 80]}\n")


def load_joined_rows(labels_tsv: Path, collapsed_tsv: Path) -> List[dict]:
    collapsed = {row["nt_scfv_id"]: row for row in read_tsv(collapsed_tsv)}
    rows = []
    for label_row in read_tsv(labels_tsv):
        seq_row = collapsed[label_row["nt_scfv_id"]]
        row = dict(label_row)
        row["nt_scfv"] = seq_row["nt_scfv"]
        row["aa_scfv"] = seq_row["aa_scfv"]
        row["nt_scfv_len"] = len(seq_row["nt_scfv"])
        row["aa_scfv_len"] = len(seq_row["aa_scfv"])
        row["barcode02"] = int(row["barcode02"])
        row["barcode03"] = int(row["barcode03"])
        row["tim3_reads"] = int(row["tim3_reads"])
        rows.append(row)
    return sorted(rows, key=lambda r: (-r["tim3_reads"], r["nt_scfv_id"]))


def make_label_summary(rows: Sequence[dict]) -> List[dict]:
    labels = sorted({row["final_label_margin05_context"] for row in rows})
    summary = []
    total_clones = len(rows)
    total_reads = sum(row["tim3_reads"] for row in rows)
    total_b02 = sum(row["barcode02"] for row in rows)
    total_b03 = sum(row["barcode03"] for row in rows)
    for label in labels:
        sub = [row for row in rows if row["final_label_margin05_context"] == label]
        b02 = sum(row["barcode02"] for row in sub)
        b03 = sum(row["barcode03"] for row in sub)
        reads = b02 + b03
        summary.append(
            {
                "label": label,
                "clones": len(sub),
                "barcode02_reads": b02,
                "barcode03_reads": b03,
                "tim3_reads": reads,
                "clone_fraction": round(len(sub) / total_clones, 6),
                "read_fraction": round(reads / total_reads, 6),
                "barcode02_read_fraction": round(b02 / total_b02, 6),
                "barcode03_read_fraction": round(b03 / total_b03, 6),
            }
        )
    for tier, tier_labels in (("strict_bind_positive", {"bind"}), ("expanded_bind_positive", {"bind", "bind_context"})):
        sub = [row for row in rows if row["final_label_margin05_context"] in tier_labels]
        summary.append(
            {
                "label": tier,
                "clones": len(sub),
                "barcode02_reads": sum(row["barcode02"] for row in sub),
                "barcode03_reads": sum(row["barcode03"] for row in sub),
                "tim3_reads": sum(row["tim3_reads"] for row in sub),
                "clone_fraction": "",
                "read_fraction": "",
                "barcode02_read_fraction": "",
                "barcode03_read_fraction": "",
            }
        )
    return summary


def make_top_clones(rows: Sequence[dict], top_n: int = 25) -> List[dict]:
    # rows are sorted by reads, so the first rows of each label are its top clones
    top_rows = []
    for label in DIRECT_LABELS:
        sub = [row for row in rows if row["final_label_margin05_context"] == label]
        for rank, row in enumerate(sub[:top_n], start=1):
            top_rows.append(
                {
                    "label": label,
                    "rank": rank,
                    "nt_scfv_id": row["nt_scfv_id"],
                    "aa_scfv_id": row["aa_scfv_id"],
                    "barcode02": row["barcode02"],
                    "barcode03": row["barcode03"],
                    "tim3_reads": row["tim3_reads"],
                    "nt_scfv_len": row["nt_scfv_len"],
                    "aa_scfv_len": row["aa_scfv_len"],
                    "best_yes_clone": row["best_yes_clone"],
                    "best_yes_identity": row["best_yes_identity"],
                    "best_no_clone": row["best_no_clone"],
                    "best_no_identity": row["best_no_identity"],
                    "full_context_mean_delta_yes_minus_no": row["full_context_mean_delta_yes_minus_no"],
                    "confidence": row["confidence"],
                }
            )
    return top_rows


def make_conflict_audit(sanger_pairs_tsv: Path) -> List[dict]:
    rows = []
    for row in read_tsv(sanger_pairs_tsv):
        if row["same_label"] != "False":
            continue
        identity = float(row["identity"])
        if identity < 98.0 or min(float(row["clone_a_cov"]), float(row["clone_b_cov"])) < 90.0:
            continue
        rows.append(
            {
                "clone_a": row["clone_a"],
                "binding_a": row["binding_a"],
                "len_a": row["len_a"],
                "clone_b": row["clone_b"],
                "binding_b": row["binding_b"],
                "len_b": row["len_b"],
                "identity": row["identity"],
                "clone_a_cov": row["clone_a_cov"],
                "clone_b_cov": row["clone_b_cov"],
                "shorter_cov": row["shorter_cov"],
                "severity": "high" if identity >= 99.0 else "moderate",
            }
        )
    return sorted(rows, key=lambda r: (-float(r["identity"]), r["clone_a"], r["clone_b"]))


def make_label_matrix(rows: Sequence[dict]) -> List[dict]:
    # Compact barcode-level matrix for plotting.
    output = []
    for row in rows:
        output.append(
            {
                "nt_scfv_id": row["nt_scfv_id"],
                "aa_scfv_id": row["aa_scfv_id"],
                "label": row["final_label_margin05_context"],
                "strict_bind_positive": int(row["final_label_margin05_context"] == "bind"),
                "expanded_bind_positive": int(row["final_label_margin05_context"] in {"bind", "bind_context"}),
                "barcode02": row["barcode02"],
                "barcode03": row["barcode03"],
                "tim3_reads": row["tim3_reads"],
                "best_yes_clone": row["best_yes_clone"],
                "best_no_clone": row["best_no_clone"],
                "confidence": row["confidence"],
            }
        )
    return output


def make_full_table(rows: Sequence[dict]) -> List[dict]:
    # One row per nt_scFv with its final label and sequences.
    return [
        {
            "nt_scfv_id": row["nt_scfv_id"],
            "aa_scfv_id": row["aa_scfv_id"],
            "final_label": row["final_label_margin05_context"],
            "confidence": row["confidence"],
            "total_tim3_reads": row["tim3_reads"],
            "barcode02_tim3_round3_reads": row["barcode02"],
            "barcode03_tim3_round4_reads": row["barcode03"],
            "best_yes_sanger_clone": row["best_yes_clone"],
            "best_no_sanger_clone": row["best_no_clone"],
            "nt_scfv_sequence": row["nt_scfv"],
            "aa_scfv_sequence": row["aa_scfv"],
        }
        for row in rows
    ]


def write_exports(rows: Sequence[dict], fasta_dir: Path) -> Dict[str, int]:
    sets = {label: [row for row in rows if row["final_label_margin05_context"] == label] for label in FASTA_LABELS}
    sets["expanded_bind_positive"] = sets["bind"] + sets["bind_context"]

    counts = {}
    for label, sub in sets.items():
        write_fasta(fasta_dir / f"tim3_{label}_nt_scfv.fasta", sub, "nt_scfv", label)
        write_fasta(fasta_dir / f"tim3_{label}_aa_scfv.fasta", sub, "aa_scfv", label)
        counts[label] = len(sub)
    return counts


def markdown_table(rows: Sequence[dict], columns: Sequence[str]) -> str:
    lines = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("|" + "|".join("---" for _ in columns) + "|")
    for row in rows:
        lines.append("| " + " | ".join(str(row[col]) for col in columns) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    summary: Sequence[dict],
    top_rows: Sequence[dict],
    conflicts: Sequence[dict],
    fasta_counts: Dict[str, int],
) -> None:
    top_bind_like = [
        row for row in top_rows if row["label"] in {"bind", "bind_context", "no_bind", "ambiguous"}
    ][:20]
    direct_summary = [row for row in summary if row["label"] in DIRECT_LABELS]
    report = f"""# TIM-3 Label Report

## Label Rules

- `bind`: direct high-confidence `nt_scFv` alignment to a validated Sanger `yes` clone.
- `bind_context`: `nt_scFv` was ambiguous (equally close to a `yes` and a `no` Sanger clone); over the full Nanopore reads, most reads match the `yes` clone better by >0.5 percentage points and the `yes` clone aligns to the full reads with >=91% identity. The competing Sanger clones have identical scFv sequences, so this call rests on sequence outside the scFv.
- `no_bind`: direct high-confidence `nt_scFv` alignment to a validated Sanger `no` clone.
- `ambiguous`: unresolved `yes`/`no` conflict.
- `unknown`: no strong validated Sanger match.

For downstream analysis, use:

- Binder evidence: `bind` only (strict bind-positive).
- Non-binders: `no_bind`.
- Not evidence either way: `ambiguous`, `unknown` and `bind_context`, and therefore also the expanded set
  (`bind + bind_context`, `expanded_bind_positive`).

## Read Summary

{markdown_table(direct_summary, ["label", "clones", "barcode02_reads", "barcode03_reads", "tim3_reads", "read_fraction"])}

## Combined Positive Tiers

{markdown_table([row for row in summary if row["label"] in {"strict_bind_positive", "expanded_bind_positive"}], ["label", "clones", "barcode02_reads", "barcode03_reads", "tim3_reads"])}

## Top Clones

{markdown_table(top_bind_like, ["label", "rank", "nt_scfv_id", "barcode02", "barcode03", "tim3_reads", "best_yes_clone", "best_no_clone", "confidence"])}

## Sanger Conflict Audit

There are {len(conflicts)} opposite-label Sanger pairs with >=98% identity and >=90% coverage. These conflicts explain why some Nanopore clones cannot be labeled from `nt_scFv` alone.

Top conflicts:

{markdown_table(conflicts[:20], ["clone_a", "binding_a", "clone_b", "binding_b", "identity", "clone_a_cov", "clone_b_cov", "severity"])}

## FASTA Exports

- `tim3_bind_nt_scfv.fasta` / `tim3_bind_aa_scfv.fasta`: {fasta_counts["bind"]} clones.
- `tim3_bind_context_nt_scfv.fasta` / `tim3_bind_context_aa_scfv.fasta`: {fasta_counts["bind_context"]} clones.
- `tim3_no_bind_nt_scfv.fasta` / `tim3_no_bind_aa_scfv.fasta`: {fasta_counts["no_bind"]} clones.
- `tim3_expanded_bind_positive_*`: {fasta_counts["expanded_bind_positive"]} clones.

## Files

- `label_read_summary.tsv`
- `top_clones_by_label.tsv`
- `sanger_conflict_audit.tsv`
- `label_matrix_for_plotting.tsv`
- `fasta/`
"""
    path.write_text(report.strip() + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="TIM-3 label report and FASTA exports.")
    parser.add_argument("--analysis-dir", type=Path, required=True, help="output dir of the matching steps")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out_dir = args.out_dir
    fasta_dir = out_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    rows = load_joined_rows(
        args.analysis_dir / "tim3_labels_with_full_context_minreads10.tsv",
        args.analysis_dir / "tim3_nanopore_nt_scfv_collapsed.tsv",
    )
    summary = make_label_summary(rows)
    top_rows = make_top_clones(rows)
    conflicts = make_conflict_audit(args.analysis_dir / "tim3_sanger_vs_sanger.tsv")
    fasta_counts = write_exports(rows, fasta_dir)

    write_tsv(out_dir / "label_read_summary.tsv", summary)
    write_tsv(out_dir / "top_clones_by_label.tsv", top_rows)
    write_tsv(out_dir / "sanger_conflict_audit.tsv", conflicts)
    write_tsv(out_dir / "label_matrix_for_plotting.tsv", make_label_matrix(rows))
    write_tsv(out_dir / "tim3_full_labeled_scfv_results_with_sequences.tsv", make_full_table(rows))
    write_report(out_dir / "TIM3_LABEL_REPORT.md", summary, top_rows, conflicts, fasta_counts)
    print(f"Wrote final TIM-3 exports to {out_dir}")


if __name__ == "__main__":
    main()
