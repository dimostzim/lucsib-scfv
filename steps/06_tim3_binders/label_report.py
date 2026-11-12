#!/usr/bin/env python3
"""Create final TIM-3 label summaries, FASTA exports, and report."""

from __future__ import annotations

import csv
import textwrap
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


ROOT = Path("/homes/dtzim01/lucsib")
ANALYSIS_DIR = ROOT / "analysis/tim3_sanger_cutoff"
LABELS_TSV = ANALYSIS_DIR / "tim3_labels_with_full_context_minreads10.tsv"
COLLAPSED_TSV = ANALYSIS_DIR / "tim3_nanopore_nt_scfv_collapsed.tsv"
SANGER_PAIRS_TSV = ANALYSIS_DIR / "tim3_sanger_vs_sanger.tsv"
SANGER_TSV = ROOT / "TIM-3_sequences_binding.tsv"
OUT_DIR = ANALYSIS_DIR / "final_exports"
FASTA_DIR = OUT_DIR / "fasta"

DIRECT_LABELS = ("bind", "bind_context", "no_bind", "ambiguous", "unknown")
FASTA_LABELS = ("bind", "bind_context", "no_bind")


def read_tsv(path: Path) -> List[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and fieldnames is None:
        path.write_text("")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_fasta(path: Path, rows: Sequence[dict], seq_field: str, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            seq = row.get(seq_field, "")
            if not seq:
                continue
            header = (
                f"{row['nt_scfv_id']}|label={label}|aa_id={row['aa_scfv_id']}|"
                f"barcode02={row['barcode02']}|barcode03={row['barcode03']}|"
                f"tim3_reads={row['tim3_reads']}|best_yes={row['best_yes_clone']}|"
                f"best_no={row['best_no_clone']}"
            )
            handle.write(f">{header}\n")
            for chunk in wrap_fasta(seq):
                handle.write(f"{chunk}\n")


def load_joined_rows() -> List[dict]:
    labels = {row["nt_scfv_id"]: row for row in read_tsv(LABELS_TSV)}
    collapsed = {row["nt_scfv_id"]: row for row in read_tsv(COLLAPSED_TSV)}
    rows = []
    for nt_id, label_row in labels.items():
        seq_row = collapsed.get(nt_id)
        if seq_row is None:
            continue
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
                "barcode02_read_fraction": round(b02 / total_b02, 6) if total_b02 else 0,
                "barcode03_read_fraction": round(b03 / total_b03, 6) if total_b03 else 0,
            }
        )
    summary.append(
        {
            "label": "strict_bind_positive",
            "clones": sum(1 for row in rows if row["final_label_margin05_context"] == "bind"),
            "barcode02_reads": sum(row["barcode02"] for row in rows if row["final_label_margin05_context"] == "bind"),
            "barcode03_reads": sum(row["barcode03"] for row in rows if row["final_label_margin05_context"] == "bind"),
            "tim3_reads": sum(row["tim3_reads"] for row in rows if row["final_label_margin05_context"] == "bind"),
            "clone_fraction": "",
            "read_fraction": "",
            "barcode02_read_fraction": "",
            "barcode03_read_fraction": "",
        }
    )
    expanded = [row for row in rows if row["final_label_margin05_context"] in {"bind", "bind_context"}]
    summary.append(
        {
            "label": "expanded_bind_positive",
            "clones": len(expanded),
            "barcode02_reads": sum(row["barcode02"] for row in expanded),
            "barcode03_reads": sum(row["barcode03"] for row in expanded),
            "tim3_reads": sum(row["tim3_reads"] for row in expanded),
            "clone_fraction": "",
            "read_fraction": "",
            "barcode02_read_fraction": "",
            "barcode03_read_fraction": "",
        }
    )
    return summary


def make_top_clones(rows: Sequence[dict], top_n: int = 25) -> List[dict]:
    top_rows = []
    for label in DIRECT_LABELS:
        sub = [row for row in rows if row["final_label_margin05_context"] == label]
        for rank, row in enumerate(sorted(sub, key=lambda r: (-r["tim3_reads"], r["nt_scfv_id"]))[:top_n], start=1):
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


def make_conflict_audit() -> List[dict]:
    sanger_labels = {
        row["clone_number"]: row["binding"]
        for row in read_tsv(SANGER_TSV)
    }
    rows = []
    for row in read_tsv(SANGER_PAIRS_TSV):
        if row["same_label"] != "False":
            continue
        identity = float(row["identity"])
        clone_a_cov = float(row["clone_a_cov"])
        clone_b_cov = float(row["clone_b_cov"])
        if identity < 98.0 or min(clone_a_cov, clone_b_cov) < 90.0:
            continue
        rows.append(
            {
                "clone_a": row["clone_a"],
                "binding_a": sanger_labels.get(row["clone_a"], row["binding_a"]),
                "len_a": row["len_a"],
                "clone_b": row["clone_b"],
                "binding_b": sanger_labels.get(row["clone_b"], row["binding_b"]),
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


def write_exports(rows: Sequence[dict]) -> Dict[str, int]:
    counts = {}
    by_label = defaultdict(list)
    for row in rows:
        by_label[row["final_label_margin05_context"]].append(row)

    for label in FASTA_LABELS:
        sub = by_label[label]
        write_fasta(FASTA_DIR / f"tim3_{label}_nt_scfv.fasta", sub, "nt_scfv", label)
        write_fasta(FASTA_DIR / f"tim3_{label}_aa_scfv.fasta", sub, "aa_scfv", label)
        counts[label] = len(sub)

    strict = by_label["bind"]
    expanded = by_label["bind"] + by_label["bind_context"]
    write_fasta(FASTA_DIR / "tim3_strict_bind_positive_nt_scfv.fasta", strict, "nt_scfv", "strict_bind_positive")
    write_fasta(FASTA_DIR / "tim3_strict_bind_positive_aa_scfv.fasta", strict, "aa_scfv", "strict_bind_positive")
    write_fasta(FASTA_DIR / "tim3_expanded_bind_positive_nt_scfv.fasta", expanded, "nt_scfv", "expanded_bind_positive")
    write_fasta(FASTA_DIR / "tim3_expanded_bind_positive_aa_scfv.fasta", expanded, "aa_scfv", "expanded_bind_positive")
    counts["strict_bind_positive"] = len(strict)
    counts["expanded_bind_positive"] = len(expanded)
    return counts


def markdown_table(rows: Sequence[dict], columns: Sequence[str], max_rows: int | None = None) -> str:
    display = rows[:max_rows] if max_rows is not None else rows
    lines = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("|" + "|".join("---" for _ in columns) + "|")
    for row in display:
        lines.append("| " + " | ".join(str(row[col]) for col in columns) + " |")
    return "\n".join(lines)


def write_report(summary: Sequence[dict], top_rows: Sequence[dict], conflicts: Sequence[dict], fasta_counts: Dict[str, int]) -> None:
    top_bind_like = [
        row for row in top_rows if row["label"] in {"bind", "bind_context", "no_bind", "ambiguous"}
    ][:20]
    direct_summary = [row for row in summary if row["label"] in DIRECT_LABELS]
    report = f"""# TIM-3 Label Report

## Label Rules

- `bind`: direct high-confidence `nt_scFv` alignment to a validated Sanger `yes` clone.
- `bind_context`: `nt_scFv` was ambiguous, but full Nanopore read context favors the Sanger `yes` reference by >0.5 percentage points in >=90% of represented reads.
- `no_bind`: direct high-confidence `nt_scFv` alignment to a validated Sanger `no` clone.
- `ambiguous`: unresolved `yes`/`no` conflict.
- `unknown`: no strong validated Sanger match.

For downstream analysis, use:

- Strict bind-positive: `bind`
- Expanded bind-positive: `bind + bind_context`
- Negative: `no_bind`
- Gray/excluded: `ambiguous + unknown`

## Read Summary

{markdown_table(direct_summary, ["label", "clones", "barcode02_reads", "barcode03_reads", "tim3_reads", "read_fraction"])}

## Combined Positive Tiers

{markdown_table([row for row in summary if row["label"] in {"strict_bind_positive", "expanded_bind_positive"}], ["label", "clones", "barcode02_reads", "barcode03_reads", "tim3_reads"])}

## Top Clones

{markdown_table(top_bind_like, ["label", "rank", "nt_scfv_id", "barcode02", "barcode03", "tim3_reads", "best_yes_clone", "best_no_clone", "confidence"], max_rows=20)}

## Sanger Conflict Audit

There are {len(conflicts)} opposite-label Sanger pairs with >=98% identity and >=90% coverage. These conflicts explain why some Nanopore clones cannot be labeled from `nt_scFv` alone.

Top conflicts:

{markdown_table(conflicts, ["clone_a", "binding_a", "clone_b", "binding_b", "identity", "clone_a_cov", "clone_b_cov", "severity"], max_rows=20)}

## FASTA Exports

- `tim3_bind_nt_scfv.fasta` / `tim3_bind_aa_scfv.fasta`: {fasta_counts.get("bind", 0)} clones.
- `tim3_bind_context_nt_scfv.fasta` / `tim3_bind_context_aa_scfv.fasta`: {fasta_counts.get("bind_context", 0)} clones.
- `tim3_no_bind_nt_scfv.fasta` / `tim3_no_bind_aa_scfv.fasta`: {fasta_counts.get("no_bind", 0)} clones.
- `tim3_strict_bind_positive_*`: {fasta_counts.get("strict_bind_positive", 0)} clones.
- `tim3_expanded_bind_positive_*`: {fasta_counts.get("expanded_bind_positive", 0)} clones.

## Files

- `label_read_summary.tsv`
- `top_clones_by_label.tsv`
- `sanger_conflict_audit.tsv`
- `label_matrix_for_plotting.tsv`
- `fasta/`
"""
    (OUT_DIR / "TIM3_LABEL_REPORT.md").write_text(textwrap.dedent(report).strip() + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_joined_rows()
    summary = make_label_summary(rows)
    top_rows = make_top_clones(rows)
    conflicts = make_conflict_audit()
    label_matrix = make_label_matrix(rows)
    fasta_counts = write_exports(rows)

    write_tsv(OUT_DIR / "label_read_summary.tsv", summary)
    write_tsv(OUT_DIR / "top_clones_by_label.tsv", top_rows)
    write_tsv(OUT_DIR / "sanger_conflict_audit.tsv", conflicts)
    write_tsv(OUT_DIR / "label_matrix_for_plotting.tsv", label_matrix)
    write_report(summary, top_rows, conflicts, fasta_counts)

    print(f"Wrote final TIM-3 exports to {OUT_DIR}")


if __name__ == "__main__":
    main()
