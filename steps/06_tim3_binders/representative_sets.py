#!/usr/bin/env python3
"""Build final TIM-3 abundance-diverse representative exports.

This script takes the outputs from tim3_abundance_diversity_sweep.py and writes
one self-contained TSV per selected mode/identity threshold. Each TSV contains
all candidates assigned to the nearest selected representative, plus cluster
summary fields and the candidate/representative sequences. It also writes one
representative-only FASTA per selected threshold.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


# Set from the command line in main()
SWEEP_DIR: Path
OUT_DIR: Path

SELECTED_EXPORTS = (
    ("protein", "95", "95%"),
    ("protein", "96", "96%"),
    ("protein", "97", "97%"),
    ("dna", "99", "99%"),
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(
            f"Missing required input: {path}\n"
            "Regenerate sweep inputs with: "
            "steps/06_tim3_binders/diversity_sweep.py"
        )
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def clean_output_dir() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path in OUT_DIR.iterdir():
        if path.is_file():
            path.unlink()


def load_summary() -> dict[tuple[str, str], dict[str, str]]:
    summary = {}
    for row in read_tsv(SWEEP_DIR / "diversity_sweep_summary.tsv"):
        threshold = row["identity_threshold"]
        if threshold.endswith(".0"):
            threshold = threshold[:-2]
        summary[(row["mode"], threshold)] = row
    return summary


def source_base(mode: str, threshold_label: str) -> Path:
    return SWEEP_DIR / f"tim3_top100_abundance_diverse_strict_binder_{mode}_identity{threshold_label}"


def write_representative_fasta(
    mode: str,
    threshold_label: str,
    threshold_display: str,
    representative_rows: list[dict[str, str]],
    clusters: dict[str, dict[str, str]],
) -> str:
    seq_field = "aa_scfv_sequence" if mode == "protein" else "nt_scfv_sequence"
    prefix = "AA" if mode == "protein" else "NT"
    out_path = OUT_DIR / f"tim3_strict_binder_{mode}_identity{threshold_label}_representatives.fasta"
    with out_path.open("w") as handle:
        for row in representative_rows:
            cluster = clusters[row["selection_rank"]]
            header = (
                f"TIM3_strict_binder_{prefix}_identity{threshold_label}_rank{int(row['selection_rank']):03d}"
                f"|mode={mode}"
                f"|identity_threshold={threshold_display}"
                f"|candidate_id={row['candidate_id']}"
                f"|abundance_rank={row['abundance_rank']}"
                f"|reads={row['total_tim3_reads']}"
                f"|barcode02={row['barcode02_tim3_round3_reads']}"
                f"|barcode03={row['barcode03_tim3_round4_reads']}"
                f"|cluster_members={cluster['cluster_member_count']}"
                f"|cluster_reads={cluster['cluster_total_reads']}"
                f"|rep_nt_id={row['representative_nt_scfv_id']}"
            )
            handle.write(f">{header}\n")
            sequence = row[seq_field]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    return out_path.name


def export_clustered_tsv(mode: str, threshold_label: str, threshold_display: str) -> dict[str, str]:
    base = source_base(mode, threshold_label)
    representative_rows = read_tsv(base.with_suffix(".tsv"))
    member_rows = read_tsv(Path(f"{base}_cluster_members.tsv"))
    cluster_rows = read_tsv(Path(f"{base}_clusters.tsv"))
    clusters = {row["cluster_rank"]: row for row in cluster_rows}

    fieldnames = [
        "mode",
        "identity_threshold",
        "cluster_rank",
        "is_representative",
        "identity_to_representative",
        "candidate_id",
        "abundance_rank",
        "total_tim3_reads",
        "barcode02_tim3_round3_reads",
        "barcode03_tim3_round4_reads",
        "aa_scfv_len",
        "nt_scfv_len",
        "nt_variant_count",
        "representative_nt_scfv_id",
        "best_yes_sanger_clones",
        "best_no_sanger_clones",
        "representative_candidate_id",
        "representative_abundance_rank",
        "representative_reads",
        "cluster_member_count",
        "cluster_total_reads",
        "cluster_barcode02_reads",
        "cluster_barcode03_reads",
        "cluster_unique_aa_sequences",
        "cluster_unique_nt_sequences",
        "min_identity_to_representative",
        "median_identity_to_representative",
        "aa_scfv_sequence",
        "nt_scfv_sequence",
        "nt_scfv_ids",
        "representative_aa_scfv_sequence",
        "representative_nt_scfv_sequence",
    ]

    output_rows = []
    for member in member_rows:
        cluster = clusters[member["cluster_rank"]]
        output_rows.append(
            {
                "mode": mode,
                "identity_threshold": threshold_display,
                "cluster_rank": member["cluster_rank"],
                "is_representative": member["is_representative"],
                "identity_to_representative": member["identity_to_representative"],
                "candidate_id": member["candidate_id"],
                "abundance_rank": member["abundance_rank"],
                "total_tim3_reads": member["total_tim3_reads"],
                "barcode02_tim3_round3_reads": member["barcode02_tim3_round3_reads"],
                "barcode03_tim3_round4_reads": member["barcode03_tim3_round4_reads"],
                "aa_scfv_len": member["aa_scfv_len"],
                "nt_scfv_len": member["nt_scfv_len"],
                "nt_variant_count": member["nt_variant_count"],
                "representative_nt_scfv_id": member["representative_nt_scfv_id"],
                "best_yes_sanger_clones": member["best_yes_sanger_clones"],
                "best_no_sanger_clones": member["best_no_sanger_clones"],
                "representative_candidate_id": cluster["representative_candidate_id"],
                "representative_abundance_rank": cluster["representative_abundance_rank"],
                "representative_reads": cluster["representative_reads"],
                "cluster_member_count": cluster["cluster_member_count"],
                "cluster_total_reads": cluster["cluster_total_reads"],
                "cluster_barcode02_reads": cluster["cluster_barcode02_reads"],
                "cluster_barcode03_reads": cluster["cluster_barcode03_reads"],
                "cluster_unique_aa_sequences": cluster["cluster_unique_aa_sequences"],
                "cluster_unique_nt_sequences": cluster["cluster_unique_nt_sequences"],
                "min_identity_to_representative": cluster["min_identity_to_representative"],
                "median_identity_to_representative": cluster["median_identity_to_representative"],
                "aa_scfv_sequence": member["aa_scfv_sequence"],
                "nt_scfv_sequence": member["nt_scfv_sequence"],
                "nt_scfv_ids": member["nt_scfv_ids"],
                "representative_aa_scfv_sequence": cluster["representative_aa_scfv_sequence"],
                "representative_nt_scfv_sequence": cluster["representative_nt_scfv_sequence"],
            }
        )

    out_path = OUT_DIR / f"tim3_strict_binder_{mode}_identity{threshold_label}_clustered_candidates.tsv"
    write_tsv(out_path, output_rows, fieldnames)
    fasta_name = write_representative_fasta(
        mode,
        threshold_label,
        threshold_display,
        representative_rows,
        clusters,
    )
    representatives = sum(1 for row in output_rows if row["is_representative"] == "1")
    return {
        "mode": mode,
        "threshold_label": threshold_label,
        "threshold_display": threshold_display,
        "file": out_path.name,
        "fasta": fasta_name,
        "rows": str(len(output_rows)),
        "representatives": str(representatives),
    }


def write_readme(exports: list[dict[str, str]], summary: dict[tuple[str, str], dict[str, str]]) -> None:
    lines = [
        "# TIM-3 Strict Binder Abundance-Diverse Clustered Sets",
        "",
        "These are the final selected TIM-3 abundance-diverse deliverables.",
        "",
        "The starting point was the conservative strict-binder label (`final_label == bind`).",
        "Those labels were assigned by globally aligning extracted Nanopore scFv nucleotide",
        "sequences against the Sanger-validated TIM-3 bind/no-bind sequences.",
        "Context-resolved probable binders (`bind_context`) are not included here.",
        "",
        "The representative-selection method is abundance-first nonredundant selection.",
        "For each selected mode/threshold, candidates were ranked by total TIM-3 read abundance.",
        "Walking down that ranked list, a candidate was kept as a representative only if its",
        "global alignment identity to every previously selected higher-abundance representative",
        "was below the threshold. After representative selection, every candidate in that mode",
        "was assigned to its nearest selected representative. This gives cluster membership and",
        "cluster-level read summaries for auditing.",
        "",
        "Protein mode collapses identical amino-acid scFv sequences and sums their nucleotide",
        "variant read counts. DNA mode uses individual nucleotide scFv sequences.",
        "The TSV files contain all candidates and cluster assignments. The FASTA files contain",
        "representatives only.",
        "",
        "## Files",
        "",
        "| TSV | FASTA | mode | identity threshold | candidate rows | representatives | selected reads | selected read fraction |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for export in exports:
        row = summary[(export["mode"], export["threshold_label"])]
        lines.append(
            "| {file} | {fasta} | {mode} | {threshold} | {rows} | {representatives} | {reads} | {fraction} |".format(
                file=f"`{export['file']}`",
                fasta=f"`{export['fasta']}`",
                mode=export["mode"],
                threshold=export["threshold_display"],
                rows=export["rows"],
                representatives=export["representatives"],
                reads=row["selected_reads"],
                fraction=row["selected_read_fraction"],
            )
        )
    lines.extend(
        [
            "",
            "## Main Columns",
            "",
            "- `is_representative`: `1` means the row is one of the selected abundance-diverse representatives.",
            "- `cluster_rank`: selected representative rank; all rows with the same value belong to that representative cluster.",
            "- `identity_to_representative`: global identity between the candidate and its assigned representative.",
            "- `cluster_member_count` and `cluster_total_reads`: size and total read support of the assigned cluster.",
            "- `aa_scfv_sequence` and `nt_scfv_sequence`: candidate sequences.",
            "- `representative_aa_scfv_sequence` and `representative_nt_scfv_sequence`: assigned representative sequences.",
            "",
            "This is not connected-component clustering. The exact connected-component sensitivity check",
            "collapsed most candidates into a few large chained clusters, so it was not used for these final",
            "deliverables.",
            "",
        ]
    )
    (OUT_DIR / "README.md").write_text("\n".join(lines))


def main() -> None:
    global SWEEP_DIR, OUT_DIR
    parser = argparse.ArgumentParser(description="Export the selected TIM-3 representative sets.")
    parser.add_argument("--sweep-dir", type=Path, required=True, help="output dir of diversity_sweep.py")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    SWEEP_DIR, OUT_DIR = args.sweep_dir, args.out_dir

    clean_output_dir()
    summary = load_summary()
    exports = [
        export_clustered_tsv(mode, threshold_label, threshold_display)
        for mode, threshold_label, threshold_display in SELECTED_EXPORTS
    ]
    write_readme(exports, summary)
    for export in exports:
        print(
            f"{export['file']}: {export['rows']} rows, "
            f"{export['representatives']} representatives",
            flush=True,
        )
    print(f"Wrote general exports to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
