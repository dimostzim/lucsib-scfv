#!/usr/bin/env python3
"""Abundance-first diversity sweeps for TIM-3 strict binder scFvs.

Two modes are produced:
  - protein: candidates are unique amino-acid scFv sequences.
  - dna: candidates are unique nucleotide scFv sequences.

For each identity threshold, candidates are sorted by abundance and greedily
selected only if they are below the threshold identity to all previously
selected higher-abundance candidates.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import edlib
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency: edlib. Use the Python environment from envs/python.yml."
    ) from exc


# Set from the command line in main()
INPUT_TSV: Path
OUT_DIR: Path
THRESHOLDS = (99.0, 98.0, 97.0, 96.0, 95.0)
TOP_N = 100


@dataclass
class Candidate:
    candidate_id: str
    mode: str
    sequence: str
    total_reads: int
    barcode02_reads: int
    barcode03_reads: int
    aa_scfv_sequence: str
    nt_scfv_sequence: str
    aa_scfv_len: int
    nt_scfv_len: int
    nt_variant_count: int
    nt_scfv_ids: List[str] = field(default_factory=list)
    representative_nt_scfv_id: str = ""
    representative_nt_reads: int = 0
    best_yes_sanger_clones: set = field(default_factory=set)
    best_no_sanger_clones: set = field(default_factory=set)
    abundance_rank: int = 0


def parse_int(value: str) -> int:
    if value is None or value == "":
        return 0
    return int(float(value))


def parse_cigar_identity(cigar: str) -> float:
    matches = mismatches = insertions = deletions = 0
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        count = int(number)
        number = ""
        if char == "=":
            matches += count
        elif char == "X":
            mismatches += count
        elif char == "I":
            insertions += count
        elif char == "D":
            deletions += count
    aligned_len = matches + mismatches + insertions + deletions
    if aligned_len == 0:
        return 0.0
    return 100.0 * matches / aligned_len


class IdentityCache:
    def __init__(self, candidates: Sequence[Candidate]):
        self.candidates = candidates
        self.cache: Dict[Tuple[int, int], float] = {}

    def identity(self, i: int, j: int) -> float:
        key = (i, j) if i < j else (j, i)
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        result = edlib.align(
            self.candidates[i].sequence,
            self.candidates[j].sequence,
            mode="NW",
            task="path",
        )
        if result["editDistance"] < 0 or not result.get("cigar"):
            ident = 0.0
        else:
            ident = parse_cigar_identity(result["cigar"])
        self.cache[key] = ident
        return ident


def load_strict_binder_rows() -> List[dict]:
    rows = []
    with INPUT_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["final_label"] == "bind":
                rows.append(row)
    if not rows:
        raise RuntimeError(f"No strict bind rows found in {INPUT_TSV}")
    return rows


def build_protein_candidates(rows: Sequence[dict]) -> List[Candidate]:
    grouped: Dict[str, Candidate] = {}
    for row in rows:
        aa = row["aa_scfv_sequence"]
        reads = parse_int(row["total_tim3_reads"])
        b02 = parse_int(row["barcode02_tim3_round3_reads"])
        b03 = parse_int(row["barcode03_tim3_round4_reads"])
        candidate = grouped.get(aa)
        if candidate is None:
            candidate = Candidate(
                candidate_id=f"aa_{len(grouped) + 1:06d}",
                mode="protein",
                sequence=aa,
                total_reads=0,
                barcode02_reads=0,
                barcode03_reads=0,
                aa_scfv_sequence=aa,
                nt_scfv_sequence=row["nt_scfv_sequence"],
                aa_scfv_len=len(aa),
                nt_scfv_len=len(row["nt_scfv_sequence"]),
                nt_variant_count=0,
            )
            grouped[aa] = candidate
        candidate.total_reads += reads
        candidate.barcode02_reads += b02
        candidate.barcode03_reads += b03
        candidate.nt_variant_count += 1
        candidate.nt_scfv_ids.append(row["nt_scfv_id"])
        if row["best_yes_sanger_clone"]:
            candidate.best_yes_sanger_clones.add(row["best_yes_sanger_clone"])
        if row["best_no_sanger_clone"]:
            candidate.best_no_sanger_clones.add(row["best_no_sanger_clone"])
        if reads > candidate.representative_nt_reads:
            candidate.representative_nt_reads = reads
            candidate.representative_nt_scfv_id = row["nt_scfv_id"]
            candidate.nt_scfv_sequence = row["nt_scfv_sequence"]
            candidate.nt_scfv_len = len(row["nt_scfv_sequence"])

    candidates = sorted(grouped.values(), key=lambda c: (-c.total_reads, c.aa_scfv_sequence))
    for rank, candidate in enumerate(candidates, start=1):
        candidate.abundance_rank = rank
    return candidates


def build_dna_candidates(rows: Sequence[dict]) -> List[Candidate]:
    candidates = []
    for row in rows:
        nt = row["nt_scfv_sequence"]
        aa = row["aa_scfv_sequence"]
        reads = parse_int(row["total_tim3_reads"])
        candidate = Candidate(
            candidate_id=row["nt_scfv_id"],
            mode="dna",
            sequence=nt,
            total_reads=reads,
            barcode02_reads=parse_int(row["barcode02_tim3_round3_reads"]),
            barcode03_reads=parse_int(row["barcode03_tim3_round4_reads"]),
            aa_scfv_sequence=aa,
            nt_scfv_sequence=nt,
            aa_scfv_len=len(aa),
            nt_scfv_len=len(nt),
            nt_variant_count=1,
            nt_scfv_ids=[row["nt_scfv_id"]],
            representative_nt_scfv_id=row["nt_scfv_id"],
            representative_nt_reads=reads,
            best_yes_sanger_clones={row["best_yes_sanger_clone"]} if row["best_yes_sanger_clone"] else set(),
            best_no_sanger_clones={row["best_no_sanger_clone"]} if row["best_no_sanger_clone"] else set(),
        )
        candidates.append(candidate)
    candidates.sort(key=lambda c: (-c.total_reads, c.nt_scfv_sequence))
    for rank, candidate in enumerate(candidates, start=1):
        candidate.abundance_rank = rank
    return candidates


def candidate_metadata(candidate: Candidate) -> dict:
    return {
        "candidate_id": candidate.candidate_id,
        "abundance_rank": candidate.abundance_rank,
        "total_tim3_reads": candidate.total_reads,
        "barcode02_tim3_round3_reads": candidate.barcode02_reads,
        "barcode03_tim3_round4_reads": candidate.barcode03_reads,
        "aa_scfv_len": candidate.aa_scfv_len,
        "nt_scfv_len": candidate.nt_scfv_len,
        "nt_variant_count": candidate.nt_variant_count,
        "representative_nt_scfv_id": candidate.representative_nt_scfv_id,
        "best_yes_sanger_clones": ";".join(sorted(candidate.best_yes_sanger_clones, key=sort_clone)),
        "best_no_sanger_clones": ";".join(sorted(candidate.best_no_sanger_clones, key=sort_clone)),
        "aa_scfv_sequence": candidate.aa_scfv_sequence,
        "nt_scfv_sequence": candidate.nt_scfv_sequence,
        "nt_scfv_ids": ";".join(candidate.nt_scfv_ids),
    }


def select_abundance_diverse(
    candidates: Sequence[Candidate],
    threshold: float,
    top_n: int = TOP_N,
) -> Tuple[List[dict], List[int], IdentityCache, dict]:
    cache = IdentityCache(candidates)
    selected_indices: List[int] = []
    selected_rows: List[dict] = []
    skipped = 0
    considered = 0

    for idx, candidate in enumerate(candidates):
        considered += 1
        nearest_identity = 0.0
        nearest_rank = ""
        nearest_candidate_id = ""
        redundant = False
        for selected_order, selected_idx in enumerate(selected_indices, start=1):
            ident = cache.identity(idx, selected_idx)
            if ident > nearest_identity:
                nearest_identity = ident
                nearest_rank = str(selected_order)
                nearest_candidate_id = candidates[selected_idx].candidate_id
            if ident >= threshold:
                redundant = True
        if redundant:
            skipped += 1
            continue

        selected_indices.append(idx)
        row = {
            "selection_rank": len(selected_indices),
            **candidate_metadata(candidate),
            "nearest_previously_selected_rank": nearest_rank,
            "nearest_previously_selected_candidate_id": nearest_candidate_id,
            "nearest_previously_selected_identity": round(nearest_identity, 4),
        }
        selected_rows.append(row)
        if len(selected_indices) >= top_n:
            break

    pairwise = []
    for i, selected_i in enumerate(selected_indices):
        for selected_j in selected_indices[i + 1 :]:
            pairwise.append(cache.identity(selected_i, selected_j))

    selected_candidates = [candidates[i] for i in selected_indices]
    total_strict_reads = sum(c.total_reads for c in candidates)
    selected_reads = sum(c.total_reads for c in selected_candidates)
    selected_aa = {c.aa_scfv_sequence for c in selected_candidates}
    selected_nt = {c.nt_scfv_sequence for c in selected_candidates}
    summary = {
        "mode": candidates[0].mode if candidates else "",
        "identity_threshold": threshold,
        "candidate_count": len(candidates),
        "selected_count": len(selected_candidates),
        "total_strict_reads": total_strict_reads,
        "selected_reads": selected_reads,
        "selected_read_fraction": round(selected_reads / total_strict_reads, 6) if total_strict_reads else 0,
        "lowest_abundance_rank_considered": considered,
        "lowest_selected_abundance_rank": max((c.abundance_rank for c in selected_candidates), default=0),
        "skipped_as_too_similar": skipped,
        "min_selected_reads": min((c.total_reads for c in selected_candidates), default=0),
        "median_selected_reads": round(statistics.median([c.total_reads for c in selected_candidates]), 4)
        if selected_candidates
        else 0,
        "unique_aa_sequences_selected": len(selected_aa),
        "unique_nt_sequences_selected": len(selected_nt),
        "min_pairwise_identity": round(min(pairwise), 4) if pairwise else "",
        "median_pairwise_identity": round(statistics.median(pairwise), 4) if pairwise else "",
        "mean_pairwise_identity": round(statistics.mean(pairwise), 4) if pairwise else "",
        "max_pairwise_identity": round(max(pairwise), 4) if pairwise else "",
    }
    return selected_rows, selected_indices, cache, summary


def assign_clusters(
    candidates: Sequence[Candidate],
    selected_indices: Sequence[int],
    selected_rows: Sequence[dict],
    cache: IdentityCache,
) -> Tuple[List[dict], List[dict], dict]:
    selected_index_to_rank = {
        selected_idx: int(row["selection_rank"])
        for selected_idx, row in zip(selected_indices, selected_rows)
    }
    rank_to_selected_idx = {
        int(row["selection_rank"]): selected_idx
        for selected_idx, row in zip(selected_indices, selected_rows)
    }

    cluster_stats = {
        rank: {
            "cluster_member_count": 0,
            "cluster_total_reads": 0,
            "cluster_barcode02_reads": 0,
            "cluster_barcode03_reads": 0,
            "cluster_unique_aa_sequences": set(),
            "cluster_unique_nt_sequences": set(),
            "member_identities_to_rep": [],
        }
        for rank in selected_index_to_rank.values()
    }

    member_rows = []
    for idx, candidate in enumerate(candidates):
        if idx in selected_index_to_rank:
            assigned_rank = selected_index_to_rank[idx]
            identity_to_rep = 100.0
        else:
            best = None
            for rank, selected_idx in rank_to_selected_idx.items():
                ident = cache.identity(idx, selected_idx)
                key = (ident, -rank)
                if best is None or key > best[0]:
                    best = (key, rank, ident)
            if best is None:
                continue
            assigned_rank = best[1]
            identity_to_rep = best[2]

        stats = cluster_stats[assigned_rank]
        stats["cluster_member_count"] += 1
        stats["cluster_total_reads"] += candidate.total_reads
        stats["cluster_barcode02_reads"] += candidate.barcode02_reads
        stats["cluster_barcode03_reads"] += candidate.barcode03_reads
        stats["cluster_unique_aa_sequences"].add(candidate.aa_scfv_sequence)
        stats["cluster_unique_nt_sequences"].add(candidate.nt_scfv_sequence)
        stats["member_identities_to_rep"].append(identity_to_rep)

        member_rows.append(
            {
                "cluster_rank": assigned_rank,
                "is_representative": int(idx == rank_to_selected_idx[assigned_rank]),
                "identity_to_representative": round(identity_to_rep, 4),
                **candidate_metadata(candidate),
            }
        )

    cluster_rows = []
    for selected_row in selected_rows:
        rank = int(selected_row["selection_rank"])
        stats = cluster_stats[rank]
        identities = stats["member_identities_to_rep"]
        cluster_rows.append(
            {
                "cluster_rank": rank,
                "representative_candidate_id": selected_row["candidate_id"],
                "representative_abundance_rank": selected_row["abundance_rank"],
                "representative_reads": selected_row["total_tim3_reads"],
                "cluster_member_count": stats["cluster_member_count"],
                "cluster_total_reads": stats["cluster_total_reads"],
                "cluster_barcode02_reads": stats["cluster_barcode02_reads"],
                "cluster_barcode03_reads": stats["cluster_barcode03_reads"],
                "cluster_unique_aa_sequences": len(stats["cluster_unique_aa_sequences"]),
                "cluster_unique_nt_sequences": len(stats["cluster_unique_nt_sequences"]),
                "min_identity_to_representative": round(min(identities), 4) if identities else "",
                "median_identity_to_representative": round(statistics.median(identities), 4)
                if identities
                else "",
                "representative_aa_scfv_sequence": selected_row["aa_scfv_sequence"],
                "representative_nt_scfv_sequence": selected_row["nt_scfv_sequence"],
            }
        )

    total_reads = sum(row["cluster_total_reads"] for row in cluster_rows)
    cluster_member_counts = [row["cluster_member_count"] for row in cluster_rows]
    cluster_read_counts = [row["cluster_total_reads"] for row in cluster_rows]
    cluster_summary = {
        "assigned_candidate_count": len(member_rows),
        "assigned_total_reads": total_reads,
        "min_cluster_size": min(cluster_member_counts) if cluster_member_counts else 0,
        "median_cluster_size": round(statistics.median(cluster_member_counts), 4)
        if cluster_member_counts
        else 0,
        "max_cluster_size": max(cluster_member_counts) if cluster_member_counts else 0,
        "min_cluster_total_reads": min(cluster_read_counts) if cluster_read_counts else 0,
        "median_cluster_total_reads": round(statistics.median(cluster_read_counts), 4)
        if cluster_read_counts
        else 0,
        "max_cluster_total_reads": max(cluster_read_counts) if cluster_read_counts else 0,
    }
    return cluster_rows, member_rows, cluster_summary


def sort_clone(value: str):
    return (0, int(value)) if value.isdigit() else (1, value)


def threshold_label(threshold: float) -> str:
    if float(threshold).is_integer():
        return str(int(threshold))
    return str(threshold).replace(".", "p")


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Optional[Sequence[str]] = None) -> None:
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


def write_fasta(path: Path, rows: Sequence[dict], mode: str, threshold: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seq_field = "aa_scfv_sequence" if mode == "protein" else "nt_scfv_sequence"
    prefix = "AA" if mode == "protein" else "NT"
    with path.open("w") as handle:
        for row in rows:
            header = (
                f"TIM3_strict_binder_{prefix}_diverse_rank{int(row['selection_rank']):03d}"
                f"|identity_threshold={threshold:g}"
                f"|abundance_rank={row['abundance_rank']}"
                f"|reads={row['total_tim3_reads']}"
                f"|barcode02={row['barcode02_tim3_round3_reads']}"
                f"|barcode03={row['barcode03_tim3_round4_reads']}"
                f"|nt_variants={row['nt_variant_count']}"
                f"|rep_nt_id={row['representative_nt_scfv_id']}"
                f"|nearest_identity={row['nearest_previously_selected_identity']}"
            )
            handle.write(f">{header}\n")
            seq = row[seq_field]
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_readme(summary_rows: Sequence[dict]) -> None:
    lines = [
        "# TIM-3 Abundance-First Diversity Sweep",
        "",
        "This folder contains greedy abundance-first diversity-filtered top-100 strict binder sets.",
        "",
        "Protein mode uses unique amino-acid scFv sequences and filters by global amino-acid identity.",
        "DNA mode uses unique nucleotide scFv sequences and filters by global nucleotide identity.",
        "",
        "Selection rule: candidates are sorted by total TIM-3 abundance. A candidate is selected only if",
        "its identity to every previously selected higher-abundance sequence is below the threshold.",
        "",
        "Summary table: `diversity_sweep_summary.tsv`.",
        "",
        "Generated thresholds: " + ", ".join(f"{t:g}%" for t in THRESHOLDS),
        "",
        "Each threshold has a representative TSV/FASTA, a cluster summary TSV, and a cluster member TSV.",
        "",
        "The representative set is abundance-first/diversity-filtered. Cluster assignment then assigns every",
        "candidate in that mode to its nearest selected representative, so cluster sizes/read totals can be audited.",
        "",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines))


def main() -> None:
    global INPUT_TSV, OUT_DIR
    parser = argparse.ArgumentParser(description="Abundance-first diversity sweep of TIM-3 strict binders.")
    parser.add_argument("--input-tsv", type=Path, required=True,
                        help="tim3_full_labeled_scfv_results_with_sequences.tsv")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    INPUT_TSV, OUT_DIR = args.input_tsv, args.out_dir

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_strict_binder_rows()
    modes = {
        "protein": build_protein_candidates(rows),
        "dna": build_dna_candidates(rows),
    }

    summary_rows = []
    for mode, candidates in modes.items():
        print(f"Mode={mode} candidates={len(candidates)}", flush=True)
        for threshold in THRESHOLDS:
            print(f"  threshold={threshold:g}", flush=True)
            selected_rows, selected_indices, cache, summary = select_abundance_diverse(candidates, threshold)
            cluster_rows, member_rows, cluster_summary = assign_clusters(
                candidates,
                selected_indices,
                selected_rows,
                cache,
            )
            summary.update(cluster_summary)
            summary_rows.append(summary)
            label = threshold_label(threshold)
            base = OUT_DIR / f"tim3_top100_abundance_diverse_strict_binder_{mode}_identity{label}"
            write_tsv(base.with_suffix(".tsv"), selected_rows)
            write_fasta(base.with_suffix(".fasta"), selected_rows, mode, threshold)
            write_tsv(OUT_DIR / f"{base.name}_clusters.tsv", cluster_rows)
            write_tsv(OUT_DIR / f"{base.name}_cluster_members.tsv", member_rows)

    write_tsv(OUT_DIR / "diversity_sweep_summary.tsv", summary_rows)
    write_readme(summary_rows)
    print(f"Wrote diversity sweep outputs to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
