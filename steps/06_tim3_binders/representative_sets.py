#!/usr/bin/env python3
"""Non-redundant sets of TIM-3 strict binders (label `bind`).

Candidates are unique amino-acid scFv sequences (protein mode) or unique nucleotide scFv sequences
(DNA mode), ranked by TIM-3 reads. Walking down the ranking, a candidate becomes a representative only
if its global identity to every representative selected so far is below the threshold, up to TOP_N
representatives. Every candidate is then assigned to its nearest representative.
"""

from __future__ import annotations

import argparse
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from common import align, read_tsv, write_tsv


# (mode, identity threshold in percent)
SETS = (("protein", 95.0), ("protein", 96.0), ("protein", 97.0), ("dna", 99.0))
TOP_N = 100


@dataclass
class Candidate:
    candidate_id: str
    sequence: str
    total_reads: int
    barcode02_reads: int
    barcode03_reads: int
    aa_scfv_sequence: str
    nt_scfv_sequence: str
    nt_variant_count: int
    nt_scfv_ids: List[str] = field(default_factory=list)
    representative_nt_scfv_id: str = ""
    representative_nt_reads: int = 0
    best_yes_sanger_clones: set = field(default_factory=set)
    best_no_sanger_clones: set = field(default_factory=set)
    abundance_rank: int = 0


class IdentityCache:
    """Global alignment identity of a candidate (query) to a representative (target)."""

    def __init__(self, candidates: Sequence[Candidate]):
        self.candidates = candidates
        self.cache: Dict[Tuple[int, int], float] = {}

    def identity(self, candidate: int, representative: int) -> float:
        # not symmetric: with query and target swapped, edlib can return another optimal path
        key = (candidate, representative)
        if key not in self.cache:
            self.cache[key] = align(
                self.candidates[candidate].sequence, self.candidates[representative].sequence, mode="NW"
            ).identity
        return self.cache[key]


def build_protein_candidates(rows: Sequence[dict]) -> List[Candidate]:
    grouped: Dict[str, Candidate] = {}
    for row in rows:
        aa = row["aa_scfv_sequence"]
        reads = int(row["total_tim3_reads"])
        candidate = grouped.get(aa)
        if candidate is None:
            candidate = Candidate(
                candidate_id=f"aa_{len(grouped) + 1:06d}",
                sequence=aa,
                total_reads=0,
                barcode02_reads=0,
                barcode03_reads=0,
                aa_scfv_sequence=aa,
                nt_scfv_sequence=row["nt_scfv_sequence"],
                nt_variant_count=0,
            )
            grouped[aa] = candidate
        candidate.total_reads += reads
        candidate.barcode02_reads += int(row["barcode02_tim3_round3_reads"])
        candidate.barcode03_reads += int(row["barcode03_tim3_round4_reads"])
        candidate.nt_variant_count += 1
        candidate.nt_scfv_ids.append(row["nt_scfv_id"])
        if row["best_yes_sanger_clone"]:
            candidate.best_yes_sanger_clones.add(row["best_yes_sanger_clone"])
        if row["best_no_sanger_clone"]:
            candidate.best_no_sanger_clones.add(row["best_no_sanger_clone"])
        # the nucleotide variant with the most reads represents the protein
        if reads > candidate.representative_nt_reads:
            candidate.representative_nt_reads = reads
            candidate.representative_nt_scfv_id = row["nt_scfv_id"]
            candidate.nt_scfv_sequence = row["nt_scfv_sequence"]

    candidates = sorted(grouped.values(), key=lambda c: (-c.total_reads, c.aa_scfv_sequence))
    for rank, candidate in enumerate(candidates, start=1):
        candidate.abundance_rank = rank
    return candidates


def build_dna_candidates(rows: Sequence[dict]) -> List[Candidate]:
    candidates = [
        Candidate(
            candidate_id=row["nt_scfv_id"],
            sequence=row["nt_scfv_sequence"],
            total_reads=int(row["total_tim3_reads"]),
            barcode02_reads=int(row["barcode02_tim3_round3_reads"]),
            barcode03_reads=int(row["barcode03_tim3_round4_reads"]),
            aa_scfv_sequence=row["aa_scfv_sequence"],
            nt_scfv_sequence=row["nt_scfv_sequence"],
            nt_variant_count=1,
            nt_scfv_ids=[row["nt_scfv_id"]],
            representative_nt_scfv_id=row["nt_scfv_id"],
            representative_nt_reads=int(row["total_tim3_reads"]),
            best_yes_sanger_clones={row["best_yes_sanger_clone"]} if row["best_yes_sanger_clone"] else set(),
            best_no_sanger_clones={row["best_no_sanger_clone"]} if row["best_no_sanger_clone"] else set(),
        )
        for row in rows
    ]
    candidates.sort(key=lambda c: (-c.total_reads, c.nt_scfv_sequence))
    for rank, candidate in enumerate(candidates, start=1):
        candidate.abundance_rank = rank
    return candidates


def select_representatives(
    candidates: Sequence[Candidate], threshold: float, cache: IdentityCache
) -> Tuple[List[int], int]:
    """Indices of the representatives, and how many candidates were considered."""
    selected: List[int] = []
    for idx in range(len(candidates)):
        if all(cache.identity(idx, rep) < threshold for rep in selected):
            selected.append(idx)
            if len(selected) == TOP_N:
                return selected, idx + 1
    return selected, len(candidates)


def assign_clusters(
    candidates: Sequence[Candidate], selected: Sequence[int], cache: IdentityCache
) -> List[Tuple[int, float]]:
    """(cluster rank, identity to its representative) of every candidate; ties go to the higher-ranked representative."""
    rank_of = {idx: rank for rank, idx in enumerate(selected, start=1)}
    assignments = []
    for idx in range(len(candidates)):
        if idx in rank_of:
            assignments.append((rank_of[idx], 100.0))
        else:
            identity, neg_rank = max(
                (cache.identity(idx, rep), -rank) for rank, rep in enumerate(selected, start=1)
            )
            assignments.append((-neg_rank, identity))
    return assignments


def describe_clusters(
    candidates: Sequence[Candidate], selected: Sequence[int], assignments: Sequence[Tuple[int, float]]
) -> List[dict]:
    """Cluster columns of the output table, one dict per representative in rank order."""
    members: List[List[Tuple[Candidate, float]]] = [[] for _ in selected]
    for candidate, (rank, identity) in zip(candidates, assignments):
        members[rank - 1].append((candidate, identity))
    clusters = []
    for rep_idx, cluster in zip(selected, members):
        rep = candidates[rep_idx]
        identities = [identity for _, identity in cluster]
        clusters.append(
            {
                "representative_candidate_id": rep.candidate_id,
                "representative_abundance_rank": rep.abundance_rank,
                "representative_reads": rep.total_reads,
                "cluster_member_count": len(cluster),
                "cluster_total_reads": sum(c.total_reads for c, _ in cluster),
                "cluster_barcode02_reads": sum(c.barcode02_reads for c, _ in cluster),
                "cluster_barcode03_reads": sum(c.barcode03_reads for c, _ in cluster),
                "cluster_unique_aa_sequences": len({c.aa_scfv_sequence for c, _ in cluster}),
                "cluster_unique_nt_sequences": len({c.nt_scfv_sequence for c, _ in cluster}),
                "min_identity_to_representative": round(min(identities), 4),
                "median_identity_to_representative": round(statistics.median(identities), 4),
            }
        )
    return clusters


def make_clustered_rows(
    mode: str,
    threshold: float,
    candidates: Sequence[Candidate],
    selected: Sequence[int],
    assignments: Sequence[Tuple[int, float]],
    clusters: Sequence[dict],
) -> List[dict]:
    rows = []
    for idx, (candidate, (rank, identity)) in enumerate(zip(candidates, assignments)):
        rep = candidates[selected[rank - 1]]
        rows.append(
            {
                "mode": mode,
                "identity_threshold": f"{threshold:g}%",
                "cluster_rank": rank,
                "is_representative": int(idx == selected[rank - 1]),
                "identity_to_representative": round(identity, 4),
                "candidate_id": candidate.candidate_id,
                "abundance_rank": candidate.abundance_rank,
                "total_tim3_reads": candidate.total_reads,
                "barcode02_tim3_round3_reads": candidate.barcode02_reads,
                "barcode03_tim3_round4_reads": candidate.barcode03_reads,
                "aa_scfv_len": len(candidate.aa_scfv_sequence),
                "nt_scfv_len": len(candidate.nt_scfv_sequence),
                "nt_variant_count": candidate.nt_variant_count,
                "representative_nt_scfv_id": candidate.representative_nt_scfv_id,
                "best_yes_sanger_clones": ";".join(sorted(candidate.best_yes_sanger_clones, key=int)),
                "best_no_sanger_clones": ";".join(sorted(candidate.best_no_sanger_clones, key=int)),
                **clusters[rank - 1],
                "aa_scfv_sequence": candidate.aa_scfv_sequence,
                "nt_scfv_sequence": candidate.nt_scfv_sequence,
                "nt_scfv_ids": ";".join(candidate.nt_scfv_ids),
                "representative_aa_scfv_sequence": rep.aa_scfv_sequence,
                "representative_nt_scfv_sequence": rep.nt_scfv_sequence,
            }
        )
    return rows


def write_representatives_fasta(
    path: Path,
    mode: str,
    threshold: float,
    candidates: Sequence[Candidate],
    selected: Sequence[int],
    clusters: Sequence[dict],
) -> None:
    prefix = "AA" if mode == "protein" else "NT"
    with path.open("w") as handle:
        for rank, (idx, cluster) in enumerate(zip(selected, clusters), start=1):
            rep = candidates[idx]
            header = (
                f"TIM3_strict_binder_{prefix}_identity{threshold:g}_rank{rank:03d}"
                f"|mode={mode}"
                f"|identity_threshold={threshold:g}%"
                f"|candidate_id={rep.candidate_id}"
                f"|abundance_rank={rep.abundance_rank}"
                f"|reads={rep.total_reads}"
                f"|barcode02={rep.barcode02_reads}"
                f"|barcode03={rep.barcode03_reads}"
                f"|cluster_members={cluster['cluster_member_count']}"
                f"|cluster_reads={cluster['cluster_total_reads']}"
                f"|rep_nt_id={rep.representative_nt_scfv_id}"
            )
            handle.write(f">{header}\n")
            for start in range(0, len(rep.sequence), 80):
                handle.write(rep.sequence[start : start + 80] + "\n")


def summarize_set(
    mode: str,
    threshold: float,
    candidates: Sequence[Candidate],
    selected: Sequence[int],
    considered: int,
    clusters: Sequence[dict],
    cache: IdentityCache,
) -> dict:
    reps = [candidates[idx] for idx in selected]
    # each representative against the representatives selected before it, as during selection
    pairwise = [cache.identity(later, earlier) for n, earlier in enumerate(selected) for later in selected[n + 1 :]]
    total_reads = sum(c.total_reads for c in candidates)
    selected_reads = sum(c.total_reads for c in reps)
    sizes = [cluster["cluster_member_count"] for cluster in clusters]
    cluster_reads = [cluster["cluster_total_reads"] for cluster in clusters]
    return {
        "mode": mode,
        "identity_threshold": threshold,
        "candidate_count": len(candidates),
        "selected_count": len(reps),
        "total_strict_reads": total_reads,
        "selected_reads": selected_reads,
        "selected_read_fraction": round(selected_reads / total_reads, 6),
        "lowest_abundance_rank_considered": considered,
        "lowest_selected_abundance_rank": max(c.abundance_rank for c in reps),
        "skipped_as_too_similar": considered - len(reps),
        "min_selected_reads": min(c.total_reads for c in reps),
        "median_selected_reads": round(statistics.median([c.total_reads for c in reps]), 4),
        "unique_aa_sequences_selected": len({c.aa_scfv_sequence for c in reps}),
        "unique_nt_sequences_selected": len({c.nt_scfv_sequence for c in reps}),
        "min_pairwise_identity": round(min(pairwise), 4),
        "median_pairwise_identity": round(statistics.median(pairwise), 4),
        "mean_pairwise_identity": round(statistics.mean(pairwise), 4),
        "max_pairwise_identity": round(max(pairwise), 4),
        "assigned_candidate_count": sum(sizes),
        "assigned_total_reads": sum(cluster_reads),
        "min_cluster_size": min(sizes),
        "median_cluster_size": round(statistics.median(sizes), 4),
        "max_cluster_size": max(sizes),
        "min_cluster_total_reads": min(cluster_reads),
        "median_cluster_total_reads": round(statistics.median(cluster_reads), 4),
        "max_cluster_total_reads": max(cluster_reads),
    }


def write_readme(path: Path, files: Sequence[Tuple[str, str, int]], summaries: Sequence[dict]) -> None:
    lines = [
        "# TIM-3 strict-binder representative sets",
        "",
        "Each set is a non-redundant selection of TIM-3 scFvs labelled `bind`. Only `bind` counts as binder",
        "evidence. `bind_context` is not binder evidence and is not used here.",
        "",
        "## Labels",
        "",
        "Each TIM-3 scFv DNA sequence (`nt_scFv`) is aligned with edlib infix alignment (mode HW: the shorter",
        "of the two sequences is aligned inside the longer one, with free end gaps in the longer one) to the",
        "5 `yes` and the 5 `no` Sanger clones sharing the most 15-letter k-mers with it (at least 3), allowing",
        "at most 100 edits. Identity is the number of identical positions divided by the alignment length, gaps",
        "included. A Sanger match counts if it has at least 99% identity and covers at least 90% of the scFv and",
        "95% of the shorter of the two sequences. The scFv is labelled `bind` when its best `yes` match counts",
        "and its best `no` match either does not count or has an identity more than 1 percentage point lower.",
        "",
        "## Selection",
        "",
        "Protein mode: each unique amino-acid scFv sequence is one candidate; its reads are the sum over its",
        "nucleotide variants, and its nucleotide sequence is the variant with the most reads. DNA mode: each",
        "unique nucleotide scFv sequence is one candidate.",
        "",
        "Candidates are ranked by TIM-3 reads (barcodes 2 and 3; ties by sequence). Walking down the ranking,",
        "a candidate becomes a representative only if its identity to every representative selected so far is",
        f"below the threshold, up to {TOP_N} representatives. Identity here comes from edlib global alignment",
        "(mode NW), with the same definition as above. Every candidate is then assigned to the representative it",
        "is most identical to (ties go to the higher-ranked representative); the candidates assigned to one",
        "representative form its cluster.",
        "",
        "## Files",
        "",
        "The TSV files hold every candidate with its cluster; the FASTA files hold the representatives only.",
        "`diversity_sweep_summary.tsv` has one row per set: candidates, representatives, their reads, the",
        "identity between representatives and the cluster sizes.",
        "",
        "| TSV | FASTA | mode | identity threshold | candidate rows | representatives | selected reads | selected read fraction |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for (tsv, fasta, rows), summary in zip(files, summaries):
        lines.append(
            f"| `{tsv}` | `{fasta}` | {summary['mode']} | {summary['identity_threshold']:g}% | {rows} "
            f"| {summary['selected_count']} | {summary['selected_reads']} | {summary['selected_read_fraction']} |"
        )
    lines += [
        "",
        "## Main columns",
        "",
        "- `is_representative`: `1` if the candidate is the representative of its cluster.",
        "- `cluster_rank`: rank of the cluster's representative; rows with the same value form one cluster.",
        "- `identity_to_representative`: global identity between the candidate and its representative.",
        "- `cluster_member_count` and `cluster_total_reads`: number of candidates and total reads of the cluster.",
        "- `aa_scfv_sequence` and `nt_scfv_sequence`: candidate sequences.",
        "- `representative_aa_scfv_sequence` and `representative_nt_scfv_sequence`: representative sequences.",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description="Non-redundant sets of TIM-3 strict binders.")
    parser.add_argument("--input-tsv", type=Path, required=True,
                        help="tim3_full_labeled_scfv_results_with_sequences.tsv")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = [row for row in read_tsv(args.input_tsv) if row["final_label"] == "bind"]
    candidates_by_mode = {"protein": build_protein_candidates(rows), "dna": build_dna_candidates(rows)}

    files = []
    summaries = []
    for mode, threshold in SETS:
        print(f"mode={mode} threshold={threshold:g}", flush=True)
        candidates = candidates_by_mode[mode]
        cache = IdentityCache(candidates)
        selected, considered = select_representatives(candidates, threshold, cache)
        assignments = assign_clusters(candidates, selected, cache)
        clusters = describe_clusters(candidates, selected, assignments)

        base = f"tim3_strict_binder_{mode}_identity{threshold:g}"
        clustered_rows = make_clustered_rows(mode, threshold, candidates, selected, assignments, clusters)
        write_tsv(out_dir / f"{base}_clustered_candidates.tsv", clustered_rows)
        write_representatives_fasta(
            out_dir / f"{base}_representatives.fasta", mode, threshold, candidates, selected, clusters
        )
        files.append((f"{base}_clustered_candidates.tsv", f"{base}_representatives.fasta", len(clustered_rows)))
        summaries.append(summarize_set(mode, threshold, candidates, selected, considered, clusters, cache))

    write_tsv(out_dir / "diversity_sweep_summary.tsv", summaries)
    write_readme(out_dir / "README.md", files, summaries)
    print(f"Wrote representative sets to {out_dir}", flush=True)


if __name__ == "__main__":
    main()
