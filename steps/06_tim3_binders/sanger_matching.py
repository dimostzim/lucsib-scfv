#!/usr/bin/env python3
"""Compare TIM-3 Nanopore scFv clones against validated Sanger clones.

The script collapses TIM-3-supported Nanopore rows by exact nt_scFv sequence,
aligns those clones to Sanger references, and writes cutoff sensitivity tables
and the first-pass labels. Alignments are edlib infix alignments (mode HW).
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from common import align, revcomp, write_tsv


KMER_SIZE = 15
MIN_SHARED_KMERS = 3
TOP_REFS_PER_LABEL = 5
MAX_EDITS_FOR_NP_MATCH = 100


@dataclass
class SangerRef:
    clone_number: str
    sequence: str
    binding: str


@dataclass
class NpClone:
    nt_scfv_id: str
    aa_scfv_id: str
    nt_scfv: str
    aa_scfv: str
    barcode02: int
    barcode03: int
    source_rows: int

    @property
    def tim3_reads(self) -> int:
        return self.barcode02 + self.barcode03


@dataclass
class Aln:
    clone_number: str
    binding: str
    strand: str
    direction: str
    identity: float
    edit_distance: int
    aligned_len: int
    sanger_cov: float
    nt_cov: float
    shorter_cov: float
    sanger_len: int
    nt_len: int


def kmers(seq: str, k: int = KMER_SIZE) -> Iterable[str]:
    seq = seq.upper()
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            yield kmer


def best_ref_alignment(
    ref: SangerRef,
    nt_seq: str,
    strands: Sequence[str] = ("+", "-"),
    max_edits: int = -1,
) -> Optional[Aln]:
    best: Optional[Aln] = None
    nt_len = len(nt_seq)
    strand_sequences = {"+": ref.sequence, "-": revcomp(ref.sequence)}
    for strand in strands:
        sanger_seq = strand_sequences[strand]
        sanger_len = len(sanger_seq)
        # the shorter sequence is aligned as an infix of the longer one
        if sanger_len <= nt_len:
            direction = "sanger_in_nt"
            raw = align(sanger_seq, nt_seq, mode="HW", max_edits=max_edits)
        else:
            direction = "nt_in_sanger"
            raw = align(nt_seq, sanger_seq, mode="HW", max_edits=max_edits)
        if raw is None:
            continue
        if direction == "sanger_in_nt":
            sanger_cov = raw.query_cov
            nt_cov = raw.target_cov
        else:
            nt_cov = raw.query_cov
            sanger_cov = raw.target_cov
        shorter_cov = min(sanger_cov, nt_cov) if sanger_len == nt_len else (
            nt_cov if nt_len <= sanger_len else sanger_cov
        )
        aln = Aln(
            clone_number=ref.clone_number,
            binding=ref.binding,
            strand=strand,
            direction=direction,
            identity=raw.identity,
            edit_distance=raw.edit_distance,
            aligned_len=raw.aligned_len,
            sanger_cov=sanger_cov,
            nt_cov=nt_cov,
            shorter_cov=shorter_cov,
            sanger_len=sanger_len,
            nt_len=nt_len,
        )
        if best is None or aln_sort_key(aln) > aln_sort_key(best):
            best = aln
    return best


def load_sanger_refs(sanger_tsv: Path) -> List[SangerRef]:
    with sanger_tsv.open() as handle:
        return [
            SangerRef(clone_number=row["clone_number"], sequence=row["sequence"], binding=row["binding"])
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def load_tim3_np_clones(counts_tsv: Path, collapsed_tsv: Path) -> List[NpClone]:
    clones: Dict[str, NpClone] = {}
    with counts_tsv.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            b02 = int(row["barcode02"])
            b03 = int(row["barcode03"])
            if b02 == 0 and b03 == 0:
                continue
            nt_scfv = row["nt_scFv"]
            current = clones.get(nt_scfv)
            if current is None:
                clones[nt_scfv] = NpClone(
                    nt_scfv_id=row["nt_scFv_id"],
                    aa_scfv_id=row["aa_scFv_id"],
                    nt_scfv=nt_scfv,
                    aa_scfv=row["aa_scFv"],
                    barcode02=b02,
                    barcode03=b03,
                    source_rows=1,
                )
            else:
                current.barcode02 += b02
                current.barcode03 += b03
                current.source_rows += 1
    print(f"Unique TIM-3 nt_scFv clones: {len(clones)}", flush=True)
    ordered = sorted(clones.values(), key=lambda c: (-c.tim3_reads, c.nt_scfv_id))
    write_tsv(
        collapsed_tsv,
        [
            {
                "nt_scfv_id": clone.nt_scfv_id,
                "aa_scfv_id": clone.aa_scfv_id,
                "nt_scfv": clone.nt_scfv,
                "aa_scfv": clone.aa_scfv,
                "barcode02": clone.barcode02,
                "barcode03": clone.barcode03,
                "tim3_reads": clone.tim3_reads,
                "source_rows": clone.source_rows,
            }
            for clone in ordered
        ],
    )
    return ordered


def build_ref_kmer_index(refs: Sequence[SangerRef]) -> Dict[str, set]:
    index: Dict[str, set] = defaultdict(set)
    for ref in refs:
        for strand, strand_seq in (("+", ref.sequence), ("-", revcomp(ref.sequence))):
            for kmer in set(kmers(strand_seq)):
                index[kmer].add((ref.clone_number, strand))
    return index


def candidate_ref_strands(
    nt_seq: str,
    refs_by_id: Dict[str, SangerRef],
    kmer_index: Dict[str, set],
) -> Dict[str, List[str]]:
    counts: Counter = Counter()
    for kmer in set(kmers(nt_seq)):
        counts.update(kmer_index.get(kmer, ()))
    by_ref: Dict[str, Tuple[int, str]] = {}
    # iterate in a fixed order so that ties do not depend on Python's random set ordering
    for (clone_number, strand), shared in sorted(counts.items(), key=lambda kv: (int(kv[0][0]), kv[0][1])):
        if shared < MIN_SHARED_KMERS:
            continue
        prev = by_ref.get(clone_number)
        if prev is None or shared > prev[0]:
            by_ref[clone_number] = (shared, strand)

    by_label: Dict[str, List[Tuple[str, int, str]]] = defaultdict(list)
    for clone_number, (shared, strand) in sorted(by_ref.items(), key=lambda kv: int(kv[0])):
        by_label[refs_by_id[clone_number].binding].append((clone_number, shared, strand))

    chosen = []
    for label_refs in by_label.values():
        chosen.extend(sorted(label_refs, key=lambda x: (-x[1], int(x[0])))[:TOP_REFS_PER_LABEL])
    # top candidates per label by shared k-mers, ties by clone number; all candidates in that order
    selected: Dict[str, List[str]] = {}
    for clone_number, _shared, strand in sorted(chosen, key=lambda x: (-x[1], int(x[0]))):
        selected.setdefault(clone_number, []).append(strand)
    return selected


def empty_match_record(clone: NpClone) -> Dict[str, object]:
    record = {
        "nt_scfv_id": clone.nt_scfv_id,
        "aa_scfv_id": clone.aa_scfv_id,
        "nt_scfv_len": len(clone.nt_scfv),
        "aa_scfv_len": len(clone.aa_scfv),
        "barcode02": clone.barcode02,
        "barcode03": clone.barcode03,
        "tim3_reads": clone.tim3_reads,
        "source_rows": clone.source_rows,
        "candidate_sanger_refs": 0,
    }
    for prefix in ("best", "best_yes", "best_no"):
        record.update(
            {
                f"{prefix}_clone": "",
                f"{prefix}_binding": "",
                f"{prefix}_identity": "",
                f"{prefix}_edit_distance": "",
                f"{prefix}_aligned_len": "",
                f"{prefix}_sanger_cov": "",
                f"{prefix}_nt_cov": "",
                f"{prefix}_shorter_cov": "",
                f"{prefix}_sanger_len": "",
                f"{prefix}_strand": "",
                f"{prefix}_direction": "",
            }
        )
    record["yes_minus_no_identity"] = ""
    return record


def fill_aln(record: Dict[str, object], prefix: str, aln: Optional[Aln]) -> None:
    if aln is None:
        return
    record.update(
        {
            f"{prefix}_clone": aln.clone_number,
            f"{prefix}_binding": aln.binding,
            f"{prefix}_identity": round(aln.identity, 4),
            f"{prefix}_edit_distance": aln.edit_distance,
            f"{prefix}_aligned_len": aln.aligned_len,
            f"{prefix}_sanger_cov": round(aln.sanger_cov, 4),
            f"{prefix}_nt_cov": round(aln.nt_cov, 4),
            f"{prefix}_shorter_cov": round(aln.shorter_cov, 4),
            f"{prefix}_sanger_len": aln.sanger_len,
            f"{prefix}_strand": aln.strand,
            f"{prefix}_direction": aln.direction,
        }
    )


def aln_sort_key(aln: Aln) -> Tuple[float, float, float, int]:
    return (aln.identity, aln.shorter_cov, aln.nt_cov, -aln.edit_distance)


def match_np_to_sanger(clones: Sequence[NpClone], refs: Sequence[SangerRef]) -> List[Dict[str, object]]:
    refs_by_id = {ref.clone_number: ref for ref in refs}
    kmer_index = build_ref_kmer_index(refs)
    records = []
    for idx, clone in enumerate(clones, start=1):
        rec = empty_match_record(clone)
        candidates = candidate_ref_strands(clone.nt_scfv, refs_by_id, kmer_index)
        rec["candidate_sanger_refs"] = len(candidates)
        # candidates come ordered by shared k-mers, then clone number; equally good alignments are
        # resolved in favour of the candidate that comes first in this order
        alignments = [
            aln
            for clone_number, strands in candidates.items()
            if (
                aln := best_ref_alignment(
                    refs_by_id[clone_number],
                    clone.nt_scfv,
                    strands=strands,
                    max_edits=MAX_EDITS_FOR_NP_MATCH,
                )
            )
            is not None
        ]
        if alignments:
            yes = [aln for aln in alignments if aln.binding == "yes"]
            no = [aln for aln in alignments if aln.binding == "no"]
            best_yes = max(yes, key=aln_sort_key) if yes else None
            best_no = max(no, key=aln_sort_key) if no else None
            fill_aln(rec, "best", max(alignments, key=aln_sort_key))
            fill_aln(rec, "best_yes", best_yes)
            fill_aln(rec, "best_no", best_no)
            if best_yes and best_no:
                rec["yes_minus_no_identity"] = round(best_yes.identity - best_no.identity, 4)
        records.append(rec)
        if idx % 10000 == 0:
            print(f"Matched {idx}/{len(clones)} clones", flush=True)
    return records


def call_record(
    record: Dict[str, object],
    identity_min: float,
    nt_cov_min: float,
    shorter_cov_min: float,
    margin_min: float,
) -> str:
    passing = {}
    for label, prefix in (("bind", "best_yes"), ("no_bind", "best_no")):
        identity = record.get(f"{prefix}_identity")
        nt_cov = record.get(f"{prefix}_nt_cov")
        shorter_cov = record.get(f"{prefix}_shorter_cov")
        if identity == "" or nt_cov == "" or shorter_cov == "":
            continue
        if (
            float(identity) >= identity_min
            and float(nt_cov) >= nt_cov_min
            and float(shorter_cov) >= shorter_cov_min
        ):
            passing[label] = float(identity)
    if not passing:
        return "unknown"
    if len(passing) == 1:
        return next(iter(passing))
    diff = passing["bind"] - passing["no_bind"]
    if abs(diff) <= margin_min:
        return "ambiguous"
    return "bind" if diff > 0 else "no_bind"


def cutoff_sensitivity(records: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    rows = []
    for identity_min, nt_cov_min, shorter_cov_min, margin_min in itertools.product(
        (95.0, 97.0, 98.0, 99.0, 99.5),
        (50.0, 70.0, 80.0, 90.0, 95.0),
        (90.0, 95.0, 99.0),
        (0.0, 1.0, 2.0),
    ):
        clone_counts = Counter()
        read_counts = Counter()
        for record in records:
            label = call_record(record, identity_min, nt_cov_min, shorter_cov_min, margin_min)
            clone_counts[label] += 1
            read_counts[label] += int(record["tim3_reads"])
        rows.append(
            {
                "identity_min": identity_min,
                "nt_cov_min": nt_cov_min,
                "shorter_cov_min": shorter_cov_min,
                "margin_min": margin_min,
                "bind_clones": clone_counts["bind"],
                "no_bind_clones": clone_counts["no_bind"],
                "ambiguous_clones": clone_counts["ambiguous"],
                "unknown_clones": clone_counts["unknown"],
                "bind_reads": read_counts["bind"],
                "no_bind_reads": read_counts["no_bind"],
                "ambiguous_reads": read_counts["ambiguous"],
                "unknown_reads": read_counts["unknown"],
            }
        )
    return rows


def sanger_vs_sanger(refs: Sequence[SangerRef]) -> List[Dict[str, object]]:
    rows = []
    for a, b in itertools.combinations(refs, 2):
        aln = best_ref_alignment(a, b.sequence, strands=("+",))
        rows.append(
            {
                "clone_a": a.clone_number,
                "binding_a": a.binding,
                "len_a": len(a.sequence),
                "clone_b": b.clone_number,
                "binding_b": b.binding,
                "len_b": len(b.sequence),
                "same_label": a.binding == b.binding,
                "identity": round(aln.identity, 4),
                "edit_distance": aln.edit_distance,
                "aligned_len": aln.aligned_len,
                "clone_a_cov": round(aln.sanger_cov, 4),
                "clone_b_cov": round(aln.nt_cov, 4),
                "shorter_cov": round(aln.shorter_cov, 4),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Match TIM-3 Nanopore scFvs to the Sanger bind/no-bind clones.")
    parser.add_argument("--sanger-tsv", type=Path, required=True)
    parser.add_argument("--counts-tsv", type=Path, required=True, help="step-03 counts table")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out_dir = args.out_dir

    out_dir.mkdir(parents=True, exist_ok=True)
    refs = load_sanger_refs(args.sanger_tsv)
    clones = load_tim3_np_clones(args.counts_tsv, out_dir / "tim3_nanopore_nt_scfv_collapsed.tsv")

    write_tsv(out_dir / "tim3_sanger_vs_sanger.tsv", sanger_vs_sanger(refs))

    match_records = match_np_to_sanger(clones, refs)
    write_tsv(out_dir / "tim3_cutoff_sensitivity.tsv", cutoff_sensitivity(match_records))

    conservative = []
    for record in match_records:
        label = call_record(
            record,
            identity_min=99.0,
            nt_cov_min=90.0,
            shorter_cov_min=95.0,
            margin_min=1.0,
        )
        out = dict(record)
        out["assigned_label"] = label
        out["cutoff"] = "identity>=99;nt_cov>=90;shorter_cov>=95;margin>1"
        conservative.append(out)
    write_tsv(out_dir / "tim3_conservative_labels.tsv", conservative)
    print(f"Wrote outputs to {out_dir}", flush=True)


if __name__ == "__main__":
    main()
