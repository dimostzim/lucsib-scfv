#!/usr/bin/env python3
"""Compare TIM-3 Nanopore scFv clones against validated Sanger clones.

The script collapses TIM-3-supported Nanopore rows by exact nt_scFv sequence,
aligns those clones to Sanger references, and writes cutoff sensitivity tables.
It uses edlib for fast semi-global edit-distance alignment.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import edlib
except ImportError as exc:  # pragma: no cover - environment check
    raise SystemExit(
        "Missing dependency: edlib. Use the Python environment from envs/python.yml."
    ) from exc


# Set from the command line in main()
SANGER_TSV: Path
COUNTS_TSV: Path
OUT_DIR: Path
COLLAPSED_TSV: Path
TIM3_BARCODES = ("barcode02", "barcode03")
KMER_SIZE = 15
MIN_SHARED_KMERS = 3
TOP_REFS_PER_LABEL = 5
MAX_EDITS_FOR_NP_MATCH = 100


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_RE = re.compile(r"(\d+)([=XID])")


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
    matches: int
    mismatches: int
    insertions: int
    deletions: int
    sanger_cov: float
    nt_cov: float
    shorter_cov: float
    sanger_len: int
    nt_len: int


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def parse_int(value: str) -> int:
    if value is None or value == "":
        return 0
    return int(float(value))


def pct(num: float, den: float) -> float:
    return 0.0 if den == 0 else 100.0 * num / den


def summarize(values: Sequence[int]) -> Dict[str, object]:
    vals = sorted(values)
    if not vals:
        return {"n": 0}
    return {
        "n": len(vals),
        "min": min(vals),
        "q1": vals[len(vals) // 4],
        "median": median(vals),
        "mean": round(mean(vals), 2),
        "q3": vals[(3 * len(vals)) // 4],
        "max": max(vals),
    }


def kmers(seq: str, k: int = KMER_SIZE) -> Iterable[str]:
    seq = seq.upper()
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            yield kmer


def parse_cigar(cigar: str) -> Tuple[int, int, int, int, int, int, int]:
    matches = mismatches = insertions = deletions = 0
    q_consumed = t_consumed = 0
    for count_s, op in CIGAR_RE.findall(cigar):
        count = int(count_s)
        if op == "=":
            matches += count
            q_consumed += count
            t_consumed += count
        elif op == "X":
            mismatches += count
            q_consumed += count
            t_consumed += count
        elif op == "I":
            insertions += count
            q_consumed += count
        elif op == "D":
            deletions += count
            t_consumed += count
    aligned_len = matches + mismatches + insertions + deletions
    return matches, mismatches, insertions, deletions, q_consumed, t_consumed, aligned_len


def align_query_to_target(query: str, target: str, max_edits: Optional[int] = None) -> Optional[dict]:
    if not query or not target:
        return None
    kwargs = {"mode": "HW", "task": "path"}
    if max_edits is not None:
        kwargs["k"] = max_edits
    result = edlib.align(query, target, **kwargs)
    if result["editDistance"] < 0 or not result.get("cigar"):
        return None
    matches, mismatches, ins, dels, q_cons, t_cons, aln_len = parse_cigar(result["cigar"])
    if aln_len == 0:
        return None
    return {
        "identity": pct(matches, aln_len),
        "edit_distance": int(result["editDistance"]),
        "aligned_len": aln_len,
        "matches": matches,
        "mismatches": mismatches,
        "insertions": ins,
        "deletions": dels,
        "q_cov": pct(q_cons, len(query)),
        "t_cov": pct(t_cons, len(target)),
    }


def best_ref_alignment(
    ref: SangerRef,
    nt_seq: str,
    strands: Sequence[str] = ("+", "-"),
    max_edits: Optional[int] = None,
) -> Optional[Aln]:
    best: Optional[Aln] = None
    nt_len = len(nt_seq)
    strand_sequences = {"+": ref.sequence, "-": revcomp(ref.sequence)}
    for strand in strands:
        sanger_seq = strand_sequences[strand]
        sanger_len = len(sanger_seq)
        if sanger_len <= nt_len:
            direction, query, target = "sanger_in_nt", sanger_seq, nt_seq
        else:
            direction, query, target = "nt_in_sanger", nt_seq, sanger_seq
        raw = align_query_to_target(query, target, max_edits=max_edits)
        if raw is None:
            continue
        if direction == "sanger_in_nt":
            sanger_cov = raw["q_cov"]
            nt_cov = raw["t_cov"]
        else:
            nt_cov = raw["q_cov"]
            sanger_cov = raw["t_cov"]
        shorter_cov = min(sanger_cov, nt_cov) if sanger_len == nt_len else (
            nt_cov if nt_len <= sanger_len else sanger_cov
        )
        aln = Aln(
            clone_number=ref.clone_number,
            binding=ref.binding,
            strand=strand,
            direction=direction,
            identity=raw["identity"],
            edit_distance=raw["edit_distance"],
            aligned_len=raw["aligned_len"],
            matches=raw["matches"],
            mismatches=raw["mismatches"],
            insertions=raw["insertions"],
            deletions=raw["deletions"],
            sanger_cov=sanger_cov,
            nt_cov=nt_cov,
            shorter_cov=shorter_cov,
            sanger_len=sanger_len,
            nt_len=nt_len,
        )
        key = (aln.identity, aln.shorter_cov, aln.nt_cov, -aln.edit_distance)
        best_key = (
            (best.identity, best.shorter_cov, best.nt_cov, -best.edit_distance)
            if best
            else None
        )
        if best is None or key > best_key:
            best = aln
    return best


def load_sanger_refs() -> List[SangerRef]:
    with SANGER_TSV.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        SangerRef(
            clone_number=row["clone_number"],
            sequence=row["sequence"].strip().upper(),
            binding=row["binding"].strip().lower(),
        )
        for row in rows
    ]


def load_tim3_np_clones() -> List[NpClone]:
    if COLLAPSED_TSV.exists():
        clones = []
        with COLLAPSED_TSV.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                clones.append(
                    NpClone(
                        nt_scfv_id=row["nt_scfv_id"],
                        aa_scfv_id=row["aa_scfv_id"],
                        nt_scfv=row["nt_scfv"],
                        aa_scfv=row["aa_scfv"],
                        barcode02=int(row["barcode02"]),
                        barcode03=int(row["barcode03"]),
                        source_rows=int(row["source_rows"]),
                    )
                )
        print(f"Loaded cached collapsed clones: unique_tim3_nt_scfv={len(clones)}", flush=True)
        return sorted(clones, key=lambda c: (-c.tim3_reads, c.nt_scfv_id))

    clones: Dict[str, NpClone] = {}
    total_rows = tim3_rows = 0
    with COUNTS_TSV.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total_rows += 1
            b02 = parse_int(row.get("barcode02", "0"))
            b03 = parse_int(row.get("barcode03", "0"))
            if b02 == 0 and b03 == 0:
                continue
            tim3_rows += 1
            nt_scfv = (row.get("nt_scFv") or "").strip().upper()
            if not nt_scfv or nt_scfv == "NA":
                continue
            key = nt_scfv
            current = clones.get(key)
            if current is None:
                clones[key] = NpClone(
                    nt_scfv_id=(row.get("nt_scFv_id") or "").strip(),
                    aa_scfv_id=(row.get("aa_scFv_id") or "").strip(),
                    nt_scfv=nt_scfv,
                    aa_scfv=(row.get("aa_scFv") or "").strip(),
                    barcode02=b02,
                    barcode03=b03,
                    source_rows=1,
                )
            else:
                current.barcode02 += b02
                current.barcode03 += b03
                current.source_rows += 1
    print(
        f"Parsed counts: total_rows={total_rows} tim3_rows={tim3_rows} "
        f"unique_tim3_nt_scfv={len(clones)}",
        flush=True,
    )
    rows = [
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
        for clone in sorted(clones.values(), key=lambda c: (-c.tim3_reads, c.nt_scfv_id))
    ]
    write_tsv(COLLAPSED_TSV, rows)
    return sorted(clones.values(), key=lambda c: (-c.tim3_reads, c.nt_scfv_id))


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
    no_candidate = 0
    for idx, clone in enumerate(clones, start=1):
        candidates = candidate_ref_strands(clone.nt_scfv, refs_by_id, kmer_index)
        if not candidates:
            no_candidate += 1
            rec = empty_match_record(clone)
            records.append(rec)
            continue

        # candidates come ordered by shared k-mers, then clone number; equally good alignments are
        # resolved in favour of the earlier candidate
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
        if not alignments:
            rec = empty_match_record(clone)
            rec["candidate_sanger_refs"] = len(candidates)
            records.append(rec)
            continue
        best = max(alignments, key=aln_sort_key)
        yes = [aln for aln in alignments if aln.binding == "yes"]
        no = [aln for aln in alignments if aln.binding == "no"]
        best_yes = max(yes, key=aln_sort_key) if yes else None
        best_no = max(no, key=aln_sort_key) if no else None

        rec = empty_match_record(clone)
        rec["candidate_sanger_refs"] = len(candidates)
        fill_aln(rec, "best", best)
        fill_aln(rec, "best_yes", best_yes)
        fill_aln(rec, "best_no", best_no)
        if best_yes and best_no:
            rec["yes_minus_no_identity"] = round(best_yes.identity - best_no.identity, 4)
        records.append(rec)
        if idx % 10000 == 0:
            print(f"Matched {idx}/{len(clones)} clones", flush=True)
    print(f"Matching complete: no_candidate={no_candidate}", flush=True)
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
        if aln is None:
            raise RuntimeError(f"No Sanger/Sanger alignment for clone {a.clone_number}/{b.clone_number}")
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


def write_tsv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    refs: Sequence[SangerRef],
    clones: Sequence[NpClone],
    match_records: Sequence[Dict[str, object]],
    sanger_pairs: Sequence[Dict[str, object]],
    sensitivity: Sequence[Dict[str, object]],
) -> None:
    sanger_lengths = [len(ref.sequence) for ref in refs]
    nt_lengths = [len(clone.nt_scfv) for clone in clones]
    binding_counts = Counter(ref.binding for ref in refs)
    no_candidate = sum(1 for r in match_records if int(r["candidate_sanger_refs"]) == 0)
    best_identities = [
        float(r["best_identity"]) for r in match_records if r["best_identity"] != ""
    ]
    opposite_pairs = [r for r in sanger_pairs if r["same_label"] is False]
    closest_opposite = sorted(opposite_pairs, key=lambda r: float(r["identity"]), reverse=True)[:10]

    chosen_rows = [
        row
        for row in sensitivity
        if row["identity_min"] == 99.0
        and row["nt_cov_min"] == 90.0
        and row["shorter_cov_min"] == 95.0
        and row["margin_min"] == 1.0
    ]
    chosen = chosen_rows[0] if chosen_rows else {}

    lines = []
    lines.append("# TIM-3 Sanger Cutoff Analysis")
    lines.append("")
    lines.append("## Inputs")
    lines.append(f"- Sanger TSV: `{SANGER_TSV}`")
    lines.append(f"- Nanopore counts TSV: `{COUNTS_TSV}`")
    lines.append(f"- TIM-3 barcodes: `{', '.join(TIM3_BARCODES)}`")
    lines.append("")
    lines.append("## Basic Counts")
    lines.append(f"- Sanger clones: {len(refs)} ({dict(binding_counts)})")
    lines.append(f"- Unique TIM-3 Nanopore nt_scFv clones: {len(clones)}")
    lines.append(f"- TIM-3 read support after nt_scFv collapse: {sum(c.tim3_reads for c in clones)}")
    lines.append(f"- Nanopore clones with no k-mer candidate Sanger reference: {no_candidate}")
    lines.append("")
    lines.append("## Lengths")
    lines.append(f"- Sanger sequence lengths: `{summarize(sanger_lengths)}`")
    lines.append(f"- Nanopore nt_scFv lengths: `{summarize(nt_lengths)}`")
    lines.append("")
    lines.append("## Match Score Overview")
    lines.append(f"- Best-match identity summary: `{summarize([round(x) for x in best_identities])}`")
    lines.append("")
    lines.append("## Closest Opposite-Label Sanger Pairs")
    lines.append("| clone_a | label_a | clone_b | label_b | identity | shorter_cov |")
    lines.append("|---|---|---|---|---:|---:|")
    for row in closest_opposite:
        lines.append(
            f"| {row['clone_a']} | {row['binding_a']} | {row['clone_b']} | "
            f"{row['binding_b']} | {row['identity']} | {row['shorter_cov']} |"
        )
    lines.append("")
    lines.append("## Example Conservative Cutoff")
    lines.append(
        "Conservative example: identity >= 99%, nt_scFv coverage >= 90%, "
        "shorter-sequence coverage >= 95%, opposite-label margin > 1%."
    )
    if chosen:
        lines.append(f"- bind clones/reads: {chosen['bind_clones']} / {chosen['bind_reads']}")
        lines.append(f"- no_bind clones/reads: {chosen['no_bind_clones']} / {chosen['no_bind_reads']}")
        lines.append(f"- ambiguous clones/reads: {chosen['ambiguous_clones']} / {chosen['ambiguous_reads']}")
        lines.append(f"- unknown clones/reads: {chosen['unknown_clones']} / {chosen['unknown_reads']}")
    lines.append("")
    lines.append("## Output Files")
    lines.append("- `tim3_sanger_vs_sanger.tsv`")
    lines.append("- `tim3_nanopore_sanger_best_matches.tsv`")
    lines.append("- `tim3_cutoff_sensitivity.tsv`")
    lines.append("- `tim3_conservative_labels.tsv`")
    lines.append("")
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    global SANGER_TSV, COUNTS_TSV, OUT_DIR, COLLAPSED_TSV
    parser = argparse.ArgumentParser(description="Match TIM-3 Nanopore scFvs to the Sanger bind/no-bind clones.")
    parser.add_argument("--sanger-tsv", type=Path, required=True)
    parser.add_argument("--counts-tsv", type=Path, required=True, help="step-03 counts table")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    SANGER_TSV, COUNTS_TSV, OUT_DIR = args.sanger_tsv, args.counts_tsv, args.out_dir
    COLLAPSED_TSV = OUT_DIR / "tim3_nanopore_nt_scfv_collapsed.tsv"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    refs = load_sanger_refs()
    clones = load_tim3_np_clones()

    sanger_pairs = sanger_vs_sanger(refs)
    write_tsv(OUT_DIR / "tim3_sanger_vs_sanger.tsv", sanger_pairs)

    match_records = match_np_to_sanger(clones, refs)
    write_tsv(OUT_DIR / "tim3_nanopore_sanger_best_matches.tsv", match_records)

    sensitivity = cutoff_sensitivity(match_records)
    write_tsv(OUT_DIR / "tim3_cutoff_sensitivity.tsv", sensitivity)

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
    write_tsv(OUT_DIR / "tim3_conservative_labels.tsv", conservative)

    write_summary(refs, clones, match_records, sanger_pairs, sensitivity)
    print(f"Wrote outputs to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
