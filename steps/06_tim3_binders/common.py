#!/usr/bin/env python3
"""Alignment and TSV helpers shared by the TIM-3 scripts."""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

import edlib

# full reads can be longer than the csv module's default field limit
csv.field_size_limit(sys.maxsize)

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_RE = re.compile(r"(\d+)([=XID])")


@dataclass
class Alignment:
    identity: float
    edit_distance: int
    aligned_len: int
    query_cov: float
    target_cov: float


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def align(query: str, target: str, mode: str, max_edits: int = -1) -> Optional[Alignment]:
    """edlib alignment; identity is matching columns over all alignment columns. None if over max_edits."""
    result = edlib.align(query, target, mode=mode, task="path", k=max_edits)
    if result["editDistance"] < 0:
        return None
    matches = aligned_len = q_consumed = t_consumed = 0
    for count_s, op in CIGAR_RE.findall(result["cigar"]):
        count = int(count_s)
        aligned_len += count
        if op == "=":
            matches += count
        if op in "=XI":
            q_consumed += count
        if op in "=XD":
            t_consumed += count
    return Alignment(
        identity=100.0 * matches / aligned_len,
        edit_distance=int(result["editDistance"]),
        aligned_len=aligned_len,
        query_cov=100.0 * q_consumed / len(query),
        target_cov=100.0 * t_consumed / len(target),
    )


def read_tsv(path: Path) -> List[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
