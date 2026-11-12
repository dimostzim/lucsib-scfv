#!/usr/bin/env python3
"""Convert FASTQ to FASTA with optional orientation and size filter."""

import argparse
from pathlib import Path
from typing import Optional, Set

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def load_id_set(path: str) -> Optional[Set[str]]:
    if not path:
        return None
    ids = set()
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                ids.add(line.split()[0])
    return ids


def main():
    parser = argparse.ArgumentParser(description="FASTQ -> FASTA with size filter and optional reverse-complement.")
    parser.add_argument("--input", required=True, help="Input FASTQ")
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument("--min-len", type=int, default=0, help="Minimum length to keep (0 = no filter)")
    parser.add_argument("--max-len", type=int, default=0, help="Maximum length to keep (0 = no filter)")
    parser.add_argument(
        "--orient",
        choices=["none", "rc", "anchor"],
        default="none",
        help="Orientation mode: none (keep), rc (reverse-complement all), anchor (orient by anchor)",
    )
    parser.add_argument(
        "--anchor",
        default="",
        help="Anchor sequence for --orient anchor (exact match, case-insensitive)",
    )
    parser.add_argument(
        "--rc-list",
        default="",
        help="Optional text file with read IDs to reverse-complement",
    )
    parser.add_argument(
        "--keep-list",
        default="",
        help="Optional text file with read IDs to keep (others are dropped)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    orient_rc = args.orient == "rc"
    orient_anchor = args.orient == "anchor"
    anchor = args.anchor.upper() if args.anchor else ""
    if orient_anchor and not anchor:
        raise SystemExit("--orient anchor requires --anchor sequence")
    rc_ids = load_id_set(args.rc_list)
    keep_ids = load_id_set(args.keep_list)

    kept = 0
    total = 0
    oriented_fwd = 0
    oriented_rc = 0
    ambiguous = 0
    no_anchor = 0
    dropped_not_kept = 0

    with input_path.open() as fh, output_path.open("w") as out:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            total += 1
            read_id = header[1:].strip().split()[0]
            if keep_ids is not None and read_id not in keep_ids:
                dropped_not_kept += 1
                continue
            seq = seq.upper()
            if orient_rc:
                seq = revcomp(seq)
                oriented_rc += 1
            elif orient_anchor:
                hit_fwd = anchor in seq
                hit_rc = anchor in revcomp(seq)
                if hit_fwd and not hit_rc:
                    oriented_fwd += 1
                elif hit_rc and not hit_fwd:
                    seq = revcomp(seq)
                    oriented_rc += 1
                elif hit_fwd and hit_rc:
                    ambiguous += 1
                else:
                    no_anchor += 1
            if rc_ids is not None and read_id in rc_ids:
                seq = revcomp(seq)
                oriented_rc += 1
            if args.min_len and len(seq) < args.min_len:
                continue
            if args.max_len and len(seq) > args.max_len:
                continue
            kept += 1
            out.write(f">{read_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")

    if orient_anchor:
        print(
            f"{input_path.name}: kept {kept}/{total} | "
            f"anchor_fwd {oriented_fwd} anchor_rc {oriented_rc} "
            f"ambig {ambiguous} no_anchor {no_anchor}"
        )
    else:
        print(f"{input_path.name}: kept {kept}/{total}")
    if keep_ids is not None:
        print(f"{input_path.name}: dropped_not_kept {dropped_not_kept}")


if __name__ == "__main__":
    main()
