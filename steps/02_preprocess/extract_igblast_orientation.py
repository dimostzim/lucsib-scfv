#!/usr/bin/env python3
"""Extract per-read orientation calls from IgBLAST AIRR output."""

import argparse
import csv
import sys
from pathlib import Path


def parse_bool(value: str):
    v = (value or "").strip().lower()
    if v in {"t", "true", "1", "yes"}:
        return True
    if v in {"f", "false", "0", "no"}:
        return False
    return None


def main():
    parser = argparse.ArgumentParser(description="Build keep/rc ID lists from IgBLAST AIRR TSV.")
    parser.add_argument("--igblast-tsv", required=True, help="IgBLAST AIRR TSV (outfmt 19)")
    parser.add_argument("--keep-out", required=True, help="Output file with oriented read IDs to keep")
    parser.add_argument("--rc-out", required=True, help="Output file with read IDs to reverse-complement")
    parser.add_argument("--ambig-out", default="", help="Optional file with ambiguous read IDs")
    args = parser.parse_args()

    # AIRR rows can contain long alignment fields; raise CSV parser limit.
    max_csv = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_csv)
            break
        except OverflowError:
            max_csv //= 10

    states = {}
    with Path(args.igblast_tsv).open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sequence_id", "rev_comp"}
        if not required.issubset(set(reader.fieldnames or [])):
            missing = sorted(required - set(reader.fieldnames or []))
            raise SystemExit(f"Missing required columns in IgBLAST TSV: {', '.join(missing)}")

        for row in reader:
            read_id = (row.get("sequence_id") or "").split("|")[0].strip()
            if not read_id:
                continue
            rev = parse_bool(row.get("rev_comp", ""))
            if rev is None:
                continue
            states.setdefault(read_id, set()).add(rev)

    keep_ids = []
    rc_ids = []
    ambig_ids = []
    for read_id, vals in states.items():
        if vals == {True}:
            keep_ids.append(read_id)
            rc_ids.append(read_id)
        elif vals == {False}:
            keep_ids.append(read_id)
        else:
            ambig_ids.append(read_id)

    keep_ids.sort()
    rc_ids.sort()
    ambig_ids.sort()

    Path(args.keep_out).write_text("".join(f"{x}\n" for x in keep_ids))
    Path(args.rc_out).write_text("".join(f"{x}\n" for x in rc_ids))
    if args.ambig_out:
        Path(args.ambig_out).write_text("".join(f"{x}\n" for x in ambig_ids))

    print(
        f"{Path(args.igblast_tsv).name}: keep={len(keep_ids)} "
        f"rc={len(rc_ids)} ambig={len(ambig_ids)}"
    )


if __name__ == "__main__":
    main()
