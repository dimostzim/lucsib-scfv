#!/usr/bin/env python3
"""Convert FASTQ to FASTA (upper case, 80 bases per line). With --igblast-tsv, reads that IgBLAST aligned in
reverse orientation are reverse-complemented."""

import argparse

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
LINE_LENGTH = 80


def _read_reversed_ids(igblast_tsv):
    with open(igblast_tsv) as fh:
        columns = fh.readline().rstrip("\n").split("\t")
        id_col = columns.index("sequence_id")
        rev_comp_col = columns.index("rev_comp")
        rows = (line.rstrip("\n").split("\t") for line in fh)
        return {row[id_col] for row in rows if row[rev_comp_col] == "T"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--igblast-tsv")
    args = parser.parse_args()

    reversed_ids = _read_reversed_ids(args.igblast_tsv) if args.igblast_tsv else set()
    with open(args.input) as fastq, open(args.output, "w") as fasta:
        for header in fastq:
            seq = next(fastq).strip().upper()
            next(fastq)  # skip the "+" and quality lines
            next(fastq)
            read_id = header[1:].split()[0]
            if read_id in reversed_ids:
                seq = seq.translate(COMPLEMENT)[::-1]
            fasta.write(f">{read_id}\n")
            for start in range(0, len(seq), LINE_LENGTH):
                fasta.write(seq[start : start + LINE_LENGTH] + "\n")


if __name__ == "__main__":
    main()
