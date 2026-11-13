#!/usr/bin/env python3
"""Download the input files from Zenodo, check their checksums and unpack the reads.

Writes to OUT_DIR:
  downloads/<file>          the files as stored on Zenodo
  reads/barcodeXX.fastq     decompressed reads, input of step 01
  samples.csv, TIM-3_sanger_binding.tsv
Files already in OUT_DIR/downloads with the expected checksum are not downloaded. Downloading needs a
Zenodo access token in the ZENODO_ACCESS_TOKEN environment variable.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import sys
import urllib.request
from pathlib import Path

API = "https://zenodo.org/api"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, token: str, dest: Path) -> None:
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(req, timeout=3600) as resp, tmp.open("wb") as out:
        shutil.copyfileobj(resp, out, length=1 << 24)
    tmp.rename(dest)


def _file_url(record_id: str, name: str) -> str:
    return f"{API}/records/{record_id}/files/{name}/content"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--record-id", required=True)
    p.add_argument("--manifest", type=Path, required=True, help="TSV with file and sha256 columns")
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    token = os.environ.get("ZENODO_ACCESS_TOKEN", "")

    with args.manifest.open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    downloads = args.out_dir / "downloads"
    reads = args.out_dir / "reads"
    downloads.mkdir(parents=True, exist_ok=True)
    reads.mkdir(parents=True, exist_ok=True)

    for i, row in enumerate(manifest, start=1):
        name, expected = row["file"], row["sha256"]
        dest = downloads / name
        if dest.exists() and _sha256(dest) == expected:
            print(f"PROGRESS {i}/{len(manifest)} {name}: already downloaded", flush=True)
        else:
            if not token:
                sys.exit(f"{name} is missing from {downloads} or has a different checksum, and ZENODO_ACCESS_TOKEN "
                         "is not set: put a Zenodo token in .env or the input files in that folder")
            print(f"PROGRESS {i}/{len(manifest)} {name}: downloading", flush=True)
            _download(_file_url(args.record_id, name), token, dest)
            got = _sha256(dest)
            if got != expected:
                sys.exit(f"{name}: sha256 {got} does not match manifest {expected}")
        if name.endswith(".fastq.gz"):
            target = reads / name[: -len(".gz")]
            tmp = target.with_name(target.name + ".part")
            with gzip.open(dest, "rb") as src, tmp.open("wb") as out:
                shutil.copyfileobj(src, out, length=1 << 24)
            tmp.rename(target)
        else:
            shutil.copyfile(dest, args.out_dir / name)
    print(f"Inputs ready in {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
