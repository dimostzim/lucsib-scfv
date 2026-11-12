#!/usr/bin/env python3
"""Create a restricted Zenodo draft and upload files to it (maintainer tool, not part of the pipeline).

The access token is read from the ZENODO_ACCESS_TOKEN environment variable.
The record is never published here; publish it from the Zenodo web page.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.request
from pathlib import Path

API = "https://zenodo.org/api"


def _request(method: str, url: str, token: str, body: bytes | None = None, headers: dict | None = None,
             data_stream=None, length: int | None = None) -> dict:
    hdrs = {"Authorization": f"Bearer {token}"}
    hdrs.update(headers or {})
    if data_stream is not None:
        hdrs["Content-Length"] = str(length)
        req = urllib.request.Request(url, data=data_stream, headers=hdrs, method=method)
    else:
        req = urllib.request.Request(url, data=body, headers=hdrs, method=method)
    with urllib.request.urlopen(req, timeout=3600) as resp:
        text = resp.read().decode()
    return json.loads(text) if text else {}


def _json(method: str, url: str, token: str, payload: dict | None = None) -> dict:
    body = json.dumps(payload).encode() if payload is not None else None
    return _request(method, url, token, body=body, headers={"Content-Type": "application/json"})


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def create(args: argparse.Namespace, token: str) -> None:
    dep = _json("POST", f"{API}/deposit/depositions", token, {})
    metadata = {
        "upload_type": "dataset",
        "title": args.title,
        "description": args.description,
        "access_right": "restricted",
        "access_conditions": args.access_conditions,
    }
    _json("PUT", f"{API}/deposit/depositions/{dep['id']}", token, {"metadata": metadata})
    print(f"Created draft {dep['id']}: {dep['links']['html']}")


def upload(args: argparse.Namespace, token: str) -> None:
    dep = _json("GET", f"{API}/deposit/depositions/{args.record_id}", token)
    bucket = dep["links"]["bucket"]
    existing = {f["filename"]: f["checksum"] for f in dep.get("files", [])}
    for path in map(Path, args.files):
        md5 = _md5(path)
        if existing.get(path.name) == md5:
            print(f"{path.name}: already uploaded, md5 matches", flush=True)
            continue
        size = path.stat().st_size
        print(f"{path.name}: uploading {size / 1e9:.2f} GB", flush=True)
        with path.open("rb") as handle:
            result = _request("PUT", f"{bucket}/{path.name}", token,
                              headers={"Content-Type": "application/octet-stream"},
                              data_stream=handle, length=size)
        remote_md5 = result["checksum"].removeprefix("md5:")
        if remote_md5 != md5:
            sys.exit(f"{path.name}: checksum mismatch after upload (local {md5}, Zenodo {remote_md5})")
        print(f"{path.name}: uploaded, md5 {md5}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("create")
    c.add_argument("--title", required=True)
    c.add_argument("--description", required=True)
    c.add_argument("--access-conditions", required=True)
    u = sub.add_parser("upload")
    u.add_argument("--record-id", required=True)
    u.add_argument("files", nargs="+")
    args = p.parse_args()
    token = os.environ["ZENODO_ACCESS_TOKEN"]
    {"create": create, "upload": upload}[args.command](args, token)


if __name__ == "__main__":
    main()
