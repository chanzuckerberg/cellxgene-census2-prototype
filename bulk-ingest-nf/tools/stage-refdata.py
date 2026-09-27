#!/usr/bin/env python3
"""Stage the reference data that `--refdata` names: the pinned gene references.

A HealthOmics task has no public internet access, so the writer never fetches a
reference. Run this once, where the internet is reachable, and pass its
destination to the pipeline as `--refdata`.

The files and their SHA-256 digests come from `gene_references` in
`census2_dataset_writer.DEFAULTS`, so this script needs the package installed. It
checks each download against the digest before it writes anything.

    tools/stage-refdata.py refdata/                     # a local directory
    tools/stage-refdata.py s3://bucket/prefix/refdata/  # an S3 prefix, via the AWS CLI
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import tempfile
import urllib.request

import census2_dataset_writer


def fetch(url: str, sha256: str) -> bytes:
    with urllib.request.urlopen(url, timeout=600) as resp:
        data = resp.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        raise SystemExit(f"{url} has digest {digest}; the pin says {sha256}")
    return data


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("dest", help="local directory or s3:// prefix")
    args = p.parse_args()

    references = census2_dataset_writer.DEFAULTS["gene_references"]
    remote = args.dest.startswith("s3://")
    with tempfile.TemporaryDirectory() as scratch:
        local = scratch if remote else args.dest
        os.makedirs(local, exist_ok=True)
        for taxon, spec in references.items():
            print(f"stage {taxon} {spec['version']}: {spec['url']}", file=sys.stderr)
            with open(os.path.join(local, spec["file"]), "wb") as fh:
                fh.write(fetch(spec["url"], spec["sha256"]))
        if remote:
            subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    "--recursive",
                    "--only-show-errors",
                    local,
                    args.dest,
                ],
                check=True,
            )
    print(f"staged {len(references)} gene references in {args.dest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
