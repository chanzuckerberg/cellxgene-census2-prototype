#!/usr/bin/env python3
"""Seed a mock S3 with the fixture datasets, fetched from the public CxG CDN.

The pipeline reads a source with `aws s3 cp`, so a test run needs each object
under the bucket and the key that the manifest names. This script puts them
there, which lets the ingest tests run with no AWS account. The CDN serves the
same objects as the source bucket, and anonymously.

Set `AWS_ENDPOINT_URL` to the mock endpoint before you run this. Any
credentials will do, because a mock accepts any value.

Uses the Python standard library and the AWS CLI, so it runs outside the
runtime container. See .github/workflows/bulk-ingest-nf.yml for how CI uses it.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request


def log(msg: str) -> None:
    print(f"seed: {msg}", file=sys.stderr, flush=True)


def bucket_of(uri: str) -> str:
    """Bucket name from an s3:// URI."""
    return uri.split("/")[2]


def fetch(url: str, dest: str) -> int:
    """Download `url` to `dest`. Return the byte count."""
    with urllib.request.urlopen(url, timeout=300) as resp, open(dest, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    return os.path.getsize(dest)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest", nargs="?", default="tests/fixtures/manifest.jsonl")
    args = p.parse_args()

    if not os.environ.get("AWS_ENDPOINT_URL"):
        raise SystemExit("set AWS_ENDPOINT_URL to the mock endpoint first")

    with open(args.manifest) as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    if not records:
        raise SystemExit(f"no records in {args.manifest}")

    for bucket in sorted({bucket_of(r["uri"]) for r in records}):
        # Creating a bucket that exists is not an error worth stopping for, so
        # the check comes from the upload instead.
        subprocess.run(
            ["aws", "s3", "mb", f"s3://{bucket}"],
            check=False,
            capture_output=True,
        )
        log(f"bucket {bucket}")

    for r in records:
        with tempfile.TemporaryDirectory() as tmp:
            local = os.path.join(tmp, os.path.basename(r["uri"]))
            got = fetch(r["url"], local)
            # The ingest step asserts this, so a mismatch here is the clearer
            # place to fail: it means the corpus moved, not that the code broke.
            if got != r["filesize"]:
                raise SystemExit(
                    f"size mismatch for {r['uri']}: CDN {got}, manifest {r['filesize']}"
                )
            subprocess.run(
                ["aws", "s3", "cp", "--only-show-errors", local, r["uri"]],
                check=True,
            )
        log(f"{r['uri']} ({got} bytes)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
