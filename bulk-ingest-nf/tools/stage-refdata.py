#!/usr/bin/env python3
"""Stage the reference data that `--refdata` names: the pinned gene references.

A HealthOmics task has no public internet access, so the writer never fetches a
reference. Run this script once, on a host with internet access. Then give its
destination to the pipeline as `--refdata`.

The files and their SHA-256 digests come from `gene_references` in
`census2_dataset_writer.DEFAULTS`, so this script needs the package installed.
`aria2c` does the download: it streams each file to disk in parallel chunks,
retries, resumes a partial file, and checks the file against its digest. A
second run checks the files that are present and downloads only what is missing
or wrong. For an S3 destination, the files go to a local scratch directory
first, so the host needs disk space for them.

A failed run exits with an error. For a local destination, it leaves each
partial file under its final name, with a `.aria2` control file beside it, and
the next run resumes the file. The writer refuses a file whose digest differs,
so a partial file cannot pass as a reference.

    tools/stage-refdata.py refdata/                     # a local directory
    tools/stage-refdata.py s3://bucket/prefix/refdata/  # an S3 prefix, via the AWS CLI
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import census2_dataset_writer

# Parallel chunks of one file, and files at the same time. A server that refuses
# range requests gets one connection for each file. aria2c splits only a file
# that is larger than two chunks.
CONNECTIONS = 4
CHUNK = "16M"
FILES = 4


def aria2_input(references: dict) -> str:
    """One entry for each file: the URL, then its options, indented."""
    return "".join(
        f"{spec['url']}\n  out={spec['file']}\n  checksum=sha-256={spec['sha256']}\n"
        for spec in references.values()
    )


def download(references: dict, dest: str) -> None:
    if shutil.which("aria2c") is None:
        raise SystemExit(
            "aria2c is not installed: `brew install aria2` or `apt-get install aria2`"
        )
    subprocess.run(
        [
            "aria2c",
            "--input-file=-",
            f"--dir={dest}",
            f"--max-concurrent-downloads={FILES}",
            f"--split={CONNECTIONS}",
            f"--max-connection-per-server={CONNECTIONS}",
            f"--min-split-size={CHUNK}",
            # Resume a partial file. Check a complete one against its digest,
            # and download it again only when the check fails.
            "--continue=true",
            "--check-integrity=true",
            "--max-tries=10",
            "--retry-wait=5",
            # Write each file under its configured name, never under a renamed copy.
            "--auto-file-renaming=false",
            "--allow-overwrite=true",
            "--console-log-level=warn",
            "--summary-interval=0",
        ],
        input=aria2_input(references),
        text=True,
        check=True,
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("dest", help="local directory or s3:// prefix")
    args = p.parse_args()

    references = census2_dataset_writer.DEFAULTS["gene_references"]
    for taxon, spec in references.items():
        print(f"stage {taxon} {spec['version']}: {spec['url']}", file=sys.stderr)
    remote = args.dest.startswith("s3://")
    with tempfile.TemporaryDirectory() as scratch:
        local = scratch if remote else args.dest
        os.makedirs(local, exist_ok=True)
        download(references, local)
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
