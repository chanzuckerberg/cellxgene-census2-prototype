#!/usr/bin/env python3
"""Consolidate one group of per-dataset stores into one consolidated store, with a receipt.

This script owns the I/O around `census2_dataset_writer.consolidation.consolidate`: where
the consolidated store goes, the check of the written store, and the receipt. The
package decides what a consolidated store holds. The group is one line of the plan
that `tools/consolidation-plan.py` writes.

The store is `<outdir>/<group>.zarr`, and the receipt is `<outdir>/<group>.json`.
With `--check`, the task compares the consolidated store with each member, value for
value. A failed check is recorded in the receipt, the task exits 0, and REPORT
fails the run; see bulk-ingest-nf/AGENTS.md for why.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

import census2_dataset_writer
import zarr
from census2_dataset_writer import consolidation


def log(msg: str) -> None:
    print(f"consolidate: {msg}", file=sys.stderr, flush=True)


def _obstore_url(path_or_uri: str) -> str:
    if census2_dataset_writer.is_remote(path_or_uri):
        return path_or_uri
    return "file://" + os.path.abspath(path_or_uri)


def store_bytes(uri: str) -> int:
    """Total bytes of a written store, local or remote."""
    import obstore
    import obstore.store

    store = obstore.store.from_url(_obstore_url(uri))
    return sum(int(item["size"]) for page in obstore.list(store) for item in page)


def write_text(uri: str, text: str) -> None:
    import obstore
    import obstore.store

    parent, name = uri.rsplit("/", 1)
    if not census2_dataset_writer.is_remote(uri):
        os.makedirs(parent, exist_ok=True)
    obstore.put(obstore.store.from_url(_obstore_url(parent)), name, text.encode())


def peak_rss_mib() -> float:
    """The peak resident memory of this process; Linux reports KiB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(peak / (2**20 if sys.platform == "darwin" else 2**10), 1)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--group", required=True, help="one plan line, as JSON")
    p.add_argument("--outdir", required=True, help="local path or s3:// prefix")
    p.add_argument("--threads", type=int, default=None, help="CPUs this task has")
    p.add_argument("--check", action="store_true", help="compare with each member")
    p.add_argument(
        "--stage-local",
        action="store_true",
        help="write into the working directory for the workflow to publish",
    )
    args = p.parse_args()

    group = json.loads(args.group)
    name = group["group"]
    uris = [m["store_uri"] for m in group["members"]]
    store = f"{args.outdir.rstrip('/')}/{name}.zarr"
    dest = f"{name}.zarr" if args.stage_local else store
    started = time.time()
    log(f"{len(uris)} members, {group['n_obs']:,} cells -> {store}")
    record = consolidation.consolidate(uris, dest, threads=args.threads)
    consolidate_sec = round(time.time() - started, 1)

    # A remote write is not a declared output, so confirm that it is complete.
    zarr.open_consolidated(census2_dataset_writer.as_store(dest), mode="r")
    if consolidation.read_consolidation(dest) != record:
        raise SystemExit(f"the store has no matching consolidation record: {dest}")

    receipt: dict = {
        **census2_dataset_writer.provenance(),
        "group": name,
        "store": store,
        "members": [m["dataset_version_id"] for m in record["members"]],
        "n_members": len(record["members"]),
        "n_obs": record["members"][-1]["row_stop"],
        "nnz": sum(m["nnz"] for m in record["members"]),
        "source_bytes": sum(int(m.get("bytes") or 0) for m in group["members"]),
        "output_bytes": store_bytes(dest),
        "consolidate_sec": consolidate_sec,
        "published_by": "workflow" if args.stage_local else "task",
    }
    if args.check:
        log("check each member")
        failures = consolidation.verify(dest, uris, progress=log)
        receipt["check"] = {"passed": not failures, "failures": failures}
        for failure in failures:
            log(f"  FAIL {failure}")
    receipt["elapsed_sec"] = round(time.time() - started, 1)
    receipt["peak_rss_mib"] = peak_rss_mib()

    text = json.dumps(receipt, indent=1)
    with open(f"{name}.json", "w") as fh:
        fh.write(text)
    if not args.stage_local:
        write_text(f"{args.outdir.rstrip('/')}/{name}.json", text)
    log(
        f"done {name} in {receipt['elapsed_sec']}s: "
        f"{receipt['source_bytes'] / 2**20:.1f} MiB in, "
        f"{receipt['output_bytes'] / 2**20:.1f} MiB out, peak {receipt['peak_rss_mib']} MiB"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
