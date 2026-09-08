#!/usr/bin/env python3
"""Ingest one CELLxGENE dataset into an AnnData-Zarr store.

This is the ingestion step of a per-dataset store model: one store per dataset,
written independently, composed at query time through a catalog. Harmonization
belongs here, during ingestion, because a store is meant to be queryable
alongside other stores and not merely the Discover h5ad in another container.

This script owns the I/O: fetch the source, write the store, publish it, and
record a receipt. The data model lives in the `census2_store` package, which
decides which slots a store carries and applies harmonization. Its harmonization
is a stub today, so the output currently does match the source.
That is a temporary state, not the design.

Memory discipline: this script never calls `read_h5ad`. `census2_store.build`
returns a whole AnnData whose small slots are in memory and whose `X` and `raw.X` are dask
arrays of sparse row blocks, so peak memory tracks the block size rather than the
file size. The `h5py` handle therefore has to stay open until `write_zarr`.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import subprocess
import sys
import time

# The data model, installed in the runtime image from census2-store/ at the
# repository root. It lives there so that other components can use it; this
# script is the pipeline's I/O around it.
import census2_store
import h5py
import numpy as np
import zarr
from anndata.io import read_elem

RECIPE = "v1"


def log(msg: str) -> None:
    print(f"ingest: {msg}", file=sys.stderr, flush=True)


def store_uri(outdir: str, name: str) -> str:
    """Final location of one artifact. Local path or s3:// URI."""
    return f"{outdir.rstrip('/')}/{RECIPE}/{name}"


def _sha256(path: str) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ingest(src: str, dest: str, *, config=None, refdata=None) -> dict:
    """Ingest one h5ad into a Zarr store at `dest`, which may be remote.

    Returns the result from `census2_store.create`, which reports whether the
    dataset was admitted. A rejection writes no store and is not an error.
    """
    return census2_store.create(src, dest, config=config, refdata=refdata)


def check(src: str, dst: str) -> dict:
    """Compare the Zarr output against the h5ad source, field by field.

    `compared` lists what the comparison covered, so a check that silently
    covered nothing is visible in the receipt instead of reading as a pass.
    """
    import hashlib

    import pandas as pd

    result: dict = {"passed": True, "compared": []}

    def fail(name, detail):
        result["passed"] = False
        result.setdefault("failures", []).append(f"{name}: {detail}")

    # Walk one segment at a time. h5py and zarr both support that, while support
    # for a nested "raw/X" lookup differs between them.
    def has(node, path):
        for part in path.split("/"):
            if part not in node:
                return False
            node = node[part]
        return True

    def get(node, path):
        for part in path.split("/"):
            node = node[part]
        return node

    def digest(arr):
        h = hashlib.sha256()
        n = arr.shape[0]
        step = max(1, census2_store.BLOCK_BYTES // max(arr.dtype.itemsize, 1))
        for i in range(0, n, step):
            block = np.ascontiguousarray(arr[i : min(i + step, n)])
            # The round trip widens sparse index arrays from int32 to int64, so
            # hashing the raw bytes would report a difference that is not one.
            # Integers compare by value; floats stay strict, because a change of
            # precision in the data is worth reporting.
            if block.dtype.kind in "iu":
                block = block.astype(np.int64, copy=False)
            h.update(block.tobytes())
        return h.hexdigest()[:16]

    with h5py.File(src, "r") as f:
        # Route through as_store: a bare s3:// string makes zarr fall back to
        # fsspec, which needs s3fs. This is why an earlier HealthOmics run failed
        # after writing the store correctly.
        zroot = zarr.open_consolidated(census2_store.as_store(dst), mode="r")
        result["consolidated"] = True

        for key in ("obs", "var"):
            a, b = read_elem(f[key]), read_elem(zroot[key])
            result["compared"].append(key)
            if not a.index.equals(b.index):
                fail(key, "index differs")
            if list(a.columns) != list(b.columns):
                fail(key, "columns differ")
            for col in a.columns:
                if str(a[col].dtype) != str(b[col].dtype):
                    fail(f"{key}.{col}", f"dtype {a[col].dtype} != {b[col].dtype}")
                elif isinstance(a[col].dtype, pd.CategoricalDtype) and (
                    list(a[col].cat.categories) != list(b[col].cat.categories)
                ):
                    fail(f"{key}.{col}", "categories differ")

        # An absent matrix is only acceptable when neither side has it. Skipping
        # a present one would let the check pass having compared nothing.
        for path in ("X", "raw/X"):
            in_src, in_out = has(f, path), has(zroot, path)
            if not in_src and not in_out:
                continue
            if not in_out:
                fail(path, "missing from the store")
                continue
            if not in_src:
                fail(path, "in the store but not in the source")
                continue

            h5o, zo = get(f, path), get(zroot, path)
            if isinstance(h5o, h5py.Group):  # sparse
                if h5o["data"].shape != zo["data"].shape:
                    fail(path, f"nnz {h5o['data'].shape} != {zo['data'].shape}")
                    continue
                for part in ("data", "indices", "indptr"):
                    if digest(h5o[part]) != digest(zo[part]):
                        fail(
                            f"{path}/{part}",
                            f"values differ (h5ad {h5o[part].dtype}, zarr {zo[part].dtype})",
                        )
                result[f"{path}.nnz"] = int(h5o["data"].shape[0])
            else:  # dense
                if tuple(h5o.shape) != tuple(zo.shape):
                    fail(path, f"shape {h5o.shape} != {zo.shape}")
                    continue
                if digest(h5o) != digest(zo):
                    fail(path, "checksum differs")
                result[f"{path}.shape"] = list(h5o.shape)
            result["compared"].append(path)

    return result


def s3_cp(uri: str, dest: str) -> None:
    """Fetch the source h5ad. The CLI does parallel multipart with retries."""
    subprocess.run(["aws", "s3", "cp", "--only-show-errors", uri, dest], check=True)


def _obstore_url(path_or_uri: str) -> str:
    """obstore addresses everything by URL, so make a local path into one."""
    if census2_store.is_remote(path_or_uri):
        return path_or_uri
    return "file://" + os.path.abspath(path_or_uri)


def verify_written(uri: str) -> None:
    """Confirm a directly written store arrived and is readable.

    A local destination is a declared Nextflow output, so the workflow fails if
    it is missing. A direct write to object storage has no such check, and an
    earlier version of this pipeline discarded stores silently. So verify here.
    """
    if store_bytes(uri) == 0:
        raise SystemExit(f"store is empty after writing: {uri}")
    zarr.open_consolidated(census2_store.as_store(uri), mode="r")


def store_bytes(uri: str) -> int:
    """Total bytes of a written store, local or remote."""
    import obstore
    import obstore.store

    store = obstore.store.from_url(_obstore_url(uri))
    return sum(int(item["size"]) for page in obstore.list(store) for item in page)


def write_text(uri: str, text: str) -> str:
    """Write a small file beside the store, local or remote."""
    import obstore
    import obstore.store

    parent, name = uri.rsplit("/", 1)
    if not census2_store.is_remote(uri):
        os.makedirs(parent, exist_ok=True)
    obstore.put(obstore.store.from_url(_obstore_url(parent)), name, text.encode())
    return uri


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--record", required=True, help="one manifest line, as JSON")
    p.add_argument("--outdir", required=True, help="local path or s3:// prefix")
    p.add_argument(
        "--all-slots",
        action="store_true",
        help="override config: carry every slot, not just the Census set",
    )
    p.add_argument(
        "--config",
        default=None,
        help="ingestion configuration file, staged by the workflow (optional)",
    )
    p.add_argument(
        "--refdata",
        default=None,
        help="directory of reference data, staged by the workflow (optional)",
    )
    p.add_argument(
        "--check",
        action="store_true",
        help="verify the output against the source; a failure lands in the receipt "
        "and REPORT fails the run",
    )
    p.add_argument(
        "--stage-local",
        action="store_true",
        help="write the store and receipt into the working directory for the "
        "workflow to publish, instead of writing them to --outdir directly",
    )
    p.add_argument(
        "--keep-source", action="store_true", help="do not delete the staged h5ad"
    )
    args = p.parse_args()

    rec = json.loads(args.record)
    version_id = rec["dataset_version_id"]
    started = time.time()

    src = f"{version_id}.h5ad"
    log(f"fetch {rec['uri']} ({rec['filesize'] / 1024**2:.1f} MiB)")
    s3_cp(rec["uri"], src)

    actual = os.path.getsize(src)
    if actual != rec["filesize"]:
        raise SystemExit(
            f"size mismatch for {version_id}: got {actual}, manifest says {rec['filesize']}"
        )

    # Policy lives in the configuration. The CLI flag is an override for one key.
    config = census2_store.load_config(args.config)
    if args.all_slots:
        config["all_slots"] = True

    # A store must be traceable to the configuration that produced it, so record
    # the effective config and not just the file it came from.
    provenance = {
        "config_file": os.path.basename(args.config) if args.config else None,
        "config_sha256": _sha256(args.config) if args.config else None,
        "config_effective": config,
        "refdata": os.path.basename(args.refdata.rstrip("/")) if args.refdata else None,
    }
    if args.config or args.refdata:
        log(f"config={provenance['config_file']} refdata={provenance['refdata']}")

    # Where the store ends up, recorded in the receipt either way.
    store = store_uri(args.outdir, f"{version_id}.zarr")
    # Where this task writes it. When staging locally the workflow publishes it,
    # so the task must not reach outside its working directory.
    dest = f"{version_id}.zarr" if args.stage_local else store

    log(f"ingest -> {store} (slots={'all' if config['all_slots'] else 'census'})")
    result = ingest(src, dest, config=config, refdata=args.refdata)

    receipt = {
        "recipe": RECIPE,
        "dataset_id": rec["dataset_id"],
        "dataset_version_id": version_id,
        "source_uri": rec["uri"],
        "source_url": rec.get("url"),
        "source_filesize": actual,
        "cell_count": rec.get("cell_count"),
        "schema_version": rec.get("schema_version"),
        "admitted": result["admitted"],
        "reason": result["reason"],
        # Which code wrote this store. The recipe still gates semantics; these
        # are the evidence for what produced a given output.
        "census2_store_version": importlib.metadata.version("census2-store"),
        "anndata_version": importlib.metadata.version("anndata"),
        "zarr_version": importlib.metadata.version("zarr"),
        **provenance,
    }

    if result["admitted"]:
        for slot in result["slots"]:
            log(f"  wrote {slot}")
        receipt.update(
            {
                "slot_policy": "all" if config["all_slots"] else "census",
                "slots_written": result["slots"],
                "zarr_format": 3,
                "store": store,
                "output_bytes": store_bytes(dest),
                "published_by": "workflow" if args.stage_local else "task",
            }
        )
        if not args.stage_local:
            # Nothing tracks a direct write, so confirm the store really arrived
            # rather than trusting a silent success.
            verify_written(dest)
        if args.check:
            log("check round trip")
            receipt["check"] = check(src, dest)
            if not receipt["check"]["passed"]:
                for fmsg in receipt["check"].get("failures", []):
                    log(f"  FAIL {fmsg}")
    else:
        # Not admitted is a normal outcome. No store is written and we succeed.
        log(f"not admitted: {result['reason']}")

    if not args.keep_source:
        os.remove(src)

    receipt["elapsed_sec"] = round(time.time() - started, 1)

    # Always written: this is the declared Nextflow output. When staging locally
    # the workflow publishes it beside the store; otherwise write it to --outdir.
    local_receipt = f"{version_id}.json"
    text = json.dumps(receipt, indent=1)
    with open(local_receipt, "w") as fh:
        fh.write(text)
    if not args.stage_local:
        write_text(store_uri(args.outdir, local_receipt), text)

    if result["admitted"]:
        log(
            f"done {version_id} in {receipt['elapsed_sec']}s, "
            f"{receipt['output_bytes'] / 1024**2:.1f} MiB out"
        )
    else:
        log(f"done {version_id} in {receipt['elapsed_sec']}s, no store written")

    # A failed check is recorded, not raised. Nextflow publishes nothing from a
    # failed task, so exiting non-zero here would hide the receipt that names the
    # failure and would stop REPORT from running. REPORT fails the run instead,
    # which also reports every failure together rather than only the first.
    return 0


if __name__ == "__main__":
    sys.exit(main())
