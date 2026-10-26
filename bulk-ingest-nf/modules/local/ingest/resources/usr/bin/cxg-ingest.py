#!/usr/bin/env python3
"""Ingest one CELLxGENE dataset into an AnnData-Zarr store.

This is the ingestion step of a per-dataset store model: one store per dataset,
written independently, composed at query time through a catalog. Harmonization
belongs here, during ingestion, because a store is meant to be queryable
alongside other stores and not merely the Discover h5ad in another container.

This script owns the I/O: fetch the source, write the store, publish it, and
record a receipt. The data model lives in the `census2_dataset_writer` package,
which decides which slots a store carries and applies harmonization: it projects
the counts onto the pinned gene reference staged in `--refdata`.

Memory discipline: this script never calls `read_h5ad`. It hands the source path
to `census2_dataset_writer.create`, which keeps the small slots in memory and
leaves `X` and `raw.X` as dask arrays of row blocks. Peak memory therefore
tracks the block size and not the file size. `create` holds the `h5py` handle
open until the write finishes.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

# The data model, installed in the runtime image from census2-dataset-writer/ at
# the repository root.
import census2_dataset_writer
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.io import read_elem, sparse_dataset


def log(msg: str) -> None:
    print(f"ingest: {msg}", file=sys.stderr, flush=True)


def store_uri(outdir: str, name: str) -> str:
    """Final location of one artifact. Local path or s3:// URI."""
    return f"{outdir.rstrip('/')}/{name}"


def _sha256(path: str) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ingest(
    src: str, dest: str, rec: dict, *, config=None, refdata=None, threads=None
) -> dict:
    """Ingest one h5ad into a Zarr store at `dest`, which may be remote.

    Returns the result from `census2_dataset_writer.create`, which reports
    whether the dataset was admitted. A rejection writes no store and is not an
    error. The identifiers come from the manifest record `rec`, because the h5ad
    does not hold them, and they go into the store's summary record.
    """
    return census2_dataset_writer.create(
        src,
        dest,
        dataset_id=rec["dataset_id"],
        dataset_version_id=rec["dataset_version_id"],
        config=config,
        refdata=refdata,
        threads=threads,
    )


def obs_column_problem(got: pd.Series, want: pd.Series, kind: str) -> str | None:
    """What is wrong with the store's obs column `got`, or None.

    `want` holds the values that the column must hold. A categorical must have
    the sorted distinct values of `want` as its string categories.
    """
    categorical = isinstance(got.dtype, pd.CategoricalDtype)
    if kind == census2_dataset_writer.CATEGORY:
        categories = sorted(set(want.astype(str)))
        if not categorical or list(got.cat.categories) != categories:
            return "not a categorical of the sorted distinct values"
    elif kind == census2_dataset_writer.BOOL:
        if got.dtype != bool:
            return f"dtype {got.dtype}, not bool"
    elif categorical:
        return "categorical, not plain strings"
    as_type = bool if kind == census2_dataset_writer.BOOL else str
    if not np.array_equal(want.astype(as_type), got.astype(as_type)):
        return "values differ"
    return None


def check(src: str, dst: str, written: dict, summary: dict) -> dict:
    """Compare the Zarr output against the h5ad source, field by field.

    `compared` lists what the comparison covered, so a check that silently
    covered nothing is visible in the receipt instead of reading as a pass.

    `summary` is the summary record that the store must hold.

    The store's obs must hold exactly the columns of `OBS_COLUMNS`, in order,
    each of its kind, and the source's index. `dataset_id` must hold the summary
    record's value, and each other column the source's values, compared as
    strings, or as bool for a bool column.

    `written` is the result of `census2_dataset_writer.create`.

    The store's `var` is the gene reference, so it is not compared with the
    source `var`. The store's `X` must hold the source counts (`raw/X`, else
    `X`) of every feature except the `dropped_var_features`, each at the
    position of its `feature_id` in the store's `var`. The check reads the store
    at those positions and compares the values with the source, one block of
    rows at a time. The writer moves each source column to its position, and
    the check reads each position back. Thus the check does not repeat the
    writer's code. Equal nonzero counts per block show that no value is outside
    those positions.
    """

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

    with h5py.File(src, "r") as f:
        # Route through as_store: a bare s3:// string makes zarr fall back to
        # fsspec, which needs s3fs.
        zroot = zarr.open_consolidated(census2_dataset_writer.as_store(dst), mode="r")
        result["consolidated"] = True

        result["compared"].append("summary")
        try:
            stored = census2_dataset_writer.read_summary(dst)
        except census2_dataset_writer.SummaryRecordError as e:
            fail("summary", str(e))
        else:
            if stored != summary:
                fail("summary", f"the store holds {stored}, expected {summary}")

        a, b = read_elem(f["obs"]), read_elem(zroot["obs"])
        result["compared"].append("obs")
        schema = census2_dataset_writer.OBS_COLUMNS
        if not a.index.equals(b.index):
            fail("obs", "index differs")
        if list(b.columns) != [name for name, _ in schema]:
            fail("obs", f"columns {list(b.columns)} are not the schema's")
        else:
            for name, kind in schema:
                # The writer adds dataset_id, with the summary record's value.
                if name == census2_dataset_writer.DATASET_ID_COLUMN:
                    want = pd.Series(summary["dataset_id"], index=b.index)
                else:
                    want = a[name]
                problem = obs_column_problem(b[name], want, kind)
                if problem:
                    fail(f"obs.{name}", problem)

        var = read_elem(zroot["var"])
        result["compared"].append("var")
        if not var.index.is_unique:
            fail("var", "feature_ids are not unique")
        if list(var.columns) != census2_dataset_writer.genes.COLUMNS:
            fail("var", f"columns {list(var.columns)} are not the reference columns")
        # anndata writes an absent raw as a null-encoded element, not a missing key.
        if has(zroot, "raw") and zroot["raw"].attrs.get("encoding-type") == "raw":
            fail("raw", "the store holds raw; the counts belong in X")

        counts = "raw/X" if has(f, "raw/X") else "X"
        source_var = read_elem(
            f["raw/var"] if counts == "raw/X" and has(f, "raw/var") else f["var"]
        )
        ids = source_var.index.astype(str)
        kept = ~ids.isin(set(written["dropped_var_features"]))
        columns = np.flatnonzero(kept)
        positions = var.index.get_indexer(ids[kept])
        h5o = get(f, counts)
        source_x = sparse_dataset(h5o) if isinstance(h5o, h5py.Group) else h5o
        n_rows = int(source_x.shape[0])
        if (positions < 0).any():
            missing = list(ids[kept][positions < 0][:5])
            fail("X", f"{int((positions < 0).sum())} features not in var: {missing}")
        elif not isinstance(zroot["X"], zarr.Group):
            fail("X", "not a sparse matrix")
        else:
            store_x = sparse_dataset(zroot["X"])
            if tuple(store_x.shape) != (n_rows, len(var)):
                fail("X", f"shape {store_x.shape} != {(n_rows, len(var))}")
            else:
                step = census2_dataset_writer.row_chunk(h5o, n_rows)
                nnz = 0
                for i in range(0, n_rows, step):
                    rows = slice(i, min(i + step, n_rows))
                    want = sp.csr_matrix(source_x[rows])[:, columns]
                    got = sp.csr_matrix(store_x[rows])
                    nnz += got.count_nonzero()
                    if got.dtype != want.dtype:
                        fail("X", f"dtype {got.dtype} != source {want.dtype}")
                        break
                    if (got[:, positions] != want).nnz or (
                        got.count_nonzero() != want.count_nonzero()
                    ):
                        fail("X", f"values differ in rows {rows.start}:{rows.stop}")
                        break
                result["X.nnz"] = int(nnz)
                result["X.source"] = counts
        result["compared"].append("X")

    return result


def s3_cp(uri: str, dest: str) -> None:
    """Fetch the source h5ad. The CLI does parallel multipart with retries."""
    subprocess.run(["aws", "s3", "cp", "--only-show-errors", uri, dest], check=True)


def _obstore_url(path_or_uri: str) -> str:
    """obstore addresses everything by URL, so make a local path into one."""
    if census2_dataset_writer.is_remote(path_or_uri):
        return path_or_uri
    return "file://" + os.path.abspath(path_or_uri)


def verify_written(uri: str) -> None:
    """Confirm a directly written store arrived, is readable, and has its summary record.

    A local destination is a declared Nextflow output, so the workflow fails if
    it is missing. A direct write to object storage has no such check, so verify
    it here.
    """
    if store_bytes(uri) == 0:
        raise SystemExit(f"store is empty after writing: {uri}")
    zarr.open_consolidated(census2_dataset_writer.as_store(uri), mode="r")
    if census2_dataset_writer.read_summary(uri) is None:
        raise SystemExit(f"store has no summary record after writing: {uri}")


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
    if not census2_dataset_writer.is_remote(uri):
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
        required=True,
        help="directory of reference data, staged by the workflow; holds the "
        "gene reference files",
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
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="CPUs this task was allocated. Bounds the dask thread pool that "
        "writes the store. Without it dask counts the host's CPUs.",
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
    config = census2_dataset_writer.load_config(args.config)
    if args.all_slots:
        config["all_slots"] = True

    # A store must be traceable to the configuration that produced it, so record
    # the effective config and not just the file it came from.
    config_provenance = {
        "config_file": os.path.basename(args.config) if args.config else None,
        "config_sha256": _sha256(args.config) if args.config else None,
        "config_effective": config,
        "refdata": os.path.basename(args.refdata.rstrip("/")),
    }
    log(
        f"config={config_provenance['config_file']} "
        f"refdata={config_provenance['refdata']}"
    )

    # Where the store ends up, recorded in the receipt either way.
    store = store_uri(args.outdir, f"{version_id}.zarr")
    # Where this task writes it. When staging locally the workflow publishes it,
    # so the task must not reach outside its working directory.
    dest = f"{version_id}.zarr" if args.stage_local else store

    log(f"ingest -> {store} (slots={'all' if config['all_slots'] else 'census'})")
    result = ingest(
        src, dest, rec, config=config, refdata=args.refdata, threads=args.threads
    )

    # The shape of a receipt is census2_dataset_writer.Receipt, so this script
    # and ingestion-report.py agree on it without importing each other. The
    # store-level facts come from the package, so they cannot disagree with the
    # code that ran. For an admitted dataset the receipt also holds a copy of the
    # store's summary record, so that the report needs to open no store.
    receipt: census2_dataset_writer.Receipt = {
        **census2_dataset_writer.provenance(config),
        "dataset_id": rec["dataset_id"],
        "dataset_version_id": version_id,
        "source_uri": rec["uri"],
        "source_url": rec.get("url"),
        "source_filesize": actual,
        "cell_count": rec.get("cell_count"),
        "schema_version": rec.get("schema_version"),
        "admitted": result["admitted"],
        "reason": result["reason"],
        **config_provenance,
    }

    if result["admitted"]:
        for slot in result["slots"]:
            log(f"  wrote {slot}")
        if result["dropped_obs_columns"]:
            dropped = result["dropped_obs_columns"]
            log(f"  dropped {len(dropped)} obs columns: {', '.join(dropped)}")
        if result["dropped_var_features"]:
            dropped = result["dropped_var_features"]
            log(f"  dropped {len(dropped)} features that are not organism genes")
        log(f"  genes projected onto {result['gene_reference_version']}")
        summary = census2_dataset_writer.summary_record(result)
        receipt.update(
            {
                "slots_written": result["slots"],
                "dropped_obs_columns": result["dropped_obs_columns"],
                "dropped_var_features": result["dropped_var_features"],
                "summary": summary,
                "store": store,
                "output_bytes": store_bytes(dest),
                "published_by": "workflow" if args.stage_local else "task",
            }
        )
        if not args.stage_local:
            verify_written(dest)
        if args.check:
            log("check round trip")
            receipt["check"] = check(src, dest, result, summary)
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
    # failed task. A non-zero exit here would therefore hide the receipt that
    # names the failure, and it would stop REPORT. REPORT fails the run instead,
    # and it reports every failure together rather than only the first.
    return 0


if __name__ == "__main__":
    sys.exit(main())
