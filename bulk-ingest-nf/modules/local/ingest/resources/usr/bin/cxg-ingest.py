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


# The slot groups that the check does not compare yet. A slot group that the
# writer wrote must be compared, or be named here, or the check fails. The
# receipt lists the groups that a store holds and the check skipped, in
# `not_compared`.
#   uns     the carried keys, such as schema_version and the organism
#   obsm, varm, obsp, varp, layers
#           the optional slots, which only `all_slots` carries
# Known gaps, not oversights. Shrink this list; do not grow it to stop a failure.
NOT_COMPARED = ("uns", "obsm", "varm", "obsp", "varp", "layers")


def slot_group(slot: str) -> str:
    """The group of a slot that the writer reports: "uns[2 keys]" is "uns",
    and "obsm/X_umap" is "obsm"."""
    return slot.split("/", 1)[0].split("[", 1)[0]


def obs_column_problem(
    got: pd.Series, want: pd.Series, col_type: census2_dataset_writer.ColType
) -> str | None:
    """What is wrong with the store's obs column `got`, or None.

    `want` holds the values that the column must hold. A categorical must have
    the sorted distinct values of `want` as its string categories.
    """
    categorical = isinstance(got.dtype, pd.CategoricalDtype)
    if col_type is census2_dataset_writer.ColType.CATEGORY:
        categories = sorted(set(want.astype(str)))
        if not categorical or list(got.cat.categories) != categories:
            return "not a categorical of the sorted distinct values"
    elif col_type is census2_dataset_writer.ColType.BOOL:
        if got.dtype != bool:
            return f"dtype {got.dtype}, not bool"
    elif categorical:
        return "categorical, not plain strings"
    as_type = bool if col_type is census2_dataset_writer.ColType.BOOL else str
    if not np.array_equal(want.astype(as_type), got.astype(as_type)):
        return "values differ"
    return None


STATS = ("raw_sum", "nnz", "raw_mean_nnz", "raw_variance_nnz", "n_measured_vars")


def row_stats(x: sp.csr_matrix) -> np.ndarray:
    """raw_sum, nnz, mean, and sample variance of the nonzero values of each row.

    These are the Census 1 definitions (its `get_obs_stats`), written here apart
    from the writer: duplicates merged and explicit zeros removed, and a variance
    of values shifted by the first value of the row, squared in float32 for float32
    counts and summed in float64. The writer must match them exactly.
    """
    x = x.copy()
    x.sum_duplicates()
    x.eliminate_zeros()
    nnz = np.diff(x.indptr)
    total = np.asarray(x.sum(axis=1, dtype=np.float64)).ravel()
    rows = np.repeat(np.arange(len(nnz)), nnz)
    data = x.data if x.data.dtype == np.float32 else x.data.astype(np.float64)
    shift = data - data[x.indptr[rows]]
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(nnz > 0, total / nnz, 0.0)
        s1 = np.bincount(rows, weights=shift, minlength=len(nnz))
        s2 = np.bincount(rows, weights=shift * shift, minlength=len(nnz))
        var = np.where(nnz > 1, (s2 - s1**2 / nnz) / (nnz - 1), 0.0)
    return np.column_stack([total, nnz, mean, var])


def check(src: str, dst: str, written: dict, summary: dict) -> dict:
    """Compare the Zarr output against the h5ad source, field by field.

    `compared` lists what the comparison covered, so a check that silently
    covered nothing is visible in the receipt instead of reading as a pass.
    The check also audits its coverage. Each slot group that the writer wrote
    (`written["slots"]`), and the summary record, must be in `compared` or in
    `NOT_COMPARED`. Otherwise the check fails. `not_compared` lists the groups
    that the store holds and the check skipped by name.

    `summary` is the summary record that the store must hold.

    The store holds the source cells that the cell filters keep, in source
    order. Its obs index must be the source index without the dropped cells,
    and the number of cells dropped must be the total of `dropped_cells`. Each
    comparison below reads the source at the kept cells. The check does not
    repeat the filters, and so it does not show that the right cells were
    dropped.

    The store's obs must hold exactly the columns of `OBS_COLUMNS`, in order,
    each of its type. `dataset_id` must hold the summary record's value, the
    `tissue_general` pair the roll-up of the source's `tissue_ontology_term_id`
    (from `census2_dataset_writer.tissues`), and each other column the source's
    values. The comparison is as strings, or as bool for a bool column. The
    statistics of the counts must equal those of the kept source counts, and
    `n_measured_vars` the number of kept source features.

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
        keep = a.index.isin(b.index)
        if not a.index[keep].equals(b.index):
            fail("obs", "index is not the source index without the dropped cells")
        if int((~keep).sum()) != sum(written["dropped_cells"].values()):
            fail(
                "obs",
                f"{int((~keep).sum())} cells dropped, "
                f"expected {sum(written['dropped_cells'].values())}",
            )
        a = a[keep]
        if list(b.columns) != [name for name, _ in schema]:
            fail("obs", f"columns {list(b.columns)} are not the schema's")
        else:
            tissue = a["tissue_ontology_term_id"].astype(str)
            general = census2_dataset_writer.tissues.tissue_general(tissue.unique())
            for name, col_type in schema:
                # The writer adds dataset_id, with the summary record's value.
                if name == census2_dataset_writer.DATASET_ID_COLUMN:
                    want = pd.Series(summary["dataset_id"], index=b.index)
                elif name == "tissue_general":
                    want = tissue.map({t: pair[1] for t, pair in general.items()})
                elif name == "tissue_general_ontology_term_id":
                    want = tissue.map({t: pair[0] for t, pair in general.items()})
                elif name in STATS:
                    # Compared with the source counts, block by block, below.
                    if b[name].dtype != col_type.value:
                        fail(
                            f"obs.{name}",
                            f"dtype {b[name].dtype}, not {col_type.value}",
                        )
                    continue
                else:
                    want = a[name]
                problem = obs_column_problem(b[name], want, col_type)
                if problem:
                    fail(f"obs.{name}", problem)

        var = read_elem(zroot["var"])
        result["compared"].append("var")
        if not var.index.is_unique:
            fail("var", "feature_ids are not unique")
        measured = census2_dataset_writer.MEASURED_COLUMN
        if list(var.columns) != [*census2_dataset_writer.genes.COLUMNS, measured]:
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
        n_kept = int(keep.sum())
        if (positions < 0).any():
            missing = list(ids[kept][positions < 0][:5])
            fail("X", f"{int((positions < 0).sum())} features not in var: {missing}")
        elif measured in var.columns and set(np.flatnonzero(var[measured])) != set(
            positions
        ):
            fail(f"var.{measured}", "does not mark exactly the kept source features")
        elif not isinstance(zroot["X"], zarr.Group):
            fail("X", "not a sparse matrix")
        else:
            store_x = sparse_dataset(zroot["X"])
            if tuple(store_x.shape) != (n_kept, len(var)):
                fail("X", f"shape {store_x.shape} != {(n_kept, len(var))}")
            else:
                step = census2_dataset_writer.row_chunk(h5o, n_rows)
                nnz = 0
                expected = []
                # The store keeps source order, so a running offset finds each block.
                start = 0
                for i in range(0, n_rows, step):
                    rows = slice(i, min(i + step, n_rows))
                    kept_rows = keep[rows]
                    stop = start + int(kept_rows.sum())
                    want = sp.csr_matrix(source_x[rows])[kept_rows][:, columns]
                    got = sp.csr_matrix(store_x[start:stop])
                    expected.append(row_stats(want))
                    start = stop
                    nnz += got.count_nonzero()
                    if got.dtype != want.dtype:
                        fail("X", f"dtype {got.dtype} != source {want.dtype}")
                        break
                    if (got[:, positions] != want).nnz or (
                        got.count_nonzero() != want.count_nonzero()
                    ):
                        fail(
                            "X",
                            f"values differ in source rows {rows.start}:{rows.stop}",
                        )
                        break
                result["X.nnz"] = int(nnz)
                result["X.source"] = counts
                # A failed block stops the loop early, so compare after a full pass.
                complete = sum(len(e) for e in expected) == len(b)
                if complete and set(STATS) <= set(b.columns):
                    want_stats = np.vstack(expected)
                    for j, name in enumerate(STATS[:4]):
                        if not np.array_equal(b[name], want_stats[:, j]):
                            fail(f"obs.{name}", "differs from the source counts")
                    if (b["n_measured_vars"] != len(columns)).any():
                        fail("obs.n_measured_vars", f"is not {len(columns)}")
                    result["compared"].append("obs statistics")
        result["compared"].append("X")

    required = {"summary"} | {slot_group(s) for s in written["slots"]}
    uncovered = sorted(required - set(result["compared"]) - set(NOT_COMPARED))
    if uncovered:
        fail("coverage", f"the store holds {uncovered}, and the check compared none")
    result["not_compared"] = sorted(required & set(NOT_COMPARED))
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
        "dropped_cells": result["dropped_cells"],
        **config_provenance,
    }
    if any(result["dropped_cells"].values()):
        log(f"  cell filters dropped {result['dropped_cells']}")

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
