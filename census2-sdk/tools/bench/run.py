"""Run one benchmark query in this process and print one JSON line.

Each run is a new process, so no connection or metadata cache carries over from
an earlier run. The shapes are:

- `single`: `get_anndata` for one dataset, the first row or `--dataset-id`. It
  opens one store at any catalog size.
- `cell_type`: `get_anndata` for one cell type across the catalog. The cell type
  is the one whose share of the source datasets is nearest to 15 percent. The
  copies are not part of that share.
- `filter`: `get_anndata` for the filter that `--filter` gives as JSON.
- `builder`: the query of `bench.filters.QUERIES` that `--builder` names, through the
  typed filter builder (`census.homo_sapiens`), by the query's mode: `get_anndata`,
  `get_obs`, `count`, or `count_by`. The builder path always opens a store with the
  `lean` open; `--reader` sets the open of the other paths.
- `open`: open every store of the catalog and read nothing more.
- `corpus_obs`: open every store and read its whole `obs`.
- `corpus_x`: open every store and read all of `X`, one shard at a time, and keep
  nothing. This measures the time of a full read of `X` without a result in memory.

`--workers` is the largest number of stores read at the same time. The record
holds the elapsed time of the catalog read and of the query, and the CPU time of
the query (user and system). It also holds the number of datasets selected, of
stores opened, and of cells, the reads and the payload, and the peak memory.

`--reader` selects how a store opens. `sdk` is `census2_sdk.catalog.open_store`:
a new obstore client for each store, and no Zarr format, so zarr looks for the
metadata of both formats. `lean` gives `zarr_format=3` and shares one obstore
client for each bucket. Both count the reads, and the SDK queries use the same
open. `--processes` splits the stores of the corpus shapes over that many
processes, with `--workers` threads in all. The processes start, and import
their modules, before the clock starts.

A catalog that names consolidated stores is read through them: a query opens each
consolidated store once for its selected datasets, and the corpus shapes open each
consolidated store in place of its datasets. `--per-dataset` reads the per-dataset
stores instead.

`--max-rss-gb` is a memory limit. When the peak memory passes it, the process
prints a record with the status `memory_limit` and stops with exit code 3.

    python -m bench.run catalog-500.parquet --shape cell_type --workers 32
"""

import argparse
import importlib.metadata
import json
import math
import multiprocessing
import os
import resource
import sys
import threading
import time
from collections.abc import AsyncIterator, Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, cast

import anndata
import pandas as pd
import zarr
import zarr.storage
from anndata.io import read_elem
from zarr.abc.store import ByteRequest, Store
from zarr.core.buffer import Buffer, BufferPrototype

import census2_sdk
from census2_sdk import catalog as _catalog

SHAPES = ("single", "cell_type", "filter", "builder", "open", "corpus_obs", "corpus_x")
READERS = ("sdk", "lean")
ORGANISM = "homo_sapiens"
MEMORY_LIMIT_EXIT = 3


class Counter:
    """The reads and the data of every counting store in this process."""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.requests = 0
        self.nbytes = 0

    def add(
        self, buffers: Iterable[Buffer | None], requests: int | None = None
    ) -> None:
        """Count the data of `buffers`, and one request each unless `requests` is given."""
        sizes = [len(b) if b is not None else 0 for b in buffers]
        with self.lock:
            self.requests += len(sizes) if requests is None else requests
            self.nbytes += sum(sizes)


COUNTER = Counter()


class CountingStore(zarr.storage.WrapperStore[Store]):
    """A store that counts each read and its data in `COUNTER`.

    The count is at the zarr store interface. A sharded array reads parts of a
    shard through `get_ranges`. The wrapped store fetches each group of ranges
    together, so one group counts as one request. The payload is the data that
    zarr keeps. A group fetch also reads the gaps between its ranges, and the
    payload does not include those gaps.
    """

    async def get(
        self,
        key: str,
        prototype: BufferPrototype,
        byte_range: ByteRequest | None = None,
    ) -> Buffer | None:
        buffer = await self._store.get(key, prototype, byte_range)
        COUNTER.add([buffer])
        return buffer

    async def get_partial_values(
        self,
        prototype: BufferPrototype,
        key_ranges: Iterable[tuple[str, ByteRequest | None]],
    ) -> list[Buffer | None]:
        buffers = await self._store.get_partial_values(prototype, key_ranges)
        COUNTER.add(buffers)
        return buffers

    async def get_ranges(
        self,
        key: str,
        byte_ranges: Sequence[ByteRequest | None],
        *,
        prototype: BufferPrototype,
        max_concurrency: int | None = None,
        max_gap_bytes: int | None = None,
        max_coalesced_bytes: int | None = None,
    ) -> AsyncIterator[Sequence[tuple[int, Buffer | None]]]:
        groups = super().get_ranges(
            key,
            byte_ranges,
            prototype=prototype,
            max_concurrency=max_concurrency,
            max_gap_bytes=max_gap_bytes,
            max_coalesced_bytes=max_coalesced_bytes,
        )
        async for group in groups:
            COUNTER.add((buffer for _, buffer in group), requests=1)
            yield group


# How `open_counted` opens a store; one of READERS. `main` sets it.
READER = "sdk"
# The lean reader's counting store over one obstore client, for each bucket.
_BUCKETS: dict[str, CountingStore] = {}
_BUCKETS_LOCK = threading.Lock()


def open_counted(uri: str | os.PathLike[str], *, lean: bool = False) -> zarr.Group:
    """Open a store through a counting store, by the open of `READER`.

    `sdk` opens as `census2_sdk.catalog.open_store` does. `lean` opens with
    `zarr_format=3`, and through one obstore client for each bucket. A caller that
    asks for the lean open, as the builder path does, gets it with either reader.
    """
    if lean or READER == "lean":
        return open_lean(uri)
    if _catalog.is_remote(uri):
        import obstore.store

        base: Store = zarr.storage.ObjectStore(
            obstore.store.from_url(str(uri)), read_only=True
        )
    else:
        base = zarr.storage.LocalStore(str(uri), read_only=True)
    return zarr.open_group(CountingStore(base), mode="r")


def open_lean(uri: str | os.PathLike[str]) -> zarr.Group:
    """Open a store with `zarr_format=3`, through the obstore client of its bucket.

    The open reads one object, the root `zarr.json` with the consolidated metadata,
    and reuses the connections of the client.
    """
    if not _catalog.is_remote(uri):
        local = CountingStore(zarr.storage.LocalStore(str(uri), read_only=True))
        return zarr.open_group(local, mode="r", zarr_format=3)
    scheme, rest = str(uri).split("://", 1)
    bucket, _, path = rest.partition("/")
    with _BUCKETS_LOCK:
        store = _BUCKETS.get(bucket)
        if store is None:
            import obstore.store

            client = obstore.store.from_url(f"{scheme}://{bucket}")
            store = CountingStore(zarr.storage.ObjectStore(client, read_only=True))
            _BUCKETS[bucket] = store
    return zarr.open_group(store, path=path.rstrip("/"), mode="r", zarr_format=3)


def read_obs(uri: str) -> int:
    """Open a store, read its whole `obs`, and return the number of cells."""
    return len(cast("pd.DataFrame", read_elem(open_counted(uri)["obs"])))


def read_x(uri: str) -> int:
    """Open a store and read all of `X`, one shard at a time. Return the cells."""
    group = open_counted(uri)
    x = group["X"]
    assert isinstance(x, zarr.Group), "X must be an encoded sparse matrix"
    rows = 0
    for name in ("indptr", "indices", "data"):
        array = x[name]
        assert isinstance(array, zarr.Array)
        step = (array.shards or array.chunks)[0]
        for start in range(0, array.shape[0], step):
            array[start : start + step]
        if name == "indptr":
            rows = array.shape[0] - 1
    return rows


def open_only(uri: str) -> int:
    """Open a store and read nothing more."""
    open_counted(uri)
    return 0


CORPUS_READS: dict[str, Callable[[str], int]] = {
    "open": open_only,
    "corpus_obs": read_obs,
    "corpus_x": read_x,
}


def cpu_seconds() -> tuple[float, float]:
    """The user and system CPU time of this process so far."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime, usage.ru_stime


def read_corpus(shape: str, uris: list[str], threads: int) -> dict[str, Any]:
    """Read the stores `uris` by `shape` with `threads` threads, in this process.

    Return the cells, and the reads, the payload, and the CPU time of this call.
    """
    requests, nbytes = COUNTER.requests, COUNTER.nbytes
    user, system = cpu_seconds()
    with ThreadPoolExecutor(max_workers=threads) as pool:
        cells = sum(pool.map(CORPUS_READS[shape], uris))
    end_user, end_system = cpu_seconds()
    return {
        "cells": cells,
        "requests": COUNTER.requests - requests,
        "nbytes": COUNTER.nbytes - nbytes,
        "cpu_user": end_user - user,
        "cpu_sys": end_system - system,
    }


def _start_child(reader: str, ready: Any) -> None:
    global READER
    READER = reader
    ready.wait()


def _child_ready(_: int) -> None:
    return None


def _child_read(task: tuple[str, list[str], int]) -> dict[str, Any]:
    return read_corpus(*task)


def read_corpus_processes(
    shape: str, uris: list[str], workers: int, processes: int
) -> tuple[dict[str, Any], float, float]:
    """`read_corpus` over `processes` processes, `workers` threads in all.

    Return the sums of the parts, and the times at which the read started and
    ended. The processes start, and wait for one another, before the read starts.
    """
    context = multiprocessing.get_context("spawn")
    ready = context.Barrier(processes)
    threads = max(1, math.ceil(workers / processes))
    with context.Pool(processes, _start_child, (READER, ready)) as pool:
        pool.map(_child_ready, range(processes), chunksize=1)
        tasks = [(shape, uris[k::processes], threads) for k in range(processes)]
        start = time.monotonic()
        parts = pool.map(_child_read, tasks, chunksize=1)
        end = time.monotonic()
    total = {key: sum(part[key] for part in parts) for key in parts[0]}
    return total, start, end


def source_rows(rows: pd.DataFrame) -> pd.DataFrame:
    """The rows of the source stores, without the copies that `catalogs` adds."""
    return rows[rows["copy"] == 0] if "copy" in rows.columns else rows


# `catalog.CONSOLIDATED[0]`, by name, so that the SDK of an earlier commit runs this
# file as well, for a baseline.
CONSOLIDATED_URI = "consolidated_store_uri"


def serving_uris(rows: pd.DataFrame, per_dataset: bool = False) -> list[str]:
    """The stores that a read of `rows` opens, each once, in catalog order.

    A row that names a consolidated store is read there, unless `per_dataset` is set.
    """
    uris = [str(u) for u in rows["store_uri"]]
    if not per_dataset and CONSOLIDATED_URI in rows.columns:
        consolidated = rows[CONSOLIDATED_URI]
        uris = [
            str(c) if pd.notna(c) else u
            for c, u in zip(consolidated, uris, strict=True)
        ]
    return list(dict.fromkeys(uris))


def pick_cell_type(rows: pd.DataFrame, share: float = 0.15) -> str:
    """The cell type whose share of the datasets is nearest to `share`."""
    counts: dict[str, int] = {}
    for values in rows["cell_type"]:
        for value in set(values):
            counts[str(value)] = counts.get(str(value), 0) + 1
    return min(counts, key=lambda c: (abs(counts[c] / len(rows) - share), c))


def peak_rss_mb() -> float:
    """The peak resident memory of this process. macOS reports bytes, Linux KiB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(peak / (2**20 if sys.platform == "darwin" else 2**10), 1)


PUBLISH_LOCK = threading.Lock()
PUBLISHED = threading.Event()


def publish(record: dict[str, Any], exit_code: int | None = None) -> bool:
    """Print the one terminal record of this process. Return False if a record exists.

    The first caller prints, and later callers print nothing. With `exit_code`,
    the process stops after the print, while it still holds the lock.
    """
    with PUBLISH_LOCK:
        if PUBLISHED.is_set():
            return False
        PUBLISHED.set()
        print(json.dumps(record), flush=True)
        if exit_code is not None:
            os._exit(exit_code)
    return True


def guard_memory(
    limit_mb: float, record: dict[str, Any], interval: float = 0.5
) -> None:
    """Stop this process when its peak memory passes `limit_mb`.

    A thread checks the peak memory at each `interval`. Above the limit, it prints
    `record` with the status `memory_limit` and stops the process at once.
    """

    def watch() -> None:
        while True:
            peak = peak_rss_mb()
            if peak > limit_mb:
                over = {**record, "status": "memory_limit", "peak_rss_mb": peak}
                publish({**over, "limit_mb": limit_mb}, MEMORY_LIMIT_EXIT)
                return
            time.sleep(interval)

    threading.Thread(target=watch, daemon=True).start()


def final_record(
    record: dict[str, Any], limit_mb: float | None
) -> tuple[dict[str, Any], int | None]:
    """The record of a finished run, and its exit code when the run passed the limit.

    The watchdog checks at intervals, so a run can pass the limit and finish between
    two checks. This check of the final peak catches that run.
    """
    if limit_mb is None or record["peak_rss_mb"] <= limit_mb:
        return record, None
    return {**record, "status": "memory_limit", "limit_mb": limit_mb}, MEMORY_LIMIT_EXIT


def run(
    catalog: str,
    shape: str,
    workers: int,
    cell_type: str | None,
    *,
    dataset_id: str | None = None,
    obs_filter: dict[str, Any] | None = None,
    per_dataset: bool = False,
    processes: int = 1,
    builder: str | None = None,
) -> dict[str, Any]:
    t0 = time.monotonic()
    # Without the keyword, so that the SDK of an earlier commit runs this as well.
    census = (
        census2_sdk.open_census([catalog], consolidated=False)
        if per_dataset
        else census2_sdk.open_census([catalog])
    )
    rows = census.catalog[census.catalog["organism"] == ORGANISM]
    t1 = time.monotonic()
    catalog_seconds = round(t1 - t0, 3)
    record: dict[str, Any] = {
        "catalog": catalog,
        "catalog_stores": len(serving_uris(rows, per_dataset)),
        "catalog_datasets": len(rows),
        "shape": shape,
        "workers": workers,
        "per_dataset": per_dataset,
        "reader": READER,
        "processes": processes,
    }
    cpu = {"cpu_user": 0.0, "cpu_sys": 0.0}

    if shape == "builder":
        if builder is None:
            raise ValueError("the builder shape needs --builder")
        from bench.filters import QUERIES, run_query

        query = QUERIES[builder]
        organism = census.homo_sapiens
        record["builder"] = builder
        record["mode"] = query.mode
        record["filter"] = repr(query.build(organism.obs))
        original = _catalog.open_store
        census2_sdk.__dict__["open_store"] = open_counted
        _catalog.open_store = open_counted
        user, system = cpu_seconds()
        try:
            cells = run_query(organism, query, workers)
        finally:
            census2_sdk.__dict__["open_store"] = original
            _catalog.open_store = original
        t2 = time.monotonic()
        end_user, end_system = cpu_seconds()
        cpu = {"cpu_user": end_user - user, "cpu_sys": end_system - system}
        stats = organism.last_query
        record["datasets"] = stats.get("candidates", 0)
        record["stores_opened"] = stats.get("stores_opened", 0)
        record["reads"] = stats.get("reads", 0)
        record["cells"] = cells
    elif shape in ("single", "cell_type", "filter"):
        if shape == "single":
            obs_filter = {"dataset_id": dataset_id or str(rows["dataset_id"].iloc[0])}
        elif shape == "cell_type":
            cell_type = cell_type or pick_cell_type(source_rows(rows))
            obs_filter = {"cell_type": cell_type}
            record["cell_type"] = cell_type
        elif obs_filter is None:
            raise ValueError("the filter shape needs a filter")
        record["filter"] = obs_filter
        # The SDK imports open_store by name, so patch both names, and restore them.
        original = _catalog.open_store
        census2_sdk.__dict__["open_store"] = open_counted
        _catalog.open_store = open_counted
        user, system = cpu_seconds()
        try:
            adata = census.get_anndata(ORGANISM, obs_filter, max_workers=workers)
        finally:
            census2_sdk.__dict__["open_store"] = original
            _catalog.open_store = original
        t2 = time.monotonic()
        end_user, end_system = cpu_seconds()
        cpu = {"cpu_user": end_user - user, "cpu_sys": end_system - system}
        selected = _catalog.exclude(
            census.catalog, ORGANISM, _catalog.canonicalize(obs_filter)
        )
        record["datasets"] = len(selected)
        record["stores_opened"] = len(serving_uris(selected, per_dataset))
        record["cells"] = int(adata.n_obs) if isinstance(adata, anndata.AnnData) else 0
    elif processes > 1:
        uris = serving_uris(rows, per_dataset)
        total, start, t2 = read_corpus_processes(shape, uris, workers, processes)
        record["process_start_seconds"] = round(start - t1, 3)
        t1 = start
        COUNTER.requests += total["requests"]
        COUNTER.nbytes += total["nbytes"]
        cpu = {k: total[k] for k in cpu}
        record["datasets"] = len(rows)
        record["stores_opened"] = len(uris)
        record["cells"] = total["cells"]
    else:
        uris = serving_uris(rows, per_dataset)
        part = read_corpus(shape, uris, workers)
        t2 = time.monotonic()
        cpu = {k: part[k] for k in cpu}
        record["datasets"] = len(rows)
        record["stores_opened"] = len(uris)
        record["cells"] = part["cells"]

    record.update(
        catalog_seconds=catalog_seconds,
        query_seconds=round(t2 - t1, 3),
        cpu_user_seconds=round(cpu["cpu_user"], 3),
        cpu_sys_seconds=round(cpu["cpu_sys"], 3),
        cpu_seconds=round(cpu["cpu_user"] + cpu["cpu_sys"], 3),
        requests=COUNTER.requests,
        mb_payload=round(COUNTER.nbytes / 2**20, 2),
        peak_rss_mb=peak_rss_mb(),
        versions={
            m: importlib.metadata.version(m) for m in ("anndata", "zarr", "pyarrow")
        },
        status="ok",
    )
    return record


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("--shape", choices=SHAPES, required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--cell-type", default=None)
    p.add_argument("--dataset-id", default=None, help="the dataset of `single`")
    p.add_argument("--filter", default=None, help="the filter of `filter`, as JSON")
    p.add_argument("--builder", default=None, help="a query of bench.filters")
    p.add_argument(
        "--per-dataset", action="store_true", help="ignore the consolidated stores"
    )
    p.add_argument("--reader", choices=READERS, default="sdk")
    p.add_argument(
        "--processes", type=int, default=1, help="processes of the corpus shapes"
    )
    p.add_argument("--max-rss-gb", type=float, default=None)
    args = p.parse_args()
    if args.processes > 1 and args.shape not in CORPUS_READS:
        p.error("--processes applies to the corpus shapes only")
    global READER
    READER = args.reader
    limit_mb = None if args.max_rss_gb is None else args.max_rss_gb * 1024
    if limit_mb is not None:
        base = {"catalog": args.catalog, "shape": args.shape, "workers": args.workers}
        guard_memory(limit_mb, base)
    record = run(
        args.catalog,
        args.shape,
        args.workers,
        args.cell_type,
        dataset_id=args.dataset_id,
        obs_filter=json.loads(args.filter) if args.filter else None,
        per_dataset=args.per_dataset,
        processes=args.processes,
        builder=args.builder,
    )
    publish(*final_record(record, limit_mb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
