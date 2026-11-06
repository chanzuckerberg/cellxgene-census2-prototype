"""Run one benchmark query in this process and print one JSON line.

Each run is a new process, so no connection or metadata cache carries over from
an earlier run. The shapes are:

- `single`: `get_anndata` for one dataset. It opens one store at any catalog size.
- `cell_type`: `get_anndata` for one cell type across the catalog. The cell type
  is the one whose share of the source datasets is nearest to 15 percent. The
  copies are not part of that share.
- `open`: open every store of the catalog and read nothing more.
- `corpus_obs`: open every store and read its whole `obs`.
- `corpus_x`: open every store and read all of `X`, one shard at a time, and keep
  nothing. This measures the time of a full read of `X` without a result in memory.
- `query`: `get_anndata` for the filter in `--query`: a mapping filter as a JSON
  object, or an expression as a JSON string (the table path only).
- `count`: the number of cells of each group of the `--by` columns (default
  `donor_id`) for the filter in `--query`, or for the `cell_type` shape's filter.

`--workers` is the largest number of stores read at the same time. The record
holds the elapsed time and the CPU time (user and system, of all threads) of the
catalog read and of the query. It also holds the number of stores and cells, the
requests and the bytes of each kind of object (`counting.py`), and the peak memory.

`--obs-table` names the `obs` table of the catalog. The query shapes then go through
the table (`census2_sdk.obs_table`): one scan selects the cells, and the SDK opens only
the stores that hold one. `corpus_obs` reads the table and no store (`--arrow-obs`
returns Arrow-backed columns), and `count` counts from the table. The record then
also holds the plan and scan times, the data files that the plan kept, and the
estimated bytes of the scan. `--x-from-indptr` reads the whole `X/indptr` of each
store, as the read without the table does, instead of the ranges in the table.
`--scan-workers` sets pyiceberg's thread pool. On Linux, the record holds the bytes
that the process read (`/proc/self/io`), which counts every kind of read.

The store reads are always counted. `--count-arrow` also counts the reads of the
table metadata and of the fragments, through a Python file system that can slow the
scan; see `counting.py`.

Without a table, `count` reads the `obs` of each store that the catalog keeps, as
`get_anndata` reads it (every column), and selects and groups the cells in pandas.
`--narrow` reads only the columns of the filter and of `--by`. An expression filter
then goes through the same catalog test as on the table path
(`expressions.can_match`), and `evaluate` applies it to each store's `obs`.

A catalog of copies (`catalogs.py`) gives each copy a new `dataset_id`. The run then
accepts a store whose `obs` names the source dataset of its copy.

`--max-rss-gb` is a memory limit. When the peak memory passes it, the process
prints a record with the status `memory_limit` and stops with exit code 3.

    python -m bench.run catalog-500.parquet --shape cell_type --workers 32
"""

import argparse
import importlib.metadata
import json
import os
import pathlib
import platform
import re
import resource
import sys
import threading
import time
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Literal, cast

import anndata
import numpy as np
import pandas as pd
import zarr
import zarr.storage
from anndata.io import read_elem
from zarr.abc.store import Store

import census2_sdk
from bench.counting import COUNTER, CountingStore
from census2_sdk import catalog as _catalog

SHAPES = (
    "single",
    "cell_type",
    "open",
    "corpus_obs",
    "corpus_x",
    "query",
    "count",
)
ORGANISM = "homo_sapiens"
MEMORY_LIMIT_EXIT = 3


def open_counted(uri: str | os.PathLike[str]) -> zarr.Group:
    """`census2_sdk.catalog.open_store`, through a counting store."""
    if _catalog.is_remote(uri):
        import obstore.store

        base: Store = zarr.storage.ObjectStore(
            obstore.store.from_url(str(uri)), read_only=True
        )
    else:
        base = zarr.storage.LocalStore(str(uri), read_only=True)
    return zarr.open_group(CountingStore(base), mode="r")


def read_obs(uri: str, columns: Sequence[str] | None = None) -> pd.DataFrame:
    """A store's `obs`, column by column, without the index, as `get_anndata` reads it.

    `columns` names the columns to read; None reads every column.
    """
    group = open_counted(uri)["obs"]
    assert isinstance(group, zarr.Group), "obs must be an encoded dataframe"
    names = cast("list[str]", group.attrs["column-order"])
    wanted = names if columns is None else [n for n in names if n in set(columns)]
    return pd.DataFrame({name: read_elem(group[name]) for name in wanted})


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


def source_rows(rows: pd.DataFrame) -> pd.DataFrame:
    """The rows of the source stores, without the copies that `catalogs` adds."""
    return rows[rows["copy"] == 0] if "copy" in rows.columns else rows


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


def cpu_seconds() -> tuple[float, float]:
    """The user and the system CPU time of this process so far, all threads."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime, usage.ru_stime


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


def accept_copies() -> None:
    """Accept a store whose `obs` names the source dataset of a copy's new `dataset_id`."""
    check = census2_sdk._check_dataset_id

    def lenient(obs: pd.DataFrame, dataset_id: str, store_uri: str) -> None:
        check(obs, re.sub(r"-c\d\d$", "", dataset_id), store_uri)

    census2_sdk.__dict__["_check_dataset_id"] = lenient


def io_read_mb() -> float | None:
    """The bytes that this process has read, from files and sockets. Linux only."""
    try:
        text = pathlib.Path("/proc/self/io").read_text()
    except OSError:
        return None
    fields = dict(line.split(": ") for line in text.splitlines())
    return round(int(fields["rchar"]) / 2**20, 2)


def table_stats(stats: list[dict[str, Any]]) -> dict[str, Any]:
    """The scan statistics of all sources, summed."""
    keys = ["plan_seconds", "scan_seconds", "files_total", "files_planned"]
    keys += ["files_after_bounds", "rows", "data_bytes"]
    out = {k: sum(s.get(k, 0) for s in stats) for k in keys}
    out["mb_scan_estimate"] = round((out.pop("data_bytes")) / 2**20, 2)
    return out


def evaluate(expr: Any, obs: pd.DataFrame) -> np.ndarray:
    """The cells of one store's `obs` that satisfy an unbound Iceberg expression.

    Apply `expressions.with_members` first: this tests strings as strings. A term on a
    column that `obs` does not hold selects no cell, as on the store path.
    """
    from pyiceberg import expressions as ex

    n = len(obs)
    if isinstance(expr, ex.AlwaysTrue):
        return np.ones(n, dtype=bool)
    if isinstance(expr, ex.AlwaysFalse):
        return np.zeros(n, dtype=bool)
    if isinstance(expr, ex.And):
        return np.asarray(evaluate(expr.left, obs) & evaluate(expr.right, obs))
    if isinstance(expr, ex.Or):
        return np.asarray(evaluate(expr.left, obs) | evaluate(expr.right, obs))
    if isinstance(expr, ex.Not):
        return np.asarray(~evaluate(expr.child, obs))
    name = str(expr.term.name)
    if name not in obs.columns:
        return np.zeros(n, dtype=bool)
    column = obs[name]
    present = column.notna().to_numpy()
    if isinstance(expr, ex.IsNull):
        return ~present
    if isinstance(expr, ex.NotNull):
        return present
    if isinstance(expr, (ex.In, ex.NotIn)):
        values = [v.value for v in expr.literals]
        hit = column.isin(values).to_numpy()
        return hit if isinstance(expr, ex.In) else ~hit & present
    value = expr.literal.value
    if isinstance(expr, (ex.StartsWith, ex.NotStartsWith)):
        hit = column.astype(str).str.startswith(value).to_numpy() & present
        return hit if isinstance(expr, ex.StartsWith) else ~hit & present
    compare: dict[type, Callable[[Any, Any], Any]] = {
        ex.EqualTo: lambda a, b: a == b,
        ex.NotEqualTo: lambda a, b: a != b,
        ex.LessThan: lambda a, b: a < b,
        ex.LessThanOrEqual: lambda a, b: a <= b,
        ex.GreaterThan: lambda a, b: a > b,
        ex.GreaterThanOrEqual: lambda a, b: a >= b,
    }
    labels: Any = column
    if isinstance(column.dtype, pd.CategoricalDtype):
        labels = column.astype(str)
    return np.asarray(compare[type(expr)](labels, value), dtype=bool) & present


def count_from_stores(
    census: census2_sdk.Census,
    obs_filter: Any,
    by: Sequence[str],
    workers: int,
    *,
    narrow: bool = False,
) -> tuple[int, int, int]:
    """The `count` shape without a table: read the `obs` of each candidate store.

    Return the number of groups, the number of stores read, and the cells.
    """
    if isinstance(obs_filter, dict):
        canonical = _catalog.canonicalize(obs_filter)
        rows = _catalog.exclude(census.catalog, ORGANISM, canonical)
        needed = {*canonical, *by} - {"dataset_id"}

        def select(obs: pd.DataFrame, row: dict[str, Any]) -> np.ndarray:
            return census2_sdk._mask(obs, canonical)

    else:
        from census2_sdk import expressions

        expr = expressions.from_filter(obs_filter)
        rows = census._candidates(ORGANISM, expr)
        needed = {*expressions.columns(expr), *by}

        def select(obs: pd.DataFrame, row: dict[str, Any]) -> np.ndarray:
            return evaluate(expressions.with_members(expr, [row]), obs)

    def one(row: dict[str, Any]) -> pd.Series:
        obs = read_obs(str(row["store_uri"]), sorted(needed) if narrow else None)
        kept = obs[select(obs, row)]
        return kept.groupby(list(by), observed=True).size()

    records = cast("list[dict[str, Any]]", rows.to_dict("records"))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(one, records))
    groups = 0
    if parts:
        groups = len(pd.concat(parts).groupby(level=list(range(len(by)))).sum())
    return groups, len(rows), int(sum(p.sum() for p in parts))


def get_anndata(
    census: census2_sdk.Census, obs_filter: Any, workers: int
) -> anndata.AnnData | None:
    """`Census.get_anndata` with every store open counted."""
    # The SDK imports open_store by name, so patch both names, and restore them.
    original = _catalog.open_store
    census2_sdk.__dict__["open_store"] = open_counted
    _catalog.open_store = open_counted
    try:
        return census.get_anndata(ORGANISM, obs_filter, max_workers=workers)
    finally:
        census2_sdk.__dict__["open_store"] = original
        _catalog.open_store = original


def run(
    catalog: str,
    shape: str,
    workers: int,
    cell_type: str | None,
    *,
    query: Any = None,
    obs_table: str | None = None,
    by: Sequence[str] = ("donor_id",),
    narrow: bool = False,
    arrow_obs: bool = False,
) -> dict[str, Any]:
    before = io_read_mb()
    c0 = cpu_seconds()
    t0 = time.monotonic()
    census = census2_sdk.open_census(
        [catalog], obs_tables=None if obs_table is None else [obs_table]
    )
    rows = census.catalog[census.catalog["organism"] == ORGANISM]
    t1 = time.monotonic()
    c1 = cpu_seconds()
    opened = COUNTER.snapshot()
    record: dict[str, Any] = {
        "catalog": catalog,
        "catalog_stores": len(rows),
        "shape": shape,
        "workers": workers,
        "via": "stores" if obs_table is None else "table",
    }
    if obs_table is not None:
        record["x_from_ranges"] = census2_sdk.X_FROM_RANGES
        record["key_first"] = census2_sdk.SCAN_KEY_FIRST
    if "copy" in rows.columns:
        accept_copies()

    def default_filter() -> dict[str, Any]:
        return {"cell_type": cell_type or pick_cell_type(source_rows(rows))}

    if shape in ("single", "cell_type", "query"):
        if shape == "single":
            obs_filter: Any = {"dataset_id": str(rows["dataset_id"].iloc[0])}
        elif shape == "query":
            if query is None:
                raise ValueError("the query shape needs --query")
            obs_filter = query
        else:
            obs_filter = default_filter()
            record["cell_type"] = obs_filter["cell_type"]
        record["query"] = obs_filter
        adata = get_anndata(census, obs_filter, workers)
        t2 = time.monotonic()
        if obs_table is not None:
            record["catalog_candidates"] = len(
                census._candidates(ORGANISM, _expression(obs_filter))
            )
            uns = adata.uns["census2"] if isinstance(adata, anndata.AnnData) else {}
            record["stores_opened"] = uns.get("stores_opened", 0)
            record.update(table_stats(uns.get("obs_table", census.last_scan)))
        else:
            candidates = _catalog.exclude(
                census.catalog, ORGANISM, _catalog.canonicalize(obs_filter)
            )
            record["catalog_candidates"] = record["stores_opened"] = len(candidates)
        record["cells"] = int(adata.n_obs) if isinstance(adata, anndata.AnnData) else 0
    elif shape == "count":
        obs_filter = query if query is not None else default_filter()
        record.update(query=obs_filter, by=list(by))
        if obs_table is None:
            groups, read, cells = count_from_stores(
                census, obs_filter, by, workers, narrow=narrow
            )
            t2 = time.monotonic()
            record.update(stores_opened=read, groups=groups, cells=cells)
            record.update(catalog_candidates=read, narrow=narrow)
        else:
            frame = census.count_obs(ORGANISM, obs_filter, by=list(by))
            t2 = time.monotonic()
            record.update(stores_opened=0, groups=len(frame))
            record["cells"] = int(frame["n_cells"].sum())
            record["catalog_candidates"] = frame.attrs["census2"]["catalog_candidates"]
            record.update(table_stats(frame.attrs["census2"]["obs_table"]))
    elif shape == "corpus_obs" and obs_table is not None:
        backend: Literal["numpy", "pyarrow"] = "pyarrow" if arrow_obs else "numpy"
        frame = census.get_obs(ORGANISM, dtype_backend=backend)
        t2 = time.monotonic()
        attrs = frame.attrs["census2"]
        record.update(stores_opened=0, cells=len(frame), dtype_backend=backend)
        record.update(table_stats(attrs["obs_table"]))
        record.update(
            table_seconds=attrs["table_seconds"], build_seconds=attrs["build_seconds"]
        )
    else:
        reader: Callable[[str], int] = {
            "open": open_only,
            "corpus_obs": lambda uri: len(read_obs(uri)),
            "corpus_x": read_x,
        }[shape]
        uris = [str(u) for u in rows["store_uri"]]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            cells = sum(pool.map(reader, uris))
        t2 = time.monotonic()
        record["stores_opened"] = len(uris)
        record["cells"] = cells

    c2 = cpu_seconds()
    done = COUNTER.snapshot()
    record.update(
        catalog_seconds=round(t1 - t0, 3),
        query_seconds=round(t2 - t1, 3),
        cpu_open_seconds=round(sum(c1) - sum(c0), 3),
        cpu_seconds=round(sum(c2) - sum(c1), 3),
        cpu_user_seconds=round(c2[0] - c1[0], 3),
        cpu_sys_seconds=round(c2[1] - c1[1], 3),
        requests_open=opened["requests"],
        requests=done["requests"] - opened["requests"],
        requests_by_kind=done["requests_by_kind"],
        heads_by_kind=done["heads_by_kind"],
        mb_by_kind=done["mb_by_kind"],
        mb_payload=round(sum(done["mb_by_kind"].values()), 2),
        peak_rss_mb=peak_rss_mb(),
        machine=platform.machine(),
        versions={
            m: importlib.metadata.version(m)
            for m in ("anndata", "zarr", "pyarrow", "pyiceberg", "pandas")
        },
        status="ok",
    )
    for table in census._tables.values():
        from census2_sdk import obs_table as _obs_table

        record.update(_obs_table.manifest_stats(table))
    after = io_read_mb()
    if before is not None and after is not None:
        record["mb_read_process"] = round(after - before, 2)
    return record


def _expression(obs_filter: Any) -> Any:
    from census2_sdk import expressions

    return expressions.from_filter(obs_filter)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("--shape", choices=SHAPES, required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--cell-type", default=None)
    p.add_argument("--max-rss-gb", type=float, default=None)
    p.add_argument("--query", default=None, help="a filter, as JSON")
    p.add_argument("--by", default="donor_id", help="the count columns, joined by +")
    p.add_argument("--label", default=None, help="a name for the query")
    p.add_argument("--obs-table", default=None)
    p.add_argument("--x-from-indptr", action="store_true")
    p.add_argument("--key-first", action="store_true")
    p.add_argument("--scan-workers", type=int, default=None)
    p.add_argument("--count-arrow", action="store_true")
    p.add_argument("--narrow", action="store_true")
    p.add_argument("--arrow-obs", action="store_true")
    args = p.parse_args()
    from bench import counting

    counting.install_read_counting()
    if args.count_arrow:
        counting.install_arrow_counting()
    if args.scan_workers is not None:
        os.environ["PYICEBERG_MAX_WORKERS"] = str(args.scan_workers)
    if args.x_from_indptr:
        census2_sdk.__dict__["X_FROM_RANGES"] = False
    if args.key_first:
        census2_sdk.__dict__["SCAN_KEY_FIRST"] = True
    if args.obs_table is not None or args.query is not None:
        # Import the table modules before the timed part.
        from census2_sdk import expressions, obs_table  # noqa: F401
    limit_mb = None if args.max_rss_gb is None else args.max_rss_gb * 1024
    if limit_mb is not None:
        base = {"catalog": args.catalog, "shape": args.shape, "workers": args.workers}
        guard_memory(limit_mb, base)
    query = None if args.query is None else json.loads(args.query)
    record = run(
        args.catalog,
        args.shape,
        args.workers,
        args.cell_type,
        query=query,
        obs_table=args.obs_table,
        by=args.by.split("+"),
        narrow=args.narrow,
        arrow_obs=args.arrow_obs,
    )
    record["count_arrow"] = args.count_arrow
    if args.label is not None:
        record["label"] = args.label
    if args.scan_workers is not None:
        record["scan_workers"] = args.scan_workers
    publish(*final_record(record, limit_mb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
