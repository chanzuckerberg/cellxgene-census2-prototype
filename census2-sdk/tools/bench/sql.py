"""Run one benchmark query through DuckDB over the obs fragments, and print one JSON line.

This is the `sql` path of the benchmark (`census2_sdk.obs_sql`). The shapes, the
arguments, and the record follow `bench.run`, and each run is a new process. The record
also holds:

- `sql_requests`: the HTTP requests of DuckDB by method, from its HTTP log, and
  `sql_mb`, the MB of their GET responses. `requests` and `mb_payload` count the zarr
  reads of the stores, as in `bench.run`, and `requests_total` adds the two.
- `cpu_seconds`: the user and system CPU time of the process during the query.
- The DuckDB settings of the run (`threads`, `memory_limit`, `files`).

`--files` selects how DuckDB finds the fragments (`obs_sql.Options.files`). `--query`
takes a mapping filter as a JSON object, or a SQL filter as a JSON string.

    python -m bench.sql catalog-103.parquet --shape query \\
        --query '"cell_type = ''neutrophil'' AND nnz >= 2000"' --workers 8 --threads 16
"""

import argparse
import importlib.metadata
import json
import platform
import resource
import time
from collections.abc import Mapping
from typing import Any

import anndata

import census2_sdk
from bench.run import (
    COUNTER,
    ORGANISM,
    final_record,
    guard_memory,
    io_read_mb,
    open_counted,
    peak_rss_mb,
    pick_cell_type,
    publish,
    source_rows,
)
from census2_sdk import catalog as _catalog
from census2_sdk import obs_sql

SHAPES = ("single", "cell_type", "query", "count_donor", "corpus_obs")


def cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime + usage.ru_stime


def sql_stats(stats: list[dict[str, Any]]) -> dict[str, Any]:
    """The DuckDB statistics of one query, for the record."""
    if not stats:
        return {}
    s = stats[0]
    return {
        "sql_requests": s.get("requests", {}),
        "sql_mb": s.get("mb_requested", 0.0),
        "sql_seconds": round(
            s.get("scan_seconds", 0.0) + s.get("arrow_seconds", 0.0), 3
        ),
        "catalog_test_seconds": s.get("catalog_seconds", 0.0),
        "pandas_seconds": s.get("pandas_seconds"),
        "files_read": s.get("files", 0),
        "scan_rows": s.get("rows", 0),
    }


def run(
    catalog: str,
    shape: str,
    workers: int,
    cell_type: str | None,
    options: obs_sql.Options,
    *,
    query: Any = None,
) -> dict[str, Any]:
    before = io_read_mb()
    t0 = time.monotonic()
    census = census2_sdk.open_census([catalog], obs_engine=options)
    rows = census.catalog[census.catalog["organism"] == ORGANISM]
    t1 = time.monotonic()
    cpu0 = cpu_seconds()
    record: dict[str, Any] = {
        "catalog": catalog,
        "catalog_stores": len(rows),
        "shape": shape,
        "workers": workers,
        "via": "sql",
        "engine_open_seconds": census.sql_engine.open_seconds,
        "threads": options.threads,
        "memory_limit": options.memory_limit,
        "files": options.files,
        "union_by_name": options.union_by_name,
    }
    stats: list[dict[str, Any]] = []
    if shape in ("single", "cell_type", "query"):
        if shape == "single":
            obs_filter: Any = {"dataset_id": str(rows["dataset_id"].iloc[0])}
        elif shape == "query":
            if query is None:
                raise ValueError("the query shape needs --query")
            obs_filter = query
            record["query"] = query
        else:
            cell_type = cell_type or pick_cell_type(source_rows(rows))
            obs_filter = {"cell_type": cell_type}
            record["cell_type"] = cell_type
        # The SDK looks up open_store by name, so patch both names, and restore them.
        original = _catalog.open_store
        census2_sdk.__dict__["open_store"] = open_counted
        _catalog.open_store = open_counted
        try:
            adata = census.get_anndata(ORGANISM, obs_filter, max_workers=workers)
        finally:
            census2_sdk.__dict__["open_store"] = original
            _catalog.open_store = original
        t2 = time.monotonic()
        stats = census.last_scan
        record["catalog_candidates"] = (
            stats[0].get("catalog_candidates") if stats else 0
        )
        if isinstance(adata, anndata.AnnData):
            record["stores_opened"] = adata.uns["census2"]["stores_opened"]
            record["cells"] = int(adata.n_obs)
        else:
            record.update(stores_opened=0, cells=0)
    elif shape == "count_donor":
        obs_filter = (
            query
            if query is not None
            else {"cell_type": cell_type or pick_cell_type(source_rows(rows))}
        )
        record["query"] = obs_filter
        frame = census.count_obs(ORGANISM, obs_filter, by=["donor_id"])
        t2 = time.monotonic()
        stats = frame.attrs["census2"]["obs_sql"]
        record.update(stores_opened=0, groups=len(frame))
        record["cells"] = int(frame["n_cells"].sum())
        record["catalog_candidates"] = frame.attrs["census2"]["catalog_candidates"]
    elif shape == "corpus_obs":
        frame = census.get_obs(ORGANISM)
        t2 = time.monotonic()
        stats = frame.attrs["census2"]["obs_sql"]
        record.update(stores_opened=0, cells=len(frame))
        record["catalog_candidates"] = frame.attrs["census2"]["catalog_candidates"]
    else:
        raise ValueError(f"no shape {shape!r} on the sql path")
    cpu1 = cpu_seconds()

    record.update(sql_stats(stats))
    sql_requests = sum(int(n) for n in record.get("sql_requests", {}).values())
    record.update(
        catalog_seconds=round(t1 - t0, 3),
        query_seconds=round(t2 - t1, 3),
        cpu_seconds=round(cpu1 - cpu0, 2),
        requests=COUNTER.requests,
        mb_payload=round(COUNTER.nbytes / 2**20, 2),
        requests_total=COUNTER.requests + sql_requests,
        peak_rss_mb=peak_rss_mb(),
        machine=platform.machine(),
        versions={
            m: importlib.metadata.version(m)
            for m in ("anndata", "zarr", "pyarrow", "duckdb")
        },
        status="ok",
    )
    after = io_read_mb()
    if before is not None and after is not None:
        record["mb_read_process"] = round(after - before, 2)
    return record


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("--shape", choices=SHAPES, required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--cell-type", default=None)
    p.add_argument("--max-rss-gb", type=float, default=None)
    p.add_argument("--query", default=None, help="a mapping or a SQL string, as JSON")
    p.add_argument("--label", default=None, help="a name for the query")
    p.add_argument("--threads", type=int, default=None, help="DuckDB's threads")
    p.add_argument("--memory-limit", default=None, help="DuckDB's memory_limit")
    p.add_argument("--files", default="catalog", help="obs_sql.Options.files")
    p.add_argument("--glob", default=None)
    p.add_argument("--iceberg", default=None, help="a metadata file")
    p.add_argument("--no-union-by-name", action="store_true")
    p.add_argument("--no-count-requests", action="store_true")
    p.add_argument("--credentials", default="auto")
    args = p.parse_args()
    options = obs_sql.Options(
        threads=args.threads,
        memory_limit=args.memory_limit,
        files=args.files,
        glob=args.glob,
        iceberg=args.iceberg,
        union_by_name=not args.no_union_by_name,
        count_requests=not args.no_count_requests,
        credentials=args.credentials,
    )
    limit_mb = None if args.max_rss_gb is None else args.max_rss_gb * 1024
    if limit_mb is not None:
        base = {"catalog": args.catalog, "shape": args.shape, "workers": args.workers}
        guard_memory(limit_mb, {**base, "via": "sql"})
    query = None if args.query is None else json.loads(args.query)
    if query is not None and not isinstance(query, (str, Mapping)):
        p.error("--query takes a JSON object or a JSON string")
    record = run(
        args.catalog, args.shape, args.workers, args.cell_type, options, query=query
    )
    if args.label is not None:
        record["label"] = args.label
    publish(*final_record(record, limit_mb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
