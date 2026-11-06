"""Run one benchmark query through DuckDB over the obs fragments, and print one JSON line.

This is the `sql` path of the benchmark (`census2_sdk.obs_sql`). The shapes, the
arguments, and the record follow `bench.run`, and each run is a new process. The
requests and the bytes of each kind of object are in the same fields as on the other
paths (`counting.py`). The zarr reads of the stores and the catalog read come from the
same counters; the fragment reads come from DuckDB's HTTP log, so they need no Python
in the read path and do not slow it. The record also holds:

- `sql_requests`: DuckDB's requests by method (GET, HEAD), and `sql_seconds`, the time
  from the query to the last Arrow batch, with the dictionary encoding of the labels.
- `catalog_test_seconds`, the time of the catalog test in DuckDB, and
  `pandas_seconds`, the time to build the pandas frame of `get_obs`.
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
import time
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import anndata

import census2_sdk
from bench import counting
from bench.counting import COUNTER
from bench.run import (
    ORGANISM,
    cpu_seconds,
    final_record,
    get_anndata,
    guard_memory,
    io_read_mb,
    peak_rss_mb,
    pick_cell_type,
    publish,
    source_rows,
)
from census2_sdk import obs_sql

SHAPES = ("single", "cell_type", "query", "count", "corpus_obs")


def sql_stats(stats: list[dict[str, Any]]) -> dict[str, Any]:
    """The DuckDB statistics of one query, for the record."""
    if not stats:
        return {}
    s = stats[0]
    return {
        "sql_requests": s.get("requests", {}),
        "sql_mb": s.get("mb_requested", 0.0),
        "sql_seconds": round(s.get("scan_seconds", 0) + s.get("arrow_seconds", 0), 3),
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
    by: Sequence[str] = ("donor_id",),
    arrow_obs: bool = False,
) -> dict[str, Any]:
    before = io_read_mb()
    c0 = cpu_seconds()
    t0 = time.monotonic()
    census = census2_sdk.open_census([catalog], obs_engine=options)
    rows = census.catalog[census.catalog["organism"] == ORGANISM]
    t1 = time.monotonic()
    c1 = cpu_seconds()
    opened = COUNTER.snapshot()
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

    def default_filter() -> dict[str, Any]:
        return {"cell_type": cell_type or pick_cell_type(source_rows(rows))}

    stats: list[dict[str, Any]] = []
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
        stats = census.last_scan
        record["catalog_candidates"] = stats[0]["catalog_candidates"] if stats else 0
        if isinstance(adata, anndata.AnnData):
            record["stores_opened"] = adata.uns["census2"]["stores_opened"]
            record["cells"] = int(adata.n_obs)
        else:
            record.update(stores_opened=0, cells=0)
    elif shape == "count":
        obs_filter = query if query is not None else default_filter()
        record.update(query=obs_filter, by=list(by))
        frame = census.count_obs(ORGANISM, obs_filter, by=list(by))
        t2 = time.monotonic()
        stats = frame.attrs["census2"]["obs_sql"]
        record.update(stores_opened=0, groups=len(frame))
        record["cells"] = int(frame["n_cells"].sum())
        record["catalog_candidates"] = frame.attrs["census2"]["catalog_candidates"]
    elif shape == "corpus_obs":
        backend: Literal["numpy", "pyarrow"] = "pyarrow" if arrow_obs else "numpy"
        frame = census.get_obs(ORGANISM, dtype_backend=backend)
        t2 = time.monotonic()
        stats = frame.attrs["census2"]["obs_sql"]
        record.update(stores_opened=0, cells=len(frame), dtype_backend=backend)
        record["catalog_candidates"] = frame.attrs["census2"]["catalog_candidates"]
    else:
        raise ValueError(f"no shape {shape!r} on the sql path")
    c2 = cpu_seconds()
    done = COUNTER.snapshot()

    record.update(sql_stats(stats))
    # The fragment reads are DuckDB's, from its HTTP log; the rest are the counters'.
    requests = dict(done["requests_by_kind"])
    heads = dict(done["heads_by_kind"])
    mb = dict(done["mb_by_kind"])
    sql_requests = record.get("sql_requests", {})
    requests["fragment"] = requests.get("fragment", 0) + sum(sql_requests.values())
    heads["fragment"] = heads.get("fragment", 0) + int(sql_requests.get("HEAD", 0))
    mb["fragment"] = round(mb.get("fragment", 0.0) + record.get("sql_mb", 0.0), 3)
    record.update(
        catalog_seconds=round(t1 - t0, 3),
        query_seconds=round(t2 - t1, 3),
        cpu_open_seconds=round(sum(c1) - sum(c0), 3),
        cpu_seconds=round(sum(c2) - sum(c1), 3),
        cpu_user_seconds=round(c2[0] - c1[0], 3),
        cpu_sys_seconds=round(c2[1] - c1[1], 3),
        requests_open=opened["requests"],
        requests=sum(requests.values()) - opened["requests"],
        requests_by_kind=requests,
        heads_by_kind=heads,
        mb_by_kind=mb,
        mb_payload=round(sum(mb.values()), 2),
        peak_rss_mb=peak_rss_mb(),
        machine=platform.machine(),
        versions={
            m: importlib.metadata.version(m)
            for m in ("anndata", "zarr", "pyarrow", "duckdb", "pandas")
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
    p.add_argument("--by", default="donor_id", help="the count columns, joined by +")
    p.add_argument("--label", default=None, help="a name for the query")
    p.add_argument("--arrow-obs", action="store_true")
    p.add_argument("--threads", type=int, default=None, help="DuckDB's threads")
    p.add_argument("--memory-limit", default=None, help="DuckDB's memory_limit")
    p.add_argument("--files", default="catalog", help="obs_sql.Options.files")
    p.add_argument("--glob", default=None)
    p.add_argument("--iceberg", default=None, help="a metadata file or table location")
    p.add_argument("--no-union-by-name", action="store_true")
    p.add_argument("--no-count-requests", action="store_true")
    p.add_argument("--credentials", default="auto")
    args = p.parse_args()
    counting.install_read_counting()
    iceberg = args.iceberg
    if iceberg is not None and not iceberg.endswith(".metadata.json"):
        # A table location: its version hint names the current metadata file.
        from census2_sdk import obs_table

        iceberg = obs_table.current_metadata(iceberg)
    options = obs_sql.Options(
        threads=args.threads,
        memory_limit=args.memory_limit,
        files=args.files,
        glob=args.glob,
        iceberg=iceberg,
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
        args.catalog,
        args.shape,
        args.workers,
        args.cell_type,
        options,
        query=query,
        by=args.by.split("+"),
        arrow_obs=args.arrow_obs,
    )
    if args.label is not None:
        record["label"] = args.label
    publish(*final_record(record, limit_mb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
