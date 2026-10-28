"""Read the obs table with one of four readers, and print one JSON record.

The readers are pyiceberg (`table.scan`, which the SDK uses), DuckDB's `iceberg`
extension, Polars' `scan_iceberg`, and pyarrow's dataset over the list of fragments in
the catalog, which needs no Iceberg metadata. DuckDB and Polars are not dependencies of
the SDK, so install them to run those readers. Each run is a new process, as in
`bench.run`, so the record holds the peak memory of one reader. The queries are:

- `select`: the rows of a mapping filter (`--query`), with the cell key and the row
  columns, in the datasets that the catalog keeps.
- `count`: the number of cells of each `donor_id`, over all cells.
- `all`: every column of every cell.

    python -m bench.readers catalog.parquet s3://bucket/census/obs/homo_sapiens \\
        --reader duckdb --query-shape select --query '{"tissue": ["lung"]}'
"""

import argparse
import importlib
import json
import os
import time
from typing import Any, cast

import census2_dataset_writer
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.dataset as ds  # type: ignore[import-untyped]
import pyarrow.fs as pafs  # type: ignore[import-untyped]

from bench.run import io_read_mb, peak_rss_mb
from census2_sdk import catalog as _catalog
from census2_sdk import expressions, obs_table

READERS = ("pyiceberg", "duckdb", "polars", "pyarrow")
KEY = ["dataset_id", "observation_joinid", *obs_table.ROW_COLUMNS]


def sql_literal(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return "'" + str(value).replace("'", "''") + "'"


def candidates(catalog_path: str, query: dict[str, Any]) -> list[dict[str, Any]]:
    """The catalog rows that the mapping filter keeps."""
    rows = _catalog.load(_catalog.read(catalog_path))
    kept = _catalog.exclude(rows, "homo_sapiens", query)
    return cast("list[dict[str, Any]]", kept.to_dict("records"))


def where_sql(expr: Any) -> str:
    """An Iceberg expression of `IN` and `AND` terms as SQL. `to_text` is close to SQL."""
    return expressions.to_text(expr)


def run_pyiceberg(location: str, shape: str, expr: Any) -> Any:
    table = obs_table.open_table(location)
    selected: tuple[str, ...] = ("*",)
    if shape == "select":
        names = [f.name for f in table.schema().fields]
        selected = tuple(names)
    elif shape == "count":
        selected = ("donor_id",)
    arrow, stats = obs_table.scan(table, expr, selected)
    if shape == "count":
        arrow = arrow.group_by(["donor_id"]).aggregate([([], "count_all")])
    return arrow.num_rows, stats.as_dict()


def run_duckdb(location: str, shape: str, expr: Any, threads: int) -> Any:
    duckdb: Any = importlib.import_module("duckdb")  # optional: not a dependency

    metadata = obs_table.open_table(location).metadata_location.removeprefix("file://")
    con = duckdb.connect()
    con.execute(f"SET threads={threads}")
    con.execute("INSTALL iceberg; INSTALL httpfs; LOAD iceberg; LOAD httpfs")
    if metadata.startswith("s3://"):
        region = os.environ.get("AWS_REGION", "us-west-2")
        con.execute(
            f"CREATE SECRET (TYPE s3, PROVIDER credential_chain, REGION '{region}')"
        )
    source = f"iceberg_scan('{metadata}')"
    if shape == "select":
        sql = f"SELECT * FROM {source} WHERE {where_sql(expr)}"
    elif shape == "count":
        sql = f"SELECT donor_id, count(*) FROM {source} GROUP BY donor_id"
    else:
        sql = f"SELECT * FROM {source}"
    arrow: Any = con.execute(sql).to_arrow_table()
    return arrow.num_rows, {}


def run_polars(location: str, shape: str, expr: Any) -> Any:
    pl: Any = importlib.import_module("polars")  # optional: not a dependency

    table = obs_table.open_table(location)
    frame = pl.scan_iceberg(table)
    if shape == "select":
        frame = frame.filter(pl.sql_expr(where_sql(expr)))
    elif shape == "count":
        frame = frame.group_by("donor_id").len()
    out = frame.collect()
    return out.height, {}


def run_pyarrow(catalog_path: str, shape: str, expr: Any, rows: Any) -> Any:
    """pyarrow's dataset over the fragments of the catalog rows, without Iceberg."""
    from pyiceberg.expressions.visitors import bind
    from pyiceberg.io.pyarrow import expression_to_pyarrow

    paths = [
        census2_dataset_writer.fragment_location(str(r["store_uri"])) for r in rows
    ]
    region = os.environ.get("AWS_REGION", "us-west-2")
    if _catalog.is_remote(paths[0]):
        fs = pafs.S3FileSystem(region=region)
        paths = [p.removeprefix("s3://") for p in paths]
    else:
        fs = pafs.LocalFileSystem()
    dataset = ds.dataset(paths, format="parquet", filesystem=fs)
    filter_ = None
    if shape == "select":
        schema = obs_table.open_table(
            os.environ["OBS_TABLE"]
        ).schema()  # to bind the expression only
        filter_ = expression_to_pyarrow(bind(schema, expr, case_sensitive=True))
    columns = ["donor_id"] if shape == "count" else None
    arrow = dataset.to_table(filter=filter_, columns=columns)
    if shape == "count":
        arrow = arrow.group_by(["donor_id"]).aggregate([([], "count_all")])
    return arrow.num_rows, {"files": len(paths)}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("table")
    p.add_argument("--reader", choices=READERS, required=True)
    p.add_argument("--query-shape", choices=("select", "count", "all"), required=True)
    p.add_argument("--query", default="{}")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 8)
    args = p.parse_args()
    os.environ["OBS_TABLE"] = args.table
    os.environ.setdefault("PYICEBERG_MAX_WORKERS", str(args.threads))
    pa.set_cpu_count(args.threads)
    pa.set_io_thread_count(max(8, args.threads))
    query = dict(json.loads(args.query))
    before = io_read_mb()
    t0 = time.monotonic()
    rows = candidates(args.catalog, query) if args.query_shape == "select" else []
    if args.query_shape == "select":
        ids = [str(r["dataset_id"]) for r in rows]
        expr = expressions.from_filter({**query, "dataset_id": ids})
        expr = expressions.with_members(expr, rows)
    else:
        expr = expressions.from_filter(None)
        rows = candidates(args.catalog, {})
    if args.reader == "pyiceberg":
        n, extra = run_pyiceberg(args.table, args.query_shape, expr)
    elif args.reader == "duckdb":
        n, extra = run_duckdb(args.table, args.query_shape, expr, args.threads)
    elif args.reader == "polars":
        n, extra = run_polars(args.table, args.query_shape, expr)
    else:
        n, extra = run_pyarrow(args.catalog, args.query_shape, expr, rows)
    seconds = time.monotonic() - t0
    record = {
        "reader": args.reader,
        "query_shape": args.query_shape,
        "query": query,
        "threads": args.threads,
        "rows": n,
        "seconds": round(seconds, 3),
        "peak_rss_mb": peak_rss_mb(),
        **{k: v for k, v in extra.items() if k != "plan"},
    }
    after = io_read_mb()
    if before is not None and after is not None:
        record["mb_read_process"] = round(after - before, 2)
    print(json.dumps(record), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
