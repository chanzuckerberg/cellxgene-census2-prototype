"""Time the catalog test of each filter language, and compare the datasets they keep.

For each catalog and each named query, the script times three catalog tests over the
rows of the organism, in this process, as the median of `--repeats`:

- `sql`: `obs_sql.catalog_test` in DuckDB, over the catalog table of an engine.
- `iceberg`: `expressions.can_match` in Python, over each row.
- `mapping`: `catalog.exclude` in Python, for a query with a mapping form.

The record holds the seconds and the number of rows that each test keeps. A query has
the forms of `bench.matrix`: a mapping, or `{"mapping": ..., "iceberg": ..., "sql": ...}`.

    python -m bench.catalog_test catalogs/catalog-103.parquet catalogs/catalog-2060.parquet \\
        --queries results/obs_duckdb/queries.json --out results/obs_duckdb/catalog_test.jsonl
"""

import argparse
import functools
import json
import pathlib
import statistics
import time
from collections.abc import Callable
from typing import Any, cast

import census2_sdk
from bench.matrix import query_form
from census2_sdk import catalog as _catalog
from census2_sdk import expressions, obs_sql

ORGANISM = "homo_sapiens"


def timed(fn: Callable[[], int], repeats: int) -> tuple[float, int]:
    seconds, kept = [], 0
    for _ in range(repeats):
        t = time.perf_counter()
        kept = fn()
        seconds.append(time.perf_counter() - t)
    return round(statistics.median(seconds), 5), kept


def kept_sql(engine: obs_sql.Engine, pred: obs_sql.Predicate) -> int:
    return len(engine.candidates(ORGANISM, pred)[0])


def kept_iceberg(records: list[dict[str, Any]], expr: Any) -> int:
    return sum(1 for r in records if expressions.can_match(expr, r))


def kept_mapping(catalog: Any, mapping: dict[str, Any]) -> int:
    return len(_catalog.exclude(catalog, ORGANISM, mapping))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalogs", nargs="+")
    p.add_argument("--queries", type=pathlib.Path, required=True)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    queries = json.loads(args.queries.read_text())
    lines = []
    for path in args.catalogs:
        census = census2_sdk.open_census([path], obs_engine="duckdb")
        engine = census.sql_engine
        scoped = census.catalog[census.catalog["organism"] == ORGANISM]
        records = cast("list[dict[str, Any]]", scoped.to_dict("records"))
        for name, query in queries.items():
            record: dict[str, Any] = {
                "catalog": path,
                "rows": len(scoped),
                "label": name,
            }
            sql = query_form(query, "sql")
            if sql is not None:
                pred = engine.predicate(sql)
                record["sql_seconds"], record["sql_kept"] = timed(
                    functools.partial(kept_sql, engine, pred), args.repeats
                )
            iceberg = query_form(query, "table")
            if iceberg is not None:
                expr = expressions.from_filter(iceberg)
                record["iceberg_seconds"], record["iceberg_kept"] = timed(
                    functools.partial(kept_iceberg, records, expr), args.repeats
                )
            mapping = query_form(query, "stores")
            if isinstance(mapping, dict) and "stores" not in query:
                record["mapping_seconds"], record["mapping_kept"] = timed(
                    functools.partial(kept_mapping, census.catalog, mapping),
                    args.repeats,
                )
            print(json.dumps(record), flush=True)
            lines.append(json.dumps(record))
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("a") as fh:
            fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
