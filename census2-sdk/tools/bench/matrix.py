"""Run the benchmark matrix: each combination in a new process.

The matrix is sizes x shapes x workers x repeats. A run that exceeds `--timeout`
gets the status `exceeded`. A run that passes `--max-rss-gb` gets the status
`memory_limit` from the run itself. A run that fails in any other way gets the
status `error` and the last lines of its output. Records append to `--out` as JSON lines,
so a stopped matrix keeps its finished runs.

    python -m bench.matrix --catalogs catalogs/ --sizes 50,100 \\
        --shapes single,cell_type,open,corpus_obs --workers 1,32 --repeats 3 \\
        --where laptop --out results/bench-laptop.jsonl

`--queries` names a JSON file of named mapping filters; a shape `query:<name>` runs
one of them, and `count_donor:<name>` counts its cells by donor. `--obs-table` gives
the `obs` table of each catalog, with `{size}` in place of the size, and `--vias`
selects `stores`, `table`, `sql`, or more than one. `--x-from-indptr` and
`--scan-workers` pass through to each run.

A query can also give one form for each via: `{"mapping": {...}, "iceberg": "...",
"sql": "..."}`. A via without a form skips the query, for example the store path for
a numeric range. The `sql` via runs `bench.sql` (DuckDB over the fragments), with
`--duckdb-threads`, `--duckdb-memory-limit`, `--files`, and `--glob`. The other vias
run `bench.counted`, which adds the CPU time to the record of `bench.run`;
`--count-arrow-requests` also counts the S3 requests of pyarrow on those runs.
"""

import argparse
import json
import pathlib
import subprocess
import sys
from typing import Any

import pandas as pd

from bench.run import MEMORY_LIMIT_EXIT, pick_cell_type, source_rows


def one(
    catalog: pathlib.Path,
    shape: str,
    workers: int,
    timeout: float,
    cell_type: str,
    max_rss_gb: float | None = None,
    extra: list[str] | None = None,
) -> dict[str, Any]:
    module = "bench.counted"
    if extra and "--sql-via" in extra:
        module = "bench.sql"
        extra = [a for a in extra if a != "--sql-via"]
    command = [
        sys.executable, "-m", module, str(catalog),
        "--shape", shape, "--workers", str(workers), "--cell-type", cell_type,
    ]  # fmt: skip
    if max_rss_gb is not None:
        command += ["--max-rss-gb", str(max_rss_gb)]
    command += extra or []
    base: dict[str, Any] = {"catalog": str(catalog), "shape": shape, "workers": workers}
    try:
        done = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return {**base, "status": "exceeded", "timeout_seconds": timeout}
    lines = done.stdout.strip().splitlines()
    if done.returncode == MEMORY_LIMIT_EXIT and lines:
        record = dict(json.loads(lines[-1]))
        if record.get("status") == "memory_limit":
            return record
    if done.returncode != 0:
        tail = (done.stderr or done.stdout).strip().splitlines()[-5:]
        return {**base, "status": "error", "error": tail}
    return dict(json.loads(done.stdout.strip().splitlines()[-1]))


def query_form(query: dict[str, Any], via: str) -> Any:
    """The form of a named query for one via, or None when the via has none."""
    forms = {"mapping", "iceberg", "sql"}
    if not set(query) <= forms:
        return query  # a mapping filter, for every via
    key = {"stores": "mapping", "table": "iceberg", "sql": "sql"}[via]
    return query.get(key, query.get("mapping"))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--catalogs", type=pathlib.Path, required=True)
    p.add_argument("--sizes", default="50,100,500,1000,2000")
    p.add_argument("--shapes", default="single,cell_type,open,corpus_obs,corpus_x")
    p.add_argument("--workers", default="1,32")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=float, default=1800)
    p.add_argument("--max-rss-gb", type=float, default=None)
    p.add_argument("--where", required=True, help="a label for the machine")
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument(
        "--cell-type", default=None, help="default: pick from the largest catalog"
    )
    p.add_argument("--queries", type=pathlib.Path, default=None)
    p.add_argument("--obs-table", default=None, help="with {size} for the size")
    p.add_argument("--vias", default="stores")
    p.add_argument("--x-from-indptr", action="store_true")
    p.add_argument("--key-first", action="store_true")
    p.add_argument("--scan-workers", type=int, default=None)
    p.add_argument("--duckdb-threads", default="", help="a list: one run each")
    p.add_argument("--duckdb-memory-limit", default=None)
    p.add_argument("--files", default="catalog", help="obs_sql.Options.files")
    p.add_argument("--glob", default=None, help="with {size} for the size")
    p.add_argument("--count-arrow-requests", action="store_true")
    args = p.parse_args()
    queries = {} if args.queries is None else json.loads(args.queries.read_text())

    # One cell type for every size, so that each size runs the same query.
    sizes = args.sizes.split(",")
    largest = args.catalogs / f"catalog-{max(sizes, key=int)}.parquet"
    cell_type = args.cell_type or pick_cell_type(source_rows(pd.read_parquet(largest)))
    print(f"cell type: {cell_type}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for size in sizes:
        catalog = args.catalogs / f"catalog-{size}.parquet"
        for name in args.shapes.split(","):
            shape, _, label = name.partition(":")
            for via in args.vias.split(","):
                extra = []
                if label:
                    form = query_form(queries[label], via)
                    if form is None:
                        continue
                    extra += ["--query", json.dumps(form), "--label", label]
                if via == "sql":
                    extra += ["--sql-via", "--files", args.files]
                    if args.duckdb_memory_limit is not None:
                        extra += ["--memory-limit", args.duckdb_memory_limit]
                    if args.glob is not None:
                        extra += ["--glob", args.glob.format(size=size)]
                elif args.count_arrow_requests:
                    extra.append("--count-arrow-requests")
                if via == "table":
                    if args.obs_table is None:
                        p.error("the table via needs --obs-table")
                    extra += ["--obs-table", args.obs_table.format(size=size)]
                    if args.x_from_indptr:
                        extra.append("--x-from-indptr")
                    if args.key_first:
                        extra.append("--key-first")
                    if args.scan_workers is not None:
                        extra += ["--scan-workers", str(args.scan_workers)]
                threads = args.duckdb_threads.split(",") if via == "sql" else [""]
                runs = [
                    (int(w), t)
                    for w in args.workers.split(",")
                    for t in (threads if any(threads) else [""])
                ]
                for workers, thread in runs:
                    for repeat in range(args.repeats):
                        more = ["--threads", thread] if thread else []
                        record = one(
                            catalog,
                            shape,
                            workers,
                            args.timeout,
                            cell_type,
                            args.max_rss_gb,
                            extra + more,
                        )
                        record.update(
                            where=args.where, repeat=repeat, size=int(size), via=via
                        )
                        if label:
                            record["label"] = label
                        with args.out.open("a") as fh:
                            fh.write(json.dumps(record) + "\n")
                        seconds = record.get("query_seconds", record["status"])
                        print(
                            f"{size:>5} {name:<22} {via:<6} w={workers:<3} r={repeat} "
                            f"{seconds}",
                            flush=True,
                        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
