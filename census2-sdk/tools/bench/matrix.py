"""Run the benchmark matrix: each combination in a new process.

The matrix is sizes x shapes x workers x repeats. A run that exceeds `--timeout`
gets the status `exceeded`. A run that passes `--max-rss-gb` gets the status
`memory_limit` from the run itself. A run that fails in any other way gets the
status `error` and the last lines of its output. Records append to `--out` as JSON lines,
so a stopped matrix keeps its finished runs.

    python -m bench.matrix --catalogs catalogs/ --sizes 50,100 \\
        --shapes single,cell_type,open,corpus_obs --workers 1,32 --repeats 3 \\
        --where laptop --out results/bench-laptop.jsonl
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
) -> dict[str, Any]:
    command = [
        sys.executable, "-m", "bench.run", str(catalog),
        "--shape", shape, "--workers", str(workers), "--cell-type", cell_type,
    ]  # fmt: skip
    if max_rss_gb is not None:
        command += ["--max-rss-gb", str(max_rss_gb)]
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
    args = p.parse_args()

    # One cell type for every size, so that each size runs the same query.
    sizes = args.sizes.split(",")
    largest = args.catalogs / f"catalog-{max(sizes, key=int)}.parquet"
    cell_type = args.cell_type or pick_cell_type(source_rows(pd.read_parquet(largest)))
    print(f"cell type: {cell_type}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for size in sizes:
        catalog = args.catalogs / f"catalog-{size}.parquet"
        for shape in args.shapes.split(","):
            for workers in (int(w) for w in args.workers.split(",")):
                for repeat in range(args.repeats):
                    record = one(
                        catalog,
                        shape,
                        workers,
                        args.timeout,
                        cell_type,
                        args.max_rss_gb,
                    )
                    record.update(where=args.where, repeat=repeat, size=int(size))
                    with args.out.open("a") as fh:
                        fh.write(json.dumps(record) + "\n")
                    seconds = record.get("query_seconds", record["status"])
                    print(f"{size:>5} {shape:<10} w={workers:<3} r={repeat} {seconds}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
