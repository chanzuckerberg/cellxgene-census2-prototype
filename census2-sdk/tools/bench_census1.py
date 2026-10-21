"""Time the benchmark questions against Census 1, for comparison (task 6.4).

The queries match the shapes of `bench.run` on the Census 2 stores:

- `single`: `get_anndata` for the dataset that `bench.run --shape single` reads.
- `cell_type`: `get_anndata` for one cell type in the source datasets of the
  benchmark. Census 1 holds the whole corpus. The dataset list scopes the query
  to the same datasets as a Census 2 catalog of 103 stores.
- `obs`: `get_obs` for all columns of the same source datasets.
- `cell_type_count`: the cells of the cell type in the whole human corpus, read
  with one column. This gives the scale of the unscoped question.

Each query runs in a new process and prints one JSON record. Run the matrix in
the oracle environment of the harness (`harness/.venv-oracle`):

    python bench_census1.py --source catalog-source.parquet --single-dataset ID \\
        --cell-type neutrophil --repeats 3 --out results/census1.jsonl

The source catalog is a Parquet file, and this script reads only its
`dataset_id` column.
"""

import argparse
import json
import pathlib
import platform
import resource
import subprocess
import sys
import time
from typing import Any

CENSUS_VERSION = "2025-11-08"
ORGANISM = "Homo sapiens"
QUERIES = ("single", "cell_type", "obs", "cell_type_count")


def peak_rss_mb() -> float:
    """The peak resident memory of this process. macOS reports bytes, Linux KiB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(peak / (2**20 if sys.platform == "darwin" else 2**10), 1)


def run_one(query: str, dataset: str, cell_type: str, ids: list[str]) -> dict[str, Any]:
    import cellxgene_census

    t0 = time.monotonic()
    census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
    t1 = time.monotonic()
    scope = f"dataset_id in {ids}"
    try:
        if query == "single":
            adata = cellxgene_census.get_anndata(
                census, ORGANISM, obs_value_filter=f"dataset_id == '{dataset}'"
            )
            cells = adata.n_obs
        elif query == "cell_type":
            adata = cellxgene_census.get_anndata(
                census,
                ORGANISM,
                obs_value_filter=f"cell_type == '{cell_type}' and {scope}",
            )
            cells = adata.n_obs
        elif query == "obs":
            cells = len(cellxgene_census.get_obs(census, ORGANISM, value_filter=scope))
        else:
            obs = cellxgene_census.get_obs(
                census,
                ORGANISM,
                value_filter=f"cell_type == '{cell_type}'",
                column_names=["soma_joinid"],
            )
            cells = len(obs)
        t2 = time.monotonic()
    finally:
        census.close()
    return {
        "stack": "census1",
        "census_version": CENSUS_VERSION,
        "query": query,
        "datasets_in_scope": {"single": 1, "cell_type_count": "release"}.get(
            query, len(ids)
        ),
        "cells": int(cells),
        "open_seconds": round(t1 - t0, 3),
        "query_seconds": round(t2 - t1, 3),
        "peak_rss_mb": peak_rss_mb(),
        "host": platform.node(),
        "status": "ok",
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=pathlib.Path, required=True)
    p.add_argument("--single-dataset", required=True)
    p.add_argument("--cell-type", required=True)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--timeout", type=float, default=1800)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--one", choices=QUERIES, help="run one query in this process")
    args = p.parse_args()

    import pandas as pd

    ids = sorted(
        str(d)
        for d in pd.read_parquet(args.source, columns=["dataset_id"])["dataset_id"]
    )
    if args.one:
        print(json.dumps(run_one(args.one, args.single_dataset, args.cell_type, ids)))
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    for query in QUERIES:
        for repeat in range(args.repeats):
            command = [sys.executable, __file__, *sys.argv[1:], "--one", query]
            base: dict[str, Any] = {"stack": "census1", "query": query}
            try:
                done = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                    check=False,
                )
                if done.returncode == 0:
                    record = dict(json.loads(done.stdout.strip().splitlines()[-1]))
                else:
                    output = f"{done.stdout}\n{done.stderr}"
                    tail = output.strip().splitlines()[-10:]
                    record = {**base, "status": "error", "error": tail}
            except subprocess.TimeoutExpired as timeout:
                # TimeoutExpired holds bytes even with text=True.
                parts = (timeout.stdout, timeout.stderr)
                output = "\n".join(
                    p.decode(errors="replace") if isinstance(p, bytes) else (p or "")
                    for p in parts
                )
                record = {
                    **base,
                    "status": "exceeded",
                    "timeout_seconds": args.timeout,
                    "output": output.strip().splitlines()[-10:],
                }
            record["repeat"] = repeat
            with args.out.open("a") as fh:
                fh.write(json.dumps(record) + "\n")
            print(
                f"{query:<16} r={repeat} {record.get('query_seconds', record['status'])}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
