"""The cell keys of each query of `bench.filters`, through the filter builder.

For each query, the script selects the cells with `census.homo_sapiens.get_obs` and
records the number of cells, the datasets, and a digest of the sorted cell keys
`(dataset_id, observation_joinid)`. `duckdb_keys.py` records the same for the SQL path
of the `spike/obs-duckdb` branch, and `compare_keys.py` compares the two.

Run it from `census2-sdk/tools`, with that directory on `PYTHONPATH`, so that
`bench.filters` imports:

    PYTHONPATH=. python ../../spikes/consolidation-filters/builder_keys.py catalog-103.parquet \\
        --out builder-keys.jsonl --export-sql queries-sql.json
"""

import argparse
import hashlib
import json
import pathlib
import time
import warnings
from typing import Any

import census2_sdk
import numpy as np
import pandas as pd
from bench.filters import QUERIES


def digest(dataset_ids: Any, joinids: Any) -> tuple[str, list[str]]:
    """SHA-256 of the sorted keys, one `dataset_id<TAB>observation_joinid` line each."""
    keys = sorted(
        f"{d}\t{j}" for d, j in zip(np.asarray(dataset_ids), np.asarray(joinids))
    )
    h = hashlib.sha256()
    for key in keys:
        h.update(key.encode())
        h.update(b"\n")
    return h.hexdigest(), keys


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--queries", default=None, help="comma-separated; default all")
    p.add_argument("--per-dataset", action="store_true")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--keys-dir", type=pathlib.Path, default=None)
    p.add_argument("--export-sql", type=pathlib.Path, default=None)
    args = p.parse_args()
    names = args.queries.split(",") if args.queries else list(QUERIES)
    names = [n for n in names if QUERIES[n].sql is not None]
    if args.export_sql is not None:
        args.export_sql.write_text(
            json.dumps({n: QUERIES[n].sql for n in names}, indent=1) + "\n"
        )
    census = census2_sdk.open_census([args.catalog], consolidated=not args.per_dataset)
    hs = census.homo_sapiens
    with args.out.open("a") as fh:
        for name in names:
            query = QUERIES[name]
            t0 = time.monotonic()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                frame = hs.get_obs(
                    query.build(hs.obs),
                    columns=["dataset_id", "observation_joinid"],
                    max_workers=args.workers,
                )
            seconds = time.monotonic() - t0
            if len(frame):
                ids, joins = (
                    frame["dataset_id"].astype(str),
                    frame["observation_joinid"],
                )
            else:
                ids, joins = pd.Series([], dtype=str), pd.Series([], dtype=str)
            value, keys = digest(ids, joins.astype(str))
            record = {
                "query": name,
                "path": "builder",
                "cells": len(keys),
                "datasets": int(pd.Series(ids).nunique()),
                "digest": value,
                "seconds": round(seconds, 3),
                "candidates": hs.last_query.get("candidates"),
                "stores_opened": hs.last_query.get("stores_opened"),
            }
            fh.write(json.dumps(record) + "\n")
            fh.flush()
            print(json.dumps(record), flush=True)
            if args.keys_dir is not None:
                args.keys_dir.mkdir(parents=True, exist_ok=True)
                (args.keys_dir / f"{name}.txt").write_text("\n".join(keys) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
