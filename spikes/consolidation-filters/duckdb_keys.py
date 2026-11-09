"""The cell keys of each SQL filter, through the DuckDB path of `spike/obs-duckdb`.

Run it with the SDK of that branch installed, not this one: it calls
`open_census(..., obs_engine="duckdb")` and `get_obs`, which only that branch has. The
filters come from `builder_keys.py --export-sql`. Each record holds the number of
cells, the datasets, and the digest of the sorted cell keys, in the form of
`builder_keys.py`.

    python duckdb_keys.py catalog-103.parquet --sql queries-sql.json --out duckdb-keys.jsonl
"""

import argparse
import hashlib
import json
import pathlib
import time

import census2_sdk  # the SDK of spike/obs-duckdb
import numpy as np


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--sql", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--keys-dir", type=pathlib.Path, default=None)
    args = p.parse_args()
    filters = json.loads(args.sql.read_text())
    census = census2_sdk.open_census([args.catalog], obs_engine="duckdb")
    with args.out.open("a") as fh:
        for name, sql in filters.items():
            t0 = time.monotonic()
            frame = census.get_obs(
                "homo_sapiens", sql, columns=["dataset_id", "observation_joinid"]
            )
            seconds = time.monotonic() - t0
            ids = np.asarray(frame["dataset_id"].astype(str)) if len(frame) else []
            joins = (
                np.asarray(frame["observation_joinid"].astype(str))
                if len(frame)
                else []
            )
            keys = sorted(f"{d}\t{j}" for d, j in zip(ids, joins))
            h = hashlib.sha256()
            for key in keys:
                h.update(key.encode())
                h.update(b"\n")
            stats = frame.attrs.get("census2", {})
            record = {
                "query": name,
                "path": "duckdb",
                "cells": len(keys),
                "datasets": len(set(ids)),
                "digest": h.hexdigest(),
                "seconds": round(seconds, 3),
                "candidates": stats.get("catalog_candidates"),
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
