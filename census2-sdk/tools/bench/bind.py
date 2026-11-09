"""Time the bind and evaluation of a filter on each store, against a full `obs` read.

For each store that a catalog names (the consolidated stores in place of their
members, unless `--per-dataset`), the tool opens the store once with the lean open,
and then, one measure at a time, in this order for each query:

- `main`: `anndata.io.read_elem(obs)`, the whole `obs` of the store, as `main`'s
  `get_anndata` reads it before its filter.
- `bind`: the builder's store side for the query over all rows of the store:
  `obs_arrow.read_columns` of the filter's columns (`read`), then
  `obs_arrow.evaluate` (`eval`), the bind and the evaluation in Arrow.

The filter is the whole normalized filter of the query, not the residual of the
catalog test, so every store does the full work. Each record holds the wall and the
CPU milliseconds of each measure (the CPU of the process, which includes zarr's
event-loop thread), and its requests and payload at the zarr store interface. Each
record is one JSON line.

    python -m bench.bind catalog.parquet --queries neutrophil,wmg --repeats 3 \\
        --where ec2 --out bind.jsonl
"""

import argparse
import functools
import json
import pathlib
import time
from typing import Any, cast

import pandas as pd
import zarr
from anndata.io import read_elem

import census2_sdk
from bench import run
from bench.filters import QUERIES
from census2_sdk import obs_arrow, obs_columns
from census2_sdk.obs_filter import as_expr, columns, normalize


def measure(fn: Any) -> tuple[Any, dict[str, float]]:
    """Run `fn()`; return its value, and its milliseconds, requests, and MB."""
    requests, nbytes = run.COUNTER.requests, run.COUNTER.nbytes
    user, system = run.cpu_seconds()
    t0 = time.perf_counter()
    value = fn()
    wall = time.perf_counter() - t0
    end_user, end_system = run.cpu_seconds()
    return value, {
        "ms": round(wall * 1e3, 3),
        "cpu_ms": round((end_user - user + end_system - system) * 1e3, 3),
        "requests": run.COUNTER.requests - requests,
        "mb": round((run.COUNTER.nbytes - nbytes) / 2**20, 3),
    }


def bench_store(uri: str, names: list[str], repeat: int) -> list[dict[str, Any]]:
    group = run.open_lean(uri)
    obs = group["obs"]
    assert isinstance(obs, zarr.Group), "obs must be an encoded dataframe"
    n_obs = int(cast("list[int]", group["X"].attrs["shape"])[0])
    record = cast("dict[str, Any]", group.attrs.get("cellxgene-census2", {}))
    layout = "consolidated" if "consolidation" in record else "per_dataset"
    namespace = obs_columns.Obs("homo_sapiens")
    out = []
    for name in names:
        query = QUERIES[name]
        built = query.build(namespace)
        if built is None:
            continue
        tree = normalize(as_expr(built))
        cols = sorted(columns(tree) | set(query.by))
        frame, main = measure(functools.partial(read_elem, obs))
        cache = obs_arrow.CategoryCache(obs)
        read_all = functools.partial(obs_arrow.read_columns, obs, cols, 0, n_obs, cache)
        data, read = measure(read_all)
        kept, evaluate = measure(
            functools.partial(obs_arrow.evaluate, tree, data, n_obs)
        )
        out.append(
            {
                "store": uri.rstrip("/").rsplit("/", 1)[-1],
                "layout": layout,
                "n_obs": n_obs,
                "obs_columns": len(cast("pd.DataFrame", frame).columns),
                "query": name,
                "columns": cols,
                "repeat": repeat,
                "kept": len(kept),
                "main": main,
                "read": read,
                "eval": evaluate,
                "bind_ms": round(read["ms"] + evaluate["ms"], 3),
                "bind_cpu_ms": round(read["cpu_ms"] + evaluate["cpu_ms"], 3),
            }
        )
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--queries", required=True, help="names in bench.filters.QUERIES")
    p.add_argument("--per-dataset", action="store_true")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--limit", type=int, default=None, help="the first N stores only")
    p.add_argument("--where", required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    names = args.queries.split(",")
    unknown = sorted(set(names) - set(QUERIES))
    if unknown:
        p.error(f"no query {unknown}")
    census = census2_sdk.open_census([args.catalog])
    rows = census.catalog[census.catalog["organism"] == "homo_sapiens"]
    uris = run.serving_uris(rows, args.per_dataset)[: args.limit]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("a") as fh:
        for repeat in range(args.repeats):
            for uri in uris:
                for record in bench_store(uri, names, repeat):
                    record["where"] = args.where
                    fh.write(json.dumps(record) + "\n")
                fh.flush()
            print(f"repeat {repeat}: {len(uris)} stores", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
