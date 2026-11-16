"""Time the read and the evaluation of a query on each store, three ways.

For each store of a catalog, the tool opens the store once and then measures, one
at a time:

- `read_elem`: `anndata.io.read_elem(obs)`, the whole `obs` of the store. It runs
  once for each repeat.
- for each query, `read`: `obs_read.read_columns` of the query's columns, the codes
  of a categorical column with the catalog's categories;
- `arrow`: `obs_read.evaluate`, the query's pyarrow expression over the
  `DictionaryArray`s, projected with Acero to one boolean for each cell, in one
  thread. This is what the SDK does.
- `filter`: the same expression through `Table.filter`, over the same columns and a
  column of row positions.
- `numpy`: the evaluation of the earlier filter spike (`spike/consolidation-filters`)
  on the same arrays: each label term is a lookup table over the categories, mapped
  by the codes with NumPy, and pyarrow evaluates only booleans and numbers. The tool
  imports that branch's `obs_filter`, `obs_columns`, and `obs_arrow` from
  `--pr43`, a checkout of its `census2-sdk/src/census2_sdk`, as a package of
  another name.

The query is the whole query, not the catalog's view of it, so every store does the
full work. Each record holds the wall and CPU milliseconds of each measure, and the
requests and bytes of each read, at the `zarr.storage.ObjectStore` interface. The
three evaluations must keep the same cells.

    python bind.py catalog-103.parquet --tools census2-sdk/tools \\
        --pr43 ../pr43/census2-sdk/src/census2_sdk --queries neutrophil,wmg \\
        --repeats 3 --where ec2 --out bind.jsonl
"""

import argparse
import functools
import importlib
import json
import pathlib
import re
import shutil
import sys
import tempfile
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import pyarrow as pa
import zarr
import zarr.storage
from anndata.io import read_elem

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from e2e import COUNTER, CountingObjectStore, usage


def measure(fn: Callable[[], Any]) -> tuple[Any, dict[str, float]]:
    requests, nbytes = COUNTER.requests, COUNTER.nbytes
    user, system = usage()
    t0 = time.perf_counter()
    value = fn()
    wall = time.perf_counter() - t0
    end_user, end_system = usage()
    return value, {
        "ms": round(wall * 1e3, 3),
        "cpu_ms": round((end_user - user + end_system - system) * 1e3, 3),
        "requests": COUNTER.requests - requests,
        "mb": round((COUNTER.nbytes - nbytes) / 2**20, 3),
    }


def vendor(src: pathlib.Path) -> Any:
    """The earlier spike's builder modules, as the package `c43`."""
    root = pathlib.Path(tempfile.mkdtemp())
    package = root / "c43"
    package.mkdir()
    (package / "__init__.py").write_text("")
    for name in ("obs_filter.py", "obs_arrow.py", "obs_columns.py"):
        text = (src / name).read_text()
        (package / name).write_text(re.sub(r"\bcensus2_sdk\.", "c43.", text))
    sys.path.insert(0, str(root))
    modules = {
        n: importlib.import_module(f"c43.{n}")
        for n in ("obs_filter", "obs_arrow", "obs_columns")
    }
    shutil.rmtree(root, ignore_errors=True)
    return modules


def filter_positions(expression: Any, data: dict[str, Any], n: int) -> np.ndarray:
    """The kept positions through `Table.filter`, with a column of row positions."""
    from census2_sdk import obs_read

    table = obs_read.table(data, n).append_column("__pos", pa.array(np.arange(n)))
    return np.asarray(table.filter(expression)["__pos"].to_numpy())


def numpy_positions(
    c43: dict[str, Any], tree: Any, data: dict[str, Any], n: int
) -> np.ndarray:
    from census2_sdk import obs_read

    arrow = c43["obs_arrow"]
    converted = {}
    for name, d in data.items():
        if isinstance(d, obs_read.Codes):
            categories = np.asarray(list(d.categories), dtype=object)
            converted[name] = arrow.Codes(d.codes, categories, d.ordered)
        else:
            converted[name] = arrow.Values(d.values, d.mask)
    return np.asarray(arrow.evaluate(tree, converted, n))


def bench_store(
    row: dict[str, Any],
    queries: list[Any],
    c43: dict[str, Any],
    clients: Any,
    repeat: int,
) -> list[dict[str, Any]]:
    import census2_sdk
    from census2_sdk import obs_read

    uri = str(row["store_uri"])
    group = census2_sdk.catalog.open_store(uri, clients)
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)
    n = census2_sdk._n_obs(group)
    listed = census2_sdk._listed(row)
    _, main = measure(lambda: read_elem(obs))
    hs = census2_sdk.Organism("homo_sapiens")
    obs43 = c43["obs_columns"].Obs("homo_sapiens")
    out = []
    for query in queries:
        q = query.build(hs.obs)
        names = sorted(q.columns)
        data, read = measure(
            functools.partial(obs_read.read_columns, obs, names, listed)
        )
        e = q.expression
        kept, arrow = measure(functools.partial(obs_read.evaluate, e, data, n))
        by_filter, filtered = measure(functools.partial(filter_positions, e, data, n))
        f43 = c43["obs_filter"]
        tree = f43.normalize(f43.as_expr(query.build(obs43)))
        by_numpy, numpy = measure(
            functools.partial(numpy_positions, c43, tree, data, n)
        )
        same = np.array_equal(kept, by_filter) and np.array_equal(kept, by_numpy)
        out.append(
            {
                "store": uri.rstrip("/").rsplit("/", 1)[-1],
                "n_obs": n,
                "query": query.name,
                "columns": names,
                "repeat": repeat,
                "kept": len(kept),
                "same": bool(same),
                "read_elem": main,
                "read": read,
                "arrow": arrow,
                "filter": filtered,
                "numpy": numpy,
            }
        )
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--tools", type=pathlib.Path, required=True)
    p.add_argument("--pr43", type=pathlib.Path, required=True)
    p.add_argument("--queries", required=True)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--where", required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    zarr.storage.ObjectStore = CountingObjectStore  # type: ignore[misc]
    pa.set_cpu_count(1)
    sys.path.insert(0, str(args.tools))
    import census2_sdk
    from bench.filters import QUERIES

    queries = [QUERIES[name] for name in args.queries.split(",")]
    c43 = vendor(args.pr43)
    census = census2_sdk.open_census([args.catalog])
    rows = census._candidates(census.homo_sapiens).to_dict("records")[: args.limit]
    clients = census2_sdk.catalog.clients(str(r["store_uri"]) for r in rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("a") as fh:
        for repeat in range(args.repeats):
            for row in rows:
                for record in bench_store(row, queries, c43, clients, repeat):
                    record["where"] = args.where
                    fh.write(json.dumps(record) + "\n")
                fh.flush()
            print(f"repeat {repeat}: {len(rows)} stores", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
