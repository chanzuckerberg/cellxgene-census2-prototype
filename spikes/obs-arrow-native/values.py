"""Compare the `get_anndata` values of a mapping query on `main` and on this branch.

`run` runs each query of `bench/filters.py` that has a mapping form, over the SDK
that this interpreter imports, with the normalized layer, and records a SHA-256 of
each part of the result: the cell keys, each `obs` column as strings, the dtypes,
`var`, `X`, and the normalized layer, with the rows in the order of `obs_names`. The
API is `mapping` (`main`) or `spike` (`Organism.where`). It also records whether
each categorical column holds only the categories in use.

`compare` compares two records part by part. The categories are the one expected
difference: this branch keeps only the categories in use.

    python values.py run --api mapping --catalog catalog-103.parquet \\
        --tools census2-sdk/tools --out values-main.jsonl
    python values.py compare values-main.jsonl values-spike.jsonl
"""

import argparse
import hashlib
import json
import pathlib
import sys
import time
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp


def sha(*arrays: Any) -> str:
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(np.asarray(a)).tobytes())
    return h.hexdigest()


def strings(values: Any) -> str:
    h = hashlib.sha256()
    for v in pd.Series(values).astype(object):
        h.update(("\x00" if pd.isna(v) else str(v)).encode())
        h.update(b"\n")
    return h.hexdigest()


def parts(adata: Any) -> dict[str, Any]:
    order = np.argsort(np.asarray(adata.obs_names))
    out: dict[str, Any] = {
        "n_obs": int(adata.n_obs),
        "keys": strings(adata.obs_names[order]),
    }
    obs = adata.obs.iloc[order]
    out["columns"] = list(obs.columns)
    out["obs"] = {c: strings(obs[c]) for c in obs.columns}
    kinds = {}
    in_use = {}
    for c in obs.columns:
        dtype = obs[c].dtype
        kinds[c] = "category" if isinstance(dtype, pd.CategoricalDtype) else str(dtype)
        if isinstance(dtype, pd.CategoricalDtype):
            used = sorted(str(v) for v in obs[c].dropna().unique())
            in_use[c] = sorted(str(v) for v in dtype.categories) == used
    out["kinds"] = kinds
    out["categories_in_use"] = in_use
    out["var"] = strings(adata.var_names) + strings(list(adata.var.columns))
    for name, matrix in [
        ("X", adata.X),
        ("normalized", adata.layers.get("normalized")),
    ]:
        if matrix is None:
            continue
        m = sp.csr_matrix(matrix)[order]
        m.sort_indices()
        out[name] = sha(m.indptr, m.indices, m.data)
        out[f"{name}_dtype"] = str(m.dtype)
    return out


def run(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(args.tools))
    import census2_sdk
    from bench.filters import QUERIES

    census = census2_sdk.open_census([args.catalog])
    names = args.queries.split(",") if args.queries else list(QUERIES)
    with args.out.open("w") as fh:
        for name in names:
            mapping = QUERIES[name].mapping
            if mapping is None:
                continue
            t0 = time.monotonic()
            if args.api == "mapping":
                adata = census.get_anndata(
                    "homo_sapiens",
                    mapping,
                    layers=["normalized"],
                    max_workers=args.workers,
                )
            else:
                adata = census.get_anndata(
                    census.homo_sapiens.where(mapping),
                    layers=["normalized"],
                    max_workers=args.workers,
                )
            record: dict[str, Any] = {"query": name, "api": args.api}
            record["seconds"] = round(time.monotonic() - t0, 2)
            record.update(parts(adata) if adata is not None else {"n_obs": 0})
            fh.write(json.dumps(record) + "\n")
            fh.flush()
            print(name, record["n_obs"], record["seconds"], flush=True)


def compare(args: argparse.Namespace) -> int:
    def load(path: pathlib.Path) -> dict[str, dict[str, Any]]:
        rows = [json.loads(x) for x in path.read_text().splitlines() if x.strip()]
        return {r["query"]: r for r in rows}

    left, right = load(args.left), load(args.right)
    rows = []
    for name, a in left.items():
        b = right.get(name)
        if b is None:
            continue
        diffs = [
            k
            for k in ("n_obs", "keys", "columns", "var", "X", "normalized", "X_dtype")
            if a.get(k) != b.get(k)
        ]
        diffs += [
            f"obs.{c}"
            for c in a.get("obs", {})
            if a["obs"][c] != b.get("obs", {}).get(c)
        ]
        diffs += [
            f"kind.{c}"
            for c in a.get("kinds", {})
            if a["kinds"][c] != b.get("kinds", {}).get(c)
        ]
        narrowed = all(b.get("categories_in_use", {}).values())
        rows.append(
            {
                "query": name,
                "n_obs": a["n_obs"],
                "differences": diffs,
                "right_in_use": narrowed,
            }
        )
        print(
            f"{name:<16} {a['n_obs']:>9} {'same' if not diffs else diffs} in_use={narrowed}"
        )
    same = sum(not r["differences"] for r in rows)
    print(f"{same} of {len(rows)} queries give the same values")
    if args.out is not None:
        args.out.write_text(json.dumps(rows, indent=1) + "\n")
    return 0 if same == len(rows) else 1


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--api", choices=["mapping", "spike"], required=True)
    r.add_argument("--catalog", required=True)
    r.add_argument("--tools", type=pathlib.Path, required=True)
    r.add_argument("--queries", default=None)
    r.add_argument("--workers", type=int, default=32)
    r.add_argument("--out", type=pathlib.Path, required=True)
    c = sub.add_parser("compare")
    c.add_argument("left", type=pathlib.Path)
    c.add_argument("right", type=pathlib.Path)
    c.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    if args.command == "run":
        run(args)
        return 0
    return compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
