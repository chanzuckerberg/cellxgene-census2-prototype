"""Catalogs of N stores for the store-count benchmark.

The source catalog has one row for each store under a source prefix. The copies
of those stores are under `<copies>/c01/` to `<copies>/cNN/`, with the same store
names. `expand` adds one block of rows for each copy and a `copy` column (0 for
the source rows). Each row of a copy gets a new `dataset_id`, because the SDK
builds `obs_names` from it, and the `store_uri` of the copy. A catalog of N
stores is the first N rows, so a smaller catalog is a prefix of a larger one.

    python -m bench.catalogs catalog-source.parquet \\
        --copies s3://bucket/census2-bench/stores --n-copies 19 \\
        --sizes 50,100,500,1000,2000 --out-dir catalogs/
"""

import argparse
import pathlib

import pandas as pd


def expand(source: pd.DataFrame, copies: str, n_copies: int) -> pd.DataFrame:
    """The source rows, then one block of rows for each copy."""
    frames = [source.assign(copy=0)]
    names = [str(uri).rstrip("/").rsplit("/", 1)[1] for uri in source["store_uri"]]
    for k in range(1, n_copies + 1):
        block = source.assign(copy=k)
        block["dataset_id"] = block["dataset_id"].astype(str) + f"-c{k:02d}"
        block["store_uri"] = [f"{copies.rstrip('/')}/c{k:02d}/{n}" for n in names]
        frames.append(block)
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("source", help="the catalog of the source stores")
    p.add_argument("--copies", required=True, help="the prefix that holds c01/...")
    p.add_argument("--n-copies", type=int, required=True)
    p.add_argument("--sizes", default="50,100,500,1000,2000")
    p.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("."))
    args = p.parse_args()

    rows = expand(pd.read_parquet(args.source), args.copies, args.n_copies)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for size in (int(s) for s in args.sizes.split(",")):
        if size > len(rows):
            p.error(f"size {size} exceeds the {len(rows)} stores")
        out = args.out_dir / f"catalog-{size}.parquet"
        rows.head(size).to_parquet(out)
        print(f"{size} rows -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
