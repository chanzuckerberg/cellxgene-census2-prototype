"""Copies of the obs fragments, and catalogs of N datasets over them.

The obs table grows with the number of datasets, and a sample of 103 datasets is far
below the corpus. A fragment is small, so copies of the source fragments are a cheap
way to measure how the plan and the scan of the table grow with N. Copy k of a fragment
holds the same rows with the `dataset_id` `<dataset_id>-c<k>`, at
`<out>/c<k>/<name>.obs.parquet`. Its catalog row is the source row with that
`dataset_id`, the copy's fragment in `obs_fragment_uri`, and the copy's size and digest.
The row keeps the source `store_uri`, so a query through the table reads the source
store for each copy. A catalog of N datasets is the first N rows, as in `catalogs`.

    python -m bench.fragments catalog-source.parquet --out s3://bucket/scale/fragments \\
        --n-copies 19 --sizes 103,500,1000,2000 --out-dir catalogs/
"""

import argparse
import pathlib
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import census2_dataset_writer
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from census2_sdk import obs_table


def copy_one(source: str, dest: str, dataset_id: str) -> dict[str, Any]:
    """Write the fragment at `source` to `dest` with `dataset_id`. Return its record."""
    table = pq.read_table(pa.BufferReader(obs_table.read_bytes(source)))
    column = pa.array([dataset_id] * table.num_rows, type=pa.string())
    table = table.set_column(0, census2_dataset_writer.DATASET_ID_COLUMN, column)
    return dict(census2_dataset_writer.write_fragment(table, dest))


def expand(
    source: pd.DataFrame, out: str, n_copies: int, workers: int = 16
) -> pd.DataFrame:
    """The source rows, then one block of rows for each copy, with the copies written."""
    out = out.rstrip("/")
    fragments = [
        census2_dataset_writer.fragment_location(str(s)) for s in source["store_uri"]
    ]
    frames = [source.assign(copy=0, obs_fragment_uri=fragments)]
    for k in range(1, n_copies + 1):
        block = source.assign(copy=k)
        block["dataset_id"] = block["dataset_id"].astype(str) + f"-c{k:02d}"
        names = [f.rsplit("/", 1)[1] for f in fragments]
        block["obs_fragment_uri"] = [f"{out}/c{k:02d}/{n}" for n in names]
        jobs = zip(
            fragments, block["obs_fragment_uri"], block["dataset_id"], strict=True
        )
        with ThreadPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(lambda job: copy_one(*job), jobs))
        block["obs_fragment_bytes"] = [r["bytes"] for r in records]
        block["obs_fragment_sha256"] = [r["sha256"] for r in records]
        frames.append(block)
        print(f"copy {k}: {len(block)} fragments")
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("source", help="the catalog of the source stores")
    p.add_argument("--out", required=True, help="the prefix of the copies")
    p.add_argument("--n-copies", type=int, required=True)
    p.add_argument("--sizes", default="103,500,1000,2000")
    p.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("."))
    args = p.parse_args()

    rows = expand(pd.read_parquet(args.source), args.out, args.n_copies)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for size in (int(s) for s in args.sizes.split(",")):
        if size > len(rows):
            p.error(f"size {size} exceeds the {len(rows)} datasets")
        path = args.out_dir / f"catalog-{size}.parquet"
        rows.head(size).to_parquet(path)
        print(f"{size} rows -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
