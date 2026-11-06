"""Variants of the obs fragments, written by DuckDB, and a catalog over them.

The writer's fragments come from pyarrow, in store row order, with Parquet dictionary
encoding and no bloom filters. This tool writes a variant of each fragment of a
catalog under another prefix, without a new ingest, so that the benchmark can compare
layouts on the same cells:

- `bloom`: the same rows in the same order, written by DuckDB, which adds a bloom
  filter to each dictionary-encoded column chunk. DuckDB tests an equality or `IN`
  filter against the bloom filter before it reads the column chunk.
- `sorted`: the rows sorted by `--sort` (default `cell_type,tissue`), then by
  `store_row`, with bloom filters. The statistics of a row group then hold a narrow
  range of labels. `store_row`, `x_start`, and `x_end` still give each cell's place in
  its store, so the SDK reads the stores as before; it sorts each dataset's rows by
  `store_row` after the scan.

The variant of a row's fragment is `<out>/<dataset_id>.obs.parquet`, so the copies of
`bench.fragments`, which share a file name, get one variant each. The output catalog is
the input catalog with each row's `obs_fragment_uri` at its variant
(`census2_sdk.fragments.URI_COLUMN`) and the variant's size in `obs_fragment_bytes`.
With `--sizes`, `--catalog-out` is a directory, and it gets `catalog-<n>.parquet`, the
first n rows, for each size, as `bench.fragments` writes them.

    python -m bench.variants catalog-103.parquet --variant sorted \\
        --out s3://my-bucket/variants/sorted --catalog-out catalog-103-sorted.parquet
"""

import argparse
import pathlib
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import duckdb
import pandas as pd

from census2_sdk import catalog as _catalog
from census2_sdk import fragments, obs_sql

VARIANTS = ("bloom", "sorted")


def copy_sql(source: str, dest: str, variant: str, sort: list[str], rows: int) -> str:
    order = ""
    if variant == "sorted":
        keys = ", ".join(obs_sql.identifier(c) for c in [*sort, fragments.STORE_ROW])
        order = f" ORDER BY {keys}"
    return (
        f"COPY (SELECT * FROM read_parquet({obs_sql.quote(source)}){order}) "
        f"TO {obs_sql.quote(dest)} "
        f"(FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE {rows})"
    )


def size_of(uri: str) -> int:
    """The size of one file, local or remote (one HEAD)."""
    if not _catalog.is_remote(uri):
        return pathlib.Path(uri).stat().st_size
    import obstore
    import obstore.store

    base, name = uri.rsplit("/", 1)
    return int(obstore.head(obstore.store.from_url(base), name)["size"])


def build(
    catalog: pd.DataFrame,
    out: str,
    variant: str,
    sort: list[str],
    rows: int,
    workers: int,
) -> pd.DataFrame:
    """Write the variant of each row's fragment; return the catalog over the variants."""
    con = duckdb.connect()
    if _catalog.is_remote(out) or any(
        _catalog.is_remote(str(u)) for u in catalog["store_uri"]
    ):
        con.install_extension("httpfs")
        con.load_extension("httpfs")
        obs_sql.s3_secret(con)
    out = out.rstrip("/")
    if not _catalog.is_remote(out):
        pathlib.Path(out).mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = catalog.to_dict("records")  # type: ignore[assignment]
    sources = [fragments.uri(r) for r in records]
    dests = [f"{out}/{r['dataset_id']}{fragments.SUFFIX}" for r in records]

    def one(pair: tuple[str, str]) -> int:
        cursor = con.cursor()
        cursor.execute(copy_sql(pair[0], pair[1], variant, sort, rows))
        return size_of(pair[1])

    t0 = time.monotonic()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        sizes = list(pool.map(one, zip(sources, dests, strict=True)))
    print(f"{len(dests)} {variant} fragments in {time.monotonic() - t0:.1f} s")
    result = catalog.copy()
    result[fragments.URI_COLUMN] = dests
    result["obs_fragment_bytes"] = sizes
    result["obs_fragment_sha256"] = None
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--variant", choices=VARIANTS, required=True)
    p.add_argument("--out", required=True, help="the prefix of the variants")
    p.add_argument("--catalog-out", type=pathlib.Path, required=True)
    p.add_argument("--sort", default="cell_type,tissue")
    p.add_argument("--row-group-rows", type=int, default=65_536)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--sizes", default=None, help="catalogs of the first n rows")
    args = p.parse_args()
    catalog = pd.read_parquet(args.catalog)
    result = build(
        catalog,
        args.out,
        args.variant,
        args.sort.split(","),
        args.row_group_rows,
        args.workers,
    )
    if args.sizes is None:
        result.to_parquet(args.catalog_out)
        print(f"{len(result)} rows -> {args.catalog_out}")
        return 0
    args.catalog_out.mkdir(parents=True, exist_ok=True)
    for size in (int(n) for n in args.sizes.split(",")):
        path = args.catalog_out / f"catalog-{size}.parquet"
        result.head(size).to_parquet(path)
        print(f"{size} rows -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
