"""The `obs` fragment next to each store: its location and the columns that it adds.

The writer puts one Parquet file next to each store, with one row for each cell: the
store's `obs` columns, then `store_row`, `x_start`, and `x_end` (int64). The fragment
of `<name>.zarr` is `<name>.obs.parquet`, in the same directory. These names agree with
`census2_dataset_writer.store_format`, which the SDK does not import.

Both `obs` query paths read the fragments: the Iceberg table (`census2_sdk.obs_table`)
and DuckDB over the files (`census2_sdk.obs_sql`). This module imports neither engine.
"""

from collections.abc import Mapping
from typing import Any

from census2_sdk.catalog import unrecorded

STORE_ROW = "store_row"
X_START = "x_start"
X_END = "x_end"
ROW_COLUMNS = (STORE_ROW, X_START, X_END)
SUFFIX = ".obs.parquet"

# A catalog row can name its fragment, for example a copy of a fragment that the
# benchmark makes (`bench.fragments`). Without it, the fragment is next to the store.
URI_COLUMN = "obs_fragment_uri"


def location(store_uri: str) -> str:
    """The fragment next to the store at `store_uri`, whose name ends with `.zarr`."""
    base = store_uri.rstrip("/")
    if not base.endswith(".zarr"):
        raise ValueError(f"a store name ends with .zarr: {store_uri!r}")
    return base.removesuffix(".zarr") + SUFFIX


def uri(row: Mapping[str, Any]) -> str:
    """The fragment of one catalog row: the row's own `obs_fragment_uri`, or the
    location next to its store."""
    named = row.get(URI_COLUMN)
    if isinstance(named, str) and not unrecorded(named):
        return named
    return location(str(row["store_uri"]))
