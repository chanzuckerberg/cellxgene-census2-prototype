"""The SDK's side of the catalog: load it, and exclude the datasets that cannot match.

The catalog is a Parquet table with one row for each store (design section 2). The
Census 2 team publishes it; `tools/build_catalog.py` is the Iteration 1 stand-in for
that builder. `exclude` applies the one correctness rule of design 2.2. A query
excludes a dataset only when the catalog proves that no cell of the dataset can match.
An absent column or an unknown value keeps the dataset.

Next to each catalog, the builder writes a cell count table (tracker 5.9). It has one
row for each dataset and each combination of the `COUNTED` values and
`is_primary_data`, with the number of cells in `n_cells`. `cell_counts_uri` gives its
location.
"""

import io
import os
from collections.abc import Iterable, Mapping
from typing import Any, cast

import numpy as np
import pandas as pd
import zarr

# obstore is the optional `s3` extra. It is imported here, at module level, when it
# is installed, so that a local catalog over local stores installs and runs
# without it. Each remote path calls `_require_obstore` first.
try:
    import obstore
    import obstore.exceptions
    import obstore.store
except ImportError:  # pragma: no cover
    # mypy strict; this repository does not use ty, so the directive is mypy's.
    obstore = None  # type: ignore[assignment]


def _require_obstore() -> None:
    if obstore is None:
        raise ImportError("a remote URI needs obstore: pip install 'census2-sdk[s3]'")


Filter = Mapping[str, str | Iterable[str]]
# A filter after `normalize`: one set of accepted values for each column.
Normalized = dict[str, frozenset[str]]

# Design 2.2: the 18 indexed columns.
PAIRS = (
    "assay",
    "cell_type",
    "development_stage",
    "disease",
    "self_reported_ethnicity",
    "sex",
    "tissue",
    "tissue_general",
)
INDEXED = tuple(c for p in PAIRS for c in (p, f"{p}_ontology_term_id")) + (
    "suspension_type",
    "tissue_type",
)
# The indexed columns that the cell count table holds (tracker 5.9).
COUNTED = tuple(
    c
    for p in ("assay", "cell_type", "disease", "sex", "tissue")
    for c in (p, f"{p}_ontology_term_id")
)
COUNTS_SUFFIX = ".cell_counts.parquet"
# Schema 2.4.0 made these two columns multi-valued. Equality is membership.
MULTIVALUED = frozenset({"disease", "disease_ontology_term_id"})
MULTIVALUE_SEP = " || "

REMOTE_SCHEMES = ("s3://", "gs://", "az://", "abfs://", "http://", "https://")


def is_remote(uri: str | os.PathLike[str]) -> bool:
    return str(uri).startswith(REMOTE_SCHEMES)


def open_store(uri: str | os.PathLike[str]) -> zarr.Group:
    """Open a store for reading. A remote URI goes through obstore."""
    if is_remote(uri):
        _require_obstore()
        store = zarr.storage.ObjectStore(
            obstore.store.from_url(str(uri)), read_only=True
        )
        return zarr.open_group(store, mode="r")
    return zarr.open_group(str(uri), mode="r")


def members(value: str, column: str) -> list[str]:
    """The values that one cell holds in `column`, one for a single-valued column."""
    if column in MULTIVALUED:
        return [m.strip() for m in value.split(MULTIVALUE_SEP)]
    return [value]


def normalize(obs_filter: Filter | None) -> Normalized:
    """One set of accepted values for each column. Refuse a column that is not indexed."""
    normalized: Normalized = {}
    for column, value in (obs_filter or {}).items():
        if column != "dataset_id" and column not in INDEXED:
            raise ValueError(
                f"{column!r} is not a catalog column; use one of {INDEXED}"
            )
        wanted = {value} if isinstance(value, str) else {str(v) for v in value}
        normalized[column] = frozenset(wanted)
    return normalized


def cell_counts_uri(catalog_uri: str | os.PathLike[str]) -> str:
    """The cell count table that the catalog builder writes next to a catalog."""
    return str(catalog_uri).removesuffix(".parquet") + COUNTS_SUFFIX


def load(uri: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a Parquet table. A missing table raises FileNotFoundError, local or remote."""
    if is_remote(uri):
        _require_obstore()
        base, _, name = str(uri).rpartition("/")
        try:
            data = bytes(obstore.get(obstore.store.from_url(base), name).bytes())
        except obstore.exceptions.NotFoundError as e:
            raise FileNotFoundError(str(uri)) from e
        return pd.read_parquet(io.BytesIO(data))
    return pd.read_parquet(uri)


def unrecorded(value: Any) -> bool:
    """True when a row records no value: the column is absent, None, or a scalar null.

    A value set (a list or an array) is recorded even when it is empty. A catalog
    without the column gives `None` through `row.get`; two sources where only one
    has the column give a float `NaN` through `pd.concat`.
    """
    return value is None or (np.ndim(value) == 0 and bool(pd.isna(value)))


def excludes(row: Mapping[str, Any], obs_filter: Normalized) -> bool:
    """True when the row proves that no cell of the dataset matches `obs_filter`.

    The proof is the design 2.2 rule. The row records the value set of each
    indexed column. The row excludes the dataset only when the set for some
    filtered indexed column is present and disjoint from the filter values. A
    column that the row does not record proves nothing, so the dataset stays
    in. For `dataset_id`, an identity and not a value set, the row excludes the
    dataset when its identity is not among the filter values. An unproven case
    costs one extra store opened; a wrong exclusion drops cells.
    """
    for column, wanted in obs_filter.items():
        if column == "dataset_id":
            if str(row["dataset_id"]) not in wanted:
                return True
            continue
        values = row.get(column)
        if unrecorded(values):
            continue  # the row records no value set; unproven keeps the dataset
        value_set = cast("Iterable[Any]", values)  # unrecorded() rules out None
        present = {m for v in value_set for m in members(str(v), column)}
        if not (present & wanted):
            return True
    return False


def exclude(
    catalog: pd.DataFrame, organism: str, obs_filter: Filter | None
) -> pd.DataFrame:
    """The catalog rows that can hold a matching cell."""
    scoped = catalog[catalog["organism"] == organism]
    normalized = normalize(obs_filter)
    if not normalized:
        return scoped
    records = cast("list[dict[str, Any]]", scoped.to_dict("records"))
    keep = [not excludes(row, normalized) for row in records]
    # A boolean array: an empty Python list would select no columns.
    return scoped[np.array(keep, dtype=bool)]
