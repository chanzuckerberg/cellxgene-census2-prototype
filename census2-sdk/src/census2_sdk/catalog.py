"""The SDK's side of the catalog: load it, and exclude the datasets that cannot match.

The catalog is a Parquet table with one row for each store (design section 2). The
Census 2 team publishes it; `tools/build_catalog.py` is the Iteration 1 stand-in for
that builder. `exclude` applies the one correctness rule of design 2.2. A query
excludes a dataset only when the catalog proves that no cell of the dataset can match.
An absent column or an unknown value keeps the dataset.
"""

import io
import os
from collections.abc import Iterable, Mapping
from typing import Any, cast

import pandas as pd
import zarr

# obstore is the optional `s3` extra. It is imported here, at module level, when it
# is installed, so that a local catalog over local stores installs and runs
# without it. Each remote path calls `_require_obstore` first.
try:
    import obstore
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


def load(uri: str | os.PathLike[str]) -> pd.DataFrame:
    if is_remote(uri):
        _require_obstore()
        base, _, name = str(uri).rpartition("/")
        data = bytes(obstore.get(obstore.store.from_url(base), name).bytes())
        return pd.read_parquet(io.BytesIO(data))
    return pd.read_parquet(uri)


def excludes(row: Mapping[str, Any], obs_filter: Normalized) -> bool:
    """True when the row proves that no cell of the dataset matches `obs_filter`."""
    for column, wanted in obs_filter.items():
        if column == "dataset_id":
            if str(row["dataset_id"]) not in wanted:
                return True
            continue
        values = row[column]
        if values is None or (not isinstance(values, list) and pd.isna(values).all()):
            continue  # the store has no such column; unproven keeps the dataset
        present = {m for v in values for m in members(str(v), column)}
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
    return scoped[keep]
