"""The SDK's side of the catalog: load it, and exclude the datasets that cannot match.

The catalog is a Parquet table with one row for each store (design section 2). The
Census 2 team publishes it; `tools/build_catalog.py` is the prototype stand-in for
that builder. `exclude` applies the one correctness rule of design 2.2. A query
excludes a dataset only when the catalog proves that no cell of the dataset can match.
An absent column or an unknown value keeps the dataset.

Each catalog row holds the cell counts of its store in the `cell_counts` column. The
column has one entry for each combination of the `COUNTED` values and `is_primary_data`,
with the number of cells in `n_cells`. `load` reads the catalog without this column, and
`load_counts` reads only the entries.
"""

import io
import os
import pathlib
from collections.abc import Iterable, Mapping
from typing import Any, cast

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]  # pyarrow has no types
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
# The catalog column of cell counts and the columns that it counts.
COUNTS = "cell_counts"
COUNTED = tuple(
    c
    for p in ("assay", "cell_type", "disease", "sex", "tissue", "tissue_general")
    for c in (p, f"{p}_ontology_term_id")
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


Source = str | os.PathLike[str] | io.BytesIO


def read(uri: str | os.PathLike[str]) -> io.BytesIO:
    """The content of a catalog. `load` and `load_counts` of it see one version."""
    if is_remote(uri):
        _require_obstore()
        base, _, name = str(uri).rpartition("/")
        return io.BytesIO(
            bytes(obstore.get(obstore.store.from_url(base), name).bytes())
        )
    return io.BytesIO(pathlib.Path(uri).read_bytes())


def _parquet(source: Source) -> tuple[io.BytesIO, list[str]]:
    """The content of `source`, at its start, and its column names."""
    content = source if isinstance(source, io.BytesIO) else read(source)
    content.seek(0)
    names = [str(n) for n in pq.read_schema(content).names]
    content.seek(0)
    return content, names


def load(source: Source) -> pd.DataFrame:
    """Read a catalog without its `cell_counts` column.

    pandas turns that column into Python dicts, which multiplies the read time of a
    catalog by about 30. Only `count_cells` needs it, through `load_counts`.
    """
    content, names = _parquet(source)
    return pd.read_parquet(content, columns=[c for c in names if c != COUNTS])


def load_counts(source: Source) -> pd.DataFrame:
    """The cell count entries of a catalog, one row each, with `dataset_version_id`.

    A catalog without the `cell_counts` column raises ValueError.
    """
    content, names = _parquet(source)
    if COUNTS not in names:
        raise ValueError(
            f"the catalog has no {COUNTS} column; "
            "build it again with tools/build_catalog.py"
        )
    frame = pd.read_parquet(
        content, columns=["dataset_version_id", COUNTS], dtype_backend="pyarrow"
    )
    # The Arrow accessors flatten the column without one Python object per entry.
    entries = frame[COUNTS]
    counts = cast("pd.DataFrame", entries.list.flatten().struct.explode())
    owners = frame["dataset_version_id"].to_numpy()
    lengths = entries.list.len().to_numpy()
    counts.insert(0, "dataset_version_id", np.repeat(owners, lengths))
    return counts.reset_index(drop=True)


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
    # pandas reads an empty Python list as a column selection.
    return scoped[np.array(keep, dtype=bool)]
