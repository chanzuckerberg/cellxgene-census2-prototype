"""Read `obs` columns of one store as Arrow arrays, and evaluate a query on them.

1. **Read.** `read_columns` reads the named columns directly with zarr. It follows
   anndata's on-disk format: it checks `encoding-type` and `encoding-version`, and an
   encoding outside `ENCODINGS` raises `ValueError`. A categorical column gives its
   codes and its categories (`Codes`), and no string for each cell. A categorical
   column of the catalog takes its categories from the catalog row, by the rule of the
   store format (see `read_column`). A numeric or boolean column gives its values
   (`Values`), and a string column, such as `observation_joinid`, its strings
   (`Strings`).
2. **Evaluate.** `evaluate` puts the columns in an Arrow table, as
   `DictionaryArray`s over the codes, and projects the query's pyarrow expression to
   one boolean for each cell with Acero. The cells where it is true are kept. pyarrow
   decodes a dictionary array to one string for each cell for `==` and `is_in`; the
   cost is small next to the reads (see `spikes/obs-arrow-native.md`).
3. **Result.** `ColumnData.take` gives the values of the kept cells as a pandas
   array. A categorical column keeps only the categories that the kept cells use.
"""

from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.compute as pc  # type: ignore[import-untyped]
import zarr
from pyarrow import acero

__all__ = [
    "ENCODINGS",
    "OBS_THREADS",
    "Arrow",
    "Codes",
    "ColumnData",
    "Strings",
    "Values",
    "evaluate",
    "mask",
    "read_column",
    "read_columns",
]

# The encodings that `read_column` reads itself, with their versions.
ENCODINGS = frozenset(
    {
        ("categorical", "0.2.0"),
        ("array", "0.2.0"),
        ("nullable-integer", "0.1.0"),
        ("nullable-boolean", "0.1.0"),
        ("string-array", "0.2.0"),
        ("nullable-string-array", "0.1.0"),
    }
)
# The columns of one store that the SDK reads at the same time.
OBS_THREADS = 8


@dataclass(frozen=True)
class Codes:
    """A categorical column: codes into `categories`, with -1 for a null."""

    codes: np.ndarray
    categories: Sequence[str]
    ordered: bool

    def arrow(self) -> Any:
        """A `DictionaryArray` over the codes, without a copy when no code is null."""
        nulls = self.codes < 0
        return pa.DictionaryArray.from_arrays(
            self.codes,
            pa.array(self.categories, type=pa.string()),
            mask=nulls if nulls.any() else None,
            ordered=self.ordered,
        )

    def take(self, keep: np.ndarray) -> Any:
        """The kept cells, with only the categories that they use."""
        codes = self.codes[keep]
        n = len(self.categories)
        used = np.flatnonzero(np.bincount(codes[codes >= 0], minlength=n))
        # A null's code -1 indexes the -1 at the end.
        renumber = np.full(n + 1, -1, dtype=np.int64)
        renumber[used] = np.arange(len(used))
        # Slice an Index of all the categories, so that an empty slice keeps the dtype.
        index = pd.Index(list(self.categories))[used]
        dtype = pd.CategoricalDtype(index, ordered=self.ordered)
        return pd.Categorical.from_codes(renumber[codes], dtype=dtype, validate=False)


@dataclass(frozen=True)
class Values:
    """Numbers or booleans, with `mask` true for a null, or None when none is null."""

    values: np.ndarray
    mask: np.ndarray | None

    def arrow(self) -> Any:
        """The values as an Arrow array; a NaN stays a NaN, not a null."""
        return pa.array(self.values, mask=self.mask)

    def take(self, keep: np.ndarray) -> Any:
        values = self.values[keep]
        if self.mask is None:
            return values
        mask = self.mask[keep]
        if values.dtype == np.bool_:
            return pd.arrays.BooleanArray(values, mask)
        return pd.arrays.IntegerArray(values, mask)


@dataclass(frozen=True)
class Strings:
    """Strings, with `mask` true for a null, or None when none is null."""

    values: np.ndarray
    mask: np.ndarray | None
    dtype: Any

    def arrow(self) -> Any:
        return pa.array(self.values, type=pa.string(), mask=self.mask)

    def take(self, keep: np.ndarray) -> Any:
        """The kept strings, with the dtype that `anndata` reads."""
        result: Any = pd.array(self.values[keep], dtype=self.dtype)
        if self.mask is not None:
            result[self.mask[keep]] = pd.NA
        return result


@dataclass(frozen=True)
class Arrow:
    """A column that is an Arrow array already, as the catalog test makes one."""

    array: Any

    def arrow(self) -> Any:
        return self.array


ColumnData = Codes | Values | Strings
FilterData = ColumnData | Arrow


def encoding(node: zarr.Array[Any] | zarr.Group) -> tuple[str, str]:
    return (
        str(node.attrs.get("encoding-type", "")),
        str(node.attrs.get("encoding-version", "")),
    )


def _array(group: zarr.Group, name: str) -> zarr.Array[Any]:
    node = group[name]
    if not isinstance(node, zarr.Array):
        raise TypeError(f"{group.store_path}/{name} must be an array")
    return node


def read_column(
    node: zarr.Array[Any] | zarr.Group, listed: Sequence[str] | None = None
) -> ColumnData:
    """One `obs` column, without a Python object for each cell.

    `listed` is the catalog row's sorted value set of a categorical column, or None.
    The store format gives a categorical column the observed values, sorted, as its
    categories, with no unused category. Thus, when the store has as many categories
    as the row lists, the categories are the row's values, and the read takes only
    the codes. Otherwise it reads the categories of the store.
    """
    kind, version = encoding(node)
    if (kind, version) not in ENCODINGS:
        raise ValueError(
            f"{node.store_path}: the SDK does not read the obs encoding "
            f"{kind or '(none)'} {version}"
        )
    if kind == "categorical" and isinstance(node, zarr.Group):
        codes, categories = _array(node, "codes"), _array(node, "categories")
        ordered = bool(node.attrs.get("ordered", False))
        if listed is not None and categories.shape == (len(listed),):
            labels: Sequence[str] = listed
        else:
            labels = [str(c) for c in np.asarray(categories[...])]
        return Codes(np.asarray(codes[...]), labels, ordered)
    if kind == "array" and isinstance(node, zarr.Array):
        values = np.asarray(node[...])
        if values.dtype.kind not in "biuf":
            raise ValueError(f"{node.store_path}: an array of {values.dtype}")
        return Values(values, None)
    if kind == "string-array" and isinstance(node, zarr.Array):
        return Strings(np.asarray(node[...]), None, object)
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path}: {kind} must be an encoded group")
    values = np.asarray(_array(node, "values")[...])
    mask = np.asarray(_array(node, "mask")[...], dtype=bool)
    if kind == "nullable-string-array":
        na_value = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
        dtype = pd.StringDtype(na_value=na_value)
        return Strings(values, mask if mask.any() else None, dtype)
    return Values(values, mask if mask.any() else None)


def read_columns(
    obs: zarr.Group,
    names: Iterable[str],
    listed: Mapping[str, Sequence[str]],
) -> dict[str, ColumnData]:
    """`read_column` for each of `names`, as many as `OBS_THREADS` at the same time.

    `listed` gives the catalog's value set of each categorical column that it lists.
    Each column is one or two requests, so the reads run at the same time. The pool
    waits on no other task, so it cannot deadlock the pool of the query.
    """
    names = list(names)

    def one(name: str) -> ColumnData:
        return read_column(obs[name], listed.get(name))

    if len(names) <= 1:
        return {name: one(name) for name in names}
    with ThreadPoolExecutor(max_workers=min(len(names), OBS_THREADS)) as pool:
        return dict(zip(names, pool.map(one, names), strict=True))


def table(data: Mapping[str, FilterData], n: int) -> Any:
    """The columns as one Arrow table of `n` rows."""
    if not data:
        return pa.table({"__row": pa.nulls(n, pa.bool_())})
    return pa.table({name: d.arrow() for name, d in data.items()})


def evaluate(expression: Any, data: Mapping[str, FilterData], n: int) -> np.ndarray:
    """The positions of the cells where `expression` is true, in order.

    A null, from a null in a column, keeps no cell.
    """
    mask = _project(expression, data, n)
    return np.asarray(pc.indices_nonzero(mask).to_numpy(), dtype=np.int64)


def mask(
    expression: Any, data: Mapping[str, FilterData], n: int
) -> tuple[np.ndarray, np.ndarray]:
    """`expression` for each row: its value, and True where it is not null."""
    result = _project(expression, data, n)
    valid = np.asarray(result.is_valid().to_numpy(zero_copy_only=False), dtype=bool)
    values = np.asarray(result.fill_null(False).to_numpy(zero_copy_only=False))
    return values.astype(bool), valid


def _project(expression: Any, data: Mapping[str, FilterData], n: int) -> Any:
    """Acero projects the expression to one boolean for each row, in this thread."""
    plan = acero.Declaration.from_sequence(
        [
            acero.Declaration(
                "table_source", acero.TableSourceNodeOptions(table(data, n))
            ),
            acero.Declaration(
                "project", acero.ProjectNodeOptions([expression], ["keep"])
            ),
        ]
    )
    result = plan.to_table(use_threads=False)["keep"].combine_chunks()
    if result.type != pa.bool_():
        raise TypeError(f"the query gives {result.type}, not bool")
    return result


def take_rows(node: zarr.Array[Any] | zarr.Group, rows: np.ndarray) -> Any:
    """The values of one string column at the row positions `rows`.

    For a nullable string column, the SDK makes strings for these rows only. The
    dtype is the one that `anndata` reads for the same column. A column in another
    encoding is read whole with `read_column`.
    """
    if node.attrs.get("encoding-type") != "nullable-string-array":
        return read_column(node).take(rows)
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded group")
    na_value = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
    dtype = pd.StringDtype(na_value=na_value)
    if not rows.size:
        return pd.array([], dtype=dtype)
    values, mask = _array(node, "values"), _array(node, "mask")
    result: Any = pd.array(np.asarray(values.oindex[rows]), dtype=dtype)
    result[np.asarray(mask.oindex[rows])] = pd.NA
    return result
