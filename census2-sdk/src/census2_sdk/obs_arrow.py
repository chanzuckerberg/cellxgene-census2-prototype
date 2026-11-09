"""Bind a resolved filter to the `obs` of one store, and evaluate it with pyarrow.

This is the store side of a query, for one read: all the rows of a per-dataset store,
or the rows of adjacent members of a consolidated store.

1. **Read.** `read_columns` reads only the columns that the residual filter names, and
   only the rows of the read, directly with zarr. It follows anndata's on-disk format:
   it checks `encoding-type` and `encoding-version`, and reads an encoding that it does
   not know with `anndata.io.read_elem`. A categorical column gives its codes and its
   categories (`Codes`), and no string for each cell. A numeric or boolean column gives
   its values (`Values`). The read does not touch `observation_joinid` or the index.
2. **Bind.** `bind` replaces each term with what this store needs. A term on labels
   (`==`, `isin`, `.str`, `has_any`, `has_all`) is evaluated once for each category,
   with the same test as the catalog (`obs_filter.label_test`), and then mapped to the
   cells by their codes: a code of -1, a null, maps to False. A term that holds for no
   category becomes False, and a term that holds for every category of a column
   without nulls becomes True, so the bound filter can fold to a constant. A numeric,
   boolean, or null term stays a term on its Arrow column.
3. **Evaluate.** `compile_expression` turns the bound filter into one
   `pyarrow.compute.Expression`, and `evaluate` projects it over an Arrow table with
   Acero, in the calling thread. The row positions of the True values are the kept
   rows. A term on a column with nulls is filled with False before `~`, `&`, and `|`
   combine it, which gives the pandas rule.

pyarrow 26 evaluates `==` and `is_in` on a dictionary array by decoding it to one
string for each cell, and it has no dictionary kernel for the string functions. The
evaluation by codes of step 2 avoids that work; see `spikes/consolidation-filters.md`.

`to_table` gives the read columns as an Arrow table: a categorical column becomes a
`DictionaryArray` over its codes without a copy, with a null for each code of -1; a
numeric column wraps its NumPy values without a copy; a boolean column packs its
values into bits.
"""

import threading
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import reduce
from typing import Any, cast

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.compute as pc  # type: ignore[import-untyped]
import zarr
from anndata.io import read_elem
from pyarrow import acero

from census2_sdk.obs_filter import (
    FALSE,
    TRUE,
    And,
    Compare,
    Const,
    Expr,
    FilterError,
    IsNull,
    Leaf,
    Not,
    Or,
    Truth,
    is_label_leaf,
    label_test,
    normalize,
    render,
)

# The encodings that `read_column` reads itself, with their versions.
ENCODINGS = frozenset(
    {
        ("categorical", "0.2.0"),
        ("array", "0.2.0"),
        ("string-array", "0.2.0"),
        ("nullable-integer", "0.1.0"),
        ("nullable-boolean", "0.1.0"),
        ("nullable-string-array", "0.1.0"),
    }
)

# The columns of one read that the SDK reads at the same time.
OBS_THREADS = 8


@dataclass(frozen=True)
class Codes:
    """A categorical column: codes into `categories`, with -1 for a null."""

    codes: np.ndarray
    categories: np.ndarray
    ordered: bool

    def __len__(self) -> int:
        return len(self.codes)

    def part(self, start: int, stop: int) -> "Codes":
        return Codes(self.codes[start:stop], self.categories, self.ordered)

    def has_null(self) -> bool:
        return bool(len(self.codes)) and int(self.codes.min()) < 0

    def arrow(self) -> Any:
        """A `DictionaryArray` over the codes, which it wraps without a copy."""
        indices = pa.array(self.codes)
        if self.has_null():
            indices = pa.array(self.codes, mask=self.codes < 0)
        dictionary = pa.array(self.categories, type=pa.string())
        return pa.DictionaryArray.from_arrays(indices, dictionary, ordered=self.ordered)


@dataclass(frozen=True)
class Values:
    """A column of values, with `mask` True for a null, or None when none is null."""

    values: np.ndarray
    mask: np.ndarray | None

    def __len__(self) -> int:
        return len(self.values)

    def part(self, start: int, stop: int) -> "Values":
        mask = None if self.mask is None else self.mask[start:stop]
        return Values(self.values[start:stop], mask)

    def kind(self) -> str:
        """bool, int, float, or string."""
        dtype = self.values.dtype
        if dtype == np.bool_:
            return "bool"
        if np.issubdtype(dtype, np.integer):
            return "int"
        if np.issubdtype(dtype, np.floating):
            return "float"
        return "string"

    def arrow(self) -> Any:
        """The values as an Arrow array: numbers without a copy, booleans as bits."""
        if self.kind() == "string":
            return pa.array(self.values, type=pa.string(), mask=self.mask)
        return pa.array(self.values, mask=self.mask)

    def as_codes(self) -> Codes:
        """A string column as codes over its distinct values, for a label term."""
        encoded = self.arrow().dictionary_encode()
        codes = encoded.indices.fill_null(-1).to_numpy(zero_copy_only=False)
        categories = encoded.dictionary.to_numpy(zero_copy_only=False)
        return Codes(np.asarray(codes, dtype=np.int64), categories, False)


ColumnData = Codes | Values


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
    node: zarr.Array[Any] | zarr.Group,
    start: int,
    stop: int,
    categories: Callable[[], np.ndarray] | None = None,
) -> ColumnData:
    """Rows `start` to `stop` of one `obs` column, without a Python object per cell.

    `categories` gives the categories of a categorical column when the caller holds
    them, as a consolidated store does for all of its members. An encoding outside
    `ENCODINGS` reads whole with `read_elem`, and then the rows are cut from it.
    """
    kind, version = encoding(node)
    if (kind, version) not in ENCODINGS:
        return _fallback(node, start, stop)
    if kind == "categorical" and isinstance(node, zarr.Group):
        codes = np.asarray(_array(node, "codes")[start:stop])
        if categories is None:
            cats = np.asarray(_array(node, "categories")[...])
        else:
            cats = categories()
        return Codes(codes, cats, bool(node.attrs.get("ordered", False)))
    if kind in ("array", "string-array") and isinstance(node, zarr.Array):
        return Values(np.asarray(node[start:stop]), None)
    if isinstance(node, zarr.Group):
        values = np.asarray(_array(node, "values")[start:stop])
        mask = np.asarray(_array(node, "mask")[start:stop], dtype=bool)
        return Values(values, mask if mask.any() else None)
    return _fallback(node, start, stop)


def _fallback(node: zarr.Array[Any] | zarr.Group, start: int, stop: int) -> ColumnData:
    """Read a column of an unknown encoding with anndata, then cut the rows."""
    value = cast("Any", read_elem(node))
    if isinstance(value, pd.Categorical):
        return Codes(
            np.asarray(value.codes[start:stop]),
            np.asarray(value.categories.astype(str)),
            bool(value.ordered),
        )
    series = pd.Series(value[start:stop])
    if isinstance(series.dtype, pd.CategoricalDtype):
        cat = series.cat
        return Codes(
            np.asarray(cat.codes), np.asarray(cat.categories.astype(str)), False
        )
    mask = series.isna().to_numpy()
    if pd.api.types.is_bool_dtype(series.dtype):
        values = series.to_numpy(dtype=bool, na_value=False)
    elif pd.api.types.is_integer_dtype(series.dtype):
        values = series.to_numpy(dtype=np.int64, na_value=0)
    elif pd.api.types.is_float_dtype(series.dtype):
        values = series.to_numpy(dtype=np.float64, na_value=np.nan)
    else:
        values = series.astype(object).where(~mask, "").to_numpy()
    return Values(np.asarray(values), mask if mask.any() else None)


class CategoryCache:
    """The categories of each categorical `obs` column of one store, read once."""

    def __init__(self, obs: zarr.Group, known: Mapping[str, pd.Index] | None = None):
        self._obs = obs
        self._lock = threading.Lock()
        self._done: dict[str, np.ndarray] = {
            k: np.asarray(v) for k, v in (known or {}).items()
        }

    def __call__(self, name: str) -> np.ndarray:
        with self._lock:
            found = self._done.get(name)
        if found is not None:
            return found
        node = self._obs[name]
        if not isinstance(node, zarr.Group):
            raise TypeError(
                f"{self._obs.store_path}/{name} must be a categorical group"
            )
        cats = np.asarray(_array(node, "categories")[...])
        with self._lock:
            self._done[name] = cats
        return cats


def read_columns(
    obs: zarr.Group,
    names: Iterable[str],
    start: int,
    stop: int,
    categories: CategoryCache | None = None,
) -> dict[str, ColumnData]:
    """`read_column` for each of `names`, up to `OBS_THREADS` at the same time.

    Each column is a read of its own, with a round trip or two, so the reads run at
    the same time. The pool waits on no other task, so it cannot deadlock the pool of
    the query.
    """
    names = list(names)

    def one(name: str) -> ColumnData:
        node = obs[name]
        cats = None
        if categories is not None and encoding(node)[0] == "categorical":
            cats = lambda: categories(name)
        return read_column(node, start, stop, cats)

    if len(names) <= 1:
        return {name: one(name) for name in names}
    with ThreadPoolExecutor(max_workers=min(len(names), OBS_THREADS)) as pool:
        return dict(zip(names, pool.map(one, names), strict=True))


def to_table(data: Mapping[str, ColumnData]) -> Any:
    """The columns as an Arrow table, without a Python object per cell for codes."""
    return pa.table({name: d.arrow() for name, d in data.items()})


# Binding.


@dataclass(frozen=True, repr=False)
class Field(Leaf):
    """A bound label term: the boolean column `name` of the bound table, with no null."""

    name: str


@dataclass(frozen=True)
class Bound:
    """A filter bound to one read: the tree, and the Arrow columns that it reads."""

    tree: Expr
    columns: dict[str, Any]


def _kind_error(leaf: Leaf, kind: str) -> FilterError:
    return FilterError(
        f"{leaf.column} holds {kind} values; {render(leaf)} does not apply"
    )


def bind(tree: Expr, data: Mapping[str, ColumnData]) -> Bound:
    """Bind `tree` to the columns of one read. See the module docstring, step 2."""
    columns: dict[str, Any] = {}
    codes_of: dict[str, Codes] = {}

    def codes(name: str) -> Codes:
        found = codes_of.get(name)
        if found is None:
            d = data[name]
            found = codes_of[name] = d if isinstance(d, Codes) else d.as_codes()
        return found

    def mask_leaf(leaf: Leaf, mask: np.ndarray) -> Expr:
        if not mask.any():
            return FALSE
        if mask.all():
            return TRUE
        name = f"__term{len(columns)}"
        columns[name] = pa.array(mask)
        return Field(leaf.organism, leaf.column, leaf.hatch, name)

    def values_column(leaf: Leaf, d: Values) -> None:
        columns.setdefault(leaf.column, d.arrow())

    def one(leaf: Leaf) -> Expr:
        d = data[leaf.column]
        if is_label_leaf(leaf):
            if isinstance(d, Values) and d.kind() != "string":
                raise _kind_error(leaf, d.kind())
            c = codes(leaf.column)
            test = label_test(leaf)
            lut = np.fromiter(
                (test(str(v)) for v in c.categories), bool, len(c.categories)
            )
            if not lut.any():
                return FALSE
            if lut.all() and not c.has_null():
                return TRUE
            return mask_leaf(leaf, np.append(lut, False)[c.codes])
        if isinstance(leaf, IsNull):
            if isinstance(d, Codes):
                return mask_leaf(leaf, d.codes < 0)
            nulls = np.zeros(len(d), bool) if d.mask is None else d.mask
            if d.kind() == "float":
                nulls = nulls | np.isnan(d.values)
            return mask_leaf(leaf, nulls)
        if isinstance(d, Codes):
            raise _kind_error(leaf, "label")
        if isinstance(leaf, Truth):
            if d.kind() != "bool":
                raise _kind_error(leaf, d.kind())
            if d.mask is None and d.values.all():
                return TRUE
            if not d.values.any():
                return FALSE
            values_column(leaf, d)
            return leaf
        if isinstance(leaf, Compare):
            if d.kind() not in ("int", "float") or isinstance(leaf.value, str | bool):
                raise _kind_error(leaf, d.kind())
            values_column(leaf, d)
            return leaf
        raise FilterError(f"{render(leaf)} does not apply to {leaf.column}")

    def walk(node: Expr) -> Expr:
        if isinstance(node, Leaf):
            return one(node)
        if isinstance(node, Not):
            return Not(node.organism, walk(node.child))
        if isinstance(node, And | Or):
            return node.__class__(node.organism, tuple(walk(c) for c in node.children))
        return node

    bound = normalize(walk(tree))
    used = {f.name for f in _fields(bound)} | {
        leaf.column for leaf in _value_leaves(bound)
    }
    return Bound(bound, {k: v for k, v in columns.items() if k in used})


def _fields(node: Expr) -> list[Field]:
    if isinstance(node, Field):
        return [node]
    if isinstance(node, Not):
        return _fields(node.child)
    if isinstance(node, And | Or):
        return [f for c in node.children for f in _fields(c)]
    return []


def _value_leaves(node: Expr) -> list[Leaf]:
    if isinstance(node, Compare | Truth):
        return [node]
    if isinstance(node, Not):
        return _value_leaves(node.child)
    if isinstance(node, And | Or):
        return [leaf for c in node.children for leaf in _value_leaves(c)]
    return []


_OPS: dict[str, Callable[[Any, Any], Any]] = {
    "==": lambda a, b: a == b,
    "<": lambda a, b: a < b,
    "<=": lambda a, b: a <= b,
    ">": lambda a, b: a > b,
    ">=": lambda a, b: a >= b,
}


def compile_expression(tree: Expr, columns: Mapping[str, Any]) -> Any:
    """`tree` as one `pyarrow.compute.Expression` over the bound columns.

    A term on a column with nulls is filled with False, so that `~` keeps the nulls.
    """

    def filled(expr: Any, column: str) -> Any:
        array = columns.get(column)
        if array is not None and array.null_count:
            return pc.coalesce(expr, pa.scalar(False))
        return expr

    def walk(node: Expr) -> Any:
        if isinstance(node, Const):
            return pc.scalar(node.value)
        if isinstance(node, Field):
            return pc.field(node.name)
        if isinstance(node, Compare):
            return filled(_OPS[node.op](pc.field(node.column), node.value), node.column)
        if isinstance(node, Truth):
            return filled(pc.field(node.column), node.column)
        if isinstance(node, Not):
            return ~walk(node.child)
        if isinstance(node, And):
            return reduce(lambda a, b: a & b, (walk(c) for c in node.children))
        if isinstance(node, Or):
            return reduce(lambda a, b: a | b, (walk(c) for c in node.children))
        raise TypeError(f"{render(node)} is not bound")

    return walk(tree)


def evaluate(tree: Expr, data: Mapping[str, ColumnData], n: int) -> np.ndarray:
    """The positions of the rows of one read that match `tree`, in order."""
    if isinstance(tree, Const):
        return np.arange(n) if tree.value else np.zeros(0, dtype=np.int64)
    bound = bind(tree, data)
    node = bound.tree
    if isinstance(node, Const):
        return np.arange(n) if node.value else np.zeros(0, dtype=np.int64)
    if isinstance(node, Field):
        return np.flatnonzero(bound.columns[node.name].to_numpy(zero_copy_only=False))
    mask = project(compile_expression(node, bound.columns), bound.columns, n)
    return np.flatnonzero(mask)


def project(expression: Any, columns: Mapping[str, Any], n: int) -> np.ndarray:
    """Evaluate `expression` over `columns` with Acero, in this thread: one bool each."""
    table = pa.table(dict(columns)) if columns else pa.table({"__n": np.zeros(n, bool)})
    plan = acero.Declaration.from_sequence(
        [
            acero.Declaration("table_source", acero.TableSourceNodeOptions(table)),
            acero.Declaration("project", acero.ProjectNodeOptions([expression], ["m"])),
        ]
    )
    result = plan.to_table(use_threads=False)["m"]
    return np.asarray(result.to_numpy(zero_copy_only=False), dtype=bool)


def to_expression(tree: Expr) -> Any:
    """`tree` as a `pyarrow.compute.Expression` over decoded columns, as a reference.

    This form evaluates each label term with pyarrow on each cell's string, and the
    pandas rule fills each term with False. The query does not use it: it binds to
    the codes instead. The tests compare the two. A multi-valued term has no form here.
    """
    from census2_sdk.obs_filter import HasAll, HasAny, IsIn, StrTest

    def fill(expr: Any) -> Any:
        return pc.coalesce(expr, pa.scalar(False))

    def walk(node: Expr) -> Any:
        if isinstance(node, Const):
            return pc.scalar(node.value)
        if isinstance(node, Not):
            return ~walk(node.child)
        if isinstance(node, And):
            return reduce(lambda a, b: a & b, (walk(c) for c in node.children))
        if isinstance(node, Or):
            return reduce(lambda a, b: a | b, (walk(c) for c in node.children))
        if not isinstance(node, Leaf):
            raise TypeError(f"not a filter node: {node!r}")
        f = pc.field(node.column)
        if isinstance(node, IsIn):
            return fill(f.isin(pa.array(list(node.values))))
        if isinstance(node, Compare):
            return fill(_OPS[node.op](f, node.value))
        if isinstance(node, Truth):
            return fill(f)
        if isinstance(node, IsNull):
            return f.is_null(nan_is_null=True)
        if isinstance(node, StrTest):
            return fill(_string_function(node, f))
        if isinstance(node, HasAny | HasAll):
            raise TypeError("a multi-valued term has no reference form")
        raise TypeError(f"no reference form for {render(node)}")

    return walk(normalize(tree))


def _string_function(node: Any, f: Any) -> Any:
    pattern = node.pattern
    if node.method in ("startswith", "endswith"):
        fn = pc.starts_with if node.method == "startswith" else pc.ends_with
        prefixes = pattern if isinstance(pattern, tuple) else (pattern,)
        return reduce(lambda a, b: a | b, (fn(f, pattern=p) for p in prefixes))
    ignore = not node.case
    if node.method == "contains":
        if node.regex:
            return pc.match_substring_regex(f, pattern=pattern, ignore_case=ignore)
        return pc.match_substring(f, pattern=pattern, ignore_case=ignore)
    if node.method == "match":
        return pc.match_substring_regex(
            f, pattern=f"^(?:{pattern})", ignore_case=ignore
        )
    return pc.match_substring_regex(f, pattern=f"^(?:{pattern})$", ignore_case=ignore)


# The result: the kept rows of the other columns.


def take(d: ColumnData, positions: np.ndarray) -> Any:
    """The values of `d` at `positions`, as a pandas array for the result `obs`.

    A categorical column keeps only the categories that the kept rows use. pandas
    does not check the codes again (`validate=False`): they come from the store.
    """
    if isinstance(d, Codes):
        codes = d.codes[positions]
        present = codes[codes >= 0]
        used = np.flatnonzero(np.bincount(present, minlength=len(d.categories)))
        renumber = np.full(len(d.categories) + 1, -1, dtype=np.int64)
        renumber[used] = np.arange(len(used))
        # Slice an Index of all the categories, so that an empty slice keeps the dtype.
        dtype = pd.CategoricalDtype(pd.Index(d.categories)[used], ordered=d.ordered)
        return pd.Categorical.from_codes(renumber[codes], dtype=dtype, validate=False)
    if d.mask is None:
        return d.values[positions]
    return _masked(d, positions)


def _masked(d: Values, positions: np.ndarray) -> Any:
    values, mask = d.values[positions], cast("np.ndarray", d.mask)[positions]
    kind = d.kind()
    if kind == "bool":
        return pd.arrays.BooleanArray(values, mask)
    if kind == "int":
        return pd.arrays.IntegerArray(values, mask)
    if kind == "float":
        return pd.arrays.FloatingArray(values, mask)
    out: Any = pd.array(values.astype(object), dtype=pd.StringDtype(na_value=np.nan))
    out[mask] = pd.NA
    return out


def read_rows(node: zarr.Array[Any] | zarr.Group, rows: np.ndarray) -> Any:
    """The values of one `obs` column at the row positions `rows`.

    For a nullable string column, the SDK makes strings for these rows only. The
    dtype is the one that `anndata` reads for the same column.
    """
    if node.attrs.get("encoding-type") != "nullable-string-array":
        return pd.array(cast("Any", read_elem(node)))[rows]
    if not isinstance(node, zarr.Group):
        raise TypeError(f"{node.store_path} must be an encoded group")
    na_value = np.nan if node.attrs.get("na-value", "NA") == "NaN" else pd.NA
    dtype = pd.StringDtype(na_value=na_value)
    if not rows.size:
        return pd.array([], dtype=dtype)
    values, mask = node["values"], node["mask"]
    if not (isinstance(values, zarr.Array) and isinstance(mask, zarr.Array)):
        raise TypeError(f"{node.store_path} must hold a values and a mask array")
    # One run of rows reads as a slice, which zarr reads faster than a selection.
    contiguous = rows[-1] - rows[0] + 1 == rows.size
    selection: Any = slice(int(rows[0]), int(rows[-1]) + 1) if contiguous else rows
    result: Any = pd.array(np.asarray(values.oindex[selection]), dtype=dtype)
    result[np.asarray(mask.oindex[selection])] = pd.NA
    return result
