"""Test a query against the catalog with pyarrow: evaluate it on the value sets.

A catalog row records the value set of each indexed column of its store. Each cell
of the store holds one value of each set, so its values are one combination of the
sets. `test` makes a table of every combination of a row's value sets, for the
query's columns, and evaluates the query's own pyarrow expression on it with
`obs_read.evaluate`. The SDK does not walk the expression.

- A column without a value set in the row, such as `nnz`, `donor_id`,
  `observation_joinid`, or a column that the row does not record, is a null in each
  combination: a value that the
  catalog does not know. pyarrow follows SQL's three-valued logic, so a term on it
  gives null, unless the other terms decide.
- The row is excluded when every combination gives false. A null, or a true,
  keeps the row.
- The row is whole when every combination gives true and no column is unknown.
  The SDK then reads no filter column of the store.

`dataset_id` holds the row's identity, and `is_primary_data` the catalog's summary
(`all_true`, `all_false`, or `mixed`). A multi-valued column holds the row's joined
values, as the store does.

The test is sound because each term of the builder gives null for a null value:
a cell that matches makes its combination true or null, never false. A row whose
table would pass `MAX_COMBINATIONS` rows keeps the hint test (`catalog.excludes` and
`catalog.covers`). The whole verdict relies on the store format's rule that an
indexed column holds no null, as `catalog.covers` does.
"""

from collections.abc import Sequence
from typing import Any

import numpy as np
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types

from census2_sdk import obs_read
from census2_sdk.catalog import INDEXED, unrecorded
from census2_sdk.query import Query
from census2_sdk.schema import OBS_COLUMNS, ColType

__all__ = ["MAX_COMBINATIONS", "test"]

# The largest table of one row. Above it the row keeps the hint test.
MAX_COMBINATIONS = 100_000
PRIMARY = {"all_true": [True], "all_false": [False], "mixed": [True, False]}
# The Arrow type of the nulls of a column that no row knows, by its `ColType`.
NULL_TYPES = {
    ColType.CATEGORY: pa.string(),
    ColType.STRING: pa.string(),
    ColType.BOOL: pa.bool_(),
    ColType.INT64: pa.int64(),
    ColType.FLOAT64: pa.float64(),
}
COL_TYPES = dict(OBS_COLUMNS)


def value_set(row: dict[str, Any], column: str) -> list[Any] | None:
    """The values that the row records for `column`, or None when it records none.

    An empty set of a store with cells is not recorded: its cells hold nulls.
    """
    if column == "dataset_id":
        return [str(row["dataset_id"])]
    found: list[Any] | None = None
    if column == "is_primary_data":
        found = PRIMARY.get(str(row.get(column)))
    elif column in INDEXED and not unrecorded(row.get(column)):
        found = [str(v) for v in row[column]]
    if found is not None and not found and int(row.get("n_obs") or 0) > 0:
        return None
    return found


def test(
    rows: Sequence[dict[str, Any]], q: Query
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each row: kept, whole, and tested (False where the row needs the hints).

    A row that is not tested is kept and not whole here; the caller applies the
    hint test to it.
    """
    n = len(rows)
    keep = np.ones(n, dtype=bool)
    whole = np.zeros(n, dtype=bool)
    tested = np.zeros(n, dtype=bool)
    if not n or q.expression is None:
        return keep, whole, tested
    columns = sorted(q.columns)
    sets = [[value_set(row, c) for c in columns] for row in rows]
    sizes = [int(np.prod([1 if s is None else len(s) for s in r])) for r in sets]
    tested = np.array([size <= MAX_COMBINATIONS for size in sizes], dtype=bool)
    chosen = np.flatnonzero(tested)
    counts = np.array([sizes[i] for i in chosen], dtype=np.int64)
    total = int(counts.sum())
    data: dict[str, obs_read.FilterData] = {}
    for j, column in enumerate(columns):
        data[column] = _column(column, [sets[i] for i in chosen], j, counts, total)
    values, valid = obs_read.mask(q.expression, data, total)
    owner = np.repeat(np.arange(len(chosen)), counts)
    true = np.bincount(owner, weights=values & valid, minlength=len(chosen))
    false = np.bincount(owner, weights=~values & valid, minlength=len(chosen))
    known = np.array([all(s is not None for s in sets[i]) for i in chosen], dtype=bool)
    keep[chosen] = false < counts
    whole[chosen] = known & (counts > 0) & (true == counts)
    return keep, whole, tested


def _column(
    column: str,
    sets: list[list[list[Any] | None]],
    j: int,
    counts: np.ndarray,
    total: int,
) -> obs_read.FilterData:
    """Column `j` of the combination tables of the rows, one after another."""
    known = [r[j] for r in sets if r[j] is not None]
    if not known:
        col_type = COL_TYPES.get(column)
        null_type = pa.null() if col_type is None else NULL_TYPES[col_type]
        return obs_read.Arrow(pa.nulls(total, null_type))
    parts: list[np.ndarray] = []
    nulls: list[np.ndarray] = []
    labels: dict[Any, int] = {}
    for r, count in zip(sets, counts, strict=True):
        if not count:
            continue
        values = r[j]
        if values is None:
            parts.append(np.zeros(count, dtype=np.int32))
            nulls.append(np.ones(count, dtype=bool))
            continue
        codes = np.array([labels.setdefault(v, len(labels)) for v in values])
        # The combinations of one row in C order: this column repeats each value
        # once for each combination of the columns after it, and tiles over the
        # columns before it.
        after = int(np.prod([1 if s is None else len(s) for s in r[j + 1 :]]))
        parts.append(np.tile(np.repeat(codes, after), count // (len(values) * after)))
        nulls.append(np.zeros(count, dtype=bool))
    codes_all = np.concatenate(parts) if parts else np.zeros(0, dtype=np.int32)
    mask = np.concatenate(nulls) if nulls else np.zeros(0, dtype=bool)
    dictionary = list(labels)
    if column == "is_primary_data":
        values_all = (
            np.array(dictionary, dtype=bool)[codes_all]
            if dictionary
            else codes_all.astype(bool)
        )
        return obs_read.Values(values_all, mask if mask.any() else None)
    return obs_read.Arrow(
        pa.DictionaryArray.from_arrays(
            codes_all.astype(np.int32),
            pa.array([str(v) for v in dictionary], type=pa.string()),
            mask=mask if mask.any() else None,
        )
    )
