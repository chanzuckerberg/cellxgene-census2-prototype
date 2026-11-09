"""Resolve a filter once for each query, and test it against the catalog.

`resolve` gives a `Plan` in four steps:

1. Check the organism of the filter, each column, and each operator. A column of the
   static schema (`obs_columns.KINDS`) takes the operators of its kind. A column that
   only `obs["name"]` names waits for the stores to give its type.
2. Check each label against the catalog's value sets, the union over the catalog rows
   of the organism. A term whose labels are all unknown raises `FilterError` with the
   nearest known labels. An unknown label with a near known label, next to known
   labels, gives a warning. Other unknown labels match no cell, as in pandas.
3. Normalize the filter (`obs_filter.normalize`).
4. Test the filter against each catalog row of the organism. The test gives one of
   three answers for each term: every cell of the dataset matches (`ALL`), no cell
   matches (`NONE`), or the catalog cannot tell (`SOME`). `&`, `|`, and `~` combine
   the answers. A row is a candidate unless its answer is `NONE`, and its residual is
   the filter with the `ALL` and `NONE` parts replaced by constants. A residual of
   `True` selects every cell of the dataset, so the stores need not read a column to
   find them.

The test is exact for one term on a column that the catalog records, and
conservative for several terms. A numeric term proves nothing, because the catalog
records no numeric range. `~` turns `NONE` into `ALL` and back: a term is `ALL` only
when every value in the row's value set passes it and no cell holds a null. The store
format holds no null in the columns of the value sets (the writer refuses one), so a
value set is complete. A row that records no value set proves nothing.
"""

import difflib
import warnings
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Literal, cast

import numpy as np
import pandas as pd

from census2_sdk import obs_columns
from census2_sdk.catalog import INDEXED, members, unrecorded
from census2_sdk.obs_filter import (
    And,
    Between,
    Compare,
    Const,
    Expr,
    FilterError,
    FilterLike,
    HasAll,
    HasAny,
    IsIn,
    IsNull,
    Leaf,
    Not,
    Or,
    StrTest,
    Truth,
    as_expr,
    is_const,
    is_label_leaf,
    label_test,
    leaves,
    normalize,
    render,
)

Proof = Literal["all", "none", "some"]
ALL: Proof = "all"
NONE: Proof = "none"
SOME: Proof = "some"

# The columns with a value set in each catalog row. `dataset_id` is the row's identity.
VALUE_SET_COLUMNS = frozenset(INDEXED)
# The operators of each kind of column.
LABEL_KINDS = frozenset({"category", "string"})


@dataclass(frozen=True)
class Candidate:
    """A catalog row that can hold a matching cell, and the filter of its cells."""

    row: dict[str, Any]
    residual: Expr


@dataclass(frozen=True)
class Plan:
    """A filter resolved against a catalog: the tree, its columns, and the candidates."""

    organism: str
    tree: Expr
    columns: frozenset[str]
    candidates: tuple[Candidate, ...]
    catalog_rows: int

    def __repr__(self) -> str:
        whole = sum(1 for c in self.candidates if is_const(c.residual, True))
        return (
            f"<Plan {self.organism}: {render(self.tree)}; {len(self.candidates)} of "
            f"{self.catalog_rows} datasets, {whole} whole>"
        )


ValueSets = dict[str, frozenset[str]]


def value_sets(rows: pd.DataFrame) -> ValueSets:
    """The known labels of each column with a value set, over `rows`.

    A multi-valued column gives its members. `dataset_id` gives the identities.
    """
    out: ValueSets = {"dataset_id": frozenset(str(v) for v in rows["dataset_id"])}
    for column in sorted(VALUE_SET_COLUMNS & set(rows.columns)):
        seen: set[str] = set()
        for values in rows[column]:
            if unrecorded(values):
                continue
            for value in cast("Any", values):
                seen.update(members(str(value), column))
        out[column] = frozenset(seen)
    return out


def resolve(
    obs_filter: FilterLike | None,
    organism: str,
    catalog: pd.DataFrame,
    known: ValueSets | None = None,
) -> Plan:
    """Check, normalize, and test `obs_filter` against the catalog rows of `organism`.

    `known` is `value_sets` of those rows, which a caller can keep for later queries.
    """
    tree = Const(organism, True) if obs_filter is None else as_expr(obs_filter)
    if tree.organism and tree.organism != organism:
        raise FilterError(
            f"the filter is over the columns of {tree.organism}, and the query is over "
            f"{organism}; build it from census.{organism}.obs"
        )
    rows = catalog[catalog["organism"] == organism]
    for leaf in leaves(tree):
        check_leaf(leaf)
    if len(rows):
        known = value_sets(rows) if known is None else known
        for leaf in leaves(normalize(tree)):
            check_labels(leaf, known)
    tree = normalize(tree)
    records = cast("list[dict[str, Any]]", rows.to_dict("records"))
    tests: dict[Leaf, Callable[[str], bool]] = {}
    candidates = []
    for row in records:
        proof, residual = prove(tree, row, tests)
        if proof != NONE:
            candidates.append(Candidate(row, residual))
    used = frozenset(leaf.column for leaf in leaves(tree))
    return Plan(organism, tree, used, tuple(candidates), len(records))


def check_leaf(leaf: Leaf) -> None:
    """Refuse an operator that the column's kind does not take.

    The typed columns allow only the right operators, so this check matters for the
    columns of `obs["name"]` that the schema names. A column that the schema does not
    name waits for the stores.
    """
    kind = obs_columns.KINDS.get(leaf.column)
    if kind is not None:
        check_kind(leaf, kind)


def check_kind(leaf: Leaf, kind: str) -> None:
    """Refuse `leaf` on a column of `kind`: category, multi, bool, string, int, float."""
    if isinstance(leaf, IsNull):
        return
    name = leaf.column
    if kind == "multi":
        if not isinstance(leaf, HasAny | HasAll):
            raise FilterError(
                f"{name} is multi-valued: a cell holds members joined by ' || '. "
                f"Use obs.{name}.has_any([...]) or obs.{name}.has_all([...]), not "
                f"{render(leaf)}"
            )
        return
    if isinstance(leaf, HasAny | HasAll):
        raise FilterError(
            f"{name} is not multi-valued; use obs.{name}.isin([...]), not {render(leaf)}"
        )
    if kind in LABEL_KINDS:
        if not is_label_leaf(leaf):
            raise FilterError(f"{name} holds labels; {render(leaf)} does not apply")
        if kind == "string" and isinstance(leaf, StrTest):
            raise FilterError(f"{name} takes == and isin, not {render(leaf)}")
        return
    if kind == "bool":
        if not isinstance(leaf, Truth):
            raise FilterError(
                f"{name} is boolean; use obs.{name} or ~obs.{name}, not {render(leaf)}"
            )
        return
    if kind in ("int", "float"):
        numeric = isinstance(leaf, Between) or (
            isinstance(leaf, Compare) and not isinstance(leaf.value, str)
        )
        if not numeric:
            raise FilterError(f"{name} holds numbers; {render(leaf)} does not apply")
        return
    raise FilterError(f"{name}: unknown column kind {kind!r}")


def suggestions(value: str, known: frozenset[str]) -> list[str]:
    """The known labels nearest to `value`: the same label in another case first."""
    folded = [k for k in known if k.casefold() == value.casefold()]
    near = difflib.get_close_matches(value, list(known), n=3, cutoff=0.75)
    return list(dict.fromkeys(folded + near))[:3]


def check_labels(leaf: Leaf, known: ValueSets) -> None:
    """Refuse a term whose labels the catalog does not know; see the module docstring."""
    if leaf.column not in known:
        return
    if isinstance(leaf, IsIn | HasAny | HasAll):
        values = [str(v) for v in leaf.values]
    else:
        return
    labels = known[leaf.column]
    unknown = [v for v in values if v not in labels]
    if not unknown:
        return
    near = {v: suggestions(v, labels) for v in unknown}
    hint = "; ".join(
        f"{v!r}: did you mean {' or '.join(repr(s) for s in near[v])}?"
        for v in unknown
        if near[v]
    )
    every = isinstance(leaf, HasAll) or len(unknown) == len(values)
    if every:
        shown = ", ".join(repr(v) for v in unknown[:5])
        raise FilterError(
            f"{leaf.column} has no value {shown} in any dataset of the catalog, so "
            f"{render(leaf)} matches no cell" + (f". {hint}" if hint else "")
        )
    if hint:
        warnings.warn(
            f"{leaf.column}: no dataset of the catalog holds {hint}",
            UserWarning,
            stacklevel=4,
        )


# The catalog test.


def prove(
    node: Expr, row: Mapping[str, Any], tests: dict[Leaf, Callable[[str], bool]]
) -> tuple[Proof, Expr]:
    """The answer of the catalog row for `node`, and the residual filter."""
    if isinstance(node, Const):
        return (ALL if node.value else NONE), node
    if isinstance(node, Not):
        proof, child = prove(node.child, row, tests)
        if proof == SOME:
            return SOME, Not(node.organism, child)
        flipped = NONE if proof == ALL else ALL
        return flipped, Const(node.organism, flipped == ALL)
    if isinstance(node, And | Or):
        decisive, neutral = (NONE, ALL) if isinstance(node, And) else (ALL, NONE)
        rest = []
        for child in node.children:
            proof, residual = prove(child, row, tests)
            if proof == decisive:
                return decisive, Const(node.organism, decisive == ALL)
            if proof == SOME:
                rest.append(residual)
        if not rest:
            return neutral, Const(node.organism, neutral == ALL)
        residual = (
            rest[0] if len(rest) == 1 else node.__class__(node.organism, tuple(rest))
        )
        return SOME, residual
    if not isinstance(node, Leaf):
        raise TypeError(f"not a filter node: {node!r}")
    proof = prove_leaf(node, row, tests)
    return proof, node if proof == SOME else Const(node.organism, proof == ALL)


def prove_leaf(
    leaf: Leaf, row: Mapping[str, Any], tests: dict[Leaf, Callable[[str], bool]]
) -> Proof:
    """The answer of one catalog row for one term."""
    column = leaf.column
    if column == "is_primary_data":
        recorded = row.get("is_primary_data")
        if isinstance(leaf, IsNull):
            return NONE if recorded in ("all_true", "all_false", "mixed") else SOME
        if isinstance(leaf, Truth):
            return {"all_true": ALL, "all_false": NONE}.get(str(recorded), SOME)
        return SOME
    if column == "dataset_id":
        values: Any = [str(row["dataset_id"])]
    elif column in VALUE_SET_COLUMNS and not unrecorded(row.get(column)):
        values = row[column]
    else:
        return SOME
    if isinstance(leaf, IsNull):
        return NONE  # the store format holds no null in this column
    if not is_label_leaf(leaf):
        return SOME
    test = tests.get(leaf)
    if test is None:
        test = tests[leaf] = label_test(leaf)
    matched = np.fromiter((test(str(v)) for v in values), dtype=bool)
    if not matched.any():
        return NONE
    return ALL if matched.all() else SOME
