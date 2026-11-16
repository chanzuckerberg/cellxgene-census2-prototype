"""A typed builder of `obs` queries, with pandas names, over pyarrow expressions.

pyarrow is the filter language. Each operator of a column gives a
`pyarrow.compute.Expression`, and `&`, `|`, and `~` combine the expressions. The SDK
evaluates the expression with pyarrow on the columns of each store. This module only
types the operators and keeps two records while it builds a query:

    hs = census.homo_sapiens
    q = hs.obs.cell_type.isin(["neutrophil"]) & (hs.obs.nnz >= 500)
    q.expression  # the pyarrow expression
    q.columns     # {"cell_type", "nnz"}: the columns to read from each store
    q.hints       # {"cell_type": {"neutrophil"}}: the values for the catalog
    q.exact       # False: the nnz term is not a hint

1. **Columns.** Each term adds its column. A pyarrow expression has no public API
   to walk it, so the builder records the columns as it goes.
2. **Hints.** A term `==` or `isin` on a catalog column (`catalog.INDEXED`,
   `dataset_id`), `has_any` on a multi-valued column, and a term on
   `is_primary_data`, give the allowed values of the column. `&` keeps the hints of
   both sides, and intersects the values when both sides give the same single-valued
   column. `|` and `~` keep no hint. `exact` is true when the expression is exactly
   the conjunction of the hints, so that the catalog value sets can prove that every
   cell of a dataset matches.

Each column has a class for its type (`CategoryColumn`, `MultiValueColumn`,
`BoolColumn`, `IntColumn`, `FloatColumn`, `StringColumn`), with only the operators
that suit the type. `census2_sdk.obs_columns` generates the attributes of `obs` from
the store format. `obs["name"]` gives an `AnyColumn` with every operator, for a column
that the static schema does not name.

The null semantics are those of pyarrow, which are SQL's three-valued logic: a
comparison with a null gives null, `~` of null is null, `&` and `|` follow Kleene
logic, and the store keeps only the cells where the expression is true. Thus
`~(obs.x == "a")` and `obs.x != "a"` drop the cells where `x` is null, unlike pandas.
pyarrow's `is_in` gives false for a null, not null, so `isin` here gives null for a
null, as `==` does. Thus every term gives null for a null, which the catalog test
needs (`census2_sdk.value_sets`). The writer stores no null in the columns of the
store format, so the rule matters only for other columns.
"""

import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from functools import reduce
from typing import TYPE_CHECKING, Any, Literal, NoReturn, overload

import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.compute as pc  # type: ignore[import-untyped]

from census2_sdk.catalog import INDEXED, MULTIVALUE_SEP, MULTIVALUED

if TYPE_CHECKING:
    from typing_extensions import deprecated
else:  # the decorator only marks a misuse for a type checker

    def deprecated(message: str) -> Callable[[object], object]:
        return lambda f: f


__all__ = [
    "AnyColumn",
    "BoolColumn",
    "CategoryColumn",
    "Column",
    "FloatColumn",
    "IntColumn",
    "MultiValueColumn",
    "ObsNamespace",
    "Query",
    "QueryLike",
    "Refused",
    "StringColumn",
    "as_query",
]

Number = int | float
Inclusive = Literal["both", "neither", "left", "right"]
# The allowed values of each hinted column: str, or bool for `is_primary_data`.
Hints = Mapping[str, frozenset[Any]]

# The columns that can give a hint.
HINT_COLUMNS = frozenset({"dataset_id", *INDEXED, "is_primary_data"})

BOOL_HINT = (
    "a query has no truth value: combine queries with & | ~, not with and, or, not, "
    "and put each comparison in parentheses, as in (obs.tissue == 'lung') & "
    "obs.is_primary_data"
)

# Display levels, for the parentheses of `repr`, from the strongest binding to the
# weakest, as Python binds them: a call or a column (and `~`), a chain of `&`, a
# chain of `|`, and a comparison.
ATOM, AND, OR, COMPARE = 0, 1, 2, 3


class Refused:
    """The static type of an operator that a column refuses.

    No value of this type is made: the operator raises `TypeError` at run time. The
    type lets a type checker report the misuse where the result is used, and the
    `deprecated` mark lets mypy report the operator itself when the error code
    `deprecated` is on.
    """


@dataclass(frozen=True, eq=False, repr=False)
class Query:
    """A query of the `obs` of one organism.

    `expression` is the pyarrow expression of the selected cells, or None to select
    every cell. `columns` names the columns that the expression reads. `hints` maps
    a column to its allowed values, and `exact` is true when the expression is
    exactly the conjunction of the hints. `text` is the pandas form, for display.
    """

    organism: str
    expression: Any
    columns: frozenset[str] = frozenset()
    hints: Hints = field(default_factory=dict)
    exact: bool = True
    text: str = ""
    level: int = ATOM

    def __and__(self, other: "QueryLike") -> "Query":
        return _and(self, as_query(other))

    def __rand__(self, other: "QueryLike") -> "Query":
        return _and(as_query(other), self)

    def __or__(self, other: "QueryLike") -> "Query":
        return _or(self, as_query(other))

    def __ror__(self, other: "QueryLike") -> "Query":
        return _or(as_query(other), self)

    def __invert__(self) -> "Query":
        if self.expression is None:
            return Query(self.organism, pc.scalar(False), exact=False, text="False")
        return Query(
            self.organism,
            ~self.expression,
            self.columns,
            exact=False,
            text=f"~{_wrap(self, ATOM)}",
        )

    def __bool__(self) -> NoReturn:
        raise TypeError(BOOL_HINT)

    def __repr__(self) -> str:
        return f"<Query of {self.organism}: {self.text or 'every cell'}>"


def _wrap(q: Query, highest: int) -> str:
    """The text of `q`, in parentheses when its level is above `highest`."""
    return f"({q.text})" if q.level > highest else q.text


def _organism(left: Query, right: Query) -> str:
    if left.organism != right.organism:
        raise ValueError(
            f"a query of {left.organism} and a query of {right.organism} do not "
            "combine; each query takes the columns of one organism"
        )
    return left.organism


def _and(left: Query, right: Query) -> Query:
    organism = _organism(left, right)
    if left.expression is None:
        return right
    if right.expression is None:
        return left
    hints = dict(left.hints)
    exact = left.exact and right.exact
    for column, values in right.hints.items():
        if column not in hints:
            hints[column] = values
        elif column in MULTIVALUED:
            # has_any(A) & has_any(B) is not has_any(A & B). Either set is a correct
            # hint alone, so keep the smaller one, and the query is no longer exact.
            hints[column] = min(hints[column], values, key=len)
            exact = False
        else:
            hints[column] = hints[column] & values
    text = f"{_wrap(left, AND)} & {_wrap(right, AND)}"
    return Query(
        organism,
        left.expression & right.expression,
        left.columns | right.columns,
        hints,
        exact,
        text,
        AND,
    )


def _or(left: Query, right: Query) -> Query:
    organism = _organism(left, right)
    if left.expression is None or right.expression is None:
        return Query(organism, None)
    text = f"{_wrap(left, OR)} | {_wrap(right, OR)}"
    return Query(
        organism,
        left.expression | right.expression,
        left.columns | right.columns,
        exact=False,
        text=text,
        level=OR,
    )


def as_query(value: object) -> Query:
    """`value` as a query: a `Query`, a boolean column, or an organism (every cell)."""
    if isinstance(value, Query):
        return value
    if isinstance(value, BoolColumn):
        return value._truth()
    if isinstance(value, HasQuery):
        return value.query()
    if isinstance(value, bool):
        raise TypeError(
            f"{value!r} is not a query; did a comparison of two columns give a bool? "
            "Compare a column with a value, as in obs.nnz >= 500"
        )
    raise TypeError(
        f"a query combines with a query, not with {type(value).__name__}; build "
        "queries from the columns of census.<organism>.obs"
    )


class HasQuery:
    """An object that stands for a query, such as an organism for all of its cells."""

    def query(self) -> Query:
        raise NotImplementedError


# The columns. Each class holds the operators of one type. A class that defines
# `__eq__` has no hash, as a pandas Series has none.


def _labels(
    values: Iterable[object], column: str, kind: str, *, empty: bool = True
) -> tuple[str, ...]:
    if isinstance(values, str):
        raise TypeError(
            f"{column}.{kind} takes a list of values, not one str; "
            f"use {column}.{kind}([{values!r}])"
        )
    out = []
    for v in values:
        if not isinstance(v, str):
            raise TypeError(f"{column}.{kind} takes str values, not {v!r}")
        out.append(v)
    if not out and not empty:
        raise ValueError(f"{column}.{kind} takes one or more values")
    return tuple(sorted(set(out)))


def _label(value: object, column: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{column} holds labels; compare it with a str, not {value!r}")
    return value


def _number(value: object, column: str) -> Number:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise TypeError(
            f"{column} holds numbers; compare it with a number, not {value!r}"
        )
    return value


def _is_in(f: Any, values: tuple[str, ...]) -> Any:
    """`is_in`, with null for a null, as `==` gives (three-valued logic)."""
    found = f.isin(pa.array(values, type=pa.string()))
    return pc.if_else(f.is_valid(), found, pa.scalar(None, type=pa.bool_()))


# RE2 takes a backslash before a punctuation character as that character.
_RE2_SPECIAL = re.compile(r"([\\.+*?()|\[\]{}^$])")


def re2_escape(text: str) -> str:
    """`text` as a literal in an RE2 regular expression."""
    return _RE2_SPECIAL.sub(r"\\\1", text)


def members_pattern(values: Iterable[str]) -> str:
    """An RE2 expression that finds one of `values` as a member of a joined value.

    The members of a cell's value are separated by `MULTIVALUE_SEP`, and the
    catalog strips the space around each member (`catalog.members`).
    """
    sep = re2_escape(MULTIVALUE_SEP)
    either = "|".join(re2_escape(v) for v in values)
    return rf"(?:^|{sep})\s*(?:{either})\s*(?:{sep}|$)"


class Column:
    """One `obs` column of one organism, the start of a term."""

    __slots__ = ("hatch", "name", "organism")

    def __init__(self, name: str, organism: str, *, hatch: bool = False) -> None:
        self.name = name
        self.organism = organism
        self.hatch = hatch

    @property
    def field(self) -> Any:
        return pc.field(self.name)

    def _term(
        self,
        expression: Any,
        text: str,
        level: int = ATOM,
        hint: frozenset[Any] | None = None,
    ) -> Query:
        """A query of one term. A hint applies only to a column of `HINT_COLUMNS`."""
        hints = {self.name: hint} if hint is not None else {}
        if self.name not in HINT_COLUMNS:
            hints = {}
        return Query(
            self.organism,
            expression,
            frozenset({self.name}),
            hints,
            bool(hints),
            text,
            level,
        )

    def __bool__(self) -> NoReturn:
        raise TypeError(BOOL_HINT)

    def __repr__(self) -> str:
        return f"obs[{self.name!r}]" if self.hatch else f"obs.{self.name}"

    # The operators that the column classes share.

    def _eq_label(self, other: object, negate: bool = False) -> Query:
        value = _label(other, self.name)
        if negate:
            return self._term(self.field != value, f"{self!r} != {value!r}", COMPARE)
        return self._term(
            self.field == value, f"{self!r} == {value!r}", COMPARE, frozenset({value})
        )

    def _isin(self, values: Iterable[str]) -> Query:
        labels = _labels(values, self.name, "isin")
        # `obs["name"]` has no type in the catalog test, where its nulls have the null
        # type; pyarrow's `is_in` then needs strings.
        f = self.field.cast(pa.string()) if self.hatch else self.field
        return self._term(
            _is_in(f, labels), f"{self!r}.isin({list(labels)!r})", ATOM,
            frozenset(labels),
        )  # fmt: skip

    def _compare(self, op: str, other: object) -> Query:
        value = _number(other, self.name)
        expression = {
            "==": lambda f: f == value,
            "!=": lambda f: f != value,
            "<": lambda f: f < value,
            "<=": lambda f: f <= value,
            ">": lambda f: f > value,
            ">=": lambda f: f >= value,
        }[op](self.field)
        return self._term(expression, f"{self!r} {op} {value!r}", COMPARE)

    def _between(self, left: object, right: object, inclusive: str) -> Query:
        if inclusive not in ("both", "neither", "left", "right"):
            raise ValueError(
                f"inclusive is both, neither, left, or right, not {inclusive!r}"
            )
        low, high = _number(left, self.name), _number(right, self.name)
        f = self.field
        lower = f >= low if inclusive in ("both", "left") else f > low
        upper = f <= high if inclusive in ("both", "right") else f < high
        args = f"{low!r}, {high!r}" + (
            "" if inclusive == "both" else f", inclusive={inclusive!r}"
        )
        return self._term(lower & upper, f"{self!r}.between({args})")

    def _has_any(self, values: Iterable[str]) -> Query:
        labels = _labels(values, self.name, "has_any", empty=False)
        expression = pc.match_substring_regex(
            self.field.cast(pa.string()), pattern=members_pattern(labels)
        )
        return self._term(
            expression, f"{self!r}.has_any({list(labels)!r})", ATOM, frozenset(labels)
        )

    def _has_all(self, values: Iterable[str]) -> Query:
        labels = _labels(values, self.name, "has_all", empty=False)
        query = reduce(_and, (self._has_any([v]) for v in labels))
        text = f"{self!r}.has_all({list(labels)!r})"
        return Query(
            query.organism,
            query.expression,
            query.columns,
            query.hints,
            query.exact,
            text,
        )


class CategoryColumn(Column):
    """A column of labels, such as `tissue`: `==`, `!=`, and `isin`."""

    __slots__ = ()

    @overload  # type: ignore[override]
    def __eq__(self, other: str) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a category column compares with a str")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Query | Refused:
        return self._eq_label(other)

    @overload  # type: ignore[override]
    def __ne__(self, other: str) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a category column compares with a str")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Query | Refused:
        return self._eq_label(other, negate=True)

    def isin(self, values: Iterable[str]) -> Query:
        """The value is one of `values`, as `Series.isin`."""
        return self._isin(values)


class StringColumn(CategoryColumn):
    """A column of distinct strings (`observation_joinid`), with the same operators."""

    __slots__ = ()


class MultiValueColumn(Column):
    """A multi-valued column (`disease`): a cell holds members joined by ' || '.

    Only `has_any` and `has_all` test it. `==`, `!=`, and `isin` raise `TypeError`,
    because a comparison of the joined value is rarely what a query means.
    """

    __slots__ = ()

    def has_any(self, values: Iterable[str]) -> Query:
        """A member of the cell's value is one of `values`."""
        return self._has_any(values)

    def has_all(self, values: Iterable[str]) -> Query:
        """Each of `values` is a member of the cell's value. Other members can occur."""
        return self._has_all(values)

    def _refuse(self, operator: str) -> NoReturn:
        raise TypeError(
            f"{self.name} is multi-valued: a cell can hold several values joined by "
            f"{MULTIVALUE_SEP!r}, so {operator} is refused. Use "
            f"obs.{self.name}.has_any([...]) for any of the values, or "
            f"obs.{self.name}.has_all([...]) for all of them"
        )

    @deprecated("a multi-valued column takes has_any([...]) or has_all([...])")
    def __eq__(self, other: object) -> Refused:  # type: ignore[override]
        self._refuse("==")

    @deprecated("a multi-valued column takes has_any([...]) or has_all([...])")
    def __ne__(self, other: object) -> Refused:  # type: ignore[override]
        self._refuse("!=")

    @deprecated("a multi-valued column takes has_any([...]) or has_all([...])")
    def isin(self, values: Iterable[str]) -> Refused:
        self._refuse("isin")


class BoolColumn(Column):
    """A boolean column (`is_primary_data`), which is a query as it is.

    `obs.is_primary_data` keeps the cells where it is true, and `~obs.is_primary_data`
    the others. It combines with `&` and `|`, and compares with a bool.
    """

    __slots__ = ()

    def _truth(self, value: bool = True) -> Query:
        if value:
            return self._term(self.field, repr(self), ATOM, frozenset({True}))
        return self._term(~self.field, f"~{self!r}", ATOM, frozenset({False}))

    def __and__(self, other: "QueryLike") -> Query:
        return _and(self._truth(), as_query(other))

    def __rand__(self, other: "QueryLike") -> Query:
        return _and(as_query(other), self._truth())

    def __or__(self, other: "QueryLike") -> Query:
        return _or(self._truth(), as_query(other))

    def __ror__(self, other: "QueryLike") -> Query:
        return _or(as_query(other), self._truth())

    def __invert__(self) -> Query:
        # Not a hint: `~` keeps none, as for every other query.
        truth = self._truth()
        return Query(
            self.organism,
            ~truth.expression,
            truth.columns,
            exact=False,
            text=f"~{self!r}",
        )

    @overload  # type: ignore[override]
    def __eq__(self, other: bool) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a boolean column compares with a bool")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Query | Refused:
        if not isinstance(other, bool):
            raise TypeError(f"{self.name} is boolean; compare it with a bool")
        return self._truth(other)

    @overload  # type: ignore[override]
    def __ne__(self, other: bool) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a boolean column compares with a bool")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Query | Refused:
        if not isinstance(other, bool):
            raise TypeError(f"{self.name} is boolean; compare it with a bool")
        return self._truth(not other)


class _Numeric(Column):
    __slots__ = ()

    def __lt__(self, other: Number) -> Query:
        return self._compare("<", other)

    def __le__(self, other: Number) -> Query:
        return self._compare("<=", other)

    def __gt__(self, other: Number) -> Query:
        return self._compare(">", other)

    def __ge__(self, other: Number) -> Query:
        return self._compare(">=", other)

    @overload  # type: ignore[override]
    def __eq__(self, other: Number) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a numeric column compares with a number")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Query | Refused:
        return self._compare("==", other)

    @overload  # type: ignore[override]
    def __ne__(self, other: Number) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a numeric column compares with a number")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Query | Refused:
        return self._compare("!=", other)

    def between(
        self, left: Number, right: Number, inclusive: Inclusive = "both"
    ) -> Query:
        """`left <= value <= right`, with the ends of `inclusive`, as `Series.between`."""
        return self._between(left, right, inclusive)


class IntColumn(_Numeric):
    """An integer column (`nnz`): comparisons and `between`."""

    __slots__ = ()


class FloatColumn(_Numeric):
    """A float column (`raw_sum`): comparisons and `between`. NaN compares false."""

    __slots__ = ()


class AnyColumn(BoolColumn):
    """A column that the static schema does not name: `obs["name"]`.

    It has the operators of every column type, and pyarrow checks the operator
    against the column's type in each store. A store without the column raises
    `ValueError`.
    """

    __slots__ = ()

    @overload  # type: ignore[override]
    def __eq__(self, other: str | Number) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a column compares with a str, a number, or a bool")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Query | Refused:
        if isinstance(other, bool):
            return BoolColumn.__eq__(self, other)
        if isinstance(other, str):
            return self._eq_label(other)
        return self._compare("==", other)

    @overload  # type: ignore[override]
    def __ne__(self, other: str | Number) -> Query: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a column compares with a str, a number, or a bool")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Query | Refused:
        if isinstance(other, bool):
            return BoolColumn.__ne__(self, other)
        if isinstance(other, str):
            return self._eq_label(other, negate=True)
        return self._compare("!=", other)

    def __lt__(self, other: Number) -> Query:
        return self._compare("<", other)

    def __le__(self, other: Number) -> Query:
        return self._compare("<=", other)

    def __gt__(self, other: Number) -> Query:
        return self._compare(">", other)

    def __ge__(self, other: Number) -> Query:
        return self._compare(">=", other)

    def between(
        self, left: Number, right: Number, inclusive: Inclusive = "both"
    ) -> Query:
        return self._between(left, right, inclusive)

    def isin(self, values: Iterable[str]) -> Query:
        return self._isin(values)

    def has_any(self, values: Iterable[str]) -> Query:
        return self._has_any(values)

    def has_all(self, values: Iterable[str]) -> Query:
        return self._has_all(values)


QueryLike = Query | BoolColumn | HasQuery


class ObsNamespace:
    """The `obs` columns of one organism. `obs_columns.Obs` adds one attribute each.

    `obs["name"]` gives an `AnyColumn` for a column that the static schema does not
    name, such as a column that a later store format adds.
    """

    def __init__(self, organism: str) -> None:
        self.organism = organism

    def __getitem__(self, name: str) -> AnyColumn:
        if not isinstance(name, str):
            raise TypeError(f"obs[...] takes a column name, not {name!r}")
        return AnyColumn(name, self.organism, hatch=True)

    def __repr__(self) -> str:
        return f"<obs of {self.organism}>"
