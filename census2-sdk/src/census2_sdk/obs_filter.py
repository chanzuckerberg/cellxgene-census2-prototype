"""A typed filter builder for `obs`, with pandas names.

A filter starts from the columns of one organism, and combines terms with `&`, `|`,
and `~`, as a pandas mask does:

    obs = census.homo_sapiens.obs
    f = obs.tissue_general.isin(["lung", "blood"]) & obs.is_primary_data
    f &= ~obs.cell_type.str.contains("macrophage", case=False)
    f &= obs.disease.has_any(["COVID-19"]) & (obs.nnz >= 500)
    adata = census.homo_sapiens.get_anndata(obs_filter=f)

Each column has a class for its type (`CategoryColumn`, `BoolColumn`, `IntColumn`,
`FloatColumn`, `StringColumn`, and `MultiValueColumn` for `disease` and
`disease_ontology_term_id`), and each class has only the operators that suit the
type. `census2_sdk.obs_columns` generates the attributes of `obs` from the store
format. `obs["name"]` gives an `AnyColumn`, with every operator, for a column that is
not in the static schema; the query checks it against the stores.

A filter is a value. It holds the terms and the organism of its columns, and `repr`
shows it in the form above. The SDK parses no text form. A filter has no truth value,
so `and`, `or`, and `not` raise `TypeError`; use `&`, `|`, and `~`.

The meaning is the meaning of a pandas mask. A null equals nothing and is in no list,
so each term gives False for a null, and `~` then gives True: `~(obs.x == "a")` keeps the
cells where `x` is null. `isna()` and `notna()` test the null itself. `x != v` is
`~(x == v)`, as in pandas.

`normalize` gives the form that the query evaluates: `~` only on terms, nested `&`
and `|` flattened, constants folded, and terms on one column merged. The catalog test
and the store evaluation (`census2_sdk.obs_plan`, `census2_sdk.obs_arrow`) read that
form.
"""

import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, NoReturn, overload

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
    "Expr",
    "FilterError",
    "FilterLike",
    "FloatColumn",
    "IntColumn",
    "MultiValueColumn",
    "ObsNamespace",
    "Refused",
    "StrMethods",
    "StringColumn",
    "normalize",
]

# The separator of the members of a multi-valued value, as in `catalog.MULTIVALUE_SEP`.
MULTIVALUE_SEP = " || "

Number = int | float
Label = str
Scalar = str | int | float | bool
Inclusive = Literal["both", "neither", "left", "right"]

BOOL_HINT = (
    "a filter has no truth value: combine filters with & | ~, not with and, or, not, "
    "and put each comparison in parentheses, as in (obs.tissue == 'lung') & "
    "obs.is_primary_data"
)


class FilterError(ValueError):
    """A filter that a query refuses: an unknown column or value, or a wrong operator."""


class Refused:
    """The static type of an operator that a column refuses.

    No value of this type is made: the operator raises `TypeError` at run time. The
    type exists so that a type checker reports the misuse where the result is used,
    and the `deprecated` mark makes mypy report the operator itself when the error
    code `deprecated` is on.
    """


# The expression tree. Every node is a frozen dataclass, so equal trees compare equal
# and hash alike. `organism` is the organism of the columns, or "" for a constant.


class Expr:
    """A filter over the `obs` of one organism."""

    organism: str

    def __and__(self, other: "FilterLike") -> "Expr":
        return _combine(And, self, other)

    def __rand__(self, other: "FilterLike") -> "Expr":
        return _combine(And, other, self)

    def __or__(self, other: "FilterLike") -> "Expr":
        return _combine(Or, self, other)

    def __ror__(self, other: "FilterLike") -> "Expr":
        return _combine(Or, other, self)

    def __invert__(self) -> "Expr":
        return Not(self.organism, self)

    def __bool__(self) -> NoReturn:
        raise TypeError(BOOL_HINT)

    def __repr__(self) -> str:
        return render(self)


@dataclass(frozen=True, repr=False)
class Leaf(Expr):
    """A term on one column. `hatch` is true when `obs["name"]` made the column."""

    organism: str
    column: str
    hatch: bool


@dataclass(frozen=True, repr=False)
class Compare(Leaf):
    """`column <op> value`, with op one of == < <= > >=. `!=` is `~(==)`."""

    op: Literal["==", "<", "<=", ">", ">="]
    value: Scalar


@dataclass(frozen=True, repr=False)
class IsIn(Leaf):
    """The value is one of `values`, which are sorted and distinct."""

    values: tuple[Scalar, ...]


@dataclass(frozen=True, repr=False)
class Between(Leaf):
    """`low <= value <= high`, with the ends that `inclusive` names, as in pandas."""

    low: Number
    high: Number
    inclusive: Inclusive


@dataclass(frozen=True, repr=False)
class IsNull(Leaf):
    """The value is null, or NaN, as `isna()` in pandas."""


@dataclass(frozen=True, repr=False)
class StrTest(Leaf):
    """A string method of pandas `.str`, applied to the value."""

    method: Literal["startswith", "endswith", "contains", "match", "fullmatch"]
    pattern: str | tuple[str, ...]
    case: bool = True
    regex: bool = True


@dataclass(frozen=True, repr=False)
class HasAny(Leaf):
    """A member of the multi-valued value is in `values`."""

    values: tuple[str, ...]


@dataclass(frozen=True, repr=False)
class HasAll(Leaf):
    """Every one of `values` is a member of the multi-valued value."""

    values: tuple[str, ...]


@dataclass(frozen=True, repr=False)
class Truth(Leaf):
    """A boolean column used as a filter: the value is true."""


@dataclass(frozen=True, repr=False)
class Const(Expr):
    """True for every cell, or for none."""

    organism: str
    value: bool


@dataclass(frozen=True, repr=False)
class And(Expr):
    organism: str
    children: tuple[Expr, ...]


@dataclass(frozen=True, repr=False)
class Or(Expr):
    organism: str
    children: tuple[Expr, ...]


@dataclass(frozen=True, repr=False)
class Not(Expr):
    organism: str
    child: Expr


TRUE = Const("", True)
FALSE = Const("", False)


def is_const(node: Expr, value: bool) -> bool:
    """True when `node` is the constant `value`, of any organism."""
    return isinstance(node, Const) and node.value is value


def _organism(left: str, right: str) -> str:
    if left and right and left != right:
        raise FilterError(
            f"a filter of {left} and a filter of {right} do not combine; "
            "each filter takes the columns of one organism"
        )
    return left or right


def as_expr(value: object) -> Expr:
    """`value` as a filter: an `Expr`, or a boolean column used as-is."""
    if isinstance(value, Expr):
        return value
    if isinstance(value, BoolColumn):
        return Truth(value.organism, value.name, value.hatch)
    if isinstance(value, bool):
        raise TypeError(
            f"{value!r} is not a filter; did a comparison of two columns give a bool? "
            "Compare a column with a value, as in obs.nnz >= 500"
        )
    raise TypeError(
        f"a filter combines with a filter, not with {type(value).__name__}; "
        "build filters from the columns of census.<organism>.obs"
    )


def _combine(kind: type["And | Or"], left: object, right: object) -> Expr:
    a, b = as_expr(left), as_expr(right)
    return kind(_organism(a.organism, b.organism), (a, b))


# The columns. Each class holds the operators of one type. A class that defines
# `__eq__` gets no hash, as a pandas Series has none.


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


class Column:
    """One `obs` column of one organism, the start of a term."""

    __slots__ = ("hatch", "name", "organism")

    def __init__(self, name: str, organism: str, *, hatch: bool = False) -> None:
        self.name = name
        self.organism = organism
        self.hatch = hatch

    def isna(self) -> Expr:
        """True where the value is null (or NaN), as `Series.isna()`."""
        return IsNull(self.organism, self.name, self.hatch)

    def notna(self) -> Expr:
        """True where the value is not null, as `Series.notna()`."""
        return ~self.isna()

    def __bool__(self) -> NoReturn:
        raise TypeError(BOOL_HINT)

    def __repr__(self) -> str:
        return _name(self.name, self.hatch)

    def _leaf_args(self) -> tuple[str, str, bool]:
        return self.organism, self.name, self.hatch


class StrMethods:
    """The pandas `.str` methods that a filter takes, for a column of labels.

    The query applies each method to each distinct label, with Python's `str` and
    `re`, as pandas does for an object column.
    """

    __slots__ = ("_column",)

    def __init__(self, column: Column) -> None:
        self._column = column

    def startswith(self, pat: str | tuple[str, ...]) -> Expr:
        return self._test("startswith", pat, True, False)

    def endswith(self, pat: str | tuple[str, ...]) -> Expr:
        return self._test("endswith", pat, True, False)

    def contains(self, pat: str, case: bool = True, regex: bool = True) -> Expr:
        """`pat` occurs in the value: a regular expression unless `regex` is false."""
        return self._test("contains", pat, case, regex)

    def match(self, pat: str, case: bool = True) -> Expr:
        """The value starts with a match of the regular expression `pat`."""
        return self._test("match", pat, case, True)

    def fullmatch(self, pat: str, case: bool = True) -> Expr:
        """The whole value matches the regular expression `pat`."""
        return self._test("fullmatch", pat, case, True)

    def _test(
        self,
        method: Literal["startswith", "endswith", "contains", "match", "fullmatch"],
        pat: str | tuple[str, ...],
        case: bool,
        regex: bool,
    ) -> Expr:
        pattern: str | tuple[str, ...]
        if isinstance(pat, tuple):
            if not all(isinstance(p, str) for p in pat):
                raise TypeError(f".str.{method} takes str patterns")
            pattern = tuple(pat)
        elif isinstance(pat, str):
            pattern = pat
        else:
            raise TypeError(f".str.{method} takes a str pattern, not {pat!r}")
        leaf = StrTest(*self._column._leaf_args(), method, pattern, case, regex)
        label_test(leaf)  # a bad regular expression fails here, not at query time
        return leaf


class CategoryColumn(Column):
    """A column of labels, such as `tissue`: `==`, `!=`, `isin`, and `.str`."""

    __slots__ = ()

    @overload  # type: ignore[override]
    def __eq__(self, other: str) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a category column compares with a str")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Expr | Refused:
        return Compare(*self._leaf_args(), "==", _label(other, self.name))

    @overload  # type: ignore[override]
    def __ne__(self, other: str) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a category column compares with a str")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Expr | Refused:
        return ~Compare(*self._leaf_args(), "==", _label(other, self.name))

    def isin(self, values: Iterable[str]) -> Expr:
        """The value is one of `values`, as `Series.isin`."""
        return IsIn(*self._leaf_args(), _labels(values, self.name, "isin"))

    @property
    def str(self) -> StrMethods:
        return StrMethods(self)


class MultiValueColumn(Column):
    """A multi-valued column (`disease`): a cell holds members joined by ' || '.

    Only `has_any` and `has_all` test it. `==`, `!=`, and `isin` raise `TypeError`,
    because a comparison of the joined value is rarely what a query means.
    """

    __slots__ = ()

    def has_any(self, values: Iterable[str]) -> Expr:
        """A member of the cell's value is one of `values`."""
        return HasAny(
            *self._leaf_args(), _labels(values, self.name, "has_any", empty=False)
        )

    def has_all(self, values: Iterable[str]) -> Expr:
        """Each of `values` is a member of the cell's value. Other members can occur."""
        return HasAll(
            *self._leaf_args(), _labels(values, self.name, "has_all", empty=False)
        )

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


class StringColumn(Column):
    """A column of distinct strings (`observation_joinid`): `==`, `!=`, and `isin`."""

    __slots__ = ()

    @overload  # type: ignore[override]
    def __eq__(self, other: str) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a string column compares with a str")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Expr | Refused:
        return Compare(*self._leaf_args(), "==", _label(other, self.name))

    @overload  # type: ignore[override]
    def __ne__(self, other: str) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a string column compares with a str")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Expr | Refused:
        return ~Compare(*self._leaf_args(), "==", _label(other, self.name))

    def isin(self, values: Iterable[str]) -> Expr:
        return IsIn(*self._leaf_args(), _labels(values, self.name, "isin"))


class BoolColumn(Column):
    """A boolean column (`is_primary_data`), which is a filter as it is.

    `obs.is_primary_data` keeps the cells where it is true, and `~obs.is_primary_data`
    the others. It combines with `&` and `|`, and compares with a bool.
    """

    __slots__ = ()

    def __and__(self, other: "FilterLike") -> Expr:
        return _combine(And, self, other)

    def __rand__(self, other: "FilterLike") -> Expr:
        return _combine(And, other, self)

    def __or__(self, other: "FilterLike") -> Expr:
        return _combine(Or, self, other)

    def __ror__(self, other: "FilterLike") -> Expr:
        return _combine(Or, other, self)

    def __invert__(self) -> Expr:
        return ~as_expr(self)

    @overload  # type: ignore[override]
    def __eq__(self, other: bool) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a boolean column compares with a bool")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Expr | Refused:
        if not isinstance(other, bool):
            raise TypeError(f"{self.name} is boolean; compare it with a bool")
        truth = as_expr(self)
        return truth if other else ~truth

    @overload  # type: ignore[override]
    def __ne__(self, other: bool) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a boolean column compares with a bool")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Expr | Refused:
        if not isinstance(other, bool):
            raise TypeError(f"{self.name} is boolean; compare it with a bool")
        return self.__eq__(not other)


def _compare(
    column: Column, op: Literal["==", "<", "<=", ">", ">="], other: object
) -> Expr:
    return Compare(*column._leaf_args(), op, _number(other, column.name))


def _between(column: Column, left: object, right: object, inclusive: str) -> Expr:
    if inclusive not in ("both", "neither", "left", "right"):
        raise ValueError(
            f"inclusive is both, neither, left, or right, not {inclusive!r}"
        )
    low, high = _number(left, column.name), _number(right, column.name)
    return Between(*column._leaf_args(), low, high, inclusive)  # type: ignore[arg-type]


class _Numeric(Column):
    __slots__ = ()

    def _compare(self, op: Literal["==", "<", "<=", ">", ">="], other: object) -> Expr:
        return _compare(self, op, other)

    def __lt__(self, other: Number) -> Expr:
        return self._compare("<", other)

    def __le__(self, other: Number) -> Expr:
        return self._compare("<=", other)

    def __gt__(self, other: Number) -> Expr:
        return self._compare(">", other)

    def __ge__(self, other: Number) -> Expr:
        return self._compare(">=", other)

    @overload  # type: ignore[override]
    def __eq__(self, other: Number) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a numeric column compares with a number")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Expr | Refused:
        return self._compare("==", other)

    @overload  # type: ignore[override]
    def __ne__(self, other: Number) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a numeric column compares with a number")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Expr | Refused:
        return ~self._compare("==", other)

    def between(
        self, left: Number, right: Number, inclusive: Inclusive = "both"
    ) -> Expr:
        """`left <= value <= right`, with the ends of `inclusive`, as `Series.between`."""
        return _between(self, left, right, inclusive)


class IntColumn(_Numeric):
    """An integer column (`nnz`): comparisons and `between`."""

    __slots__ = ()


class FloatColumn(_Numeric):
    """A float column (`raw_sum`): comparisons and `between`. NaN compares false."""

    __slots__ = ()


class AnyColumn(BoolColumn):
    """A column that the static schema does not name: `obs["name"]`.

    It has the operators of every column type. The query checks the column and the
    operator: against the schema when the schema names the column, and otherwise
    against each store, which refuses an operator that does not suit the column's
    type. A store without the column gives no cells, with a warning.
    """

    __slots__ = ()

    @overload  # type: ignore[override]
    def __eq__(self, other: Scalar) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a column compares with a str, a number, or a bool")
    def __eq__(self, other: object) -> Refused: ...
    def __eq__(self, other: object) -> Expr | Refused:
        if isinstance(other, bool):
            return BoolColumn.__eq__(self, other)
        return Compare(*self._leaf_args(), "==", _scalar(other, self.name))

    @overload  # type: ignore[override]
    def __ne__(self, other: Scalar) -> Expr: ...  # type: ignore[overload-overlap]
    @overload
    @deprecated("a column compares with a str, a number, or a bool")
    def __ne__(self, other: object) -> Refused: ...
    def __ne__(self, other: object) -> Expr | Refused:
        if isinstance(other, bool):
            return BoolColumn.__ne__(self, other)
        return ~Compare(*self._leaf_args(), "==", _scalar(other, self.name))

    def __lt__(self, other: Number) -> Expr:
        return _compare(self, "<", other)

    def __le__(self, other: Number) -> Expr:
        return _compare(self, "<=", other)

    def __gt__(self, other: Number) -> Expr:
        return _compare(self, ">", other)

    def __ge__(self, other: Number) -> Expr:
        return _compare(self, ">=", other)

    def between(
        self, left: Number, right: Number, inclusive: Inclusive = "both"
    ) -> Expr:
        return _between(self, left, right, inclusive)

    def isin(self, values: Iterable[str]) -> Expr:
        return IsIn(*self._leaf_args(), _labels(values, self.name, "isin"))

    def has_any(self, values: Iterable[str]) -> Expr:
        return HasAny(
            *self._leaf_args(), _labels(values, self.name, "has_any", empty=False)
        )

    def has_all(self, values: Iterable[str]) -> Expr:
        return HasAll(
            *self._leaf_args(), _labels(values, self.name, "has_all", empty=False)
        )

    @property
    def str(self) -> StrMethods:
        return StrMethods(self)


def _scalar(value: object, column: str) -> Scalar:
    if isinstance(value, str | int | float):
        return value
    raise TypeError(f"{column} compares with a str, a number, or a bool, not {value!r}")


FilterLike = Expr | BoolColumn


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


# Display. `repr` is for a person; the SDK reads no text form of a filter.


def _name(column: str, hatch: bool) -> str:
    return f"obs[{column!r}]" if hatch else f"obs.{column}"


def render(node: Expr) -> str:
    """The pandas form of a filter, for display."""
    if isinstance(node, Const):
        return repr(node.value)
    if isinstance(node, And | Or):
        sep = " & " if isinstance(node, And) else " | "
        return sep.join(_wrapped(c) for c in node.children)
    if isinstance(node, Not):
        child = node.child
        if isinstance(child, IsNull):
            return f"{_name(child.column, child.hatch)}.notna()"
        if isinstance(child, Truth):
            return f"~{_name(child.column, child.hatch)}"
        if isinstance(child, IsIn) and len(child.values) == 1:
            return f"{_name(child.column, child.hatch)} != {child.values[0]!r}"
        return f"~{_wrapped(child)}"
    if not isinstance(node, Leaf):
        raise TypeError(f"not a filter node: {node!r}")
    col = _name(node.column, node.hatch)
    if isinstance(node, Compare):
        return f"{col} {node.op} {node.value!r}"
    if isinstance(node, IsIn):
        if len(node.values) == 1:
            return f"{col} == {node.values[0]!r}"
        return f"{col}.isin({list(node.values)!r})"
    if isinstance(node, Between):
        tail = "" if node.inclusive == "both" else f", inclusive={node.inclusive!r}"
        return f"{col}.between({node.low!r}, {node.high!r}{tail})"
    if isinstance(node, IsNull):
        return f"{col}.isna()"
    if isinstance(node, StrTest):
        args = [repr(node.pattern)]
        if not node.case:
            args.append("case=False")
        if node.method == "contains" and not node.regex:
            args.append("regex=False")
        return f"{col}.str.{node.method}({', '.join(args)})"
    if isinstance(node, HasAny | HasAll):
        kind = "has_any" if isinstance(node, HasAny) else "has_all"
        return f"{col}.{kind}({list(node.values)!r})"
    if isinstance(node, Truth):
        return col
    return f"<{type(node).__name__} of {col}>"


def _wrapped(node: Expr) -> str:
    """`render(node)`, in parentheses unless `&`, `|`, and `~` bind less tightly."""
    text = render(node)
    if isinstance(node, Not):
        if text.startswith("~"):
            return text
        node = node.child
    one = isinstance(node, IsIn) and len(node.values) == 1
    simple = isinstance(node, Truth | Const) or (
        isinstance(node, Leaf) and not isinstance(node, Compare) and not one
    )
    return text if simple else f"({text})"


# The label tests. One function serves the catalog test and the store evaluation, so
# the two cannot disagree on a value.


def members(value: str) -> list[str]:
    """The members of a multi-valued value, as `catalog.members` splits them."""
    return [m.strip() for m in value.split(MULTIVALUE_SEP)]


def label_test(leaf: Leaf) -> Callable[[str], bool]:
    """The test of one label value for a label term. A null is never passed.

    Raises `TypeError` for a term that does not test labels.
    """
    if isinstance(leaf, Compare) and leaf.op == "==" and isinstance(leaf.value, str):
        wanted = leaf.value
        return lambda v: v == wanted
    if isinstance(leaf, IsIn):
        values = frozenset(leaf.values)
        return lambda v: v in values
    if isinstance(leaf, HasAny):
        any_of = frozenset(leaf.values)
        return lambda v: not any_of.isdisjoint(members(v))
    if isinstance(leaf, HasAll):
        all_of = frozenset(leaf.values)
        return lambda v: all_of.issubset(members(v))
    if isinstance(leaf, StrTest):
        return _str_test(leaf)
    raise TypeError(f"{render(leaf)} does not test labels")


def _str_test(leaf: StrTest) -> Callable[[str], bool]:
    pattern = leaf.pattern
    if leaf.method in ("startswith", "endswith"):
        prefixes = pattern if isinstance(pattern, tuple) else (pattern,)
        if leaf.method == "startswith":
            return lambda v: v.startswith(prefixes)
        return lambda v: v.endswith(prefixes)
    if not isinstance(pattern, str):
        raise TypeError(f".str.{leaf.method} takes one str pattern")
    if leaf.method == "contains" and not leaf.regex:
        if leaf.case:
            return lambda v: pattern in v
        upper = pattern.upper()
        return lambda v: upper in v.upper()
    try:
        compiled = re.compile(pattern, 0 if leaf.case else re.IGNORECASE)
    except re.error as error:
        raise FilterError(f"{render(leaf)}: {error}") from error
    if leaf.method == "contains":
        return lambda v: compiled.search(v) is not None
    if leaf.method == "match":
        return lambda v: compiled.match(v) is not None
    return lambda v: compiled.fullmatch(v) is not None


def is_label_leaf(leaf: Leaf) -> bool:
    """True for a term that `label_test` evaluates."""
    if isinstance(leaf, Compare):
        return leaf.op == "==" and isinstance(leaf.value, str)
    return isinstance(leaf, IsIn | HasAny | HasAll | StrTest)


# Normal form.


def normalize(node: Expr) -> Expr:
    """`node` with `~` on terms only, flat `&` and `|`, constants folded, terms merged.

    The pandas rule makes each term two-valued, so De Morgan's laws hold. `between`
    becomes two comparisons. Equality and `isin` terms on one label column merge
    under `|` (a union) and their negations under `&`; `has_any` terms on one column
    merge under `|`, and `has_all` terms under `&`.
    """
    return _simplify(_push(node, negate=False))


def _push(node: Expr, negate: bool) -> Expr:
    if isinstance(node, Not):
        return _push(node.child, not negate)
    if isinstance(node, And | Or):
        kind: type[And | Or] = node.__class__
        if negate:
            kind = Or if isinstance(node, And) else And
        return kind(node.organism, tuple(_push(c, negate) for c in node.children))
    if isinstance(node, Const):
        return Const(node.organism, node.value != negate)
    if isinstance(node, Between):
        lo_op: Literal[">=", ">"] = ">=" if node.inclusive in ("both", "left") else ">"
        hi_op: Literal["<=", "<"] = "<=" if node.inclusive in ("both", "right") else "<"
        args = (node.organism, node.column, node.hatch)
        both = And(
            node.organism,
            (Compare(*args, lo_op, node.low), Compare(*args, hi_op, node.high)),
        )
        return _push(both, negate)
    if isinstance(node, IsIn) and not node.values:
        return Const(node.organism, negate)
    if isinstance(node, HasAny) and not node.values:
        return Const(node.organism, negate)
    if isinstance(node, HasAll) and not node.values:
        return Const(node.organism, not negate)
    if isinstance(node, Compare) and node.op == "==" and isinstance(node.value, str):
        node = IsIn(node.organism, node.column, node.hatch, (node.value,))
    return Not(node.organism, node) if negate else node


def _simplify(node: Expr) -> Expr:
    if not isinstance(node, And | Or):
        return node
    is_and = isinstance(node, And)
    flat: list[Expr] = []
    for child in (_simplify(c) for c in node.children):
        if isinstance(child, node.__class__):
            flat.extend(child.children)
        elif isinstance(child, Const):
            if child.value != is_and:  # False in an And, True in an Or
                return Const(node.organism, child.value)
        else:
            flat.append(child)
    flat = _merge(flat, is_and)
    unique = list(dict.fromkeys(flat))
    if not unique:
        return Const(node.organism, is_and)
    if len(unique) == 1:
        return unique[0]
    return node.__class__(node.organism, tuple(unique))


def _merge(children: list[Expr], is_and: bool) -> list[Expr]:
    """Merge the label terms on one column, keeping the first position of each."""
    out: list[Expr] = []
    slots: dict[tuple[str, str, bool, bool], int] = {}
    for child in children:
        negated = isinstance(child, Not)
        leaf = child.child if isinstance(child, Not) else child
        key: tuple[str, str, bool, bool] | None = None
        if isinstance(leaf, IsIn) and negated == is_and:
            key = ("isin", leaf.column, leaf.hatch, negated)
        elif isinstance(leaf, HasAny) and not negated and not is_and:
            key = ("has_any", leaf.column, leaf.hatch, False)
        elif isinstance(leaf, HasAll) and not negated and is_and:
            key = ("has_all", leaf.column, leaf.hatch, False)
        if key is None:
            out.append(child)
            continue
        if key not in slots:
            slots[key] = len(out)
            out.append(child)
            continue
        previous = out[slots[key]]
        base = previous.child if isinstance(previous, Not) else previous
        assert isinstance(base, IsIn | HasAny | HasAll) and type(base) is type(leaf)
        values = tuple(sorted(set(base.values) | set(leaf.values), key=_order))
        merged: Expr = type(base)(base.organism, base.column, base.hatch, values)  # type: ignore[arg-type]
        out[slots[key]] = Not(merged.organism, merged) if negated else merged
    return out


def _order(value: Scalar) -> tuple[str, str]:
    return type(value).__name__, str(value)


def leaves(node: Expr) -> list[Leaf]:
    """The terms of a filter, in order."""
    if isinstance(node, Leaf):
        return [node]
    if isinstance(node, Not):
        return leaves(node.child)
    if isinstance(node, And | Or):
        return [leaf for c in node.children for leaf in leaves(c)]
    return []


def columns(node: Expr) -> frozenset[str]:
    """The columns that a filter reads."""
    return frozenset(leaf.column for leaf in leaves(node))


def from_mapping(
    obs: ObsNamespace,
    mapping: dict[str, str | Sequence[str]],
    multivalued: Iterable[str],
) -> Expr:
    """The filter of a mapping filter of `get_anndata`, for the columns of `obs`.

    Each column must hold one of its values; a multi-valued column must hold one of
    them as a member. This is the meaning of `catalog.canonicalize` and `_mask`.
    """
    multi = frozenset(multivalued)
    terms: list[Expr] = []
    for column, value in mapping.items():
        values = [value] if isinstance(value, str) else list(value)
        col = getattr(obs, column, None)
        if column in multi:
            multi_col = col if isinstance(col, MultiValueColumn) else obs[column]
            terms.append(multi_col.has_any(values))
        else:
            label_col = col if isinstance(col, CategoryColumn) else obs[column]
            terms.append(label_col.isin(values))
    result: Expr = Const(obs.organism, True)
    for term in terms:
        result = result & term
    return result
