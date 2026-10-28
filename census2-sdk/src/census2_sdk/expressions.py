"""A query filter as an Iceberg expression, and its push-down to the catalog.

The `obs` table (`census2_sdk.obs_table`) takes a filter in the Iceberg expression
language of pyiceberg: an expression object, or a string that pyiceberg parses, for
example

    tissue IN ('lung', 'lung parenchyma') AND cell_type = 'macrophage'
        AND is_primary_data = true

The language has `AND`, `OR`, `NOT`, comparisons, `IN`, `NOT IN`, `IS NULL`, `IS NOT
NULL`, and `LIKE 'prefix%'`. `from_filter` also takes the SDK's mapping filter and gives
the same expression as a conjunction of `IN` terms.

The query runs in two steps. First `can_match` evaluates the expression against each
catalog row, so the catalog excludes the datasets that cannot hold a matching cell. The
rule for a term is: the dataset can match when some value of the row's value set
satisfies the term. A cell holds one value of each column, so the rule never excludes a
dataset that has a match. `AND` and `OR` combine the results of their terms, and `NOT`
moves to the terms first. A term on a column that the row does not record proves nothing,
so the dataset stays in. Then the scan of the remaining fragments applies the expression
to each cell.

The `disease` columns hold `" || "`-joined values. For these columns, `=`, `!=`, `IN`, and
`NOT IN` mean membership: `disease = 'COVID-19'` selects a cell when one of its values is
`COVID-19`. A literal that holds the separator, as a Census 1 filter gives, selects the
cells with exactly its values. Iceberg has no membership term for a string, so `with_members` rewrites each
such term into `IN` or `NOT IN` over the joined values that the catalog records for the
remaining datasets. Every other term on these columns compares the joined string.
"""

import ast
import math
from collections.abc import Callable, Iterable, Mapping
from typing import Any, cast

from pyiceberg.expressions import (
    AlwaysFalse,
    AlwaysTrue,
    And,
    BooleanExpression,
    EqualTo,
    GreaterThan,
    GreaterThanOrEqual,
    In,
    IsNaN,
    IsNull,
    LessThan,
    LessThanOrEqual,
    Not,
    NotEqualTo,
    NotIn,
    NotNaN,
    NotNull,
    NotStartsWith,
    Or,
    StartsWith,
    UnboundPredicate,
)
from pyiceberg.expressions import parser as _parser
from pyiceberg.expressions.visitors import rewrite_not

from census2_sdk.catalog import (
    INDEXED,
    MULTIVALUE_SEP,
    MULTIVALUED,
    Filter,
    members,
    normalize,
    unrecorded,
)
from census2_sdk.catalog import PRIMARY as PRIMARY_COLUMN

# What `from_filter` takes: the SDK's mapping, an expression string, or an expression.
AnyFilter = Filter | str | BooleanExpression

# The catalog records `is_primary_data` as a tri-state, which is a set of values.
PRIMARY = {"all_true": {True}, "all_false": {False}, "mixed": {True, False}}


def from_filter(obs_filter: AnyFilter | None) -> BooleanExpression:
    """The filter as an unbound Iceberg expression.

    A mapping gives one `IN` term for each column, joined with `AND`. It takes the
    columns that `catalog.normalize` takes.
    """
    if obs_filter is None:
        return AlwaysTrue()
    if isinstance(obs_filter, BooleanExpression):
        return obs_filter
    if isinstance(obs_filter, str):
        return _parser.parse(obs_filter)
    terms: list[BooleanExpression] = []
    for column, wanted in normalize(obs_filter).items():
        if column == PRIMARY_COLUMN:
            terms.append(is_in(column, sorted(w == "True" for w in wanted)))
        else:
            terms.append(is_in(column, sorted(wanted)))
    return _all(terms)


# pyiceberg's terms take a column name and literals at run time, but their pydantic
# signatures name other fields, so a type checker refuses the call. These functions
# state the run-time form once.
def is_in(column: str, values: Iterable[Any]) -> BooleanExpression:
    """`column IN (values)`. One value gives `=`, and no value gives `false`."""
    return cast("BooleanExpression", cast("Any", In)(column, list(values)))


def not_in(column: str, values: Iterable[Any]) -> BooleanExpression:
    """`column NOT IN (values)`."""
    return cast("BooleanExpression", cast("Any", NotIn)(column, list(values)))


def not_null(column: str) -> BooleanExpression:
    """`column IS NOT NULL`."""
    return cast("BooleanExpression", cast("Any", NotNull)(column))


def _text_literal(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    return str(value)


_OPERATORS = {
    EqualTo: "=",
    NotEqualTo: "!=",
    LessThan: "<",
    LessThanOrEqual: "<=",
    GreaterThan: ">",
    GreaterThanOrEqual: ">=",
}


def to_text(expr: BooleanExpression) -> str:
    """The expression in the string form that `from_filter` parses."""
    if isinstance(expr, AlwaysTrue):
        return "true"
    if isinstance(expr, AlwaysFalse):
        return "false"
    if isinstance(expr, And):
        return f"({to_text(expr.left)} AND {to_text(expr.right)})"
    if isinstance(expr, Or):
        return f"({to_text(expr.left)} OR {to_text(expr.right)})"
    if isinstance(expr, Not):
        return f"NOT {to_text(expr.child)}"
    if not isinstance(expr, UnboundPredicate):
        raise TypeError(f"no text form for {expr!r}")
    column = str(expr.term.name)  # type: ignore[attr-defined]
    if isinstance(expr, (In, NotIn)):
        values = ", ".join(_text_literal(v) for v in sorted(_literals(expr)))
        return f"{column} {'IN' if isinstance(expr, In) else 'NOT IN'} ({values})"
    if isinstance(expr, (IsNull, NotNull)):
        return f"{column} IS {'' if isinstance(expr, IsNull) else 'NOT '}NULL"
    if isinstance(expr, (StartsWith, NotStartsWith)):
        prefix = _text_literal(str(_literal(expr)) + "%")
        return f"{column} {'' if isinstance(expr, StartsWith) else 'NOT '}LIKE {prefix}"
    operator = _OPERATORS.get(type(expr))
    if operator is None:
        raise ValueError(f"no text form for {expr!r}")
    return f"{column} {operator} {_text_literal(_literal(expr))}"


def from_value_filter(value_filter: str) -> BooleanExpression:
    """The Iceberg expression of a TileDB-SOMA `value_filter` string of Census 1.

    The translator takes `and`, `or`, `not`, the six comparisons, `in [...]`, and
    `not in [...]`, with string, number, and boolean literals. It raises ValueError
    for any other construct, for example a comparison of two columns. A migration
    error can then name the equivalent new filter.

    The translation keeps the meaning of each term except one: on a multi-valued
    column, Census 1 compares the joined string, and this SDK tests membership (see
    the module docstring).
    """
    try:
        tree = ast.parse(value_filter.strip(), mode="eval").body
    except SyntaxError as error:
        raise ValueError(f"not a value filter: {value_filter!r}") from error
    return _from_node(tree)


_COMPARISONS: dict[type, type] = {
    ast.Eq: EqualTo,
    ast.NotEq: NotEqualTo,
    ast.Lt: LessThan,
    ast.LtE: LessThanOrEqual,
    ast.Gt: GreaterThan,
    ast.GtE: GreaterThanOrEqual,
}


def _from_node(node: ast.expr) -> BooleanExpression:
    if isinstance(node, ast.BoolOp):
        parts = [_from_node(v) for v in node.values]
        joined = And if isinstance(node.op, ast.And) else Or
        return cast("BooleanExpression", joined(*parts))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return Not(_from_node(node.operand))
    if isinstance(node, ast.Compare):
        terms: list[BooleanExpression] = []
        left = node.left
        for op, right in zip(node.ops, node.comparators, strict=True):
            terms.append(_comparison(left, op, right))
            left = right
        return _all(terms)
    raise ValueError(f"no Iceberg term for {ast.unparse(node)!r}")


def _comparison(left: ast.expr, op: ast.cmpop, right: ast.expr) -> BooleanExpression:
    # Each refusal is a ValueError: the input is a filter string, whatever the node.
    column = left.id if isinstance(left, ast.Name) else None
    if column is None:
        raise ValueError(f"a term starts with a column: {ast.unparse(left)!r}")
    try:
        value = ast.literal_eval(right)
    except ValueError as error:
        raise ValueError(f"no Iceberg literal for {ast.unparse(right)!r}") from error
    if isinstance(op, (ast.In, ast.NotIn)):
        values = value if isinstance(value, (list, tuple, set)) else None
        if values is None:
            raise ValueError(f"`in` takes a list: {ast.unparse(right)!r}")
        maker = is_in if isinstance(op, ast.In) else not_in
        return maker(column, values)
    kind = _COMPARISONS.get(type(op))
    if kind is None or isinstance(value, (list, tuple, set, dict)):
        raise ValueError(f"no Iceberg term for {ast.unparse(op)!r}")
    return cast("BooleanExpression", cast("Any", kind)(column, value))


def _all(terms: list[BooleanExpression]) -> BooleanExpression:
    if not terms:
        return AlwaysTrue()
    return terms[0] if len(terms) == 1 else And(*terms)


def columns(expr: BooleanExpression) -> set[str]:
    """The columns that the expression names."""
    if isinstance(expr, (And, Or)):
        return columns(expr.left) | columns(expr.right)
    if isinstance(expr, Not):
        return columns(expr.child)
    if isinstance(expr, UnboundPredicate):
        return {str(expr.term.name)}  # type: ignore[attr-defined]
    return set()


def prefix_columns(expr: BooleanExpression) -> set[str]:
    """The columns of the `LIKE 'prefix%'` terms of the expression."""
    if isinstance(expr, (And, Or)):
        return prefix_columns(expr.left) | prefix_columns(expr.right)
    if isinstance(expr, Not):
        return prefix_columns(expr.child)
    if isinstance(expr, (StartsWith, NotStartsWith)):
        return {str(expr.term.name)}  # type: ignore[attr-defined]
    return set()


def _literal(pred: Any) -> Any:
    return pred.literal.value


def _literals(pred: Any) -> set[Any]:
    return {lit.value for lit in pred.literals}


# The rule of each term for one value. A comparison of values of two types raises
# TypeError, and the caller then keeps the dataset.
_HOLDS: dict[type, Callable[[Any, Any], bool]] = {
    EqualTo: lambda p, v: bool(v == _literal(p)),
    NotEqualTo: lambda p, v: bool(v != _literal(p)),
    In: lambda p, v: v in _literals(p),
    NotIn: lambda p, v: v not in _literals(p),
    LessThan: lambda p, v: bool(v < _literal(p)),
    LessThanOrEqual: lambda p, v: bool(v <= _literal(p)),
    GreaterThan: lambda p, v: bool(v > _literal(p)),
    GreaterThanOrEqual: lambda p, v: bool(v >= _literal(p)),
    StartsWith: lambda p, v: str(v).startswith(str(_literal(p))),
    NotStartsWith: lambda p, v: not str(v).startswith(str(_literal(p))),
    NotNull: lambda p, v: v is not None,
    NotNaN: lambda p, v: not (isinstance(v, float) and math.isnan(v)),
}


# The membership rule of each term on a multi-valued column, for one joined value.
def _member(literal: Any, held: set[str]) -> bool:
    """True when a cell with the members `held` matches one literal.

    A bare term matches when it is a member. A joined literal, as a Census 1 filter
    gives, matches the cells with exactly its members.
    """
    text = str(literal)
    separator = MULTIVALUE_SEP.strip()
    if separator in text:
        return held == {m.strip() for m in text.split(separator)}
    return text in held


_HOLDS_MEMBER: dict[type, Callable[[Any, set[str]], bool]] = {
    EqualTo: lambda p, m: _member(_literal(p), m),
    NotEqualTo: lambda p, m: not _member(_literal(p), m),
    In: lambda p, m: any(_member(v, m) for v in _literals(p)),
    NotIn: lambda p, m: not any(_member(v, m) for v in _literals(p)),
}


def holds(pred: UnboundPredicate, value: Any, column: str) -> bool:
    """True when one cell value satisfies the term. A null satisfies no term here."""
    if value is None:
        return False
    member_rule = _HOLDS_MEMBER.get(type(pred)) if column in MULTIVALUED else None
    if member_rule is not None:
        return member_rule(pred, set(members(str(value), column)))
    rule = _HOLDS.get(type(pred))
    if rule is None:
        return True  # a term that this module does not know proves nothing
    return rule(pred, value)


def value_set(row: Mapping[str, Any], column: str) -> set[Any] | None:
    """The values of `column` that the catalog row records, or None when it records none."""
    if column == "dataset_id":
        return {str(row["dataset_id"])}
    if column == PRIMARY_COLUMN:
        state = row.get(PRIMARY_COLUMN)
        return None if unrecorded(state) else PRIMARY.get(str(state))
    if column not in INDEXED:
        return None
    values = row.get(column)
    if unrecorded(values):
        return None
    return {str(v) for v in cast("Iterable[Any]", values)}


def can_match(expr: BooleanExpression, row: Mapping[str, Any]) -> bool:
    """False when the catalog row proves that no cell of the dataset satisfies `expr`.

    The row records the value set of a column, without nulls. Thus `IS NULL` and
    `IS NAN` prove nothing, and a term on an unrecorded column proves nothing.
    """
    return _can_match(rewrite_not(expr), row)


def _can_match(expr: BooleanExpression, row: Mapping[str, Any]) -> bool:
    if isinstance(expr, AlwaysFalse):
        return False
    if isinstance(expr, And):
        return _can_match(expr.left, row) and _can_match(expr.right, row)
    if isinstance(expr, Or):
        return _can_match(expr.left, row) or _can_match(expr.right, row)
    if not isinstance(expr, UnboundPredicate) or isinstance(expr, (IsNull, IsNaN)):
        return True
    column = str(expr.term.name)  # type: ignore[attr-defined]
    values = value_set(row, column)
    if values is None:
        return True
    try:
        return any(holds(expr, v, column) for v in values)
    except TypeError:
        return True  # the term compares another type; it proves nothing


def with_members(
    expr: BooleanExpression, rows: Iterable[Mapping[str, Any]]
) -> BooleanExpression:
    """`expr` with each membership term on a multi-valued column as a string term.

    `rows` are the catalog rows that the scan reads. The joined values of a column are
    the values that those rows record. A term becomes `IN` or `NOT IN` over the joined
    values that hold, or do not hold, a wanted member.
    """
    expr = rewrite_not(expr)
    named = columns(expr) & MULTIVALUED
    if not named:
        return expr
    joined: dict[str, set[str]] = {c: set() for c in named}
    for row in rows:
        for column in named:
            joined[column] |= value_set(row, column) or set()
    return _with_members(expr, joined)


def _with_members(
    expr: BooleanExpression, joined: Mapping[str, set[str]]
) -> BooleanExpression:
    if isinstance(expr, And):
        return And(_with_members(expr.left, joined), _with_members(expr.right, joined))
    if isinstance(expr, Or):
        return Or(_with_members(expr.left, joined), _with_members(expr.right, joined))
    if not isinstance(expr, UnboundPredicate):
        return expr
    column = str(expr.term.name)  # type: ignore[attr-defined]
    if column not in joined:
        return expr
    if isinstance(expr, (EqualTo, In)):
        hits = sorted(j for j in joined[column] if holds(expr, j, column))
        return is_in(column, hits)
    if isinstance(expr, (NotEqualTo, NotIn)):
        # The joined values that hold a member that the term refuses.
        wanted = _literals(expr) if isinstance(expr, NotIn) else {_literal(expr)}
        positive = cast("UnboundPredicate", is_in(column, sorted(wanted)))
        refused = sorted(j for j in joined[column] if holds(positive, j, column))
        return not_in(column, refused) if refused else not_null(column)
    return expr
