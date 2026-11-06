"""The `obs` query in SQL: DuckDB over the `obs` fragments, without Iceberg.

The writer puts one Parquet fragment next to each store (`census2_sdk.fragments`). The
catalog lists the stores of a release, so it also lists their fragments. On this path,
DuckDB reads the fragments directly:

1. The catalog excludes the datasets that cannot match (`catalog_test`).
2. One DuckDB query reads the fragments of the remaining datasets and selects the cells
   of the filter. It returns the cell key, the `obs` columns, and the `X` range of each
   selected cell.
3. The SDK opens only the stores that hold a selected cell, and reads only the `X`
   values of those cells, as on the Iceberg path.

    census = census2_sdk.open_census(["catalog.parquet"], obs_engine="duckdb")
    adata = census.get_anndata(
        "homo_sapiens", "cell_type = 'neutrophil' AND nnz >= 500", max_workers=8
    )

There is no table metadata. The catalog pins the dataset versions of a release, because
each catalog row names one store and thus one fragment.

The filter language
-------------------
A filter is one SQL boolean expression over the `obs` columns: the part of a query after
`WHERE`. For example:

    tissue IN ('lung', 'lung parenchyma') AND cell_type = 'macrophage'
        AND is_primary_data AND nnz >= 500

`parse` reads the text with DuckDB's own parser (`json_serialize_sql`). That function
parses a statement and does not run it. `parse` then accepts only these parts of the
parse tree:

- `AND`, `OR`, `NOT`, the comparisons, `[NOT] BETWEEN`, `[NOT] IN (...)`,
  `IS [NOT] NULL`, `IS [NOT] DISTINCT FROM`, `[NOT] LIKE`, `[NOT] ILIKE`, `GLOB`,
  `SIMILAR TO`, and arithmetic;
- string, number, boolean, and null literals, list literals, and `CAST` to a scalar
  type;
- column names without a table name;
- the scalar functions in `FUNCTIONS`: string, list, and number functions that read no
  state;
- one form of subquery, for ontology closure (see below).

Thus a filter cannot name a table, a table function, a file, a setting, a sequence, or
an environment variable, and it cannot add a clause to the query. `render` writes the
SQL again from the accepted tree, and the SDK runs that text, not the text of the user.
As a second barrier, the DuckDB connection of the SDK refuses file access outside the
directories of the fragments (`enable_external_access = false` with
`allowed_directories`), and it locks its configuration before the first query.

A filter means what SQL means. In particular, `disease = 'COVID-19'` compares the whole
joined value, as Census 1 does. Membership in a multi-valued column is explicit:

    list_contains(string_split(disease, ' || '), 'COVID-19')
    list_has_any(string_split(disease_ontology_term_id, ' || '), ['MONDO:0100096'])

`from_mapping` gives the SQL of the SDK's mapping filter. On a multi-valued column it
gives this explicit membership.

Ontology closure
----------------
A filter can select the descendants of an ontology term with a semi-join on the
`descendants` table, which has the columns `ontology`, `ancestor`, and `descendant`
(each term is its own descendant):

    cell_type_ontology_term_id IN (
        SELECT descendant FROM descendants WHERE ancestor = 'CL:0000084')

This is the one subquery that the language accepts: `<expression> [NOT] IN (SELECT
descendant FROM descendants [WHERE ...])`, where the `WHERE` part follows the same rules
over the three columns. `Options.descendants` names the Parquet file of the table
(`tools/build_descendants.py`). Before the query, the SDK runs each such subquery and
puts the list of terms in its place, so the catalog test and the Parquet statistics see
a plain `IN` list.

Push-down
---------
1. The catalog. `catalog_test` gives a SQL test over the catalog rows. The filter goes
   to negation normal form first, so that each `NOT` applies to one term. A term on one
   column with a value set in the catalog becomes
   `list_bool_or(list_transform(<values>, lambda v: <term on v>))`. DuckDB thus
   evaluates the term with the same functions as the scan, and the test needs no
   second implementation of each function. The values include a null, because the
   value set records no nulls. A term on two columns, or on a column without a value
   set, proves nothing. The test is exact for one term and conservative for a
   conjunction, as `expressions.can_match` is.
2. Parquet. DuckDB pushes comparisons, `IN`, and `IS NULL` into the scan of each
   fragment, and skips the row groups whose statistics exclude them.
3. The rows. DuckDB applies the whole filter to each row.

`duckdb` is the optional `duckdb` extra. A remote fragment needs DuckDB's `httpfs`
extension, which DuckDB downloads on first use unless it is already installed.
"""

import ast
import json
import os
import pathlib
import re
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, cast

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.compute as pc  # type: ignore[import-untyped]

from census2_sdk import catalog as _catalog
from census2_sdk import fragments
from census2_sdk.catalog import INDEXED, MULTIVALUE_SEP, MULTIVALUED, PRIMARY, Filter

if TYPE_CHECKING:
    import anndata

    from census2_sdk import Census

Node = dict[str, Any]


class FilterError(ValueError):
    """A filter that the SQL filter language does not accept."""


# The scalar functions that a filter can call. Each one reads only its arguments: no
# file, setting, sequence, clock, random state, or environment variable.
FUNCTIONS = frozenset(
    {
        # strings
        "lower",
        "upper",
        "lcase",
        "ucase",
        "trim",
        "ltrim",
        "rtrim",
        "length",
        "strlen",
        "starts_with",
        "ends_with",
        "prefix",
        "suffix",
        "contains",
        "strpos",
        "instr",
        "substring",
        "substr",
        "left",
        "right",
        "replace",
        "concat",
        "split_part",
        "string_split",
        "str_split",
        "string_split_regex",
        "regexp_split_to_array",
        "regexp_matches",
        "regexp_full_match",
        "regexp_extract",
        # lists
        "list_value",
        "list_contains",
        "list_has",
        "array_contains",
        "array_has",
        "list_has_any",
        "array_has_any",
        "list_has_all",
        "array_has_all",
        "len",
        "array_length",
        # numbers
        "abs",
        "round",
        "floor",
        "ceil",
        "ln",
        "log",
        "log2",
        "log10",
        "sqrt",
        "exp",
        "pow",
        "power",
        "greatest",
        "least",
        "isnan",
        "isinf",
        "isfinite",
        # nulls
        "coalesce",
        "ifnull",
        "nullif",
    }
)

# DuckDB's parser gives these operators as functions.
_INFIX = {
    "~~": "LIKE",
    "!~~": "NOT LIKE",
    "~~*": "ILIKE",
    "!~~*": "NOT ILIKE",
    "~~~": "GLOB",
    "+": "+",
    "-": "-",
    "*": "*",
    "/": "/",
    "//": "//",
    "%": "%",
}

_COMPARISONS = {
    "COMPARE_EQUAL": "=",
    "COMPARE_NOTEQUAL": "<>",
    "COMPARE_LESSTHAN": "<",
    "COMPARE_GREATERTHAN": ">",
    "COMPARE_LESSTHANOREQUALTO": "<=",
    "COMPARE_GREATERTHANOREQUALTO": ">=",
    "COMPARE_DISTINCT_FROM": "IS DISTINCT FROM",
    "COMPARE_NOT_DISTINCT_FROM": "IS NOT DISTINCT FROM",
}

_CASTS = frozenset(
    {"VARCHAR", "BOOLEAN", "INTEGER", "BIGINT", "DOUBLE", "FLOAT", "REAL"}
)

_NUMBER = re.compile(r"-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?")

# The table of ontology closure, and its columns.
DESCENDANTS = "descendants"
DESCENDANT_COLUMNS = ("ontology", "ancestor", "descendant")

# The longest filter that `parse` reads. A longer one is not a hand-written filter.
MAX_FILTER_CHARS = 1_000_000

# The columns that the scan always returns: the cell key and the fragment's row
# columns. `get_anndata` needs them to find each cell in its store.
KEY = ("dataset_id", "observation_joinid")


@dataclass(frozen=True)
class Predicate:
    """An accepted filter: its parse tree, and its SQL as the SDK runs it."""

    tree: Node

    @property
    def text(self) -> str:
        return render(self.tree)

    @property
    def columns(self) -> frozenset[str]:
        return frozenset(_columns(self.tree))

    def has_subquery(self) -> bool:
        return _any(self.tree, lambda n: n.get("class") == "SUBQUERY")


def parse(text: str) -> Predicate:
    """Read a filter, and refuse it unless the language accepts every part of it.

    Raise `FilterError` with the reason.
    """
    if len(text) > MAX_FILTER_CHARS:
        raise FilterError(f"a filter has at most {MAX_FILTER_CHARS} characters")
    statement = f"SELECT 1 FROM obs WHERE {text}"
    con = _parser()
    found = con.execute("SELECT json_serialize_sql(?)", [statement]).fetchone()
    tree = json.loads(str(found[0] if found else "{}"))
    if tree.get("error"):
        raise FilterError(f"not a SQL expression: {tree.get('error_message')}")
    statements = tree.get("statements", [])
    if len(statements) != 1:
        raise FilterError("a filter is one expression")
    node = statements[0]["node"]
    _check_select(node, table="obs", select=None)
    where = node["where_clause"]
    _check(where, set())
    return Predicate(where)


_PARSER: duckdb.DuckDBPyConnection | None = None


def _parser() -> duckdb.DuckDBPyConnection:
    """A connection that only parses. It reads no file and holds no table."""
    global _PARSER
    if _PARSER is None:
        _PARSER = duckdb.connect(config={"enable_external_access": "false"})
    return _PARSER


def _check_select(node: Node, *, table: str, select: str | None) -> None:
    """The statement is `SELECT 1 FROM <table> WHERE ...` and nothing more.

    With `select`, the select list is that one column instead of the constant 1. Any
    other clause means that the filter text broke out of the `WHERE` part.
    """
    if node.get("type") != "SELECT_NODE":
        raise FilterError("a filter is one expression, not a query")
    empty = ("modifiers", "group_expressions", "group_sets")
    if any(node.get(k) for k in empty) or node.get("cte_map", {}).get("map"):
        raise FilterError("a filter cannot add a clause to the query")
    if any(node.get(k) is not None for k in ("having", "qualify", "sample")):
        raise FilterError("a filter cannot add a clause to the query")
    if node.get("aggregate_handling") != "STANDARD_HANDLING":
        raise FilterError("a filter cannot add a clause to the query")
    source = node.get("from_table") or {}
    plain = (
        source.get("type") == "BASE_TABLE"
        and source.get("table_name") == table
        and not source.get("schema_name")
        and not source.get("catalog_name")
        and not source.get("column_name_alias")
        and source.get("sample") is None
        and source.get("at_clause") is None
    )
    if not plain:
        raise FilterError(f"a filter reads the {table} table only")
    items = node.get("select_list") or []
    item = items[0] if len(items) == 1 else {}
    if select is None:
        wanted = item.get("class") == "CONSTANT"
        wanted = wanted and item.get("literal") == {"kind": "INTEGER", "text": "1"}
    else:
        wanted = item.get("class") == "COLUMN_REF"
        wanted = wanted and item.get("column_names") == [select]
    if not wanted or item.get("alias"):
        raise FilterError("a filter cannot change the select list")
    if node.get("where_clause") is None and select is None:
        raise FilterError("a filter is an expression")


def _check(node: Node, lambdas: set[str]) -> None:
    """Refuse any part of the tree that the language does not accept."""
    if not isinstance(node, dict):
        raise FilterError("not a SQL expression")
    if node.get("alias"):
        raise FilterError("a filter has no alias")
    kind, cls = node.get("type"), node.get("class")
    if (
        cls == "CONJUNCTION"
        and kind in ("CONJUNCTION_AND", "CONJUNCTION_OR")
        or cls == "OPERATOR"
        and kind
        in (
            "OPERATOR_NOT",
            "COMPARE_IN",
            "COMPARE_NOT_IN",
            "OPERATOR_IS_NULL",
            "OPERATOR_IS_NOT_NULL",
            "OPERATOR_COALESCE",
        )
    ):
        for child in node["children"]:
            _check(child, lambdas)
    elif cls == "COMPARISON" and kind in _COMPARISONS:
        _check(node["left"], lambdas)
        _check(node["right"], lambdas)
    elif cls == "BETWEEN":
        for key in ("input", "lower", "upper"):
            _check(node[key], lambdas)
    elif cls == "CONSTANT":
        literal = node.get("literal") or {}
        kind = literal.get("kind")
        text = str(literal.get("text", ""))
        if kind in ("INTEGER", "NUMERIC") and not _NUMBER.fullmatch(text):
            raise FilterError(f"not a number: {text!r}")
        if kind == "BOOLEAN" and text not in ("true", "false"):
            raise FilterError(f"not a boolean: {text!r}")
        if kind not in ("STRING", "INTEGER", "NUMERIC", "BOOLEAN", "NULL_LITERAL"):
            raise FilterError(f"no literal of kind {kind}")
    elif cls == "COLUMN_REF":
        names = node.get("column_names") or []
        if len(names) != 1:
            raise FilterError(f"a column name has no table name: {'.'.join(names)}")
    elif cls == "FUNCTION":
        _check_function(node, lambdas)
    elif cls == "CAST":
        type_expr = node.get("type_expr") or {}
        name = str(type_expr.get("type_name", "")).upper()
        if name not in _CASTS or type_expr.get("children"):
            raise FilterError(f"no cast to {name or 'that type'}")
        _check(node["child"], lambdas)
    elif cls == "SUBQUERY":
        _check_descendants(node, lambdas)
    else:
        raise FilterError(f"the filter language has no {str(kind or cls).lower()}")


def _check_function(node: Node, lambdas: set[str]) -> None:
    name = str(node.get("function_name", ""))
    if name not in FUNCTIONS and name not in _INFIX:
        raise FilterError(f"the filter language has no function {name}")
    if node.get("schema") or node.get("catalog"):
        raise FilterError(f"a function name has no schema: {name}")
    if node.get("distinct") or node.get("filter") or node.get("export_state"):
        raise FilterError(f"{name} is not an aggregate")
    if (node.get("order_bys") or {}).get("orders"):
        raise FilterError(f"{name} is not an aggregate")
    arguments = node.get("arguments") or []
    if name in _INFIX and len(arguments) not in (1, 2):
        raise FilterError(f"{name} takes two operands")
    if name in _INFIX and len(arguments) == 1 and name != "-":
        raise FilterError(f"{name} takes two operands")
    for argument in arguments:
        if argument.get("name"):
            raise FilterError(f"{name} takes no named argument")
        _check(argument["expression"], lambdas)


def _check_descendants(node: Node, lambdas: set[str]) -> None:
    """`<expression> IN (SELECT descendant FROM descendants [WHERE ...])` only."""
    shape = node.get("subquery_type") == "ANY"
    shape = shape and node.get("comparison_type") == "COMPARE_EQUAL"
    if not shape or node.get("child") is None:
        raise FilterError(
            "the one subquery is <column> IN (SELECT descendant FROM descendants ...)"
        )
    _check(node["child"], lambdas)
    inner = node["subquery"]
    if inner.get("named_param_map"):
        raise FilterError("a filter has no parameters")
    select = inner["node"]
    _check_select(select, table=DESCENDANTS, select="descendant")
    where = select.get("where_clause")
    if where is not None:
        if _any(where, lambda n: n.get("class") == "SUBQUERY"):
            raise FilterError("a subquery of descendants has no subquery")
        _check(where, lambdas)
        unknown = _columns(where) - set(DESCENDANT_COLUMNS)
        if unknown:
            raise FilterError(f"descendants has no column {sorted(unknown)}")


def _any(node: Any, test: Callable[[Node], bool]) -> bool:
    """True when `test` holds for some node of the tree."""
    if isinstance(node, dict):
        if test(node):
            return True
        return any(_any(v, test) for v in node.values())
    if isinstance(node, list):
        return any(_any(v, test) for v in node)
    return False


def _children(node: Node) -> list[Node]:
    """The expression children of one accepted node, outside subqueries."""
    cls = node.get("class")
    if cls in ("CONJUNCTION", "OPERATOR"):
        return list(node["children"])
    if cls == "COMPARISON":
        return [node["left"], node["right"]]
    if cls == "BETWEEN":
        return [node["input"], node["lower"], node["upper"]]
    if cls == "FUNCTION":
        return [a["expression"] for a in node.get("arguments") or []]
    if cls == "CAST":
        return [node["child"]]
    if cls == "SUBQUERY":
        return [node["child"]]
    return []


def _columns(node: Node) -> set[str]:
    """The columns that an accepted expression names, outside subqueries."""
    if node.get("class") == "COLUMN_REF":
        return {str(node["column_names"][0])}
    return set().union(*(_columns(c) for c in _children(node)))


def quote(value: str) -> str:
    """A SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def identifier(name: str) -> str:
    """A SQL identifier, quoted."""
    return '"' + name.replace('"', '""') + '"'


def render(node: Node, column: Callable[[str], str] = identifier) -> str:
    """The SQL of an accepted expression. Every operation has its own parentheses.

    `column` gives the SQL of each column name; the catalog test puts a lambda
    parameter in place of one column.
    """

    def r(n: Node) -> str:
        return render(n, column)

    cls, kind = node["class"], node["type"]
    if cls == "CONJUNCTION":
        joiner = " AND " if kind == "CONJUNCTION_AND" else " OR "
        return "(" + joiner.join(r(c) for c in node["children"]) + ")"
    if cls == "OPERATOR":
        children = node["children"]
        if kind == "OPERATOR_NOT":
            return f"(NOT {r(children[0])})"
        if kind == "OPERATOR_COALESCE":
            return f"coalesce({', '.join(r(c) for c in children)})"
        if kind in ("COMPARE_IN", "COMPARE_NOT_IN"):
            word = "IN" if kind == "COMPARE_IN" else "NOT IN"
            values = ", ".join(r(c) for c in children[1:])
            return f"({r(children[0])} {word} ({values}))"
        word = "IS NULL" if kind == "OPERATOR_IS_NULL" else "IS NOT NULL"
        return f"({r(children[0])} {word})"
    if cls == "COMPARISON":
        return f"({r(node['left'])} {_COMPARISONS[kind]} {r(node['right'])})"
    if cls == "BETWEEN":
        return f"({r(node['input'])} BETWEEN {r(node['lower'])} AND {r(node['upper'])})"
    if cls == "CONSTANT":
        literal = node["literal"]
        if literal["kind"] == "STRING":
            return quote(str(literal["text"]))
        if literal["kind"] == "NULL_LITERAL":
            return "NULL"
        return str(literal["text"])
    if cls == "COLUMN_REF":
        return column(str(node["column_names"][0]))
    if cls == "FUNCTION":
        name = str(node["function_name"])
        arguments = [r(a["expression"]) for a in node.get("arguments") or []]
        if name in _INFIX:
            if len(arguments) == 1:
                return f"(- {arguments[0]})"
            return f"({arguments[0]} {_INFIX[name]} {arguments[1]})"
        if name == "list_value":
            return "[" + ", ".join(arguments) + "]"
        return f"{name}({', '.join(arguments)})"
    if cls == "CAST":
        word = "TRY_CAST" if node.get("try_cast") else "CAST"
        return f"{word}({r(node['child'])} AS {node['type_expr']['type_name']})"
    if cls == "SUBQUERY":
        select = node["subquery"]["node"]
        where = select.get("where_clause")
        clause = "" if where is None else f" WHERE {render(where)}"
        inner = f"SELECT descendant FROM {DESCENDANTS}{clause}"
        return f"({r(node['child'])} IN ({inner}))"
    raise FilterError(f"no SQL for {kind}")


def _constant(value: Any) -> Node:
    """A literal node."""
    if value is None:
        return {"class": "CONSTANT", "type": "VALUE_CONSTANT", "literal": {"kind": "NULL_LITERAL", "text": ""}}  # fmt: skip
    if isinstance(value, bool):
        literal = {"kind": "BOOLEAN", "text": "true" if value else "false"}
    elif isinstance(value, str):
        literal = {"kind": "STRING", "text": value}
    else:
        raise TypeError(f"no literal for {value!r}")
    return {"class": "CONSTANT", "type": "VALUE_CONSTANT", "literal": literal}


def expand_descendants(pred: Predicate, run: Callable[[str], list[str]]) -> Predicate:
    """`pred` with each subquery of descendants as the `IN` list of its terms.

    `run` runs the SQL of one subquery and returns its terms. An empty list gives
    `false`, as `IN` of an empty subquery does.
    """

    def walk(node: Any) -> Any:
        if isinstance(node, list):
            return [walk(v) for v in node]
        if not isinstance(node, dict):
            return node
        if node.get("class") == "SUBQUERY":
            select = node["subquery"]["node"]
            where = select.get("where_clause")
            clause = "" if where is None else f" WHERE {render(where)}"
            terms = sorted(set(run(f"SELECT descendant FROM {DESCENDANTS}{clause}")))
            if not terms:
                return _constant(False)
            children = [walk(node["child"]), *(_constant(t) for t in terms)]
            return {"class": "OPERATOR", "type": "COMPARE_IN", "children": children}
        return {k: walk(v) for k, v in node.items()}

    return Predicate(walk(pred.tree))


def from_mapping(obs_filter: Filter | None) -> str:
    """The SQL of the SDK's mapping filter.

    Each column gives one `IN` term, and the terms join with `AND`. A multi-valued
    column gives explicit membership. The SQL splits on the exact separator of the CxG
    schema; `catalog.members` also strips spaces around each member, so the two agree
    for values that follow the schema.
    """
    terms = []
    for column, wanted in _catalog.canonicalize(obs_filter).items():
        name = identifier(column)
        if column == PRIMARY:
            flags = sorted({w == "True" for w in wanted})
            terms.append(f"{name} IN ({', '.join(str(f).lower() for f in flags)})")
            continue
        values = ", ".join(quote(v) for v in sorted(wanted))
        if column in MULTIVALUED:
            split = f"string_split({name}, {quote(MULTIVALUE_SEP)})"
            terms.append(f"list_has_any({split}, [{values}])")
        else:
            terms.append(f"{name} IN ({values})")
    return " AND ".join(terms) if terms else "true"


def from_value_filter(value_filter: str) -> str:
    """The SQL of a TileDB-SOMA `value_filter` string of Census 1.

    The translator takes `and`, `or`, `not`, the six comparisons, `in [...]`, and
    `not in [...]`, with string, number, and boolean literals. It raises `FilterError`
    for any other construct, for example a comparison of two columns, so that a
    migration error can name the equivalent filter.

    Each term keeps its Census 1 meaning, also on a multi-valued column: Census 1
    compares the joined string, and so does SQL.
    """
    try:
        tree = ast.parse(value_filter.strip(), mode="eval").body
    except SyntaxError as error:
        raise FilterError(f"not a value filter: {value_filter!r}") from error
    return _sql_of(tree)


_SQL_COMPARISONS: dict[type, str] = {
    ast.Eq: "=",
    ast.NotEq: "<>",
    ast.Lt: "<",
    ast.LtE: "<=",
    ast.Gt: ">",
    ast.GtE: ">=",
}


def _sql_literal(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return quote(value)
    if isinstance(value, (int, float)):
        return repr(value)
    raise FilterError(f"no SQL literal for {value!r}")


def _sql_of(node: ast.expr) -> str:
    if isinstance(node, ast.BoolOp):
        joiner = " AND " if isinstance(node.op, ast.And) else " OR "
        return "(" + joiner.join(_sql_of(v) for v in node.values) + ")"
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return f"(NOT {_sql_of(node.operand)})"
    if isinstance(node, ast.Compare):
        terms, left = [], node.left
        for op, right in zip(node.ops, node.comparators, strict=True):
            if not isinstance(left, ast.Name):
                raise FilterError(f"a term starts with a column: {ast.unparse(left)!r}")
            try:
                value = ast.literal_eval(right)
            except ValueError as error:
                raise FilterError(f"not a literal: {ast.unparse(right)!r}") from error
            column = identifier(left.id)
            if isinstance(op, (ast.In, ast.NotIn)):
                if not isinstance(value, (list, tuple, set)):
                    raise FilterError(f"`in` takes a list: {ast.unparse(right)!r}")
                word = "IN" if isinstance(op, ast.In) else "NOT IN"
                listed = ", ".join(_sql_literal(v) for v in value)
                terms.append(f"({column} {word} ({listed}))" if value else "false")
            elif type(op) in _SQL_COMPARISONS:
                terms.append(
                    f"({column} {_SQL_COMPARISONS[type(op)]} {_sql_literal(value)})"
                )
            else:
                raise FilterError(f"no SQL term for {ast.unparse(node)!r}")
            left = right
        return terms[0] if len(terms) == 1 else "(" + " AND ".join(terms) + ")"
    raise FilterError(f"no SQL term for {ast.unparse(node)!r}")


def from_filter(obs_filter: object) -> Predicate:
    """A mapping filter, a SQL filter, or None (every cell), as a `Predicate`."""
    if obs_filter is None:
        return parse("true")
    if isinstance(obs_filter, str):
        return parse(obs_filter)
    if isinstance(obs_filter, Mapping):
        return parse(from_mapping(obs_filter))
    raise TypeError(f"a SQL census takes a mapping or a SQL filter, not {obs_filter!r}")


# The catalog side of each column with a value set: a list of values, or NULL when
# the row records none.
def _value_set(column: str) -> str | None:
    if column == "dataset_id":
        return "[dataset_id]"
    if column == PRIMARY:
        return (
            f"CASE {identifier(PRIMARY)} WHEN 'all_true' THEN [true] "
            "WHEN 'all_false' THEN [false] WHEN 'mixed' THEN [true, false] END"
        )
    if column in INDEXED:
        return identifier(column)
    return None


_LAMBDA = "census2_v"


def catalog_test(pred: Predicate) -> str:
    """A SQL test over the catalog: false when the row proves that no cell matches.

    See "Push-down" in the module docstring.
    """
    return _catalog_test(_nnf(pred.tree, negate=False))


def _nnf(node: Node, *, negate: bool) -> Node:
    """Negation normal form: each NOT applies to one term."""
    cls, kind = node.get("class"), node.get("type")
    if cls == "CONJUNCTION":
        flipped = {"CONJUNCTION_AND": "CONJUNCTION_OR", "CONJUNCTION_OR": "CONJUNCTION_AND"}  # fmt: skip
        children = [_nnf(c, negate=negate) for c in node["children"]]
        kind = flipped[str(kind)] if negate else kind
        return {**node, "type": kind, "children": children}
    if cls == "OPERATOR" and kind == "OPERATOR_NOT":
        return _nnf(node["children"][0], negate=not negate)
    if negate:
        return {"class": "OPERATOR", "type": "OPERATOR_NOT", "children": [node]}
    return node


def _catalog_test(node: Node) -> str:
    if node.get("class") == "CONJUNCTION":
        joiner = " AND " if node["type"] == "CONJUNCTION_AND" else " OR "
        return "(" + joiner.join(_catalog_test(c) for c in node["children"]) + ")"
    named = _columns(node)
    if not named:
        return f"coalesce({render(node)}, false)"
    if len(named) > 1:
        return "true"  # a term on two columns proves nothing here
    (column,) = named
    values = _value_set(column)
    if values is None:
        return "true"  # the catalog records no value set of this column
    body = render(node, lambda name: _LAMBDA if name == column else identifier(name))
    # The value set has no null; a cell can hold one.
    some = (
        f"list_bool_or(list_transform(list_append({values}, NULL), "
        f"lambda {_LAMBDA}: {body}))"
    )
    return f"(CASE WHEN ({values}) IS NULL THEN true ELSE coalesce({some}, false) END)"


@dataclass
class Options:
    """How the SDK runs DuckDB.

    The defaults suit a query through the catalog. The other values of `files` are
    for the comparison in the benchmark:

    - `catalog`: the fragments of the catalog rows that the catalog test keeps.
    - `all`: the fragments of every catalog row of the organism. The scan then relies
      on the Parquet statistics alone.
    - `glob`: `glob`, a pattern over the fragment directory, and a test of each file
      name against the fragments of the catalog. DuckDB lists the directory first.
    - `iceberg`: `iceberg_scan` of the metadata file in `iceberg`, with a
      `dataset_id IN (...)` term for the datasets that the catalog keeps.
    """

    # DuckDB's worker threads. None keeps DuckDB's default, the number of cores.
    threads: int | None = None
    # DuckDB's memory limit, for example "8GB". None keeps DuckDB's default.
    memory_limit: str | None = None
    files: str = "catalog"
    glob: str | None = None
    iceberg: str | None = None
    # Unify the columns of the fragments by name. A fragment without a column, for
    # example one of an older writer, then gives nulls in it.
    union_by_name: bool = True
    # Count the HTTP requests of each query in its statistics (DuckDB's HTTP log).
    count_requests: bool = False
    # S3 credentials: "env" takes the AWS_* variables, as obstore does; "chain" is
    # DuckDB's credential chain through its `aws` extension; "auto" is "env" when
    # AWS_ACCESS_KEY_ID is set and "chain" otherwise.
    credentials: str = "auto"
    # A Parquet file of the descendants table, for ontology closure.
    descendants: str | None = None
    # Rows per Arrow batch of a result.
    batch_rows: int = 1_000_000
    # Other DuckDB settings, set before the configuration locks.
    settings: dict[str, str] = field(default_factory=dict)


# The CA bundles that DuckDB's httpfs looks for. With external access off, the bundle
# must be an allowed path, or the TLS connection cannot start.
_CA_BUNDLES = (
    "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem",
    "/etc/pki/tls/certs/ca-bundle.crt",
    "/etc/ssl/ca-bundle.pem",
    "/etc/ssl/cert.pem",
)


@dataclass
class QueryStats:
    """What one DuckDB query did."""

    catalog_seconds: float = 0.0
    scan_seconds: float = 0.0
    arrow_seconds: float = 0.0
    files: int = 0
    rows: int = 0
    requests: dict[str, int] = field(default_factory=dict)
    mb_requested: float = 0.0
    extra: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        out = {k: v for k, v in self.__dict__.items() if k != "extra"}
        return {**out, **self.extra}


class Engine:
    """One DuckDB connection over the catalog and the fragments of a census.

    The connection holds the catalog as a table. It may read only the directories of
    the fragments, and its configuration is locked; see the module docstring.
    """

    def __init__(self, catalog: pd.DataFrame, options: Options | None = None) -> None:
        self.options = options or Options()
        config: dict[str, Any] = {}
        if self.options.threads is not None:
            config["threads"] = self.options.threads
        if self.options.memory_limit is not None:
            config["memory_limit"] = self.options.memory_limit
        self.catalog = catalog.reset_index(drop=True)
        locations = [fragments.uri(r) for r in _records(self.catalog)]
        self.fragment_uris = locations
        t0 = time.monotonic()
        con = duckdb.connect(config=config)
        for key, value in self.options.settings.items():
            con.execute(f"SET {key} = {quote(str(value))}")
        remote = any(_catalog.is_remote(u) for u in locations)
        if self.options.iceberg is not None:
            remote = remote or _catalog.is_remote(self.options.iceberg)
        if remote:
            con.install_extension("httpfs")
            con.load_extension("httpfs")
            self._secret(con)
        if self.options.files == "iceberg":
            con.install_extension("iceberg")
            con.load_extension("iceberg")
        if self.options.count_requests:
            con.execute("CALL enable_logging('HTTP')")
        con.register("census2_catalog_arrow", catalog_arrow(self.catalog))
        con.execute("CREATE TABLE catalog AS SELECT * FROM census2_catalog_arrow")
        con.unregister("census2_catalog_arrow")
        if self.options.descendants is not None:
            path = self.options.descendants
            con.execute(
                f"CREATE TABLE {DESCENDANTS} AS SELECT * FROM read_parquet({quote(path)})"
            )
        directories = sorted({_parent(u) for u in locations})
        if self.options.iceberg is not None:
            directories.append(_parent(_parent(self.options.iceberg)))
        if self.options.glob is not None:
            directories.append(_parent(self.options.glob))
        paths = [p for p in _CA_BUNDLES if os.path.exists(p)]
        con.execute(f"SET allowed_directories = [{', '.join(map(quote, directories))}]")
        con.execute(f"SET allowed_paths = [{', '.join(map(quote, paths))}]")
        con.execute("SET enable_external_access = false")
        con.execute("SET lock_configuration = true")
        self.con = con
        self.open_seconds = round(time.monotonic() - t0, 3)

    def _secret(self, con: duckdb.DuckDBPyConnection) -> None:
        """An S3 secret: the AWS_* variables, or DuckDB's credential chain."""
        mode = self.options.credentials
        if mode == "auto":
            mode = "env" if os.environ.get("AWS_ACCESS_KEY_ID") else "chain"
        region = os.environ.get("AWS_REGION") or os.environ.get("AWS_DEFAULT_REGION")
        parts = ["TYPE s3"]
        if region:
            parts.append(f"REGION {quote(region)}")
        if mode == "chain":
            con.install_extension("aws")
            con.load_extension("aws")
            parts.append("PROVIDER credential_chain")
        elif mode == "env":
            parts.append(f"KEY_ID {quote(os.environ['AWS_ACCESS_KEY_ID'])}")
            parts.append(f"SECRET {quote(os.environ['AWS_SECRET_ACCESS_KEY'])}")
            token = os.environ.get("AWS_SESSION_TOKEN")
            if token:
                parts.append(f"SESSION_TOKEN {quote(token)}")
        else:
            raise ValueError(f"no credentials mode {mode!r}")
        con.execute(f"CREATE SECRET census2_s3 ({', '.join(parts)})")

    def terms(self, sql: str) -> list[str]:
        """The terms of one subquery of the descendants table."""
        if self.options.descendants is None:
            raise FilterError("a filter on descendants needs Options.descendants")
        return [str(r[0]) for r in self.con.execute(sql).fetchall()]

    def predicate(self, obs_filter: object) -> Predicate:
        pred = from_filter(obs_filter)
        if pred.has_subquery():
            pred = expand_descendants(pred, self.terms)
        return pred

    def candidates(self, organism: str, pred: Predicate) -> tuple[pd.DataFrame, float]:
        """The catalog rows of `organism` that the catalog test keeps, and its time."""
        t0 = time.monotonic()
        sql = (
            "SELECT census2_row FROM catalog WHERE organism = ? AND "
            f"{catalog_test(pred)} ORDER BY census2_row"
        )
        try:
            keep = [int(r[0]) for r in self.con.execute(sql, [organism]).fetchall()]
        except duckdb.Error as error:
            raise FilterError(f"the catalog test failed: {error}") from error
        return self.catalog.iloc[keep], round(time.monotonic() - t0, 4)

    def _source(self, rows: pd.DataFrame) -> tuple[str, str, int]:
        """The FROM clause of a scan, a term that keeps only the rows' datasets, and
        the number of files that the clause names."""
        by_name = "true" if self.options.union_by_name else "false"
        mode = self.options.files
        if mode in ("catalog", "all"):
            scoped = rows if mode == "catalog" else self._organism(rows)
            files = [fragments.uri(r) for r in _records(scoped)]
            listed = ", ".join(quote(f) for f in files)
            return (
                f"read_parquet([{listed}], union_by_name = {by_name})",
                "true",
                len(files),
            )
        if mode == "glob":
            if self.options.glob is None:
                raise ValueError("files='glob' needs Options.glob")
            files = [fragments.uri(r) for r in _records(rows)]
            listed = ", ".join(quote(f) for f in files) or "NULL"
            source = (
                f"read_parquet({quote(self.options.glob)}, union_by_name = {by_name}, "
                "filename = true)"
            )
            return source, f"(filename IN ({listed}))", len(files)
        if mode == "iceberg":
            if self.options.iceberg is None:
                raise ValueError("files='iceberg' needs Options.iceberg")
            ids = ", ".join(quote(str(d)) for d in rows["dataset_id"]) or "NULL"
            source = f"iceberg_scan({quote(self.options.iceberg)})"
            return source, f"(dataset_id IN ({ids}))", len(rows)
        raise ValueError(f"no files mode {mode!r}")

    def _organism(self, rows: pd.DataFrame) -> pd.DataFrame:
        organisms = set(rows["organism"]) if len(rows) else set()
        return self.catalog[self.catalog["organism"].isin(organisms)]

    def _http(self) -> tuple[dict[str, int], float]:
        """The HTTP requests since the last call, by method, and the MB of the GETs."""
        if not self.options.count_requests:
            return {}, 0.0
        found = self.con.execute(
            "SELECT request.type, count(*), "
            "sum(TRY_CAST(response.headers['Content-Length'] AS BIGINT)) "
            "FROM duckdb_logs_parsed('HTTP') GROUP BY ALL"
        ).fetchall()
        self.con.execute("CALL truncate_duckdb_logs()")
        counts = {str(t): int(n) for t, n, _ in found}
        nbytes = sum(int(b or 0) for t, _, b in found if t == "GET")
        return counts, round(nbytes / 2**20, 2)

    def select(
        self,
        rows: pd.DataFrame,
        pred: Predicate,
        columns: Sequence[str] | None,
        *,
        dictionary: bool = True,
    ) -> tuple[pa.Table, QueryStats]:
        """The selected cells of the fragments of `rows`, as one Arrow table.

        `columns` names the `obs` columns; None gives all. The cell key and the row
        columns are always in the table. With `dictionary`, each string column except
        `observation_joinid` is dictionary-encoded, as a pandas categorical needs.
        """
        stats = QueryStats()
        self._http()
        if not len(rows):
            return pa.table({}), stats
        source, scope, stats.files = self._source(rows)
        wanted = (
            "*"
            if columns is None
            else ", ".join(
                identifier(c)
                for c in dict.fromkeys([*KEY, *columns, *fragments.ROW_COLUMNS])
            )
        )
        if columns is None and self.options.files == "glob":
            wanted = "* EXCLUDE (filename)"
        sql = f"SELECT {wanted} FROM {source} WHERE {scope} AND {pred.text}"
        t0 = time.monotonic()
        try:
            reader = self.con.sql(sql).to_arrow_reader(self.options.batch_rows)
            t1 = time.monotonic()
            batches = list(_encoded(reader, dictionary))
        except duckdb.BinderException as error:
            raise FilterError(str(error).splitlines()[0]) from error
        t2 = time.monotonic()
        schema = reader.schema
        if dictionary:
            schema = pa.schema([_dictionary_field(f) for f in schema])
        table = pa.Table.from_batches(batches, schema=schema)
        stats.scan_seconds = round(t1 - t0, 3)
        stats.arrow_seconds = round(t2 - t1, 3)
        stats.rows = int(table.num_rows)
        stats.requests, stats.mb_requested = self._http()
        return table, stats

    def count(
        self, rows: pd.DataFrame, pred: Predicate, by: Sequence[str]
    ) -> tuple[pd.DataFrame, QueryStats]:
        """The number of selected cells in each group of `by`, counted by DuckDB."""
        stats = QueryStats()
        self._http()
        keys = ", ".join(identifier(c) for c in by)
        if not len(rows):
            return pd.DataFrame({**{c: [] for c in by}, "n_cells": []}), stats
        source, scope, stats.files = self._source(rows)
        select = f"{keys}, count(*) AS n_cells" if by else "count(*) AS n_cells"
        group = f" GROUP BY {keys}" if by else ""
        sql = f"SELECT {select} FROM {source} WHERE {scope} AND {pred.text}{group}"
        t0 = time.monotonic()
        try:
            table = self.con.sql(sql).to_arrow_table()
        except duckdb.BinderException as error:
            raise FilterError(str(error).splitlines()[0]) from error
        stats.scan_seconds = round(time.monotonic() - t0, 3)
        stats.rows = int(table.num_rows)
        stats.requests, stats.mb_requested = self._http()
        frame = cast("pd.DataFrame", table.to_pandas())
        return frame[[*by, "n_cells"]], stats


def _dictionary_field(f: Any) -> Any:
    if pa.types.is_string(f.type) and f.name != "observation_joinid":
        return pa.field(f.name, pa.dictionary(pa.int32(), pa.string()))
    return f


def _encoded(reader: Any, dictionary: bool) -> Iterable[Any]:
    """The batches of `reader`, with the string columns dictionary-encoded.

    The encoding of one batch runs while DuckDB produces the next one, in a pool of
    threads; pyarrow does not hold the GIL while it encodes.
    """
    if not dictionary:
        yield from reader
        return

    def encode(batch: Any) -> Any:
        arrays = []
        for f, array in zip(batch.schema, batch.columns, strict=True):
            if pa.types.is_string(f.type) and f.name != "observation_joinid":
                array = pc.dictionary_encode(array).cast(
                    pa.dictionary(pa.int32(), pa.string())
                )
            elif pa.types.is_large_string(f.type) or pa.types.is_string_view(f.type):
                array = array.cast(pa.string())
            arrays.append(array)
        names = list(batch.schema.names)
        return pa.RecordBatch.from_arrays(arrays, names=names)

    with ThreadPoolExecutor(max_workers=4) as pool:
        yield from pool.map(encode, reader)


def catalog_arrow(catalog: pd.DataFrame) -> pa.Table:
    """The catalog columns that the catalog test reads, as an Arrow table.

    Each indexed column is a list of strings, null where the row records no value
    set, also when the catalog has no such column. `census2_row` is the row's position.
    """
    n = len(catalog)
    arrays: dict[str, Any] = {
        "census2_row": pa.array(np.arange(n, dtype=np.int64)),
        "dataset_id": pa.array(catalog["dataset_id"].astype(str).tolist()),
        "organism": pa.array(catalog["organism"].astype(str).tolist()),
    }
    primary = catalog.get(PRIMARY, pd.Series([None] * n))
    arrays[PRIMARY] = pa.array(
        [None if _catalog.unrecorded(v) else str(v) for v in primary], pa.string()
    )
    for column in INDEXED:
        values = catalog.get(column, pd.Series([None] * n))
        arrays[column] = pa.array(
            [None if _catalog.unrecorded(v) else [str(x) for x in v] for v in values],
            pa.list_(pa.string()),
        )
    return pa.table(arrays)


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return cast("list[dict[str, Any]]", frame.to_dict("records"))


def _parent(uri: str) -> str:
    """The directory of a file, with a trailing slash."""
    if _catalog.is_remote(uri):
        return uri.rstrip("/").rsplit("/", 1)[0] + "/"
    return str(pathlib.Path(uri).resolve().parent) + os.sep


def by_dataset(arrow: pa.Table) -> dict[str, pa.Table]:
    """The rows of each dataset, in store row order.

    DuckDB keeps the order of the files and of their rows, so the rows of each
    dataset are together; the SDK checks this and sorts when they are not.
    """
    if not arrow.num_rows:
        return {}
    ids = arrow.column("dataset_id").combine_chunks()
    if not pa.types.is_dictionary(ids.type):
        ids = pc.dictionary_encode(ids)
    codes = ids.indices.to_numpy(zero_copy_only=False)
    rows_of = arrow.column(fragments.STORE_ROW).to_numpy()
    change = np.flatnonzero(codes[1:] != codes[:-1]) + 1
    grouped = len(change) + 1 == len(np.unique(codes)) and not bool(
        np.any((np.diff(rows_of) <= 0) & (codes[1:] == codes[:-1]))
    )
    if not grouped:
        order = np.lexsort((rows_of, codes))
        arrow, codes = arrow.take(order), codes[order]
        change = np.flatnonzero(codes[1:] != codes[:-1]) + 1
    names = ids.dictionary.to_pylist()
    starts = [0, *change.tolist()]
    stops = [*change.tolist(), len(codes)]
    return {
        str(names[codes[a]]): arrow.slice(a, b - a)
        for a, b in zip(starts, stops, strict=True)
    }


def get_anndata(
    census: "Census",
    organism: str,
    obs_filter: object,
    layers: Sequence[str],
    max_workers: int,
) -> "anndata.AnnData | None":
    """`Census.get_anndata` through DuckDB. See the module docstring."""
    import anndata
    from anndata.io import read_elem

    import census2_sdk as sdk
    from census2_sdk.gene_reference import MEASURED_COLUMN, ORGANISMS

    engine = census.sql_engine
    pred = engine.predicate(obs_filter)
    rows, catalog_seconds = engine.candidates(organism, pred)
    sdk._check_composition(rows, ORGANISMS[organism]["version"])
    arrow, stats = engine.select(rows, pred, None)
    stats.catalog_seconds = catalog_seconds
    stats.extra["catalog_candidates"] = len(rows)
    census.last_scan = [stats.as_dict()]
    cells = by_dataset(arrow)
    scope = engine.catalog if engine.options.files != "catalog" else rows
    records = [r for r in _records(scope) if str(r["dataset_id"]) in cells]
    if engine.options.files != "catalog":
        # The scan read files that the catalog test did not choose; the stores of
        # the hits must still compose.
        sdk._check_composition(
            scope[scope["dataset_id"].astype(str).isin(cells)],
            ORGANISMS[organism]["version"],
        )

    def part(row: dict[str, Any]) -> "anndata.AnnData":
        return sdk._read_cells(row, cells[str(row["dataset_id"])], layers)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        parts = list(pool.map(part, records))
    if not sum(p.n_obs for p in parts):
        return None
    uris = [str(r["store_uri"]) for r in records]
    # The benchmark counts the reads of `catalog.open_store`; look it up at call time.
    var = cast("pd.DataFrame", read_elem(_catalog.open_store(uris[0])["var"]))
    for p, uri in zip(parts, uris, strict=True):
        if p.n_vars != len(var):
            raise sdk.SourceConflict(
                f"{uri}: X has {p.n_vars} genes, and {uris[0]} has {len(var)}; "
                "the stores are not on one gene reference"
            )
        p.var_names = cast("Sequence[str]", var.index)
    result = anndata.concat(parts, join="inner")
    if not result.obs_names.is_unique:
        raise sdk.SourceConflict("the cell key is not unique across the result")
    result.var = var.drop(columns=MEASURED_COLUMN, errors="ignore")
    result.uns["census2"] = {
        "sources": census.sources,
        "stores_in_catalog": int((census.catalog["organism"] == organism).sum()),
        "catalog_candidates": len(rows),
        "stores_opened": len(parts),
        "obs_sql": [stats.as_dict()],
    }
    return result


def get_obs(
    census: "Census",
    organism: str,
    obs_filter: object,
    columns: Sequence[str] | None,
) -> pd.DataFrame:
    """`Census.get_obs` through DuckDB. The SDK opens no store."""
    import census2_sdk as sdk

    engine = census.sql_engine
    pred = engine.predicate(obs_filter)
    rows, catalog_seconds = engine.candidates(organism, pred)
    arrow, stats = engine.select(rows, pred, columns)
    stats.catalog_seconds = catalog_seconds
    t0 = time.monotonic()
    frame = pd.DataFrame()
    if arrow.num_rows:
        arrow = arrow.unify_dictionaries()
        frame = sdk._obs_frame(arrow, keep_rows=True)
        key = pc.binary_join_element_wise(
            pc.cast(arrow.column("dataset_id"), pa.string()),
            arrow.column("observation_joinid"),
            "_",
        )
        frame.index = pd.Index(key.to_pandas())
    stats.extra["pandas_seconds"] = round(time.monotonic() - t0, 3)
    frame.attrs["census2"] = {
        "catalog_candidates": len(rows),
        "obs_sql": [stats.as_dict()],
    }
    return frame


def count_obs(
    census: "Census",
    organism: str,
    obs_filter: object,
    by: Sequence[str],
) -> pd.DataFrame:
    """`Census.count_obs` through DuckDB. The SDK opens no store."""
    engine = census.sql_engine
    pred = engine.predicate(obs_filter)
    rows, catalog_seconds = engine.candidates(organism, pred)
    frame, stats = engine.count(rows, pred, by)
    stats.catalog_seconds = catalog_seconds
    frame.attrs["census2"] = {
        "catalog_candidates": len(rows),
        "obs_sql": [stats.as_dict()],
    }
    return frame
