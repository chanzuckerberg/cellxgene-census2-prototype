"""Tabulate the query builder over the harness corpus and the probes of the SQL spike.

Two parts, in the form of `spike/obs-duckdb`'s `filter_language_sql.py`, so that the
two tables compare row for row:

1. The corpus. For each bound query of `corpus.expand()` with an `obs` filter, the
   script writes the Census 1 `value_filter` in builder form (`builder_of`), by a
   translation of its Python syntax tree. This is a tool of the harness: the SDK
   parses no text form of a query.
2. The probes: one query for each kind of term, with its builder form, or the reason
   that the builder has none.

With `--catalog`, the script also tests each query against the catalog: the
candidates, the candidates that the catalog proves whole, and the `obs` columns that
the SDK reads from the other candidates. With `--count`, it also selects the cells
over the stores with `get_obs`, which reads no `X`.

    python -m harness.filter_language_builder --catalog catalog-103.parquet \\
        --out results/obs_arrow_native/filter_language.json
"""

import argparse
import ast
import collections
import json
import pathlib
import time
from collections.abc import Callable
from typing import Any

import census2_sdk
from census2_sdk import obs_columns
from census2_sdk.catalog import MULTIVALUE_SEP, MULTIVALUED
from census2_sdk.query import Query, QueryLike, as_query
from harness.corpus import BoundQuery, expand, load_values


class Untranslated(ValueError):
    """A value filter with no builder form."""


def builder_of(text: str, obs: obs_columns.Obs) -> tuple[Query, list[str]]:
    """The builder form of a Census 1 `value_filter`, and notes on its meaning.

    `==` and `in` on a multi-valued column name a joined value in Census 1. The
    builder has no equality on such a column, so each joined value becomes `has_all`
    of its members, which also selects a cell that holds more members. The note says so.
    """
    notes: list[str] = []
    tree = ast.parse(text, mode="eval").body

    def column(node: ast.expr) -> Any:
        if not isinstance(node, ast.Name):
            raise Untranslated(f"not a column: {ast.unparse(node)}")
        found = getattr(obs, node.id, None)
        return found if found is not None else obs[node.id]

    def literal(node: ast.expr) -> Any:
        return ast.literal_eval(node)

    def multi(col: Any, values: list[str], negate: bool) -> Query:
        notes.append(
            f"{col.name}: a joined value of Census 1 becomes has_all of its members"
        )
        terms: list[Query] = [
            col.has_all([m.strip() for m in v.split(MULTIVALUE_SEP)]) for v in values
        ]
        out = terms[0]
        for t in terms[1:]:
            out = out | t
        return ~out if negate else out

    def walk(node: ast.expr) -> QueryLike:
        if isinstance(node, ast.BoolOp):
            parts = [as_query(walk(v)) for v in node.values]
            out = parts[0]
            for p in parts[1:]:
                out = out & p if isinstance(node.op, ast.And) else out | p
            return out
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            return ~as_query(walk(node.operand))
        if isinstance(node, ast.Name):
            return as_query(column(node))
        if not isinstance(node, ast.Compare) or len(node.ops) != 1:
            raise Untranslated(f"no builder form for {ast.unparse(node)}")
        col, op, right = column(node.left), node.ops[0], literal(node.comparators[0])
        name = col.name
        if name in MULTIVALUED and isinstance(
            op, ast.Eq | ast.NotEq | ast.In | ast.NotIn
        ):
            values = [right] if isinstance(op, ast.Eq | ast.NotEq) else list(right)
            return multi(col, values, isinstance(op, ast.NotEq | ast.NotIn))
        result: Any
        if isinstance(op, ast.Eq):
            result = col == right
        elif isinstance(op, ast.NotEq):
            result = col != right
        elif isinstance(op, ast.In):
            result = col.isin(list(right))
        elif isinstance(op, ast.NotIn):
            result = ~col.isin(list(right))
        elif isinstance(op, ast.Lt):
            result = col < right
        elif isinstance(op, ast.LtE):
            result = col <= right
        elif isinstance(op, ast.Gt):
            result = col > right
        elif isinstance(op, ast.GtE):
            result = col >= right
        else:
            raise Untranslated(f"no builder form for {ast.unparse(node)}")
        return as_query(result)

    return as_query(walk(tree)), notes


Probe = tuple[str, Callable[[obs_columns.Obs], QueryLike] | None, str]

# The probes of `filter_language_sql.py`, in builder form. None: no builder form.
PROBES: list[Probe] = [
    (
        "compound boolean",
        lambda o: (
            ((o.tissue == "lung") | (o.tissue_general == "lung"))
            & ~(o.cell_type == "macrophage")
        ),
        "",
    ),
    ("IN", lambda o: o.cell_type.isin(["neutrophil", "monocyte"]), ""),
    (
        "NOT IN",
        lambda o: (
            (o.tissue_general == "lung")
            & ~o.cell_type.isin(["alveolar macrophage", "type II pneumocyte"])
        ),
        "",
    ),
    ("nulls", lambda o: o.raw_variance_nnz.isna(), ""),
    ("membership, multi-valued disease", lambda o: o.disease.has_any(["COVID-19"]), ""),
    (
        "joined value, as Census 1",
        lambda o: o.disease.has_all(["COVID-19", "normal"]),
        "has_all also selects a cell with more members; no exact form of the joined value",
    ),
    (
        "numeric range",
        lambda o: (o.cell_type == "neutrophil") & (o.nnz >= 2000) & (o.raw_sum < 20000),
        "",
    ),
    ("numeric range only", lambda o: o.nnz >= 8000, ""),
    ("prefix", lambda o: o.cell_type.str.startswith("CD8"), ""),
    (
        "substring, case-insensitive",
        lambda o: o.cell_type.str.contains("macrophage", case=False, regex=False),
        "",
    ),
    ("regular expression", lambda o: o.cell_type.str.match(r"CD(4|8)-positive"), ""),
    (
        "ontology closure (T cell)",
        None,
        "out of scope by decision; a caller passes isin of the descendant terms",
    ),
    (
        "arithmetic",
        None,
        (
            "Organism.expression takes a pyarrow expression, such as "
            "raw_sum / nnz > 3, with the list of its columns"
        ),
    ),
]


def classify(query: BoundQuery) -> dict[str, Any]:
    out: dict[str, Any] = {"label": query.label, "key": query.key}
    if query.accessor in ("get_var", "get_presence_matrix"):
        out["status"] = "no obs filter"
        return out
    if not query.obs_filter:
        out["status"] = "var axis"
        return out
    out["value_filter"] = query.obs_filter
    obs = obs_columns.Obs(query.experiment)
    try:
        q, notes = builder_of(query.obs_filter, obs)
    except (Untranslated, TypeError, ValueError) as error:
        out.update(status="refused", reason=str(error))
        return out
    out.update(status="expressible", builder=q.text, notes=notes)
    out["organism"] = query.experiment
    return out


def test_catalog(census: Any, q: Query, count: bool) -> dict[str, Any]:
    """The candidates, the whole ones, the column reads, and the cells of `q`."""
    t0 = time.monotonic()
    rows = census.candidates(q)
    seconds = time.monotonic() - t0
    whole = int(rows["whole"].sum())
    out: dict[str, Any] = {
        "catalog_candidates": len(rows),
        "catalog_whole": whole,
        "catalog_seconds": round(seconds, 4),
        "columns": sorted(q.columns),
        "column_reads": (len(rows) - whole) * len(q.columns),
        "hints": sorted(q.hints),
        "exact": q.exact,
    }
    if count:
        frame = census.get_obs(q, columns=["dataset_id"], max_workers=8)
        ids = frame["dataset_id"].astype(str) if len(frame) else []
        out.update(cells=len(frame), datasets_with_cells=len(set(ids)))
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=pathlib.Path, default=None)
    p.add_argument("--catalog", default=None, help="test each query on the catalog")
    p.add_argument("--count", action="store_true", help="select the cells as well")
    args = p.parse_args()
    values = load_values()
    queries = expand(values)
    corpus = [classify(q) for q in queries]
    census = census2_sdk.open_census([args.catalog]) if args.catalog else None
    if census is not None:
        for bound, row in zip(queries, corpus, strict=True):
            if row["status"] == "expressible" and bound.experiment == "homo_sapiens":
                obs = obs_columns.Obs(bound.experiment)
                q, _ = builder_of(str(bound.obs_filter), obs)
                row.update(test_catalog(census, q, args.count))
    probes = []
    hs = obs_columns.Obs("homo_sapiens")
    for kind, build, note in PROBES:
        probe: dict[str, Any] = {"kind": kind, "note": note}
        if build is None:
            probe["status"] = "no builder form"
        else:
            q = as_query(build(hs))
            probe.update(status="expressible", builder=q.text)
            if census is not None:
                probe.update(test_catalog(census, q, args.count))
        probes.append(probe)

    by_status = collections.Counter(r["status"] for r in corpus)
    with_filter = [r for r in corpus if "value_filter" in r]
    ok = sum(1 for r in with_filter if r["status"] == "expressible")
    print(f"{len(corpus)} bound queries: {dict(by_status)}")
    print(f"of {len(with_filter)} with an obs filter: builder {ok}")
    for r in with_filter:
        tail = f"{r.get('catalog_candidates', '')} {r.get('cells', '')}".strip()
        print(f"| {r['label']} | {r.get('builder', r.get('reason'))} | {tail} |")
    for r in probes:
        tail = f"{r.get('catalog_candidates', '')} {r.get('cells', '')}".strip()
        print(f"| {r['kind']} | {r.get('builder', r['note'])} | {tail} |")
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = {"corpus": corpus, "probes": probes}
        args.out.write_text(json.dumps(payload, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
