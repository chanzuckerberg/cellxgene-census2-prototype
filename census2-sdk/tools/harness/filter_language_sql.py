"""Tabulate the filter languages: the SDK mapping, Iceberg expressions, and SQL.

Two parts:

1. The corpus. For each bound query of `corpus.expand()`, the script translates the
   Census 1 `value_filter` into SQL (`obs_sql.from_value_filter`) and into an Iceberg
   expression (`expressions.from_value_filter`), and binds each one to the fragment
   schema of the writer (`census2_dataset_writer.store_format.OBS_COLUMNS` and the row
   columns). The harness adapter's `translate` gives the mapping filter.
2. The probes: one filter for each kind of term that the corpus does not exercise
   (compound boolean logic, `IN` and `NOT IN`, nulls, membership in a multi-valued
   column, numeric ranges, patterns, ontology closure, arithmetic). Each probe gives
   its form in each language, or none where the language has none.

With `--catalog`, the script also runs each probe over the fragments of a catalog: it
records the datasets that each catalog test keeps and the datasets and cells that the
SQL filter selects. This measures how exact each push-down is on real data.

    python -m harness.filter_language_sql --out results/obs_duckdb/filter_language.json \\
        --catalog catalog-103.parquet --descendants descendants.parquet
"""

import argparse
import collections
import json
import pathlib
from typing import Any, cast

import duckdb
from census2_dataset_writer.store_format import OBS_COLUMNS, ColType

from census2_sdk import catalog as _catalog
from census2_sdk import fragments, obs_sql
from harness.corpus import BoundQuery, expand, load_values
from harness.sdk_adapter import Unsupported, translate

SQL_TYPES = {
    ColType.CATEGORY: "VARCHAR",
    ColType.STRING: "VARCHAR",
    ColType.BOOL: "BOOLEAN",
    ColType.INT64: "BIGINT",
    ColType.FLOAT64: "DOUBLE",
}

# One probe for each kind of term. None: the language has no form for it.
PROBES: list[dict[str, Any]] = [
    {
        "kind": "compound boolean",
        "mapping": None,
        "iceberg": "(tissue = 'lung' OR tissue_general = 'lung') AND NOT cell_type = 'macrophage'",
        "sql": "(tissue = 'lung' OR tissue_general = 'lung') AND NOT cell_type = 'macrophage'",
    },
    {
        "kind": "IN",
        "mapping": {"cell_type": ["neutrophil", "monocyte"]},
        "iceberg": "cell_type IN ('neutrophil', 'monocyte')",
        "sql": "cell_type IN ('neutrophil', 'monocyte')",
    },
    {
        "kind": "NOT IN",
        "mapping": None,
        "iceberg": "tissue_general = 'lung' AND cell_type NOT IN ('alveolar macrophage', 'type II pneumocyte')",
        "sql": "tissue_general = 'lung' AND cell_type NOT IN ('alveolar macrophage', 'type II pneumocyte')",
    },
    {
        "kind": "nulls",
        "mapping": None,
        "iceberg": "raw_variance_nnz IS NULL",
        "sql": "raw_variance_nnz IS NULL",
    },
    {
        "kind": "membership, multi-valued disease",
        # The mapping and Iceberg give `=` a membership meaning on this column.
        "mapping": {"disease": "COVID-19"},
        "iceberg": "disease = 'COVID-19'",
        "sql": "list_contains(string_split(disease, ' || '), 'COVID-19')",
    },
    {
        "kind": "joined value, as Census 1",
        "mapping": None,
        "iceberg": "disease = 'COVID-19 || normal'",
        "sql": "disease = 'COVID-19 || normal'",
    },
    {
        "kind": "numeric range",
        "mapping": None,
        "iceberg": "cell_type = 'neutrophil' AND nnz >= 2000 AND raw_sum < 20000",
        "sql": "cell_type = 'neutrophil' AND nnz >= 2000 AND raw_sum < 20000",
    },
    {
        "kind": "numeric range only",
        "mapping": None,
        "iceberg": "nnz >= 8000",
        "sql": "nnz >= 8000",
    },
    {
        "kind": "prefix",
        "mapping": None,
        "iceberg": "cell_type LIKE 'CD8%'",
        "sql": "cell_type LIKE 'CD8%'",
    },
    {
        "kind": "substring, case-insensitive",
        "mapping": None,
        "iceberg": None,
        "sql": "cell_type ILIKE '%macrophage%'",
    },
    {
        "kind": "regular expression",
        "mapping": None,
        "iceberg": None,
        "sql": "regexp_matches(cell_type, '^CD(4|8)-positive')",
    },
    {
        "kind": "ontology closure (T cell)",
        "mapping": None,
        "iceberg": None,
        "sql": (
            "cell_type_ontology_term_id IN (SELECT descendant FROM descendants "
            "WHERE ancestor = 'CL:0000084')"
        ),
    },
    {
        "kind": "arithmetic",
        "mapping": None,
        "iceberg": None,
        "sql": "raw_sum / nnz > 3.0 AND tissue_general = 'blood'",
    },
]  # fmt: skip


def schema_connection(descendants: str | None = None) -> duckdb.DuckDBPyConnection:
    """An empty `obs` table with the fragment schema, to bind a filter against."""
    columns = [(n, SQL_TYPES[t]) for n, t in OBS_COLUMNS]
    columns += [(c, "BIGINT") for c in fragments.ROW_COLUMNS]
    con = duckdb.connect()
    body = ", ".join(f"{obs_sql.identifier(n)} {t}" for n, t in columns)
    con.execute(f"CREATE TABLE obs ({body})")
    con.execute(
        f"CREATE TABLE {obs_sql.DESCENDANTS} (ontology VARCHAR, ancestor VARCHAR, "
        "descendant VARCHAR)"
    )
    if descendants is not None:
        con.execute(
            f"INSERT INTO {obs_sql.DESCENDANTS} SELECT ontology, ancestor, descendant "
            f"FROM read_parquet({obs_sql.quote(descendants)})"
        )
    return con


def sql_status(text: str, con: duckdb.DuckDBPyConnection) -> tuple[str, str]:
    """`accepted` when the filter parses, passes the check, and binds; else the reason."""
    try:
        pred = obs_sql.parse(text)
        con.execute(f"SELECT count(*) FROM obs WHERE {pred.text}").fetchall()
    except (obs_sql.FilterError, duckdb.Error) as error:
        return "refused", str(error).splitlines()[0]
    return "accepted", pred.text


def iceberg_status(text: str) -> tuple[str, str]:
    """`accepted` when pyiceberg parses the expression and binds it to the schema."""
    import pyarrow as pa  # type: ignore[import-untyped]
    from pyiceberg.catalog import Catalog
    from pyiceberg.expressions.visitors import bind
    from pyiceberg.schema import assign_fresh_schema_ids
    from pyparsing.exceptions import ParseBaseException

    from census2_sdk import expressions

    types = {"VARCHAR": pa.string(), "BOOLEAN": pa.bool_(), "BIGINT": pa.int64()}
    types["DOUBLE"] = pa.float64()
    fields = [pa.field(n, types[SQL_TYPES[t]]) for n, t in OBS_COLUMNS]
    schema = assign_fresh_schema_ids(
        Catalog._convert_schema_if_needed(pa.schema(fields))
    )
    try:
        expr = expressions.from_filter(text)
        bind(schema, expr, case_sensitive=True)
    except (ValueError, TypeError, ParseBaseException) as error:
        return "refused", str(error).splitlines()[0]
    return "accepted", expressions.to_text(expr)


def classify(query: BoundQuery, pins: dict[str, Any], con: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"label": query.label, "key": query.key}
    if query.accessor in ("get_var", "get_presence_matrix"):
        out["status"] = "no obs filter"
        return out
    if not query.obs_filter:
        out["status"] = "var axis"
        return out
    out["value_filter"] = query.obs_filter
    try:
        out["sql"] = obs_sql.from_value_filter(query.obs_filter)
        out["sql_status"], detail = sql_status(out["sql"], con)
        if out["sql_status"] != "accepted":
            out["sql_reason"] = detail
    except obs_sql.FilterError as error:
        out.update(sql_status="refused", sql_reason=str(error))
    try:
        from census2_sdk import expressions

        text = expressions.to_text(expressions.from_value_filter(query.obs_filter))
        out["iceberg_status"], out["iceberg"] = iceberg_status(text)
    except ValueError as error:
        out.update(iceberg_status="refused", iceberg=str(error))
    try:
        out["mapping"] = translate(query.obs_filter, pins)
    except Unsupported as error:
        out.update(mapping=None, mapping_reason=str(error))
    out["status"] = "expressible" if out["sql_status"] == "accepted" else "refused"
    return out


def probe(item: dict[str, Any], con: Any) -> dict[str, Any]:
    out = dict(item)
    out["sql_status"], out["sql_rendered"] = sql_status(item["sql"], con)
    if item["iceberg"] is not None:
        out["iceberg_status"], out["iceberg_rendered"] = iceberg_status(item["iceberg"])
    if item["mapping"] is not None:
        _catalog.canonicalize(item["mapping"])  # raises on a column it refuses
    return out


def run_probes(
    rows: list[dict[str, Any]], catalog: str, descendants: str | None
) -> None:
    """Add the catalog candidates and the selected datasets of each probe, in place."""
    import census2_sdk
    from census2_sdk import expressions

    census = census2_sdk.open_census(
        [catalog], obs_engine=obs_sql.Options(descendants=descendants)
    )
    engine = census.sql_engine
    organism = "homo_sapiens"
    scoped = census.catalog[census.catalog["organism"] == organism]
    for row in rows:
        pred = engine.predicate(row["sql"])
        candidates, seconds = engine.candidates(organism, pred)
        counts, _ = engine.count(candidates, pred, ["dataset_id"])
        row["sql_catalog_candidates"] = len(candidates)
        row["sql_catalog_seconds"] = seconds
        row["datasets_with_cells"] = int((counts["n_cells"] > 0).sum())
        row["cells"] = int(counts["n_cells"].sum())
        if row["iceberg"] is not None:
            expr = expressions.from_filter(row["iceberg"])
            records = cast("list[dict[str, Any]]", scoped.to_dict("records"))
            kept = sum(1 for r in records if expressions.can_match(expr, r))
            row["iceberg_catalog_candidates"] = kept
        if row["mapping"] is not None:
            kept_rows = _catalog.exclude(census.catalog, organism, row["mapping"])
            row["mapping_catalog_candidates"] = len(kept_rows)
    row_count = len(scoped)
    for row in rows:
        row["catalog_rows"] = row_count


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=pathlib.Path, default=None)
    p.add_argument("--catalog", default=None, help="run the probes over its fragments")
    p.add_argument("--descendants", default=None, help="the descendants Parquet")
    args = p.parse_args()
    con = schema_connection(args.descendants)
    values = load_values()
    corpus = [
        classify(q, values.get(f"{q.modality}/{q.experiment}", {}), con)
        for q in expand(values)
    ]
    probes = [probe(item, con) for item in PROBES]
    if args.catalog is not None:
        run_probes(probes, args.catalog, args.descendants)

    by_status = collections.Counter(r["status"] for r in corpus)
    with_filter = [r for r in corpus if "value_filter" in r]
    sql_ok = sum(1 for r in with_filter if r["sql_status"] == "accepted")
    ice_ok = sum(1 for r in with_filter if r.get("iceberg_status") == "accepted")
    mapping_ok = sum(1 for r in with_filter if r.get("mapping") is not None)
    print(f"{len(corpus)} bound queries: {dict(by_status)}")
    print(
        f"of {len(with_filter)} with an obs filter: SQL {sql_ok}, Iceberg {ice_ok}, "
        f"mapping {mapping_ok}"
    )
    print("| key | queries | SQL | Iceberg | mapping | SQL (first) |")
    print("|---|---|---|---|---|---|")
    for key in dict.fromkeys(r["key"] for r in corpus):
        mine = [r for r in corpus if r["key"] == key]
        if "value_filter" not in mine[0]:
            print(f"| {key} | {len(mine)} | {mine[0]['status']} | | | |")
            continue
        n_sql = sum(1 for r in mine if r["sql_status"] == "accepted")
        n_ice = sum(1 for r in mine if r.get("iceberg_status") == "accepted")
        n_map = sum(1 for r in mine if r.get("mapping") is not None)
        text = mine[0].get("sql", "")
        text = text if len(text) < 80 else text[:77] + "..."
        print(
            f"| {key} | {len(mine)} | {n_sql}/{len(mine)} | {n_ice}/{len(mine)} | "
            f"{n_map}/{len(mine)} | `{text}` |"
        )
    print()
    print(
        "| kind | mapping | Iceberg | SQL | candidates (map/ice/sql) | datasets | cells |"
    )
    print("|---|---|---|---|---|---|---|")
    for r in probes:
        cand = "/".join(
            str(r.get(f"{k}_catalog_candidates", "-"))
            for k in ("mapping", "iceberg", "sql")
        )
        print(
            f"| {r['kind']} | {'yes' if r['mapping'] else '-'} | "
            f"{r.get('iceberg_status', '-')} | {r['sql_status']} | {cand} | "
            f"{r.get('datasets_with_cells', '')} | {r.get('cells', '')} |"
        )
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        payload = {"corpus": corpus, "probes": probes}
        args.out.write_text(json.dumps(payload, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
