"""Tabulate the query corpus in the Iceberg expression language of the obs table.

For each of the bound queries of `corpus.expand()`, this script translates the Census 1
`value_filter` with `census2_sdk.expressions.from_value_filter`, binds the expression to
the schema of the obs fragments, and records one status:

- `expressible`: the expression binds to the columns that a store holds now.
- `expressible, membership`: as above, and a term is on a multi-valued column. The SDK
  rewrites it (`expressions.with_members`); Iceberg has no membership term.
- `expressible, later column`: the expression needs a column that the store does not
  hold yet, such as `tissue_general` or `nnz`.
- `var axis`: the query filters genes only, which is outside the obs table.
- `no obs filter`: the accessor reads no cells (`get_var`, `get_presence_matrix`).

It also records whether the mapping filter of the SDK expresses the query, through the
harness adapter's `translate`.

    python -m harness.filter_language --out results/filter_language.json
"""

import argparse
import collections
import json
import pathlib
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
from pyiceberg.catalog import Catalog
from pyiceberg.expressions.visitors import bind
from pyiceberg.schema import assign_fresh_schema_ids

from census2_sdk import expressions
from census2_sdk.catalog import MULTIVALUED
from harness.corpus import BoundQuery, expand, load_values
from harness.sdk_adapter import Unsupported, translate

# The obs columns of a store now, as in census2_dataset_writer.store_format.OBS_COLUMNS.
STORE_COLUMNS: dict[str, Any] = {
    name: pa.bool_() if name == "is_primary_data" else pa.string()
    for name in [
        "dataset_id",
        "assay",
        "assay_ontology_term_id",
        "cell_type",
        "cell_type_ontology_term_id",
        "development_stage",
        "development_stage_ontology_term_id",
        "disease",
        "disease_ontology_term_id",
        "donor_id",
        "is_primary_data",
        "observation_joinid",
        "self_reported_ethnicity",
        "self_reported_ethnicity_ontology_term_id",
        "sex",
        "sex_ontology_term_id",
        "suspension_type",
        "tissue",
        "tissue_ontology_term_id",
        "tissue_type",
    ]
}
# The columns that Census 1 also holds and that a store does not hold yet.
LATER_COLUMNS: dict[str, Any] = {
    "tissue_general": pa.string(),
    "tissue_general_ontology_term_id": pa.string(),
    "raw_sum": pa.float64(),
    "nnz": pa.int64(),
    "raw_mean_nnz": pa.float64(),
    "raw_variance_nnz": pa.float64(),
    "n_measured_vars": pa.int64(),
}


def schema() -> Any:
    """The Iceberg schema of a fragment with the later columns."""
    fields = [pa.field(n, t) for n, t in {**STORE_COLUMNS, **LATER_COLUMNS}.items()]
    return assign_fresh_schema_ids(Catalog._convert_schema_if_needed(pa.schema(fields)))


def classify(query: BoundQuery, pins: dict[str, Any], table: Any) -> dict[str, Any]:
    """The status of one bound query in the obs table language and in the SDK mapping."""
    out: dict[str, Any] = {"label": query.label, "key": query.key}
    if query.accessor in ("get_var", "get_presence_matrix"):
        out["status"] = "no obs filter"
    elif not query.obs_filter:
        out["status"] = "var axis"
    else:
        try:
            expr = expressions.from_value_filter(query.obs_filter)
        except ValueError as error:
            out.update(status="not expressible", reason=str(error))
            return out
        out["iceberg"] = expressions.to_text(expr)
        bind(table, expr, case_sensitive=True)  # raises on a column or type error
        named = expressions.columns(expr)
        if named & set(LATER_COLUMNS):
            out["status"] = "expressible, later column"
            out["later_columns"] = sorted(named & set(LATER_COLUMNS))
        elif named & MULTIVALUED:
            out["status"] = "expressible, membership"
        else:
            out["status"] = "expressible"
    if query.obs_filter:
        try:
            out["sdk_mapping"] = translate(query.obs_filter, pins)
        except Unsupported as error:
            out["sdk_mapping"] = None
            out["sdk_reason"] = str(error)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    values = load_values()
    table = schema()
    rows = [
        classify(q, values.get(f"{q.modality}/{q.experiment}", {}), table)
        for q in expand(values)
    ]
    by_status = collections.Counter(r["status"] for r in rows)
    mapping = sum(1 for r in rows if r.get("sdk_mapping") is not None)
    with_filter = sum(1 for r in rows if "sdk_mapping" in r)
    print(f"{len(rows)} bound queries: {dict(by_status)}")
    print(f"SDK mapping filter: {mapping} of {with_filter} queries with an obs filter")
    print("| key | queries | status | Iceberg expression (first) | SDK mapping |")
    print("|---|---|---|---|---|")
    for key in dict.fromkeys(r["key"] for r in rows):
        mine = [r for r in rows if r["key"] == key]
        first = mine[0]
        expressed = sum(1 for r in mine if r.get("sdk_mapping") is not None)
        text = first.get("iceberg", "")
        text = text if len(text) < 90 else text[:87] + "..."
        sdk = f"{expressed}/{len(mine)}" if "sdk_mapping" in first else "n/a"
        print(f"| {key} | {len(mine)} | {first['status']} | `{text}` | {sdk} |")
    if args.out is not None:
        args.out.write_text(json.dumps(rows, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
