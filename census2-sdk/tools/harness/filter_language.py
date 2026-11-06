"""Tabulate the query corpus in the Iceberg expression language of the obs table.

For each of the bound queries of `corpus.expand()`, this script translates the Census 1
`value_filter` with `census2_sdk.expressions.from_value_filter`, binds the expression to
the schema of the obs fragments, and records one status:

- `expressible`: the expression binds to the columns of a store's `obs`.
- `expressible, membership`: as above, and a term is on a multi-valued column. The SDK
  rewrites it (`expressions.with_members`); Iceberg has no membership term.
- `var axis`: the query filters genes only, which is outside the obs table.
- `no obs filter`: the accessor reads no cells (`get_var`, `get_presence_matrix`).

The schema is the fragment schema: the columns of
`census2_dataset_writer.OBS_COLUMNS`, with the Iceberg type of each column type. The
script also records the columns of the expression that the catalog does not index
(`table_only_columns`). Only the table tests a term on such a column, for example a
range on `nnz`; the catalog keeps every dataset for it. And it records whether the
mapping filter of the SDK expresses the query, through the harness adapter's
`translate`.

    python -m harness.filter_language --out results/filter_language.json
"""

import argparse
import collections
import json
import pathlib
from typing import Any

import census2_dataset_writer
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
from pyiceberg.catalog import Catalog
from pyiceberg.expressions.visitors import bind
from pyiceberg.schema import assign_fresh_schema_ids

from census2_sdk import expressions
from census2_sdk.catalog import INDEXED, MULTIVALUED, PRIMARY
from harness.corpus import BoundQuery, expand, load_values
from harness.sdk_adapter import Unsupported, translate

# The Arrow type of each store column type, as the fragment holds it.
ARROW_TYPES: dict[census2_dataset_writer.ColType, Any] = {
    census2_dataset_writer.ColType.CATEGORY: pa.string(),
    census2_dataset_writer.ColType.STRING: pa.string(),
    census2_dataset_writer.ColType.BOOL: pa.bool_(),
    census2_dataset_writer.ColType.INT64: pa.int64(),
    census2_dataset_writer.ColType.FLOAT64: pa.float64(),
}
STORE_COLUMNS: dict[str, Any] = {
    name: ARROW_TYPES[col_type] for name, col_type in census2_dataset_writer.OBS_COLUMNS
}


def schema() -> Any:
    """The Iceberg schema of the `obs` columns of a fragment."""
    fields = [pa.field(n, t) for n, t in STORE_COLUMNS.items()]
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
        out["status"] = (
            "expressible, membership" if named & MULTIVALUED else "expressible"
        )
        catalog_columns = {*INDEXED, "dataset_id", PRIMARY}
        out["table_only_columns"] = sorted(named - catalog_columns)
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
    table_only = collections.Counter(
        c for r in rows for c in r.get("table_only_columns", [])
    )
    print(f"terms that only the table tests, by column: {dict(table_only)}")
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
