"""Estimate a catalog test by pyarrow itself: evaluate a query on the value sets.

The SDK of this branch does not do this; the script estimates it for the open
questions. For each catalog row and each query whose columns all have a value set
(the indexed columns, `dataset_id`, and `is_primary_data` from the catalog's
`all_true`, `all_false`, or `mixed`), the script makes a table of every combination of
the row's values and evaluates the query's pyarrow expression on it with
`obs_read.evaluate`. Each cell of the store holds one of the combinations, because
these columns hold no null. Thus no true combination proves that no cell matches,
and all true combinations prove that every cell matches. A query with a column
without a value set, such as `nnz`, gets no estimate.

    python valuesets.py catalog-103.parquet --tools census2-sdk/tools --out valuesets.json
"""

import argparse
import itertools
import json
import pathlib
import sys
import time
from typing import Any

import numpy as np
import pandas as pd

PRIMARY = {"all_true": [True], "all_false": [False], "mixed": [True, False]}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--tools", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    sys.path.insert(0, str(args.tools))
    import census2_sdk
    from bench.filters import QUERIES
    from census2_sdk import obs_read
    from census2_sdk.catalog import INDEXED, unrecorded
    from census2_sdk.query import as_query
    from harness import filter_language_builder as flb
    from harness.corpus import expand, load_values

    census = census2_sdk.open_census([args.catalog])
    hs = census.homo_sapiens
    catalog = census.catalog
    rows = catalog[catalog["organism"] == "homo_sapiens"].to_dict("records")

    def values(row: dict[str, Any], column: str) -> list[Any] | None:
        if column == "dataset_id":
            return [str(row["dataset_id"])]
        if column == "is_primary_data":
            return PRIMARY.get(str(row.get("is_primary_data")))
        if column in INDEXED and not unrecorded(row.get(column)):
            return [str(v) for v in row[column]]
        return None

    def estimate(q: Any) -> dict[str, Any] | None:
        columns = sorted(q.columns)
        keep = whole = largest = 0
        for row in rows:
            sets = [values(row, c) for c in columns]
            if any(s is None for s in sets):
                return None
            combos = list(itertools.product(*sets))  # type: ignore[arg-type]
            largest = max(largest, len(combos))
            if not combos:
                continue
            data: dict[str, Any] = {}
            for i, c in enumerate(columns):
                column = [t[i] for t in combos]
                if c == "is_primary_data":
                    data[c] = obs_read.Values(np.array(column, dtype=bool), None)
                else:
                    data[c] = obs_read.Elem(pd.Categorical(column))
            kept = len(obs_read.evaluate(q.expression, data, len(combos)))
            keep += kept > 0
            whole += kept == len(combos)
        return {"candidates": keep, "whole": whole, "largest_table": largest}

    items = [(f"filters:{n}", q.build(hs.obs)) for n, q in QUERIES.items()]
    for bound in expand(load_values()):
        if bound.experiment != "homo_sapiens" or not bound.obs_filter:
            continue
        if bound.accessor in ("get_var", "get_presence_matrix"):
            continue
        built, _ = flb.builder_of(str(bound.obs_filter), hs.obs)
        items.append((f"corpus:{bound.label}", built))
    items += [(f"probe:{k}", b(hs.obs)) for k, b, _ in flb.PROBES if b is not None]
    out = {}
    for name, built in items:
        if built is None:
            continue
        q = as_query(built)
        rows_h = census.candidates(q)
        t0 = time.monotonic()
        result = estimate(q)
        record: dict[str, Any] = {
            "hint_candidates": len(rows_h),
            "hint_whole": int(rows_h["whole"].sum()),
        }
        if result is not None:
            record.update(result, ms=round((time.monotonic() - t0) * 1e3, 1))
        out[name] = record
        print(name, record, flush=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
