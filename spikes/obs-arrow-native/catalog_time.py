"""Time the catalog test of each query, with and without pyarrow on the value sets.

For each query of the corpus of `cost.py` (the named queries of `bench/filters.py`,
the human queries of the harness corpus, and the probes), the script times
`Census.candidates` as the SDK runs it, and with the hint test only (the value-set
test replaced by a function that tests no row). It records the medians of
`--repeats` runs, the candidates and whole stores of each, the largest table of one
row, the rows of all tables, and the rows that kept the hint test because their
table passed `value_sets.MAX_COMBINATIONS`. It reads only the catalog.

    python catalog_time.py catalog-2060.parquet --tools census2-sdk/tools --out t.json
"""

import argparse
import json
import pathlib
import statistics
import sys
import time
from typing import Any

import numpy as np


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--tools", type=pathlib.Path, required=True)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--out", type=pathlib.Path, required=True)
    args = p.parse_args()
    sys.path.insert(0, str(args.tools))
    import warnings

    import census2_sdk
    from bench.filters import QUERIES
    from census2_sdk import value_sets
    from census2_sdk.query import as_query
    from harness import filter_language_builder as flb
    from harness.corpus import expand, load_values

    census = census2_sdk.open_census([args.catalog])
    hs = census.homo_sapiens
    items = [(f"filters:{n}", q.build(hs.obs)) for n, q in QUERIES.items()]
    for bound in expand(load_values()):
        if bound.experiment != "homo_sapiens" or not bound.obs_filter:
            continue
        if bound.accessor in ("get_var", "get_presence_matrix"):
            continue
        built, _ = flb.builder_of(str(bound.obs_filter), hs.obs)
        items.append((f"corpus:{bound.label}", built))
    items += [(f"probe:{k}", b(hs.obs)) for k, b, _ in flb.PROBES if b is not None]

    real = value_sets.test

    def hints_only(rows: Any, q: Any) -> Any:
        n = len(rows)
        return np.ones(n, bool), np.zeros(n, bool), np.zeros(n, bool)

    def timed(q: Any) -> tuple[float, Any]:
        times, rows = [], None
        for _ in range(args.repeats):
            t0 = time.perf_counter()
            rows = census.candidates(q)
            times.append(time.perf_counter() - t0)
        return round(statistics.median(times) * 1e3, 2), rows

    out: dict[str, Any] = {"catalog_rows": len(census.catalog)}
    for name, built in items:
        q = as_query(hs if built is None else built)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value_sets.test = real
            ms, rows = timed(q)
            value_sets.test = hints_only
            hint_ms, hint_rows = timed(q)
            value_sets.test = real
        records = census.catalog[census.catalog["organism"] == "homo_sapiens"]
        sizes = []
        for row in records.to_dict("records"):
            sets = [value_sets.value_set(row, c) for c in sorted(q.columns)]
            sizes.append(int(np.prod([1 if s is None else len(s) for s in sets])))
        record = {
            "ms": ms,
            "hint_ms": hint_ms,
            "candidates": len(rows),
            "whole": int(rows["whole"].sum()),
            "hint_candidates": len(hint_rows),
            "hint_whole": int(hint_rows["whole"].sum()),
            "null_safe": q.null_safe,
            "largest_table": max(sizes),
            "table_rows": int(
                sum(s for s in sizes if s <= value_sets.MAX_COMBINATIONS)
            ),
            "rows_over_limit": sum(s > value_sets.MAX_COMBINATIONS for s in sizes),
        }
        out[name] = record
        print(name, json.dumps(record), flush=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
