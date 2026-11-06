"""Count the datasets that pass a two-term catalog test but hold no matching cell.

The catalog records the value set of each column for each dataset, so it tests each
term of a query alone. A query `tissue = T AND cell_type = C` keeps a dataset that holds
some cell of tissue T and some cell of type C, also when no cell is both. The SDK then
opens that store for nothing; the obs table does not. This script reads the `tissue`
and `cell_type` columns of the obs table and counts, for each pair (T, C) that the
catalog test can pass, the datasets that pass it (`passing`) and the datasets that
hold a cell of the pair (`joint`).

    python -m bench.false_positives catalog.parquet s3://bucket/census/obs/homo_sapiens \\
        --out results/false_positives.json
"""

import argparse
import collections
import json
import pathlib
from typing import Any

import pandas as pd

from census2_sdk import catalog as _catalog
from census2_sdk import expressions, obs_table

ORGANISM = "homo_sapiens"


def pairs(catalog_path: str, table_location: str, a: str, b: str) -> pd.DataFrame:
    """One row for each pair of values of columns `a` and `b`: passing and joint."""
    rows = _catalog.load(_catalog.read(catalog_path))
    rows = rows[rows["organism"] == ORGANISM]
    table = obs_table.open_table(table_location)
    arrow, _ = obs_table.scan(
        table,
        expressions.from_filter(None),
        ("dataset_id", a, b),
        dictionary=("dataset_id", a, b),
    )
    cells = arrow.to_pandas()
    joint = cells.groupby([a, b], observed=True)["dataset_id"].nunique()
    sets = {
        str(r["dataset_id"]): (set(r[a]), set(r[b]))
        for r in rows.to_dict("records")
        if r[a] is not None and r[b] is not None
    }
    # Index the datasets by each value of `a`, so that a pair counts its passing
    # datasets without a scan of all datasets.
    by_a: dict[str, set[str]] = collections.defaultdict(set)
    by_b: dict[str, set[str]] = collections.defaultdict(set)
    for dataset, (values_a, values_b) in sets.items():
        for v in values_a:
            by_a[v].add(dataset)
        for v in values_b:
            by_b[v].add(dataset)
    out = []
    for value_a, datasets_a in by_a.items():
        for value_b, datasets_b in by_b.items():
            passing = len(datasets_a & datasets_b)
            if passing:
                found = int(joint.get((value_a, value_b), 0))
                out.append((value_a, value_b, passing, found))
    frame = pd.DataFrame(out, columns=[a, b, "passing", "joint"])
    frame["false_positive"] = frame["passing"] - frame["joint"]
    return frame


def summary(frame: pd.DataFrame) -> dict[str, Any]:
    """The totals over all pairs, and over the pairs that a cell holds."""
    held = frame[frame["joint"] > 0]
    worst = frame.sort_values(["false_positive", "passing"], ascending=False).head(10)
    return {
        "pairs": len(frame),
        "pairs_without_a_cell": int((frame["joint"] == 0).sum()),
        "passing": int(frame["passing"].sum()),
        "false_positive": int(frame["false_positive"].sum()),
        "false_positive_share": round(
            float(frame["false_positive"].sum() / frame["passing"].sum()), 3
        ),
        "held_pairs": len(held),
        "held_passing": int(held["passing"].sum()),
        "held_false_positive": int(held["false_positive"].sum()),
        "held_false_positive_share": round(
            float(held["false_positive"].sum() / max(1, held["passing"].sum())), 3
        ),
        "median_false_positive_share_of_held_pairs": round(
            float((held["false_positive"] / held["passing"]).median()), 3
        ),
        "worst": worst.to_dict("records"),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("table")
    p.add_argument("--columns", default="tissue,cell_type")
    p.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    a, b = args.columns.split(",")
    result = summary(pairs(args.catalog, args.table, a, b))
    print(json.dumps(result, indent=1))
    if args.out is not None:
        args.out.write_text(json.dumps(result, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
