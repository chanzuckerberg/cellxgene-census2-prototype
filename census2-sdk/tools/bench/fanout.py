"""Count the datasets and the stores that each query reads, from a catalog alone.

The count needs no store, so it runs on a catalog of any size. A consolidation plan
from `consolidate-nf/tools/consolidation-plan.py` stands in for its consolidated stores: each
member of a group reads from one store, the group. A plan thus shows the store
count that a consolidation would give before any store is written.

For each query and each plan, the record holds the datasets that the catalog keeps,
the stores that a read opens, and the datasets that hold a matching cell, when the
cell counts of the catalog can tell (a filter on the counted columns alone). A kept
dataset without a matching cell is a false positive of the catalog.

    python -m bench.fanout catalog.parquet --queries queries.json \\
        --plan size=plan-size.jsonl --plan tissue=plan-tissue.jsonl
"""

import argparse
import json
import pathlib
from typing import Any

import numpy as np
import pandas as pd

from bench.run import ORGANISM, pick_cell_type, serving_uris, source_rows
from census2_sdk import catalog as _catalog


def with_plan(rows: pd.DataFrame, groups: list[dict[str, Any]]) -> pd.DataFrame:
    """The catalog rows, each member of a group read from the group's store."""
    store = {
        m["dataset_version_id"]: f"plan:{g['group']}"
        for g in groups
        for m in g["members"]
    }
    out = rows.copy()
    out[_catalog.CONSOLIDATED[0]] = [
        store.get(str(v)) for v in out["dataset_version_id"]
    ]
    return out


def query_filter(query: dict[str, Any], rows: pd.DataFrame) -> dict[str, Any]:
    if query["shape"] == "single":
        return {"dataset_id": query.get("dataset_id", str(rows["dataset_id"].iloc[0]))}
    if query["shape"] == "cell_type":
        return {
            "cell_type": query.get("cell_type") or pick_cell_type(source_rows(rows))
        }
    if query["shape"] == "filter":
        return dict(query["filter"])
    return {}


def matching(counts: pd.DataFrame, versions: pd.Series, canonical: Any) -> int | None:
    """The kept datasets that hold a matching cell, or None if the counts cannot tell."""
    if set(canonical) - {"dataset_id", *_catalog.COUNTED}:
        return None
    rows = counts[counts["dataset_version_id"].isin(versions)]
    keep = np.ones(len(rows), dtype=bool)
    for column, wanted in canonical.items():
        if column == "dataset_id":
            continue
        values = rows[column].astype(object)
        keep &= np.array(
            [
                not pd.isna(v) and bool(wanted & set(_catalog.members(str(v), column)))
                for v in values
            ],
            dtype=bool,
        )
    return int(rows.loc[keep, "dataset_version_id"].nunique())


def fanout(
    path: pathlib.Path,
    queries: list[dict[str, Any]],
    plans: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    content = _catalog.read(path)
    rows = _catalog.load(content)
    rows = rows[rows["organism"] == ORGANISM]
    try:
        counts = _catalog.load_counts(content)
    except ValueError:
        counts = None
    records = []
    for query in queries:
        obs_filter = query_filter(query, rows)
        canonical = _catalog.canonicalize(obs_filter)
        kept = _catalog.exclude(rows, ORGANISM, canonical) if canonical else rows
        hits = (
            None
            if counts is None
            else matching(counts, kept["dataset_version_id"], canonical)
        )
        record: dict[str, Any] = {
            "query": query["query"],
            "filter": obs_filter,
            "datasets_in_catalog": len(rows),
            "datasets": len(kept),
            "datasets_with_cells": hits,
            "stores": {"per_dataset": len(serving_uris(kept, per_dataset=True))},
        }
        for name, groups in plans.items():
            planned = with_plan(kept, groups)
            record["stores"][name] = len(serving_uris(planned))
        record["stores_in_catalog"] = {"per_dataset": len(rows)} | {
            name: len(serving_uris(with_plan(rows, groups)))
            for name, groups in plans.items()
        }
        records.append(record)
    return records


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("catalog", type=pathlib.Path)
    p.add_argument("--queries", type=pathlib.Path, required=True)
    p.add_argument("--plan", action="append", default=[], help="NAME=plan.jsonl")
    p.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    plans = {}
    for item in args.plan:
        name, _, path = item.partition("=")
        lines = pathlib.Path(path).read_text().splitlines()
        plans[name] = [json.loads(line) for line in lines if line.strip()]
    records = fanout(args.catalog, json.loads(args.queries.read_text()), plans)
    for r in records:
        stores = " ".join(f"{k}={v}" for k, v in r["stores"].items())
        print(
            f"{r['query']:<16} datasets {r['datasets']:>5} "
            f"(with cells {r['datasets_with_cells']}) stores {stores}"
        )
    if args.out:
        args.out.write_text("".join(json.dumps(r) + "\n" for r in records))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
