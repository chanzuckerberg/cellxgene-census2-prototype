#!/usr/bin/env python3
"""Plan a consolidation: put the small per-dataset stores of a catalog into groups.

The plan is JSONL, one group per line, and each group becomes one consolidated store.
The planner reads a catalog that `census2-sdk/tools/build_catalog.py` wrote, and the
size of each store from a listing of its objects. It runs outside the pipeline.

A store is small when it has fewer than `--small-mib` MiB. The planner packs the
small stores into groups of at most `--target-mib` MiB. A large store stays alone,
so no consolidated store is larger than a large per-dataset store, and a query can
still read many stores at the same time. A group with fewer than `--min-members`
members is dropped, and its stores stay alone.

The members of a group must compose: one organism, one gene reference version, and
all spatial or all non-spatial. The planner packs each such partition apart. They
must also share one access grant, because a consolidated store has one. A catalog is
one access domain, so plan only from a public catalog.

Two strategies order the packing:

- `size` packs first-fit decreasing by store size, which gives the fewest groups.
- `tissue` sorts the stores by their main tissue (the tissue with the most cells)
  and fills one group at a time in that order. Stores of one tissue then share a
  consolidated store, so a tissue query opens fewer stores.

The name of a group is `c-` and the first 16 hexadecimal digits of the SHA-256 of
its sorted member dataset versions, so one set of members has one name.

    tools/consolidation-plan.py catalog.parquet --small-mib 256 --target-mib 1024 \\
        --strategy size --out plan.jsonl

The listing needs the `s3` extra of census2-sdk, for obstore, and credentials in the
environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import pyarrow.parquet as pq

MiB = 2**20


def store_bytes(uri: str) -> int:
    """Total bytes of a store, from a listing of its objects."""
    import obstore
    import obstore.store

    url = uri if "://" in uri else "file://" + os.path.abspath(uri)
    store = obstore.store.from_url(url)
    return sum(int(item["size"]) for page in obstore.list(store) for item in page)


def main_tissues(path: str) -> dict[str, str]:
    """The tissue with the most cells in each dataset version, from the cell counts."""
    table = pq.read_table(path, columns=["dataset_version_id", "cell_counts"])
    frame = table.to_pandas()
    rows = []
    for version, entries in zip(frame["dataset_version_id"], frame["cell_counts"]):
        for entry in entries:
            rows.append((version, entry["tissue_ontology_term_id"], entry["n_cells"]))
    counts = pd.DataFrame(rows, columns=["version", "tissue", "n"])
    totals = (
        counts.groupby(["version", "tissue"], dropna=False)["n"].sum().reset_index()
    )
    best = totals.sort_values(["version", "n", "tissue"], ascending=[True, False, True])
    return {
        str(v): str(t)
        for v, t in best.drop_duplicates("version")[["version", "tissue"]].itertuples(
            index=False
        )
    }


def first_fit_decreasing(stores: list[dict], target: int) -> list[list[dict]]:
    groups: list[list[dict]] = []
    sizes: list[int] = []
    for store in sorted(stores, key=lambda s: (-s["bytes"], s["dataset_version_id"])):
        for k, size in enumerate(sizes):
            if size + store["bytes"] <= target:
                groups[k].append(store)
                sizes[k] += store["bytes"]
                break
        else:
            groups.append([store])
            sizes.append(store["bytes"])
    return groups


def by_tissue(stores: list[dict], target: int) -> list[list[dict]]:
    groups: list[list[dict]] = []
    size = target + 1
    for store in sorted(stores, key=lambda s: (s["tissue"], s["dataset_version_id"])):
        if size + store["bytes"] > target:
            groups.append([])
            size = 0
        groups[-1].append(store)
        size += store["bytes"]
    return groups


STRATEGIES: dict[str, Callable[[list[dict], int], list[list[dict]]]] = {
    "size": first_fit_decreasing,
    "tissue": by_tissue,
}


def group_name(members: list[dict]) -> str:
    versions = "\n".join(sorted(m["dataset_version_id"] for m in members))
    return "c-" + hashlib.sha256(versions.encode()).hexdigest()[:16]


def plan(
    catalog: pd.DataFrame,
    tissues: dict[str, str],
    *,
    small: int,
    target: int,
    min_members: int,
    strategy: str,
) -> list[dict]:
    """The groups of the plan, from catalog rows that hold `bytes`."""
    out = []
    keys = ["organism", "gene_reference_version", "spatial"]
    for key, rows in catalog.groupby(keys, sort=True):
        organism, reference, spatial = key
        stores = [
            {
                "dataset_id": str(r.dataset_id),
                "dataset_version_id": str(r.dataset_version_id),
                "store_uri": str(r.store_uri),
                "n_obs": int(r.n_obs),
                "bytes": int(r.bytes),
                "tissue": tissues.get(str(r.dataset_version_id), ""),
            }
            for r in rows.itertuples()
            if int(r.bytes) < small
        ]
        for members in STRATEGIES[strategy](stores, target):
            if len(members) < min_members:
                continue
            out.append(
                {
                    "group": group_name(members),
                    "strategy": strategy,
                    "organism": str(organism),
                    "gene_reference_version": str(reference),
                    "spatial": bool(spatial),
                    "n_obs": sum(m["n_obs"] for m in members),
                    "bytes": sum(m["bytes"] for m in members),
                    "members": members,
                }
            )
    return out


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("catalog", help="a local catalog from build_catalog.py")
    p.add_argument("--small-mib", type=float, default=256)
    p.add_argument("--target-mib", type=float, default=1024)
    p.add_argument("--min-members", type=int, default=2)
    p.add_argument("--strategy", choices=sorted(STRATEGIES), default="size")
    p.add_argument(
        "--sizes",
        default=None,
        help="a JSON file of store bytes by store_uri; the listing writes it when absent",
    )
    p.add_argument(
        "--workers", type=int, default=16, help="stores listed at the same time"
    )
    p.add_argument("--out", required=True)
    args = p.parse_args()

    catalog = pd.read_parquet(
        args.catalog,
        columns=[
            "dataset_id",
            "dataset_version_id",
            "store_uri",
            "n_obs",
            "organism",
            "gene_reference_version",
            "spatial",
        ],
    )
    sizes: dict[str, int] = {}
    if args.sizes and os.path.exists(args.sizes):
        with open(args.sizes) as fh:
            sizes = json.load(fh)
    missing = [u for u in catalog["store_uri"] if u not in sizes]
    if missing:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            sizes.update(zip(missing, pool.map(store_bytes, missing), strict=True))
        if args.sizes:
            with open(args.sizes, "w") as fh:
                json.dump(sizes, fh, indent=1)
    catalog["bytes"] = [sizes[u] for u in catalog["store_uri"]]

    groups = plan(
        catalog,
        main_tissues(args.catalog),
        small=int(args.small_mib * MiB),
        target=int(args.target_mib * MiB),
        min_members=args.min_members,
        strategy=args.strategy,
    )
    with open(args.out, "w") as fh:
        fh.writelines(json.dumps(group) + "\n" for group in groups)

    members = sum(len(g["members"]) for g in groups)
    after = len(catalog) - members + len(groups)
    print(
        f"{len(catalog)} stores -> {after} stores: {members} small stores in "
        f"{len(groups)} consolidated stores ({args.strategy}), {args.out}",
        file=sys.stderr,
    )
    for g in groups:
        print(
            f"  {g['group']}: {len(g['members']):3d} members, {g['n_obs']:>9,} cells, "
            f"{g['bytes'] / MiB:7.1f} MiB",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
