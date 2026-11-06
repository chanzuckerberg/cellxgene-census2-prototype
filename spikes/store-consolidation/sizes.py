"""Size distribution before and after consolidation, and where the cells of a cell type are.

Usage: python sizes.py CATALOG_CONSOLIDATED SIZES_JSON RECEIPTS_DIR [CELL_TYPE]

CATALOG_CONSOLIDATED is a catalog from build_catalog.py --consolidated. SIZES_JSON
is the store bytes by store_uri that the planner wrote. RECEIPTS_DIR holds the
receipts of the consolidation run (one JSON per group).
"""

import glob
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MiB = 2**20


def pct(values: np.ndarray, qs: list[float]) -> list[float]:
    return [float(np.percentile(values, q)) for q in qs]


def rank(value: float, others: np.ndarray) -> float:
    """The percentile of `value` among `others`: the share at or below it."""
    return 100.0 * float(np.mean(others <= value))


def main() -> None:
    catalog_path, sizes_path, receipts_dir = sys.argv[1:4]
    cell_type = sys.argv[4] if len(sys.argv) > 4 else "neutrophil"
    rows = pd.read_parquet(
        catalog_path,
        columns=[
            "dataset_id",
            "dataset_version_id",
            "store_uri",
            "n_obs",
            "consolidated_store_uri",
        ],
    )
    sizes = json.loads(pathlib.Path(sizes_path).read_text())
    rows["bytes"] = [sizes[u] for u in rows["store_uri"]]
    receipts = [
        json.loads(pathlib.Path(p).read_text())
        for p in sorted(glob.glob(f"{receipts_dir}/*.json"))
    ]
    receipts = [r for r in receipts if "group" in r and "output_bytes" in r]
    by_store = {r["store"].rstrip("/"): r for r in receipts}

    alone = rows[rows["consolidated_store_uri"].isna()]
    cons = pd.DataFrame(
        [
            {
                "store": s,
                "bytes": r["output_bytes"],
                "n_obs": r["n_obs"],
                "members": r["n_members"],
                "source_bytes": r["source_bytes"],
            }
            for s, r in by_store.items()
        ]
    )
    qs = [0, 25, 50, 75, 90, 100]
    out: dict = {
        "datasets": len(rows),
        "stores_before": len(rows),
        "stores_after": len(alone) + len(cons),
        "consolidated_stores": len(cons),
        "members": int(rows["consolidated_store_uri"].notna().sum()),
        "bytes_total_gib": round(rows["bytes"].sum() / 2**30, 2),
        "cells_total": int(rows["n_obs"].sum()),
        "second_copy_share_of_bytes": round(
            cons["bytes"].sum() / rows["bytes"].sum(), 4
        ),
    }
    before_b = rows["bytes"].to_numpy() / MiB
    after_b = (
        np.concatenate([alone["bytes"].to_numpy(), cons["bytes"].to_numpy()]) / MiB
    )
    before_c = rows["n_obs"].to_numpy()
    after_c = np.concatenate([alone["n_obs"].to_numpy(), cons["n_obs"].to_numpy()])
    out["mib_before_p0_25_50_75_90_100"] = [round(v, 1) for v in pct(before_b, qs)]
    out["mib_after_p0_25_50_75_90_100"] = [round(v, 1) for v in pct(after_b, qs)]
    out["cells_before_p0_25_50_75_90_100"] = [int(v) for v in pct(before_c, qs)]
    out["cells_after_p0_25_50_75_90_100"] = [int(v) for v in pct(after_c, qs)]
    out["small_share_of_stores"] = round(float(np.mean(before_b < 256)), 3)
    out["small_share_of_bytes"] = round(
        float(before_b[before_b < 256].sum() / before_b.sum()), 4
    )
    alone_b = alone["bytes"].to_numpy() / MiB
    alone_c = alone["n_obs"].to_numpy()
    out["consolidated"] = [
        {
            "store": str(r.store).rsplit("/", 1)[1],
            "members": int(r.members),
            "mib": round(r.bytes / MiB, 1),
            "cells": int(r.n_obs),
            "pct_rank_bytes_among_alone": round(rank(r.bytes / MiB, alone_b), 1),
            "pct_rank_cells_among_alone": round(rank(r.n_obs, alone_c), 1),
            "out_over_in": round(r.bytes / r.source_bytes, 3),
        }
        for r in cons.sort_values("bytes").itertuples()
    ]

    # Where the cells of one cell type are, from the catalog's cell counts.
    counts = pq.read_table(
        catalog_path, columns=["dataset_version_id", "cell_counts"]
    ).to_pandas()
    per = {}
    for version, entries in zip(counts["dataset_version_id"], counts["cell_counts"]):
        n = sum(int(e["n_cells"]) for e in entries if e["cell_type"] == cell_type)
        if n:
            per[str(version)] = n
    hit = rows[rows["dataset_version_id"].isin(per)].copy()
    hit["cells"] = [per[v] for v in hit["dataset_version_id"]]
    in_cons = hit[hit["consolidated_store_uri"].notna()]
    out["cell_type"] = {
        "name": cell_type,
        "datasets_with_cells": len(hit),
        "cells": int(hit["cells"].sum()),
        "datasets_alone": int(hit["consolidated_store_uri"].isna().sum()),
        "datasets_in_consolidated": len(in_cons),
        "consolidated_stores": int(in_cons["consolidated_store_uri"].nunique()),
        "share_of_cells_in_consolidated": round(
            float(in_cons["cells"].sum() / hit["cells"].sum()), 4
        ),
        "alone_top4_share": round(
            float(
                hit[hit["consolidated_store_uri"].isna()]["cells"].nlargest(4).sum()
                / hit["cells"].sum()
            ),
            4,
        ),
        "stores_per_dataset_path": len(hit),
        "stores_consolidated_path": int(
            hit["consolidated_store_uri"].isna().sum()
            + in_cons["consolidated_store_uri"].nunique()
        ),
    }
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
