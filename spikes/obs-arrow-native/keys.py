"""The cell keys of each query of `bench/filters.py`, through this branch's SDK.

For each query, the script selects the cells with `census.get_obs` and records the
number of cells, the datasets, the catalog candidates and the whole ones, and a
SHA-256 of the sorted cell keys, one `dataset_id<TAB>observation_joinid` line each.
The digest is the one of `builder_keys.py` on `spike/consolidation-filters`, so the
records compare with that branch's records. `--compare` compares with such a file.

    python keys.py catalog-103.parquet --tools census2-sdk/tools --out keys.jsonl \\
        --compare builder-keys-iceberg-stores.jsonl
"""

import argparse
import hashlib
import json
import pathlib
import sys
import time
from typing import Any

import numpy as np


def digest(dataset_ids: Any, joinids: Any) -> tuple[str, list[str]]:
    """SHA-256 of the sorted keys, one `dataset_id<TAB>observation_joinid` line each."""
    keys = sorted(
        f"{d}\t{j}" for d, j in zip(np.asarray(dataset_ids), np.asarray(joinids))
    )
    h = hashlib.sha256()
    for key in keys:
        h.update(key.encode())
        h.update(b"\n")
    return h.hexdigest(), keys


def compare(ours: pathlib.Path, theirs: pathlib.Path) -> list[dict[str, Any]]:
    def load(path: pathlib.Path) -> dict[str, dict[str, Any]]:
        lines = path.read_text().splitlines()
        return {r["query"]: r for r in (json.loads(x) for x in lines if x.strip())}

    left, right = load(ours), load(theirs)
    rows = []
    for name, a in left.items():
        b = right.get(name)
        if b is None:
            continue
        same = a["cells"] == b["cells"] and a["digest"] == b["digest"]
        rows.append(
            {
                "query": name,
                "same": same,
                "cells": a["cells"],
                "their_cells": b["cells"],
                "candidates": a["candidates"],
                "their_candidates": b.get("candidates"),
            }
        )
        print(
            f"{name:<24} {a['cells']:>10} {b['cells']:>10} "
            f"{a['candidates']:>4} {b.get('candidates'):>4} "
            f"{'same' if same else 'DIFFERENT'}"
        )
    print(f"{sum(r['same'] for r in rows)} of {len(rows)} queries give the same keys")
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("--tools", type=pathlib.Path, required=True)
    p.add_argument("--queries", default=None, help="comma-separated; default all")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--keys-dir", type=pathlib.Path, default=None)
    p.add_argument("--compare", type=pathlib.Path, default=None)
    p.add_argument("--compare-out", type=pathlib.Path, default=None)
    args = p.parse_args()
    sys.path.insert(0, str(args.tools))
    import census2_sdk
    from bench.filters import QUERIES

    names = args.queries.split(",") if args.queries else list(QUERIES)
    census = census2_sdk.open_census([args.catalog])
    hs = census.homo_sapiens
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for name in names:
            built = QUERIES[name].build(hs.obs)
            q = hs if built is None else built
            t0 = time.monotonic()
            frame = census.get_obs(q, columns=["dataset_id"], max_workers=args.workers)
            seconds = time.monotonic() - t0
            ids = frame["dataset_id"].astype(str) if len(frame) else []
            joins = frame["observation_joinid"].astype(str) if len(frame) else []
            value, keys = digest(ids, joins)
            rows = census.candidates(q)
            record = {
                "query": name,
                "path": "spike",
                "cells": len(keys),
                "datasets": len(set(np.asarray(ids))),
                "digest": value,
                "seconds": round(seconds, 3),
                "candidates": len(rows),
                "whole": int(rows["whole"].sum()),
                "columns": sorted(as_query(q).columns),
                "exact": as_query(q).exact,
            }
            fh.write(json.dumps(record) + "\n")
            fh.flush()
            print(json.dumps(record), flush=True)
            if args.keys_dir is not None:
                args.keys_dir.mkdir(parents=True, exist_ok=True)
                (args.keys_dir / f"{name}.txt").write_text("\n".join(keys) + "\n")
    if args.compare is not None:
        rows = compare(args.out, args.compare)
        if args.compare_out is not None:
            args.compare_out.write_text(json.dumps(rows, indent=1) + "\n")
    return 0


def as_query(q: Any) -> Any:
    from census2_sdk.query import as_query as convert

    return convert(q)


if __name__ == "__main__":
    raise SystemExit(main())
