"""Compare the cell keys of two paths, query by query.

Each input is a JSON-lines file of `builder_keys.py` or `duckdb_keys.py`. A query is
the same when the number of cells and the digest of the sorted keys agree.

    python compare_keys.py builder-keys.jsonl duckdb-keys.jsonl --out crosscheck.json
"""

import argparse
import json
import pathlib
from typing import Any


def load(path: pathlib.Path) -> dict[str, dict[str, Any]]:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return {r["query"]: r for r in rows}  # the last record of a query wins


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("left", type=pathlib.Path)
    p.add_argument("right", type=pathlib.Path)
    p.add_argument("--out", type=pathlib.Path, default=None)
    args = p.parse_args()
    left, right = load(args.left), load(args.right)
    rows = []
    for name in left:
        if name not in right:
            continue
        a, b = left[name], right[name]
        same = a["cells"] == b["cells"] and a["digest"] == b["digest"]
        rows.append(
            {
                "query": name,
                "cells": a["cells"],
                "datasets": a["datasets"],
                "same": same,
                f"{a['path']}_cells": a["cells"],
                f"{b['path']}_cells": b["cells"],
                f"{a['path']}_candidates": a.get("candidates"),
                f"{b['path']}_candidates": b.get("candidates"),
            }
        )
        mark = "same" if same else "DIFFERENT"
        print(f"{name:<26} {a['cells']:>10} {b['cells']:>10} {mark}")
    print(f"{sum(r['same'] for r in rows)} of {len(rows)} queries give the same keys")
    if args.out is not None:
        args.out.write_text(json.dumps(rows, indent=1) + "\n")
    return 0 if all(r["same"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
