"""Render the tables of `spikes/obs-arrow-native.md` from the result files.

python summarize.py census2-sdk/tools/results/obs_arrow_native
"""

import argparse
import collections
import json
import pathlib
import statistics
from typing import Any


def jsonl(path: pathlib.Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()]


def e2e(results: pathlib.Path) -> None:
    rows = [r for r in jsonl(results / "e2e-ec2.jsonl") if r.get("status") == "ok"]
    groups: dict[tuple[str, int, str], list[dict[str, Any]]] = collections.defaultdict(
        list
    )
    for r in rows:
        groups[(r["query"], r["workers"], r["arm"])].append(r)
    arms = ["main", "pr45", "pr43", "spike"]
    queries = list(dict.fromkeys(r["query"] for r in rows))

    def med(key: tuple[str, int, str], field: str) -> Any:
        values = [r[field] for r in groups.get(key, [])]
        return statistics.median(values) if values else None

    for workers in (8, 32):
        print(f"\n### {workers} workers: seconds (CPU seconds), median of 3\n")
        print(
            "| Query | Cells | Stores opened: mapping / #43 / spike (whole) | "
            + " | ".join(arms)
            + " |"
        )
        print("|---|---|---|" + "---|" * len(arms))
        for q in queries:
            cells = med((q, workers, "spike"), "cells")
            opened = [
                med((q, workers, a), "stores_opened") for a in ("main", "pr43", "spike")
            ]
            whole = med((q, workers, "spike"), "stores_whole")
            cols = []
            for a in arms:
                s, c = (
                    med((q, workers, a), "seconds"),
                    med((q, workers, a), "cpu_seconds"),
                )
                cols.append("–" if s is None else f"{s:.2f} ({c:.1f})")
            op = " / ".join("–" if o is None else str(int(o)) for o in opened)
            print(
                f"| {q} | {int(cells or 0):,} | {op} ({int(whole or 0)}) | "
                + " | ".join(cols)
                + " |"
            )
        print(f"\n### {workers} workers: requests / MB / peak RSS GiB, median of 3\n")
        print("| Query | " + " | ".join(arms) + " |")
        print("|---|" + "---|" * len(arms))
        for q in queries:
            cols = []
            for a in arms:
                req = med((q, workers, a), "requests")
                mb = med((q, workers, a), "mb")
                rss = med((q, workers, a), "peak_rss_mb")
                cols.append(
                    "–"
                    if req is None
                    else f"{int(req):,} / {mb:,.0f} / {rss / 1024:.1f}"
                )
            print(f"| {q} | " + " | ".join(cols) + " |")


def bind(results: pathlib.Path) -> None:
    rows = jsonl(results / "bind-ec2.jsonl")
    if not rows:
        return
    by: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    for r in rows:
        by[(r["query"], r["store"])].append(r)
    print(
        "\n### Per store: medians over the stores of the per-store medians of 3 repeats\n"
    )
    print(
        "| Query | Columns | read_elem(obs) ms | read ms | arrow ms | Table.filter ms | "
        "numpy (#43) ms | arrow share of read + arrow | Sum over stores, s: read_elem / read + arrow |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    queries = list(dict.fromkeys(r["query"] for r in rows))
    for q in queries:
        stores = [s for (qq, s) in by if qq == q]
        per: dict[str, list[float]] = collections.defaultdict(list)
        for s in stores:
            rs = by[(q, s)]
            for m in ("read_elem", "read", "arrow", "filter", "numpy"):
                per[m].append(statistics.median(r[m]["ms"] for r in rs))
        share = [a / (r + a) for r, a in zip(per["read"], per["arrow"], strict=True)]
        same = all(r["same"] for s in stores for r in by[(q, s)])
        ncols = len(by[(q, stores[0])][0]["columns"])
        print(
            f"| {q} | {ncols} | {statistics.median(per['read_elem']):,.0f} | "
            f"{statistics.median(per['read']):.1f} | {statistics.median(per['arrow']):.2f} | "
            f"{statistics.median(per['filter']):.2f} | {statistics.median(per['numpy']):.2f} | "
            f"{100 * statistics.median(share):.1f}% | "
            f"{sum(per['read_elem']) / 1e3:.1f} / "
            f"{(sum(per['read']) + sum(per['arrow'])) / 1e3:.2f} |"
            + ("" if same else " (DIFFERENT CELLS)")
        )
    print("\n### Per store, CPU ms (medians) and the largest store\n")
    print(
        "| Query | read_elem CPU | read CPU | arrow CPU | numpy CPU | largest store: n_obs, arrow ms, numpy ms, filter ms |"
    )
    print("|---|---|---|---|---|---|")
    for q in queries:
        stores = [s for (qq, s) in by if qq == q]
        cpu: dict[str, list[float]] = collections.defaultdict(list)
        for s in stores:
            for m in ("read_elem", "read", "arrow", "numpy"):
                cpu[m].append(statistics.median(r[m]["cpu_ms"] for r in by[(q, s)]))
        big = max(stores, key=lambda s: by[(q, s)][0]["n_obs"])
        rs = by[(q, big)]
        print(
            f"| {q} | {statistics.median(cpu['read_elem']):.1f} | {statistics.median(cpu['read']):.1f} | "
            f"{statistics.median(cpu['arrow']):.2f} | {statistics.median(cpu['numpy']):.2f} | "
            f"{rs[0]['n_obs']:,}, {statistics.median(r['arrow']['ms'] for r in rs):.1f}, "
            f"{statistics.median(r['numpy']['ms'] for r in rs):.1f}, "
            f"{statistics.median(r['filter']['ms'] for r in rs):.1f} |"
        )


def dictionary(results: pathlib.Path) -> None:
    rows = jsonl(results / "dictionary-ec2.jsonl") + jsonl(
        results / "dictionary-laptop.jsonl"
    )
    if not rows:
        return
    cells = {(r["where"], r["categories"], r["measure"]): r["ms"] for r in rows}
    measures = list(dict.fromkeys(r["measure"] for r in rows))
    cols = sorted({(r["where"], r["categories"]) for r in rows})
    print("\n### Dictionary terms, 2 million cells, one thread, best of 7, ms\n")
    print("| Measure | " + " | ".join(f"{w}, {k:,}" for w, k in cols) + " |")
    print("|---|" + "---|" * len(cols))
    for m in measures:
        print(
            f"| {m} | "
            + " | ".join(str(cells.get((w, k, m), "")) for w, k in cols)
            + " |"
        )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", type=pathlib.Path)
    args = p.parse_args()
    dictionary(args.results)
    bind(args.results)
    e2e(args.results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
