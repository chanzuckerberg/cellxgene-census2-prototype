"""Tables and a cost model of file planning, from `bench.sql` records at several sizes.

Each argument is `<mode>=<records.jsonl>`, one file for each way of finding the
fragments. For each query and size, the script prints the median seconds, the
requests and MB of every object read except the catalog, and the files that the scan
read. The cost model is the least-squares slope of the median seconds and of the
requests against the number of catalog rows N, per 1,000 rows, with the median over
the queries.

    python -m bench.plan_cost catalog=ec2-plan-catalog.jsonl footers=ec2-plan-footers.jsonl
"""

import argparse
import json
import statistics
from collections import defaultdict
from typing import Any


def name(record: dict[str, Any]) -> str:
    label = str(record.get("label") or record["shape"])
    return f"count {label}" if record["shape"] == "count" else label


def summary(records: list[dict[str, Any]]) -> dict[str, float] | None:
    ok = [r for r in records if r.get("status") == "ok"]
    if not ok:
        return None
    first = ok[0]
    return {
        "seconds": statistics.median(float(r["query_seconds"]) for r in ok),
        "requests": sum(
            v for k, v in first["requests_by_kind"].items() if k != "catalog"
        ),
        "mb": sum(v for k, v in first["mb_by_kind"].items() if k != "catalog"),
        "files": float(first.get("files_read", 0)),
        "cpu": statistics.median(float(r["cpu_seconds"]) for r in ok),
        "peak_gib": statistics.median(float(r["peak_rss_mb"]) for r in ok) / 1024,
    }


def slope(points: list[tuple[float, float]]) -> float | None:
    if len(points) < 2:
        return None
    mx = statistics.mean(x for x, _ in points)
    my = statistics.mean(y for _, y in points)
    den = sum((x - mx) ** 2 for x, _ in points)
    return sum((x - mx) * (y - my) for x, y in points) / den if den else None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("files", nargs="+", help="<mode>=<records.jsonl>")
    args = p.parse_args()
    groups: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    modes: list[str] = []
    for arg in args.files:
        mode, path = arg.split("=", 1)
        modes.append(mode)
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    record = json.loads(line)
                    groups[(mode, int(record["size"]), name(record))].append(record)
    sizes = sorted({k[1] for k in groups})
    queries = list(dict.fromkeys(k[2] for k in groups))
    for query in queries:
        print(f"\n{query}")
        print("| N | " + " | ".join(modes) + " |")
        print("|---|" + "---|" * len(modes))
        for n in sizes:
            cells = []
            for mode in modes:
                s = summary(groups.get((mode, n, query), []))
                cells.append(
                    f"{s['seconds']:.2f} s, {s['requests']:.0f} req, {s['mb']:.0f} MB"
                    if s
                    else "–"
                )
            print(f"| {n} | " + " | ".join(cells) + " |")
    print("\nper 1,000 catalog rows: seconds, requests (median slope over the queries)")
    for mode in modes:
        seconds, requests = [], []
        for query in queries:
            stats = [(n, summary(groups.get((mode, n, query), []))) for n in sizes]
            have = [(float(n), s) for n, s in stats if s is not None]
            a = slope([(n, s["seconds"]) for n, s in have])
            b = slope([(n, s["requests"]) for n, s in have])
            if a is not None and b is not None:
                seconds.append(a * 1000)
                requests.append(b * 1000)
        if seconds:
            print(
                f"{mode}: {statistics.median(seconds):.2f} s, "
                f"{statistics.median(requests):.0f} requests"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
