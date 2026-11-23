"""Medians of the bench records, by configuration, open, processes, query, and threads.

Usage: python medians.py RUNS.jsonl [WHERE]

A record with `"sdk_categories": "serial"` ran before the SDK read the categories of
a consolidated store at the same time; such records form groups of their own, with
`*` after the configuration. An in-region record shows its batch after a `/`. The
in-region selective queries ran in both batches; batch 2 ran them after the change.
"""

import json
import pathlib
import statistics
import sys
from collections import defaultdict
from typing import Any

path = sys.argv[1]
where = sys.argv[2] if len(sys.argv) > 2 else None
groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
bad = []
for line in pathlib.Path(path).read_text().splitlines():
    r = json.loads(line)
    if where and r.get("where") != where:
        continue
    config = str(r.get("config")) + ("*" if r.get("sdk_categories") else "")
    if "batch" in r:
        config += f"/b{r['batch']}"
    key = (
        config,
        r.get("reader", "sdk"),
        r.get("processes", 1),
        r.get("query"),
        r.get("workers"),
    )
    if r.get("status") != "ok":
        bad.append((key, r.get("status"), r.get("peak_rss_mb")))
        continue
    groups[key].append(r)


def med(rs: list[dict[str, Any]], k: str) -> float:
    return statistics.median(float(r[k]) for r in rs if k in r)


print(
    f"{'config':<18}{'reader':<6}{'P':>2} {'query':<16}{'w':>3} {'n':>2} "
    f"{'wall':>7} {'min':>7} {'max':>7} {'cpu':>7} {'req':>6} {'MB':>7} {'rssGB':>6} "
    f"{'stores':>6} {'cells':>9}"
)
for key in sorted(groups, key=lambda k: (k[3] or "", k[4] or 0, k[0], k[1], k[2])):
    rs = groups[key]
    walls = [float(r["query_seconds"]) for r in rs]
    print(
        f"{key[0]:<18}{key[1]:<6}{key[2]:>2} {key[3]:<16}{key[4]:>3} {len(rs):>2} "
        f"{med(rs, 'query_seconds'):>7.2f} {min(walls):>7.2f} {max(walls):>7.2f} "
        f"{med(rs, 'cpu_seconds'):>7.2f} {med(rs, 'requests'):>6.0f} "
        f"{med(rs, 'mb_payload'):>7.1f} {med(rs, 'peak_rss_mb') / 1024:>6.2f} "
        f"{rs[0]['stores_opened']:>6} {rs[0]['cells']:>9}"
    )
for b in bad:
    print("not ok:", b)
