"""Tabulate the whole-obs reads: medians by configuration, organism, and workers.

    python spikes/store-consolidation/obs_summary.py \\
        census2-sdk/tools/results/consolidation-obs-runs.jsonl

The records with `repeat` 0 carry the digest of the frame and are not timed runs.
The script also checks that every Census 2 configuration gives the same digest.
"""

import json
import statistics
import sys
from collections import defaultdict


def main(path: str) -> None:
    with open(path) as fh:
        rows = [json.loads(line) for line in fh if line.strip()]
    timed = defaultdict(list)
    digests = defaultdict(dict)
    for r in rows:
        if r.get("repeat") == 0:
            digests[r["organism"]][r["config"]] = r.get("digest")
            continue
        timed[(r["config"], r["organism"], r.get("workers", "-"))].append(r)

    print(
        "| config | organism | workers | stores | cells | s | CPU s | concat s | requests | peak GiB |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|")
    for (config, organism, workers), runs in sorted(timed.items(), key=str):

        def med(key: str, runs: list[dict] = runs) -> float | None:
            values = [r[key] for r in runs if r.get(key) is not None]
            return statistics.median(values) if values else None

        peak = med("peak_rss_gib") or (med("peak_rss_mb") or 0) / 1024
        print(
            f"| {config} | {organism} | {workers} | {runs[0].get('stores_opened', '-')} "
            f"| {runs[0].get('rows', runs[0].get('cells'))} | {med('query_seconds'):.1f} "
            f"| {med('cpu_seconds'):.0f} | {med('concat_seconds') or 0:.1f} "
            f"| {med('requests') or '-'} | {peak:.1f} |"
        )
    for organism, by_config in digests.items():
        values = {json.dumps(d, sort_keys=True) for d in by_config.values()}
        same = "the same" if len(values) == 1 else "DIFFERENT"
        print(f"{organism}: digests of {sorted(by_config)} are {same}")


if __name__ == "__main__":
    main(sys.argv[1])
