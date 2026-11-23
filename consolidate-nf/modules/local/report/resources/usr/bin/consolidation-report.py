#!/usr/bin/env python3
"""Read the consolidation receipts and report on the run as a whole.

Reads every receipt JSON named on the command line and writes one combined report.
It opens no store. Exits non-zero when a receipt records a failed check, which is
what fails the pipeline run. Uses only the standard library.
"""

from __future__ import annotations

import argparse
import json
import sys


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("receipts", nargs="+")
    p.add_argument("--out", default="consolidation_report.json")
    args = p.parse_args()

    receipts = []
    for path in args.receipts:
        with open(path) as fh:
            receipts.append(json.load(fh))
    checked = [r for r in receipts if "check" in r]
    failed = [r for r in checked if not r["check"]["passed"]]
    report = {
        "groups": len(receipts),
        # Each group replaces its members with one store for a query.
        "members": sum(r["n_members"] for r in receipts),
        "cells": sum(r["n_obs"] for r in receipts),
        "source_bytes": sum(r["source_bytes"] for r in receipts),
        "output_bytes": sum(r["output_bytes"] for r in receipts),
        "consolidate_sec": round(sum(r["consolidate_sec"] for r in receipts), 1),
        "max_elapsed_sec": max(r["elapsed_sec"] for r in receipts),
        "max_peak_rss_mib": max(r.get("peak_rss_mib") or 0 for r in receipts),
        "census2_dataset_writer_version": sorted(
            {r["census2_dataset_writer_version"] for r in receipts}
        ),
        "checked": len(checked),
        "check_failures": [r["group"] for r in failed],
        "receipts": receipts,
    }
    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)

    print(
        f"report: {report['groups']} consolidated stores from {report['members']} "
        f"datasets, {report['cells']:,} cells",
        file=sys.stderr,
    )
    print(
        f"report: {report['source_bytes'] / 2**20:.1f} MiB of member stores -> "
        f"{report['output_bytes'] / 2**20:.1f} MiB consolidated",
        file=sys.stderr,
    )
    for group in report["check_failures"]:
        print(f"report: CHECK FAILED {group}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
