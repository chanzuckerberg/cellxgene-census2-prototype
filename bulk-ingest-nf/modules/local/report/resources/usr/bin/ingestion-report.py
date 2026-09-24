#!/usr/bin/env python3
"""Read the ingestion receipts and report on the run as a whole.

Reads every receipt JSON named on the command line, writes one combined report,
and prints a short version to stderr. It opens no store and no h5ad: every number
comes from the receipts.

Exits non-zero when a receipt records a failed round-trip check, which is what
fails the pipeline run. Uses only the standard library.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # For the shape only. This script stays runnable with the standard library
    # alone, so that a receipt can be read anywhere, not just in the image.
    from census2_dataset_writer import Receipt


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("receipts", nargs="+")
    p.add_argument("--out", default="ingestion_report.json")
    args = p.parse_args()

    receipts: list[Receipt] = []
    for path in args.receipts:
        with open(path) as fh:
            receipts.append(json.load(fh))

    # Not admitted is a normal outcome, so it is reported and not a failure.
    admitted = [r for r in receipts if r.get("admitted", True)]
    rejected = [r for r in receipts if not r.get("admitted", True)]
    checked = [r for r in admitted if "check" in r]
    failed = [r for r in checked if not r["check"].get("passed")]
    dropped = [r for r in admitted if r.get("dropped_obs_columns")]

    report = {
        "recipe": receipts[0]["recipe"] if receipts else None,
        "datasets": len(receipts),
        "admitted": len(admitted),
        "not_admitted": len(rejected),
        "not_admitted_reasons": sorted(
            {r.get("reason") or "unstated" for r in rejected}
        ),
        "source_bytes": sum(r["source_filesize"] for r in receipts),
        "output_bytes": sum(r.get("output_bytes") or 0 for r in admitted),
        "cells": sum(r.get("cell_count") or 0 for r in admitted),
        "elapsed_sec": round(sum(r.get("elapsed_sec") or 0 for r in receipts), 1),
        "slot_policy": sorted(
            {r["slot_policy"] for r in admitted if "slot_policy" in r}
        ),
        # More than one value in either list means the tasks ran different code.
        "census2_dataset_writer_version": sorted(
            {r["census2_dataset_writer_version"] for r in receipts}
        ),
        "anndata_version": sorted({r["anndata_version"] for r in receipts}),
        "zarr_version": sorted({r["zarr_version"] for r in receipts}),
        "checked": len(checked),
        "check_failures": [r["dataset_version_id"] for r in failed],
        # Harmonization drops an obs column that Zarr cannot name. Reported,
        # because a store that quietly holds less than its dataset should not
        # be invisible in the run summary.
        "datasets_with_dropped_obs_columns": len(dropped),
        "dropped_obs_columns": sorted(
            {c for r in dropped for c in r["dropped_obs_columns"]}
        ),
        "receipts": receipts,
    }

    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=1)

    src, out = report["source_bytes"], report["output_bytes"]
    ratio = (out / src) if src else 0
    print(
        f"report: {report['datasets']} datasets seen, "
        f"{report['admitted']} admitted, {report['not_admitted']} not admitted",
        file=sys.stderr,
    )
    for reason in report["not_admitted_reasons"]:
        n = sum(1 for r in rejected if (r.get("reason") or "unstated") == reason)
        print(f"report:   not admitted x{n}: {reason}", file=sys.stderr)
    print(f"report: {report['cells']:,} cells admitted", file=sys.stderr)
    print(
        f"report: {src / 1024**2:.1f} MiB h5ad -> {out / 1024**2:.1f} MiB zarr ({ratio:.2f}x)",
        file=sys.stderr,
    )
    print(f"report: {report['elapsed_sec']}s total task time", file=sys.stderr)
    if checked:
        print(
            f"report: round-trip checked {len(checked)}, failed {len(failed)}",
            file=sys.stderr,
        )
    if dropped:
        print(
            f"report: {len(dropped)} datasets had obs columns dropped: "
            f"{', '.join(report['dropped_obs_columns'])}",
            file=sys.stderr,
        )
    for vid in report["check_failures"]:
        print(f"report: CHECK FAILED {vid}", file=sys.stderr)

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
