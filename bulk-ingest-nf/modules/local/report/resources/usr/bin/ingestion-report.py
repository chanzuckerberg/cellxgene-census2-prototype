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
    dropped_var = [r for r in admitted if r.get("dropped_var_features")]

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
        # cell_count is the source count from the manifest, so subtract the drops.
        "cells": sum(
            r["cell_count"] - sum(r.get("dropped_cells", {}).values())
            for r in admitted
            if r.get("cell_count") is not None
        ),
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
        # The obs projection drops every source column outside the store's
        # schema, so most datasets drop some. The report gives the count; each
        # receipt names its dataset's columns.
        "datasets_with_dropped_obs_columns": len(dropped),
        # Unlike `cells`, this total includes the rejected datasets.
        "dropped_cells": {
            name: sum(r.get("dropped_cells", {}).get(name, 0) for r in receipts)
            for name in sorted(
                {k for r in receipts for k in r.get("dropped_cells", {})}
            )
        },
        # More than one value means the stores do not concatenate by position.
        "gene_reference_version": sorted(
            {str(r["summary"]["gene_reference_version"]) for r in admitted}
        ),
        # The gene projection drops the features that are not organism genes,
        # such as ERCC spike-ins. The report gives a count and not the list,
        # because one dataset can drop many features.
        "datasets_with_dropped_var_features": len(dropped_var),
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
            f"report: {len(dropped)} datasets had obs columns outside the schema",
            file=sys.stderr,
        )
    if dropped_var:
        print(
            f"report: {len(dropped_var)} datasets had non-gene features dropped",
            file=sys.stderr,
        )
    for vid in report["check_failures"]:
        print(f"report: CHECK FAILED {vid}", file=sys.stderr)

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
