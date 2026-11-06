"""How far the column statistics prune a scan of the obs table, for each filter.

For each filter, the tool prints one JSON line with:

- `catalog_keeps`: the datasets that the catalog keeps (`expressions.can_match`).
- `files` and `files_after_bounds`: the data files of the table, and the files that
  Iceberg's plan keeps by the column bounds in the manifests.
- `row_groups` and `row_groups_after_stats`: the row groups of the files after the
  bounds, and the row groups whose Parquet statistics can match the filter (pyarrow's
  `ParquetFileFragment.subset`).
- `rows`: the rows of those files that match.

A numeric column has exact bounds, so a range on it can exclude a file or a row group.
A label column has bounds that are the first 16 characters of its lowest and highest
label, which exclude a file only when the wanted label sorts outside them.

    python -m bench.pruning catalog-103.parquet s3://my-bucket/spike/tables/obs-103 \\
        --filter "nnz >= 500" --filter "cell_type = 'neutrophil'"
"""

import argparse
import json
import time
from typing import Any, cast

import pandas as pd
import pyarrow.dataset as ds  # type: ignore[import-untyped]  # pyarrow has no types
from pyiceberg.expressions import parser
from pyiceberg.expressions.visitors import bind
from pyiceberg.io.pyarrow import PyArrowFileIO, expression_to_pyarrow

from census2_sdk import catalog as _catalog
from census2_sdk import expressions, obs_table


def prune(catalog: pd.DataFrame, table: Any, text: str) -> dict[str, Any]:
    expr = parser.parse(text)
    records = cast("list[dict[str, Any]]", catalog.to_dict("records"))
    keeps = sum(1 for row in records if expressions.can_match(expr, row))
    t0 = time.monotonic()
    tasks = list(table.scan(row_filter=expr).plan_files())
    plan_seconds = time.monotonic() - t0
    bound = bind(table.schema(), expressions.with_members(expr, records), True)
    filter_ = expression_to_pyarrow(bound)
    io = cast("PyArrowFileIO", table.io)
    groups = kept = rows = 0
    for task in tasks:
        scheme, netloc, path = io.parse_location(task.file.file_path)
        fs = io.fs_by_scheme(scheme, netloc)
        fragment = ds.ParquetFileFormat().make_fragment(path, filesystem=fs)
        fragment.ensure_complete_metadata()
        groups += len(fragment.row_groups)
        subset = fragment.subset(filter=filter_)
        kept += len(subset.row_groups)
        rows += int(subset.count_rows(filter=filter_))
    return {
        "filter": text,
        "catalog_rows": len(records),
        "catalog_keeps": keeps,
        "files": int(table.inspect.files().num_rows),
        "files_after_bounds": len(tasks),
        "plan_seconds": round(plan_seconds, 3),
        "row_groups": groups,
        "row_groups_after_stats": kept,
        "rows": rows,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("catalog")
    p.add_argument("table")
    p.add_argument("--filter", action="append", required=True)
    args = p.parse_args()
    frame = _catalog.load(_catalog.read(args.catalog))
    table = obs_table.open_table(args.table)
    for text in args.filter:
        print(json.dumps(prune(frame, table, text)), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
