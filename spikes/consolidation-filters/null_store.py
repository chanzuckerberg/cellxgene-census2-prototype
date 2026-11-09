"""Check the pandas null rule of the filter builder on a real store with nulls.

No column of the store format holds a null, so the script makes a test store: a copy
of one store with an extra `obs` column, `author_cell_type`, which is the cell type
with every fifth cell null. It writes the copy and a one-row catalog under `OUT`,
then compares the builder's cells with pandas over the copy's whole `obs`, for terms
that a null meets: a negation, `!=`, `isna`, `notna`, and a negated `isin`. Only
`obs["author_cell_type"]` names the column, because the static schema does not.

Run it from `census2-sdk/tools`, with that directory on `PYTHONPATH`:

    PYTHONPATH=. python ../../spikes/consolidation-filters/null_store.py \\
        s3://my-bucket/stores/<dataset_version_id>.zarr null-store --out-json null-rule.json
"""

import argparse
import json
import pathlib
from typing import Any, cast

import anndata
import build_catalog
import census2_sdk
import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem
from census2_sdk.catalog import open_store

COLUMN = "author_cell_type"


def make(source: str, out: pathlib.Path) -> tuple[pathlib.Path, pd.DataFrame]:
    group = open_store(source)
    adata = cast("anndata.AnnData", read_elem(group))
    obs = cast("pd.DataFrame", adata.obs)
    labels = obs["cell_type"].astype(str).to_numpy(dtype=object)
    labels[::5] = None
    obs[COLUMN] = pd.Categorical(labels)
    prefix = out / "out"
    prefix.mkdir(parents=True, exist_ok=True)
    name = source.rstrip("/").rsplit("/", 1)[-1]
    path = prefix / name
    adata.write_zarr(path)
    root = zarr.open_group(path, mode="r+")
    root.update_attributes(dict(group.attrs))
    catalog = out / "catalog.parquet"
    build_catalog.build([str(prefix)], catalog)
    frame = cast("pd.DataFrame", read_elem(zarr.open_group(path, mode="r")["obs"]))
    return catalog, frame


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("source")
    p.add_argument("out", type=pathlib.Path)
    p.add_argument("--out-json", type=pathlib.Path, default=None)
    args = p.parse_args()
    catalog, frame = make(args.source, args.out)
    census = census2_sdk.open_census([catalog])
    hs = census.homo_sapiens
    col = hs.obs[COLUMN]
    label = frame[COLUMN]
    common = str(label.value_counts().index[0])
    pairs: list[tuple[Any, Any]] = [
        (~(col == common), ~(label == common)),
        (col != common, label != common),
        (col == common, label == common),
        (col.isna(), label.isna()),
        (col.notna(), label.notna()),
        (~col.isin([common]), ~label.isin([common])),
        (
            ~((col == common) & hs.obs.is_primary_data),
            ~((label == common) & frame["is_primary_data"]),
        ),
    ]
    rows = []
    for built, want in pairs:
        got = hs.get_obs(built, columns=["observation_joinid"])
        expected = frame["observation_joinid"][np.asarray(want)].astype(str)
        same = list(got["observation_joinid"].astype(str)) == list(expected)
        rows.append(
            {
                "filter": repr(built),
                "cells": len(got),
                "pandas": int(np.sum(want)),
                "same": same,
            }
        )
        print(rows[-1], flush=True)
    result = {"cells": len(frame), "nulls": int(label.isna().sum()), "terms": rows}
    if args.out_json is not None:
        args.out_json.write_text(json.dumps(result, indent=1) + "\n")
    return 0 if all(r["same"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
