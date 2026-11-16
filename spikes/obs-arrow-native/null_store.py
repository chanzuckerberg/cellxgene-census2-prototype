"""Show the null semantics of the query builder on a real store with nulls.

No column of the store format holds a null, so the script makes a test store, as
`null_store.py` of the earlier filter spike does: a copy of one store with an extra
`obs` column, `author_cell_type`, which is the cell type with every fifth cell null.
It writes the copy and a one-row catalog under `OUT`. For terms that a null meets, it
compares the SDK's cells with two rules over the copy's whole `obs`: SQL's
three-valued logic, computed with pandas' nullable booleans, and the pandas rule of a
mask, which the earlier spike follows. Only `obs["author_cell_type"]` names the
column, because the static schema does not.

    python null_store.py s3://my-bucket/stores/<dataset_version_id>.zarr null-store \\
        --tools census2-sdk/tools --out-json null-rule.json
"""

import argparse
import json
import pathlib
import sys
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem

COLUMN = "author_cell_type"


def make(source: str, out: pathlib.Path) -> tuple[pathlib.Path, pd.DataFrame]:
    import build_catalog
    from census2_sdk.catalog import open_store

    group = open_store(source)
    adata = cast("anndata.AnnData", read_elem(group))
    obs = cast("pd.DataFrame", adata.obs)
    labels = obs["cell_type"].astype(str).to_numpy(dtype=object)
    labels[::5] = None
    obs[COLUMN] = pd.Categorical(labels)
    prefix = out / "out"
    prefix.mkdir(parents=True, exist_ok=True)
    path = prefix / source.rstrip("/").rsplit("/", 1)[-1]
    adata.write_zarr(path)
    zarr.open_group(path, mode="r+").update_attributes(dict(group.attrs))
    catalog = out / "catalog.parquet"
    build_catalog.build(str(prefix), catalog)
    frame = cast("pd.DataFrame", read_elem(zarr.open_group(path, mode="r")["obs"]))
    return catalog, frame


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("source")
    p.add_argument("out", type=pathlib.Path)
    p.add_argument("--tools", type=pathlib.Path, required=True)
    p.add_argument("--out-json", type=pathlib.Path, default=None)
    args = p.parse_args()
    sys.path.insert(0, str(args.tools))
    import census2_sdk

    catalog, frame = make(args.source, args.out)
    census = census2_sdk.open_census([catalog])
    hs = census.homo_sapiens
    col = hs.obs[COLUMN]
    label = frame[COLUMN]
    common = str(label.value_counts().index[0])
    # SQL's logic: pandas' nullable booleans follow Kleene logic, and a null keeps
    # no cell. The pandas rule of a mask: a comparison with a null is False.
    sql_eq = label.astype("string") == common
    primary = frame["is_primary_data"].astype("boolean")
    cases: list[tuple[Any, Any, Any]] = [
        (~(col == common), ~sql_eq, ~(label == common)),
        (col != common, sql_eq.map(lambda v: pd.NA if v is pd.NA else not v).astype("boolean"), label != common),
        (col == common, sql_eq, label == common),
        (~col.isin([common]), ~sql_eq, ~label.isin([common])),
        (~((col == common) & hs.obs.is_primary_data), ~(sql_eq & primary), ~((label == common) & frame["is_primary_data"])),
        ((col == common) | hs.obs.is_primary_data, sql_eq | primary, (label == common) | frame["is_primary_data"]),
    ]  # fmt: skip
    rows = []
    for built, sql, pandas_rule in cases:
        got = census._get_obs(built, columns=["dataset_id"])
        keys = list(got["observation_joinid"].astype(str))
        sql_mask = np.asarray(pd.array(sql, dtype="boolean").fillna(False), dtype=bool)
        want = list(frame["observation_joinid"][sql_mask].astype(str))
        rows.append(
            {
                "query": built.text,
                "cells": len(keys),
                "sql_rule": int(sql_mask.sum()),
                "pandas_rule": int(np.asarray(pandas_rule, dtype=bool).sum()),
                "same_as_sql": keys == want,
            }
        )
        print(rows[-1], flush=True)
    result = {
        "cells": len(frame),
        "nulls": int(label.isna().sum()),
        "is_primary_data": sorted({bool(v) for v in frame["is_primary_data"]}),
        "terms": rows,
    }
    if args.out_json is not None:
        args.out_json.write_text(json.dumps(result, indent=1) + "\n")
    return 0 if all(r["same_as_sql"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
