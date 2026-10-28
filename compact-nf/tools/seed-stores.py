#!/usr/bin/env python3
"""Write small per-dataset stores and a plan of one group, for the pipeline tests.

Each store is in the store format: the `obs` schema of `OBS_COLUMNS`, a summary
record, and a shared six-gene `var` in place of the pinned reference. The writer's
`write` checks the record and the `dataset_id` column, as for an ingested store. No
source h5ad and no gene reference are necessary, so the tests run offline.

    tools/seed-stores.py tests/stores

It needs the census2-dataset-writer package, so run it in the runtime image.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import anndata
import census2_dataset_writer
import numpy as np
import pandas as pd
import scipy.sparse as sp

GENES = [f"ENSG{j:011d}" for j in range(6)]


def store(n_obs: int, k: int) -> anndata.AnnData:
    rng = np.random.default_rng(k)
    columns = {}
    for name, col_type in census2_dataset_writer.OBS_COLUMNS:
        if name == census2_dataset_writer.DATASET_ID_COLUMN:
            continue
        if col_type is census2_dataset_writer.ColType.BOOL:
            columns[name] = np.ones(n_obs, dtype=bool)
        elif col_type is census2_dataset_writer.ColType.STRING:
            columns[name] = np.array([f"j{i}" for i in range(n_obs)], dtype=object)
        else:
            values = [f"{name}-{v}" for v in rng.choice(["a", "b", f"only-{k}"], n_obs)]
            columns[name] = pd.Categorical(values, categories=sorted(set(values)))
    obs = pd.DataFrame(columns, index=[f"cell{i}" for i in range(n_obs)])
    x = sp.random(
        n_obs, len(GENES), density=0.5, format="csr", dtype="float32", random_state=k
    )
    var = pd.DataFrame(
        {"feature_name": [f"G{j}" for j in range(len(GENES))]}, index=GENES
    )
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.0.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    census2_dataset_writer.add_dataset_id(adata, f"dataset-{k}")
    return adata


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=pathlib.Path)
    p.add_argument("--sizes", default="40,25,60")
    args = p.parse_args()

    members = []
    for k, n_obs in enumerate(int(s) for s in args.sizes.split(",")):
        dest = args.out / f"version-{k}.zarr"
        summary: census2_dataset_writer.SummaryRecord = {
            "version": census2_dataset_writer.SUMMARY_VERSION,
            "dataset_id": f"dataset-{k}",
            "dataset_version_id": f"version-{k}",
            "gene_reference_version": "TEST v1",
        }
        census2_dataset_writer.write(store(n_obs, k), dest, summary=summary)
        members.append(
            {
                "dataset_id": f"dataset-{k}",
                "dataset_version_id": f"version-{k}",
                "store_uri": str(dest.resolve()),
                "n_obs": n_obs,
                "bytes": 0,
            }
        )
    group = {
        "group": "c-test",
        "strategy": "size",
        "organism": "homo_sapiens",
        "gene_reference_version": "TEST v1",
        "spatial": False,
        "n_obs": sum(m["n_obs"] for m in members),
        "bytes": 0,
        "members": members,
    }
    (args.out / "plan.jsonl").write_text(json.dumps(group) + "\n")
    print(f"{len(members)} stores and plan.jsonl -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
