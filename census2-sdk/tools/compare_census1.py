#!/usr/bin/env python3
"""Compare one query through the SDK with the same query through Census 1.

Census 1 is the oracle: release 2025-11-08 through `cellxgene-census` 1.18.0. The
validation harness uses the same release and SDK version, but it also pins
`tiledbsoma` 1.15.3 and `anndata` 0.13.2, and this command does not. A store written from a
re-curated dataset version can hold other cells, so a cell-set difference is a finding
to record and not a failure. Values are compared on the cells and genes that both sides
hold. The exit code is 0 when the values on the intersection are equal.

    uv run --python 3.12 --with cellxgene-census==1.18.0 --with './census2-sdk[s3]' \\
        census2-sdk/tools/compare_census1.py --catalog catalog.parquet \\
        --dataset-id <dataset_id> [--tissue lung] [--cell-type macrophage] \\
        [--out census2-sdk/tools/results/compare_<dataset_id>.json]

`--out` writes the figures as JSON. `tools/results/` holds the runs that the PR cites.

Needs AWS credentials in the environment for the stores (obstore reads no profile):

    eval "$(aws configure export-credentials --profile <profile> --format env)"
"""

import argparse
import json
import sys
import time

import anndata
import cellxgene_census
import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma

import census2_sdk

CENSUS_VERSION = "2025-11-08"


def census1(dataset_id: str, obs_filter: dict[str, str]) -> anndata.AnnData:
    # repr() gives a Python literal, which is what a SOMA value_filter parses.
    clauses = [f"dataset_id == {dataset_id!r}"]
    clauses += [f"{column} == {value!r}" for column, value in obs_filter.items()]
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        return cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=" and ".join(clauses),
            X_name="raw",
        )


def keyed(adata: anndata.AnnData) -> tuple[pd.Index, pd.Index, sp.csr_matrix]:
    obs = pd.DataFrame(adata.obs)
    cells = pd.Index(
        obs["dataset_id"].astype(str) + "_" + obs["observation_joinid"].astype(str)
    )
    var = pd.DataFrame(adata.var)
    genes = pd.Index(
        var["feature_id"].astype(str) if "feature_id" in var else adata.var_names
    )
    return cells, genes, sp.csr_matrix(adata.X)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--tissue")
    parser.add_argument("--cell-type", dest="cell_type")
    parser.add_argument("--out", help="write the figures to this JSON file")
    args = parser.parse_args()
    obs_filter = {
        k: v for k, v in (("tissue", args.tissue), ("cell_type", args.cell_type)) if v
    }

    started = time.time()
    census = census2_sdk.open_census([args.catalog])
    ours = census.get_anndata(
        "homo_sapiens", {"dataset_id": args.dataset_id, **obs_filter}
    )
    sdk_seconds = round(time.time() - started, 1)
    if ours is None:
        raise SystemExit(f"the SDK selects no cell of {args.dataset_id}")
    print(f"sdk: {ours.shape} in {sdk_seconds}s")

    started = time.time()
    theirs = census1(args.dataset_id, obs_filter)
    census1_seconds = round(time.time() - started, 1)
    print(
        f"census 1: {theirs.shape} in {census1_seconds}s "
        f"(cellxgene-census {cellxgene_census.__version__}, "
        f"tiledbsoma {tiledbsoma.__version__})"
    )

    our_cells, our_genes, our_x = keyed(ours)
    their_cells, their_genes, their_x = keyed(theirs)
    cells = our_cells.intersection(their_cells)
    genes = our_genes.intersection(their_genes)
    print(
        f"cells: sdk {len(our_cells)}, census 1 {len(their_cells)}, both {len(cells)}"
    )
    print(
        f"genes: sdk {len(our_genes)}, census 1 {len(their_genes)}, both {len(genes)}"
    )

    if not len(cells) or not len(genes):
        print("nothing to compare")
        return 1
    a = our_x[our_cells.get_indexer(cells)][:, our_genes.get_indexer(genes)]
    b = their_x[their_cells.get_indexer(cells)][:, their_genes.get_indexer(genes)]
    diff = (a != b).nnz
    union = (abs(a) + abs(b)).nnz
    only_ours = our_x[our_cells.get_indexer(cells)][
        :, our_genes.get_indexer(our_genes.difference(genes))
    ].nnz
    print(
        f"values differing on the intersection: {diff} of {union} positions "
        f"that are nonzero on either side (sdk nnz {a.nnz}, census 1 nnz {b.nnz})"
    )
    print(f"nonzeros in genes the SDK has and Census 1 lacks: {only_ours}")
    print(f"sdk counts total {a.sum():.0f}, census 1 counts total {b.sum():.0f}")
    max_abs = float(np.abs((a - b).data).max()) if diff else 0.0
    print("max abs difference:", max_abs)
    if args.out:
        summary = {
            "dataset_id": args.dataset_id,
            "obs_filter": obs_filter,
            "census_version": CENSUS_VERSION,
            "cellxgene_census": cellxgene_census.__version__,
            "tiledbsoma": tiledbsoma.__version__,
            "cells": {
                "sdk": len(our_cells),
                "census1": len(their_cells),
                "both": len(cells),
            },
            "genes": {
                "sdk": len(our_genes),
                "census1": len(their_genes),
                "both": len(genes),
            },
            "nonzero_positions_either_side": int(union),
            "values_differing": int(diff),
            "nonzeros_only_in_sdk_genes": int(only_ours),
            "counts_total": {"sdk": float(a.sum()), "census1": float(b.sum())},
            "max_abs_difference": max_abs,
            "seconds": {"sdk": sdk_seconds, "census1": census1_seconds},
        }
        with open(args.out, "w") as handle:
            json.dump(summary, handle, indent=2)
            handle.write("\n")
        print(f"wrote {args.out}")
    return 0 if diff == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
