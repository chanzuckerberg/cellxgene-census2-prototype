#!/usr/bin/env python3
"""Compare one query through the SDK with the same query through Census 1.

Census 1 is the oracle: release 2025-11-08 through `cellxgene-census` 1.18.0. The
validation harness uses the same release and SDK version, but it also pins
`tiledbsoma` 1.15.3 and `anndata` 0.13.2, and this command does not. The Iteration 1
stores come from re-curated dataset versions, so a cell-set difference is a finding to
record and not a failure. Values are compared on the cells and genes that both sides
hold. The exit code is 0 when the values on the intersection are equal.

    uv run --python 3.12 --with cellxgene-census==1.18.0 --with './census2-sdk[s3]' \\
        census2-sdk/tools/compare_census1.py --catalog catalog.parquet \\
        --dataset-id <dataset_id> [--tissue lung] [--cell-type macrophage]

Needs AWS credentials in the environment for the stores (obstore reads no profile):

    eval "$(aws configure export-credentials --profile single-cell-dev --format env)"
"""

from __future__ import annotations

import argparse
import sys
import time

import anndata
import cellxgene_census
import census2_sdk
import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma

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
    args = parser.parse_args()
    obs_filter = {
        k: v for k, v in (("tissue", args.tissue), ("cell_type", args.cell_type)) if v
    }

    started = time.time()
    census = census2_sdk.open_census([args.catalog])
    ours = census.get_anndata(
        "homo_sapiens", {"dataset_id": args.dataset_id, **obs_filter}
    )
    print(f"sdk: {ours.shape} in {time.time() - started:.1f}s")

    started = time.time()
    theirs = census1(args.dataset_id, obs_filter)
    print(
        f"census 1: {theirs.shape} in {time.time() - started:.1f}s "
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
    print(
        "max abs difference:",
        np.abs((a - b).data).max() if diff else 0,
    )
    return 0 if diff == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
