"""Census 1 get_anndata for a cell_type list, every obs column, one JSON record.

python census1_anndata.py homo_sapiens queries/rare_19.json --repeat 1
"""

import argparse
import importlib.metadata
import json
import resource
import time


def cpu() -> float:
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_utime + u.ru_stime


p = argparse.ArgumentParser()
p.add_argument("organism")
p.add_argument("query")
p.add_argument("--census-version", default="2025-11-08")
p.add_argument("--repeat", type=int, default=0)
args = p.parse_args()

import cellxgene_census

with open(args.query) as fh:
    types = json.load(fh)["cell_type"]
value_filter = "cell_type in [" + ", ".join(json.dumps(t) for t in types) + "]"
t0 = time.monotonic()
census = cellxgene_census.open_soma(census_version=args.census_version)
t1, c1 = time.monotonic(), cpu()
adata = cellxgene_census.get_anndata(
    census, args.organism, obs_value_filter=value_filter
)
t2, c2 = time.monotonic(), cpu()
print(
    json.dumps(
        {
            "config": "census1",
            "shape": "filter",
            "organism": args.organism,
            "repeat": args.repeat,
            "open_seconds": round(t1 - t0, 3),
            "query_seconds": round(t2 - t1, 3),
            "cpu_seconds": round(c2 - c1, 3),
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "nnz": int(adata.X.nnz),
            "obs_columns": len(adata.obs.columns),
            "peak_rss_gib": round(
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 2
            ),
            "versions": {
                m: importlib.metadata.version(m)
                for m in ("cellxgene-census", "tiledbsoma", "pandas")
            },
        }
    ),
    flush=True,
)
