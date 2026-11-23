"""Digest the result of each query, so that results from two SDKs can be compared.

Usage: python digest.py CATALOG QUERIES OUT [--per-dataset] [--cases-from CATALOG2]

QUERIES is the bench queries file. The cases are its selective queries, and one
case for each consolidated store of CATALOG2 (or of CATALOG): all of its members,
read whole. Each case reads with `layers=["normalized"]`. The digest covers the
`obs_names`, every `obs` and `var` column with its dtype (and the categories of a
categorical), `X`, and `layers["normalized"]`, but not `uns`, which names the stores
read. The SDK of `main` has no `consolidated` argument, so a run without
`--per-dataset` calls `open_census` without it.
"""

import hashlib
import json
import pathlib
import sys
import time
from typing import Any

import census2_sdk
import numpy as np
import pandas as pd
import scipy.sparse as sp


def h(*parts: bytes) -> str:
    d = hashlib.sha256()
    for p in parts:
        d.update(p)
        d.update(b"\0")
    return d.hexdigest()[:16]


def frame_digest(frame: pd.DataFrame) -> dict[str, str]:
    out = {"index": h("\n".join(map(str, frame.index)).encode())}
    for name in frame.columns:
        col = frame[name]
        parts = [str(col.dtype).encode()]
        if isinstance(col.dtype, pd.CategoricalDtype):
            parts.append("\n".join(map(str, col.cat.categories)).encode())
            parts.append(np.asarray(col.cat.codes, dtype=np.int64).tobytes())
        else:
            hashed = pd.util.hash_pandas_object(col, index=False)
            parts.append(hashed.to_numpy().tobytes())
        out[str(name)] = h(*parts)
    return out


def matrix_digest(m: Any) -> str:
    m = sp.csr_matrix(m)
    m.sort_indices()
    return h(
        str(m.dtype).encode(),
        str(m.shape).encode(),
        np.asarray(m.indptr, dtype=np.int64).tobytes(),
        np.asarray(m.indices, dtype=np.int64).tobytes(),
        m.data.tobytes(),
    )


def cases(
    queries: list[dict[str, Any]], catalog: pd.DataFrame
) -> list[tuple[str, Any]]:
    out = []
    for q in queries:
        if q["shape"] == "single":
            out.append((q["query"], {"dataset_id": q["dataset_id"]}))
        elif q["shape"] == "cell_type":
            out.append((q["query"], {"cell_type": q["cell_type"]}))
        elif q["shape"] == "filter":
            out.append((q["query"], q["filter"]))
    if "consolidated_store_uri" in catalog.columns:
        rows = catalog[catalog["consolidated_store_uri"].notna()]
        for uri, group in rows.groupby("consolidated_store_uri"):
            name = "members:" + str(uri).rsplit("/", 1)[1]
            out.append((name, {"dataset_id": sorted(group["dataset_id"])}))
    return out


def main() -> None:
    catalog, queries_path, out_path = sys.argv[1:4]
    per_dataset = "--per-dataset" in sys.argv
    source = catalog
    if "--cases-from" in sys.argv:
        source = sys.argv[sys.argv.index("--cases-from") + 1]
    census = (
        census2_sdk.open_census([catalog], consolidated=False)
        if per_dataset
        else census2_sdk.open_census([catalog])
    )
    from census2_sdk import catalog as _catalog

    queries = json.loads(pathlib.Path(queries_path).read_text())
    results = {}
    for name, obs_filter in cases(queries, _catalog.load(source)):
        t = time.time()
        a = census.get_anndata(
            "homo_sapiens", obs_filter, layers=["normalized"], max_workers=16
        )
        if a is None:
            results[name] = None
        else:
            results[name] = {
                "n_obs": int(a.n_obs),
                "obs": frame_digest(pd.DataFrame(a.obs)),
                "var": frame_digest(pd.DataFrame(a.var)),
                "X": matrix_digest(a.X),
                "normalized": matrix_digest(a.layers["normalized"]),
                "stores_opened": a.uns["census2"]["stores_opened"],
            }
        print(
            name,
            results[name] and results[name]["n_obs"],
            round(time.time() - t, 1),
            flush=True,
        )
        del a
    pathlib.Path(out_path).write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
