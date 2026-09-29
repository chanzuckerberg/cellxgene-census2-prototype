"""The canonical result form (task 2.3): a deterministic serialization of a query result.

All normalization happens here, in the serialization, and none in the comparator. The
rules:

- `obs` rows are ordered by identity, `(dataset_id, observation_joinid)`, and `var`
  rows by `feature_id`. Position plays no part.
- A categorical column is decoded to its values. The values keep their dtype when it
  has a NumPy form and no value is null. In every other case they become strings.
  Codes do not leave the stack.
- Each X layer is reordered to the sorted obs x sorted var coordinate space, as CSR
  with duplicate coordinates summed, explicit zeros removed, and sorted indices, dtype
  preserved. Raw counts are exact.
- Every array is written as `.npy`, whose bytes are a function of shape, dtype, and
  content only. A string column becomes two arrays, uint64 offsets and UTF-8 bytes.
  Zip, gzip, and Parquet embed timestamps or writer metadata, so they are not used.
- The manifest is JSON with sorted keys: one SHA-256 per file, and a `result_digest`
  over the sorted file digests. Two runs agree when their result digests agree.
"""

import hashlib
import io
import json
import pathlib
from collections.abc import Iterable
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

Arrays = dict[str, np.ndarray]


def _npy_bytes(arr: np.ndarray) -> bytes:
    buf = io.BytesIO()
    np.save(buf, arr, allow_pickle=False)
    return buf.getvalue()


def _string_arrays(values: Iterable[Any]) -> Arrays:
    encoded = [str(v).encode("utf-8") for v in values]
    offsets = np.zeros(len(encoded) + 1, dtype=np.uint64)
    np.cumsum([len(b) for b in encoded], out=offsets[1:])
    return {
        "offsets": offsets,
        "bytes": np.frombuffer(b"".join(encoded), dtype=np.uint8),
    }


def _column_arrays(name: str, col: pd.Series) -> Arrays:
    if isinstance(col.dtype, pd.CategoricalDtype):
        # decode to the categories' dtype, so a numeric categorical stays numeric
        cats = col.cat.categories.dtype
        col = col.astype(object if col.isna().any() else cats)
    dtype = str(col.dtype)
    if col.dtype == object or dtype in ("string", "str"):
        return _string_column(name, col)
    if dtype.startswith(
        ("Int", "UInt", "Float", "boolean")
    ):  # nullable extension dtypes
        assert not col.isna().any(), (
            f"column {name} has nulls; a null policy must be decided"
        )
        col = col.astype(cast("Any", dtype.lower()))
    arr = np.ascontiguousarray(col.to_numpy())
    if (
        arr.dtype == object
    ):  # an extension dtype with no NumPy form, such as tz datetimes
        return _string_column(name, col)
    return {name: arr}


def _string_column(name: str, col: pd.Series) -> Arrays:
    parts = _string_arrays(col.astype(str))
    return {f"{name}.offsets": parts["offsets"], f"{name}.bytes": parts["bytes"]}


def frame_arrays(kind: str, frame: pd.DataFrame) -> Arrays:
    out: Arrays = {}
    for name in sorted(frame.columns):
        out.update(_column_arrays(f"{kind}.{name}", frame[name]))
    return out


def layer_arrays(name: str, matrix: sparse.csr_matrix) -> Arrays:
    m = matrix.tocsr().copy()
    m.sum_duplicates()
    m.eliminate_zeros()
    m.sort_indices()
    return {
        f"X.{name}.indptr": m.indptr,
        f"X.{name}.indices": m.indices,
        f"X.{name}.data": m.data,
        f"X.{name}.shape": np.asarray(m.shape, dtype=np.int64),
    }


def decode_strings(offsets: np.ndarray, data: np.ndarray) -> list[str]:
    """The inverse of the string encoding, to name identities in a comparison report."""
    raw = data.tobytes()
    return [
        raw[offsets[i] : offsets[i + 1]].decode("utf-8")
        for i in range(len(offsets) - 1)
    ]


def load_result(result_dir: pathlib.Path) -> Arrays:
    """Load a serialized result back into its arrays.

    Only the files the manifest lists are loaded, so a rewrite that dropped a column
    does not read the stale `.npy` left behind.
    """
    manifest = result_dir / "manifest.json"
    if manifest.exists():
        names = sorted(json.loads(manifest.read_text())["files"])
    else:
        names = sorted(p.name.removesuffix(".npy") for p in result_dir.glob("*.npy"))
    return {n: np.load(result_dir / f"{n}.npy", allow_pickle=False) for n in names}


def result_digest(arrays: Arrays) -> str:
    """The digest that `write_result` records, computed without a write."""
    digests = {
        name: hashlib.sha256(_npy_bytes(arr)).hexdigest()
        for name, arr in sorted(arrays.items())
    }
    return hashlib.sha256(json.dumps(digests, sort_keys=True).encode()).hexdigest()


def write_result(out_dir: pathlib.Path, arrays: Arrays, meta: dict[str, Any]) -> str:
    """Write the arrays and the manifest. Return the result digest."""
    out_dir.mkdir(parents=True, exist_ok=True)
    digests = {}
    for name, arr in sorted(arrays.items()):
        blob = _npy_bytes(arr)
        (out_dir / f"{name}.npy").write_bytes(blob)
        digests[name] = hashlib.sha256(blob).hexdigest()
    digest = hashlib.sha256(json.dumps(digests, sort_keys=True).encode()).hexdigest()
    manifest = {"files": digests, "result_digest": digest, **meta}
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=1, sort_keys=True) + "\n"
    )
    return digest


def canonicalize_anndata(adata: anndata.AnnData) -> Arrays:
    """The canonical arrays of an AnnData that holds `obs`, `var`, and layers.

    `obs` must hold `dataset_id` and `observation_joinid`, and `var` must hold
    `feature_id`. `uns["x_layer_name"]` names the layer in `X`; the default is `raw`.
    """
    obs = cast("pd.DataFrame", adata.obs).reset_index(drop=True)
    var = cast("pd.DataFrame", adata.var).reset_index(drop=True)
    assert {"dataset_id", "observation_joinid"} <= set(obs.columns), (
        "identity columns missing from obs"
    )
    assert "feature_id" in var.columns, "feature_id missing from var"

    obs_order = np.lexsort(
        (
            obs["observation_joinid"].astype(str).to_numpy(),
            obs["dataset_id"].astype(str).to_numpy(),
        )
    )
    var_order = np.argsort(var["feature_id"].astype(str).to_numpy(), kind="stable")
    obs = obs.iloc[obs_order].reset_index(drop=True)
    var = var.iloc[var_order].reset_index(drop=True)

    arrays = frame_arrays("obs", obs)
    arrays.update(frame_arrays("var", var))
    first_layer = str(adata.uns.get("x_layer_name", "raw"))
    if adata.X is not None:
        x = sparse.csr_matrix(cast("Any", adata.X))[obs_order][:, var_order]
        arrays.update(layer_arrays(first_layer, x))
    for lname, layer in adata.layers.items():
        # anndata 0.13 exposes X itself under the key None; it is serialized above.
        if lname is None or (lname == first_layer and adata.X is not None):
            continue
        arrays.update(
            layer_arrays(lname, sparse.csr_matrix(layer)[obs_order][:, var_order])
        )
    return arrays
