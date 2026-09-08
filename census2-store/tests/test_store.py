"""Unit tests for the store builder. No network, no AWS, no container.

These cover the properties the module docstring promises, because those are the
ones a caller relies on and the ones that regress quietly: the round trip, the
slot policy, admission as a normal outcome, and the memory contract.
"""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path
from typing import Any

import anndata
import census2_store
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr
from anndata.io import read_elem


def dataset(
    n_obs: int = 40,
    n_var: int = 12,
    *,
    dense: bool = False,
    raw: bool = True,
    extras: bool = False,
    seed: int = 0,
) -> anndata.AnnData:
    """A small AnnData shaped like a Discover dataset."""
    rng = np.random.default_rng(seed)
    x: Any = (
        rng.random((n_obs, n_var), dtype="float32")
        if dense
        else sp.random(
            n_obs, n_var, density=0.3, format="csr", dtype="float32", random_state=seed
        )
    )
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(rng.choice(["a", "b"], n_obs))},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"feature_name": [f"gene{j}" for j in range(n_var)]},
        index=[f"ENSG{j:05d}" for j in range(n_var)],
    )
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
        "title": "not a Census uns key",
    }
    if extras:
        adata.obsm["X_pca"] = rng.random((n_obs, 3), dtype="float32")
        adata.varm["gene_stuff"] = rng.random((n_var, 2), dtype="float32")
        adata.layers["counts"] = x
    if raw:
        adata.raw = adata
    return adata


def written(tmp_path: Path, **kwargs: Any) -> tuple[Path, census2_store.Result]:
    """Write a dataset to h5ad, build a store from it, return both."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(**kwargs).write_h5ad(source)
    result = census2_store.create(source, dest)
    return dest, result


def values(array: Any) -> np.ndarray[Any, Any]:
    """A whole array as values, so a comparison ignores the stored width.

    Typed as Any at both ends: a zarr group member and a read_elem result are
    each a wide union, so this is where the dynamic boundary is named.
    """
    return np.asarray(array)


# --- the round trip -------------------------------------------------------


@pytest.mark.parametrize("dense", [False, True])
def test_round_trip(tmp_path: Path, dense: bool) -> None:
    dest, result = written(tmp_path, dense=dense, raw=not dense)
    assert result["admitted"] is True

    source: Any = anndata.read_h5ad(tmp_path / "source.h5ad")
    store: Any = zarr.open_consolidated(dest, mode="r")

    # Compare through read_elem, the same surface a reader uses. Reaching into
    # the layout instead is misleading: obs/_index reads back as a group.
    obs: Any = read_elem(store["obs"])
    assert list(obs.index) == list(source.obs.index)
    assert list(obs.columns) == list(source.obs.columns)
    assert obs["cell_type"].dtype == source.obs["cell_type"].dtype

    x: Any = read_elem(store["X"])
    if dense:
        assert np.allclose(x, source.X)
    else:
        assert x.nnz == source.X.nnz
        assert np.allclose(x.toarray(), source.X.toarray())


def test_sparse_index_values_survive_the_round_trip(tmp_path: Path) -> None:
    """Compare the index arrays by value, not by stored width.

    The round trip widens them from int32 to int64, which is why the pipeline's
    own round-trip check hashes values rather than bytes.
    """
    dest, _ = written(tmp_path)
    source: Any = anndata.read_h5ad(tmp_path / "source.h5ad")
    store: Any = zarr.open_consolidated(dest, mode="r")

    for part in ("data", "indices", "indptr"):
        assert np.array_equal(values(store["X"][part]), values(getattr(source.X, part)))


def test_raw_is_carried(tmp_path: Path) -> None:
    dest, result = written(tmp_path, raw=True)
    store: Any = zarr.open_consolidated(dest, mode="r")
    assert "raw" in result["slots"]
    assert values(store["raw"]["X"]["data"]).size > 0


def test_consolidated_metadata_is_present(tmp_path: Path) -> None:
    """Written last, so opening this way is also the completeness check."""
    dest, _ = written(tmp_path)
    assert zarr.open_consolidated(dest, mode="r") is not None


# --- slot policy ----------------------------------------------------------


def test_census_slot_policy_drops_extras(tmp_path: Path) -> None:
    dest, result = written(tmp_path, extras=True)
    store: Any = zarr.open_consolidated(dest, mode="r")

    assert result["slots"] == ["obs", "var", "X", "raw", "uns[2 keys]"]
    assert set(store["uns"]) == set(census2_store.CENSUS_UNS_KEYS)
    for slot in census2_store.OPTIONAL_SLOTS:
        assert slot not in store or len(list(store[slot])) == 0


def test_all_slots_carries_everything(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(extras=True).write_h5ad(source)
    result = census2_store.create(source, dest, config={"all_slots": True})

    store: Any = zarr.open_consolidated(dest, mode="r")
    assert "obsm/X_pca" in result["slots"]
    assert "varm/gene_stuff" in result["slots"]
    assert "layers/counts" in result["slots"]
    # The full uns, so the key the Census set drops survives here.
    assert "title" in set(store["uns"])


# --- admission ------------------------------------------------------------


def test_not_admitted_writes_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A rejection is a normal outcome: no store, a reason, and no exception."""
    monkeypatch.setattr(
        census2_store, "admit", lambda adata, config=None, refdata=None: "no thanks"
    )
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)

    result = census2_store.create(source, dest)

    assert result == {"admitted": False, "reason": "no thanks", "slots": []}
    assert not dest.exists()


# --- the memory contract --------------------------------------------------


def test_build_leaves_the_matrix_lazy(tmp_path: Path) -> None:
    """`build` must not materialize X. A large dataset depends on this."""
    import dask.array as da

    source = tmp_path / "source.h5ad"
    dataset(400, 20).write_h5ad(source)

    with h5py.File(source) as f:
        adata, _ = census2_store.build(f)
        assert adata is not None
        assert isinstance(adata.X, da.Array), type(adata.X)
        # Blocked by rows, whole rows per block.
        assert adata.X.chunks[1] == (20,)


def test_row_chunk_is_bounded_by_the_row_count(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset(40, 12).write_h5ad(source)
    with h5py.File(source) as f:
        assert census2_store.row_chunk(f["X"], 40) == 40


# --- provenance -----------------------------------------------------------


def test_provenance_reports_the_code_that_wrote_a_store() -> None:
    p = census2_store.provenance()
    assert p["recipe"] == census2_store.RECIPE
    assert p["zarr_format"] == census2_store.ZARR_FORMAT
    # Read from the installed distributions, so a receipt cannot disagree with
    # the code that produced it.
    assert p["census2_store_version"] == version("census2-store")
    assert p["anndata_version"] == version("anndata")
    assert p["zarr_version"] == version("zarr")


def test_provenance_follows_the_slot_policy() -> None:
    assert census2_store.provenance()["slot_policy"] == "census"
    assert census2_store.provenance({"all_slots": True})["slot_policy"] == "all"


# --- destinations ---------------------------------------------------------


def test_is_remote_covers_every_scheme() -> None:
    for scheme in census2_store.REMOTE_SCHEMES:
        assert census2_store.is_remote(f"{scheme}bucket/key")
    assert not census2_store.is_remote("/tmp/local.zarr")
    assert not census2_store.is_remote(Path("/tmp/local.zarr"))
    assert not census2_store.is_remote(None)


def test_as_store_passes_a_local_path_through(tmp_path: Path) -> None:
    assert census2_store.as_store(tmp_path) is tmp_path


def test_as_store_adapts_a_remote_url() -> None:
    """A bare s3:// string would make zarr fall back to fsspec and fail."""
    pytest.importorskip("obstore")
    store = census2_store.as_store("s3://bucket/prefix.zarr")
    assert isinstance(store, zarr.storage.ObjectStore)
