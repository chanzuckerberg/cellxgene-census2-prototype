"""Unit tests for the store builder. No network, no AWS, no container.

These cover the properties the module docstring promises, because a caller
relies on them: the round trip, the slot policy, admission as a normal outcome,
the memory contract, provenance, and destination handling.
"""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path
from typing import Any

import anndata
import census2_dataset_writer
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


def written(
    tmp_path: Path, **kwargs: Any
) -> tuple[Path, census2_dataset_writer.Result]:
    """Write a dataset to h5ad, build a store from it, return both."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(**kwargs).write_h5ad(source)
    result = census2_dataset_writer.create(source, dest)
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
    assert set(store["uns"]) == set(census2_dataset_writer.DEFAULTS["uns_keys"])
    for slot in census2_dataset_writer.OPTIONAL_SLOTS:
        assert slot not in store or len(list(store[slot])) == 0


def test_uns_keys_is_configurable(tmp_path: Path) -> None:
    """The carried uns keys are policy, not a fixed set in the code."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)
    census2_dataset_writer.create(
        source, dest, config={"all_slots": False, "uns_keys": ["title"]}
    )

    store: Any = zarr.open_consolidated(dest, mode="r")
    assert set(store["uns"]) == {"title"}


def test_unwanted_uns_keys_are_never_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A spatial dataset holds a large image under uns; the policy drops it, so
    reading it would be paid for and thrown away."""
    source = tmp_path / "source.h5ad"
    adata = dataset()
    adata.uns["spatial"] = {"image": np.zeros((64, 64, 3), dtype="uint8")}
    adata.write_h5ad(source)

    seen: list[str] = []

    def spy(elem: Any, *args: Any, **kwargs: Any) -> Any:
        seen.append(getattr(elem, "name", "?"))
        return read_elem(elem, *args, **kwargs)

    # By name: read_elem is imported by the module, not exported from it.
    monkeypatch.setattr("census2_dataset_writer.read_elem", spy)
    with h5py.File(source) as f:
        built, _ = census2_dataset_writer.build(f)

    assert built is not None
    assert "spatial" not in built.uns
    # Reading the group reads everything under it, so the group itself must
    # never be read: only the keys the policy keeps.
    assert "/uns" not in seen, seen
    assert "/uns/schema_version" in seen, seen


def test_all_slots_carries_everything(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(extras=True).write_h5ad(source)
    result = census2_dataset_writer.create(source, dest, config={"all_slots": True})

    store: Any = zarr.open_consolidated(dest, mode="r")
    assert "obsm/X_pca" in result["slots"]
    assert "varm/gene_stuff" in result["slots"]
    assert "layers/counts" in result["slots"]
    # The full uns, so the key the Census set drops survives here.
    assert "title" in set(store["uns"])


def test_slots_report_what_harmonize_left(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The list describes the written object, not the one build assembled."""

    def rename_and_drop(
        adata: anndata.AnnData, config: Any = None, refdata: Any = None
    ) -> anndata.AnnData:
        assert adata.X is not None
        x: Any = adata.X  # the getter's union is wider than layers accepts
        adata.layers["harmonized"] = x
        del adata.raw
        return adata

    monkeypatch.setattr(census2_dataset_writer, "harmonize", rename_and_drop)
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)

    with h5py.File(source) as f:
        adata, result = census2_dataset_writer.build(f)

    assert adata is not None
    assert "layers/harmonized" in result["slots"]
    assert "raw" not in result["slots"]


# --- admission ------------------------------------------------------------


def test_not_admitted_writes_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A rejection is a normal outcome: no store, a reason, and no exception."""
    monkeypatch.setattr(
        census2_dataset_writer,
        "admit",
        lambda adata, config=None, refdata=None: "no thanks",
    )
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)

    result = census2_dataset_writer.create(source, dest)

    assert result == {
        "admitted": False,
        "reason": "no thanks",
        "slots": [],
        "dropped_obs_columns": [],
    }
    assert not dest.exists()


# --- the assay allowlist --------------------------------------------------

ON_LIST = "EFO:0009922"  # 10x 3' v3
OFF_LIST = "EFO:0030027"  # snmC-Seq2, a DNA methylation assay


def with_assay(adata: anndata.AnnData, *terms: str) -> anndata.AnnData:
    """Give each cell an assay term, cycling through `terms`."""
    n = adata.n_obs
    adata.obs["assay_ontology_term_id"] = pd.Categorical(
        [terms[i % len(terms)] for i in range(n)]
    )
    return adata


def test_off_list_assay_is_rejected(tmp_path: Path) -> None:
    """The criterion that keeps dense methylation matrices out of the workflow."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), OFF_LIST).write_h5ad(source)

    result = census2_dataset_writer.create(source, dest)

    assert result["admitted"] is False
    assert OFF_LIST in str(result["reason"])
    assert not dest.exists()


def test_one_on_list_cell_admits_the_dataset(tmp_path: Path) -> None:
    """Census 1 masks the off-list cells. This does not; it takes them all."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), ON_LIST, OFF_LIST, OFF_LIST).write_h5ad(source)

    result = census2_dataset_writer.create(source, dest)

    assert result["admitted"] is True
    assert zarr.open_consolidated(str(dest))["obs"] is not None


def test_missing_assay_column_is_admitted(tmp_path: Path) -> None:
    """An AnnData without the column is not a Discover dataset, so the
    criterion does not apply. Every other test here relies on this."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)
    assert "assay_ontology_term_id" not in dataset().obs

    assert census2_dataset_writer.create(source, dest)["admitted"] is True


def test_empty_allowlist_turns_the_criterion_off(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), OFF_LIST).write_h5ad(source)

    config = census2_dataset_writer.load_config()
    config["assay_allowlist"] = []

    assert (
        census2_dataset_writer.create(source, dest, config=config)["admitted"] is True
    )


# --- the obs projection ---------------------------------------------------

# Both names are taken from the corpus: a cell2location factor on a Visium
# dataset, and a clinical unit on a 10x one. The second holds three slashes.
SLASH_COLS = ["T/NK_cycling", "eGFR_(ml/min/1.73m2)"]


def with_obs_columns(adata: anndata.AnnData, *names: str) -> anndata.AnnData:
    for i, name in enumerate(names):
        adata.obs[name] = np.full(adata.n_obs, float(i), dtype="float32")
    return adata


def write_source(adata: anndata.AnnData, path: Path) -> None:
    """Write a fixture h5ad that holds a slash in an obs column name.

    anndata 0.12 added `disallow_forward_slash_in_h5ad`, on by default, so
    anndata now refuses to write such a name even though HDF5 permits it. The
    corpus datasets that carry one predate the setting. Turn it off to build
    the fixture; the reader never needs it.
    """
    with anndata.settings.override(disallow_forward_slash_in_h5ad=False):
        adata.write_h5ad(path)


def test_project_obs_removes_only_the_slash_columns() -> None:
    adata = with_obs_columns(dataset(), *SLASH_COLS)
    kept = list(census2_dataset_writer.project_obs(adata).obs.columns)
    assert kept == ["cell_type"]


def test_project_obs_keeps_a_clean_frame_untouched() -> None:
    adata = dataset()
    before = list(adata.obs.columns)
    assert list(census2_dataset_writer.project_obs(adata).obs.columns) == before


def test_harmonize_projects_obs() -> None:
    """harmonize is the seam; project_obs is what it currently does."""
    adata = with_obs_columns(dataset(), *SLASH_COLS)
    assert list(census2_dataset_writer.harmonize(adata).obs.columns) == ["cell_type"]


def test_a_slash_column_no_longer_fails_the_write(tmp_path: Path) -> None:
    """Zarr gives `/` its path meaning, so this write used to fail outright."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    write_source(with_obs_columns(dataset(), *SLASH_COLS), source)

    result = census2_dataset_writer.create(source, dest)

    assert result["admitted"] is True
    assert result["dropped_obs_columns"] == SLASH_COLS  # source order
    store = zarr.open_consolidated(str(dest))
    obs: Any = read_elem(store["obs"])
    assert list(obs.columns) == ["cell_type"]


def test_dropped_obs_columns_is_empty_when_nothing_is_dropped(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    with h5py.File(source) as f:
        _, result = census2_dataset_writer.build(f)
    assert result["dropped_obs_columns"] == []


def test_a_rejected_dataset_reports_no_dropped_columns(tmp_path: Path) -> None:
    """A rejection short-circuits before harmonize, so there is nothing to drop."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    write_source(with_assay(with_obs_columns(dataset(), *SLASH_COLS), OFF_LIST), source)

    result = census2_dataset_writer.create(source, dest)

    assert result["admitted"] is False
    assert result["dropped_obs_columns"] == []


def test_default_allowlist_holds_rna_and_spatial() -> None:
    """The default is Census 1's two lists, merged."""
    allowlist = set(census2_dataset_writer.DEFAULTS["assay_allowlist"])
    assert ON_LIST in allowlist
    assert OFF_LIST not in allowlist
    assert "EFO:0022857" in allowlist  # Visium, from ALLOWED_SPATIAL_ASSAYS
    assert allowlist == set(census2_dataset_writer.assays.DEFAULT_ALLOWLIST)


# --- the write thread bound -----------------------------------------------


def test_threads_bounds_the_write_and_restores_the_setting(tmp_path: Path) -> None:
    """dask counts the host's CPUs unless told otherwise, and a container with
    a CPU share sets no cgroup quota for it to read."""
    import dask

    before = dask.config.get("num_workers", None)
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)

    assert census2_dataset_writer.create(source, dest, threads=2)["admitted"] is True
    assert dask.config.get("num_workers", None) == before


# --- the memory contract --------------------------------------------------


def test_build_leaves_the_matrix_lazy(tmp_path: Path) -> None:
    """`build` must not materialize X. A large dataset does not fit in memory."""
    import dask.array as da

    source = tmp_path / "source.h5ad"
    dataset(400, 20).write_h5ad(source)

    with h5py.File(source) as f:
        adata, _ = census2_dataset_writer.build(f)
        assert adata is not None
        assert isinstance(adata.X, da.Array), type(adata.X)
        # Blocked by rows, whole rows per block.
        assert adata.X.chunks[1] == (20,)


def test_row_chunk_is_bounded_by_the_row_count(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset(40, 12).write_h5ad(source)
    with h5py.File(source) as f:
        assert census2_dataset_writer.row_chunk(f["X"], 40) == 40


# --- provenance -----------------------------------------------------------


def test_provenance_reports_the_code_that_wrote_a_store() -> None:
    p = census2_dataset_writer.provenance()
    assert p["recipe"] == census2_dataset_writer.RECIPE
    assert p["zarr_format"] == census2_dataset_writer.ZARR_FORMAT
    # Reported from the installed distributions.
    assert p["census2_dataset_writer_version"] == version("census2-dataset-writer")
    assert p["anndata_version"] == version("anndata")
    assert p["zarr_version"] == version("zarr")


def test_provenance_follows_the_slot_policy() -> None:
    assert census2_dataset_writer.provenance()["slot_policy"] == "census"
    assert (
        census2_dataset_writer.provenance({"all_slots": True})["slot_policy"] == "all"
    )


# --- destinations ---------------------------------------------------------


def test_is_remote_covers_every_scheme() -> None:
    for scheme in census2_dataset_writer.REMOTE_SCHEMES:
        assert census2_dataset_writer.is_remote(f"{scheme}bucket/key")
    assert not census2_dataset_writer.is_remote("/tmp/local.zarr")
    assert not census2_dataset_writer.is_remote(Path("/tmp/local.zarr"))
    assert not census2_dataset_writer.is_remote(None)


def test_as_store_passes_a_local_path_through(tmp_path: Path) -> None:
    assert census2_dataset_writer.as_store(tmp_path) is tmp_path


def test_as_store_adapts_a_remote_url() -> None:
    """A bare s3:// string would make zarr fall back to fsspec and fail."""
    pytest.importorskip("obstore")
    store = census2_dataset_writer.as_store("s3://bucket/prefix.zarr")
    assert isinstance(store, zarr.storage.ObjectStore)
