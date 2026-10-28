"""Unit tests for compaction: one store from several per-dataset stores.

The per-dataset stores come from `create`, over the small datasets and the test
gene reference of `test_writer`. Their `obs_names` are the same in every dataset,
as author cell names often are, so a test can see that the compacted index is the
cell key.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import anndata
import census2_dataset_writer
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr
from anndata.io import read_elem
from census2_dataset_writer import compaction
from test_writer import dataset, refdata  # noqa: F401  # refdata is an autouse fixture

SIZES = (40, 25, 60)


def per_dataset_stores(tmp_path: Path, sizes: tuple[int, ...] = SIZES) -> list[Path]:
    """One store for each size, each with its own identifiers and cell types."""
    stores = []
    for k, n_obs in enumerate(sizes):
        adata = dataset(n_obs=n_obs, seed=k)
        # Each dataset has a cell type of its own, so the union has more categories
        # than any member.
        cell_type = adata.obs["cell_type"].astype(str).to_numpy(dtype=object)
        cell_type[0] = f"only-in-{k}"
        adata.obs["cell_type"] = pd.Categorical(cell_type)
        source = tmp_path / f"source-{k}.h5ad"
        adata.write_h5ad(source)
        dest = tmp_path / "stores" / f"version-{k}.zarr"
        census2_dataset_writer.create(
            source,
            dest,
            dataset_id=f"dataset-{k}",
            dataset_version_id=f"version-{k}",
            refdata=tmp_path / "refdata",
        )
        stores.append(dest)
    return stores


def compacted(
    tmp_path: Path,
) -> tuple[list[Path], Path, census2_dataset_writer.CompactionRecord]:
    stores = per_dataset_stores(tmp_path)
    dest = tmp_path / "compacted" / "c-test.zarr"
    record = compaction.compact(stores, dest)
    return stores, dest, record


def x_of(path: Path) -> sp.csr_matrix:
    x: Any = read_elem(zarr.open_group(path, mode="r")["X"])
    return sp.csr_matrix(x)


def obs_of(path: Path) -> pd.DataFrame:
    return cast("pd.DataFrame", read_elem(zarr.open_group(path, mode="r")["obs"]))


def test_members_follow_one_another_in_source_order(tmp_path: Path) -> None:
    _, dest, record = compacted(tmp_path)
    starts = np.cumsum((0, *SIZES))
    assert [(m["row_start"], m["row_stop"]) for m in record["members"]] == list(
        zip(starts[:-1].tolist(), starts[1:].tolist(), strict=True)
    )
    assert [m["dataset_version_id"] for m in record["members"]] == [
        f"version-{k}" for k in range(len(SIZES))
    ]
    assert compaction.read_compaction(dest) == record
    # A compacted store has no summary record, so a reader that wants one dataset
    # refuses it.
    assert census2_dataset_writer.read_summary(dest) is None


def test_a_member_holds_the_cells_and_counts_of_its_store(tmp_path: Path) -> None:
    stores, dest, record = compacted(tmp_path)
    x, obs = x_of(dest), obs_of(dest)
    for member, store in zip(record["members"], stores, strict=True):
        rows = slice(member["row_start"], member["row_stop"])
        assert np.array_equal(x[rows].toarray(), x_of(store).toarray())
        assert member["nnz"] == x_of(store).nnz
        want = obs_of(store)
        got = obs.iloc[rows]
        for column in want.columns:
            part = got[column]
            if isinstance(part.dtype, pd.CategoricalDtype):
                # Without the unused categories, a member's column is its store's.
                part = part.cat.remove_unused_categories()
            pd.testing.assert_series_equal(
                part.reset_index(drop=True), want[column].reset_index(drop=True)
            )
        assert member["obs_digest"] == census2_dataset_writer.obs_digest(
            want["observation_joinid"].astype(str)
        )
    assert compaction.verify(dest, stores) == []


def test_categories_are_the_sorted_union(tmp_path: Path) -> None:
    stores, dest, _ = compacted(tmp_path)
    obs = obs_of(dest)
    union = sorted(
        {str(c) for s in stores for c in obs_of(s)["cell_type"].cat.categories}
    )
    assert list(obs["cell_type"].cat.categories) == union
    assert list(obs["dataset_id"].cat.categories) == [
        f"dataset-{k}" for k in range(len(SIZES))
    ]
    assert not obs["cell_type"].cat.ordered


def test_the_index_is_the_cell_key(tmp_path: Path) -> None:
    stores, dest, _ = compacted(tmp_path)
    # The author cell names collide across the members.
    assert obs_of(stores[0]).index[0] == obs_of(stores[1]).index[0]
    obs = obs_of(dest)
    assert obs.index.is_unique
    want = obs["dataset_id"].astype(str) + "_" + obs["observation_joinid"].astype(str)
    assert list(obs.index) == list(want)


def test_a_plain_reader_opens_the_store(tmp_path: Path) -> None:
    stores, dest, _ = compacted(tmp_path)
    adata = anndata.read_zarr(dest)
    assert adata.n_obs == sum(SIZES)
    var = cast("pd.DataFrame", anndata.read_zarr(stores[0]).var)
    assert cast("pd.DataFrame", adata.var).equals(var)
    assert adata.uns["organism_ontology_term_id"] == "NCBITaxon:9606"
    # The consolidated metadata survives the record, which is the last write.
    assert zarr.open_consolidated(dest, mode="r") is not None


def test_obs_and_x_have_the_fixed_layout(tmp_path: Path) -> None:
    _, dest, _ = compacted(tmp_path)
    root = zarr.open_group(dest, mode="r")
    codes = root["obs/cell_type/codes"]
    assert isinstance(codes, zarr.Array)
    assert codes.chunks == (compaction.OBS_CHUNK_ROWS,)
    assert codes.shards == (compaction.OBS_SHARD_ROWS,)
    data = root["X/data"]
    assert isinstance(data, zarr.Array)
    assert data.chunks == (compaction.X_CHUNK_VALUES,)


def test_x_copies_across_shard_boundaries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Small shards, so the copy writes many whole shards and one partial shard."""
    monkeypatch.setattr(compaction, "X_CHUNK_VALUES", 8)
    monkeypatch.setattr(compaction, "X_SHARD_VALUES", 32)
    stores, dest, _ = compacted(tmp_path)
    assert compaction.verify(dest, stores) == []
    x = x_of(dest)
    want = sp.vstack([x_of(s) for s in stores])
    assert np.array_equal(x.toarray(), want.toarray())


def test_verify_finds_a_changed_value(tmp_path: Path) -> None:
    stores, dest, record = compacted(tmp_path)
    data = zarr.open_group(dest, mode="r+")["X/data"]
    assert isinstance(data, zarr.Array)
    data[0:1] = np.asarray(data[0:1]) + 1
    failures = compaction.verify(dest, stores)
    assert failures == [
        f"{record['members'][0]['dataset_version_id']}: X data differs at 0"
    ]


def test_a_store_without_a_summary_record_is_refused(tmp_path: Path) -> None:
    stores = per_dataset_stores(tmp_path, (10, 12))
    root = zarr.open_group(stores[1], mode="r+")
    root.update_attributes({census2_dataset_writer.SUMMARY_NAMESPACE: {}})
    with pytest.raises(compaction.CompactionError, match="no summary record"):
        compaction.compact(stores, tmp_path / "out.zarr")


def test_a_dataset_version_twice_is_refused(tmp_path: Path) -> None:
    stores = per_dataset_stores(tmp_path, (10,))
    with pytest.raises(compaction.CompactionError, match="occurs twice"):
        compaction.compact([stores[0], stores[0]], tmp_path / "out.zarr")


def test_another_gene_reference_is_refused(tmp_path: Path) -> None:
    stores = per_dataset_stores(tmp_path, (10, 12))
    record = dict(census2_dataset_writer.read_summary(stores[1]) or {})
    record["gene_reference_version"] = "TEST v2"
    root = zarr.open_group(stores[1], mode="r+")
    root.update_attributes(
        {
            census2_dataset_writer.SUMMARY_NAMESPACE: {
                census2_dataset_writer.SUMMARY_KEY: record
            }
        }
    )
    with pytest.raises(compaction.CompactionError, match="gene_reference_version"):
        compaction.compact(stores, tmp_path / "out.zarr")


def member(start: int, stop: int, name: str) -> dict[str, Any]:
    return {
        "dataset_id": f"d-{name}",
        "dataset_version_id": f"v-{name}",
        "row_start": start,
        "row_stop": stop,
        "nnz": 1,
        "obs_digest": "0" * 64,
        "uns": {},
    }


@pytest.mark.parametrize(
    ("members", "message"),
    [
        ([member(1, 4, "a")], "the next row is 0"),
        ([member(0, 4, "a"), member(5, 6, "b")], "the next row is 4"),
        ([member(0, 4, "a"), member(4, 4, "b")], "has rows 4:4"),
        ([member(0, 4, "a"), member(4, 6, "a")], "occurs twice"),
        ([], "one or more members"),
    ],
)
def test_validate_compaction_refuses_bad_rows(
    members: list[dict[str, Any]], message: str
) -> None:
    record = {
        "version": census2_dataset_writer.COMPACTION_VERSION,
        "gene_reference_version": "TEST v1",
        "members": members,
    }
    with pytest.raises(census2_dataset_writer.SummaryRecordError, match=message):
        census2_dataset_writer.validate_compaction(record)


def test_obs_digest_depends_on_the_order() -> None:
    digest = census2_dataset_writer.obs_digest
    assert digest(["a", "b"]) != digest(["b", "a"])
    # The separator keeps the boundaries, so a join of the values is not enough.
    assert digest(["ab", "c"]) != digest(["a", "bc"])
