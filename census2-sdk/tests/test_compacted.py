"""Tests for queries that read compacted stores. No network, no AWS.

The corpus is the three stores of `test_sdk`. `ds-a` and `ds-b` also have a second
copy in one compacted store, which `census2_dataset_writer.compaction` writes. The
catalog builder then names that store in their rows. A query must give the same
result from the compacted store as from the per-dataset stores.
"""

import shutil
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import census2_dataset_writer
import numpy as np
import pandas as pd
import pytest
import zarr
from census2_dataset_writer import compaction
from test_sdk import STORES, X, expected, make_store, same, write_one

import census2_sdk
from census2_sdk import catalog

MEMBERS = ("ds-a", "ds-b")


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("compacted")
    prefix, compacted = root / "out", root / "compacted"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_one(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
    sources = [prefix / f"v-{d}.zarr" for d in MEMBERS]
    compaction.compact(sources, compacted / "c-ab.zarr")
    build_catalog.build(
        [str(prefix)], root / "catalog.parquet", compacted=[str(compacted)]
    )
    return {
        "root": root,
        "prefix": prefix,
        "compacted": compacted,
        "catalog": root / "catalog.parquet",
    }


FILTERS: list[dict[str, Any]] = [
    {},
    {"tissue": "lung"},
    {"disease": "flu"},
    {"tissue": ["lung", "kidney"], "disease": "normal"},
    {"dataset_id": "ds-b"},
    {"dataset_id": ["ds-b", "ds-c"], "tissue": "lung"},
]


def test_catalog_names_the_compacted_store(corpus: dict[str, Path]) -> None:
    frame = catalog.load(corpus["catalog"]).set_index("dataset_id")
    uri = str(corpus["compacted"] / "c-ab.zarr")
    assert list(frame.loc[list(MEMBERS), "compacted_store_uri"]) == [uri, uri]
    assert frame.loc["ds-a", "compacted_row_start"] == 0
    assert frame.loc["ds-a", "compacted_row_stop"] == STORES["ds-a"][0]
    assert frame.loc["ds-b", "compacted_row_start"] == STORES["ds-a"][0]
    assert pd.isna(frame.loc["ds-c", "compacted_store_uri"])
    # The per-dataset store stays the store of each row.
    assert frame.loc["ds-a", "store_uri"] == str(corpus["prefix"] / "v-ds-a.zarr")


@pytest.mark.parametrize("obs_filter", FILTERS)
@pytest.mark.parametrize("max_workers", [1, 4])
def test_compacted_reads_equal_per_dataset_reads(
    corpus: dict[str, Path], obs_filter: dict[str, Any], max_workers: int
) -> None:
    fast = census2_sdk.open_census([corpus["catalog"]])
    slow = census2_sdk.open_census([corpus["catalog"]], compacted=False)
    got = fast.get_anndata("homo_sapiens", obs_filter, max_workers=max_workers)
    want = slow.get_anndata("homo_sapiens", obs_filter, max_workers=max_workers)
    assert got is not None and want is not None
    assert list(got.obs_names) == list(want.obs_names)
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", got.obs), cast("pd.DataFrame", want.obs)
    )
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", got.var), cast("pd.DataFrame", want.var)
    )
    assert X(got).dtype == X(want).dtype
    assert np.array_equal(X(got).toarray(), X(want).toarray())
    same(got, expected(corpus, obs_filter))
    # The per-dataset path opens one store for each dataset.
    assert want.uns["census2"]["stores_opened"] == want.uns["census2"]["datasets_read"]


def test_one_open_serves_every_member(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens")
    assert result is not None
    assert result.uns["census2"]["datasets_read"] == 3
    assert result.uns["census2"]["stores_opened"] == 2
    lung = census.get_anndata("homo_sapiens", {"tissue": "lung"})
    assert lung is not None
    assert lung.uns["census2"]["stores_opened"] == 1


def test_a_member_keeps_the_categories_of_its_store(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", {"dataset_id": "ds-c"})
    assert result is not None
    one = census.get_anndata("homo_sapiens", {"dataset_id": "ds-a"})
    assert one is not None
    stored = anndata.read_zarr(corpus["prefix"] / "v-ds-a.zarr").obs["tissue"]
    assert list(one.obs["tissue"].cat.categories) == list(stored.cat.categories)


def test_count_cells_is_the_same(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    for obs_filter in FILTERS:
        result = census.get_anndata("homo_sapiens", obs_filter)
        assert census.count_cells("homo_sapiens", obs_filter) == result.n_obs  # type: ignore[union-attr]


def test_rows_that_the_record_does_not_list_are_refused(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    frame = pd.read_parquet(corpus["catalog"])
    frame.loc[frame["dataset_id"] == "ds-b", "compacted_row_start"] = 1
    frame.to_parquet(tmp_path / "catalog.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"])
    with pytest.raises(ValueError, match="does not list v-ds-b at rows"):
        census.get_anndata("homo_sapiens", {"dataset_id": "ds-b"})


def test_a_per_dataset_store_named_as_compacted_is_refused(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    frame = pd.read_parquet(corpus["catalog"])
    frame.loc[frame["dataset_id"] == "ds-c", "compacted_store_uri"] = str(
        corpus["prefix"] / "v-ds-c.zarr"
    )
    frame.loc[frame["dataset_id"] == "ds-c", "compacted_row_start"] = 0
    frame.loc[frame["dataset_id"] == "ds-c", "compacted_row_stop"] = 20
    frame.to_parquet(tmp_path / "catalog.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"])
    with pytest.raises(ValueError, match="no compaction record"):
        census.get_anndata("homo_sapiens", {"dataset_id": "ds-c"})


# --- the catalog builder checks -------------------------------------------------


def copy_corpus(corpus: dict[str, Path], tmp_path: Path) -> tuple[Path, Path]:
    prefix, compacted = tmp_path / "out", tmp_path / "compacted"
    shutil.copytree(corpus["prefix"], prefix)
    shutil.copytree(corpus["compacted"], compacted)
    return prefix, compacted


def test_builder_refuses_a_member_with_other_cells(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix, compacted = copy_corpus(corpus, tmp_path)
    # The per-dataset store of ds-b gets other cells: a store of another version
    # written at the same address.
    shutil.rmtree(prefix / "v-ds-b.zarr")
    adata = make_store("ds-b", *STORES["ds-b"], seed=1)
    adata.obs["observation_joinid"] = pd.array(
        [f"other-{i}" for i in range(adata.n_obs)]
    )
    write_one(prefix, adata, "ds-b")
    with pytest.raises(ValueError, match=r"v-ds-b differs .*\['obs_digest'\]"):
        build_catalog.build([str(prefix)], compacted=[str(compacted)])


def test_builder_refuses_a_member_without_a_store(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix, compacted = copy_corpus(corpus, tmp_path)
    shutil.rmtree(prefix / "v-ds-a.zarr")
    with pytest.raises(ValueError, match="v-ds-a has no per-dataset store"):
        build_catalog.build([str(prefix)], compacted=[str(compacted)])


def test_builder_refuses_a_dataset_in_two_compacted_stores(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix, compacted = copy_corpus(corpus, tmp_path)
    compaction.compact(
        [prefix / "v-ds-b.zarr", prefix / "v-ds-c.zarr"], compacted / "c-bc.zarr"
    )
    with pytest.raises(ValueError, match="v-ds-b is in two compacted stores"):
        build_catalog.build([str(prefix)], compacted=[str(compacted)])


def test_builder_refuses_a_store_that_differs_from_its_record(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix, compacted = copy_corpus(corpus, tmp_path)
    store = compacted / "c-ab.zarr"
    record = compaction.read_compaction(store)
    assert record is not None
    first = record["members"][0]
    first["nnz"] = first["nnz"] + 1
    root = zarr.open_group(store, mode="r+")
    root.update_attributes(
        {
            census2_dataset_writer.SUMMARY_NAMESPACE: {
                census2_dataset_writer.COMPACTION_KEY: record
            }
        }
    )
    with pytest.raises(ValueError, match=r"v-ds-a differs .*\['nnz'\]"):
        build_catalog.build([str(prefix)], compacted=[str(compacted)])


def test_builder_refuses_members_of_another_spatial_type(
    corpus: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The compactor refuses such a group, so the test changes the rule afterwards."""
    prefix, compacted = copy_corpus(corpus, tmp_path)
    summarize = build_catalog.summarize

    def spatial_a(group: zarr.Group) -> tuple[dict[str, Any], pd.DataFrame]:
        row, cells = summarize(group)
        row["spatial"] = bool(np.array_equal(row["tissue"], ["liver", "lung"]))
        return row, cells

    monkeypatch.setattr(build_catalog, "summarize", spatial_a)
    with pytest.raises(ValueError, match="the members differ in spatial"):
        build_catalog.build([str(prefix)], compacted=[str(compacted)])


def test_a_catalog_without_compacted_stores_has_null_columns(tmp_path: Path) -> None:
    prefix = tmp_path / "out"
    prefix.mkdir()
    write_one(
        prefix,
        make_store("ds-x", 4, ["ENSG00000000000"], ["lung"], ["normal"], 3),
        "ds-x",
    )
    build_catalog.build([str(prefix)], tmp_path / "catalog.parquet")
    frame = catalog.load(tmp_path / "catalog.parquet")
    assert frame[list(catalog.COMPACTED)].isna().all().all()
    result = census2_sdk.open_census([tmp_path / "catalog.parquet"]).get_anndata(
        "homo_sapiens"
    )
    assert result is not None and result.n_obs == 4
