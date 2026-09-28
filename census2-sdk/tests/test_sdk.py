"""Tests for the catalog and the query path. No network, no AWS.

The fixture writes three small stores in the form that the writer produces: the
counts in `X`, a shared 20-gene reference as `var`, `observation_joinid` in `obs`,
and no `dataset_id`. The receipt beside each store records `gene_reference_version`. Each store measured a
different subset of the reference genes, so the others are zero columns.
"""

import json
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import census2_sdk
from census2_sdk import catalog

REFERENCE_IDS = [f"ENSG{j:011d}" for j in range(20)]

STORES = {
    # dataset_id: (n_obs, measured gene ids, tissues, diseases)
    "ds-a": (30, REFERENCE_IDS[:12], ["lung", "liver"], ["normal"]),
    "ds-b": (25, REFERENCE_IDS[5:20], ["lung"], ["normal ||  flu", "flu"]),
    "ds-c": (20, REFERENCE_IDS[::2], ["kidney"], ["normal"]),
}


def reference_var(ids: list[str]) -> pd.DataFrame:
    """A `var` in the writer's form: the reference, indexed by `feature_id`."""
    return pd.DataFrame(
        {
            "feature_id": ids,
            "feature_name": [f"GENE{j}" for j in range(len(ids))],
            "feature_version": ["1"] * len(ids),
            "feature_length": [1000] * len(ids),
            "feature_type": ["gene"] * len(ids),
        },
        index=ids,
    )


def make_store(
    dataset_id: str,
    n_obs: int,
    measured: list[str],
    tissues: list[str],
    diseases: list[str],
    seed: int,
    *,
    var_ids: list[str] = REFERENCE_IDS,
) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    x = sp.random(
        n_obs,
        len(var_ids),
        density=0.4,
        format="csr",
        dtype="float32",
        random_state=seed,
    )
    x.data = np.floor(x.data * 10) + 1
    # The genes that the dataset did not measure are zero columns.
    kept: Any = x.multiply(np.isin(var_ids, measured)[None, :])
    x = sp.csr_matrix(kept)
    obs = pd.DataFrame(
        {
            "tissue": pd.Categorical(rng.choice(tissues, n_obs)),
            "disease": pd.Categorical(rng.choice(diseases, n_obs)),
            "assay_ontology_term_id": pd.Categorical(["EFO:0009922"] * n_obs),
            "is_primary_data": np.ones(n_obs, dtype=bool),
            "observation_joinid": pd.array([f"{dataset_id}-{i}" for i in range(n_obs)]),
        },
        index=[f"cell{i}" for i in range(n_obs)],
    )
    adata = anndata.AnnData(X=x, obs=obs, var=reference_var(var_ids))
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    return adata


def write_one(
    prefix: Path,
    adata: anndata.AnnData,
    dataset_id: str,
    version: str | None = "GENCODE v48",
) -> None:
    """One store and its receipt, as the ingestion pipeline writes them."""
    version_id = f"v-{dataset_id}"
    adata.write_zarr(prefix / f"{version_id}.zarr")
    receipt = {
        "admitted": True,
        "dataset_id": dataset_id,
        "dataset_version_id": version_id,
        "schema_version": "7.1.0",
    }
    if version is not None:
        receipt["gene_reference_version"] = version
    (prefix / f"{version_id}.json").write_text(json.dumps(receipt))


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("corpus")
    prefix = root / "out"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_one(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
    (prefix / "v-rejected.json").write_text(json.dumps({"admitted": False}))
    frame = build_catalog.build(str(prefix), root / "catalog.parquet")
    assert len(frame) == 3
    return {"root": root, "prefix": prefix, "catalog": root / "catalog.parquet"}


def one_store_catalog(
    tmp_path: Path,
    adata: anndata.AnnData,
    dataset_id: str,
    version: str | None = "GENCODE v48",
) -> Path:
    prefix = tmp_path / "out"
    prefix.mkdir()
    write_one(prefix, adata, dataset_id, version)
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    return tmp_path / "catalog.parquet"


def expected(
    corpus: dict[str, Path], obs_filter: dict[str, str | list[str]]
) -> anndata.AnnData:
    """The answer by brute force: read every store whole and filter in pandas."""
    parts = []
    for dataset_id in STORES:
        adata = anndata.read_zarr(corpus["prefix"] / f"v-{dataset_id}.zarr")
        mask = np.ones(adata.n_obs, dtype=bool)
        for column, value in obs_filter.items():
            wanted = {value} if isinstance(value, str) else set(value)
            if column == "dataset_id":
                mask &= dataset_id in wanted
                continue
            series = adata.obs[column].astype(str)
            if column == "disease":
                mask &= np.array(
                    [
                        bool(wanted & {p.strip() for p in s.split(" || ")})
                        for s in series
                    ]
                )
            else:
                mask &= series.isin(wanted).to_numpy()
        sub = adata[mask].copy()
        obs = cast("pd.DataFrame", sub.obs).copy()
        obs.insert(0, "dataset_id", dataset_id)
        obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
        var = cast("pd.DataFrame", sub.var)
        parts.append(anndata.AnnData(X=X(sub), obs=obs, var=var))
    return anndata.concat(parts, join="inner")


def X(adata: anndata.AnnData) -> sp.csr_matrix:
    return sp.csr_matrix(cast("Any", adata.X))


def same(left: anndata.AnnData, right: anndata.AnnData) -> None:
    assert set(left.obs_names) == set(right.obs_names)
    right = right[left.obs_names].copy()
    assert list(left.var_names) == list(right.var_names)
    assert np.array_equal(X(left).toarray(), X(right).toarray())
    for column in ("dataset_id", "observation_joinid", "tissue", "disease"):
        assert list(left.obs[column].astype(str)) == list(right.obs[column].astype(str))


def test_catalog_rows(corpus: dict[str, Path]) -> None:
    frame = catalog.load(corpus["catalog"])
    assert list(frame["dataset_id"]) == ["ds-a", "ds-b", "ds-c"]
    assert list(frame["n_obs"]) == [30, 25, 20]
    assert set(frame["organism"]) == {"homo_sapiens"}
    assert not frame["spatial"].any()
    assert frame["tissue_general"].isna().all()  # absent from the stores
    assert set(frame["gene_reference_version"]) == {"GENCODE v48"}


def test_exclusion_is_conservative(corpus: dict[str, Path]) -> None:
    frame = catalog.load(corpus["catalog"])
    kept = catalog.exclude(frame, "homo_sapiens", {"tissue": "lung"})
    assert list(kept["dataset_id"]) == ["ds-a", "ds-b"]
    kept = catalog.exclude(frame, "homo_sapiens", {"disease": "flu"})
    assert list(kept["dataset_id"]) == ["ds-b"]  # membership in "normal ||  flu"
    kept = catalog.exclude(frame, "homo_sapiens", {"tissue_general": "lung"})
    assert len(kept) == 3  # an absent column proves nothing
    kept = catalog.exclude(frame, "homo_sapiens", {"dataset_id": ["ds-c"]})
    assert list(kept["dataset_id"]) == ["ds-c"]
    with pytest.raises(ValueError, match="not a catalog column"):
        catalog.exclude(frame, "homo_sapiens", {"donor_id": "x"})


@pytest.mark.parametrize(
    "obs_filter",
    [
        {"tissue": "lung"},
        {"disease": "flu"},
        {"tissue": ["lung", "kidney"], "disease": "normal"},
        {"dataset_id": "ds-b", "tissue": "lung"},
    ],
)
def test_get_anndata_matches_brute_force(
    corpus: dict[str, Path], obs_filter: dict[str, str | list[str]]
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", obs_filter)
    assert list(result.var_names) == REFERENCE_IDS
    assert result.uns["census2"]["stores_in_catalog"] == 3
    same(result, expected(corpus, obs_filter))


def test_result_shape_and_provenance(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", {"tissue": "lung"})
    assert result.uns["census2"]["stores_opened"] == 2
    assert result.uns["census2"]["sources"] == [str(corpus["catalog"])]
    assert result.obs.index.is_unique
    assert X(result).dtype == np.float32
    # var is the stores' var, with its reference columns.
    assert list(result.var.columns) == list(reference_var(REFERENCE_IDS).columns)


def test_empty_result(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", {"tissue": "brain"})
    # No store is opened, so no store supplies a var.
    assert result.shape == (0, 0)
    assert result.uns["census2"]["stores_opened"] == 0
    assert list(result.obs.columns) == ["dataset_id", "observation_joinid"]


def test_over_exclusion_is_visible(corpus: dict[str, Path]) -> None:
    """The check the harness relies on: drop a true candidate, the result shrinks."""
    census = census2_sdk.open_census([corpus["catalog"]])
    full = census.get_anndata("homo_sapiens", {"tissue": "lung"})
    census.catalog = census.catalog[census.catalog["dataset_id"] != "ds-a"]
    mutated = census.get_anndata("homo_sapiens", {"tissue": "lung"})
    assert set(mutated.obs_names) < set(full.obs_names)


def test_dataset_in_two_sources(corpus: dict[str, Path]) -> None:
    with pytest.raises(census2_sdk.SourceConflict, match="more than one source"):
        census2_sdk.open_census([corpus["catalog"], corpus["catalog"]])


def test_spatial_mix_refused(corpus: dict[str, Path], tmp_path: Path) -> None:
    frame = catalog.load(corpus["catalog"])
    frame.loc[frame["dataset_id"] == "ds-a", "spatial"] = True
    frame.to_parquet(tmp_path / "catalog.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match="spatial"):
        census.get_anndata("homo_sapiens", {"tissue": "lung"})
    # Excluding the spatial store makes the rest compose.
    assert census.get_anndata("homo_sapiens", {"dataset_id": "ds-b"}).n_obs == 25


def test_store_with_another_var_refused(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    """A store whose genes are in another order cannot concatenate by position."""
    prefix = tmp_path / "out"
    prefix.mkdir()
    for seed, dataset_id in enumerate(["ds-a", "ds-x"]):
        var_ids = REFERENCE_IDS if dataset_id == "ds-a" else REFERENCE_IDS[::-1]
        adata = make_store(
            dataset_id, 5, REFERENCE_IDS, ["lung"], ["normal"], seed, var_ids=var_ids
        )
        write_one(prefix, adata, dataset_id)
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match="v-ds-x.zarr: var differs"):
        census.get_anndata("homo_sapiens")


def test_row_without_a_recorded_reference_refused(tmp_path: Path) -> None:
    """A receipt that names no reference gives a catalog row that the SDK refuses."""
    adata = make_store("ds-u", 5, REFERENCE_IDS, ["lung"], ["normal"], 7)
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-u", None)])
    with pytest.raises(census2_sdk.SourceConflict, match="no gene_reference_version"):
        census.get_anndata("homo_sapiens")


def test_absent_column_selects_no_cells(corpus: dict[str, Path]) -> None:
    """The catalog keeps a store whose obs lacks the column; the read then matches none."""
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", {"tissue_general": "lung"})
    assert result.shape == (0, 20)
    assert result.uns["census2"]["stores_opened"] == 3
    # The opened stores still give the result its obs schema.
    assert "tissue" in result.obs.columns and "dataset_id" in result.obs.columns


def test_filter_values_can_be_any_iterable(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    generator = (t for t in ["lung", "kidney"])
    result = census.get_anndata("homo_sapiens", {"tissue": generator})
    same(result, expected(corpus, {"tissue": ["lung", "kidney"]}))


def test_duplicate_cell_key_refused(corpus: dict[str, Path], tmp_path: Path) -> None:
    adata = make_store("ds-d", 4, REFERENCE_IDS[:3], ["lung"], ["normal"], 3)
    adata.obs["observation_joinid"] = pd.array(["k", "k", "l", "m"])
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-d")])
    with pytest.raises(ValueError, match="not unique"):
        census.get_anndata("homo_sapiens")


def test_null_is_not_a_value(corpus: dict[str, Path], tmp_path: Path) -> None:
    """A null obs value is absent from the catalog and matches no filter."""
    adata = make_store("ds-n", 6, REFERENCE_IDS[:3], ["lung"], ["normal"], 5)
    tissue = adata.obs["tissue"].cat.add_categories(["liver"])
    tissue.iloc[:2] = None
    adata.obs["tissue"] = tissue
    path = one_store_catalog(tmp_path, adata, "ds-n")
    assert list(catalog.load(path)["tissue"][0]) == ["lung"]
    census = census2_sdk.open_census([path])
    assert census.get_anndata("homo_sapiens", {"tissue": "nan"}).n_obs == 0
    assert census.get_anndata("homo_sapiens", {"tissue": "lung"}).n_obs == 4


def test_recorded_reference_must_match_the_pin(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    frame = catalog.load(corpus["catalog"])
    frame["gene_reference_version"] = "GENCODE v47"
    frame.to_parquet(tmp_path / "v47.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "v47.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match="v47.*pinned.*v48"):
        census.get_anndata("homo_sapiens", {"tissue": "lung"})
    frame["gene_reference_version"] = "GENCODE v48"
    frame.to_parquet(tmp_path / "v48.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "v48.parquet"])
    assert census.get_anndata("homo_sapiens", {"tissue": "lung"}).n_obs > 0
