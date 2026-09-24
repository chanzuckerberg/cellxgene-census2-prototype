"""Tests for the catalog and the query path. No network, no AWS.

The fixture writes three small stores in the Iteration 1 form: source `var` with a
subset of a 20-gene reference in a store-specific order, `observation_joinid` in
`obs`, no `dataset_id`, and one store with a `raw` group. The reference itself is a
20-row file in a `refdata` directory, so nothing is downloaded.
"""

import gzip
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
from census2_sdk.gene_reference import ReferencePinMismatch

REFERENCE_IDS = [f"ENSG{j:011d}" for j in range(20)]

STORES = {
    # dataset_id: (n_obs, gene ids in store order, tissues, diseases, raw)
    "ds-a": (30, REFERENCE_IDS[:12], ["lung", "liver"], ["normal"], False),
    "ds-b": (25, REFERENCE_IDS[5:20][::-1], ["lung"], ["normal ||  flu", "flu"], True),
    "ds-c": (20, REFERENCE_IDS[::2], ["kidney"], ["normal"], False),
}


def write_reference(refdata: Path) -> None:
    refdata.mkdir()
    rows = [f"{g},GENE{j},1,1000,gene" for j, g in enumerate(REFERENCE_IDS)]
    (refdata / "genes_homo_sapiens.csv.gz").write_bytes(
        gzip.compress(("\n".join(rows) + "\n").encode())
    )


def make_store(
    dataset_id: str,
    n_obs: int,
    genes: list[str],
    tissues: list[str],
    diseases: list[str],
    raw: bool,
    seed: int,
) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    x = sp.random(
        n_obs, len(genes), density=0.4, format="csr", dtype="float32", random_state=seed
    )
    x.data = np.floor(x.data * 10) + 1
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
    var = pd.DataFrame(
        {
            "feature_biotype": pd.Categorical(["gene"] * len(genes)),
            "feature_reference": pd.Categorical(["NCBITaxon:9606"] * len(genes)),
        },
        index=genes,
    )
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    if raw:
        # X is a transformed value; raw/X holds the counts and one gene more.
        adata.raw = anndata.AnnData(
            X=sp.csr_matrix(
                np.hstack([x.toarray(), np.zeros((n_obs, 1), dtype="float32")])
            ),
            var=pd.DataFrame(index=[*genes, REFERENCE_IDS[0]]),
        )
        adata.X = x * 0.5
    return adata


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("corpus")
    prefix = root / "out"
    prefix.mkdir()
    write_reference(root / "refdata")
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        version_id = f"v-{dataset_id}"
        adata = make_store(dataset_id, *spec, seed=seed)
        adata.write_zarr(prefix / f"{version_id}.zarr")
        (prefix / f"{version_id}.json").write_text(
            json.dumps(
                {
                    "admitted": True,
                    "dataset_id": dataset_id,
                    "dataset_version_id": version_id,
                    "schema_version": "7.1.0",
                }
            )
        )
    (prefix / "v-rejected.json").write_text(json.dumps({"admitted": False}))
    frame = build_catalog.build(str(prefix), root / "catalog.parquet")
    assert len(frame) == 3
    return {"root": root, "prefix": prefix, "catalog": root / "catalog.parquet"}


def expected(
    corpus: dict[str, Path], obs_filter: dict[str, str | list[str]]
) -> anndata.AnnData:
    """The answer by brute force: read every store whole, filter in pandas, align genes."""
    parts = []
    for dataset_id, (_, _, _, _, raw) in STORES.items():
        adata = anndata.read_zarr(corpus["prefix"] / f"v-{dataset_id}.zarr")
        counts = adata.raw.to_adata() if raw else adata
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
        sub = counts[mask].copy()
        dense = pd.DataFrame(X(sub).toarray(), columns=sub.var_names).reindex(
            columns=REFERENCE_IDS, fill_value=0
        )
        obs = cast("pd.DataFrame", sub.obs).copy()
        obs.insert(0, "dataset_id", dataset_id)
        obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
        parts.append(
            anndata.AnnData(
                X=sp.csr_matrix(dense.to_numpy(dtype=np.float32)),
                obs=obs,
                var=pd.DataFrame(index=REFERENCE_IDS),
            )
        )
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
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
    result = census.get_anndata("homo_sapiens", obs_filter)
    assert list(result.var_names) == REFERENCE_IDS
    assert result.uns["census2"]["stores_in_catalog"] == 3
    same(result, expected(corpus, obs_filter))


def test_result_shape_and_provenance(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
    result = census.get_anndata("homo_sapiens", {"tissue": "lung"})
    assert result.uns["census2"]["stores_opened"] == 2
    assert result.uns["census2"]["sources"] == [str(corpus["catalog"])]
    assert result.obs.index.is_unique
    assert X(result).dtype == np.float32
    # ds-b counts come from raw/X, not from the halved X.
    b = result[result.obs["dataset_id"] == "ds-b"]
    assert (X(b).data % 1 == 0).all()


def test_empty_result(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
    result = census.get_anndata("homo_sapiens", {"tissue": "brain"})
    assert result.shape == (0, 20)
    assert result.uns["census2"]["stores_opened"] == 0
    assert list(result.obs.columns) == ["dataset_id", "observation_joinid"]


def test_over_exclusion_is_visible(corpus: dict[str, Path]) -> None:
    """The check the harness relies on: drop a true candidate, the result shrinks."""
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
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
    census = census2_sdk.open_census(
        [tmp_path / "catalog.parquet"], refdata=corpus["root"] / "refdata"
    )
    with pytest.raises(census2_sdk.SourceConflict, match="spatial"):
        census.get_anndata("homo_sapiens", {"tissue": "lung"})
    # Excluding the spatial store makes the rest compose.
    assert census.get_anndata("homo_sapiens", {"dataset_id": "ds-b"}).n_obs == 25


def test_gene_outside_reference_refused(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix = tmp_path / "out"
    prefix.mkdir()
    adata = make_store(
        "ds-x",
        5,
        ["ENSG99999999999", *REFERENCE_IDS[:3]],
        ["lung"],
        ["normal"],
        False,
        9,
    )
    adata.write_zarr(prefix / "v-x.zarr")
    (prefix / "v-x.json").write_text(
        json.dumps(
            {"admitted": True, "dataset_id": "ds-x", "dataset_version_id": "v-x"}
        )
    )
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    census = census2_sdk.open_census(
        [tmp_path / "catalog.parquet"], refdata=corpus["root"] / "refdata"
    )
    with pytest.raises(ReferencePinMismatch, match="1 of 4"):
        census.get_anndata("homo_sapiens")


def test_absent_column_selects_no_cells(corpus: dict[str, Path]) -> None:
    """The catalog keeps a store whose obs lacks the column; the read then matches none."""
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
    result = census.get_anndata("homo_sapiens", {"tissue_general": "lung"})
    assert result.shape == (0, 20)
    assert result.uns["census2"]["stores_opened"] == 3
    # The opened stores still give the result its obs schema.
    assert "tissue" in result.obs.columns and "dataset_id" in result.obs.columns


def test_filter_values_can_be_any_iterable(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census(
        [corpus["catalog"]], refdata=corpus["root"] / "refdata"
    )
    generator = (t for t in ["lung", "kidney"])
    result = census.get_anndata("homo_sapiens", {"tissue": generator})
    same(result, expected(corpus, {"tissue": ["lung", "kidney"]}))


def test_duplicate_cell_key_refused(corpus: dict[str, Path], tmp_path: Path) -> None:
    prefix = tmp_path / "out"
    prefix.mkdir()
    adata = make_store("ds-d", 4, REFERENCE_IDS[:3], ["lung"], ["normal"], False, 3)
    adata.obs["observation_joinid"] = pd.array(["k", "k", "l", "m"])
    adata.write_zarr(prefix / "v-d.zarr")
    (prefix / "v-d.json").write_text(
        json.dumps(
            {"admitted": True, "dataset_id": "ds-d", "dataset_version_id": "v-d"}
        )
    )
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    census = census2_sdk.open_census(
        [tmp_path / "catalog.parquet"], refdata=corpus["root"] / "refdata"
    )
    with pytest.raises(ValueError, match="not unique"):
        census.get_anndata("homo_sapiens")
