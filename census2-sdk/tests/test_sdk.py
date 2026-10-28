"""Tests for the catalog and the query path. No network, no AWS.

The fixture writes three small stores in the form that the writer produces: the
counts in `X`, a shared 20-gene reference as `var`, `observation_joinid` in `obs`,
and `dataset_id` first in `obs`. Each store's summary record, in the attributes of its root group, holds the
dataset's identifiers and `gene_reference_version`. Each store measured a different subset of the
reference genes, so the others are zero columns.
"""

import shutil
import threading
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import census2_dataset_writer
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr

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
            "dataset_id": pd.Categorical([dataset_id] * n_obs),
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
) -> None:
    """One store and its summary record, as the writer writes them."""
    version_id = f"v-{dataset_id}"
    path = prefix / f"{version_id}.zarr"
    adata.write_zarr(path)
    record: census2_dataset_writer.SummaryRecord = {
        "version": census2_dataset_writer.SUMMARY_VERSION,
        "dataset_id": dataset_id,
        "dataset_version_id": version_id,
        "gene_reference_version": "GENCODE v48",
    }
    set_summary(path, dict(record))


def set_summary(path: Path, record: dict[str, Any] | None) -> None:
    """Set or remove the summary record in the store's root attributes.

    It writes the record as given, unchecked, so that a test can write a malformed one.
    """
    namespace = census2_dataset_writer.SUMMARY_NAMESPACE
    root = zarr.open_group(path, mode="r+")
    attrs = dict(root.attrs)
    attrs.pop(namespace, None)
    if record is not None:
        attrs[namespace] = {census2_dataset_writer.SUMMARY_KEY: record}
    root.attrs.clear()
    root.update_attributes(attrs)


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("corpus")
    prefix = root / "out"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_one(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
    frame = build_catalog.build([str(prefix)], root / "catalog.parquet")
    assert len(frame) == 3
    return {"root": root, "prefix": prefix, "catalog": root / "catalog.parquet"}


def one_store_catalog(
    tmp_path: Path,
    adata: anndata.AnnData,
    dataset_id: str,
) -> Path:
    prefix = tmp_path / "out"
    prefix.mkdir()
    write_one(prefix, adata, dataset_id)
    build_catalog.build([str(prefix)], tmp_path / "catalog.parquet")
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
        obs.index = pd.Index(dataset_id + "_" + obs["observation_joinid"].astype(str))
        var = cast("pd.DataFrame", sub.var)
        parts.append(anndata.AnnData(X=X(sub), obs=obs, var=var))
    return anndata.concat(parts, join="inner")


def query(census: census2_sdk.Census, obs_filter: Any) -> anndata.AnnData:
    """A query that must select cells: `get_anndata` returns None when none match."""
    result = census.get_anndata("homo_sapiens", obs_filter)
    assert result is not None
    return result


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
    result = query(census, obs_filter)
    assert list(result.var_names) == REFERENCE_IDS
    assert result.uns["census2"]["datasets_in_catalog"] == 3
    same(result, expected(corpus, obs_filter))


def test_result_shape_and_provenance(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = query(census, {"tissue": "lung"})
    assert result.uns["census2"]["stores_opened"] == 2
    assert result.uns["census2"]["sources"] == [str(corpus["catalog"])]
    assert result.obs.index.is_unique
    assert X(result).dtype == np.float32
    # var is the stores' var, with its reference columns.
    assert list(result.var.columns) == list(reference_var(REFERENCE_IDS).columns)


def test_concurrent_reads_match_one_at_a_time(
    corpus: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    obs_filter = {"tissue": ["lung", "kidney"]}  # selects all three stores
    census = census2_sdk.open_census([corpus["catalog"]])
    one = query(census, obs_filter)

    # Each read waits for all three, so a serial read times out. Each store then
    # waits for the next one in catalog order, so the first store finishes last.
    read_rows = census2_sdk._read_rows
    order = list(STORES)
    all_running = threading.Barrier(len(order), timeout=10)
    done = {dataset_id: threading.Event() for dataset_id in order}
    finished: list[str] = []

    def overlapping(opened: Any, read: Any, normalized: Any) -> anndata.AnnData:
        all_running.wait()
        part = read_rows(opened, read, normalized)
        position = order.index(read.row["dataset_id"])
        if position + 1 < len(order):
            assert done[order[position + 1]].wait(timeout=10)
        finished.append(str(read.row["dataset_id"]))
        done[str(read.row["dataset_id"])].set()
        return part

    monkeypatch.setattr(census2_sdk, "_read_rows", overlapping)
    many = census.get_anndata("homo_sapiens", obs_filter, max_workers=len(STORES))
    assert many is not None
    assert finished == list(reversed(STORES))
    assert list(many.obs_names) == list(one.obs_names)  # catalog order
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", many.obs), cast("pd.DataFrame", one.obs)
    )
    assert np.array_equal(X(many).toarray(), X(one).toarray())
    assert many.uns["census2"] == one.uns["census2"]


@pytest.mark.parametrize("max_workers", [0, -1])
def test_max_workers_must_be_positive(
    corpus: dict[str, Path], max_workers: int
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    with pytest.raises(ValueError, match="max_workers"):
        census.get_anndata("homo_sapiens", {"tissue": "lung"}, max_workers=max_workers)


def test_empty_result(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    assert census.get_anndata("homo_sapiens", {"tissue": "brain"}) is None


def test_over_exclusion_is_visible(corpus: dict[str, Path]) -> None:
    """The check the harness relies on: drop a true candidate, the result shrinks."""
    census = census2_sdk.open_census([corpus["catalog"]])
    full = query(census, {"tissue": "lung"})
    census.catalog = census.catalog[census.catalog["dataset_id"] != "ds-a"]
    mutated = query(census, {"tissue": "lung"})
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
        query(census, {"tissue": "lung"})
    # Excluding the spatial store makes the rest compose.
    assert query(census, {"dataset_id": "ds-b"}).n_obs == 25


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
    build_catalog.build([str(prefix)], tmp_path / "catalog.parquet")
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match="v-ds-x.zarr: var differs"):
        census.get_anndata("homo_sapiens")


def test_row_without_a_recorded_reference_refused(tmp_path: Path) -> None:
    """A catalog row that names no reference is refused.

    The builder refuses a summary record without a reference, so the test removes the
    value from the catalog row.
    """
    adata = make_store("ds-u", 5, REFERENCE_IDS, ["lung"], ["normal"], 7)
    path = one_store_catalog(tmp_path, adata, "ds-u")
    frame = catalog.load(path)
    frame["gene_reference_version"] = None
    frame.to_parquet(path, index=False)
    census = census2_sdk.open_census([path])
    with pytest.raises(census2_sdk.SourceConflict, match="no gene_reference_version"):
        census.get_anndata("homo_sapiens")


def test_absent_column_selects_no_cells(corpus: dict[str, Path]) -> None:
    """The catalog keeps a store whose obs lacks the column; the read then matches none.

    The result is None, as when the catalog excludes every store.
    """
    census = census2_sdk.open_census([corpus["catalog"]])
    assert census.get_anndata("homo_sapiens", {"tissue_general": "lung"}) is None


def test_catalog_without_a_column_keeps_the_dataset(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    """A catalog that lacks an indexed column proves nothing about it (design 2.2).

    Alone, and mixed with a catalog that has the column, the dataset stays in and
    the read decides. An unrecorded `spatial` value refuses a result of two or
    more stores. An unrecorded `gene_reference_version` refuses the row.
    """
    frame = catalog.load(corpus["catalog"])
    full = frame[frame["dataset_id"] != "ds-c"]
    partial = frame[frame["dataset_id"] == "ds-c"]
    full.to_parquet(tmp_path / "full.parquet", index=False)
    partial.drop(columns=["tissue"]).to_parquet(
        tmp_path / "no-tissue.parquet", index=False
    )
    partial.drop(columns=["spatial"]).to_parquet(
        tmp_path / "no-spatial.parquet", index=False
    )
    partial.drop(columns=["gene_reference_version"]).to_parquet(
        tmp_path / "no-version.parquet", index=False
    )

    # The catalog keeps the store; its cells hold no lung, so the result is None.
    alone = census2_sdk.open_census([tmp_path / "no-tissue.parquet"])
    lung = catalog.normalize({"tissue": "lung"})
    assert len(catalog.exclude(alone.catalog, "homo_sapiens", lung)) == 1
    assert alone.get_anndata("homo_sapiens", {"tissue": "lung"}) is None

    mixed = census2_sdk.open_census(
        [tmp_path / "full.parquet", tmp_path / "no-tissue.parquet"]
    )
    result = query(mixed, {"tissue": "lung"})
    assert result.uns["census2"]["stores_opened"] == 3
    assert result.n_obs == expected(corpus, {"tissue": "lung"}).n_obs
    assert query(mixed, {"tissue": "kidney"}).n_obs == 20

    unknown = census2_sdk.open_census(
        [tmp_path / "full.parquet", tmp_path / "no-spatial.parquet"]
    )
    assert query(unknown, {"tissue": "kidney"}).n_obs == 20
    with pytest.raises(census2_sdk.SourceConflict, match=r"no spatial.*ds-c"):
        unknown.get_anndata("homo_sapiens", {"dataset_id": ["ds-a", "ds-c"]})

    unpinned = census2_sdk.open_census([tmp_path / "no-version.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match=r"no gene_reference.*ds-c"):
        unpinned.get_anndata("homo_sapiens", {"tissue": "kidney"})


def test_filter_values_can_be_any_iterable(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    generator = (t for t in ["lung", "kidney"])
    result = query(census, {"tissue": generator})
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
    assert census.get_anndata("homo_sapiens", {"tissue": "nan"}) is None
    assert query(census, {"tissue": "lung"}).n_obs == 4


def test_mask_by_category_equals_mask_by_cell() -> None:
    """The categorical path and the per-cell fallback apply the same rule."""
    disease = ["normal || flu", "flu", None, "normal", "cold"]
    obs = pd.DataFrame({"disease": pd.Categorical(disease)})
    wanted = {"disease": frozenset({"flu"})}
    expected = np.array([True, True, False, False, False])
    assert (census2_sdk._mask(obs, wanted) == expected).all()
    obs["disease"] = obs["disease"].astype(object)
    assert (census2_sdk._mask(obs, wanted) == expected).all()
    empty = pd.DataFrame({"disease": pd.Categorical([None, None])})
    assert not census2_sdk._mask(empty, wanted).any()


def test_recorded_reference_must_match_the_pin(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    frame = catalog.load(corpus["catalog"])
    frame["gene_reference_version"] = "GENCODE v47"
    frame.to_parquet(tmp_path / "v47.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "v47.parquet"])
    with pytest.raises(census2_sdk.SourceConflict, match="v47.*pinned.*v48"):
        query(census, {"tissue": "lung"})
    frame["gene_reference_version"] = "GENCODE v48"
    frame.to_parquet(tmp_path / "v48.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "v48.parquet"])
    assert query(census, {"tissue": "lung"}).n_obs > 0


# --- the cell counts and count_cells (tracker 5.9) -----------------------------

COUNT_FILTERS: list[dict[str, Any]] = [
    {},
    {"tissue": "lung"},
    {"disease": "flu"},  # a member of a multi-valued value
    {"tissue": ["lung", "kidney"], "disease": "normal"},
    {"dataset_id": "ds-b", "tissue": "lung"},
    {"assay_ontology_term_id": "EFO:0009922"},
    {"tissue": "brain"},  # no store holds the value
    {"cell_type": "T cell"},  # no store holds the column
]


@pytest.mark.parametrize("obs_filter", COUNT_FILTERS)
def test_count_cells_equals_the_query(
    corpus: dict[str, Path], obs_filter: dict[str, Any]
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata("homo_sapiens", obs_filter)
    want = 0 if result is None else result.n_obs
    assert census.count_cells("homo_sapiens", obs_filter) == want


def test_cell_counts_add_up_to_n_obs(corpus: dict[str, Path]) -> None:
    rows = catalog.load(corpus["catalog"])
    assert catalog.COUNTS not in rows.columns
    counts = catalog.load_counts(corpus["catalog"])
    assert list(counts.columns) == [
        "dataset_version_id",
        *catalog.COUNTED,
        "is_primary_data",
        "n_cells",
    ]
    # No store holds cell_type, and the column keeps its declared string type.
    assert counts["cell_type"].isna().all()
    assert str(counts["cell_type"].dtype).startswith("string")
    totals = counts.groupby("dataset_version_id")["n_cells"].sum().to_dict()
    assert totals == dict(zip(rows["dataset_version_id"], rows["n_obs"], strict=True))


def test_count_cells_with_nulls_and_primary_only(tmp_path: Path) -> None:
    adata = make_store("ds-p", 6, REFERENCE_IDS[:3], ["lung"], ["normal"], 7)
    tissue = adata.obs["tissue"].cat.add_categories(["liver"])
    tissue.iloc[:2] = None
    adata.obs["tissue"] = tissue
    adata.obs["is_primary_data"] = np.array([True, False, True, False, True, True])
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-p")])
    assert census.count_cells("homo_sapiens") == 6
    assert census.count_cells("homo_sapiens", {"tissue": "lung"}) == 4
    assert census.count_cells("homo_sapiens", primary_only=True) == 4
    assert (
        census.count_cells("homo_sapiens", {"tissue": "lung"}, primary_only=True) == 3
    )
    assert census.count_cells("mus_musculus") == 0
    # For an organism with no store, a query with a filter gives zero and None.
    assert census.count_cells("mus_musculus", {"tissue": "lung"}) == 0
    assert census.get_anndata("mus_musculus", {"tissue": "lung"}) is None


def test_count_cells_reads_the_catalog_that_open_census_read(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    """A catalog that changes after `open_census` does not change the counts."""
    want = census2_sdk.open_census([corpus["catalog"]]).count_cells("homo_sapiens")
    path = tmp_path / "catalog.parquet"
    shutil.copy(corpus["catalog"], path)
    census = census2_sdk.open_census([path])
    # The builder writes a catalog of one small store over the same path.
    one_store_catalog(
        tmp_path, make_store("ds-z", 2, REFERENCE_IDS, ["lung"], ["normal"], 15), "ds-z"
    )
    assert census2_sdk.open_census([path]).count_cells("homo_sapiens") == 2
    assert census.count_cells("homo_sapiens") == want


def test_count_cells_refuses_an_uncounted_column(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    with pytest.raises(ValueError, match="development_stage"):
        census.count_cells("homo_sapiens", {"development_stage": "adult"})


def test_count_cells_needs_the_counts_column(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    catalog.load(corpus["catalog"]).to_parquet(tmp_path / "bare.parquet", index=False)
    census = census2_sdk.open_census([tmp_path / "bare.parquet"])
    with pytest.raises(ValueError, match="no cell_counts column"):
        census.count_cells("homo_sapiens")


# --- the catalog builder and the summary record -------------------------------


def builder_prefix(tmp_path: Path) -> tuple[Path, Path]:
    adata = make_store("ds-s", 5, REFERENCE_IDS, ["lung"], ["normal"], 9)
    prefix = tmp_path / "out"
    prefix.mkdir()
    write_one(prefix, adata, "ds-s")
    return prefix, prefix / "v-ds-s.zarr"


def test_builder_reads_the_summary_record(tmp_path: Path) -> None:
    prefix, _ = builder_prefix(tmp_path)
    row = build_catalog.build([str(prefix)]).iloc[0]
    assert row["dataset_id"] == "ds-s"
    assert row["dataset_version_id"] == "v-ds-s"
    assert row["gene_reference_version"] == "GENCODE v48"
    assert row["curation_schema_version"] == "7.1.0"
    assert not list(prefix.glob("*.json"))


def test_builder_refuses_a_store_without_a_summary_record(tmp_path: Path) -> None:
    prefix, store = builder_prefix(tmp_path)
    set_summary(store, None)
    with pytest.raises(ValueError, match="v-ds-s.zarr: no summary record"):
        build_catalog.build([str(prefix)])


def test_builder_refuses_an_unknown_record_version(tmp_path: Path) -> None:
    prefix, store = builder_prefix(tmp_path)
    record = {
        "version": 2,
        "dataset_id": "ds-s",
        "dataset_version_id": "v-ds-s",
        "gene_reference_version": "GENCODE v48",
    }
    set_summary(store, record)
    with pytest.raises(ValueError, match="summary record version 2"):
        build_catalog.build([str(prefix)])


def test_builder_refuses_a_record_for_another_store(tmp_path: Path) -> None:
    prefix, store = builder_prefix(tmp_path)
    record = {
        "version": 1,
        "dataset_id": "ds-s",
        "dataset_version_id": "v-other",
        "gene_reference_version": "GENCODE v48",
    }
    set_summary(store, record)
    with pytest.raises(ValueError, match="names dataset version 'v-other'"):
        build_catalog.build([str(prefix)])


def test_builder_refuses_a_malformed_record(tmp_path: Path) -> None:
    prefix, store = builder_prefix(tmp_path)
    record = {
        "version": 1,
        "dataset_id": 7,
        "dataset_version_id": "v-ds-s",
        "gene_reference_version": "GENCODE v48",
    }
    set_summary(store, record)
    with pytest.raises(ValueError, match="v-ds-s.zarr: summary record dataset_id is 7"):
        build_catalog.build([str(prefix)])


def test_builder_marks_a_visium_store_spatial(tmp_path: Path) -> None:
    """The spatial flag follows the writer's list of spatial assays."""
    adata = make_store("ds-v", 5, REFERENCE_IDS, ["lung"], ["normal"], 11)
    visium = "EFO:0022857"  # Visium Spatial Gene Expression V1
    assert visium in census2_dataset_writer.assays.SPATIAL
    adata.obs["assay_ontology_term_id"] = pd.Categorical([visium] * adata.n_obs)
    prefix = tmp_path / "out"
    prefix.mkdir()
    write_one(prefix, adata, "ds-v")
    assert bool(build_catalog.build([str(prefix)]).iloc[0]["spatial"]) is True


# --- the dataset_id column ------------------------------------------------------


def test_store_without_a_dataset_id_column_is_refused(tmp_path: Path) -> None:
    adata = make_store("ds-n", 5, REFERENCE_IDS, ["lung"], ["normal"], 13)
    adata.obs = cast("pd.DataFrame", adata.obs).drop(columns="dataset_id")
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-n")])
    with pytest.raises(ValueError, match="v-ds-n.zarr: obs has no dataset_id column"):
        census.get_anndata("homo_sapiens")


def test_store_for_another_dataset_is_refused(tmp_path: Path) -> None:
    adata = make_store("ds-m", 5, REFERENCE_IDS, ["lung"], ["normal"], 14)
    adata.obs["dataset_id"] = pd.Categorical(["ds-other"] * adata.n_obs)
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-m")])
    with pytest.raises(ValueError, match="holds \\['ds-other'\\].*names ds-m"):
        census.get_anndata("homo_sapiens")
