"""Tests for the `obs` table path. No network, no AWS.

The fixture writes three small stores with `census2_dataset_writer.write`, so each store
has its `obs` fragment and fragment record. It builds the catalog and the Iceberg table
over the fragments with the two tools. A query through the table must give the result
of the same query through the stores, value for value.
"""

import shutil
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import build_obs_table
import census2_dataset_writer
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from pyiceberg.expressions import parser

import census2_sdk
from census2_sdk import catalog, expressions, obs_table

REFERENCE_IDS = [f"ENSG{j:011d}" for j in range(20)]

STORES = {
    # dataset_id: (n_obs, tissues, cell types, diseases, assay)
    "ds-a": (30, ["lung", "liver"], ["macrophage", "T cell"], ["normal"], "EFO:0009922"),
    "ds-b": (25, ["lung"], ["macrophage", "B cell"], ["normal ||  flu", "flu"], "EFO:0009922"),
    # Smart-seq2, a full-gene assay: the normalized layer divides by the gene length.
    "ds-c": (20, ["kidney"], ["T cell"], ["normal"], "EFO:0008931"),
}  # fmt: skip

# The tissue_general of each tissue, as the writer rolls them up.
GENERAL = {"lung": "respiratory system", "liver": "liver", "kidney": "kidney"}


def reference_var(ids: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature_id": ids,
            "feature_name": [f"GENE{j}" for j in range(len(ids))],
            "feature_version": ["1"] * len(ids),
            "feature_length": [100 * (j + 1) for j in range(len(ids))],
            "feature_type": ["gene"] * len(ids),
        },
        index=ids,
    )


def make_store(
    dataset_id: str,
    n_obs: int,
    tissues: list[str],
    cell_types: list[str],
    diseases: list[str],
    assay: str,
    seed: int,
) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    x = sp.random(
        n_obs, len(REFERENCE_IDS), density=0.4, format="csr", dtype="float32",
        random_state=seed,
    )  # fmt: skip
    x.data = np.floor(x.data * 10) + 1

    def labels(values: list[str]) -> Any:
        return pd.Categorical(rng.choice(values, n_obs))

    tissue = labels(tissues)
    obs = pd.DataFrame(
        {
            "dataset_id": pd.Categorical([dataset_id] * n_obs),
            "tissue": tissue,
            "cell_type": labels(cell_types),
            "disease": labels(diseases),
            "donor_id": labels([f"{dataset_id}-d1", f"{dataset_id}-d2"]),
            "assay_ontology_term_id": pd.Categorical([assay] * n_obs),
            "is_primary_data": rng.random(n_obs) < 0.7,
            "observation_joinid": [f"j{i}" for i in range(n_obs)],
            "tissue_general": tissue.map(GENERAL),
            # Two of the count statistics, a float64 and an int64 column.
            "raw_sum": np.asarray(x.sum(axis=1), dtype=np.float64).ravel(),
            "nnz": np.diff(x.indptr).astype(np.int64),
        },
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = reference_var(REFERENCE_IDS)
    # The writer marks the dataset's assayed genes; here every gene.
    var["measured"] = np.ones(len(var), dtype=bool)
    adata = anndata.AnnData(X=sp.csr_matrix(x), obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    return adata


def write_store(prefix: Path, adata: anndata.AnnData, dataset_id: str) -> None:
    """One store with its fragment and records, as the writer writes them."""
    version_id = f"v-{dataset_id}"
    path = prefix / f"{version_id}.zarr"
    census2_dataset_writer.write(
        adata,
        path,
        summary={
            "version": census2_dataset_writer.SUMMARY_VERSION,
            "dataset_id": dataset_id,
            "dataset_version_id": version_id,
            "gene_reference_version": "GENCODE v48",
        },
        fragment=census2_dataset_writer.fragment_location(str(path)),
    )


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("corpus")
    prefix = root / "out"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_store(prefix, make_store(dataset_id, *spec, seed=seed), dataset_id)
    build_catalog.build(str(prefix), root / "catalog.parquet")
    summary = build_obs_table.build(str(root / "catalog.parquet"), str(root / "obs"))
    assert summary["added"] == 3
    return {"root": root, "prefix": prefix, "catalog": root / "catalog.parquet"}


@pytest.fixture(scope="module")
def stores(corpus: dict[str, Path]) -> census2_sdk.Census:
    return census2_sdk.open_census([corpus["catalog"]])


@pytest.fixture(scope="module")
def table(corpus: dict[str, Path]) -> census2_sdk.Census:
    return census2_sdk.open_census(
        [corpus["catalog"]], obs_tables=[corpus["root"] / "obs"]
    )


def layer_names(adata: anndata.AnnData) -> list[str]:
    """The names of the layers. anndata 0.13 also lists `X`, as the layer None."""
    layers = cast("dict[str | None, Any]", adata.layers)
    return sorted(k for k in layers if k is not None)


def assert_same(a: anndata.AnnData | None, b: anndata.AnnData | None) -> None:
    """Equal results, value for value. The categories can differ; see the SDK."""
    assert (a is None) == (b is None)
    if a is None or b is None:
        return
    assert a.obs_names.equals(b.obs_names)
    assert list(a.obs.columns) == list(b.obs.columns)
    for column in a.obs.columns:
        assert np.array_equal(
            a.obs[column].astype(str).to_numpy(), b.obs[column].astype(str).to_numpy()
        ), column
        assert a.obs[column].dtype.kind == b.obs[column].dtype.kind, column
    assert cast("pd.DataFrame", a.var).equals(cast("pd.DataFrame", b.var))
    assert "measured" not in a.var.columns
    names = layer_names(a)
    assert names == layer_names(b)
    pairs = [(a.X, b.X), *((a.layers[k], b.layers[k]) for k in names)]
    for one, other in pairs:
        x_a, x_b = sp.csr_matrix(cast("Any", one)), sp.csr_matrix(cast("Any", other))
        assert x_a.shape == x_b.shape and cast("Any", x_a != x_b).nnz == 0
        assert x_a.dtype == x_b.dtype


MAPPINGS: list[dict[str, Any] | None] = [
    None,
    {"tissue": "lung"},
    {"cell_type": ["B cell", "T cell"]},
    {"tissue": "lung", "cell_type": "macrophage", "is_primary_data": True},
    {"disease": "flu"},
    {"dataset_id": ["ds-c", "ds-a"]},
    {"tissue": "brain"},
]


@pytest.mark.parametrize("obs_filter", MAPPINGS)
def test_the_table_gives_the_result_of_the_stores(
    stores: census2_sdk.Census,
    table: census2_sdk.Census,
    obs_filter: dict[str, Any] | None,
) -> None:
    assert_same(
        stores.get_anndata("homo_sapiens", obs_filter),
        table.get_anndata("homo_sapiens", obs_filter, max_workers=3),
    )


@pytest.mark.parametrize("obs_filter", [None, {"tissue": "lung"}, {"tissue": "kidney"}])
def test_the_table_gives_the_normalized_layer_of_the_stores(
    stores: census2_sdk.Census,
    table: census2_sdk.Census,
    obs_filter: dict[str, Any] | None,
) -> None:
    """The layer is calculated on read on both paths; ds-c is a full-gene assay."""
    want = stores.get_anndata("homo_sapiens", obs_filter, layers=["normalized"])
    got = table.get_anndata(
        "homo_sapiens", obs_filter, layers=["normalized"], max_workers=3
    )
    assert got is not None and "normalized" in got.layers
    assert_same(want, got)


def test_a_layer_that_does_not_exist_is_refused(table: census2_sdk.Census) -> None:
    with pytest.raises(ValueError, match="no layer"):
        table.get_anndata("homo_sapiens", layers=["log1p"])


def test_x_from_the_ranges_equals_x_from_indptr(
    table: census2_sdk.Census, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = {"cell_type": "macrophage"}
    from_ranges = table.get_anndata("homo_sapiens", query)
    monkeypatch.setattr(census2_sdk, "X_FROM_RANGES", False)
    assert_same(from_ranges, table.get_anndata("homo_sapiens", query))


@pytest.mark.parametrize("obs_filter", MAPPINGS)
def test_a_key_first_scan_gives_the_same_result(
    stores: census2_sdk.Census,
    table: census2_sdk.Census,
    obs_filter: dict[str, Any] | None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(census2_sdk, "SCAN_KEY_FIRST", True)
    assert_same(
        stores.get_anndata("homo_sapiens", obs_filter),
        table.get_anndata("homo_sapiens", obs_filter),
    )


def test_a_store_without_a_selected_cell_is_not_opened(
    table: census2_sdk.Census, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ds-c passes the catalog test of each term, but holds no lung T cell."""
    opened: list[str] = []
    original = catalog.open_store

    def spy(uri: Any) -> Any:
        opened.append(Path(str(uri)).name)
        return original(uri)

    monkeypatch.setattr(census2_sdk, "open_store", spy)
    query = "tissue IN ('lung', 'kidney') AND cell_type = 'B cell'"
    adata = table.get_anndata("homo_sapiens", query)
    assert adata is not None
    assert adata.uns["census2"]["catalog_candidates"] == 1
    assert {o for o in opened} == {"v-ds-b.zarr"}
    assert set(adata.obs["dataset_id"].astype(str)) == {"ds-b"}


def all_cells(stores: census2_sdk.Census) -> pd.DataFrame:
    adata = stores.get_anndata("homo_sapiens")
    assert adata is not None
    return cast("pd.DataFrame", adata.obs)


def members_of(value: str) -> set[str]:
    return set(catalog.members(value, "disease"))


# Each expression with the same test in pandas, over the obs of all cells.
EXPRESSIONS: list[tuple[str, Any]] = [
    (
        "cell_type = 'T cell' OR tissue = 'liver'",
        lambda o: (o.cell_type == "T cell") | (o.tissue == "liver"),
    ),
    (
        "tissue = 'lung' AND NOT (cell_type = 'macrophage')",
        lambda o: (o.tissue == "lung") & (o.cell_type != "macrophage"),
    ),
    ("cell_type NOT IN ('T cell', 'B cell')", lambda o: ~o.cell_type.isin(["T cell", "B cell"])),
    ("donor_id LIKE 'ds-a%'", lambda o: o.donor_id.astype(str).str.startswith("ds-a")),
    ("donor_id = 'ds-b-d2' AND is_primary_data = false", lambda o: (o.donor_id == "ds-b-d2") & ~o.is_primary_data),
    ("disease = 'flu'", lambda o: o.disease.astype(str).map(lambda v: "flu" in members_of(v))),
    ("disease != 'flu'", lambda o: o.disease.astype(str).map(lambda v: "flu" not in members_of(v))),
    ("disease IN ('normal')", lambda o: o.disease.astype(str).map(lambda v: "normal" in members_of(v))),
    ("disease NOT IN ('normal', 'flu')", lambda o: o.disease.astype(str).map(lambda v: not members_of(v) & {"normal", "flu"})),
    ("tissue_general = 'respiratory system'", lambda o: o.tissue_general == "respiratory system"),
    ("nnz >= 9", lambda o: o.nnz >= 9),
    ("raw_sum > 40.5 AND nnz < 12", lambda o: (o.raw_sum > 40.5) & (o.nnz < 12)),
    ("tissue_general = 'kidney' OR nnz <= 4", lambda o: (o.tissue_general == "kidney") | (o.nnz <= 4)),
]  # fmt: skip


@pytest.mark.parametrize(("text", "expected"), EXPRESSIONS)
def test_an_expression_selects_the_cells_that_satisfy_it(
    stores: census2_sdk.Census, table: census2_sdk.Census, text: str, expected: Any
) -> None:
    obs = all_cells(stores)
    want = obs.index[np.asarray(expected(obs), dtype=bool)]
    adata = table.get_anndata("homo_sapiens", text)
    got = pd.Index([]) if adata is None else adata.obs_names
    assert list(got) == list(want)


@pytest.mark.parametrize(("text", "expected"), EXPRESSIONS)
def test_the_catalog_excludes_only_datasets_without_a_match(
    stores: census2_sdk.Census, text: str, expected: Any
) -> None:
    obs = all_cells(stores)
    matching = set(obs["dataset_id"][np.asarray(expected(obs), dtype=bool)].astype(str))
    expr = parser.parse(text)
    rows = cast("list[dict[str, Any]]", stores.catalog.to_dict("records"))
    kept = {str(r["dataset_id"]) for r in rows if expressions.can_match(expr, r)}
    assert matching <= kept


@pytest.mark.parametrize("obs_filter", MAPPINGS)
def test_can_match_agrees_with_the_catalog_rule_for_a_mapping(
    stores: census2_sdk.Census, obs_filter: dict[str, Any] | None
) -> None:
    expr = expressions.from_filter(obs_filter)
    rows = cast("list[dict[str, Any]]", stores.catalog.to_dict("records"))
    kept = [str(r["dataset_id"]) for r in rows if expressions.can_match(expr, r)]
    excluded = catalog.exclude(stores.catalog, "homo_sapiens", obs_filter)
    assert kept == [str(d) for d in excluded["dataset_id"]]


def test_get_obs_reads_the_table_and_no_store(
    stores: census2_sdk.Census,
    table: census2_sdk.Census,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obs = all_cells(stores)

    def refuse(uri: Any) -> Any:
        raise AssertionError(f"opened {uri}")

    monkeypatch.setattr(census2_sdk, "open_store", refuse)
    frame = table.get_obs("homo_sapiens", "tissue = 'lung'", columns=["donor_id"])
    assert list(frame.index) == list(obs.index[obs.tissue == "lung"])
    # The cell key is always there. The columns keep the order of the table schema.
    names = ["dataset_id", "donor_id", "observation_joinid"]
    assert list(frame.columns) == [*names, *obs_table.ROW_COLUMNS]
    assert isinstance(frame["donor_id"].dtype, pd.CategoricalDtype)
    assert frame.attrs["census2"]["obs_table"][0]["files_planned"] == 2


def test_count_obs_counts_on_any_column(
    stores: census2_sdk.Census, table: census2_sdk.Census
) -> None:
    obs = all_cells(stores)
    lung = obs[obs.tissue == "lung"]
    want = lung.groupby("donor_id", observed=True).size()
    frame = table.count_obs("homo_sapiens", "tissue = 'lung'", by=["donor_id"])
    got = frame.set_index("donor_id")["n_cells"]
    assert got.sort_index().to_dict() == want.sort_index().to_dict()
    total = table.count_obs("homo_sapiens", {"tissue": "lung"})
    assert int(total["n_cells"].iloc[0]) == len(lung)


def test_an_expression_needs_a_table(stores: census2_sdk.Census) -> None:
    with pytest.raises(ValueError, match="needs an obs table"):
        stores.get_anndata("homo_sapiens", "tissue = 'lung'")


def test_a_table_of_another_catalog_is_refused(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    other = tmp_path / "catalog.parquet"
    frame = pd.read_parquet(corpus["catalog"]).iloc[:2]
    frame.to_parquet(other)
    with pytest.raises(ValueError, match="built from another catalog"):
        census2_sdk.open_census([other], obs_tables=[corpus["root"] / "obs"])


def test_the_table_opens_from_its_metadata_file(corpus: dict[str, Path]) -> None:
    metadata = obs_table.current_metadata(str(corpus["root"] / "obs"))
    assert metadata is not None and metadata.endswith(".metadata.json")
    census = census2_sdk.open_census([corpus["catalog"]], obs_tables=[metadata])
    adata = census.get_anndata("homo_sapiens", {"tissue": "kidney"})
    assert adata is not None and set(adata.obs["dataset_id"].astype(str)) == {"ds-c"}


def test_the_builder_commits_the_change_of_the_catalog(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    full = pd.read_parquet(corpus["catalog"])
    location = str(tmp_path / "obs")
    first = tmp_path / "first.parquet"
    full.iloc[:2].to_parquet(first)
    assert build_obs_table.build(str(first), location)["added"] == 2

    second = tmp_path / "second.parquet"
    full.iloc[1:].to_parquet(second)
    summary = build_obs_table.build(str(second), location)
    assert (summary["added"], summary["deleted"]) == (1, 1)
    census = census2_sdk.open_census([second], obs_tables=[location])
    frame = census.get_obs("homo_sapiens", columns=["dataset_id"])
    assert sorted(set(frame["dataset_id"].astype(str))) == ["ds-b", "ds-c"]
    # The metadata of the first commit stays, so a release that pins it still reads.
    table = obs_table.open_table(location)
    assert len(table.metadata.snapshots) >= 3


def test_the_builder_refuses_a_fragment_of_another_size(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    prefix = tmp_path / "out"
    shutil.copytree(corpus["prefix"], prefix)
    (prefix / "v-ds-a.obs.parquet").write_bytes(b"not the fragment")
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    with pytest.raises(ValueError, match="size 16, the record gives"):
        build_obs_table.build(str(tmp_path / "catalog.parquet"), str(tmp_path / "obs"))


def test_the_builder_refuses_a_store_without_a_fragment(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    frame = pd.read_parquet(corpus["catalog"])
    frame.loc[0, "obs_fragment_bytes"] = None
    frame.to_parquet(tmp_path / "catalog.parquet")
    with pytest.raises(ValueError, match="no obs fragment for \\['ds-a'\\]"):
        build_obs_table.build(str(tmp_path / "catalog.parquet"), str(tmp_path / "obs"))


def test_with_members_rewrites_membership_into_joined_values() -> None:
    rows = [{"dataset_id": "d", "disease": ["flu", "normal ||  flu", "normal"]}]
    expr = expressions.with_members(parser.parse("disease = 'flu'"), rows)
    assert expr == parser.parse("disease IN ('flu', 'normal ||  flu')")
    expr = expressions.with_members(parser.parse("disease != 'flu'"), rows)
    assert expr == parser.parse("disease NOT IN ('flu', 'normal ||  flu')")
    expr = expressions.with_members(parser.parse("disease = 'cold'"), rows)
    assert expr == parser.parse("false")


def test_a_joined_literal_selects_exactly_its_values() -> None:
    """A Census 1 filter names the joined value; its translation keeps that meaning."""
    rows = [{"dataset_id": "d", "disease": ["flu", "normal ||  flu", "flu || normal"]}]
    census1 = expressions.from_value_filter("disease in ['normal || flu']")
    expr = expressions.with_members(census1, rows)
    assert expr == parser.parse("disease IN ('normal ||  flu', 'flu || normal')")


@pytest.mark.parametrize(
    ("value_filter", "text"),
    [
        ("tissue == 'lung' and cell_type != 'T cell'", "(tissue = 'lung' AND cell_type != 'T cell')"),
        ("nnz >= 500 or not (sex in ['male'])", "(nnz >= 500 OR NOT sex = 'male')"),
        ("is_primary_data == True", "is_primary_data = true"),
    ],
)  # fmt: skip
def test_a_value_filter_translates_to_an_expression(
    value_filter: str, text: str
) -> None:
    expr = expressions.from_value_filter(value_filter)
    assert expressions.to_text(expr) == text
    assert parser.parse(text) == expr


@pytest.mark.parametrize(
    "value_filter", ["a == b", "1 < nnz", "tissue in 'lung'", "f(x)"]
)
def test_a_value_filter_without_an_iceberg_term_is_refused(value_filter: str) -> None:
    with pytest.raises(ValueError, match="Iceberg|column|list"):
        expressions.from_value_filter(value_filter)


def test_get_obs_with_arrow_columns_gives_the_same_values(
    table: census2_sdk.Census,
) -> None:
    query = "tissue = 'lung' AND nnz >= 6"
    numpy_frame = table.get_obs("homo_sapiens", query)
    arrow_frame = table.get_obs("homo_sapiens", query, dtype_backend="pyarrow")
    assert len(numpy_frame) and list(arrow_frame.index) == list(numpy_frame.index)
    assert list(arrow_frame.columns) == list(numpy_frame.columns)
    for column in numpy_frame.columns:
        assert isinstance(arrow_frame[column].dtype, pd.ArrowDtype), column
        assert list(arrow_frame[column].astype(str)) == list(
            numpy_frame[column].astype(str)
        ), column
    timing = arrow_frame.attrs["census2"]
    assert timing["table_seconds"] >= 0 and timing["build_seconds"] >= 0
    with pytest.raises(ValueError, match="dtype_backend"):
        table.get_obs("homo_sapiens", query, dtype_backend=cast("Any", "polars"))
