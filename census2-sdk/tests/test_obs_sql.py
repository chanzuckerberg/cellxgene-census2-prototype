"""Tests for the SQL path: DuckDB over the `obs` fragments. No network, no AWS.

The fixture writes three small stores with `census2_dataset_writer.write`, so each store
has its `obs` fragment, and builds the catalog. A query through DuckDB must give the
result of the same query through the stores, value for value, and a SQL filter must
select the cells that the same test in pandas selects.
"""

import shutil
from pathlib import Path
from typing import Any, cast

import anndata
import build_catalog
import census2_dataset_writer
import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
import scipy.sparse as sp

import census2_sdk
from census2_sdk import catalog, fragments, obs_sql

REFERENCE_IDS = [f"ENSG{j:011d}" for j in range(20)]

STORES = {
    # dataset_id: (n_obs, tissues, cell types, diseases, assay)
    "ds-a": (30, ["lung", "liver"], ["macrophage", "T cell"], ["normal"], "EFO:0009922"),
    "ds-b": (25, ["lung"], ["macrophage", "B cell"], ["flu || normal", "flu"], "EFO:0009922"),
    # Smart-seq2, a full-gene assay: the normalized layer divides by the gene length.
    "ds-c": (20, ["kidney"], ["T cell", "CD8 T cell"], ["normal"], "EFO:0008931"),
}  # fmt: skip

GENERAL = {"lung": "respiratory system", "liver": "liver", "kidney": "kidney"}
CELL_TYPE_IDS = {
    "macrophage": "CL:0000235",
    "T cell": "CL:0000084",
    "CD8 T cell": "CL:0000625",
    "B cell": "CL:0000236",
}


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
    cell_type = labels(cell_types)
    nnz = np.diff(x.indptr).astype(np.int64)
    variance = np.where(nnz > 7, rng.random(n_obs), np.nan)
    obs = pd.DataFrame(
        {
            "dataset_id": pd.Categorical([dataset_id] * n_obs),
            "tissue": tissue,
            "cell_type": cell_type,
            "cell_type_ontology_term_id": cell_type.map(CELL_TYPE_IDS),
            "disease": labels(diseases),
            "donor_id": labels([f"{dataset_id}-d1", f"{dataset_id}-d2"]),
            "assay_ontology_term_id": pd.Categorical([assay] * n_obs),
            "is_primary_data": rng.random(n_obs) < 0.7,
            "observation_joinid": [f"j{i}" for i in range(n_obs)],
            "tissue_general": tissue.map(GENERAL),
            "raw_sum": np.asarray(x.sum(axis=1), dtype=np.float64).ravel(),
            "nnz": nnz,
            # A float column with NaN, which the fragment holds as a null.
            "raw_variance_nnz": variance,
        },
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = reference_var(REFERENCE_IDS)
    var["measured"] = np.ones(len(var), dtype=bool)
    adata = anndata.AnnData(X=sp.csr_matrix(x), obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    return adata


def write_store(prefix: Path, adata: anndata.AnnData, dataset_id: str) -> None:
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
    # A small descendants table: T cell and its one descendant here.
    closure = pa.table(
        {
            "ontology": ["CL"] * 5,
            "ancestor": ["CL:0000084", "CL:0000084", "CL:0000625", "CL:0000235", "CL:0000236"],
            "descendant": ["CL:0000084", "CL:0000625", "CL:0000625", "CL:0000235", "CL:0000236"],
        }
    )  # fmt: skip
    pq.write_table(closure, root / "descendants.parquet")
    return {"root": root, "prefix": prefix, "catalog": root / "catalog.parquet"}


@pytest.fixture(scope="module")
def stores(corpus: dict[str, Path]) -> census2_sdk.Census:
    return census2_sdk.open_census([corpus["catalog"]])


@pytest.fixture(scope="module")
def sql(corpus: dict[str, Path]) -> census2_sdk.Census:
    options = obs_sql.Options(
        threads=2, descendants=str(corpus["root"] / "descendants.parquet")
    )
    return census2_sdk.open_census([corpus["catalog"]], obs_engine=options)


def layer_names(adata: anndata.AnnData) -> list[str]:
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
        # As strings, with a null equal to a null.
        pd.testing.assert_series_equal(
            a.obs[column].astype(str), b.obs[column].astype(str), obj=column
        )
        assert a.obs[column].dtype.kind == b.obs[column].dtype.kind, column
        is_category = isinstance(a.obs[column].dtype, pd.CategoricalDtype)
        assert is_category == isinstance(b.obs[column].dtype, pd.CategoricalDtype)
    assert cast("pd.DataFrame", a.var).equals(cast("pd.DataFrame", b.var))
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
def test_duckdb_gives_the_result_of_the_stores(
    stores: census2_sdk.Census,
    sql: census2_sdk.Census,
    obs_filter: dict[str, Any] | None,
) -> None:
    """Value for value, with the normalized layer; ds-c is a full-gene assay."""
    assert_same(
        stores.get_anndata("homo_sapiens", obs_filter, layers=["normalized"]),
        sql.get_anndata(
            "homo_sapiens", obs_filter, layers=["normalized"], max_workers=3
        ),
    )


@pytest.mark.parametrize("files", ["all", "glob"])
@pytest.mark.parametrize("obs_filter", MAPPINGS[:5])
def test_each_way_to_find_the_files_gives_the_same_result(
    corpus: dict[str, Path],
    sql: census2_sdk.Census,
    files: str,
    obs_filter: dict[str, Any] | None,
) -> None:
    pattern = str(corpus["prefix"].resolve() / "*.obs.parquet")
    options = obs_sql.Options(files=files, glob=pattern)
    other = census2_sdk.open_census([corpus["catalog"]], obs_engine=options)
    assert_same(
        sql.get_anndata("homo_sapiens", obs_filter),
        other.get_anndata("homo_sapiens", obs_filter),
    )


def all_cells(stores: census2_sdk.Census) -> pd.DataFrame:
    adata = stores.get_anndata("homo_sapiens")
    assert adata is not None
    return cast("pd.DataFrame", adata.obs)


def members_of(value: str) -> set[str]:
    return set(value.split(" || "))


# Each filter with the same test in pandas, over the obs of all cells.
FILTERS: list[tuple[str, Any]] = [
    ("cell_type = 'T cell' OR tissue = 'liver'", lambda o: (o.cell_type == "T cell") | (o.tissue == "liver")),
    ("tissue = 'lung' AND NOT (cell_type = 'macrophage')", lambda o: (o.tissue == "lung") & (o.cell_type != "macrophage")),
    ("cell_type NOT IN ('T cell', 'B cell')", lambda o: ~o.cell_type.isin(["T cell", "B cell"])),
    ("(tissue = 'lung' OR tissue = 'kidney') AND NOT (is_primary_data AND donor_id LIKE '%d1')", lambda o: o.tissue.isin(["lung", "kidney"]) & ~(o.is_primary_data & o.donor_id.astype(str).str.endswith("d1"))),
    ("donor_id LIKE 'ds-a%'", lambda o: o.donor_id.astype(str).str.startswith("ds-a")),
    ("regexp_matches(donor_id, '-d[2-9]$') AND is_primary_data = false", lambda o: o.donor_id.astype(str).str.endswith("-d2") & ~o.is_primary_data),
    ("lower(cell_type) ILIKE '%T CELL'", lambda o: o.cell_type.astype(str).str.lower().str.endswith("t cell")),
    # The disease column: string semantics, and explicit membership.
    ("disease = 'flu'", lambda o: o.disease == "flu"),
    ("list_contains(string_split(disease, ' || '), 'flu')", lambda o: o.disease.astype(str).map(lambda v: "flu" in members_of(v))),
    ("NOT list_has_any(string_split(disease, ' || '), ['flu'])", lambda o: o.disease.astype(str).map(lambda v: "flu" not in members_of(v))),
    ("list_has_all(string_split(disease, ' || '), ['flu', 'normal'])", lambda o: o.disease.astype(str).map(lambda v: {"flu", "normal"} <= members_of(v))),
    # Numbers and nulls.
    ("nnz >= 9 AND raw_sum < 60", lambda o: (o.nnz >= 9) & (o.raw_sum < 60)),
    ("nnz BETWEEN 6 AND 8 OR tissue_general = 'kidney'", lambda o: o.nnz.between(6, 8) | (o.tissue_general == "kidney")),
    ("raw_sum / nnz > 6.0", lambda o: o.raw_sum / o.nnz > 6.0),
    ("raw_variance_nnz IS NULL", lambda o: o.raw_variance_nnz.isna()),
    ("raw_variance_nnz > 0.5", lambda o: o.raw_variance_nnz > 0.5),
    ("NOT (raw_variance_nnz > 0.5)", lambda o: o.raw_variance_nnz <= 0.5),
    ("coalesce(raw_variance_nnz, 0) < 0.5", lambda o: o.raw_variance_nnz.fillna(0) < 0.5),
    # Ontology closure through the descendants table.
    ("cell_type_ontology_term_id IN (SELECT descendant FROM descendants WHERE ancestor = 'CL:0000084')", lambda o: o.cell_type.isin(["T cell", "CD8 T cell"])),
    ("cell_type_ontology_term_id NOT IN (SELECT descendant FROM descendants WHERE ancestor = 'CL:0000084') AND tissue = 'lung'", lambda o: ~o.cell_type.isin(["T cell", "CD8 T cell"]) & (o.tissue == "lung")),
    ("cell_type_ontology_term_id IN (SELECT descendant FROM descendants WHERE ancestor = 'CL:9999999')", lambda o: o.tissue == "nowhere"),
]  # fmt: skip


@pytest.mark.parametrize(("text", "expected"), FILTERS)
def test_a_filter_selects_the_cells_that_satisfy_it(
    stores: census2_sdk.Census, sql: census2_sdk.Census, text: str, expected: Any
) -> None:
    obs = all_cells(stores)
    want = obs.index[np.asarray(expected(obs), dtype=bool)]
    adata = sql.get_anndata("homo_sapiens", text)
    got = pd.Index([]) if adata is None else adata.obs_names
    assert list(got) == list(want)


@pytest.mark.parametrize(("text", "expected"), FILTERS)
def test_the_catalog_excludes_only_datasets_without_a_match(
    stores: census2_sdk.Census, sql: census2_sdk.Census, text: str, expected: Any
) -> None:
    obs = all_cells(stores)
    matching = set(obs["dataset_id"][np.asarray(expected(obs), dtype=bool)].astype(str))
    engine = sql.sql_engine
    rows, _ = engine.candidates("homo_sapiens", engine.predicate(text))
    assert matching <= set(rows["dataset_id"].astype(str))


@pytest.mark.parametrize("obs_filter", MAPPINGS)
def test_the_catalog_test_of_a_mapping_is_the_catalog_rule(
    stores: census2_sdk.Census,
    sql: census2_sdk.Census,
    obs_filter: dict[str, Any] | None,
) -> None:
    engine = sql.sql_engine
    rows, _ = engine.candidates("homo_sapiens", engine.predicate(obs_filter))
    excluded = catalog.exclude(stores.catalog, "homo_sapiens", obs_filter)
    assert list(rows["dataset_id"]) == list(excluded["dataset_id"])


def test_the_catalog_test_excludes_by_one_term_and_keeps_a_conjunction(
    sql: census2_sdk.Census,
) -> None:
    """ds-c holds T cells and only kidney: one term excludes it, a pair cannot."""
    engine = sql.sql_engine

    def kept(text: str) -> list[str]:
        rows, _ = engine.candidates("homo_sapiens", engine.predicate(text))
        return sorted(rows["dataset_id"].astype(str))

    assert kept("tissue = 'kidney'") == ["ds-c"]
    assert kept("NOT (tissue = 'kidney')") == ["ds-a", "ds-b"]
    assert kept("tissue <> 'kidney' AND cell_type = 'CD8 T cell'") == []
    assert kept("list_contains(string_split(disease, ' || '), 'flu')") == ["ds-b"]
    assert kept("nnz > 1000") == ["ds-a", "ds-b", "ds-c"]  # no numeric range yet
    assert kept("raw_sum / nnz > 1 AND tissue LIKE 'kid%'") == ["ds-c"]
    assert kept("tissue IS NULL") == ["ds-a", "ds-b", "ds-c"]


def test_a_store_without_a_selected_cell_is_not_opened(
    sql: census2_sdk.Census, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ds-a passes the catalog test of each term, but holds no liver B cell."""
    opened: list[str] = []
    original = catalog.open_store

    def spy(uri: Any) -> Any:
        opened.append(Path(str(uri)).name)
        return original(uri)

    monkeypatch.setattr(census2_sdk, "open_store", spy)
    monkeypatch.setattr(catalog, "open_store", spy)
    adata = sql.get_anndata(
        "homo_sapiens", "tissue IN ('lung', 'liver') AND cell_type = 'B cell'"
    )
    assert adata is not None
    assert adata.uns["census2"]["catalog_candidates"] == 1
    assert set(opened) == {"v-ds-b.zarr"}
    none = sql.get_anndata("homo_sapiens", "tissue = 'liver' AND cell_type = 'B cell'")
    assert none is None
    assert sql.last_scan[0]["catalog_candidates"] == 0


def test_get_obs_and_count_obs_open_no_store(
    stores: census2_sdk.Census,
    sql: census2_sdk.Census,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obs = all_cells(stores)

    def refuse(uri: Any) -> Any:
        raise AssertionError(f"opened {uri}")

    monkeypatch.setattr(census2_sdk, "open_store", refuse)
    monkeypatch.setattr(catalog, "open_store", refuse)
    frame = sql.get_obs("homo_sapiens", "tissue = 'lung'", columns=["donor_id"])
    assert list(frame.index) == list(obs.index[obs.tissue == "lung"])
    names = ["dataset_id", "observation_joinid", "donor_id"]
    assert list(frame.columns) == [*names, *fragments.ROW_COLUMNS]
    assert isinstance(frame["donor_id"].dtype, pd.CategoricalDtype)
    every = sql.get_obs("homo_sapiens")
    assert every.index.equals(obs.index)
    for column in obs.columns:
        pd.testing.assert_series_equal(
            every[column].astype(str), obs[column].astype(str), check_names=False
        )

    lung = obs[obs.tissue == "lung"]
    want = lung.groupby("donor_id", observed=True).size()
    counted = sql.count_obs("homo_sapiens", "tissue = 'lung'", by=["donor_id"])
    got = counted.set_index("donor_id")["n_cells"]
    assert got.sort_index().to_dict() == want.sort_index().to_dict()
    total = sql.count_obs("homo_sapiens", {"tissue": "lung"})
    assert int(total["n_cells"].iloc[0]) == len(lung)


REFUSED = [
    "true; DROP TABLE catalog",
    "true UNION SELECT 1",
    "true ORDER BY 1",
    "true LIMIT 1",
    "true GROUP BY 1",
    "getenv('HOME') = ''",
    "current_setting('threads') = 1",
    "read_text('/etc/passwd') = ''",
    "x IN (SELECT 1)",
    "x = (SELECT max(nnz) FROM obs)",
    "x IN (SELECT descendant FROM descendants, catalog)",
    "x IN (SELECT descendant FROM descendants WHERE ancestor IN (SELECT 1))",
    "x IN (SELECT ancestor FROM descendants)",
    "obs.tissue = 'lung'",
    "count(*) > 0",
    "nextval('s') > 0",
    "random() < 0.5",
    "tissue = ?",
    "CAST(tissue AS STRUCT(a INT)) IS NULL",
    "list_transform([1], lambda v: v + 1) = [2]",
]


@pytest.mark.parametrize("text", REFUSED)
def test_the_language_refuses_what_is_not_a_filter(text: str) -> None:
    with pytest.raises(obs_sql.FilterError):
        obs_sql.parse(text)


def test_a_filter_on_an_unknown_column_is_a_clear_error(
    sql: census2_sdk.Census,
) -> None:
    with pytest.raises(obs_sql.FilterError, match="cel_type"):
        sql.get_anndata("homo_sapiens", "cel_type = 'T cell'")


def test_render_gives_the_tree_that_parse_reads() -> None:
    for text, _ in FILTERS:
        first = obs_sql.parse(text)
        again = obs_sql.parse(first.text)
        assert again.text == first.text


def test_the_connection_reads_no_other_file(
    corpus: dict[str, Path], sql: census2_sdk.Census, tmp_path: Path
) -> None:
    """The second barrier, under the filter language: DuckDB's own file access rules.

    The rules allow reads and writes in the fragment directories. Only the filter
    language keeps a `COPY` out of a query.
    """
    con = sql.sql_engine.con
    with pytest.raises(duckdb.PermissionException):
        con.execute(f"SELECT * FROM read_csv('{corpus['root'] / 'catalog.parquet'}')")
    with pytest.raises(duckdb.Error):
        con.execute("SET enable_external_access = true")
    with pytest.raises(duckdb.PermissionException):
        con.execute(f"COPY catalog TO '{tmp_path / 'copy.parquet'}'")


def test_a_fragment_of_an_older_writer_gives_nulls(
    corpus: dict[str, Path], tmp_path: Path
) -> None:
    """union_by_name: a fragment without `nnz` gives nulls; without it, DuckDB refuses."""
    prefix = tmp_path / "out"
    shutil.copytree(corpus["prefix"], prefix)
    path = prefix / "v-ds-c.obs.parquet"
    pq.write_table(pq.read_table(path).drop_columns(["nnz"]), path)
    build_catalog.build(str(prefix), tmp_path / "catalog.parquet")
    census = census2_sdk.open_census(
        [tmp_path / "catalog.parquet"], obs_engine="duckdb"
    )
    frame = census.get_obs("homo_sapiens", "nnz IS NULL", columns=["nnz"])
    assert set(frame["dataset_id"].astype(str)) == {"ds-c"}
    strict = obs_sql.Options(union_by_name=False)
    census = census2_sdk.open_census([tmp_path / "catalog.parquet"], obs_engine=strict)
    with pytest.raises((obs_sql.FilterError, duckdb.Error)):
        census.get_obs("homo_sapiens", "nnz IS NULL", columns=["nnz"])


def test_from_mapping_states_membership_explicitly() -> None:
    text = obs_sql.from_mapping(
        {"disease": ["flu"], "tissue": "lung", "is_primary_data": True}
    )
    assert text == (
        "list_has_any(string_split(\"disease\", ' || '), ['flu']) AND "
        '"tissue" IN (\'lung\') AND "is_primary_data" IN (true)'
    )


def test_a_census_takes_one_obs_engine(corpus: dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="no obs engine"):
        census2_sdk.open_census([corpus["catalog"]], obs_engine="sqlite")
