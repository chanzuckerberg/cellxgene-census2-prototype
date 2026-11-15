"""Tests for the query builder, its hints, and its evaluation with pyarrow.

mypy checks this file under `--strict`, with the error code `deprecated` on (see
`pyproject.toml`). Each `# type: ignore[deprecated]` proves that mypy reports a
refused operator: mypy reports an ignore that it does not need.
"""

import itertools
import random
from pathlib import Path
from typing import Any, cast

import anndata
import generate_obs_columns
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]  # pyarrow has no types
import pyarrow.compute as pc  # type: ignore[import-untyped]
import pytest
from test_sdk import STORES, X, corpus, expected, make_store, one_store_catalog

import census2_sdk
from census2_sdk import obs_read, value_sets
from census2_sdk.catalog import open_store
from census2_sdk.obs_columns import KINDS
from census2_sdk.query import Query, as_query, members_pattern

__all__ = ["corpus"]  # the fixture of test_sdk

hs = census2_sdk.Organism("homo_sapiens")
obs = hs.obs


def test_obs_columns_match_the_store_format() -> None:
    """The generated copy of the schema has each store column, with its class."""
    assert generate_obs_columns.current() == generate_obs_columns.expected()


# --- hints and columns ----------------------------------------------------------


def test_a_term_records_its_column_and_its_hint() -> None:
    q = obs.cell_type.isin(["neutrophil", "B cell"]) & (obs.nnz >= 500)
    assert q.columns == {"cell_type", "nnz"}
    assert q.hints == {"cell_type": {"B cell", "neutrophil"}}
    assert not q.exact  # the numeric term is not a hint
    assert isinstance(q.expression, pc.Expression)
    exact = (obs.tissue == "lung") & obs.disease.has_any(["COVID-19"])
    assert exact.exact
    assert exact.hints == {"tissue": {"lung"}, "disease": {"COVID-19"}}


def test_and_intersects_the_values_of_one_column() -> None:
    q = obs.tissue.isin(["lung", "liver"]) & (obs.tissue == "lung")
    assert q.hints == {"tissue": {"lung"}} and q.exact
    none = (obs.tissue == "lung") & (obs.tissue == "liver")
    assert none.hints == {"tissue": frozenset()} and none.exact


def test_and_keeps_one_set_of_a_multi_valued_column() -> None:
    """has_any(A) & has_any(B) is not has_any(A & B): keep one, and it is not exact."""
    q = obs.disease.has_any(["a", "b", "c"]) & obs.disease.has_any(["d"])
    assert q.hints == {"disease": {"d"}} and not q.exact
    every = obs.disease.has_all(["a", "b"])
    assert every.hints == {"disease": {"a"}} and not every.exact
    assert obs.disease.has_all(["a"]).exact


def test_or_and_not_keep_the_columns_and_drop_the_hints() -> None:
    for q in [
        (obs.tissue == "lung") | (obs.cell_type == "neutrophil"),
        ~(obs.tissue == "lung"),
        obs.tissue != "lung",
        ~obs.is_primary_data,
    ]:
        assert not q.hints and not q.exact
        assert q.columns
    q = (obs.tissue == "lung") | (obs.cell_type == "neutrophil")
    assert q.columns == {"tissue", "cell_type"}


def test_only_catalog_columns_give_hints() -> None:
    assert (obs.dataset_id == "d").hints == {"dataset_id": {"d"}}
    assert obs.is_primary_data._truth().hints == {"is_primary_data": {True}}
    assert (obs.is_primary_data == False).hints == {"is_primary_data": {False}}
    assert not (obs.donor_id == "x").hints
    assert not (obs.donor_id == "x").exact
    assert (obs["tissue"] == "lung").hints == {"tissue": {"lung"}}
    assert not obs.cell_type.str.startswith("CD8").hints


def test_where_builds_the_mapping() -> None:
    q = hs.where({"tissue": ["lung", "liver"], "disease": "flu", "dataset_id": "d"})
    assert q.exact
    assert q.hints == {
        "tissue": {"lung", "liver"},
        "disease": {"flu"},
        "dataset_id": {"d"},
    }
    assert as_query(hs.where({})).expression is None
    with pytest.raises(ValueError, match="not a catalog column"):
        hs.where({"donor_id": "x"})


def test_an_organism_is_the_query_of_every_cell() -> None:
    every = as_query(hs)
    assert every.expression is None and every.exact and not every.columns
    q = obs.tissue == "lung"
    assert (hs.query() & q) is q
    assert (q & hs.query()) is q


def test_queries_of_two_organisms_do_not_combine() -> None:
    mouse = census2_sdk.Organism("mus_musculus").obs
    with pytest.raises(ValueError, match="do not combine"):
        (obs.tissue == "lung") & (mouse.tissue == "lung")
    with pytest.raises(ValueError, match="no organism"):
        census2_sdk.Organism("danio_rerio")


def test_a_query_has_no_truth_value() -> None:
    q = obs.tissue == "lung"
    with pytest.raises(TypeError, match="no truth value"):
        bool(q)
    with pytest.raises(TypeError, match="no truth value"):
        q and q  # noqa: B018
    with pytest.raises(TypeError, match="no truth value"):
        bool(obs.tissue)


def test_refused_operators() -> None:
    with pytest.raises(TypeError, match="multi-valued"):
        obs.disease == "COVID-19"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="multi-valued"):
        obs.disease != "COVID-19"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="multi-valued"):
        obs.disease_ontology_term_id.isin(["MONDO:0100096"])  # type: ignore[deprecated]
    with pytest.raises(TypeError, match="holds labels"):
        obs.tissue == 3  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="holds numbers"):
        obs.nnz == "3"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="compare it with a bool"):
        obs.is_primary_data == 1  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="not one str"):
        obs.tissue.isin("lung")
    with pytest.raises(ValueError, match="RE2 refuses"):
        obs.cell_type.str.match("(unclosed")
    with pytest.raises(TypeError, match="not with int"):
        (obs.tissue == "lung") & cast("Any", 3)


def test_repr_is_the_pandas_form() -> None:
    q = ((obs.tissue == "lung") | (obs.tissue_general == "lung")) & ~(
        obs.cell_type == "macrophage"
    )
    assert repr(q) == (
        "<Query of homo_sapiens: ((obs.tissue == 'lung') | "
        "(obs.tissue_general == 'lung')) & ~(obs.cell_type == 'macrophage')>"
    )
    q = obs.cell_type.isin(["a"]) & (obs.nnz >= 500) & obs.is_primary_data
    assert (
        q.text == "obs.cell_type.isin(['a']) & (obs.nnz >= 500) & obs.is_primary_data"
    )
    assert obs["x"].isna().text == "obs['x'].isna()"
    assert obs.nnz.between(1, 5, inclusive="left").text == (
        "obs.nnz.between(1, 5, inclusive='left')"
    )


# --- evaluation with pyarrow ----------------------------------------------------

DISEASE = [
    "normal || flu",
    "flu",
    None,
    "normal",
    "cold",
    "a.b || c (d)",
    "normal ||  flu",
]
_DISEASE = pd.Categorical(DISEASE)
TABLE: dict[str, obs_read.FilterData] = {
    "disease": obs_read.Codes(
        np.asarray(_DISEASE.codes, dtype=np.int8), list(_DISEASE.categories), False
    ),
    "nnz": obs_read.Values(np.array([1, 5, 9, 3, 7, 2, 4]), None),
    "raw_variance_nnz": obs_read.Values(
        np.array([0.5, np.nan, 1.0, 2.0, np.nan, 3.0, 4.0]), None
    ),
    "flag": obs_read.Values(
        np.array([True, False, True, True, False, True, False]),
        np.array([False, False, True, False, False, False, False]),
    ),
}


def kept(q: Query) -> list[int]:
    return list(obs_read.evaluate(q.expression, TABLE, len(DISEASE)))


def test_has_any_finds_members() -> None:
    """A member matches; the catalog's rule strips the space around each member."""
    assert kept(obs.disease.has_any(["flu"])) == [0, 1, 6]
    assert kept(obs.disease.has_any(["normal"])) == [0, 3, 6]
    assert kept(obs.disease.has_all(["normal", "flu"])) == [0, 6]
    # Regular expression characters in a value are literal.
    assert kept(obs.disease.has_any(["a.b"])) == [5]
    assert kept(obs.disease.has_any(["c (d)"])) == [5]
    assert kept(obs.disease.has_any(["a?b", "c"])) == []
    assert members_pattern(["a|b"]) == r"(?:^| \|\| )\s*(?:a\|b)\s*(?: \|\| |$)"


def test_nulls_follow_three_valued_logic() -> None:
    """`~` of a null is null, so a negated term drops the null cells, unlike pandas."""
    flag = obs["flag"]
    assert kept(flag == True) == [0, 3, 5]
    assert kept(~flag) == [1, 4, 6]  # the null cell 2 is in neither
    assert kept(flag.isna()) == [2]
    assert kept(flag.notna()) == [0, 1, 3, 4, 5, 6]
    has = obs.disease.has_any(["flu"])
    assert kept(~has) == [3, 4, 5]  # the null disease of cell 2 drops
    isin = obs["disease"].isin(["flu"])
    assert kept(isin) == [1]
    assert kept(~isin) == [0, 3, 4, 5, 6]  # isin gives null for a null, as ==
    assert kept(obs["disease"] != "flu") == [0, 3, 4, 5, 6]
    # Kleene logic: null | true is true, and null & false is false.
    assert 2 in kept(flag | (obs.nnz > 8))
    assert kept(~(flag & (obs.nnz > 100))) == [0, 1, 2, 3, 4, 5, 6]


def test_numeric_terms_and_nan() -> None:
    assert kept(obs.nnz.between(2, 5)) == [1, 3, 5, 6]
    assert kept(obs.nnz.between(2, 5, inclusive="neither")) == [3, 6]
    assert kept(obs.raw_variance_nnz.isna()) == [1, 4]  # NaN, as in pandas
    assert kept(obs.raw_variance_nnz > 1.5) == [3, 5, 6]  # NaN compares false
    assert kept(~(obs.raw_variance_nnz > 1.5)) == [0, 1, 2, 4]


def test_a_missing_column_is_null() -> None:
    data = {"tissue": obs_read.Missing(3, "category")}
    q = obs.tissue == "lung"
    assert list(obs_read.evaluate(q.expression, data, 3)) == []
    assert list(obs_read.evaluate((~q).expression, data, 3)) == []
    assert list(obs_read.evaluate(obs.tissue.isna().expression, data, 3)) == [0, 1, 2]
    unknown = {"x": obs_read.Missing(2, "")}
    assert list(obs_read.evaluate((obs["x"] > 3).expression, unknown, 2)) == []


def test_string_methods() -> None:
    labels = ["CD8-positive T cell", "cd4 T cell", "B cell", "CD4-positive x"]
    data = {
        "cell_type": obs_read.Codes(
            np.array([0, 1, 2, 3], np.int8), sorted(labels), False
        )
    }
    by_label = {label: i for i, label in enumerate(sorted(labels))}

    def labels_of(q: Query) -> set[str]:
        found = obs_read.evaluate(q.expression, data, 4)
        return {sorted(labels)[int(i)] for i in found} & set(by_label)

    ct = obs.cell_type
    assert labels_of(ct.str.startswith("CD")) == {
        "CD8-positive T cell",
        "CD4-positive x",
    }
    assert labels_of(ct.str.startswith("cd", case=False)) == {
        "CD8-positive T cell",
        "cd4 T cell",
        "CD4-positive x",
    }
    assert labels_of(ct.str.endswith(("x", "B cell"))) == {"CD4-positive x", "B cell"}
    assert labels_of(ct.str.contains("T c", regex=False)) == {
        "CD8-positive T cell",
        "cd4 T cell",
    }
    assert labels_of(ct.str.match(r"CD(4|8)-positive")) == {
        "CD8-positive T cell",
        "CD4-positive x",
    }
    assert labels_of(ct.str.fullmatch("B cell")) == {"B cell"}


def test_expression_takes_a_pyarrow_expression() -> None:
    ratio = pc.divide(pc.field("nnz").cast(pa.float64()), pc.scalar(2.0)) > 2
    q = hs.expression(ratio, columns=["nnz"]) & obs.disease.has_any(["flu"])
    assert q.columns == {"nnz", "disease"}
    assert q.hints == {"disease": {"flu"}} and not q.exact
    assert kept(q) == [1]
    with pytest.raises(TypeError, match="pyarrow Expression"):
        hs.expression("nnz > 3", columns=["nnz"])


# --- queries over the test stores -------------------------------------------------


def brute_force(corpus: dict[str, Path], test: Any) -> set[str]:
    """The cell keys of `test`, a function of one store's whole obs, in pandas."""
    keys: set[str] = set()
    for dataset_id in STORES:
        adata = anndata.read_zarr(corpus["prefix"] / f"v-{dataset_id}.zarr")
        frame = cast("pd.DataFrame", adata.obs)
        mask = np.asarray(test(frame), dtype=bool)
        joinids = frame["observation_joinid"].astype(str)[mask]
        keys |= {f"{dataset_id}_{j}" for j in joinids}
    return keys


def _members(series: pd.Series) -> pd.Series:
    return series.astype(str).map(lambda v: {m.strip() for m in v.split(" || ")})


CASES: list[tuple[str, Any, Any]] = [
    ("isin", lambda o: o.tissue.isin(["lung", "kidney"]), lambda f: f["tissue"].isin(["lung", "kidney"])),
    ("ne", lambda o: o.tissue != "lung", lambda f: f["tissue"] != "lung"),
    ("or", lambda o: (o.tissue == "liver") | o.disease.has_any(["flu"]),
     lambda f: (f["tissue"] == "liver") | _members(f["disease"]).map(lambda m: "flu" in m)),
    ("not and", lambda o: ~((o.tissue == "lung") & o.disease.has_any(["normal"])),
     lambda f: ~((f["tissue"] == "lung") & _members(f["disease"]).map(lambda m: "normal" in m))),
    ("has_all", lambda o: o.disease.has_all(["normal", "flu"]),
     lambda f: _members(f["disease"]).map(lambda m: {"normal", "flu"} <= m)),
    ("primary", lambda o: o.is_primary_data & (o.dataset_id == "ds-b"),
     lambda f: f["is_primary_data"] & (f["dataset_id"] == "ds-b")),
    ("prefix", lambda o: o.tissue.str.startswith("l"), lambda f: f["tissue"].astype(str).str.startswith("l")),
]  # fmt: skip


@pytest.mark.parametrize(("name", "build", "test"), CASES, ids=[c[0] for c in CASES])
def test_queries_match_pandas(
    corpus: dict[str, Path], name: str, build: Any, test: Any
) -> None:
    """The test stores hold no null, so the pandas rule and SQL's rule agree."""
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata(build(census.homo_sapiens.obs))
    got = set() if result is None else set(result.obs_names)
    assert got == brute_force(corpus, test)


MAPPINGS: list[dict[str, Any]] = [
    {"tissue": "lung"},
    {"disease": "flu"},
    {"tissue": ["lung", "kidney"], "disease": "normal"},
    {"dataset_id": "ds-b", "tissue": "lung"},
    {"tissue": ["lung", "liver"]},
]


@pytest.mark.parametrize("mapping", MAPPINGS)
def test_mapping_equals_the_builder(
    corpus: dict[str, Path], mapping: dict[str, Any]
) -> None:
    """`where` and the builder form give one result, with the values of pandas."""
    census = census2_sdk.open_census([corpus["catalog"]])
    o = census.homo_sapiens.obs
    terms = [
        o.disease.has_any([v] if isinstance(v, str) else v)
        if column == "disease"
        else getattr(o, column).isin([v] if isinstance(v, str) else v)
        for column, v in mapping.items()
    ]
    built = terms[0]
    for term in terms[1:]:
        built = built & term
    left = census.get_anndata(census.homo_sapiens.where(mapping))
    right = census.get_anndata(built & (o.tissue.notna()))  # evaluated, not whole
    assert left is not None and right is not None
    assert list(left.obs_names) == list(right.obs_names)
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", left.obs), cast("pd.DataFrame", right.obs)
    )
    assert np.array_equal(X(left).toarray(), X(right).toarray())
    want = expected(corpus, mapping)
    assert set(left.obs_names) == set(want.obs_names)


def test_whole_stores_read_no_filter_column(corpus: dict[str, Path]) -> None:
    """A store whose value sets lie within the hints keeps all of its cells."""
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    rows = census.candidates(hs.where({"tissue": ["lung", "liver"]}))
    assert dict(zip(rows["dataset_id"], rows["whole"], strict=True)) == {
        "ds-a": True,
        "ds-b": True,
    }
    result = census.get_anndata(hs.where({"tissue": ["lung", "liver"]}))
    assert result is not None
    assert result.uns["census2"]["stores_whole"] == 2
    assert result.n_obs == 30 + 25
    # A column without a value set (nnz), or that the catalog does not record
    # (tissue_general), is never whole. The catalog's is_primary_data summary is a
    # value set: every cell of the test stores is primary.
    for q in [
        hs.where({"tissue": ["lung", "liver"]}) & (hs.obs.nnz > 0),
        hs.obs.tissue_general == "lung",
    ]:
        assert not census.candidates(q)["whole"].any()
    primary = hs.where({"tissue": ["lung", "liver"]}) & hs.obs.is_primary_data
    assert census.candidates(primary)["whole"].all()
    assert census.candidates(~hs.obs.is_primary_data).empty
    assert census.candidates(hs)["whole"].all()


def test_a_store_that_selects_no_cell_reads_nothing_more(
    corpus: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    read: list[str] = []
    original = obs_read.read_columns

    def spy(group: Any, names: Any, listed: Any) -> Any:
        read.extend(names)
        return original(group, names, listed)

    monkeypatch.setattr(obs_read, "read_columns", spy)
    # The catalog keeps the two lung stores. The stores hold no nnz column, so each
    # cell holds a null, and no cell matches.
    q = (hs.obs.tissue == "lung") & (hs.obs.nnz > 0)
    assert census.get_anndata(q) is None
    assert read == ["tissue", "tissue"]


def test_get_obs_reads_no_x(
    corpus: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    monkeypatch.setattr(census2_sdk, "_rows", None)
    frame = census.get_obs(hs.where({"tissue": "lung"}), columns=["dataset_id"])
    full = census.get_anndata
    monkeypatch.undo()
    result = full(hs.where({"tissue": "lung"}))
    assert result is not None
    assert list(frame.columns) == ["dataset_id", "observation_joinid"]
    assert list(frame.index) == list(result.obs_names)
    assert census.get_obs(hs.where({"tissue": "brain"})).empty


def test_categories_are_the_labels_in_use(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    result = census.get_anndata(census.homo_sapiens.where({"tissue": "lung"}))
    assert result is not None
    assert list(result.obs["tissue"].cat.categories) == ["lung"]


def test_count_cells_takes_only_an_exact_query(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    o = hs.obs
    for q in [o.tissue != "lung", (o.tissue == "lung") | (o.tissue == "kidney")]:
        with pytest.raises(ValueError, match="count_cells sums"):
            census.count_cells(q)
    with pytest.raises(ValueError, match="counts no"):
        census.count_cells(hs.where({"development_stage": "adult"}))
    q = o.disease.has_any(["flu"]) & o.is_primary_data
    result = census.get_anndata(q)
    assert result is not None
    assert census.count_cells(q) == result.n_obs


def test_random_conjunctions_count_as_they_read(corpus: dict[str, Path]) -> None:
    """For exact queries, the count, the whole proof, and the read agree."""
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    rng = random.Random(0)
    tissues = ["lung", "liver", "kidney", "brain"]
    diseases = ["normal", "flu"]
    for _ in range(12):
        mapping: dict[str, Any] = {}
        if rng.random() < 0.8:
            mapping["tissue"] = rng.sample(tissues, rng.randint(1, 3))
        if rng.random() < 0.5:
            mapping["disease"] = rng.sample(diseases, rng.randint(1, 2))
        q = hs.where(mapping)
        result = census.get_anndata(q)
        n = 0 if result is None else result.n_obs
        assert census.count_cells(q) == n
        assert n == len(expected(corpus, mapping)) if mapping else n == 75


def test_store_without_a_listed_category_count(tmp_path: Path) -> None:
    """A store whose categories differ in number from the row reads its own."""
    adata = make_store("ds-z", 6, ["ENSG00000000000"], ["lung", "liver"], ["normal"], 1)
    adata.obs["tissue"] = adata.obs["tissue"].cat.add_categories(["brain"])
    census = census2_sdk.open_census([one_store_catalog(tmp_path, adata, "ds-z")])
    hs = census.homo_sapiens
    q = hs.obs.tissue.isin(["lung"]) & hs.obs.is_primary_data
    result = census.get_anndata(q)
    want = int((adata.obs["tissue"] == "lung").sum())
    assert (0 if result is None else result.n_obs) == want


def test_every_combination_of_two_terms(corpus: dict[str, Path]) -> None:
    """`&`, `|`, and `~` of two terms equal pandas on the test stores."""
    census = census2_sdk.open_census([corpus["catalog"]])
    o = census.homo_sapiens.obs
    terms = {
        "lung": (o.tissue == "lung", lambda f: f["tissue"] == "lung"),
        "flu": (
            o.disease.has_any(["flu"]),
            lambda f: _members(f["disease"]).map(lambda m: "flu" in m),
        ),
    }
    for (a, (qa, fa)), (b, (qb, fb)) in itertools.product(terms.items(), repeat=2):
        for op in ("&", "|"):
            q = qa & ~qb if op == "&" else ~qa | qb

            def test(f: pd.DataFrame, fa: Any = fa, fb: Any = fb, op: str = op) -> Any:
                left, right = np.asarray(fa(f)), np.asarray(fb(f))
                return left & ~right if op == "&" else ~left | right

            result = census.get_anndata(q)
            got = set() if result is None else set(result.obs_names)
            assert got == brute_force(corpus, test), (a, op, b)


# --- the catalog test by pyarrow on the value sets --------------------------------


def scan(census: census2_sdk.Census, q: Query) -> set[str]:
    """The cell keys of `q` from every store of the organism, with no catalog test."""
    keys: set[str] = set()
    rows = census.catalog[census.catalog["organism"] == q.organism]
    for row in cast("list[dict[str, Any]]", rows.to_dict("records")):
        group = open_store(str(row["store_uri"]))
        obs = cast("Any", group["obs"])
        n = census2_sdk._n_obs(group)
        held = set(obs.attrs["column-order"])
        data: dict[str, obs_read.FilterData] = dict(
            obs_read.read_columns(obs, [c for c in q.columns if c in held], {})
        )
        for c in q.columns - held:
            data[c] = obs_read.Missing(n, KINDS.get(c, ""))
        kept = obs_read.evaluate(q.expression, data, n)
        joinids = np.asarray(obs_read.take_rows(obs["observation_joinid"], kept))
        keys |= {f"{row['dataset_id']}_{j}" for j in joinids}
    return keys


def random_query(o: Any, rng: random.Random) -> Query:
    terms = [
        o.tissue == "lung",
        o.tissue.isin(["liver", "kidney"]),
        o.tissue.str.startswith("l"),
        o.disease.has_any(["flu"]),
        o.disease.has_all(["normal", "flu"]),
        as_query(o.is_primary_data),
        o.nnz > 3,  # no store holds nnz: unknown in the catalog, null in the stores
        o.dataset_id == "ds-b",
        o.donor_id == "d",  # no value set, and absent from the stores
    ]
    q: Query = rng.choice(terms)
    for _ in range(rng.randint(0, 3)):
        other: Query = rng.choice(terms)
        if rng.random() < 0.3:
            other = ~other
        q = q & other if rng.random() < 0.5 else q | other
    return ~q if rng.random() < 0.3 else q


def test_the_catalog_test_keeps_every_store_with_a_cell(
    corpus: dict[str, Path],
) -> None:
    """For random queries, the result through the catalog equals a scan of all stores.

    A whole store keeps all of its cells, so a wrong whole verdict shows too.
    """
    census = census2_sdk.open_census([corpus["catalog"]])
    o = census.homo_sapiens.obs
    rng = random.Random(1)
    for _ in range(40):
        q = random_query(o, rng)
        result = census.get_anndata(q)
        got = set() if result is None else set(result.obs_names)
        assert got == scan(census, q), q


def test_the_value_sets_prune_or_not_and_string_methods(
    corpus: dict[str, Path],
) -> None:
    census = census2_sdk.open_census([corpus["catalog"]])
    o = census.homo_sapiens.obs

    def ids(q: Query) -> list[str]:
        return list(census.candidates(q)["dataset_id"])

    assert ids((o.tissue == "kidney") | o.disease.has_any(["flu"])) == ["ds-b", "ds-c"]
    assert ids(~o.tissue.isin(["lung", "liver"])) == ["ds-c"]
    assert ids(o.tissue.str.startswith("kid")) == ["ds-c"]
    assert ids(o.disease.has_all(["normal", "flu"])) == ["ds-b"]
    assert ids((o.nnz > 3) & (o.tissue == "kidney")) == ["ds-c"]
    rows = census.candidates(o.tissue.str.startswith("kid") | (o.nnz > 3))
    assert len(rows) == 3  # null | false is null, so each row stays
    # True for every combination of ds-c, but nnz is unknown: not whole.
    assert not rows["whole"].any()


def test_the_value_sets_give_way_to_the_hints(
    corpus: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Above the size limit, or for a query that is not null-safe, the hints decide."""
    census = census2_sdk.open_census([corpus["catalog"]])
    hs = census.homo_sapiens
    prefix = hs.obs.tissue.str.startswith("kid")
    assert len(census.candidates(prefix)) == 1
    assert len(census.candidates(prefix | hs.obs.tissue.isna())) == 3
    raw = hs.expression(pc.field("tissue") == "kidney", columns=["tissue"])
    assert len(census.candidates(raw)) == 3
    monkeypatch.setattr(value_sets, "MAX_COMBINATIONS", 0)
    assert len(census.candidates(prefix)) == 3
    row: dict[str, Any] = {"dataset_id": "d", "tissue": [], "n_obs": 4}
    assert value_sets.value_set(row, "tissue") is None  # the cells hold nulls
    assert value_sets.value_set({**row, "n_obs": 0}, "tissue") == []
