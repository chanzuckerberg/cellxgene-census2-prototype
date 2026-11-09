"""Tests for the builder path over small local stores. No network, no AWS.

The corpus is `obs_stores`: four stores with every column of the store format's
`obs`, a catalog, and a second catalog in which three of them are members of one
consolidated store. The tests cover the resolution of a filter, the catalog test, the
Arrow compiler, the pandas null rule, the equality with `main`'s path for mapping
filters, and the consolidated reads.
"""

import random
import warnings
from pathlib import Path
from typing import Any, cast

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pytest
import scipy.sparse as sp
import zarr
from anndata.io import read_elem
from obs_stores import SPECS, make_corpus

import census2_sdk
from census2_sdk import _mask, obs_arrow, obs_plan
from census2_sdk.catalog import MULTIVALUED, canonicalize
from census2_sdk.obs_arrow import Codes, Values
from census2_sdk.obs_filter import Const, FilterError, from_mapping, normalize

ORG = "homo_sapiens"


@pytest.fixture(scope="module")
def corpus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return make_corpus(tmp_path_factory.mktemp("obs"))


@pytest.fixture(scope="module")
def census(corpus: dict[str, Path]) -> census2_sdk.Census:
    return census2_sdk.open_census([corpus["catalog"]])


def store_obs(corpus: dict[str, Path], dataset_id: str) -> pd.DataFrame:
    group = zarr.open_group(corpus["prefix"] / f"v-{dataset_id}.zarr", mode="r")
    return cast("pd.DataFrame", read_elem(group["obs"]))


# Resolution.


def test_an_unknown_label_fails_with_a_suggestion(census: census2_sdk.Census) -> None:
    hs = census.homo_sapiens
    with pytest.raises(FilterError, match="did you mean 'lung'"):
        hs.plan(hs.obs.tissue == "lugn")
    with pytest.raises(FilterError, match="did you mean 'neutrophil'"):
        hs.plan(hs.obs.cell_type.isin(["Neutrophil"]))
    # has_all with one unknown member can match no cell.
    with pytest.raises(FilterError, match="matches no cell"):
        hs.plan(hs.obs.disease.has_all(["COVID-19", "influenza"]))
    with pytest.raises(FilterError, match="dataset_id"):
        hs.plan(hs.obs.dataset_id == "ds-9")


def test_a_near_unknown_label_beside_known_labels_warns(
    census: census2_sdk.Census,
) -> None:
    hs = census.homo_sapiens
    with pytest.warns(UserWarning, match="did you mean 'liver'"):
        plan = hs.plan(hs.obs.tissue.isin(["lung", "livr"]))
    assert {c.row["dataset_id"] for c in plan.candidates} == {"ds-1", "ds-2"}
    # A label with no near known label passes in silence: it matches no cell.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        hs.plan(hs.obs.cell_type.isin(["neutrophil", "Purkinje cell of cerebellum"]))


def test_the_escape_hatch(corpus: dict[str, Path], census: census2_sdk.Census) -> None:
    hs = census.homo_sapiens
    # A column of the schema keeps the rules of its kind.
    with pytest.raises(FilterError, match="multi-valued"):
        hs.plan(hs.obs["disease"] == "normal")
    with pytest.raises(FilterError, match="holds numbers"):
        hs.plan(hs.obs["nnz"] == "a")
    with pytest.raises(FilterError, match="boolean"):
        hs.plan(hs.obs["is_primary_data"] == "a")
    assert hs.count(hs.obs["tissue"] == "kidney") == SPECS["ds-3"]["n"]
    # A column that the schema does not name waits for the stores.
    with pytest.raises(FilterError, match="did you mean 'author_label'"):
        hs.count(hs.obs["author_labl"] == "a")
    with pytest.warns(UserWarning, match="3 of the stores read have no obs column"):
        n = hs.count(hs.obs["author_label"] == "a")
    labels = store_obs(corpus, "ds-4")
    assert n == int((labels["author_label"] == "a").sum())
    with pytest.raises(FilterError, match="label values"):
        hs.count(hs.obs["author_label"] > 1)


def test_a_filter_of_another_organism_is_refused(census: census2_sdk.Census) -> None:
    with pytest.raises(FilterError, match="census.homo_sapiens.obs"):
        census.homo_sapiens.plan(census.mus_musculus.obs.tissue == "lung")


# The catalog test.


def candidates(census: census2_sdk.Census, f: Any) -> dict[str, str]:
    plan = census.homo_sapiens.plan(f)
    return {str(c.row["dataset_id"]): repr(c.residual) for c in plan.candidates}


def test_one_term_is_exact(census: census2_sdk.Census) -> None:
    obs = census.homo_sapiens.obs
    assert candidates(census, obs.tissue == "kidney") == {"ds-3": "True"}
    assert candidates(census, obs.tissue == "lung") == {
        "ds-1": "obs.tissue == 'lung'",
        "ds-2": "True",
    }
    assert candidates(census, obs.is_primary_data) == {
        "ds-1": "True",
        "ds-2": "obs.is_primary_data",
        "ds-4": "obs.is_primary_data",
    }
    assert candidates(census, obs.dataset_id.isin(["ds-2", "ds-4"])) == {
        "ds-2": "True",
        "ds-4": "True",
    }
    assert candidates(census, obs.disease.has_any(["COVID-19"])) == {
        "ds-2": "obs.disease.has_any(['COVID-19'])"
    }
    assert candidates(census, obs.cell_type.str.endswith("T cell")) == {
        "ds-1": "obs.cell_type.str.endswith('T cell')",
        "ds-4": "obs.cell_type.str.endswith('T cell')",
    }


def test_negation_meets_the_pandas_null_rule(census: census2_sdk.Census) -> None:
    obs = census.homo_sapiens.obs
    # The value set of ds-3 is {kidney}, and the column holds no null, so no cell
    # of ds-3 is outside kidney.
    kept = candidates(census, ~(obs.tissue == "kidney"))
    assert "ds-3" not in kept and kept["ds-4"] == "True"
    # A column that the catalog does not record proves nothing, so a negated term
    # on it keeps every dataset: a null there matches the negation.
    assert len(candidates(census, ~(obs["author_label"] == "a"))) == 4
    assert len(candidates(census, obs["author_label"].isna())) == 4
    # isna of a column without nulls excludes every dataset.
    assert candidates(census, obs.tissue.isna()) == {}
    assert len(candidates(census, obs.tissue.notna())) == 4


def test_a_conjunction_is_conservative(census: census2_sdk.Census) -> None:
    obs = census.homo_sapiens.obs
    f = (obs.tissue == "lung") & (obs.cell_type == "alveolar macrophage")
    # ds-2 holds lung and no alveolar macrophage: excluded by its second term.
    assert set(candidates(census, f)) == {"ds-1"}
    # A numeric term proves nothing.
    assert len(candidates(census, obs.nnz > 10**9)) == 4
    g = (obs.tissue == "kidney") | (obs.nnz > 10**9)
    assert candidates(census, g)["ds-3"] == "True"


# The Arrow compiler.


def test_codes_become_a_dictionary_array_with_nulls() -> None:
    codes = np.array([0, 1, -1, 1], dtype=np.int8)
    arr = Codes(codes, np.array(["a", "b"], dtype=object), False).arrow()
    assert arr.type == pa.dictionary(pa.int8(), pa.string())
    assert arr.to_pylist() == ["a", "b", None, "b"]
    clean = np.array([0, 1, 1], dtype=np.int8)
    arr = Codes(clean, np.array(["a", "b"], dtype=object), False).arrow()
    # Without a null, the indices wrap the codes without a copy.
    assert arr.indices.buffers()[1].address == clean.ctypes.data
    values = np.arange(5, dtype=np.int64)
    assert Values(values, None).arrow().buffers()[1].address == values.ctypes.data


def reference(tree: Any, data: dict[str, Any]) -> np.ndarray:
    """The rows by pyarrow on decoded strings, with Table.filter: the reference."""
    table = pa.table(
        {
            k: (v.arrow().cast(pa.string()) if isinstance(v, Codes) else v.arrow())
            for k, v in data.items()
        }
    ).append_column("__row__", pa.array(np.arange(len(next(iter(data.values()))))))
    expression = obs_arrow.to_expression(tree)
    return np.asarray(table.filter(expression)["__row__"].to_numpy())


@pytest.mark.parametrize("seed", range(40))
def test_the_bound_filter_equals_pyarrow_on_strings(seed: int) -> None:
    rng = np.random.default_rng(seed)
    n = 200
    cats = np.array(["lung", "liver", "lung epithelium", "kidney"], dtype=object)
    codes = rng.integers(-1, len(cats), n).astype(np.int8)  # -1 is a null
    nnz = rng.integers(0, 100, n).astype(np.int64)
    flag = rng.random(n) < 0.5
    flag_mask = rng.random(n) < 0.2
    data: dict[str, Any] = {
        "tissue": Codes(codes, cats, False),
        "nnz": Values(nnz, None),
        "flag": Values(flag, flag_mask),
    }
    obs = census2_sdk.obs_columns.Obs(ORG)
    tissue, number, truth = obs["tissue"], obs["nnz"], obs["flag"]
    terms: list[Any] = [
        tissue == "lung",
        tissue.isin(["liver", "kidney"]),
        tissue.str.startswith("lung"),
        tissue.str.contains("EPI", case=False),
        tissue.isna(),
        number >= int(rng.integers(0, 100)),
        number.between(10, 60),
        truth,
        truth.isna(),
    ]
    random.seed(seed)
    tree: Any = random.choice(terms)
    for _ in range(int(rng.integers(1, 5))):
        other = random.choice(terms)
        if rng.random() < 0.3:
            other = ~other
        tree = tree & other if rng.random() < 0.5 else tree | other
    if rng.random() < 0.5:
        tree = ~tree
    tree = normalize(tree)
    got = obs_arrow.evaluate(tree, data, n)
    assert list(got) == list(reference(tree, data))


def test_pyarrow_decodes_a_dictionary_for_equality() -> None:
    """The finding of the spike: pyarrow 26 has no dictionary kernel for these."""
    d = pa.DictionaryArray.from_arrays(pa.array([0, 1, 0], pa.int8()), ["a", "b"])
    with pytest.raises(pa.ArrowNotImplementedError):
        pa.compute.starts_with(d, "a")
    assert pa.compute.equal(d, "a").to_pylist() == [True, False, True]


# The null rule on a store.


def test_the_null_rule_matches_pandas(
    corpus: dict[str, Path], census: census2_sdk.Census
) -> None:
    hs = census.homo_sapiens
    obs = hs.obs
    frame = store_obs(corpus, "ds-4")
    label = frame["author_label"]
    assert label.isna().any()
    cases = [
        (~(obs["author_label"] == "a"), ~(label == "a")),
        (obs["author_label"] != "a", label != "a"),
        (obs["author_label"].isna(), label.isna()),
        (obs["author_label"].notna(), label.notna()),
        (~obs["author_label"].isin(["a", "b"]), ~label.isin(["a", "b"])),
        (
            ~((obs["author_label"] == "a") | (obs.tissue == "lung")),
            ~((label == "a") | (frame["tissue"] == "lung")),
        ),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        for f, want in cases:
            got = hs.get_obs(
                f & (obs.dataset_id == "ds-4"), columns=["observation_joinid"]
            )
            expected = frame["observation_joinid"][want.to_numpy()].astype(str)
            assert list(got["observation_joinid"].astype(str)) == list(expected), f


# Equality with main's path.


def same_values(got: anndata.AnnData, want: anndata.AnnData) -> None:
    assert list(got.obs_names) == list(want.obs_names)
    assert list(got.obs.columns) == list(want.obs.columns)
    for name in want.obs.columns:
        g, w = got.obs[name], want.obs[name]
        if isinstance(w.dtype, pd.CategoricalDtype):
            assert isinstance(g.dtype, pd.CategoricalDtype), name
            assert list(g.astype(str)) == list(w.astype(str)), name
            # The builder keeps the labels that the result uses.
            assert set(g.cat.categories) == set(w.dropna().astype(str)), name
        else:
            assert g.dtype == w.dtype, name
            assert list(g) == list(w), name
    pd.testing.assert_frame_equal(
        cast("pd.DataFrame", got.var), cast("pd.DataFrame", want.var)
    )
    gx, wx = sp.csr_matrix(cast("Any", got.X)), sp.csr_matrix(cast("Any", want.X))
    assert gx.dtype == wx.dtype
    assert np.array_equal(gx.toarray(), wx.toarray())
    for layer in list(want.layers.keys()):
        assert np.array_equal(
            sp.csr_matrix(got.layers[layer]).toarray(),
            sp.csr_matrix(want.layers[layer]).toarray(),
        )


def random_mapping(rng: random.Random, known: dict[str, list[str]]) -> dict[str, Any]:
    columns = rng.sample(sorted(known), rng.randint(1, 3))
    out: dict[str, Any] = {}
    for column in columns:
        values = rng.sample(known[column], rng.randint(1, min(3, len(known[column]))))
        out[column] = values[0] if len(values) == 1 and rng.random() < 0.5 else values
    return out


def known_values(census: census2_sdk.Census) -> dict[str, list[str]]:
    sets = obs_plan.value_sets(census.catalog)
    columns = ["tissue", "cell_type", "disease", "sex", "dataset_id", "tissue_general"]
    return {c: sorted(sets[c]) for c in columns}


@pytest.mark.parametrize("seed", range(30))
def test_the_store_evaluation_equals_main_mask(
    corpus: dict[str, Path], census: census2_sdk.Census, seed: int
) -> None:
    rng = random.Random(seed)
    mapping = random_mapping(rng, known_values(census))
    tree = normalize(from_mapping(census.homo_sapiens.obs, mapping, MULTIVALUED))
    canonical = canonicalize(mapping)
    for dataset_id in SPECS:
        group = zarr.open_group(corpus["prefix"] / f"v-{dataset_id}.zarr", mode="r")
        obs_group = group["obs"]
        assert isinstance(obs_group, zarr.Group)
        frame = cast("pd.DataFrame", read_elem(obs_group))
        data = obs_arrow.read_columns(obs_group, sorted(mapping), 0, len(frame))
        got = obs_arrow.evaluate(tree, data, len(frame))
        # `_mask` leaves dataset_id to the catalog; the store holds one dataset.
        want = np.flatnonzero(_mask(frame, canonical))
        if "dataset_id" in canonical and dataset_id not in canonical["dataset_id"]:
            want = want[:0]
        assert list(got) == list(want), (mapping, dataset_id)


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("key", ["catalog", "catalog_c"])
def test_get_anndata_equals_main(corpus: dict[str, Path], seed: int, key: str) -> None:
    census = census2_sdk.open_census([corpus[key]])
    rng = random.Random(seed)
    mapping = random_mapping(rng, known_values(census))
    tree = from_mapping(census.homo_sapiens.obs, mapping, MULTIVALUED)
    want = census.get_anndata(ORG, mapping, layers=["normalized"], max_workers=2)
    got = census.homo_sapiens.get_anndata(tree, layers=["normalized"], max_workers=2)
    if want is None:
        assert got is None
        return
    assert got is not None
    same_values(got, want)


def test_a_filter_beyond_the_mapping(corpus: dict[str, Path]) -> None:
    """The same cells as pandas over the whole `obs` of every store."""
    census = census2_sdk.open_census([corpus["catalog_c"]])
    hs = census.homo_sapiens
    obs = hs.obs
    f = (
        ((obs.cell_type == "neutrophil") | (obs.tissue_general == "kidney"))
        & ~obs.disease.has_any(["COVID-19"])
        & (obs.nnz.between(100, 2500) | obs.is_primary_data)
    )
    got = hs.get_anndata(f, max_workers=3)
    assert got is not None
    parts = []
    for dataset_id in SPECS:
        frame = store_obs(corpus, dataset_id)
        covid = frame["disease"].astype(str).str.split(" || ", regex=False)
        mask = (
            (
                (frame["cell_type"] == "neutrophil")
                | (frame["tissue_general"] == "kidney")
            )
            & ~covid.map(lambda m: "COVID-19" in m)
            & (frame["nnz"].between(100, 2500) | frame["is_primary_data"])
        )
        parts += [f"{dataset_id}_{j}" for j in frame["observation_joinid"][mask]]
    assert list(got.obs_names) == parts
    assert hs.count(f) == len(parts)


# Consolidated reads.


def test_adjacent_members_read_once(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog_c"]])
    hs = census.homo_sapiens
    hs.count(hs.obs.cell_type == "neutrophil")
    # ds-1 and ds-2 are adjacent members of one store; ds-4 is alone.
    assert hs.last_query == {"candidates": 3, "stores_opened": 2, "reads": 2}
    flat = census2_sdk.open_census([corpus["catalog_c"]], consolidated=False)
    flat.homo_sapiens.count(flat.homo_sapiens.obs.cell_type == "neutrophil")
    assert flat.homo_sapiens.last_query["reads"] == 3


def test_a_whole_dataset_counts_with_no_store_open(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog_c"]])
    hs = census.homo_sapiens
    n = hs.count(hs.obs.dataset_id.isin(["ds-2", "ds-4"]) | (hs.obs.tissue == "kidney"))
    assert n == SPECS["ds-2"]["n"] + SPECS["ds-3"]["n"] + SPECS["ds-4"]["n"]
    assert hs.last_query["stores_opened"] == 0


@pytest.mark.parametrize("key", ["catalog", "catalog_c"])
def test_count_by_equals_pandas(corpus: dict[str, Path], key: str) -> None:
    census = census2_sdk.open_census([corpus[key]])
    hs = census.homo_sapiens
    f = hs.obs.cell_type.isin(["neutrophil", "monocyte"])
    got = hs.count_by(["tissue", "sex"], f, max_workers=2)
    frames = [store_obs(corpus, d) for d in SPECS]
    frame = pd.concat([f[["tissue", "sex", "cell_type"]].astype(str) for f in frames])
    frame = frame[frame["cell_type"].isin(["neutrophil", "monocyte"])]
    want = frame.groupby(["tissue", "sex"]).size()
    assert got.to_dict() == want.to_dict()
    assert hs.count(f) == int(want.sum())


def test_get_obs_columns_and_cell_key(corpus: dict[str, Path]) -> None:
    census = census2_sdk.open_census([corpus["catalog_c"]])
    hs = census.homo_sapiens
    frame = hs.get_obs(hs.obs.tissue == "lung", columns=["dataset_id", "cell_type"])
    assert list(frame.columns) == ["dataset_id", "cell_type"]
    assert isinstance(frame["cell_type"].dtype, pd.CategoricalDtype)
    assert frame.index.is_unique and frame.index[0].startswith("ds-1_ds-1-j")
    assert set(frame["dataset_id"]) == {"ds-1", "ds-2"}
    assert len(frame) == hs.count(hs.obs.tissue == "lung")
    assert hs.get_obs(hs.obs.tissue == "kidney", columns=[]).shape == (30, 0)
    assert census.homo_sapiens.plan(None).tree == Const(ORG, True)


@pytest.mark.parametrize("key", ["catalog", "catalog_c"])
def test_get_obs_with_a_read_that_keeps_no_cell(
    corpus: dict[str, Path], key: str
) -> None:
    census = census2_sdk.open_census([corpus[key]])
    hs = census.homo_sapiens
    f = (hs.obs.cell_type == "neutrophil") & (hs.obs.nnz > 2900)
    frame = hs.get_obs(f, max_workers=2)
    frames = [store_obs(corpus, d) for d in SPECS]
    keys = [
        f"{d}_{j}"
        for d, fr in zip(SPECS, frames, strict=True)
        for j in fr["observation_joinid"][
            (fr["cell_type"] == "neutrophil") & (fr["nnz"] > 2900)
        ]
    ]
    assert list(frame.index) == keys
    assert hs.last_query["candidates"] == 3
    # At least one candidate keeps no cell, and one keeps some.
    per_store = {
        d: sum(k.startswith(f"{d}_") for k in keys) for d in ("ds-1", "ds-2", "ds-4")
    }
    assert min(per_store.values()) == 0 < max(per_store.values()), per_store
    assert isinstance(frame["tissue"].dtype, pd.CategoricalDtype)
