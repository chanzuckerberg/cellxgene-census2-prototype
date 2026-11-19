"""The SDK adapter's filter translation and the driver's gene restriction and owners."""

import pathlib
from collections.abc import Callable, Sequence
from dataclasses import replace
from typing import Any, cast

import anndata as ad
import build_catalog
import numpy as np
import pandas as pd
import pytest
from harness.canonical import (
    Arrays,
    canonicalize_anndata,
    canonicalize_frame,
    frame_arrays,
    layer_arrays,
    presence_arrays,
    write_result,
)
from harness.comparator import Finding
from harness.compare_stacks import (
    SDK_ONLY_GENES,
    expected,
    owner,
    restrict_genes,
    run_one,
)
from harness.corpus import BoundQuery
from harness.sdk_adapter import NoCells, Unsupported, answer, clauses, translate
from scipy import sparse
from test_sdk import REFERENCE_IDS, STORES, make_store, reference_var, write_one

import census2_sdk
from census2_sdk import SourceConflict

PINS = {
    "disease_ontology_term_id": "MONDO:0005109",
    "disease_ontology_term_id_matches": ["MONDO:0005109 || MONDO:0005445"],
}


def test_translate_takes_equality_membership_and_conjunction() -> None:
    assert translate("cell_type_ontology_term_id == 'CL:0000897'", PINS) == {
        "cell_type_ontology_term_id": ["CL:0000897"]
    }
    assert translate(
        "tissue_general in ['brain', 'lung'] and sex == 'female'", PINS
    ) == {"tissue_general": ["brain", "lung"], "sex": ["female"]}
    both = "sex in ['female', 'male'] and sex == 'female'"
    assert translate(both, PINS) == {"sex": ["female"]}
    assert translate("sex == 'female' and sex == 'male'", PINS) == {"sex": []}


def test_translate_gives_a_multivalued_column_the_bare_term() -> None:
    joined = "disease_ontology_term_id in ['MONDO:0005109 || MONDO:0005445']"
    assert translate(joined, PINS) == {"disease_ontology_term_id": ["MONDO:0005109"]}
    with pytest.raises(Unsupported):
        translate("disease_ontology_term_id == 'MONDO:0001627 || MONDO:0004975'", PINS)


@pytest.mark.parametrize(
    "value_filter",
    [
        "nnz >= 500",
        "cell_type_ontology_term_id != 'unknown'",
        "is_primary_data == True and nnz >= 500",
    ],
)
def test_translate_refuses_what_the_sdk_cannot_express(value_filter: str) -> None:
    with pytest.raises(Unsupported):
        translate(value_filter, PINS)


def test_restrict_genes_keeps_var_and_layer_columns_aligned() -> None:
    var = pd.DataFrame(
        {"feature_id": ["g1", "g2", "g3"], "feature_length": [10, 20, 30]}
    )
    x = sparse.csr_matrix(np.array([[1, 0, 3], [0, 5, 6]], dtype=np.float32))
    arrays = frame_arrays("var", var)
    arrays.update(layer_arrays("raw", x))
    out = restrict_genes(arrays, {"g1", "g3"})
    assert list(out["var.feature_length"]) == [10, 30]
    rows, cols = (int(n) for n in out["X.raw.shape"])
    kept = sparse.csr_matrix(
        (out["X.raw.data"], out["X.raw.indices"], out["X.raw.indptr"]),
        shape=(rows, cols),
    )
    assert kept.toarray().tolist() == [[1, 3], [0, 6]]


def test_owner_names_the_task_that_removes_a_finding() -> None:
    assert (
        owner(Finding("column-set", "obs.raw_sum", "present on left only"), set())[0]
        == "4.5"
    )
    assert (
        owner(Finding("column-set", "obs.leiden", "present on right only"), set())[0]
        == "4.4"
    )
    assert (
        owner(Finding("column-set", "X.normalized", "present on left only"), set())[0]
        == "4.13"
    )
    tissue = Finding("strict-subset", "obs", "right is a strict subset")
    assert owner(tissue, {"tissue_general_ontology_term_id", "dataset_id"})[0] == "4.5"
    assert owner(Finding("value", "X.raw", "cell x gene"), set())[0] == "6.2"
    assert owner(Finding("value", "obs.donor_id", "identity x"), set())[0] == "6.2"


def test_clauses_split_outside_quotes_and_brackets() -> None:
    assert clauses("a == 'x' and b in ['p and q', 'r']") == [
        "a == 'x'",
        "b in ['p and q', 'r']",
    ]
    assert clauses("tissue == 'blood and lymph'") == ["tissue == 'blood and lymph'"]
    assert clauses(r"a == 'it\'s and more' and b == 'y'") == [
        r"a == 'it\'s and more'",
        "b == 'y'",
    ]


def test_translate_reads_quoted_values_and_booleans() -> None:
    assert translate("tissue == 'blood and lymph'", PINS) == {
        "tissue": ["blood and lymph"]
    }
    assert translate("is_primary_data == True", PINS) == {"is_primary_data": ["True"]}
    with pytest.raises(Unsupported):
        translate("tissue == blood", PINS)


QUERY = BoundQuery(
    key="cell_type_slice",
    accessor="get_anndata",
    modality="census_data",
    experiment="homo_sapiens",
    relied_on_by=("6",),
    obs_filter="cell_type_ontology_term_id == 'CL:0000897'",
    var_filter=None,
    layers=("raw",),
    obs_columns=(),
    label="cell_type_slice:census_data/homo_sapiens",
)


class Empty:
    def get_anndata(
        self, organism: str, obs_filter: object, layers: object = ()
    ) -> None:
        return None


def test_answer_lets_an_sdk_defect_fail_the_run() -> None:
    class Broken:
        def get_anndata(
            self, organism: str, obs_filter: object, layers: object = ()
        ) -> None:
            raise SourceConflict("the cell key is not unique across the result")

    with pytest.raises(SourceConflict):
        answer(QUERY, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
    unindexed = replace(QUERY, obs_filter="is_primary_data == True")
    with pytest.raises(Unsupported, match="not a catalog column"):
        answer(unindexed, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
    with pytest.raises(NoCells):
        answer(QUERY, cast("census2_sdk.Census", Empty()), ["d1"], PINS)


def _run_one(
    tmp_path: pathlib.Path,
    query: BoundQuery,
    census: object,
    census1: Arrays,
    clamp: Sequence[str] = ("d1",),
) -> dict[str, Any]:
    """`run_one` over a Census 1 result of `census1` and a clamp of `clamp`."""
    meta = {"clamp": list(clamp), "env": {"corpus_version": "1.1"}}
    write_result(tmp_path / query.label.replace("/", "__"), census1, meta)
    return run_one(
        query,
        tmp_path,
        cast("census2_sdk.Census", census),
        set(clamp),
        {"census_data/homo_sapiens": PINS},
        {"census_data/homo_sapiens": {"dataset_ids": list(clamp)}},
    )


@pytest.mark.parametrize(("census1_rows", "status"), [(2, "fail"), (0, "vacuous")])
def test_run_one_records_an_empty_sdk_result(
    tmp_path: pathlib.Path, census1_rows: int, status: str
) -> None:
    query = replace(QUERY, obs_filter="tissue_general == 'lung'")
    census1 = {"obs.dataset_id.offsets": np.arange(census1_rows + 1, dtype=np.uint64)}
    record = _run_one(tmp_path, query, Empty(), census1)
    assert record["status"] == status
    assert record["cells"] == {"census1": census1_rows, "sdk": 0}
    assert {"seconds", "notes", "findings"} <= record.keys()
    assert [f["owner"] for f in record["findings"]] == (["4.5"] if census1_rows else [])


def _adata(genes: list[str]) -> ad.AnnData:
    obs = pd.DataFrame(
        {"dataset_id": ["d1", "d1"], "observation_joinid": [0, 1]}, index=["0", "1"]
    )
    var = pd.DataFrame({"feature_id": genes}, index=genes)
    return ad.AnnData(sparse.csr_matrix(np.ones((2, len(genes)), np.float32)), obs, var)


@pytest.mark.parametrize(
    ("sdk_genes", "status"), [(["g1", "g2", "g3"], "pass"), (["g1"], "fail")]
)
def test_run_one_fails_when_the_sdk_lacks_a_census_1_gene(
    tmp_path: pathlib.Path, sdk_genes: list[str], status: str
) -> None:
    class Stores:
        def get_anndata(
            self, organism: str, obs_filter: object, layers: object = ()
        ) -> ad.AnnData:
            return _adata(sdk_genes)

    census1 = canonicalize_anndata(_adata(["g1", "g2"]))
    record = _run_one(tmp_path, QUERY, Stores(), census1)
    assert record["status"] == status


@pytest.mark.parametrize(
    ("layers", "asked"), [(("raw",), []), (("raw", "normalized"), ["normalized"])]
)
def test_answer_asks_the_sdk_for_the_query_layers(
    layers: tuple[str, ...], asked: list[str]
) -> None:
    """The adapter asks for every layer but `raw`, and the layer reaches the arrays."""
    seen: list[list[str]] = []

    class Stores:
        def get_anndata(
            self, organism: str, obs_filter: object, layers: Sequence[str] = ()
        ) -> ad.AnnData:
            seen.append(list(layers))
            adata = _adata(["g1", "g2"])
            for name in layers:
                adata.layers[name] = sparse.csr_matrix(np.full((2, 2), 0.5, np.float32))
            return adata

    query = replace(QUERY, layers=layers)
    arrays, _ = answer(query, cast("census2_sdk.Census", Stores()), ["d1"], PINS)
    assert seen == [asked]
    if asked:
        assert np.array_equal(arrays["X.normalized.data"], np.full(4, 0.5))
    else:
        assert "X.normalized.data" not in arrays


def test_an_expected_difference_holds_in_its_declared_direction() -> None:
    def why(finding: Finding) -> str | None:
        return expected(finding, {}, {})

    assert why(Finding("column-set", "var.feature_version", "present on right only"))
    assert not why(Finding("column-set", "var.feature_version", "present on left only"))
    assert why(Finding("column-set", "var.nnz", "present on left only"))
    assert not why(Finding("column-set", "var.nnz", "present on right only"))
    assert not why(Finding("column-set", "X.normalized", "present on left only"))


S, M, N = "raw_sum", "raw_mean_nnz", "nnz"


@pytest.mark.parametrize(
    ("census1_obs", "sdk_obs", "status"),
    [
        ({S: 16_777_216.0}, {S: 16_777_217.0}, "pass"),  # Census 1 rounds the sum
        ({S: 1_000.0}, {S: 1_001.0}, "fail"),  # below 2**24, float32 sums are exact
        ({S: 16_777_216.0}, {S: 16_800_000.0}, "fail"),  # too large for a rounding
        (
            {S: 2.0**24, M: 2.0**23, N: 2},
            {S: 2.0**24 + 1, M: 2.0**23 + 0.5, N: 2},
            "pass",
        ),
        ({S: 2.0**25, M: 2.0**24, N: 2}, {S: 2.0**25, M: 2.0**24 + 1, N: 2}, "fail"),
        (
            {S: 2.0**24, M: 2.0**23, N: 2},
            {M: 2.0**23 + 0.5, N: 2},
            "fail",
        ),  # no raw_sum
        # the SDK mean does not follow from the SDK sum
        (
            {S: 2.0**24, M: 2.0**23, N: 2},
            {S: 2.0**24 + 1, M: 2.0**23 - 0.5, N: 2},
            "fail",
        ),
        # a later cell differs below 2**24
        ({S: [2.0**24, 5.0]}, {S: [2.0**24 + 1, 6.0]}, "fail"),
        # a later mean differs at an equal sum
        (
            {S: [2.0**24, 10.0], M: [2.0**23, 5.0], N: [2, 2]},
            {S: [2.0**24 + 1, 10.0], M: [2.0**23 + 0.5, 6.0], N: [2, 2]},
            "fail",
        ),
    ],
)
def test_run_one_accepts_only_the_float32_rounding_of_raw_sum(
    tmp_path: pathlib.Path,
    census1_obs: dict[str, float | list[float]],
    sdk_obs: dict[str, float | list[float]],
    status: str,
) -> None:
    def cells(values: dict[str, float | list[float]]) -> ad.AnnData:
        """Two cells. A single value is the first cell's, and the second cell holds 5."""
        adata = _adata(["g1"])
        for column, value in values.items():
            adata.obs[column] = value if isinstance(value, list) else [value, 5.0]
        return adata

    class Stores:
        def get_anndata(
            self, organism: str, obs_filter: object, layers: object = ()
        ) -> ad.AnnData:
            return cells(sdk_obs)

    census1 = canonicalize_anndata(cells(census1_obs))
    record = _run_one(tmp_path, QUERY, Stores(), census1)
    assert record["status"] == status


def test_run_one_names_the_owner_of_an_unsupported_query(
    tmp_path: pathlib.Path,
) -> None:
    """A numeric clause waits for the filter interface. Triage takes any other reason."""
    numeric = replace(QUERY, obs_filter="nnz >= 500")
    record = _run_one(tmp_path, numeric, Empty(), {})
    assert (record["status"], record["owner"]) == ("unsupported", "6.7")
    spatial = replace(QUERY, modality="census_spatial_sequencing")
    record = run_one(
        spatial,
        tmp_path,
        cast("census2_sdk.Census", Empty()),
        {"d1"},
        {"census_spatial_sequencing/homo_sapiens": PINS},
        {"census_spatial_sequencing/homo_sapiens": {"dataset_ids": ["d1"]}},
    )
    assert (record["status"], record["owner"]) == ("unsupported", "6.2")
    assert "no spatial experiment" in record["reason"]


LUNG = ("ds-a", "ds-b")  # the stores of `test_sdk` that hold lung cells
OBS_QUERY = replace(
    QUERY,
    key="obs_metadata_only",
    accessor="get_obs",
    obs_filter="tissue == 'lung'",
    obs_columns=("tissue", "disease"),
    label="obs_metadata_only:census_data/homo_sapiens",
)


@pytest.fixture(scope="module")
def stores(tmp_path_factory: pytest.TempPathFactory) -> pathlib.Path:
    """The three stores of `test_sdk` in `out/`, and their `catalog.parquet`."""
    root = tmp_path_factory.mktemp("stores")
    (root / "out").mkdir()
    for seed, (dataset_id, spec) in enumerate(STORES.items()):
        write_one(root / "out", make_store(dataset_id, *spec, seed=seed), dataset_id)
    build_catalog.build(str(root / "out"), root / "catalog.parquet")
    return root


def _lung_obs(stores: pathlib.Path) -> pd.DataFrame:
    """The lung cells of the `LUNG` stores, as the Census 1 side of `OBS_QUERY`."""
    frames = [ad.read_zarr(stores / "out" / f"v-{d}.zarr").obs for d in LUNG]
    obs = pd.concat(cast("list[pd.DataFrame]", frames))
    columns = ["dataset_id", "observation_joinid", *OBS_QUERY.obs_columns]
    return obs.loc[obs["tissue"].astype(str) == "lung", columns]


def test_answer_reads_obs_through_get_obs(stores: pathlib.Path) -> None:
    census = census2_sdk.open_census([stores / "catalog.parquet"])
    arrays, notes = answer(OBS_QUERY, census, LUNG, PINS)
    want = canonicalize_frame("obs", _lung_obs(stores))
    assert notes == [] and len(want["obs.dataset_id.offsets"]) > 1
    assert arrays.keys() == want.keys()
    for name, array in want.items():
        assert np.array_equal(arrays[name], array), name
    with pytest.raises(NoCells):
        answer(replace(OBS_QUERY, obs_filter="tissue == 'brain'"), census, LUNG, PINS)


@pytest.mark.parametrize(("dropped", "status"), [(0, "pass"), (1, "fail")])
def test_run_one_compares_a_get_obs_result(
    tmp_path: pathlib.Path, stores: pathlib.Path, dropped: int, status: str
) -> None:
    """A result without `var` compares. A cell that Census 1 lacks fails the query."""
    census = census2_sdk.open_census([stores / "catalog.parquet"])
    obs = _lung_obs(stores)
    census1 = canonicalize_frame("obs", obs.iloc[dropped:])
    record = _run_one(tmp_path, OBS_QUERY, census, census1, LUNG)
    assert record["status"] == status
    assert record["cells"] == {"census1": len(obs) - dropped, "sdk": len(obs)}
    assert [f["kind"] for f in record["findings"]] == ["row-set"] * dropped


VAR_QUERY = replace(
    QUERY,
    key="var_all_genes",
    accessor="get_var",
    obs_filter=None,
    label="var_all_genes:census_data/homo_sapiens",
)


def test_answer_reads_var_through_get_var(stores: pathlib.Path) -> None:
    """The clamp does not apply. A var filter selects genes after the read."""
    census = census2_sdk.open_census([stores / "catalog.parquet"])
    reference = reference_var(REFERENCE_IDS)
    arrays, notes = answer(VAR_QUERY, census, LUNG[:1], PINS)
    want = canonicalize_frame("var", reference)
    assert notes == [] and arrays.keys() == want.keys()
    for name, array in want.items():
        assert np.array_equal(arrays[name], array), name
    one_name = replace(VAR_QUERY, var_filter="feature_name in ['GENE3', 'GENE7']")
    arrays, notes = answer(one_name, census, LUNG, PINS)
    assert list(arrays["var.feature_length"]) == [1000, 1000] and len(notes) == 1
    for var_filter in ("feature_length >= 5", "symbol == 'GENE3'"):
        with pytest.raises(Unsupported):
            answer(replace(VAR_QUERY, var_filter=var_filter), census, LUNG, PINS)


@pytest.mark.parametrize(
    ("var_filter", "census1_name", "status"),
    [
        (None, "GENE3", "pass"),
        (None, "OTHER", "fail"),  # a value differs
        ("feature_name == 'GENE3'", "GENE3", "pass"),
        ("feature_name == 'GENE3'", "OTHER", "fail"),
        ("feature_name == 'NONE'", "GENE3", "vacuous"),
    ],
)
def test_run_one_compares_a_get_var_result(
    tmp_path: pathlib.Path,
    stores: pathlib.Path,
    var_filter: str | None,
    census1_name: str,
    status: str,
) -> None:
    """A result without `obs` compares, with the expected differences of `var`."""
    census = census2_sdk.open_census([stores / "catalog.parquet"])
    var = reference_var(REFERENCE_IDS).drop(columns="feature_version")
    var["nnz"] = 5  # a corpus aggregate, which only Census 1 holds
    if var_filter:
        var = var[var["feature_name"] == var_filter.split("'")[1]]
    var = var.replace({"feature_name": {"GENE3": census1_name}})
    query = replace(VAR_QUERY, var_filter=var_filter)
    record = _run_one(tmp_path, query, census, canonicalize_frame("var", var), LUNG)
    assert record["status"] == status
    assert record["genes"]["census1"] == len(var)
    blocking = [f["kind"] for f in record["findings"] if "expected" not in f]
    assert blocking == (["value"] if status == "fail" else [])


PRESENCE_QUERY = replace(
    QUERY,
    key="presence_matrix",
    accessor="get_presence_matrix",
    obs_filter=None,
    label="presence_matrix:census_data/homo_sapiens",
)
PRESENCE_STORES = ("ds-a", "ds-c")  # no store of the two measured an odd gene above 12


@pytest.fixture(scope="module")
def presence(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, pd.DataFrame]:
    """A census of the `PRESENCE_STORES` of `test_sdk`, and its presence matrix."""
    root = tmp_path_factory.mktemp("presence")
    (root / "out").mkdir()
    for seed, dataset_id in enumerate(PRESENCE_STORES):
        adata = make_store(dataset_id, *STORES[dataset_id], seed=seed)
        write_one(root / "out", adata, dataset_id)
    build_catalog.build(str(root / "out"), root / "catalog.parquet")
    census = census2_sdk.open_census([root / "catalog.parquet"])
    return census, census.get_presence_matrix("homo_sapiens")


def _census1_presence(frame: pd.DataFrame) -> Arrays:
    """`frame` in the form of the Census 1 recording.

    The rows and the columns are in the opposite order. One more row holds no gene, as
    the row of a dataset of another experiment does.
    """
    frame = frame.iloc[::-1, ::-1].copy()
    frame.loc["ds-other"] = False
    matrix = sparse.csr_matrix(frame.to_numpy(dtype=np.uint8))
    return presence_arrays(matrix, frame.index, frame.columns)


def _measured_only(frame: pd.DataFrame) -> pd.DataFrame:
    """As the Census 1 build: only the genes that a dataset measured."""
    return frame.loc[:, frame.any()]


def _flipped(frame: pd.DataFrame) -> pd.DataFrame:
    frame = _measured_only(frame)
    frame.iloc[1, 3] = not frame.iloc[1, 3]
    return frame


def _one_more_dataset(frame: pd.DataFrame) -> pd.DataFrame:
    frame = _measured_only(frame)
    frame.loc["ds-z"] = True
    return frame


def _dataset_without_a_gene(frame: pd.DataFrame) -> pd.DataFrame:
    frame = _measured_only(frame)
    frame.loc["ds-c"] = False
    return frame


def _without_a_measured_gene(frame: pd.DataFrame) -> pd.DataFrame:
    return _measured_only(frame).drop(columns=REFERENCE_IDS[0])


def _one_more_gene(frame: pd.DataFrame) -> pd.DataFrame:
    return _measured_only(frame).assign(ENSG99999999999=True)


@pytest.mark.parametrize(
    ("census1_form", "status", "blocking"),
    [
        (pd.DataFrame.copy, "pass", []),
        (_measured_only, "pass", []),
        (_flipped, "fail", ["X.presence.indptr"]),
        (_one_more_dataset, "fail", ["presence datasets, Census 1 only"]),
        (_dataset_without_a_gene, "fail", ["presence datasets, SDK only"]),
        (_without_a_measured_gene, "fail", ["X.presence"]),
        (_one_more_gene, "fail", ["var genes, Census 1 only"]),
    ],
)
def test_run_one_compares_a_presence_result_by_label(
    tmp_path: pathlib.Path,
    presence: tuple[Any, pd.DataFrame],
    census1_form: Callable[[pd.DataFrame], pd.DataFrame],
    status: str,
    blocking: list[str],
) -> None:
    """The order of a side plays no part. A dataset on one side only fails the query.

    A gene that only the SDK holds is an expected difference when no dataset measured it.
    """
    census, sdk = presence
    census1 = census1_form(sdk)
    arrays = _census1_presence(census1)
    record = _run_one(tmp_path, PRESENCE_QUERY, census, arrays, PRESENCE_STORES)
    assert record["status"] == status
    assert record["datasets"] == {"census1": int(census1.any(axis=1).sum()), "sdk": 2}
    assert record["genes"] == {"census1": census1.shape[1], "sdk": len(REFERENCE_IDS)}
    found = {f["where"]: "expected" in f for f in record["findings"]}
    assert [where for where, known in found.items() if not known] == blocking
    assert found.get(SDK_ONLY_GENES, True)
    assert (SDK_ONLY_GENES in found) == (census1_form is not pd.DataFrame.copy)


def test_run_one_refuses_a_presence_result_without_labels(
    tmp_path: pathlib.Path, presence: tuple[Any, pd.DataFrame]
) -> None:
    census, sdk = presence
    arrays = layer_arrays("presence", sparse.csr_matrix(sdk.to_numpy()))
    with pytest.raises(SystemExit, match="names no dataset"):
        _run_one(tmp_path, PRESENCE_QUERY, census, arrays, PRESENCE_STORES)
