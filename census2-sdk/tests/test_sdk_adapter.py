"""The SDK adapter's filter translation and the driver's gene restriction and owners."""

import pathlib
from dataclasses import replace
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from harness.canonical import (
    Arrays,
    canonicalize_anndata,
    frame_arrays,
    layer_arrays,
    write_result,
)
from harness.comparator import Finding
from harness.compare_stacks import expected, owner, restrict_genes, run_one
from harness.corpus import BoundQuery
from harness.sdk_adapter import NoCells, Unsupported, answer, clauses, translate
from scipy import sparse

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
    def get_anndata(self, organism: str, obs_filter: object) -> None:
        return None


def test_answer_lets_an_sdk_defect_fail_the_run() -> None:
    class Broken:
        def get_anndata(self, organism: str, obs_filter: object) -> None:
            raise SourceConflict("the cell key is not unique across the result")

    with pytest.raises(SourceConflict):
        answer(QUERY, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
    unindexed = replace(QUERY, obs_filter="is_primary_data == True")
    with pytest.raises(Unsupported, match="not a catalog column"):
        answer(unindexed, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
    with pytest.raises(NoCells):
        answer(QUERY, cast("census2_sdk.Census", Empty()), ["d1"], PINS)


def _run_one(
    tmp_path: pathlib.Path, query: BoundQuery, census: object, census1: Arrays
) -> dict[str, Any]:
    """`run_one` over a Census 1 result of `census1` and a clamp of one dataset."""
    meta = {"clamp": ["d1"], "env": {"corpus_version": "1.1"}}
    write_result(tmp_path / query.label.replace("/", "__"), census1, meta)
    return run_one(
        query,
        tmp_path,
        cast("census2_sdk.Census", census),
        {"d1"},
        {"census_data/homo_sapiens": PINS},
        {"census_data/homo_sapiens": {"dataset_ids": ["d1"]}},
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
        def get_anndata(self, organism: str, obs_filter: object) -> ad.AnnData:
            return _adata(sdk_genes)

    census1 = canonicalize_anndata(_adata(["g1", "g2"]))
    record = _run_one(tmp_path, QUERY, Stores(), census1)
    assert record["status"] == status


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
        def get_anndata(self, organism: str, obs_filter: object) -> ad.AnnData:
            return cells(sdk_obs)

    census1 = canonicalize_anndata(cells(census1_obs))
    record = _run_one(tmp_path, QUERY, Stores(), census1)
    assert record["status"] == status


def test_run_one_names_the_task_that_adds_an_unsupported_accessor(
    tmp_path: pathlib.Path,
) -> None:
    query = replace(QUERY, accessor="get_obs")
    record = _run_one(tmp_path, query, Empty(), {})
    assert (record["status"], record["owner"]) == ("unsupported", "6.7")
