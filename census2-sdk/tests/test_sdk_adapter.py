"""The 6.2 adapter's filter translation and the driver's gene restriction and owners."""

from dataclasses import replace
from typing import cast

import numpy as np
import pandas as pd
import pytest
from harness.canonical import frame_arrays, layer_arrays
from harness.comparator import Finding
from harness.compare_stacks import owner, relabel, restrict_genes
from harness.corpus import BoundQuery
from harness.sdk_adapter import Unsupported, answer, clauses, translate
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
        == "4.5"
    )
    tissue = Finding("strict-subset", "obs", "right is a strict subset")
    assert owner(tissue, {"tissue_general_ontology_term_id", "dataset_id"})[0] == "4.5"
    assert owner(Finding("value", "X.raw", "cell x gene"), set())[0] == "6.2"
    assert owner(Finding("value", "obs.donor_id", "identity x"), set())[0] == "6.2"


def test_relabel_gives_a_label_difference_with_agreeing_terms_to_recuration() -> None:
    def finding(kind: str, where: str) -> dict[str, str]:
        return {"kind": kind, "where": where, "owner": "4.4", "reason": ""}

    compared = {
        "obs.cell_type",
        "obs.cell_type_ontology_term_id",
        "obs.tissue",
        "obs.tissue_ontology_term_id",
        "obs.sex",
        "obs.sex_ontology_term_id",
    }
    findings = [
        finding("value", "obs.cell_type"),  # its term agrees
        finding("value", "obs.tissue"),  # its term differs
        finding("value", "obs.tissue_ontology_term_id"),
        finding("value", "obs.sex"),  # its term has another kind of finding
        finding("dtype", "obs.sex_ontology_term_id"),
        finding("value", "obs.disease"),  # its term is not on both sides
    ]
    relabel(findings, compared)
    assert [f["owner"] for f in findings] == [
        "re-curation",
        "4.4",
        "4.4",
        "4.4",
        "4.4",
        "4.4",
    ]


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


def test_answer_lets_an_sdk_defect_fail_the_run() -> None:
    class Broken:
        def get_anndata(self, organism: str, obs_filter: object) -> None:
            raise SourceConflict("the cell key is not unique across the result")

    query = BoundQuery(
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
    with pytest.raises(SourceConflict):
        answer(query, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
    unindexed = replace(query, obs_filter="is_primary_data == True")
    with pytest.raises(Unsupported, match="not a catalog column"):
        answer(unindexed, cast("census2_sdk.Census", Broken()), ["d1"], PINS)
