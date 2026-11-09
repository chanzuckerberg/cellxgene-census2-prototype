"""Tests for the harness translation of Census 1 value filters into builder form."""

import pytest
from harness import filter_language_builder as flb
from harness.corpus import expand, load_values

from census2_sdk import obs_columns
from census2_sdk.obs_filter import render

obs = obs_columns.Obs("homo_sapiens")


@pytest.mark.parametrize(
    ("text", "builder"),
    [
        ("tissue_general == 'vasculature'", "obs.tissue_general == 'vasculature'"),
        (
            "tissue_general in ['brain', 'lung']",
            "obs.tissue_general.isin(['brain', 'lung'])",
        ),
        (
            "is_primary_data == True and nnz >= 500 and cell_type_ontology_term_id != 'x'",
            "obs.is_primary_data & (obs.nnz >= 500) & (obs.cell_type_ontology_term_id != 'x')",
        ),
        ("not (sex == 'male') or nnz < 3", "(obs.sex != 'male') | (obs.nnz < 3)"),
        (
            "disease_ontology_term_id in ['MONDO:1 || MONDO:2']",
            "obs.disease_ontology_term_id.has_all(['MONDO:1', 'MONDO:2'])",
        ),
        ("donor_id not in ['a', 'b']", "~obs.donor_id.isin(['a', 'b'])"),
    ],
)
def test_builder_of(text: str, builder: str) -> None:
    expr, _ = flb.builder_of(text, obs)
    assert render(expr) == builder


def test_a_joined_value_has_a_note() -> None:
    _, notes = flb.builder_of("disease == 'a || b'", obs)
    assert notes and "has_all" in notes[0]


def test_every_corpus_filter_has_a_builder_form() -> None:
    rows = [flb.classify(q) for q in expand(load_values())]
    with_filter = [r for r in rows if "value_filter" in r]
    assert len(with_filter) == 50
    assert all(r["status"] == "expressible" for r in with_filter)
