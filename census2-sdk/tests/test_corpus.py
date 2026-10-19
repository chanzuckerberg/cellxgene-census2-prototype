"""The bound corpus: its size and the human pins of corpus v1.1."""

import json

import pytest
from harness.corpus import CORPUS_VERSION, VALUES_PATH, expand, parse_picks

KEY = "census_data/homo_sapiens"


def test_parse_picks_sorts_the_ids_of_each_experiment() -> None:
    assert parse_picks([f"{KEY}=c,a,b"], 3) == {KEY: ["a", "b", "c"]}
    assert parse_picks([], 3) == {}


@pytest.mark.parametrize(
    "items",
    [
        [KEY],  # no "="
        ["homo_sapiens=a,b,c"],  # no modality
        [f"{KEY}="],  # empty
        [f"{KEY}=a,b"],  # too few
        [f"{KEY}=a,a,b"],  # a repeated ID
        [f"{KEY}=a,,b"],  # an empty ID
        [f"{KEY}=a,b,c", f"{KEY}=d,e,f"],  # the experiment twice
    ],
)
def test_parse_picks_refuses_a_malformed_pick(items: list[str]) -> None:
    with pytest.raises(ValueError):
        parse_picks(items, 3)


def test_the_corpus_binds_74_queries() -> None:
    assert CORPUS_VERSION == "1.1"
    assert len(expand()) == 74


def test_the_human_disease_filter_lists_the_joined_values_of_the_term() -> None:
    values = json.loads(VALUES_PATH.read_text())["census_data/homo_sapiens"]
    term = values["disease_ontology_term_id"]
    matches = values["disease_ontology_term_id_matches"]
    assert " || " not in term
    assert matches and all(term in m.split(" || ") for m in matches)
    (query,) = [
        q
        for q in expand()
        if q.key == "multivalued_disease" and q.experiment == "homo_sapiens"
    ]
    assert query.obs_filter == f"disease_ontology_term_id in {matches}"
