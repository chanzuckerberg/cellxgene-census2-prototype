"""Tests for the typed filter builder, without stores. No network, no AWS.

mypy checks this file under `--strict`, with the error code `deprecated` on (see
`pyproject.toml`). `assert_type` checks the static type of each operator. A line that
misuses an operator carries `# type: ignore[<code>]`: mypy then reports an unused
ignore when it does not see the misuse, so these lines test the static types too. At
run time, each misuse raises `TypeError` with a hint.
"""

from typing import assert_type

import generate_obs_columns
import pytest

from census2_sdk import obs_columns
from census2_sdk.obs_filter import (
    And,
    AnyColumn,
    BoolColumn,
    CategoryColumn,
    Compare,
    Const,
    Expr,
    FilterError,
    FloatColumn,
    IntColumn,
    IsIn,
    MultiValueColumn,
    Not,
    Or,
    StringColumn,
    Truth,
    as_expr,
    from_mapping,
    label_test,
    leaves,
    normalize,
    render,
)

obs = obs_columns.Obs("homo_sapiens")
mouse = obs_columns.Obs("mus_musculus")


def test_the_columns_are_the_store_format() -> None:
    assert generate_obs_columns.current() == generate_obs_columns.expected()
    assert isinstance(obs.disease, MultiValueColumn)
    assert isinstance(obs.disease_ontology_term_id, MultiValueColumn)
    assert isinstance(obs.is_primary_data, BoolColumn)
    assert isinstance(obs.observation_joinid, StringColumn)
    assert isinstance(obs.nnz, IntColumn)
    assert isinstance(obs.raw_sum, FloatColumn)
    assert isinstance(obs.tissue, CategoryColumn)


def test_static_types() -> None:
    assert_type(obs.tissue, CategoryColumn)
    assert_type(obs.tissue == "lung", Expr)
    assert_type(obs.tissue != "lung", Expr)
    assert_type(obs.tissue.isin(["lung", "liver"]), Expr)
    assert_type(obs.cell_type.str.startswith("CD8"), Expr)
    assert_type(obs.cell_type.str.contains("macro", case=False), Expr)
    assert_type(obs.cell_type.str.match(r"CD(4|8)"), Expr)
    assert_type(obs.disease.has_any(["COVID-19"]), Expr)
    assert_type(obs.disease.has_all(["COVID-19", "normal"]), Expr)
    assert_type(obs.nnz >= 500, Expr)
    assert_type(obs.raw_sum < 2.5, Expr)
    assert_type(obs.nnz == 3, Expr)
    assert_type(obs.nnz.between(1, 2), Expr)
    assert_type(obs.observation_joinid == "x", Expr)
    assert_type(obs.observation_joinid.isin(["x"]), Expr)
    assert_type(~obs.is_primary_data, Expr)
    assert_type(obs.is_primary_data & (obs.tissue == "lung"), Expr)
    assert_type((obs.tissue == "lung") & obs.is_primary_data, Expr)
    assert_type((obs.tissue == "lung") | (obs.nnz > 1), Expr)
    assert_type(~(obs.tissue == "lung"), Expr)
    assert_type(obs.tissue.isna(), Expr)
    assert_type(obs["author_label"], AnyColumn)
    assert_type(obs["author_label"] == "a", Expr)
    assert_type(obs["author_label"] > 1, Expr)


def test_a_multi_valued_column_refuses_equality_and_isin() -> None:
    with pytest.raises(TypeError, match="has_any"):
        obs.disease == "COVID-19"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="has_all"):
        obs.disease != "COVID-19"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="multi-valued"):
        obs.disease_ontology_term_id.isin(["MONDO:0100096"])  # type: ignore[deprecated]


def test_an_operator_of_another_type_is_refused() -> None:
    with pytest.raises(TypeError, match="compare it with a str"):
        obs.tissue == 3  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="compare it with a number"):
        obs.nnz == "3"  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError, match="compare it with a bool"):
        obs.is_primary_data == 1  # type: ignore[deprecated]  # noqa: B015
    with pytest.raises(TypeError):
        obs.tissue > "a"  # type: ignore[operator]  # noqa: B015
    with pytest.raises(TypeError, match="holds numbers"):
        obs.nnz < "a"  # type: ignore[operator]  # noqa: B015
    with pytest.raises(TypeError, match="list of values"):
        obs.tissue.isin("lung")
    with pytest.raises(ValueError, match="one or more"):
        obs.disease.has_any([])
    with pytest.raises(FilterError, match="unterminated"):
        obs.cell_type.str.match("CD(4")


def test_and_or_not_raise_with_a_hint() -> None:
    lung = obs.tissue == "lung"
    with pytest.raises(TypeError, match=r"& \| ~"):
        bool(lung)
    with pytest.raises(TypeError, match="parentheses"):
        _ = lung and obs.is_primary_data
    with pytest.raises(TypeError, match="truth value"):
        _ = not obs.is_primary_data
    with pytest.raises(TypeError, match="truth value"):
        if obs.tissue:
            pass


def test_a_filter_remembers_its_organism() -> None:
    f = (obs.tissue == "lung") & obs.is_primary_data
    assert f.organism == "homo_sapiens"
    with pytest.raises(FilterError, match="do not combine"):
        f & (mouse.tissue == "lung")
    with pytest.raises(TypeError, match="combines with a filter"):
        f & "lung"  # type: ignore[operator]


def test_repr_is_the_pandas_form() -> None:
    f = (obs.tissue == "lung") & ~obs.cell_type.isin(["a", "b"]) | obs.is_primary_data
    assert repr(f) == (
        "((obs.tissue == 'lung') & ~obs.cell_type.isin(['a', 'b'])) | obs.is_primary_data"
    )
    assert repr(obs.disease.has_any(["x"])) == "obs.disease.has_any(['x'])"
    assert repr(~obs.is_primary_data) == "~obs.is_primary_data"
    assert repr(obs["x"].notna()) == "obs['x'].notna()"
    assert repr(obs.nnz.between(1, 5, inclusive="left")) == (
        "obs.nnz.between(1, 5, inclusive='left')"
    )
    assert repr(obs.cell_type.str.contains("m", case=False, regex=False)) == (
        "obs.cell_type.str.contains('m', case=False, regex=False)"
    )
    assert render(normalize(obs.tissue != "lung")) == "obs.tissue != 'lung'"


def test_normalize_pushes_not_to_the_terms() -> None:
    a, b = obs.tissue == "lung", obs.nnz > 5
    tree = normalize(~(a & b))
    assert isinstance(tree, Or)
    assert all(isinstance(c, Not) for c in tree.children)
    assert normalize(~~a) == normalize(a)
    assert normalize(obs.tissue != "lung") == normalize(~(obs.tissue == "lung"))
    # A constant folds away.
    assert normalize(a & Const("", True)) == normalize(a)
    assert normalize(a | Const("", True)) == Const("homo_sapiens", True)
    assert normalize(obs.tissue.isin([])) == Const("homo_sapiens", False)
    assert normalize(~obs.tissue.isin([])) == Const("homo_sapiens", True)


def test_normalize_merges_terms_on_one_column() -> None:
    tree = normalize((obs.tissue == "lung") | (obs.tissue == "liver") | (obs.nnz > 1))
    assert isinstance(tree, Or)
    assert IsIn("homo_sapiens", "tissue", False, ("liver", "lung")) in tree.children
    tree = normalize(~(obs.tissue == "lung") & ~obs.tissue.isin(["liver"]))
    assert tree == Not(
        "homo_sapiens", IsIn("homo_sapiens", "tissue", False, ("liver", "lung"))
    )
    tree = normalize(obs.disease.has_all(["a"]) & obs.disease.has_all(["b"]))
    assert render(tree) == "obs.disease.has_all(['a', 'b'])"


def test_between_is_two_comparisons() -> None:
    tree = normalize(obs.nnz.between(1, 5, inclusive="right"))
    assert tree == And(
        "homo_sapiens",
        (
            Compare("homo_sapiens", "nnz", False, ">", 1),
            Compare("homo_sapiens", "nnz", False, "<=", 5),
        ),
    )


def test_a_boolean_column_is_a_filter() -> None:
    assert as_expr(obs.is_primary_data) == Truth(
        "homo_sapiens", "is_primary_data", False
    )
    assert normalize(obs.is_primary_data == False) == normalize(~obs.is_primary_data)
    assert normalize(obs.is_primary_data != True) == normalize(~obs.is_primary_data)


@pytest.mark.parametrize(
    ("term", "value", "hit"),
    [
        (obs.tissue == "lung", "lung", True),
        (obs.tissue.isin(["a", "b"]), "c", False),
        (obs.disease.has_any(["COVID-19"]), "COVID-19 || normal", True),
        (obs.disease.has_any(["COVID-19"]), "COVID-19-like", False),
        (obs.disease.has_all(["a", "b"]), "b || a", True),
        (obs.disease.has_all(["a", "b"]), "a", False),
        (obs.cell_type.str.startswith(("CD4", "CD8")), "CD8-positive T cell", True),
        (obs.cell_type.str.contains("MACRO", case=False), "lung macrophage", True),
        (obs.cell_type.str.contains("a.c", regex=False), "abc", False),
        (obs.cell_type.str.match("CD(4|8)-pos"), "CD4-positive", True),
        (obs.cell_type.str.match("positive"), "CD4-positive", False),
        (obs.cell_type.str.fullmatch("CD4"), "CD4-positive", False),
    ],
)
def test_label_test(term: Expr, value: str, hit: bool) -> None:
    (leaf,) = leaves(normalize(term))
    assert label_test(leaf)(value) is hit


def test_from_mapping_is_the_mapping_meaning() -> None:
    tree = from_mapping(
        obs, {"tissue": ["lung", "liver"], "disease": "flu"}, ["disease"]
    )
    assert render(normalize(tree)) == (
        "obs.tissue.isin(['liver', 'lung']) & obs.disease.has_any(['flu'])"
    )
