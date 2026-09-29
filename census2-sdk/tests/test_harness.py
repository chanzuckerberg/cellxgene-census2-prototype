"""The corruption suite of the comparator: one deliberate mutation per test.

Each corruption is applied to a copy of one AnnData and canonicalized again, because
that is how a broken stack presents it. Every test asserts that the comparator finds
the corruption and that the report names the right thing.
"""

from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from harness.canonical import canonicalize_anndata
from harness.comparator import Finding, compare

RNG = np.random.default_rng(7)  # fixed seed: the fixture must be reproducible


def make_adata(n_obs: int = 40, n_var: int = 15) -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "dataset_id": pd.Categorical([f"ds{i % 3}" for i in range(n_obs)]),
            "observation_joinid": [f"cell{i:04d}" for i in range(n_obs)],
            "cell_type": pd.Categorical(
                RNG.choice(["neuron", "astrocyte", "microglia"], n_obs)
            ),
            "n_genes": RNG.integers(100, 5000, n_obs).astype(np.int32),
        }
    )
    var = pd.DataFrame(
        {
            "feature_id": [f"ENSG{i:011d}" for i in range(n_var)],
            "feature_name": [f"GENE{i}" for i in range(n_var)],
        }
    )
    raw = sp.random(
        n_obs,
        n_var,
        density=0.4,
        random_state=7,
        dtype=np.float32,
        data_rvs=lambda n: RNG.integers(1, 100, n).astype(np.float32),
    ).tocsr()
    adata = ad.AnnData(X=raw.copy(), obs=obs, var=var)
    adata.layers["normalized"] = (
        raw.multiply(1.0 / raw.sum(axis=1)).tocsr().astype(np.float32)
    )
    adata.uns["x_layer_name"] = "raw"
    return adata


def x_of(adata: ad.AnnData) -> sp.csr_matrix:
    return cast("sp.csr_matrix", adata.X).tocsr()


def layer_of(adata: ad.AnnData, name: str) -> sp.csr_matrix:
    return cast("sp.csr_matrix", adata.layers[name]).tocsr()


@pytest.fixture
def truth() -> ad.AnnData:
    return make_adata()


def kinds(findings: list[Finding]) -> set[str]:
    return {f.kind for f in findings}


def test_equal_to_itself(truth: ad.AnnData) -> None:
    assert compare(canonicalize_anndata(truth), canonicalize_anndata(truth)) == []


def test_row_order_is_invisible(truth: ad.AnnData) -> None:
    shuffled = truth[RNG.permutation(truth.n_obs)].copy()
    shuffled.uns["x_layer_name"] = "raw"
    assert compare(canonicalize_anndata(truth), canonicalize_anndata(shuffled)) == []


def test_dropped_row(truth: ad.AnnData) -> None:
    corrupt = truth[1:].copy()
    corrupt.uns["x_layer_name"] = "raw"
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert kinds(findings) == {"strict-subset"}
    assert "cell0000" in findings[0].detail


def test_duplicated_row(truth: ad.AnnData) -> None:
    corrupt = ad.concat([truth, truth[:1]], merge="same")
    corrupt.uns["x_layer_name"] = "raw"
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert any(
        f.kind == "row-set" and "duplicated" in f.detail and "cell0000" in f.detail
        for f in findings
    )


def test_single_changed_expression_value(truth: ad.AnnData) -> None:
    corrupt = truth.copy()
    x = x_of(corrupt)
    x.data[5] += 1.0
    corrupt.X = x
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert [f.kind for f in findings] == ["value"]
    assert findings[0].where == "X.raw"
    assert "cell" in findings[0].detail and "ENSG" in findings[0].detail


def test_permuted_column(truth: ad.AnnData) -> None:
    corrupt = truth.copy()
    corrupt.obs["cell_type"] = np.roll(corrupt.obs["cell_type"].to_numpy(), 1)
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert all(f.kind == "value" and f.where == "obs.cell_type" for f in findings)
    assert "|cell" in findings[0].detail  # named by identity, not by position


def test_dropped_metadata_column(truth: ad.AnnData) -> None:
    corrupt = truth.copy()
    corrupt.obs = cast("pd.DataFrame", corrupt.obs).drop(columns=["n_genes"])
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert findings == [Finding("column-set", "obs.n_genes", "present on left only")]


def test_remapped_categorical(truth: ad.AnnData) -> None:
    corrupt = truth.copy()
    cell_type = cast("Any", corrupt.obs["cell_type"])
    corrupt.obs["cell_type"] = cell_type.cat.rename_categories({"astrocyte": "neuron2"})
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert findings[0].kind == "value" and findings[0].where == "obs.cell_type"
    assert "'astrocyte'" in findings[0].detail and "'neuron2'" in findings[0].detail


def test_truncated_strict_subset(truth: ad.AnnData) -> None:
    corrupt = truth[: truth.n_obs // 2].copy()
    corrupt.uns["x_layer_name"] = "raw"
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert kinds(findings) == {"strict-subset"}
    assert f"{truth.n_obs - truth.n_obs // 2} rows missing" in findings[0].detail


def test_normalized_tolerance_is_bounded(truth: ad.AnnData) -> None:
    """The stated rtol forgives last-ulp noise and nothing more."""
    corrupt = truth.copy()
    layer = layer_of(corrupt, "normalized")
    within = layer.copy()
    within.data[0] *= 1 + 1e-7  # inside rtol=1e-6: the same answer
    corrupt.layers["normalized"] = within
    assert compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt)) == []
    beyond = layer.copy()
    beyond.data[0] *= 1 + 1e-4  # beyond rtol: a finding
    corrupt.layers["normalized"] = beyond
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert [f.where for f in findings] == ["X.normalized"]


def test_raw_is_exact(truth: ad.AnnData) -> None:
    """Raw counts get no tolerance."""
    corrupt = truth.copy()
    x = x_of(corrupt)
    x.data[0] += np.float32(0.001)
    corrupt.X = x
    findings = compare(canonicalize_anndata(truth), canonicalize_anndata(corrupt))
    assert [f.where for f in findings] == ["X.raw"]
