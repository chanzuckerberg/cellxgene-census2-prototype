"""Small local stores with every column of the store format's `obs`, for the builder tests.

`make_corpus` writes four stores and a catalog, and a consolidated store of two of them
with a second catalog. The values are random, but each store holds few labels, so a
filter selects some cells of some stores and none of others:

- `ds-1` and `ds-2` share tissues, `ds-3` holds only kidney, and `ds-4` holds blood.
- `disease` holds joined values in `ds-2` and `ds-4`.
- `is_primary_data` is all true in `ds-1`, all false in `ds-3`, and mixed elsewhere.
- `ds-4` holds an extra column, `author_label`, with nulls, which only `obs["..."]` names.
"""

from pathlib import Path
from typing import Any

import anndata
import build_catalog
import numpy as np
import pandas as pd
import scipy.sparse as sp
from census2_dataset_writer import consolidation
from census2_dataset_writer.store_format import OBS_COLUMNS, ColType
from test_sdk import REFERENCE_IDS, reference_var, write_one

SPECS: dict[str, dict[str, Any]] = {
    "ds-1": {
        "n": 60,
        "tissue": ["lung", "liver"],
        "cell_type": ["neutrophil", "alveolar macrophage", "CD8-positive T cell"],
        "disease": ["normal"],
        "primary": "all",
    },
    "ds-2": {
        "n": 45,
        "tissue": ["lung"],
        "cell_type": ["neutrophil", "type II pneumocyte"],
        "disease": ["normal", "COVID-19 || normal", "COVID-19"],
        "primary": "mixed",
    },
    "ds-3": {
        "n": 30,
        "tissue": ["kidney"],
        "cell_type": ["kidney resident macrophage", "podocyte"],
        "disease": ["normal"],
        "primary": "none",
    },
    "ds-4": {
        "n": 50,
        "tissue": ["blood"],
        "cell_type": ["neutrophil", "CD4-positive T cell", "monocyte"],
        "disease": ["HIV infectious disease || visceral leishmaniasis", "normal"],
        "primary": "mixed",
        "author_label": ["a", "b", None],
    },
}
TERM = {
    "lung": "UBERON:0002048",
    "liver": "UBERON:0002107",
    "kidney": "UBERON:0002113",
    "blood": "UBERON:0000178",
}
GENERAL = {"lung": "lung", "liver": "liver", "kidney": "kidney", "blood": "blood"}


def _labels(rng: np.random.Generator, values: list[Any], n: int) -> pd.Categorical:
    return pd.Categorical(rng.choice(np.array(values, dtype=object), n))


def make_store(dataset_id: str, spec: dict[str, Any], seed: int) -> anndata.AnnData:
    """One store with every column of `OBS_COLUMNS`, in the writer's types."""
    rng = np.random.default_rng(seed)
    n = int(spec["n"])
    tissue = _labels(rng, spec["tissue"], n)
    cell_type = _labels(rng, spec["cell_type"], n)
    disease = _labels(rng, spec["disease"], n)
    columns: dict[str, Any] = {
        "dataset_id": pd.Categorical([dataset_id] * n),
        "assay": pd.Categorical(["10x 3' v3"] * n),
        "assay_ontology_term_id": pd.Categorical(["EFO:0009922"] * n),
        "cell_type": cell_type,
        "cell_type_ontology_term_id": pd.Categorical(
            [f"CL:{abs(hash(str(v))) % 10_000_000:07d}" for v in cell_type]
        ),
        "development_stage": _labels(rng, ["adult", "child"], n),
        "development_stage_ontology_term_id": pd.Categorical(["HsapDv:0000087"] * n),
        "disease": disease,
        "disease_ontology_term_id": pd.Categorical(
            [str(v).replace("COVID-19", "MONDO:0100096") for v in disease]
        ),
        "donor_id": _labels(rng, [f"{dataset_id}-d{i}" for i in range(4)], n),
        "is_primary_data": {
            "all": np.ones(n, bool),
            "none": np.zeros(n, bool),
            "mixed": rng.random(n) < 0.5,
        }[spec["primary"]],
        "observation_joinid": pd.array([f"{dataset_id}-j{i:04d}" for i in range(n)]),
        "self_reported_ethnicity": pd.Categorical(["unknown"] * n),
        "self_reported_ethnicity_ontology_term_id": pd.Categorical(["unknown"] * n),
        "sex": _labels(rng, ["female", "male"], n),
        "sex_ontology_term_id": _labels(rng, ["PATO:0000383", "PATO:0000384"], n),
        "suspension_type": pd.Categorical(["cell"] * n),
        "tissue": tissue,
        "tissue_ontology_term_id": pd.Categorical([TERM[str(v)] for v in tissue]),
        "tissue_type": pd.Categorical(["tissue"] * n),
        "tissue_general": pd.Categorical([GENERAL[str(v)] for v in tissue]),
        "tissue_general_ontology_term_id": pd.Categorical(
            [TERM[str(v)] for v in tissue]
        ),
        "raw_sum": rng.integers(100, 5000, n).astype(np.float64),
        "nnz": rng.integers(50, 3000, n).astype(np.int64),
        "raw_mean_nnz": rng.random(n) * 3,
        "raw_variance_nnz": rng.random(n) * 9,
        "n_measured_vars": np.full(n, len(REFERENCE_IDS), dtype=np.int64),
    }
    assert [name for name, _ in OBS_COLUMNS] == list(columns)
    for name, col_type in OBS_COLUMNS:
        if col_type is ColType.CATEGORY:
            columns[name] = pd.Categorical(columns[name])
    obs = pd.DataFrame(columns, index=[f"cell{i}" for i in range(n)])
    if "author_label" in spec:
        labels = rng.choice(np.array(spec["author_label"], dtype=object), n)
        obs["author_label"] = pd.Categorical(labels)
    x = sp.random(
        n, len(REFERENCE_IDS), density=0.4, format="csr", random_state=seed
    ).astype(np.float32)
    x.data = np.floor(x.data * 10) + 1
    var = reference_var(REFERENCE_IDS)
    var["measured"] = True
    adata = anndata.AnnData(X=sp.csr_matrix(x, dtype=np.float32), obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": "NCBITaxon:9606",
    }
    return adata


def make_corpus(root: Path) -> dict[str, Path]:
    """Write the stores, a catalog, and a catalog with a consolidated store of two."""
    prefix, consolidated = root / "out", root / "consolidated"
    prefix.mkdir()
    for seed, (dataset_id, spec) in enumerate(SPECS.items()):
        write_one(prefix, make_store(dataset_id, spec, seed), dataset_id)
    build_catalog.build([str(prefix)], root / "catalog.parquet")
    members = [prefix / f"v-{d}.zarr" for d in ("ds-1", "ds-2", "ds-3")]
    consolidation.consolidate(members, consolidated / "c-123.zarr")
    build_catalog.build(
        [str(prefix)], root / "catalog-c.parquet", consolidated=[str(consolidated)]
    )
    return {
        "root": root,
        "prefix": prefix,
        "catalog": root / "catalog.parquet",
        "catalog_c": root / "catalog-c.parquet",
    }
