"""Unit tests for the store builder. No network, no AWS, no container.

These cover the properties the module docstring promises, because a caller
relies on them: the round trip, the gene projection, the slot policy, admission
as a normal outcome, the memory contract, provenance, and destination handling.

The gene reference is a 24-row file in a `refdata` directory, staged for every
test, and the default configuration names it. Its order is the reverse of the
fixture's gene order, so a projection that keeps the source order fails.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, cast

import anndata
import census2_dataset_writer
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr
from anndata.io import read_elem

TAXON = "NCBITaxon:9606"
REF_IDS = [f"ENSG{j:05d}" for j in reversed(range(24))]

# A Discover h5ad holds no identifiers, so every store needs these from its caller.
IDS = {"dataset_id": "dataset-1", "dataset_version_id": "version-1"}


def build(
    f: h5py.File, **kwargs: Any
) -> tuple[anndata.AnnData | None, census2_dataset_writer.Result]:
    """`census2_dataset_writer.build`, with the test identifiers."""
    return census2_dataset_writer.build(
        f,
        dataset_id=IDS["dataset_id"],
        dataset_version_id=IDS["dataset_version_id"],
        **kwargs,
    )


def create(source: Any, dest: Any, **kwargs: Any) -> census2_dataset_writer.Result:
    """`census2_dataset_writer.create`, with the test identifiers."""
    return census2_dataset_writer.create(
        source,
        dest,
        dataset_id=IDS["dataset_id"],
        dataset_version_id=IDS["dataset_version_id"],
        **kwargs,
    )


@pytest.fixture(autouse=True)
def refdata(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Stage the test reference, and make the default configuration name it."""
    root = tmp_path / "refdata"
    root.mkdir()
    rows = "".join(f"{g},name{g},1,100,protein_coding\n" for g in REF_IDS)
    data = gzip.compress(rows.encode())
    (root / "genes_test.csv.gz").write_bytes(data)
    spec = {
        "file": "genes_test.csv.gz",
        "version": "TEST v1",
        "sha256": hashlib.sha256(data).hexdigest(),
        "url": "https://example.invalid/genes_test.csv.gz",
    }
    monkeypatch.setitem(
        census2_dataset_writer.DEFAULTS, "gene_references", {TAXON: spec}
    )
    return root


def dataset(
    n_obs: int = 40,
    n_var: int = 12,
    *,
    dense: bool = False,
    raw: bool = True,
    extras: bool = False,
    seed: int = 0,
) -> anndata.AnnData:
    """A small AnnData shaped like a Discover dataset.

    With `raw`, `raw.X` holds the counts and `X` holds half of each count, so a
    test can tell which matrix a store took.
    """
    rng = np.random.default_rng(seed)
    x: Any = (
        rng.random((n_obs, n_var), dtype="float32")
        if dense
        else sp.random(
            n_obs, n_var, density=0.3, format="csr", dtype="float32", random_state=seed
        )
    )
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(rng.choice(["a", "b"], n_obs))},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {
            "feature_name": [f"gene{j}" for j in range(n_var)],
            "feature_biotype": pd.Categorical(["gene"] * n_var),
            "feature_reference": pd.Categorical([TAXON] * n_var),
        },
        index=[f"ENSG{j:05d}" for j in range(n_var)],
    )
    adata = anndata.AnnData(X=x, obs=obs, var=var)
    adata.uns = {
        "schema_version": "7.1.0",
        "organism_ontology_term_id": TAXON,
        "title": "not a Census uns key",
    }
    if extras:
        adata.obsm["X_pca"] = rng.random((n_obs, 3), dtype="float32")
        adata.varm["gene_stuff"] = rng.random((n_var, 2), dtype="float32")
        adata.layers["counts"] = x
    if raw:
        adata.raw = adata
        adata.X = x * 0.5
    return adata


def written(
    tmp_path: Path, **kwargs: Any
) -> tuple[Path, census2_dataset_writer.Result]:
    """Write a dataset to h5ad, build a store from it, return both."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(**kwargs).write_h5ad(source)
    result = create(source, dest, refdata=tmp_path / "refdata")
    return dest, result


def counts(adata: Any) -> tuple[np.ndarray[Any, Any], list[str]]:
    """The source counts, dense, and their feature_ids: raw when it exists."""
    part = adata.raw if adata.raw is not None else adata
    x = part.X.toarray() if hasattr(part.X, "toarray") else np.asarray(part.X)
    return x, [str(i) for i in part.var.index]


def projected(adata: Any) -> np.ndarray[Any, Any]:
    """The source counts on the reference genes, by brute force."""
    x, ids = counts(adata)
    out = np.zeros((x.shape[0], len(REF_IDS)), dtype=x.dtype)
    for j, gene in enumerate(ids):
        if gene in REF_IDS:
            out[:, REF_IDS.index(gene)] = x[:, j]
    return out


def stored_x(dest: Path) -> sp.csr_matrix:
    x: Any = read_elem(zarr.open_consolidated(dest, mode="r")["X"])
    assert isinstance(x, (sp.csr_matrix, sp.csr_array)), type(x)
    return sp.csr_matrix(x)


def values(array: Any) -> np.ndarray[Any, Any]:
    """A whole array as values, so a comparison ignores the stored width.

    Typed as Any at both ends: a zarr group member and a read_elem result are
    each a wide union, so this is where the dynamic boundary is named.
    """
    return np.asarray(array)


# --- the round trip -------------------------------------------------------


@pytest.mark.parametrize(
    ("dense", "raw"), [(False, True), (False, False), (True, False), (True, True)]
)
def test_round_trip(tmp_path: Path, dense: bool, raw: bool) -> None:
    dest, result = written(tmp_path, dense=dense, raw=raw)
    assert result["admitted"] is True

    source: Any = anndata.read_h5ad(tmp_path / "source.h5ad")
    store: Any = zarr.open_consolidated(dest, mode="r")

    # Compare through read_elem, the same surface a reader uses. Reaching into
    # the layout instead is misleading: obs/_index reads back as a group.
    obs: Any = read_elem(store["obs"])
    assert list(obs.index) == list(source.obs.index)
    # The store adds `dataset_id`, first, and keeps every source column.
    assert list(obs.columns) == ["dataset_id", *source.obs.columns]
    assert set(obs["dataset_id"]) == {IDS["dataset_id"]}
    assert obs["cell_type"].dtype == source.obs["cell_type"].dtype

    # X is the counts on the reference genes, and sparse even from a dense source.
    x = stored_x(dest)
    assert np.array_equal(x.toarray(), projected(source))


def test_the_store_holds_counts_on_the_reference(tmp_path: Path) -> None:
    dest, result = written(tmp_path, raw=True)
    store: Any = zarr.open_consolidated(dest, mode="r")
    source: Any = anndata.read_h5ad(tmp_path / "source.h5ad")

    var: Any = read_elem(store["var"])
    assert list(var.index) == REF_IDS
    assert list(var.columns) == census2_dataset_writer.genes.COLUMNS
    assert result["gene_reference_version"] == "TEST v1"
    # raw.X holds the counts, and X the halved values: the store took raw.X.
    assert "raw" not in result["slots"]
    assert "raw" not in store or store["raw"].attrs.get("encoding-type") != "raw"
    assert np.array_equal(stored_x(dest).toarray(), projected(source))
    assert not np.array_equal(stored_x(dest).toarray(), 0.5 * projected(source))


def test_the_projection_keeps_every_stored_value(tmp_path: Path) -> None:
    """The remap moves values and drops none, and leaves the CSR canonical."""
    dest, _ = written(tmp_path, raw=False)
    source: Any = anndata.read_h5ad(tmp_path / "source.h5ad")
    x = stored_x(dest)
    assert x.nnz == source.X.nnz
    assert x.has_canonical_format
    assert np.array_equal(values(x.indptr), values(source.X.indptr))


def test_non_genes_and_other_organisms_are_dropped(tmp_path: Path) -> None:
    """Disposition row 12: a spike-in and a viral gene have no reference position."""
    adata = dataset(n_var=4, raw=False)
    extra = pd.DataFrame(
        {
            "feature_name": ["ERCC-00002", "ORF1ab"],
            "feature_biotype": pd.Categorical(["spike-in", "gene"]),
            "feature_reference": pd.Categorical(
                ["NCBITaxon:32630", "NCBITaxon:2697049"]
            ),
        },
        index=["ERCC-00002", "ENSSASG00005000002"],
    )
    rng = np.random.default_rng(1)
    extra_x = rng.integers(1, 5, (adata.n_obs, 2)).astype("float32")
    x, _ = counts(adata)
    obs: Any = adata.obs
    source_var: Any = adata.var
    var = pd.concat([source_var, extra])
    adata = anndata.AnnData(
        X=sp.csr_matrix(np.hstack([x, extra_x])), obs=obs, var=var, uns=dict(adata.uns)
    )
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    adata.write_h5ad(source)

    result = create(source, dest, refdata=tmp_path / "refdata")

    assert result["dropped_var_features"] == ["ERCC-00002", "ENSSASG00005000002"]
    assert np.array_equal(stored_x(dest).toarray(), projected(adata))


def test_a_gene_outside_the_reference_is_refused(tmp_path: Path) -> None:
    adata = dataset(n_var=4, raw=False)
    adata.var.index = pd.Index(["ENSG00000", "ENSG00001", "ENSG00002", "ENSG99999"])
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    adata.write_h5ad(source)

    with pytest.raises(
        census2_dataset_writer.genes.GeneReferenceMismatch, match="ENSG99999"
    ):
        create(source, dest, refdata=tmp_path / "refdata")
    assert not dest.exists()


@pytest.mark.parametrize("column", ["feature_biotype", "feature_reference"])
def test_a_var_without_a_feature_column_is_refused(refdata: Path, column: str) -> None:
    adata = dataset(raw=False)
    var: Any = adata.var  # the getter's union includes Dataset2D
    adata.var = var.drop(columns=[column])
    with pytest.raises(ValueError, match=f"var lacks {column}"):
        census2_dataset_writer.harmonize(adata, refdata=refdata)


def test_a_reference_with_another_digest_is_refused(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    config = census2_dataset_writer.load_config()
    config["gene_references"] = {
        TAXON: {**config["gene_references"][TAXON], "sha256": "0" * 64}
    }
    with pytest.raises(ValueError, match="digest"):
        create(
            source, tmp_path / "store.zarr", config=config, refdata=tmp_path / "refdata"
        )


def test_an_admitted_dataset_needs_refdata(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    with pytest.raises(ValueError, match="refdata"):
        create(source, tmp_path / "store.zarr")


def test_select_counts_takes_raw_and_drops_it() -> None:
    adata = dataset(extras=True, raw=True)
    raw_x, raw_ids = counts(adata)
    selected = census2_dataset_writer.select_counts(adata)
    x: Any = selected.X
    assert selected.raw is None
    assert np.array_equal(x.toarray(), raw_x)
    assert list(selected.var.index) == raw_ids
    # These index the genes of the normalized X, not of raw.X. anndata 0.13
    # keeps X in layers[None].
    assert [k for k in selected.layers if k is not None] == []
    assert list(selected.varm) == []
    assert list(selected.obsm) == ["X_pca"]


def test_select_counts_keeps_a_dataset_without_raw() -> None:
    adata = dataset(raw=False)
    assert census2_dataset_writer.select_counts(adata) is adata


def test_project_genes_refuses_a_dataset_with_raw(refdata: Path) -> None:
    with pytest.raises(ValueError, match="select_counts"):
        census2_dataset_writer.project_genes(dataset(raw=True), refdata=refdata)


def test_default_references_pin_human_and_mouse() -> None:
    references = census2_dataset_writer.genes.REFERENCES
    assert set(references) == {"NCBITaxon:9606", "NCBITaxon:10090"}
    for spec in references.values():
        assert census2_dataset_writer.genes.CURATION_TAG in spec["url"]
        assert spec["url"].endswith(spec["file"])


def test_consolidated_metadata_is_present(tmp_path: Path) -> None:
    dest, _ = written(tmp_path)
    assert zarr.open_consolidated(dest, mode="r") is not None


# --- the summary record -----------------------------------------------------


def test_summary_record_is_in_the_root_attributes(tmp_path: Path) -> None:
    dest, result = written(tmp_path)
    expected = {
        "version": census2_dataset_writer.store_format.SUMMARY_VERSION,
        **IDS,
        "gene_reference_version": result["gene_reference_version"],
    }
    assert census2_dataset_writer.read_summary(dest) == expected
    root = json.loads((dest / "zarr.json").read_text())
    assert root["attributes"]["cellxgene-census2"]["summary"] == expected
    # anndata's own attributes stay beside the record.
    assert root["attributes"]["encoding-type"] == "anndata"


def test_summary_record_keeps_the_consolidated_metadata(tmp_path: Path) -> None:
    """The record rewrites the root zarr.json, which must keep every consolidated node."""
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    with h5py.File(source) as f:
        adata, result = build(f, refdata=tmp_path / "refdata")
        assert adata is not None
        summary = census2_dataset_writer.summary_record(result)
        census2_dataset_writer.write(adata, tmp_path / "store.zarr", summary=summary)
        adata.write_zarr(tmp_path / "plain.zarr")

    def nodes(path: Path) -> set[str]:
        root = json.loads((path / "zarr.json").read_text())
        return set(root["consolidated_metadata"]["metadata"])

    assert nodes(tmp_path / "store.zarr") == nodes(tmp_path / "plain.zarr")
    assert zarr.open_consolidated(tmp_path / "store.zarr", mode="r") is not None
    assert anndata.read_zarr(tmp_path / "store.zarr").n_obs == adata.n_obs


def test_a_store_without_a_summary_record_reads_as_none(tmp_path: Path) -> None:
    anndata.AnnData(X=np.zeros((2, 2))).write_zarr(tmp_path / "plain.zarr")
    assert census2_dataset_writer.read_summary(tmp_path / "plain.zarr") is None


GOOD_RECORD = {
    "version": 1,
    "dataset_id": "d",
    "dataset_version_id": "v",
    "gene_reference_version": "TEST v1",
}


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"version": 2}, "version 2"),
        ({"version": True}, "version True"),
        ({"dataset_id": 7}, "dataset_id is 7"),
        ({"gene_reference_version": ""}, "gene_reference_version is ''"),
        ({"gene_reference_version": None}, "gene_reference_version is None"),
        ({"extra": 1}, "unknown \\['extra'\\]"),
    ],
)
def test_validate_summary_refuses_a_malformed_record(
    change: dict[str, Any], message: str
) -> None:
    with pytest.raises(
        census2_dataset_writer.store_format.SummaryRecordError, match=message
    ):
        census2_dataset_writer.store_format.validate_summary({**GOOD_RECORD, **change})


def test_validate_summary_refuses_a_missing_field() -> None:
    record = {k: v for k, v in GOOD_RECORD.items() if k != "dataset_id"}
    with pytest.raises(
        census2_dataset_writer.store_format.SummaryRecordError,
        match="missing \\['dataset_id'\\]",
    ):
        census2_dataset_writer.store_format.validate_summary(record)


def test_read_summary_refuses_a_malformed_record(tmp_path: Path) -> None:
    dest, _ = written(tmp_path)
    root = zarr.open_group(dest, mode="r+")
    root.update_attributes(
        {"cellxgene-census2": {"summary": {**GOOD_RECORD, "version": 9}}}
    )
    with pytest.raises(
        census2_dataset_writer.store_format.SummaryRecordError, match="version 9"
    ):
        census2_dataset_writer.read_summary(dest)
    # An open root group reads the same way.
    with pytest.raises(
        census2_dataset_writer.store_format.SummaryRecordError, match="version 9"
    ):
        census2_dataset_writer.read_summary(zarr.open_group(dest, mode="r"))


def test_write_refuses_a_malformed_record_and_writes_nothing(tmp_path: Path) -> None:
    write: Any = census2_dataset_writer.write
    bad = {**GOOD_RECORD, "dataset_id": ""}
    with pytest.raises(
        census2_dataset_writer.store_format.SummaryRecordError, match="dataset_id"
    ):
        write(anndata.AnnData(X=np.zeros((2, 2))), tmp_path / "store.zarr", summary=bad)
    assert not (tmp_path / "store.zarr").exists()


def test_write_needs_a_summary_record(tmp_path: Path) -> None:
    write: Any = census2_dataset_writer.write
    with pytest.raises(TypeError, match="summary"):
        write(anndata.AnnData(X=np.zeros((2, 2))), tmp_path / "store.zarr")


def test_a_rejected_dataset_has_no_summary_record() -> None:
    rejected: census2_dataset_writer.Result = {
        "dataset_id": "d",
        "dataset_version_id": "v",
        "admitted": False,
        "reason": "no cells on the assay allowlist",
        "slots": [],
        "dropped_obs_columns": [],
        "dropped_var_features": [],
        "gene_reference_version": None,
    }
    with pytest.raises(ValueError, match="not admitted"):
        census2_dataset_writer.summary_record(rejected)


# --- slot policy ----------------------------------------------------------


def test_census_slot_policy_drops_extras(tmp_path: Path) -> None:
    dest, result = written(tmp_path, extras=True)
    store: Any = zarr.open_consolidated(dest, mode="r")

    assert result["slots"] == ["obs", "var", "X", "uns[2 keys]"]
    assert set(store["uns"]) == set(census2_dataset_writer.DEFAULTS["uns_keys"])
    for slot in census2_dataset_writer.OPTIONAL_SLOTS:
        assert slot not in store or len(list(store[slot])) == 0


def test_uns_keys_is_configurable(tmp_path: Path) -> None:
    """The carried uns keys are policy, not a fixed set in the code.

    The organism key stays, because the gene projection reads it.
    """
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)
    keys = ["organism_ontology_term_id", "title"]
    create(
        source,
        dest,
        config={"all_slots": False, "uns_keys": keys},
        refdata=tmp_path / "refdata",
    )

    store: Any = zarr.open_consolidated(dest, mode="r")
    assert set(store["uns"]) == set(keys)


def test_unwanted_uns_keys_are_never_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A spatial dataset holds a large image under uns; the policy drops it, so
    reading it would be paid for and thrown away."""
    source = tmp_path / "source.h5ad"
    adata = dataset()
    adata.uns["spatial"] = {"image": np.zeros((64, 64, 3), dtype="uint8")}
    adata.write_h5ad(source)

    seen: list[str] = []

    def spy(elem: Any, *args: Any, **kwargs: Any) -> Any:
        seen.append(getattr(elem, "name", "?"))
        return read_elem(elem, *args, **kwargs)

    # By name: read_elem is imported by the module, not exported from it.
    monkeypatch.setattr("census2_dataset_writer.read_elem", spy)
    with h5py.File(source) as f:
        built, _ = build(f, refdata=tmp_path / "refdata")

    assert built is not None
    assert "spatial" not in built.uns
    # Reading the group reads everything under it, so the group itself must
    # never be read: only the keys the policy keeps.
    assert "/uns" not in seen, seen
    assert "/uns/schema_version" in seen, seen


def test_all_slots_carries_the_cell_axis_slots(tmp_path: Path) -> None:
    """The gene projection drops the slots that index the source genes."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset(extras=True).write_h5ad(source)
    result = create(
        source, dest, config={"all_slots": True}, refdata=tmp_path / "refdata"
    )

    store: Any = zarr.open_consolidated(dest, mode="r")
    assert "obsm/X_pca" in result["slots"]
    assert "varm/gene_stuff" not in result["slots"]
    assert "layers/counts" not in result["slots"]
    # The full uns, so the key the Census set drops survives here.
    assert "title" in set(store["uns"])


def test_slots_report_what_harmonize_left(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The list describes the written object, not the one build assembled."""

    def rename_and_drop(
        adata: anndata.AnnData, config: Any = None, refdata: Any = None
    ) -> anndata.AnnData:
        assert adata.X is not None
        x: Any = adata.X  # the getter's union is wider than layers accepts
        adata.layers["harmonized"] = x
        del adata.raw
        return adata

    monkeypatch.setattr(census2_dataset_writer, "harmonize", rename_and_drop)
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)

    with h5py.File(source) as f:
        adata, result = build(f)

    assert adata is not None
    assert "layers/harmonized" in result["slots"]
    assert "raw" not in result["slots"]


# --- admission ------------------------------------------------------------


def test_not_admitted_writes_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A rejection is a normal outcome: no store, a reason, and no exception."""
    monkeypatch.setattr(
        census2_dataset_writer,
        "admit",
        lambda adata, config=None, refdata=None: "no thanks",
    )
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)

    result = create(source, dest)

    assert result == {
        **IDS,
        "admitted": False,
        "reason": "no thanks",
        "slots": [],
        "dropped_obs_columns": [],
        "dropped_var_features": [],
        "gene_reference_version": None,
    }
    assert not dest.exists()


def test_an_organism_without_a_reference_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    adata = dataset()
    adata.uns["organism_ontology_term_id"] = "NCBITaxon:7955"
    adata.write_h5ad(source)

    result = create(source, dest, refdata=tmp_path / "refdata")

    assert result["admitted"] is False
    assert "NCBITaxon:7955" in str(result["reason"])
    assert not dest.exists()


# --- the assay allowlist --------------------------------------------------

ON_LIST = "EFO:0009922"  # 10x 3' v3
OFF_LIST = "EFO:0030027"  # snmC-Seq2, a DNA methylation assay


def with_assay(adata: anndata.AnnData, *terms: str) -> anndata.AnnData:
    """Give each cell an assay term, cycling through `terms`."""
    n = adata.n_obs
    adata.obs["assay_ontology_term_id"] = pd.Categorical(
        [terms[i % len(terms)] for i in range(n)]
    )
    return adata


def test_off_list_assay_is_rejected(tmp_path: Path) -> None:
    """The criterion that keeps dense methylation matrices out of the workflow."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), OFF_LIST).write_h5ad(source)

    result = create(source, dest)

    assert result["admitted"] is False
    assert OFF_LIST in str(result["reason"])
    assert not dest.exists()


def test_one_on_list_cell_admits_the_dataset(tmp_path: Path) -> None:
    """Census 1 masks the off-list cells. This does not; it takes them all."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), ON_LIST, OFF_LIST, OFF_LIST).write_h5ad(source)

    result = create(source, dest, refdata=tmp_path / "refdata")

    assert result["admitted"] is True
    assert zarr.open_consolidated(str(dest))["obs"] is not None


def test_missing_assay_column_is_admitted(tmp_path: Path) -> None:
    """An AnnData without the column is not a Discover dataset, so the
    criterion does not apply. Every other test here relies on this."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)
    assert "assay_ontology_term_id" not in dataset().obs

    assert create(source, dest, refdata=tmp_path / "refdata")["admitted"] is True


def test_empty_allowlist_turns_the_criterion_off(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    with_assay(dataset(), OFF_LIST).write_h5ad(source)

    config = census2_dataset_writer.load_config()
    config["assay_allowlist"] = []

    result = create(source, dest, config=config, refdata=tmp_path / "refdata")
    assert result["admitted"] is True


# --- the obs projection ---------------------------------------------------

# Both names are taken from the corpus: a cell2location factor on a Visium
# dataset, and a clinical unit on a 10x one. The second holds three slashes.
SLASH_COLS = ["T/NK_cycling", "eGFR_(ml/min/1.73m2)"]


def with_obs_columns(adata: anndata.AnnData, *names: str) -> anndata.AnnData:
    for i, name in enumerate(names):
        adata.obs[name] = np.full(adata.n_obs, float(i), dtype="float32")
    return adata


def write_source(adata: anndata.AnnData, path: Path) -> None:
    """Write a fixture h5ad that holds a slash in an obs column name.

    anndata 0.12 added `disallow_forward_slash_in_h5ad`, on by default, so
    anndata refuses to write such a name even though HDF5 permits it. The
    corpus datasets that carry one predate the setting. Turn it off to build
    the fixture; the reader never needs it.
    """
    with anndata.settings.override(disallow_forward_slash_in_h5ad=False):
        adata.write_h5ad(path)


def test_project_obs_removes_only_the_slash_columns() -> None:
    adata = with_obs_columns(dataset(), *SLASH_COLS)
    kept = list(census2_dataset_writer.project_obs(adata).obs.columns)
    assert kept == ["cell_type"]


def test_project_obs_keeps_a_clean_frame_untouched() -> None:
    adata = dataset()
    before = list(adata.obs.columns)
    assert list(census2_dataset_writer.project_obs(adata).obs.columns) == before


def test_harmonize_projects_obs_and_genes(refdata: Path) -> None:
    """harmonize is the seam; project_obs, select_counts, and project_genes are
    what it does."""
    adata = with_obs_columns(dataset(), *SLASH_COLS)
    expected = projected(adata)
    harmonized = census2_dataset_writer.harmonize(adata, refdata=refdata)
    assert list(harmonized.obs.columns) == ["cell_type"]
    assert list(harmonized.var.index) == REF_IDS
    assert harmonized.raw is None
    x: Any = harmonized.X
    assert np.array_equal(x.toarray(), expected)


def test_a_slash_column_no_longer_fails_the_write(tmp_path: Path) -> None:
    """Zarr gives `/` its path meaning, so without `project_obs` this write fails."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    write_source(with_obs_columns(dataset(), *SLASH_COLS), source)

    result = create(source, dest, refdata=tmp_path / "refdata")

    assert result["admitted"] is True
    assert result["dropped_obs_columns"] == SLASH_COLS  # source order
    store = zarr.open_consolidated(str(dest))
    obs: Any = read_elem(store["obs"])
    assert list(obs.columns) == ["dataset_id", "cell_type"]


def test_dropped_obs_columns_is_empty_when_nothing_is_dropped(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    with h5py.File(source) as f:
        _, result = build(f, refdata=tmp_path / "refdata")
    assert result["dropped_obs_columns"] == []
    assert result["dropped_var_features"] == []


def test_a_rejected_dataset_reports_no_dropped_columns(tmp_path: Path) -> None:
    """A rejection short-circuits before harmonize, so there is nothing to drop."""
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    write_source(with_assay(with_obs_columns(dataset(), *SLASH_COLS), OFF_LIST), source)

    result = create(source, dest)

    assert result["admitted"] is False
    assert result["dropped_obs_columns"] == []


def test_default_allowlist_holds_rna_and_spatial() -> None:
    """The default is Census 1's two lists, merged."""
    allowlist = set(census2_dataset_writer.DEFAULTS["assay_allowlist"])
    assert ON_LIST in allowlist
    assert OFF_LIST not in allowlist
    assert "EFO:0022857" in allowlist  # Visium, from ALLOWED_SPATIAL_ASSAYS
    assert allowlist == set(census2_dataset_writer.assays.DEFAULT_ALLOWLIST)


# --- the write thread bound -----------------------------------------------


def test_threads_bounds_the_write_and_restores_the_setting(tmp_path: Path) -> None:
    """dask counts the host's CPUs unless told otherwise, and a container with
    a CPU share sets no cgroup quota for it to read."""
    import dask

    before = dask.config.get("num_workers", None)
    source = tmp_path / "source.h5ad"
    dest = tmp_path / "store.zarr"
    dataset().write_h5ad(source)

    result = create(source, dest, threads=2, refdata=tmp_path / "refdata")
    assert result["admitted"] is True
    assert dask.config.get("num_workers", None) == before


# --- the memory contract --------------------------------------------------


def test_build_leaves_the_matrix_lazy(tmp_path: Path) -> None:
    """`build` must not materialize X. A large dataset does not fit in memory."""
    import dask.array as da

    source = tmp_path / "source.h5ad"
    dataset(400, 20).write_h5ad(source)

    with h5py.File(source) as f:
        adata, _ = build(f, refdata=tmp_path / "refdata")
        assert adata is not None
        # The projected X too: blocked by rows, whole reference rows per block.
        assert isinstance(adata.X, da.Array), type(adata.X)
        assert adata.X.chunks[1] == (len(REF_IDS),)


def test_row_chunk_is_bounded_by_the_row_count(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset(40, 12).write_h5ad(source)
    with h5py.File(source) as f:
        assert census2_dataset_writer.row_chunk(f["X"], 40) == 40


# --- provenance -----------------------------------------------------------


def test_provenance_reports_the_code_that_wrote_a_store() -> None:
    p = census2_dataset_writer.provenance()
    assert p["recipe"] == census2_dataset_writer.RECIPE
    assert p["zarr_format"] == census2_dataset_writer.ZARR_FORMAT
    # Reported from the installed distributions.
    assert p["census2_dataset_writer_version"] == version("census2-dataset-writer")
    assert p["anndata_version"] == version("anndata")
    assert p["zarr_version"] == version("zarr")


def test_provenance_follows_the_slot_policy() -> None:
    assert census2_dataset_writer.provenance()["slot_policy"] == "census"
    assert (
        census2_dataset_writer.provenance({"all_slots": True})["slot_policy"] == "all"
    )


# --- destinations ---------------------------------------------------------


def test_is_remote_covers_every_scheme() -> None:
    for scheme in census2_dataset_writer.REMOTE_SCHEMES:
        assert census2_dataset_writer.is_remote(f"{scheme}bucket/key")
    assert not census2_dataset_writer.is_remote("/tmp/local.zarr")
    assert not census2_dataset_writer.is_remote(Path("/tmp/local.zarr"))
    assert not census2_dataset_writer.is_remote(None)


def test_as_store_passes_a_local_path_through(tmp_path: Path) -> None:
    assert census2_dataset_writer.as_store(tmp_path) is tmp_path


def test_as_store_adapts_a_remote_url() -> None:
    """A bare s3:// string would make zarr fall back to fsspec and fail."""
    pytest.importorskip("obstore")
    store = census2_dataset_writer.as_store("s3://bucket/prefix.zarr")
    assert isinstance(store, zarr.storage.ObjectStore)


# --- the dataset_id column ------------------------------------------------------


def test_store_obs_names_its_dataset(tmp_path: Path) -> None:
    dest, _ = written(tmp_path)
    obs = cast("pd.DataFrame", read_elem(zarr.open_consolidated(dest, mode="r")["obs"]))
    assert obs.columns[0] == "dataset_id"
    assert isinstance(obs["dataset_id"].dtype, pd.CategoricalDtype)
    assert list(obs["dataset_id"].cat.categories) == [IDS["dataset_id"]]
    assert (obs["dataset_id"] == IDS["dataset_id"]).all()


def test_a_source_dataset_id_column_fails_the_build(tmp_path: Path) -> None:
    """The name belongs to the store, so the writer refuses to guess."""
    source = tmp_path / "source.h5ad"
    adata = dataset()
    adata.obs["dataset_id"] = "curator-label"
    adata.write_h5ad(source)
    with pytest.raises(
        ValueError, match="unexpected dataset_id column in the source obs"
    ):
        create(source, tmp_path / "store.zarr", refdata=tmp_path / "refdata")
    assert not (tmp_path / "store.zarr").exists()


def test_add_dataset_id_refuses_an_existing_column() -> None:
    adata = anndata.AnnData(
        X=np.zeros((2, 2)), obs=pd.DataFrame({"dataset_id": ["a", "b"]})
    )
    with pytest.raises(ValueError, match="already has a dataset_id column"):
        census2_dataset_writer.add_dataset_id(adata, "d")


def test_build_returns_the_object_that_the_store_holds(tmp_path: Path) -> None:
    source = tmp_path / "source.h5ad"
    dataset().write_h5ad(source)
    with h5py.File(source) as f:
        adata, result = build(f, refdata=tmp_path / "refdata")
    assert adata is not None
    obs = cast("pd.DataFrame", adata.obs)
    assert obs.columns[0] == "dataset_id"
    assert set(obs["dataset_id"]) == {IDS["dataset_id"]}
    assert (result["dataset_id"], result["dataset_version_id"]) == (
        IDS["dataset_id"],
        IDS["dataset_version_id"],
    )


def test_write_refuses_obs_for_another_dataset(tmp_path: Path) -> None:
    adata = anndata.AnnData(
        X=np.zeros((2, 2)), obs=pd.DataFrame({"dataset_id": ["other", "other"]})
    )
    write: Any = census2_dataset_writer.write
    with pytest.raises(ValueError, match="must hold only the summary record's d"):
        write(adata, tmp_path / "store.zarr", summary=dict(GOOD_RECORD))
    assert not (tmp_path / "store.zarr").exists()


# --- the store format module -------------------------------------------------


def test_store_format_imports_only_the_standard_library() -> None:
    """A reader can take the store format without the writer's dependencies.

    The test loads the module file alone, in a new interpreter, so that the
    package's own imports do not count.
    """
    import subprocess
    import sys

    path = Path(census2_dataset_writer.store_format.__file__)
    code = (
        "import importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('store_format', {str(path)!r})\n"
        "module = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(module)\n"
        "heavy = {'anndata', 'numpy', 'pandas', 'zarr', 'dask', 'h5py', 'scipy'}\n"
        "print(sorted(heavy & set(sys.modules)))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]"


def test_gene_references_use_the_store_format_labels() -> None:
    labels = census2_dataset_writer.store_format.GENE_REFERENCE_VERSIONS
    references = census2_dataset_writer.genes.REFERENCES
    assert {taxon: spec["version"] for taxon, spec in references.items()} == labels
