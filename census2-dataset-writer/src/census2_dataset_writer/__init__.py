"""Write one Census 2 store from one Discover dataset.

The writer reads one published Discover `.h5ad` and writes one AnnData-Zarr
store. The bulk ingestion pipeline (`bulk-ingest-nf`) calls this module once
for each dataset, and a developer can call it from Python. The module reads
one dataset and not the corpus, so a step that needs a corpus-wide
aggregate cannot run here.

Public calls
------------
`create(source, dest, *, dataset_id, dataset_version_id, ...)` builds the
store and writes it. It writes nothing when the dataset is not admitted, and
returns a `Result` in both cases. A Discover h5ad does not hold its own
identifiers, so the caller gives them. `build` and `write` split the same work
into two calls, so a caller can inspect the AnnData before the write:

    with h5py.File("dataset.h5ad") as f:
        adata, result = census2_dataset_writer.build(
            f, dataset_id=dataset_id, dataset_version_id=dataset_version_id
        )
        if adata is not None:
            summary = census2_dataset_writer.summary_record(result)
            census2_dataset_writer.write(adata, "store.zarr", summary=summary)

Keep `f` open until `write` returns, because `X` is a dask graph until then
(see Memory). `create` does this itself when it receives a path.

`slots_of(adata)` lists what a store written from an AnnData will hold.
`read_summary(dest)` returns the summary record of a store. `load_config(path)`
returns the configuration mapping. `provenance(config)` returns the part of a
receipt that this module knows. `as_store(dest)` adapts a remote URL for zarr,
which ships no S3 client. `census2_dataset_writer.store_format` names what a
store holds beyond AnnData-Zarr, for the writer and for a reader of a store.

Decisions
---------
`build` runs four steps in order: assemble the core slots (`obs`, `var`,
`X`, `raw`, the carried `uns` keys), admit and filter the cells, apply the
slot policy, harmonize.
It then adds the `dataset_id` column to `obs`. A source `obs` that already has
a `dataset_id` column is an error, because that name belongs to the store.

1. Admission (`admit`) decides whether the dataset belongs in Census 2. A
   rejection is a normal result: `build` returns no AnnData and a reason,
   and the caller records the reason and succeeds. This decision is largely
   open. The stub applies one criterion: the organism must have a gene
   reference, so that harmonization can project the genes. `cell_mask` then
   applies the Census 1 cell filters: tissue type, assay, and primary data
   on a spatial assay. The store keeps only the cells that pass. A dataset
   with no cell left is rejected with a reason that names the filters. The
   result counts the cells that each filter drops. The assay filter also
   keeps the workflow away from modalities a Census is unlikely to carry.
2. The slot policy decides which AnnData slots the store carries. This
   decision is made: the set that the Census 1 builder reads (`obs`, `var`,
   `X`, `raw`, and the `uns` keys that the `uns_keys` configuration lists),
   or every slot when `all_slots` is set. Harmonization then puts the
   counts from `raw` in `X`, and removes the slots that index the source
   genes.
3. Harmonization (`harmonize`) changes values and annotations so that a
   query can combine stores. It runs five steps in order. `project_obs`
   projects `obs` onto `OBS_COLUMNS`, the CxG columns that Census 1 copies,
   each with a fixed type; it drops every other source column.
   `add_tissue_general` adds the Census 1 tissue roll-up (see
   `census2_dataset_writer.tissues`). `select_counts` puts the counts in `X`:
   `raw.X` when the dataset has `raw`, else `X`. The store has no `raw`.
   `add_obs_stats` adds the Census 1 statistics of the counts to `obs`.
   `project_genes` projects the counts onto the pinned gene reference of the
   dataset's organism. Thus every store of one organism has the same genes in
   `var`, in the same order. `var` also marks the dataset's assayed genes in
   `MEASURED_COLUMN`. See `census2_dataset_writer.genes`.

Admission runs before the optional slots are read, so a rejected dataset
does not read them. The `uns` policy is the exception: it applies at the read,
so `admit` sees the carried keys and not the whole group. A criterion that
needs another key adds it to `uns_keys`. `admit` and `harmonize` both receive
an AnnData, so a step author works with one kind of object.

Configuration and reference data
--------------------------------
`config` is a small mapping of policy: thresholds, criteria, the names of
lists. It lives in version control next to a release. `refdata` is a
directory of bulk reference data that the caller staged before the run:
ontology term lists, gene references, tissue roll-ups. The configuration
names a reference file, and a step resolves the name against `refdata`.
This module does not fetch a file at run time, because a HealthOmics task
has no public internet access.

An admitted dataset needs `refdata`. `gene_references` names one gene file
for each organism, with its SHA-256 digest. `project_genes` refuses a file
that is absent or that has another digest.

`load_config` is a stub: no file format is defined, and it returns
`DEFAULTS`. Do not define a configuration schema before the step that needs
it exists.

Summary record
--------------
Every store holds a summary record: a small JSON object in the attributes of
the store's root Zarr group. It holds the facts about the store that a catalog
needs and that the store does not otherwise hold: the dataset's identifiers and
the version of the gene reference. `store_format` defines the record, where it
is, and its check.
`build` also adds `dataset_id` to `obs`, as the first column, so that a cell
keeps its dataset when stores are concatenated.

`write` adds the record after anndata writes the store and its consolidated
metadata, so the record is the last write.

Receipt
-------
A receipt is the JSON a pipeline writes beside a store, about one run of the
pipeline. `Receipt` declares its shape, and `provenance` fills the part this
package knows; see their docstrings. A receipt holds a copy of the store's
summary record, so that a report on a run needs to open no store. Nothing
outside the pipeline reads a receipt, because the facts about a store are in
the store.

Memory
------
`build` returns a complete AnnData but reads only `obs`, `var`, and `uns`
into memory. The source `X` and `raw.X` are dask arrays of row blocks of
about `BLOCK_BYTES`, sparse when the source is sparse. The projected `X` is
a dask array of sparse blocks over the same rows. Each block holds at most
the stored values of its source block. Memory use follows the block
size and not the file size; the corpus p99 is about 26 GiB and the maximum
is 50 GiB. Two rules follow:

- An `h5py.File` given to `build` must stay open until `write` returns.
- A step that needs the whole matrix at one time forces a compute, and so
  does a mask computed from `X`. That compute can fail on a large dataset.
  Row-wise work maps over blocks, and a subset of the AnnData keeps `obs`,
  `var`, `X`, and the other slots aligned.

`add_obs_stats` maps over the blocks and computes one row of statistics for
each cell, so `build` reads `X` once and the write reads it again. The
statistics are `obs` columns, and `write_zarr` writes `obs` and computes `X`
in one call, so the statistics must exist before the write starts. One
shared pass needs the writer to write `X` itself, so that one dask compute
feeds the store and the statistics. A load of the whole `X` can fail, as the
rules above state. Revisit this choice when a third pass over `X` appears or
the ingest cost matters. The per-gene nnz file of tracker 4.14 can share the
pass of `add_obs_stats`.

`README.md` in this directory covers installation, use, and development.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import os
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, TypedDict, assert_never, cast

import anndata
import dask
import dask.array as da
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem
from zarr.abc.store import Store

from . import assays, genes, tissues
from . import store_format as store_format  # a public module; callers name it

# The store format, re-exported, so a caller writes census2_dataset_writer.X.
from .store_format import DATASET_ID_COLUMN as DATASET_ID_COLUMN
from .store_format import DERIVED_OBS_COLUMNS as DERIVED_OBS_COLUMNS
from .store_format import MEASURED_COLUMN as MEASURED_COLUMN
from .store_format import OBS_COLUMNS as OBS_COLUMNS
from .store_format import SUMMARY_KEY as SUMMARY_KEY
from .store_format import SUMMARY_NAMESPACE as SUMMARY_NAMESPACE
from .store_format import SUMMARY_VERSION as SUMMARY_VERSION
from .store_format import ColType as ColType
from .store_format import SummaryRecord as SummaryRecord
from .store_format import SummaryRecordError as SummaryRecordError
from .store_format import validate_summary as validate_summary

if TYPE_CHECKING:  # these are only ever types here, never constructed
    import numpy.typing as npt

# What the public functions accept. The prose above describes these unions; the
# aliases are the same statement in a form a type checker and an editor can use.
PathLike: TypeAlias = str | os.PathLike[str]
StoreLike: TypeAlias = "PathLike | Store"
# h5py ships no type information, so anything reached through it is Any. The
# names still record what a caller has to pass:
#   H5Matrix  a group of data/indices/indptr when sparse, a dataset when dense
#   Source    a path to an h5ad, or an already open h5py.File
H5Matrix: TypeAlias = Any
Source: TypeAlias = Any
# Configuration is open: its schema is not decided, so the values stay Any.
Config: TypeAlias = dict[str, Any]


class Result(TypedDict):
    """What `build` and `create` report about one dataset.

    A dict, so that it serializes into a receipt unchanged. It is declared, so
    that a read of an undeclared key fails the type check.
    """

    # The identifiers that the caller gave `build`. A Discover h5ad holds none.
    dataset_id: str
    dataset_version_id: str
    admitted: bool
    reason: str | None
    slots: list[str]
    # The cells that each filter drops (see `cell_mask`). Empty when `admit` rejects.
    dropped_cells: dict[str, int]
    # The `obs` columns that harmonization removed, in source order. Empty for
    # a rejected dataset, which is never harmonized.
    dropped_obs_columns: list[str]
    # The features of the source counts matrix that are not in the store's
    # `var`, in source order: the non-gene features and the genes of another
    # organism. Empty for a rejected dataset.
    dropped_var_features: list[str]
    # The version of the gene reference in the store's `var`. None for a
    # rejected dataset.
    gene_reference_version: str | None


# Slots that only `all_slots` carries. The Census builder reads none of them.
OPTIONAL_SLOTS = ("obsm", "varm", "obsp", "varp", "layers")

# Target bytes per row block of the expression matrix.
BLOCK_BYTES = 64 * 1024**2

# The recipe names the store format: the anndata version, the Zarr format, the X
# encoding, the slot policy, and harmonization. Raise it when any of those
# change, so that a receipt says which format its store holds. Nothing is keyed
# on it: a store's location is its dataset version alone. v2 adds the fixed `obs`
# schema, the Census 1 cell filters, `tissue_general`, the count statistics, and
# the `measured` column of `var`.
RECIPE = "v2"

# anndata 0.13 writes Zarr v3 by default. Recorded, so that a reader does not
# have to infer it.
ZARR_FORMAT = 3


def row_chunk(h5_obj: H5Matrix, n_rows: int) -> int:
    """Rows per dask block, sized so one block is about BLOCK_BYTES."""
    if isinstance(h5_obj, h5py.Group):  # sparse: data + indices per stored value
        nnz = int(h5_obj["data"].shape[0])
        per_value = int(h5_obj["data"].dtype.itemsize) + int(
            h5_obj["indices"].dtype.itemsize
        )
        row_bytes = max(1, (nnz // max(n_rows, 1)) * per_value)
    else:  # dense
        row_bytes = max(1, int(np.prod(h5_obj.shape[1:])) * int(h5_obj.dtype.itemsize))
    return max(1, min(n_rows, BLOCK_BYTES // row_bytes))


def lazy_matrix(h5_obj: H5Matrix, n_rows: int) -> da.Array:
    """A dask array of row blocks. The blocks are sparse if the source is."""
    # `read_elem_lazy` is typed for every element kind it can return. For a
    # matrix it returns a dask array, and the memory contract depends on that,
    # so state it here rather than widening the signature.
    return cast(
        "da.Array", read_elem_lazy(h5_obj, chunks=(row_chunk(h5_obj, n_rows), -1))
    )


# Policy that is already decided. Anything a caller can choose about what a store
# contains belongs here, not in a separate function argument.
DEFAULTS: Config = {
    # Carry every AnnData slot instead of the Census set. See build().
    "all_slots": False,
    # The uns keys a store carries, when `all_slots` is false. The default is
    # the set the CELLxGENE Census builder reads: see AnnDataProxy in the
    # chanzuckerberg/cellxgene-census repository, at
    # tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/anndata.py
    "uns_keys": ["schema_version", "organism_ontology_term_id"],
    # The assay ontology terms whose cells a store keeps; see `cell_mask`. An
    # empty list turns the assay filter off. A sorted list rather than the
    # frozenset, because a config file is JSON.
    # See census2_dataset_writer.assays for where the terms come from.
    "assay_allowlist": sorted(assays.DEFAULT_ALLOWLIST),
    # Census 1 keeps these `tissue_type` values and drops cell culture.
    "tissue_types": ["organoid", "tissue"],
    # Census 1 requires `is_primary_data` on these assays, in its spatial experiments.
    "primary_data_assays": sorted(assays.SPATIAL),
    # The gene reference of each organism, keyed by organism ontology term: the
    # file name in `refdata`, the annotation version, the SHA-256 digest, and
    # the URL to stage it from. An organism with no entry is not admitted.
    # See census2_dataset_writer.genes for the pin.
    "gene_references": {taxon: dict(spec) for taxon, spec in genes.REFERENCES.items()},
}


def _gene_references(config: Config) -> dict[str, genes.Reference]:
    """The `gene_references` of `config`, or the defaults when it names none."""
    # A config value is Any by design, so name its shape where it is read.
    return cast(
        "dict[str, genes.Reference]",
        config.get("gene_references", DEFAULTS["gene_references"]),
    )


def load_config(path: PathLike | None = None) -> Config:
    """Read ingestion configuration from `path`, over `DEFAULTS`.

    STUB: ignores `path` and returns the defaults. No file format is defined yet.

    Keep this the only place that parses configuration, so the format stays one
    decision. See the module docstring for how config and reference data differ
    and why reference data cannot be fetched at run time.
    """
    return dict(DEFAULTS)


class StoreProvenance(TypedDict):
    """What produced a store. Always present in a receipt, because this package
    supplies all of it.
    """

    recipe: str
    slot_policy: Literal["all", "census"]
    zarr_format: int
    census2_dataset_writer_version: str
    anndata_version: str
    zarr_version: str
    # The package version decides the `tissue_general` mapping.
    cellxgene_ontology_guide_version: str


class Receipt(StoreProvenance, total=False):
    """The JSON a pipeline writes beside a store, about one run.

    The pipeline writes it, and its report step reads it. The shape lives here
    because the script that writes one and the script that reads one cannot
    import each other, and because `provenance` fills the part that only this
    package knows. For an admitted dataset, `summary` is a copy of the store's
    summary record, so that the report reads no store.

    The `StoreProvenance` keys are always present. The rest is what a writer
    knows about where the dataset came from and what it did, so it is optional
    here. `bulk-ingest-nf` fills all of them for an admitted dataset. A rejected
    dataset gets no `slots_written`, `summary`, `store`, `output_bytes`, or
    `published_by`, and `check` appears only when the task runs with `--check`.
    """

    dataset_id: str
    dataset_version_id: str
    source_uri: str
    source_url: str | None
    source_filesize: int
    cell_count: int | None
    schema_version: str | None
    admitted: bool
    reason: str | None
    slots_written: list[str]
    dropped_cells: dict[str, int]
    dropped_obs_columns: list[str]
    dropped_var_features: list[str]
    summary: SummaryRecord
    store: str
    output_bytes: int
    published_by: str
    elapsed_sec: float
    check: dict[str, Any]
    config_file: str | None
    config_sha256: str | None
    config_effective: Config
    refdata: str | None


def provenance(config: Config | None = None) -> StoreProvenance:
    """The store-level facts for a receipt: what wrote a store, and to what format.

    Read from the installed distributions rather than passed in, so that a
    receipt cannot disagree with the code that produced it.
    """
    config = DEFAULTS if config is None else config
    return {
        "recipe": RECIPE,
        "slot_policy": "all" if config.get("all_slots", False) else "census",
        "zarr_format": ZARR_FORMAT,
        "census2_dataset_writer_version": importlib.metadata.version(
            "census2-dataset-writer"
        ),
        "anndata_version": importlib.metadata.version("anndata"),
        "zarr_version": importlib.metadata.version("zarr"),
        "cellxgene_ontology_guide_version": importlib.metadata.version(
            "cellxgene-ontology-guide"
        ),
    }


def admit(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> str | None:
    """Decide whether Census 2 admits this dataset.

    THIS IS A MINIMAL STUB. It applies one criterion: the organism must have
    an entry in `gene_references`, so that `project_genes` can project the
    dataset. It is not production-ready admission logic. `admit` holds the
    criteria on the whole dataset, and `cell_mask` decides which cells a store
    keeps, assays included. Real admission criteria are open.

    Return `None` to admit the dataset. Return a short reason to reject it.
    A rejection is a normal result. The caller records the reason and exits
    with success. Do not raise an exception for a dataset that you reject.
    Raise only for an input that this function cannot read.

    `adata` holds `obs`, `var`, `X`, `raw`, and the `uns` keys that the slot
    policy carries. A criterion that needs another key adds it to `uns_keys`.
    `X` and `raw.X` are not loaded. A read of `X` or `raw.X` loads the full
    matrix.

    `refdata` is the directory that the caller staged. The configuration
    names a reference file. It does not give the path. Resolve each name
    against `refdata`.

    A `uns` without `organism_ontology_term_id` is admitted. Such an object is
    not a Discover dataset, so a Discover criterion does not apply to it.
    `project_genes` cannot project it, and raises an error.
    """
    config = DEFAULTS if config is None else config
    references = _gene_references(config)
    organism = adata.uns.get("organism_ontology_term_id")
    if organism is not None and str(organism) not in references:
        return f"no gene reference for organism {organism}"
    return None


# A cell that fails several filters counts under the first one in this order.
CELL_FILTERS = {
    "tissue_type": ("tissue_type",),
    "assay": ("assay_ontology_term_id",),
    "primary_data": ("assay_ontology_term_id", "is_primary_data"),
}


def cell_mask(
    obs: pd.DataFrame, config: Config | None = None
) -> tuple[npt.NDArray[np.bool_], dict[str, int]]:
    """The cells that a store keeps, and the number of cells that each filter drops.

    These are the Census 1 cell filters (its `make_anndata_cell_filter`). A
    cell passes when its `tissue_type` is in `tissue_types`, when its assay is
    in `assay_allowlist`, and, for an assay in `primary_data_assays`, when
    `is_primary_data` is true. Census 1 writes RNA and spatial cells to
    different experiments and requires primary data in the spatial ones.
    Census 2 merges the two assay lists, so it applies that rule to each cell.
    An empty `assay_allowlist` turns the assay filter off.

    The counts follow `CELL_FILTERS`. A cell counts under the first filter that
    it fails, so the counts add up to the cells dropped.

    Raises `ValueError` when `obs` lacks a column that a filter reads. The CxG
    schema requires each of them.
    """
    config = DEFAULTS if config is None else config
    missing = sorted({c for cols in CELL_FILTERS.values() for c in cols} - set(obs))
    if missing:
        raise ValueError(
            f"the source obs lacks the columns {missing} of the cell filters"
        )
    tissue_types = cast(
        "list[str]", config.get("tissue_types", DEFAULTS["tissue_types"])
    )
    allowlist = cast("list[str]", config.get("assay_allowlist") or [])
    primary_assays = cast(
        "list[str]", config.get("primary_data_assays", DEFAULTS["primary_data_assays"])
    )
    assay = obs["assay_ontology_term_id"]
    passes = {
        "tissue_type": obs["tissue_type"].isin(tissue_types).to_numpy(),
        "assay": assay.isin(allowlist).to_numpy()
        if allowlist
        else np.ones(len(obs), bool),
        "primary_data": ~assay.isin(primary_assays).to_numpy()
        | obs["is_primary_data"].astype(bool).to_numpy(),
    }
    keep = np.ones(len(obs), dtype=bool)
    dropped: dict[str, int] = {}
    for name in CELL_FILTERS:
        dropped[name] = int((keep & ~passes[name]).sum())
        keep &= passes[name]
    return keep, dropped


def _no_cells_reason(obs: pd.DataFrame, dropped: dict[str, int]) -> str:
    """The rejection reason of a dataset that loses every cell to `cell_mask`.

    It names the filters that dropped cells and, for the assay filter, the
    assays present. It holds no counts, because `dropped_cells` holds them and
    a run report groups rejections by their exact text.
    """
    parts = []
    for name in (n for n in CELL_FILTERS if dropped.get(n)):
        if name == "assay":
            present = sorted({str(t) for t in obs["assay_ontology_term_id"].unique()})
            parts.append(f"assay ({', '.join(present[:3])})")
        else:
            parts.append(name)
    return f"no cells pass the cell filters: {', '.join(parts) or 'no cells'}"


def project_obs(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Project `obs` onto the columns of `OBS_COLUMNS`.

    Census 1 projects `obs` onto a closed schema, with fixed types: the columns
    that the CxG schema requires in every dataset, plus columns that it adds.
    This step copies the CxG columns, every column of `OBS_COLUMNS` except
    `DERIVED_OBS_COLUMNS`, in schema order, and gives each one its type: a
    categorical with the sorted observed values as string categories, NumPy
    `bool` for `is_primary_data`, and plain strings for `observation_joinid`.
    The source index stays. Every other source column goes: the result reports
    them in `dropped_obs_columns`.

    It raises `ValueError` when the source lacks a schema column, when a column
    holds a null, when `is_primary_data` is not boolean, or when
    `observation_joinid` is not unique. The CxG schema requires each of these,
    so the writer does not choose a repair.

    `add_tissue_general` adds the `tissue_general` pair. `build` adds
    `dataset_id` after harmonization.

    This step consumes `adata`, as `harmonize` states.
    """
    # `AnnData.obs` is typed `DataFrame | Dataset2D`, the second being what
    # `read_lazy` produces. `build` always assembles obs with `read_elem`, so
    # narrow to the pandas frame here rather than widen the signature.
    obs = cast("pd.DataFrame", adata.obs)
    columns = [
        (name, col_type)
        for name, col_type in OBS_COLUMNS
        if name not in DERIVED_OBS_COLUMNS
    ]
    missing = [name for name, _ in columns if name not in obs.columns]
    if missing:
        raise ValueError(f"the source obs lacks the schema columns {missing}")
    adata.obs = pd.DataFrame(
        {
            name: _project_obs_column(obs[name], name, col_type)
            for name, col_type in columns
        },
        index=obs.index,
    )
    return adata


def _project_obs_column(column: pd.Series, name: str, col_type: ColType) -> pd.Series:
    """One `obs` column with the type that `OBS_COLUMNS` gives it."""
    if bool(column.isna().any()):
        raise ValueError(f"obs {name} holds a null")
    if col_type is ColType.CATEGORY:
        if isinstance(column.dtype, pd.CategoricalDtype):
            # Work on the categories, not on one value per cell.
            column = column.cat.remove_unused_categories()
            names = column.cat.categories.astype(str)
            if not names.is_unique:
                raise ValueError(f"obs {name} has categories that are equal as strings")
            column = column.cat.rename_categories(list(names))
        else:
            column = column.astype(str).astype("category")
        categories = sorted(column.cat.categories)
        return column.cat.reorder_categories(categories).cat.as_unordered()
    if col_type is ColType.BOOL:
        if column.dtype == bool:
            return column
        if str(column.dtype) == "boolean":
            return column.astype(bool)
        raise ValueError(f"obs {name} is {column.dtype}, not boolean")
    if col_type is ColType.STRING:
        values = pd.Series(
            column.astype(str).to_numpy(dtype=object), index=column.index
        )
        _require_unique(values, name)
        return values
    if col_type is ColType.INT64 or col_type is ColType.FLOAT64:
        return column.astype(col_type.value)
    assert_never(col_type)


def add_tissue_general(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Add the `tissue_general` pair to `obs`, after `tissue_type`.

    Each cell's `tissue_ontology_term_id` maps to a high-level tissue term and
    its label, as in Census 1; see `census2_dataset_writer.tissues`. Run
    `project_obs` first. Raises `ValueError` for a tissue term that the ontology
    does not hold. This step consumes `adata`, as `harmonize` states.
    """
    obs = cast("pd.DataFrame", adata.obs)
    terms = obs["tissue_ontology_term_id"].astype(str)
    mapping = tissues.tissue_general(terms.unique())
    labels = {term: label for term, (_, label) in mapping.items()}
    general_ids = {term: general for term, (general, _) in mapping.items()}
    obs["tissue_general"] = _project_obs_column(
        terms.map(labels), "tissue_general", ColType.CATEGORY
    )
    obs["tissue_general_ontology_term_id"] = _project_obs_column(
        terms.map(general_ids), "tissue_general_ontology_term_id", ColType.CATEGORY
    )
    return adata


def _require_unique(values: pd.Series, name: str) -> None:
    """Raise `ValueError` for a null in `obs` column `name`, or name its repeated values."""
    if bool(values.isna().any()):
        raise ValueError(f"obs {name} holds a null")
    repeated = sorted({str(v) for v in values[values.duplicated()]})
    if repeated:
        raise ValueError(
            f"obs {name} is not unique: {len(repeated)} values repeat, "
            f"such as {', '.join(repeated[:5])}"
        )


def select_counts(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Put the counts in `X`, and remove `raw`.

    A Discover dataset holds the counts in `raw.X` when `X` holds normalized
    values, and in `X` when it has no `raw`. With `raw`, return a new AnnData
    whose `X` and `var` are `raw.X` and `raw.var`. `layers`, `varm`, and `varp`
    index the genes of the normalized `X`, so they are dropped. Without `raw`,
    return `adata` unchanged. `X` stays lazy. This step consumes `adata`, as
    `harmonize` states.
    """
    if adata.raw is None:
        return adata
    return genes.with_gene_axis(adata, adata.raw.X, adata.raw.var)


def project_genes(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Project the counts in `X` onto the organism's pinned gene reference.

    Return a new AnnData. `X` has its columns moved to reference positions, and
    `var` is the reference. `varm`, `varp`, and `layers` index the source genes,
    so they are dropped. `X` stays lazy. Run `select_counts` first: this step
    refuses an AnnData that has `raw`. See `census2_dataset_writer.genes` for the
    rules.

    `gene_references` in `config` names the reference file for each organism.
    The file must be staged in `refdata`. This step consumes `adata`, as
    `harmonize` states.
    """
    config = DEFAULTS if config is None else config
    return genes.project(adata, _gene_references(config), refdata)


_ROW_STATS = ("raw_sum", "nnz", "raw_mean_nnz", "raw_variance_nnz")


def _row_stats(block: Any) -> npt.NDArray[np.float64]:
    """The Census 1 statistics of each row of one block, in `_ROW_STATS` order."""
    x = sp.csr_matrix(block, copy=True)
    # As in Census 1, merge duplicates (this sorts each row), then drop zeros.
    x.sum_duplicates()
    x.eliminate_zeros()
    nnz = np.diff(x.indptr)
    raw_sum = np.asarray(x.sum(axis=1, dtype=np.float64)).ravel()
    rows = np.repeat(np.arange(len(nnz)), nnz)
    # Census 1's `_var_ndarray` shifts by the row's first value, squares in float32,
    # and sums in float64 in storage order. bincount keeps that order, so the variance
    # is bit-identical. Counts outside the CxG float32 type keep their precision.
    data = x.data if x.data.dtype == np.float32 else x.data.astype(np.float64)
    shifted = data - data[x.indptr[rows]]
    ex = np.bincount(rows, weights=shifted, minlength=len(nnz))
    ex2 = np.bincount(rows, weights=shifted * shifted, minlength=len(nnz))
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = raw_sum / nnz
        variance = (ex2 - ex**2 / nnz) / (nnz - 1)
    mean[nnz == 0] = 0.0
    variance[nnz < 2] = 0.0
    return np.column_stack([raw_sum, nnz, mean, variance])


def add_obs_stats(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Add the five Census 1 statistics of the counts to `obs`, after `tissue_general`.

    `raw_sum`, `nnz`, `raw_mean_nnz`, and `raw_variance_nnz` describe the stored
    values of each cell, as Census 1 computes them (its `get_obs_stats`, after it
    removes explicit zeros). `n_measured_vars` is the number of genes of the
    organism in the source, the same for every cell; `project_genes` marks the
    same genes in `MEASURED_COLUMN`.

    Run `select_counts` first, and `project_genes` after. The step reads the
    organism's genes in source order, the matrix that Census 1 reads, so that
    the variance sums in the same order. It reads `X` once, one row block at a
    time, and the write reads it again. This step consumes `adata`, as
    `harmonize` states.
    """
    taxon = str(adata.uns.get("organism_ontology_term_id"))
    kept = genes.organism_genes(cast("pd.DataFrame", adata.var), taxon)
    # `AnnData.X` is typed as any matrix or None. Here it is dask or a matrix in memory.
    counts: Any = adata.X
    x = counts[:, kept]
    if isinstance(x, da.Array):
        # dask types `map_blocks` for same-shape blocks. This gives 4 values a cell.
        map_blocks: Any = x.map_blocks
        stats = np.asarray(
            map_blocks(
                _row_stats,
                chunks=(x.chunks[0], (len(_ROW_STATS),)),
                dtype=np.float64,
                meta=np.empty((0, len(_ROW_STATS)), dtype=np.float64),
            ).compute()
        )
    else:
        stats = _row_stats(x)
    obs = cast("pd.DataFrame", adata.obs)
    for i, name in enumerate(_ROW_STATS):
        obs[name] = stats[:, i].astype(np.int64 if name == "nnz" else np.float64)
    obs["n_measured_vars"] = np.full(len(obs), int(kept.sum()), dtype=np.int64)
    return adata


def harmonize(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Apply harmonization: `project_obs`, `add_tissue_general`, `select_counts`,
    `add_obs_stats`, then `project_genes`.

    Harmonization changes values and annotations so that a query can combine
    stores from different datasets. The other steps are not specified yet.
    See the module docstring for the constraints on whatever lands here.

    `harmonize` and each step consume the AnnData that they receive. A step can
    change that object and return it, or return a new object that shares its
    parts. Do not use the input after the call, and do not assume that the
    result is the same object. A caller that needs the input afterwards passes a
    copy.
    """
    adata = project_obs(adata, config, refdata)
    adata = add_tissue_general(adata, config, refdata)
    adata = select_counts(adata, config, refdata)
    adata = add_obs_stats(adata, config, refdata)
    return project_genes(adata, config, refdata)


def read_uns(f: h5py.File, keys: Iterable[str] | None = None) -> Config:
    """The `uns` mapping, restricted to `keys` unless `keys` is None.

    Reading the whole group is not safe on every dataset. A spatial dataset can
    hold a multi-GB image under `uns/spatial`, and the slot policy discards it
    a moment later. Read what the policy keeps, and nothing else.
    """
    if "uns" not in f:
        return {}
    if keys is None:
        return cast("Config", read_elem(f["uns"]))
    return {k: cast("Any", read_elem(f["uns"][k])) for k in keys if k in f["uns"]}


def slots_of(adata: anndata.AnnData) -> list[str]:
    """What a store written from this AnnData will hold.

    Read from the object, so it reports what is there rather than what a caller
    meant to put there.
    """
    slots = ["obs", "var"]
    if adata.X is not None:
        slots.append("X")
    if adata.raw is not None:
        slots.append("raw")
    slots.append(f"uns[{len(adata.uns)} keys]")
    for name in OPTIONAL_SLOTS:
        # anndata 0.13 keeps `X` in `layers[None]` in memory. It is reported
        # above, and a None key does not sort against the rest.
        keys = [k for k in getattr(adata, name) if k is not None]
        slots += [f"{name}/{key}" for key in sorted(keys)]
    return slots


def build(
    f: h5py.File,
    *,
    dataset_id: str,
    dataset_version_id: str,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> tuple[anndata.AnnData | None, Result]:
    """Read one Discover dataset and return (AnnData to store, result).

    The AnnData is `None` when the dataset is not admitted, which is a normal
    outcome. `result` always carries the identifiers, `admitted`, `reason`, and
    `slots`. `slots` describes the object after `harmonize`, so it reports what
    a store will hold. After `harmonize`, `build` adds the `dataset_id` column
    to `obs`; see `add_dataset_id`. The AnnData is thus what the store holds.

    Raises `ValueError` when the source `obs` already has a `dataset_id` column.
    That name belongs to the store, and the writer does not decide what to do
    with a source column of that name.

    `f` must stay open until the returned object is written.
    """
    config = dict(DEFAULTS) if config is None else config
    all_slots = bool(config.get("all_slots", False))
    # A config value is Any by design, so name its shape where it is read.
    uns_keys = set(cast("list[str]", config.get("uns_keys", DEFAULTS["uns_keys"])))

    # Assemble the core object first: obs, var, uns are small and X is lazy, so
    # this is cheap even for a dataset that turns out not to be admitted.
    # `read_elem` returns a wide union, so name what these are.
    obs = cast("pd.DataFrame", read_elem(f["obs"]))
    if DATASET_ID_COLUMN in obs.columns:
        raise ValueError(
            f"unexpected {DATASET_ID_COLUMN} column in the source obs: "
            "a CxG h5ad does not hold one, and in a store it is the column that the "
            "writer adds. The input is not a Discover dataset as published"
        )
    var = cast("pd.DataFrame", read_elem(f["var"]))
    n_rows = int(obs.shape[0])

    adata = anndata.AnnData(X=lazy_matrix(f["X"], n_rows), obs=obs, var=var)
    adata.uns = read_uns(f, None if all_slots else uns_keys)
    if "raw" in f and "X" in f["raw"]:
        raw_var = (
            cast("pd.DataFrame", read_elem(f["raw/var"])) if "var" in f["raw"] else var
        )
        adata.raw = anndata.AnnData(
            X=lazy_matrix(f["raw/X"], n_rows), obs=obs[[]], var=raw_var
        )

    # Admission and the cell filters run before the optional slots are read.
    reason = admit(adata, config, refdata)
    dropped_cells: dict[str, int] = {}
    if reason is None:
        keep, dropped_cells = cell_mask(obs, config)
        if not keep.any():
            reason = _no_cells_reason(obs, dropped_cells)
    if reason is not None:
        return None, {
            "dataset_id": dataset_id,
            "dataset_version_id": dataset_version_id,
            "admitted": False,
            "reason": reason,
            "slots": [],
            "dropped_cells": dropped_cells,
            "dropped_obs_columns": [],
            "dropped_var_features": [],
            "gene_reference_version": None,
        }

    if all_slots:
        for key in OPTIONAL_SLOTS:
            if key not in f:
                continue
            slot = getattr(adata, key)
            for sub in f[key]:
                slot[sub] = read_elem(f[key][sub])

    # This runs before harmonize, so `project_obs` sees only the kept cells.
    if not keep.all():
        adata = adata[keep].copy()

    # Report the object that harmonize returns, not the one assembled above, so
    # that a step which renames, adds, or drops a slot is reflected. The same
    # rule gives the dropped columns and features: diff `obs` and the counts
    # `var` across the call rather than have harmonize report them, so a drop is
    # reported whatever its reason.
    before = list(adata.obs.columns)
    counts_var = adata.raw.var if adata.raw is not None else adata.var
    features = [str(i) for i in counts_var.index]
    adata = harmonize(adata, config, refdata)
    kept = set(adata.obs.columns)
    stored = set(adata.var.index.astype(str))
    add_dataset_id(adata, dataset_id)
    # harmonize projects onto the reference that `gene_references` names for
    # the organism, which admission requires.
    spec = _gene_references(config).get(str(adata.uns.get("organism_ontology_term_id")))
    return adata, {
        "dataset_id": dataset_id,
        "dataset_version_id": dataset_version_id,
        "admitted": True,
        "reason": None,
        "slots": slots_of(adata),
        "dropped_cells": dropped_cells,
        "dropped_obs_columns": [c for c in before if c not in kept],
        "dropped_var_features": [i for i in features if i not in stored],
        "gene_reference_version": None if spec is None else spec["version"],
    }


REMOTE_SCHEMES = ("s3://", "gs://", "az://", "abfs://", "http://", "https://")


def is_remote(dest: object) -> bool:
    """Whether `dest` is a remote URL rather than a local path."""
    return isinstance(dest, (str, os.PathLike)) and str(dest).startswith(REMOTE_SCHEMES)


def as_store(dest: StoreLike) -> StoreLike:
    """Adapt a destination to something zarr can write to.

    zarr ships no S3 client, so a remote URL becomes a `zarr.storage.ObjectStore`
    over `obstore`. Local paths and Store instances pass through untouched.
    """
    if is_remote(dest):
        import obstore.store

        return zarr.storage.ObjectStore(obstore.store.from_url(str(dest)))
    return dest


def _thread_limit(threads: int | None) -> contextlib.AbstractContextManager[Any]:
    """Bound the dask threads to `threads`, or leave the default when None.

    `dask.config.set` changes the setting when it is called, not when its block
    starts, so call this in the `with` statement that it bounds.
    """
    if threads:
        return dask.config.set(num_workers=threads)
    return contextlib.nullcontext()


def write(
    adata: anndata.AnnData,
    dest: StoreLike,
    *,
    summary: SummaryRecord,
    consolidate: bool = True,
    threads: int | None = None,
) -> None:
    """Write an AnnData to a Zarr store at `dest`, then its summary record.

    anndata consolidates the metadata itself, so this passes the choice through.
    Leave it on. Consolidated metadata is not part of the Zarr v3 specification,
    so zarr warns about it. Readers such as `anndata.experimental.read_lazy`
    still need it for tolerable remote latency.

    `adata.obs` must hold the `dataset_id` column that `build` adds, with the
    summary record's `dataset_id`, and an `observation_joinid` with no null and
    no repeated value. The pair is the cell key (design 1.6), so every write
    checks it, and no configuration turns the check off. The summary record
    goes into the root group's attributes after anndata finishes, so it is the
    last write.
    The root group opens with its consolidated metadata, and the attribute
    update rewrites the root `zarr.json` with that metadata intact.

    Writing computes the lazy matrix, which anndata does with the dask threaded
    scheduler. `threads` is how many threads it may use. Give it the CPU count
    the caller was allocated. The default asks dask, and dask counts the CPUs of
    the host unless a cgroup sets a quota. A container with a CPU *share* sets
    no quota, so the default over-counts: a HealthOmics task allocated 2 vCPU
    was measured using 4.67.
    """
    # Refuse a malformed record before anything is written, so that a caller with
    # no type checker cannot leave a store with one.
    validate_summary(summary)
    obs = cast("pd.DataFrame", adata.obs)
    # unique() works on the categorical codes, not on one string per cell.
    if DATASET_ID_COLUMN not in obs.columns or {
        str(v) for v in obs[DATASET_ID_COLUMN].unique()
    } != {summary["dataset_id"]}:
        raise ValueError(
            f"obs {DATASET_ID_COLUMN} must hold only the summary record's "
            f"{summary['dataset_id']}; `build` adds that column"
        )
    if "observation_joinid" not in obs.columns:
        raise ValueError(
            "obs lacks observation_joinid, the second half of the cell key"
        )
    _require_unique(obs["observation_joinid"], "observation_joinid")
    store = as_store(dest)
    with _thread_limit(threads):
        adata.write_zarr(store, consolidate_metadata=consolidate)
    root = zarr.open_group(store, mode="r+")
    root.update_attributes({SUMMARY_NAMESPACE: {SUMMARY_KEY: dict(summary)}})


def add_dataset_id(adata: anndata.AnnData, dataset_id: str) -> None:
    """Add the `dataset_id` column to `adata.obs`, as its first column.

    `build` calls this after `harmonize`. The column is categorical with one
    category, so it costs one code per cell. `obs` must not hold a column of
    that name already.
    """
    obs = cast("pd.DataFrame", adata.obs)
    if DATASET_ID_COLUMN in obs.columns:
        raise ValueError(f"obs already has a {DATASET_ID_COLUMN} column")
    # From the codes, so that no list of one string per cell is made.
    codes = np.zeros(len(obs), dtype=np.int8)
    categorical = pd.Categorical.from_codes(codes, categories=pd.Index([dataset_id]))
    column = pd.Series(categorical, index=obs.index)
    obs.insert(0, DATASET_ID_COLUMN, column)


def summary_record(result: Result) -> SummaryRecord:
    """The summary record of an admitted dataset, from its result."""
    if not result["admitted"] or result["gene_reference_version"] is None:
        raise ValueError("a dataset that is not admitted has no store and no record")
    return {
        "version": SUMMARY_VERSION,
        "dataset_id": result["dataset_id"],
        "dataset_version_id": result["dataset_version_id"],
        "gene_reference_version": result["gene_reference_version"],
    }


def read_summary(dest: StoreLike | zarr.Group) -> SummaryRecord | None:
    """The summary record of a store, or None when the store has none.

    `dest` is a store's location, or its root group already open. A record that is
    present but malformed raises `SummaryRecordError`; see
    `validate_summary`.
    """
    root = (
        dest
        if isinstance(dest, zarr.Group)
        else zarr.open_group(as_store(dest), mode="r")
    )
    namespace = root.attrs.get(SUMMARY_NAMESPACE)
    if namespace is None:
        return None
    if not isinstance(namespace, dict):
        raise SummaryRecordError(f"{SUMMARY_NAMESPACE} attribute is not a JSON object")
    if SUMMARY_KEY not in namespace:
        return None
    return validate_summary(namespace[SUMMARY_KEY])


@contextlib.contextmanager
def _opened(source: Source) -> Iterator[h5py.File]:
    """Yield an open h5py.File, whether given a path or an already-open file."""
    if isinstance(source, (str, os.PathLike)):
        with h5py.File(source, "r") as f:
            yield f
    else:
        yield source


def create(
    source: Source,
    dest: StoreLike,
    *,
    dataset_id: str,
    dataset_version_id: str,
    consolidate: bool = True,
    config: Config | None = None,
    refdata: PathLike | None = None,
    threads: int | None = None,
) -> Result:
    """Build a store from `source` and write it to `dest`. Return the result.

    Writes nothing and reports `admitted: False` with a reason when the dataset
    is not admitted. That is a normal outcome, not an error.

    `source` is a path to an h5ad or an open `h5py.File`. `dest` is a path, a
    remote URL, or a zarr Store. `dataset_id` and `dataset_version_id` go into
    the store's summary record. The source stays open across the write, because
    the expression matrix is read lazily. `threads` bounds the write, and the
    read of `X` in `add_obs_stats`; see `write`.
    """
    with _opened(source) as f:
        with _thread_limit(threads):
            adata, result = build(
                f,
                dataset_id=dataset_id,
                dataset_version_id=dataset_version_id,
                config=config,
                refdata=refdata,
            )
        if adata is None:
            return result
        summary = summary_record(result)
        write(adata, dest, summary=summary, consolidate=consolidate, threads=threads)
    return result
