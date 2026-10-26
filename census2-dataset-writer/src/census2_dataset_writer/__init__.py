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
`X`, `raw`, the carried `uns` keys), admit, apply the slot policy, harmonize.
It then adds the `dataset_id` column to `obs`. A source `obs` that already has
a `dataset_id` column is an error, because that name belongs to the store.

1. Admission (`admit`) decides whether the dataset belongs in Census 2. A
   rejection is a normal result: `build` returns no AnnData and a reason,
   and the caller records the reason and succeeds. This decision is largely
   open. The stub applies two criteria. The assay allowlist keeps the
   workflow away from modalities a Census is unlikely to carry. The
   organism must have a gene reference, so that harmonization can project
   the genes.
2. The slot policy decides which AnnData slots the store carries. This
   decision is made: the set that the Census 1 builder reads (`obs`, `var`,
   `X`, `raw`, and the `uns` keys that the `uns_keys` configuration lists),
   or every slot when `all_slots` is set. Harmonization then puts the
   counts from `raw` in `X`, and removes the slots that index the source
   genes.
3. Harmonization (`harmonize`) changes values and annotations so that a
   query can combine stores. It runs three steps in order. `project_obs`
   projects `obs` onto `OBS_COLUMNS`, the CxG columns that Census
   1 copies, each with a fixed kind; it drops every other source column. `select_counts` puts the counts in `X`: `raw.X`
   when the dataset has `raw`, else `X`. The store has no `raw`.
   `project_genes` projects the counts onto the pinned gene reference of the
   dataset's organism. Thus every store of one organism has the same `var`,
   in the same order. See `census2_dataset_writer.genes`.

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

`README.md` in this directory covers installation, use, and development.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import os
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, TypedDict, cast

import anndata
import dask
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem
from zarr.abc.store import Store

from . import assays, genes
from . import store_format as store_format  # a public module; callers name it

# The store format, re-exported, so a caller writes census2_dataset_writer.X.
from .store_format import BOOL as BOOL
from .store_format import CATEGORY as CATEGORY
from .store_format import DATASET_ID_COLUMN as DATASET_ID_COLUMN
from .store_format import OBS_COLUMNS as OBS_COLUMNS
from .store_format import STRING as STRING
from .store_format import SUMMARY_KEY as SUMMARY_KEY
from .store_format import SUMMARY_NAMESPACE as SUMMARY_NAMESPACE
from .store_format import SUMMARY_VERSION as SUMMARY_VERSION
from .store_format import SummaryRecord as SummaryRecord
from .store_format import SummaryRecordError as SummaryRecordError
from .store_format import validate_summary as validate_summary

if TYPE_CHECKING:  # these are only ever types here, never constructed
    import dask.array as da

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
# on it: a store's location is its dataset version alone.
RECIPE = "v1"

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
    # The assay ontology terms that admission accepts. A dataset needs at least
    # one cell on one of them. An empty list turns the criterion off. A sorted
    # list rather than the frozenset, because a config file is JSON.
    # See census2_dataset_writer.assays for where the terms come from.
    "assay_allowlist": sorted(assays.DEFAULT_ALLOWLIST),
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
    }


def admit(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> str | None:
    """Decide whether Census 2 admits this dataset.

    THIS IS A MINIMAL STUB. It applies two criteria. The dataset must hold at
    least one cell on the assay allowlist, and its organism must have an entry
    in `gene_references`. It is not production-ready admission logic. The
    assay criterion keeps the bulk ingestion workflow away from data shapes
    that a Census is unlikely to carry, such as the dense matrices of a DNA
    methylation dataset, so that the workflow is not contorted to hold them.
    The organism criterion admits only the datasets that `project_genes` can
    project. Real admission criteria are open.

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

    Two limits of the assay criterion, both deliberate:

    - One cell on the allowlist admits every cell. Census 1 masks the cells
      that do not match. This does not, so a mixed-assay dataset carries its
      off-list cells as well.
    - An `obs` without `assay_ontology_term_id` is admitted. Such an object is
      not a Discover dataset, so a Discover criterion does not apply to it.

    A `uns` without `organism_ontology_term_id` is admitted for the same
    reason. `project_genes` cannot project it, and raises an error.
    """
    config = DEFAULTS if config is None else config
    # A config value is Any by design, so name its shape where it is read.
    allowlist = set(cast("list[str]", config.get("assay_allowlist") or []))
    terms = adata.obs.get("assay_ontology_term_id")
    if allowlist and terms is not None and not bool(terms.isin(allowlist).any()):
        present = sorted({str(t) for t in terms.unique()})
        return f"no cells on the assay allowlist: {', '.join(present[:3])}"

    references = _gene_references(config)
    organism = adata.uns.get("organism_ontology_term_id")
    if organism is not None and str(organism) not in references:
        return f"no gene reference for organism {organism}"
    return None


def project_obs(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Project `obs` onto the columns of `OBS_COLUMNS`.

    Census 1 projects `obs` onto a closed schema, with fixed types: the columns
    that the CxG schema requires in every dataset, plus columns that it adds.
    This step copies the CxG columns, every column of `OBS_COLUMNS` except
    `dataset_id`, in schema order, and gives each one its kind: a categorical with the sorted
    observed values as string categories, NumPy `bool` for `is_primary_data`,
    and plain strings for `observation_joinid`. The source index stays. Every
    other source column goes: the result reports them in `dropped_obs_columns`.

    It raises `ValueError` when the source lacks a schema column, when a column
    holds a null, when `is_primary_data` is not boolean, or when
    `observation_joinid` is not unique. The CxG schema requires each of these,
    so the writer does not choose a repair.

    It does not yet add Census 1's `tissue_general` pair or the statistics of
    the counts. `build` adds `dataset_id` after harmonization.

    This step consumes `adata`, as `harmonize` states.
    """
    # `AnnData.obs` is typed `DataFrame | Dataset2D`, the second being what
    # `read_lazy` produces. `build` always assembles obs with `read_elem`, so
    # narrow to the pandas frame here rather than widen the signature.
    obs = cast("pd.DataFrame", adata.obs)
    columns = [(name, kind) for name, kind in OBS_COLUMNS if name != DATASET_ID_COLUMN]
    missing = [name for name, _ in columns if name not in obs.columns]
    if missing:
        raise ValueError(f"the source obs lacks the schema columns {missing}")
    adata.obs = pd.DataFrame(
        {name: _project_obs_column(obs[name], name, kind) for name, kind in columns},
        index=obs.index,
    )
    return adata


def _project_obs_column(column: pd.Series, name: str, kind: str) -> pd.Series:
    """One `obs` column with the kind that `OBS_COLUMNS` gives it."""
    if bool(column.isna().any()):
        raise ValueError(f"obs {name} holds a null")
    if kind == CATEGORY:
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
    if kind == BOOL:
        if column.dtype == bool:
            return column
        if str(column.dtype) == "boolean":
            return column.astype(bool)
        raise ValueError(f"obs {name} is {column.dtype}, not boolean")
    if kind == STRING:
        values = pd.Series(
            column.astype(str).to_numpy(dtype=object), index=column.index
        )
        if not values.is_unique:
            raise ValueError(f"obs {name} is not unique")
        return values
    raise ValueError(f"unknown obs column kind {kind!r}")


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


def harmonize(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Apply harmonization: `project_obs`, `select_counts`, then `project_genes`.

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
    adata = select_counts(adata, config, refdata)
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

    # Admission runs before the optional slots are read.
    reason = admit(adata, config, refdata)
    if reason is not None:
        return None, {
            "dataset_id": dataset_id,
            "dataset_version_id": dataset_version_id,
            "admitted": False,
            "reason": reason,
            "slots": [],
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
    summary record's `dataset_id`. The summary record goes into the root group's
    attributes after anndata finishes, so it is the last write.
    The root group opens with its consolidated metadata, and the attribute
    update rewrites the root `zarr.json` with that metadata intact.

    Writing computes the lazy matrix, which anndata does with the dask threaded
    scheduler. `threads` is how many threads it may use. Give it the CPU count
    the caller was allocated. The default asks dask, and dask counts the CPUs of
    the host unless a cgroup sets a quota. A container with a CPU *share* sets
    no quota, so the default over-counts: a HealthOmics task allocated 2 vCPU
    was measured using 4.67.
    """
    limit = (
        dask.config.set(num_workers=threads) if threads else contextlib.nullcontext()
    )
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
    store = as_store(dest)
    with limit:
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
    the expression matrix is read lazily. `threads` bounds the write; see
    `write`.
    """
    with _opened(source) as f:
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
