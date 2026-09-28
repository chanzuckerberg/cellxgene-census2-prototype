"""Write one Census 2 store from one Discover dataset.

The writer reads one published Discover `.h5ad` and writes one AnnData-Zarr
store. The bulk ingestion pipeline (`bulk-ingest-nf`) calls this module once
for each dataset, and a developer can call it from Python. The module reads
one dataset and not the corpus, so a step that needs a corpus-wide
aggregate cannot run here.

Public calls
------------
`create(source, dest, ...)` builds the store and writes it. It writes
nothing when the dataset is not admitted, and returns a `Result` in both
cases. `build` and `write` split the same work into two calls, so a caller
can inspect the AnnData before the write:

    with h5py.File("dataset.h5ad") as f:
        adata, result = census2_dataset_writer.build(f)
        if adata is not None:
            census2_dataset_writer.write(adata, "store.zarr")

Keep `f` open until `write` returns, because `X` is a dask graph until then
(see Memory). `create` does this itself when it receives a path.

`slots_of(adata)` lists what a store written from an AnnData will hold.
`load_config(path)` returns the configuration mapping. `provenance(config)`
returns the part of a receipt that this module knows. `as_store(dest)`
adapts a remote URL for zarr, which ships no S3 client.

Decisions
---------
`build` runs four steps in order: assemble the core slots (`obs`, `var`,
`X`, `raw`, the carried `uns` keys), admit, apply the slot policy, harmonize.

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
   decides which `obs` columns a store carries, and removes only the columns
   that Zarr cannot name. `select_counts` puts the counts in `X`: `raw.X`
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

Receipt
-------
A receipt is the JSON a pipeline writes beside a store. `Receipt` declares
its shape, and `provenance` fills the part this package knows; see their
docstrings.

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
import zarr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem
from zarr.abc.store import Store

from . import assays, genes

if TYPE_CHECKING:  # these are only ever types here, never constructed
    import dask.array as da
    import pandas as pd

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
    """The JSON a pipeline writes beside a store.

    A private log today, read only by the pipeline's own report step. The shape
    lives here because the script that writes one and the script that reads one
    cannot import each other, and because `provenance` fills the part that only
    this package knows.

    The `StoreProvenance` keys are always present. The rest is what a writer
    knows about where the dataset came from and what it did, so it is optional
    here. `bulk-ingest-nf` fills all of them for an admitted dataset. A rejected
    dataset gets no `slots_written`, `store`, `output_bytes`, or `published_by`,
    and `check` appears only when the task runs with `--check`.
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
    gene_reference_version: str | None
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

    references = cast(
        "dict[str, Any]",
        config.get("gene_references", DEFAULTS["gene_references"]),
    )
    organism = adata.uns.get("organism_ontology_term_id")
    if organism is not None and str(organism) not in references:
        return f"no gene reference for organism {organism}"
    return None


def project_obs(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Project `obs` onto the columns a store carries.

    THIS IS A MINIMAL STUB. It removes the `obs` columns whose names hold a
    `/`, and keeps every other column. It is not production-ready projection
    logic. Zarr gives `/` its path meaning, so a column named that way fails
    the write. Removing it is what lets the dataset be stored at all.

    Census 1 projects `obs` onto a closed schema: about 20 named CxG fields,
    plus columns it computes. It carries no author-chosen column, because its
    `obs` is one dataframe joined across every dataset and so must have one
    column set. Census 2 writes one store per dataset and has no such
    constraint, so what a store carries is open. Until that is decided, this
    keeps every column that Zarr can hold.

    `var` is not projected. The names seen with a `/` were all in `obs`.
    """
    # `AnnData.obs` is typed `DataFrame | Dataset2D`, the second being what
    # `read_lazy` produces. `build` always assembles obs with `read_elem`, so
    # narrow to the pandas frame here rather than widen the signature.
    obs = cast("pd.DataFrame", adata.obs)
    dropped = [c for c in obs.columns if "/" in str(c)]
    if dropped:
        adata.obs = obs.drop(columns=dropped)
    return adata


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
    return `adata` unchanged. `X` stays lazy.
    """
    if adata.raw is None:
        return adata
    obsm: dict[str, Any] = {key: adata.obsm[key] for key in adata.obsm}
    obsp: dict[str, Any] = {key: adata.obsp[key] for key in adata.obsp}
    # `AnnData.obs` is typed `DataFrame | Dataset2D`, the second being what
    # `read_lazy` produces. `build` assembles obs with `read_elem`.
    return anndata.AnnData(
        X=adata.raw.X,
        obs=cast("pd.DataFrame", adata.obs),
        var=adata.raw.var,
        uns=dict(adata.uns),
        obsm=obsm,
        obsp=obsp,
    )


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
    The file must be staged in `refdata`.
    """
    config = DEFAULTS if config is None else config
    references = cast(
        "dict[str, genes.Reference]",
        config.get("gene_references", DEFAULTS["gene_references"]),
    )
    return genes.project(adata, references, refdata)


def harmonize(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Apply harmonization: `project_obs`, `select_counts`, then `project_genes`.

    Harmonization changes values and annotations so that a query can combine
    stores from different datasets. The other steps are not specified yet.
    See the module docstring for the constraints on whatever lands here.
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
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> tuple[anndata.AnnData | None, Result]:
    """Read one Discover dataset and return (AnnData to store, result).

    The AnnData is `None` when the dataset is not admitted, which is a normal
    outcome. `result` always carries `admitted`, `reason`, and `slots`. `slots`
    describes the object after `harmonize`, so it reports what a store will
    hold.

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
    # harmonize projects onto the reference that `gene_references` names for
    # the organism, which admission requires.
    references = cast(
        "dict[str, genes.Reference]",
        config.get("gene_references", DEFAULTS["gene_references"]),
    )
    spec = references.get(str(adata.uns.get("organism_ontology_term_id")))
    return adata, {
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
    consolidate: bool = True,
    threads: int | None = None,
) -> None:
    """Write an AnnData to a Zarr store at `dest`.

    anndata consolidates the metadata itself, so this passes the choice through.
    Leave it on. Consolidated metadata is not part of the Zarr v3 specification,
    so zarr warns about it. Readers such as `anndata.experimental.read_lazy`
    still need it for tolerable remote latency. anndata also writes it last, so a
    reader that opens with `zarr.open_consolidated` fails on a partially written
    store rather than reading part of one.

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
    with limit:
        adata.write_zarr(as_store(dest), consolidate_metadata=consolidate)


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
    consolidate: bool = True,
    config: Config | None = None,
    refdata: PathLike | None = None,
    threads: int | None = None,
) -> Result:
    """Build a store from `source` and write it to `dest`. Return the result.

    Writes nothing and reports `admitted: False` with a reason when the dataset
    is not admitted. That is a normal outcome, not an error.

    `source` is a path to an h5ad or an open `h5py.File`. `dest` is a path, a
    remote URL, or a zarr Store. The source stays open across the write, because
    the expression matrix is read lazily. `threads` bounds the write; see
    `write`.
    """
    with _opened(source) as f:
        adata, result = build(f, config=config, refdata=refdata)
        if adata is None:
            return result
        write(adata, dest, consolidate=consolidate, threads=threads)
    return result
