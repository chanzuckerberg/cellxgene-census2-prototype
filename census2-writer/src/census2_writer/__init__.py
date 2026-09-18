"""Write one Census 2 store from one Discover dataset.

The census writer reads one published Discover `.h5ad` and writes one
AnnData-Zarr store. The bulk ingestion pipeline (`bulk-ingest-nf`) calls
this module once for each dataset. A developer can call it from Python.
The module reads one dataset. It does not read the corpus, so a step that
needs a corpus-wide aggregate cannot run here.

Public calls
------------
`create(source, dest, *, config=None, refdata=None, consolidate=True)`
builds the store and writes it. It writes nothing when the writer does not
admit the dataset. It returns a `Result` with `admitted`, `reason`, and
`slots`.

    import census2_writer

    result = census2_writer.create("dataset.h5ad", "store.zarr")
    if not result["admitted"]:
        print("skipped:", result["reason"])

`build(f, *, config=None, refdata=None)` returns the AnnData and the
`Result`. It writes nothing. `write(adata, dest, *, consolidate=True)` writes
what `build` returned. Keep `f` open until `write` returns, because `X` is a
dask graph until then (see Memory).

    with h5py.File("dataset.h5ad") as f:
        adata, result = census2_writer.build(f)
        if adata is not None:
            census2_writer.write(adata, "store.zarr")

`load_config(path)` returns the configuration mapping. `provenance(config)`
returns the part of the receipt that this module knows.

`dest` is a local path, a remote URL such as `s3://bucket/prefix.zarr`, or a
`zarr.abc.store.Store`. zarr has no S3 client, so `as_store` adapts a remote
URL to `zarr.storage.ObjectStore` over `obstore`. Credentials come from the
standard AWS chain.

Decisions
---------
`build` runs four steps in this order: assemble the core slots (`obs`,
`var`, `X`, `raw`, full `uns`), admit, apply the slot policy, harmonize.

1. Admission decides whether the dataset belongs in Census 2. A rejection
   is a normal result. `build` returns no AnnData and a reason, and the
   caller records the reason and exits with success. This decision is open.
   `admit` is a stub that admits every dataset.
2. The slot policy decides which AnnData slots the store carries. This
   decision is made. The default set is the set that the Census 1 builder
   reads: `obs`, `var`, `X`, `raw`, and the two `uns` keys in
   `CENSUS_UNS_KEYS`. The `all_slots` configuration key carries every slot.
3. Harmonization changes values and annotations so that a query can combine
   stores. This decision is open. `harmonize` is a stub that returns its
   input.

Admission runs before the slot policy so that a criterion can read any
`uns` key. It runs before the optional slots are read so that a rejected
dataset does not read `obsm` or `obsp`. `admit` and `harmonize` both
receive an AnnData, so a step author works with one kind of object.

Configuration and reference data
--------------------------------
`config` is a mapping that the caller loads and passes in. It holds policy:
thresholds, criteria, and the names of lists. It is small and lives in
version control next to a release. `refdata` is the path of a directory
that the caller staged before the run. It holds bulk reference data:
ontology term lists, gene references, tissue roll-ups. The configuration
names a reference file and does not give its path. A step resolves each
name against `refdata`, because the staged location is known only inside
the task. This module does not fetch a file at run time. A HealthOmics task
has no public internet access.

Today `load_config` ignores `path` and returns `DEFAULTS`. No file format is
defined. `all_slots` is the only key that `build` and `provenance` read.
`admit` and `harmonize` receive `config` and `refdata` so that a step can
use them. Do not define a configuration schema before the step that needs
it exists.

Receipt
-------
A receipt is the JSON file that a pipeline writes next to a store. Today it
is a private log. `bulk-ingest-nf` writes one for each dataset, and its
report step reads them to summarize the run. Nothing outside a run reads
one. `Receipt` declares its keys, so the script that writes a receipt and
the script that reads it agree on the shape; the two cannot import each
other. `provenance` supplies the recipe (`RECIPE`), the slot policy, the
Zarr format, and the versions of this package, anndata, and zarr. It reads
the versions from the installed distributions. The pipeline adds the
dataset source and the run facts.

Registering a store, with metadata that a reader outside the run can rely
on, is a separate job and is not implemented yet.

Memory
------
`build` returns a complete AnnData but reads only `obs`, `var`, and `uns`
into memory. `X` and `raw.X` are dask arrays of row blocks of about
`BLOCK_BYTES`: sparse blocks for a sparse dataset, dense blocks for a dense
dataset. Memory use follows the block size and not the file size. The
corpus p99 is about 26 GiB and the maximum is 50 GiB.

Two rules follow:

- An `h5py.File` given to `build` must stay open until `write` returns,
  because `X` is an unevaluated dask graph. `create` keeps the file open
  when it receives a path.
- A step that needs the complete matrix at one time forces a compute. So
  does a mask computed from `X`. That compute can fail on a large dataset.
  Row-wise work maps over blocks.

A subset of the AnnData keeps `obs`, `var`, `X`, and the other slots
aligned.

The README in this directory covers installation, use, and the development
workflow.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import os
from collections.abc import Iterator
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, TypedDict, cast

import anndata
import h5py
import numpy as np
import zarr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem
from zarr.abc.store import Store

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


# uns keys the CELLxGENE Census builder reads. See AnnDataProxy in the
# chanzuckerberg/cellxgene-census repository, at
# tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/anndata.py
CENSUS_UNS_KEYS = ("schema_version", "organism_ontology_term_id")
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
    census2_writer_version: str
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
        "census2_writer_version": importlib.metadata.version("census2-writer"),
        "anndata_version": importlib.metadata.version("anndata"),
        "zarr_version": importlib.metadata.version("zarr"),
    }


def admit(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> str | None:
    """Decide whether Census 2 admits this dataset.

    Return `None` to admit the dataset. Return a short reason to reject it.
    A rejection is a normal result. The caller records the reason and exits
    with success. Do not raise an exception for a dataset that you reject.
    Raise only for an input that this function cannot read.

    `adata` holds `obs`, `var`, `X`, `raw`, and the full `uns`. The slot
    policy runs after admission. All `uns` keys are present, and a criterion
    can read any of them. `X` and `raw.X` are not loaded. A read of `X` or
    `raw.X` loads the full matrix.

    `refdata` is the directory that the caller staged. The configuration
    names a reference file. It does not give the path. Resolve each name
    against `refdata`.

    This function is a stub that admits every dataset. The admission
    criteria are open.
    """
    return None


def harmonize(
    adata: anndata.AnnData,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> anndata.AnnData:
    """Apply harmonization. STUB: returns its input unchanged.

    Nothing is specified yet. See the module docstring for the constraints on
    whatever lands here.
    """
    return adata


def build(
    f: h5py.File,
    *,
    config: Config | None = None,
    refdata: PathLike | None = None,
) -> tuple[anndata.AnnData | None, Result]:
    """Read one Discover dataset and return (AnnData to store, result).

    The AnnData is `None` when the dataset is not admitted, which is a normal
    outcome. `result` always carries `admitted`, `reason`, and `slots`.

    `f` must stay open until the returned object is written.
    """
    config = dict(DEFAULTS) if config is None else config
    all_slots = bool(config.get("all_slots", False))

    # Assemble the core object first: obs, var, uns are small and X is lazy, so
    # this is cheap even for a dataset that turns out not to be admitted.
    # `read_elem` returns a wide union, so name what these are.
    obs = cast("pd.DataFrame", read_elem(f["obs"]))
    var = cast("pd.DataFrame", read_elem(f["var"]))
    n_rows = int(obs.shape[0])

    adata = anndata.AnnData(X=lazy_matrix(f["X"], n_rows), obs=obs, var=var)
    adata.uns = cast("Config", read_elem(f["uns"])) if "uns" in f else {}
    if "raw" in f and "X" in f["raw"]:
        raw_var = (
            cast("pd.DataFrame", read_elem(f["raw/var"])) if "var" in f["raw"] else var
        )
        adata.raw = anndata.AnnData(
            X=lazy_matrix(f["raw/X"], n_rows), obs=obs[[]], var=raw_var
        )

    # Admission runs before slot policy, and before the optional slots are read.
    reason = admit(adata, config, refdata)
    if reason is not None:
        return None, {"admitted": False, "reason": reason, "slots": []}

    carried: list[str] = ["obs", "var", "X"]
    if adata.raw is not None:
        carried.append("raw")

    if not all_slots:
        adata.uns = {k: v for k, v in adata.uns.items() if k in CENSUS_UNS_KEYS}
    carried.append(f"uns[{len(adata.uns)} keys]")

    if all_slots:
        for key in OPTIONAL_SLOTS:
            if key not in f:
                continue
            slot = getattr(adata, key)
            for sub in f[key]:
                slot[sub] = read_elem(f[key][sub])
                carried.append(f"{key}/{sub}")

    return harmonize(adata, config, refdata), {
        "admitted": True,
        "reason": None,
        "slots": carried,
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


def write(adata: anndata.AnnData, dest: StoreLike, *, consolidate: bool = True) -> None:
    """Write an AnnData to a Zarr store at `dest`.

    anndata consolidates the metadata itself, so this passes the choice through.
    Leave it on. Consolidated metadata is not part of the Zarr v3 specification,
    so zarr warns about it. Readers such as `anndata.experimental.read_lazy`
    still need it for tolerable remote latency. anndata also writes it last, so a
    reader that opens with `zarr.open_consolidated` fails on a partially written
    store rather than reading part of one.
    """
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
) -> Result:
    """Build a store from `source` and write it to `dest`. Return the result.

    Writes nothing and reports `admitted: False` with a reason when the dataset
    is not admitted. That is a normal outcome, not an error.

    `source` is a path to an h5ad or an open `h5py.File`. `dest` is a path, a
    remote URL, or a zarr Store. The source stays open across the write, because
    the expression matrix is read lazily.
    """
    with _opened(source) as f:
        adata, result = build(f, config=config, refdata=refdata)
        if adata is None:
            return result
        write(adata, dest, consolidate=consolidate)
    return result
