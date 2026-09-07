"""Build one Census 2 store from one Discover dataset.

This module owns the data model, meaning what a store contains, and it can write
the result. It is usable standalone:

    import census2_store

    # both steps; nothing is written if the dataset is not admitted
    result = census2_store.create("dataset.h5ad", "store.zarr")

    # policy comes from configuration
    cfg = census2_store.load_config("ingest.json") | {"all_slots": True}
    census2_store.create("dataset.h5ad", "store.zarr", config=cfg)
    if not result["admitted"]:
        print("skipped:", result["reason"])

    # or the processing alone, to inspect or alter before writing
    with h5py.File("dataset.h5ad") as f:
        adata, result = census2_store.build(f)
        if adata is not None:
            census2_store.write(adata, "store.zarr")

`dest` is a local path, a remote URL such as `s3://bucket/prefix.zarr`, or a
`zarr.abc.store.Store` instance. zarr has no S3 client of its own, so a remote
URL is adapted to `zarr.storage.ObjectStore` over `obstore`. Credentials come
from the standard AWS chain.

Three kinds of decision live here.

1. **Admission.** Whether the dataset belongs in Census 2 at all. Not every
   Discover dataset is admitted. A rejection is a **normal outcome**, not an
   error: `build` returns no AnnData and a reason, and the caller records that
   and succeeds. Not decided. `admit` is a stub that admits everything.
2. **Slot selection.** Which AnnData slots a store carries. Decided: by default
   the set the Census 1 builder reads, which is `obs`, `var`, `X`, `raw`, and two
   `uns` keys. Set `all_slots` in the configuration to carry everything instead.
3. **Harmonization.** How values and annotations are changed so that stores are
   queryable together. Not decided. `harmonize` is a stub that returns its input.

Census 2 writes one store per dataset and composes stores at query time, so
cross-dataset queries need harmonized data and harmonization has to happen while
writing a single store. A store is therefore not the Discover h5ad in another
container. The constraint on whatever lands in `harmonize`: ingestion reads
exactly one dataset, so a step needing a corpus-wide aggregate has to be
redesigned or eliminated.

Configuration and reference data
--------------------------------
STUB. `load_config` returns an empty mapping and nothing reads it yet.

Admission and harmonization will need both, and they are different things:

- **Configuration** is policy: thresholds, which criteria apply, which lists to
  use. Small, and it belongs in version control next to a release.
- **Reference data** is bulk: ontology term lists, gene references, tissue
  roll-ups. Too big for config, and it changes on its own schedule.

The constraint that shapes how they arrive: a HealthOmics task has no public
internet access, so **nothing here may fetch a reference file at run time.** It
has to be staged into the task as a workflow input. So `config` is a mapping the
caller loads and passes in, and `refdata` is the path to a directory the caller
staged. Config should name reference files rather than give paths to them,
because the staged location is only known inside the task; resolve names against
`refdata`.

Both are threaded through `admit` and `harmonize` now so the plumbing exists.
`DEFAULTS` holds the one policy item that is already decided, `all_slots`, which
makes configuration the single place policy lives rather than a second mechanism
beside it. Everything else is undecided, so do not invent a schema here before
the steps that need it.

Both `admit` and `harmonize` receive an AnnData, so a step author works with one
kind of object rather than two. `build` runs them in this order:

    assemble core (obs, var, X, raw, full uns)  ->  admit  ->  slot policy  ->  harmonize

Admission comes before slot policy so a criterion can read any `uns` key, and
before the optional slots are read so a rejected dataset never pays for `obsm`
or `obsp`.

Memory contract
---------------
`build` returns a whole AnnData, but not a whole read. `obs`, `var`, and
`uns` are in memory, because they are small. `X` and `raw.X` are dask arrays of
sparse row blocks, so memory tracks the block size and not the file size. That
matters at the corpus p99 of about 26 GiB and its 50 GiB maximum.

Two consequences the caller and any future step must respect:

- A `h5py.File` handed to `build` **must stay open until the store is written**,
  because `X` is still an unevaluated dask graph. `create` handles that for you
  when given a path; it is only a trap when composing `build` and `write`
  yourself.
- A step that needs the whole matrix at once, or a mask computed from `X`, forces
  a compute that a large dataset may not survive. Row-wise work maps over blocks
  and stays cheap.

Subsetting the AnnData keeps `obs`, `var`, `X`, and the other slots aligned, so
consistency is the object's job rather than a step author's.
"""

from __future__ import annotations

import contextlib
import os

import anndata
import h5py
import numpy as np
import zarr
from anndata.experimental import read_elem_lazy
from anndata.io import read_elem

# Slots the CELLxGENE Census builder reads. See AnnDataProxy in the
# chanzuckerberg/cellxgene-census repository, at
# tools/cellxgene_census_builder/src/cellxgene_census_builder/build_soma/anndata.py
CENSUS_UNS_KEYS = ("schema_version", "organism_ontology_term_id")
OPTIONAL_SLOTS = ("obsm", "varm", "obsp", "varp", "layers")

# Target bytes per row block of the expression matrix.
BLOCK_BYTES = 64 * 1024**2


def row_chunk(h5_obj, n_rows: int) -> int:
    """Rows per dask block, sized so one block is about BLOCK_BYTES."""
    if isinstance(h5_obj, h5py.Group):  # sparse: data + indices per stored value
        nnz = h5_obj["data"].shape[0]
        per_value = h5_obj["data"].dtype.itemsize + h5_obj["indices"].dtype.itemsize
        row_bytes = max(1, (nnz // max(n_rows, 1)) * per_value)
    else:  # dense
        row_bytes = max(1, int(np.prod(h5_obj.shape[1:])) * h5_obj.dtype.itemsize)
    return max(1, min(n_rows, BLOCK_BYTES // row_bytes))


def lazy_matrix(h5_obj, n_rows: int):
    """A dask array of sparse row blocks."""
    return read_elem_lazy(h5_obj, chunks=(row_chunk(h5_obj, n_rows), -1))


# Policy that is already decided. Anything a caller can choose about what a store
# contains belongs here, not in a separate function argument.
DEFAULTS: dict = {
    # Carry every AnnData slot instead of the Census set. See build().
    "all_slots": False,
}


def load_config(path=None) -> dict:
    """Read ingestion configuration from `path`, over `DEFAULTS`.

    STUB: ignores `path` and returns the defaults. No file format is defined yet.

    Keep this the only place that parses configuration, so the format stays one
    decision. See the module docstring for how config and reference data differ
    and why reference data cannot be fetched at run time.
    """
    return dict(DEFAULTS)


def admit(
    adata: anndata.AnnData, config: dict | None = None, refdata=None
) -> str | None:
    """Decide whether this dataset belongs in Census 2.

    Return `None` to admit it, or a short human-readable reason to reject it.
    A rejection is a normal outcome that the caller records and succeeds on, so
    do not raise here for an unsuitable dataset. Raise only for a genuinely
    broken input.

    `adata` is the core object with its **full** `uns`, before slot policy trims
    it, so a criterion may read any `uns` key. `X` and `raw.X` are still lazy, so
    reading either forces a compute.

    STUB: admits everything. Which datasets Census 2 admits is not decided yet.
    """
    return None


def harmonize(
    adata: anndata.AnnData, config: dict | None = None, refdata=None
) -> anndata.AnnData:
    """Apply harmonization. STUB: returns its input unchanged.

    Nothing is specified yet. See the module docstring for the constraints on
    whatever lands here.
    """
    return adata


def build(
    f: h5py.File, *, config: dict | None = None, refdata=None
) -> tuple[anndata.AnnData | None, dict]:
    """Read one Discover dataset and return (AnnData to store, result).

    The AnnData is `None` when the dataset is not admitted, which is a normal
    outcome. `result` always carries `admitted`, `reason`, and `slots`.

    `f` must stay open until the returned object is written.
    """
    config = dict(DEFAULTS) if config is None else config
    all_slots = bool(config.get("all_slots", False))

    # Assemble the core object first: obs, var, uns are small and X is lazy, so
    # this is cheap even for a dataset that turns out not to be admitted.
    obs = read_elem(f["obs"])
    var = read_elem(f["var"])
    n_rows = obs.shape[0]

    adata = anndata.AnnData(X=lazy_matrix(f["X"], n_rows), obs=obs, var=var)
    adata.uns = read_elem(f["uns"]) if "uns" in f else {}
    if "raw" in f and "X" in f["raw"]:
        raw_var = read_elem(f["raw/var"]) if "var" in f["raw"] else var
        adata.raw = anndata.AnnData(
            X=lazy_matrix(f["raw/X"], n_rows), obs=obs[[]], var=raw_var
        )

    # Admission runs before slot policy, so a criterion can read any uns key.
    # It also runs before the optional slots are read, so a rejected dataset
    # never pays for obsm or obsp.
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


def is_remote(dest) -> bool:
    return isinstance(dest, (str, os.PathLike)) and str(dest).startswith(REMOTE_SCHEMES)


def as_store(dest):
    """Adapt a destination to something zarr can write to.

    zarr ships no S3 client, so a remote URL becomes a `zarr.storage.ObjectStore`
    over `obstore`. Local paths and Store instances pass through untouched.
    """
    if is_remote(dest):
        import obstore.store

        return zarr.storage.ObjectStore(obstore.store.from_url(str(dest)))
    return dest


def write(adata: anndata.AnnData, dest, *, consolidate: bool = True) -> None:
    """Write an AnnData to a Zarr store at `dest`.

    anndata consolidates the metadata itself, so this passes the choice through
    rather than consolidating again. Leave it on: consolidated metadata is not
    part of the Zarr v3 specification, so zarr warns about it, but readers such
    as `anndata.experimental.read_lazy` need it for tolerable remote latency.
    It is also written last, so a reader that opens with `zarr.open_consolidated`
    fails on a partially written store rather than reading part of one.
    """
    adata.write_zarr(as_store(dest), consolidate_metadata=consolidate)


@contextlib.contextmanager
def _opened(source):
    """Yield an open h5py.File, whether given a path or an already-open file."""
    if isinstance(source, (str, os.PathLike)):
        with h5py.File(source, "r") as f:
            yield f
    else:
        yield source


def create(
    source,
    dest,
    *,
    consolidate: bool = True,
    config: dict | None = None,
    refdata=None,
) -> dict:
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
