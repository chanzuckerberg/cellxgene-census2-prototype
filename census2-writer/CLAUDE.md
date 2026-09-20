# census2-writer

The data model: what a Census 2 store contains. `README.md` covers installing and using it. The
module docstring in `src/census2_writer/__init__.py` is the design record: the three decisions,
the memory contract, and the constraints on anything added to `harmonize`. The repository root
`CLAUDE.md` gives the Census 2 tentpoles and the code and documentation policy.

This file records the sharp edges: what cost us effort to find, and what is easy to catch again.

## Sharp edges: anndata and zarr

zarr has no S3 client of its own. It defines a store interface, and two adapters exist:
`FsspecStore` over `s3fs`, and `ObjectStore` over `obstore`. This component uses `obstore` for its
smaller dependency footprint and its object listing, so sizing a store needs no second dependency.
Every zarr call that takes a destination must go through `census2_writer.as_store`. A
bare `s3://` string makes zarr fall back to fsspec and fail with `Install s3fs to access S3`, and
a local path hides the mistake because it works.

`obstore` reads credentials from the environment or from instance metadata. It does not resolve
AWS profiles, so `AWS_PROFILE` alone is not enough on a workstation.

anndata consolidates Zarr metadata itself, through `write_zarr(consolidate_metadata=...)`; a
second consolidation is redundant. Consolidated metadata is written last, so a reader that opens
with `zarr.open_consolidated` fails on a partially written store instead of reading part of one.

anndata 0.13 made Zarr v3 the default write format with automatic sharding, so `pyproject.toml`
declares `anndata>=0.13,<0.14`. The upper bound is the constraint that matters, and this is the
only place any component states it: the runtime image resolves its dependencies from this file.
The same release moved `X` into `layers[None]` in memory; on disk `X` stays at `/X`.

`anndata.experimental.read_lazy` would give a lazy whole object, but as of 0.13 it does not
round-trip: anndata cannot write an AnnData whose `obs` is a `Dataset2D`, and `read_lazy` takes no
`chunks` argument. So `build` reads the small slots eagerly and the matrices through
`read_elem_lazy`. Because `X` is then an unevaluated dask graph, the `h5py` handle stays open
until the store is written.

The round trip widens the sparse index arrays of `X` from int32 to int64; the values do not
change. `--check` therefore compares integer arrays by value and not by stored bytes.

A gene-axis query against a CSR matrix costs close to a full scan: selecting any set of genes
reads the whole `indices` array, and no chunk size, shard size, or codec changes that. Serving
those queries efficiently needs a second orientation or a dense two-dimensional chunk grid, which
is a change of output format rather than tuning. Census 1's TileDB-SOMA arrays are configured with
equal 2048 tile extents on roughly 159 million cells by 61,500 genes, which makes that build
obs-major in practice; the format supports other layouts, so that is a property of the build
rather than of TileDB-SOMA. Efficient gene-axis queries are unsolved on both sides, and neither
corpus can be presented as the answer for them.

## Sharp edges: types

`read_elem` and `read_elem_lazy` are typed for every element kind they can return, so both return
a wide union. The annotations here narrow at the boundary with `cast`, which records what we
know that the checker cannot: for a matrix, `read_elem_lazy` returns a dask array, and the memory
contract depends on that.

h5py ships no type information, so `H5Matrix` and `Source` are `Any` with a comment naming what a
caller must pass, and h5py is the only entry in the `ignore_missing_imports` override. Keep it the
only entry: adding a package there discards types we can have. Everything else is typed: dask,
obstore, anndata, and zarr ship `py.typed`, and `pandas-stubs` and `scipy-stubs` are dev
dependencies.

`Result` and `Receipt` are `TypedDict`s rather than named tuples, because both serialize into JSON
and callers index them. A `TypedDict` is closed, which is why `Config` is not one: the config file
format is undecided and must carry keys we have not named. When that format is designed it needs
run-time validation, not annotations, because it is external input.
