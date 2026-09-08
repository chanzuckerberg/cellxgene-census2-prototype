# census2-store

The data model: what a Census 2 store contains. `README.md` covers installing and using it, and
the module docstring in `src/census2_store/__init__.py` is the design record — the three decisions,
the memory contract, and the constraints on anything added to `harmonize`.

This file records the rakes: what cost us effort to find, and what is easy to step on again. The
repository root `CLAUDE.md` gives the Census 2 tentpoles and the code and documentation policy.

## Rakes: anndata and zarr

zarr has no S3 client of its own. It defines a store interface, and two adapters exist:
`FsspecStore` over `s3fs`, and `ObjectStore` over `obstore`. This component uses `obstore` for its
smaller dependency footprint and its object listing, which saves a second dependency for sizing a
store. Every zarr call that takes a destination has to go through `census2_store.as_store`. A bare
`s3://` string makes zarr fall back to fsspec and fail with `Install s3fs to access S3`, and a
local path hides the mistake because it works.

`obstore` reads credentials from the environment or from instance metadata. It does not resolve
AWS profiles, so `AWS_PROFILE` alone is not enough on a workstation.

anndata consolidates Zarr metadata itself, through `write_zarr(consolidate_metadata=...)`. A second
consolidation is redundant. Consolidated metadata is written last, so a reader that opens with
`zarr.open_consolidated` fails on a partially written store instead of reading part of one.

anndata 0.13 made Zarr v3 the default write format with automatic sharding, so `pyproject.toml`
declares `anndata>=0.13,<0.14`. The upper bound is the constraint that matters, and it is the only
place any component states it: the runtime image resolves its dependencies from this file. The same release moved `X` into `layers[None]` in memory; on disk `X` stays at
`/X`.

`anndata.experimental.read_lazy` looks like the natural way to get a lazy whole object, but as of
0.13 it does not round-trip: anndata cannot write an AnnData whose `obs` is a `Dataset2D`, and
`read_lazy` takes no `chunks` argument. Hence the hybrid of eager small slots and
`read_elem_lazy` matrices.

Because `X` arrives as an unevaluated dask graph, the `h5py` handle has to stay open until the
store is written.

The round trip widens the sparse index arrays of `X` from int32 to int64. The values do not change.
`--check` therefore compares integer arrays by value and not by stored bytes.

A CSR matrix cannot serve gene-axis queries. Selecting genes requires reading the whole `indices`
array, and no chunk size, shard size, or codec changes that. Serving those queries needs a second
orientation or a dense two-dimensional chunk grid, which is a change of output format rather than
tuning.

Census 1's TileDB-SOMA arrays are configured with equal 2048 tile extents on an array of about
159 million cells by 61,500 genes, which makes that build obs-major in practice: selecting a few
genes reads close to the volume of selecting all of them. The format supports other layouts, so
this is a property of that build rather than of TileDB-SOMA. It matters here only because it means
gene-axis queries are unsolved on both sides, and neither corpus can be presented as the answer
for them.

## Rakes: types

`read_elem` and `read_elem_lazy` are typed for every element kind they can return, so both hand
back a wide union. The annotations here narrow at the boundary with `cast`, which records what we
know that the checker cannot: for a matrix, `read_elem_lazy` returns a dask array, and the memory
contract depends on that.

h5py ships no type information, so `H5Matrix` and `Source` are `Any` with a comment naming what a
caller must pass. `pandas-stubs` is a dev dependency, and without it `obs` and `var` check as `Any`
rather than as DataFrames.
