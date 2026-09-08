# census2-store

Build one Census 2 store from one CELLxGENE Discover dataset.

A Census 2 store is one AnnData-Zarr store per dataset, written independently and composed at query
time through a catalog. This package holds the data model: what a store contains, and the three
decisions behind it. It is the same code the bulk ingestion pipeline runs, so a store you build
here matches a store the Census publishes.

Two of the three decisions are stubs today. The module docstring is the design record — read
`help(census2_store)` or `src/census2_store/__init__.py` for the reasoning, the memory contract,
and the constraints on anything added to `harmonize`.

## Install

```bash
pip install ./census2-store          # local h5ad to a local store
pip install './census2-store[s3]'    # also write to object storage
```

zarr ships no S3 client, and `as_store` imports `obstore` only when the destination is remote, so
the `s3` extra is unnecessary for a local conversion.

## Use

```python
import census2_store

# build and write in one call; writes nothing if the dataset is not admitted
result = census2_store.create("dataset.h5ad", "store.zarr")
if not result["admitted"]:
    print("skipped:", result["reason"])
```

To inspect or alter the object before it is written, run the two steps yourself. The source must
stay open until the write, because the expression matrix is read lazily:

```python
import h5py

with h5py.File("dataset.h5ad") as f:
    adata, result = census2_store.build(f)
    if adata is not None:
        census2_store.write(adata, "store.zarr")
```

`dest` is a local path, a remote URL such as `s3://bucket/prefix.zarr`, or a
`zarr.abc.store.Store`. A remote URL always goes through `as_store`: a bare `s3://` string makes
zarr fall back to fsspec and fail, while a local path hides the mistake because it works.

## The three decisions

| Decision | State | Where |
|---|---|---|
| Admission — whether a dataset belongs in Census 2 | Stub, admits everything | `admit` |
| Slot policy — which AnnData slots a store carries | Decided: the Census set, or `all_slots` | `build` |
| Harmonization — what makes stores queryable together | Stub, returns its input | `harmonize` |

A rejection is a normal outcome and not an error: `build` returns no AnnData and a reason, and the
caller records that and succeeds.

Ingestion reads exactly one dataset, so a harmonization step that needs a corpus-wide aggregate has
to be redesigned or eliminated. That constraint is what lets many Census releases reuse one store.

## Develop

```bash
pip install './census2-store[dev,s3]'
pytest census2-store
mypy census2-store/src
```

The tests need no network, no credentials, and no container. See the repository root `CLAUDE.md`
for the lint, format, and type-check policy, and `CLAUDE.md` here for the anndata and zarr
behaviour that shaped this code.
