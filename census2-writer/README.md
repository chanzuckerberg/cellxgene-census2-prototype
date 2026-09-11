# census2-writer

Convert one CELLxGENE Discover dataset into one Census 2 store.

A Census 2 store is an AnnData-Zarr store holding a single dataset, written independently and
composed at query time through a catalog. This package is the data model: what a store contains,
and the [three decisions](#the-three-decisions) behind it. The bulk ingestion pipeline runs this
same code, so what you build here is what that pipeline writes.

Two of those three decisions are stubs today, so a store is currently the dataset's Census slot set
written to Zarr, and no more than that. Both stubs will change what a store contains.

It reads one dataset at a time and never the corpus. That is a constraint the data model lives
with, not a feature it lacks. The module docstring in
[`src/census2_writer/__init__.py`](src/census2_writer/__init__.py) is the design record: the
reasoning, the memory contract, and what a new harmonization step has to respect.

## Install

```bash
pip install ./census2-writer          # local h5ad to a local store
pip install './census2-writer[s3]'    # also write to object storage
```

zarr ships no S3 client, and `as_store` imports `obstore` only when the destination is remote, so
the `s3` extra is unnecessary for a local conversion.

## Use

```python
import census2_writer

# build and write in one call; writes nothing if the dataset is not admitted
result = census2_writer.create("dataset.h5ad", "store.zarr")
if not result["admitted"]:
    print("skipped:", result["reason"])
```

To inspect or alter the object before it is written, run the two steps yourself. The source must
stay open until the write, because the expression matrix is read lazily:

```python
import h5py

with h5py.File("dataset.h5ad") as f:
    adata, result = census2_writer.build(f)
    if adata is not None:
        census2_writer.write(adata, "store.zarr")
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

## Receipts

A receipt is the JSON a writer publishes beside a store, and `Receipt` types it here rather than in
any one pipeline, so that a reader parses the same keys whichever pipeline wrote the store.
`provenance()` supplies the part this package knows — the recipe, the slot policy, the Zarr format,
and the versions that ran — read from the installed distributions, so a receipt cannot disagree
with the code that produced it. A pipeline adds where the dataset came from and what it did.

```python
census2_writer.provenance()
# {'recipe': 'v1', 'slot_policy': 'census', 'zarr_format': 3,
#  'census2_writer_version': ..., 'anndata_version': ..., 'zarr_version': ...}
```

## Develop

```bash
cd census2-writer
pip install -e '.[dev,s3]'
pytest
mypy
```

Run both from this directory. `mypy` reads its settings from `pyproject.toml` here, so invoking it
from the repository root silently drops `strict` and reports an `h5py` error it should ignore.
Install editable, or the tests import the copy in site-packages and keep passing against code you
have already changed. CI installs non-editable on purpose, because there it should test what the
image ships.

The tests need no network, no credentials, and no container. See the repository root
[`CLAUDE.md`](../CLAUDE.md) for the lint, format, and type-check policy, and
[`CLAUDE.md`](CLAUDE.md) here for the anndata and zarr behaviour that shaped this code.
