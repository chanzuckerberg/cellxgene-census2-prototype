# census2-dataset-writer

Convert one CELLxGENE Discover dataset into one Census 2 store.

A Census 2 store is an AnnData-Zarr store that holds a single dataset. Stores are written
independently and composed at query time through a catalog. This package is the data model: it
decides what a store contains. The bulk ingestion pipeline
([`bulk-ingest-nf`](../bulk-ingest-nf/README.md)) runs this same code, so a store built here
matches a store that pipeline writes. The module docstring in
[`src/census2_dataset_writer/__init__.py`](src/census2_dataset_writer/__init__.py) is the design
record: the reasoning, the memory contract, and the constraints on new harmonization steps.

## Install

```bash
pip install ./census2-dataset-writer          # local h5ad to a local store
pip install './census2-dataset-writer[s3]'    # also write to object storage
```

The `s3` extra is unnecessary for a local conversion.

## Use

```python
import census2_dataset_writer

# build and write in one call; writes nothing if the dataset is not admitted
result = census2_dataset_writer.create("dataset.h5ad", "store.zarr", refdata="refdata/")
if not result["admitted"]:
    print("skipped:", result["reason"])
```

`build` and `write` can also be called separately; the module docstring covers that.

`refdata` is a directory that holds the pinned gene reference files, `genes_homo_sapiens.csv.gz`
and `genes_mus_musculus.csv.gz` from the `single-cell-curation` tag `v7.0.0`. The writer does not
download them. Stage them once with
[`bulk-ingest-nf/tools/stage-refdata.py`](../bulk-ingest-nf/tools/stage-refdata.py), which checks
each file against the SHA-256 digest in `DEFAULTS["gene_references"]`.

The destination is a local path, a remote URL such as `s3://bucket/prefix.zarr`, or a
`zarr.abc.store.Store`. For a remote destination, export AWS credentials into the environment
first; `AWS_PROFILE` alone does not work.

## The three decisions

| Decision | State | Where |
|---|---|---|
| Admission — does the dataset belong in Census 2 | Minimal stub: the `assay_allowlist`, and an organism with a gene reference | `admit` |
| Slot policy — which AnnData slots a store carries | Decided: the Census set, or `all_slots`; the carried `uns` keys come from the `uns_keys` config | `build` |
| Harmonization — what makes stores queryable together | Partly made: the `obs` columns that Zarr cannot name, and the gene projection | `harmonize` |

A store today holds the dataset's cells, with the counts on the organism's pinned gene reference:

- `X` holds the counts: `raw.X` when the source has `raw`, else `X`. The store has no `raw`.
- `var` is the reference (GENCODE v48 for human, 78,894 rows; GENCODE M37 for mouse), in reference
  order. Every store of one organism has the same `var`, so a reader concatenates stores by
  position. A gene that the dataset did not measure is a zero column.
- The features that are not genes of the organism, such as ERCC spike-ins, are dropped. A gene
  that is not in the reference fails the dataset; the writer does not map genes across reference
  versions.
- `uns["gene_reference_version"]` records the reference, and the result reports it with the dropped
  features.

The assay criterion rejects a dataset with no cell on the assay allowlist, which defaults to the
Census 1 assay list; see `census2_dataset_writer.assays`. It is there to keep the bulk ingestion
workflow away from modalities a Census is unlikely to carry, not to settle what Census 2 admits.
The organism criterion rejects a dataset whose organism has no gene reference.

A rejection is a normal outcome, not an error: `build` returns no AnnData and a reason, and the
caller records that and succeeds.

Ingestion reads one dataset at a time and not the corpus. A harmonization step that needs a
corpus-wide aggregate must be redesigned or eliminated. This constraint lets many Census releases
reuse one store.

## Receipts

A receipt is the JSON a pipeline writes beside a store. `Receipt` types its shape, `provenance()`
supplies the part this package knows, and the pipeline adds where the dataset came from and what
it did.

```python
census2_dataset_writer.provenance()
# {'recipe': 'v2', 'slot_policy': 'census', 'zarr_format': 3,
#  'census2_dataset_writer_version': ..., 'anndata_version': ...,
#  'zarr_version': ...}
```

Registering a store in a catalog, with metadata a reader outside the run can rely on, is a
separate job and is not implemented yet.

## Develop

```bash
cd census2-dataset-writer
pip install -e '.[dev,s3]'
pytest
mypy
```

Run pytest and mypy from this directory; mypy finds its configuration only here. Install
editable, or the tests import the copy in site-packages. The tests need no network, no
credentials, and no container.

See the repository root [`CLAUDE.md`](../CLAUDE.md) for the lint, format, and type-check policy,
and [`CLAUDE.md`](CLAUDE.md) here for the anndata and zarr behavior behind this code.
