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
result = census2_dataset_writer.create(
    "dataset.h5ad",
    "store.zarr",
    dataset_id="...",
    dataset_version_id="...",
    refdata="refdata/",
)
if not result["admitted"]:
    print("skipped:", result["reason"])
```

A Discover h5ad does not hold its own identifiers, so the caller gives `dataset_id` and
`dataset_version_id`. They go into the store's summary record, and `dataset_id` also goes into
`obs`. `build` and `write` can also be called separately; the module docstring covers that.

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
| Harmonization — what makes stores queryable together | Partly made: the `obs` columns that Zarr cannot name, the counts matrix, and the gene projection | `harmonize` |

A store holds the dataset's cells, with the counts on the organism's pinned gene reference:

- `X` holds the counts: `raw.X` when the source has `raw`, else `X`. The store has no `raw`.
- `var` is the reference (GENCODE v48 for human, 78,894 rows; GENCODE M37 for mouse), in reference
  order. Every store of one organism has the same `var`, so a reader concatenates stores by
  position. A gene that the dataset did not measure is a zero column.
- The writer drops the features that are not genes of the organism, such as ERCC spike-ins. A
  gene that is not in the reference fails the dataset. The writer does not map genes across
  reference versions.
- The result reports the reference version and the dropped features. The store's summary record
  holds the reference version.

The assay criterion rejects a dataset with no cell on the assay allowlist, which defaults to the
Census 1 assay list; see `census2_dataset_writer.assays`. It is there to keep the bulk ingestion
workflow away from modalities a Census is unlikely to carry, not to settle what Census 2 admits.
The organism criterion rejects a dataset whose organism has no gene reference.

A rejection is a normal outcome, not an error: `build` returns no AnnData and a reason, and the
caller records that and succeeds.

Ingestion reads one dataset at a time and not the corpus. A harmonization step that needs a
corpus-wide aggregate must be redesigned or eliminated. This constraint lets many Census releases
reuse one store.

## Summary record

Every store holds a summary record: the facts about the store that a catalog needs and that the
store does not otherwise hold. The record is a JSON object in the attributes of the store's root
Zarr group, so it is in the root `zarr.json`, beside anndata's own attributes:

```json
"cellxgene-census2": {
  "summary": {
    "version": 1,
    "dataset_id": "...",
    "dataset_version_id": "...",
    "gene_reference_version": "GENCODE v48"
  }
}
```

- `version` gives the format of the record. Its value is `SUMMARY_VERSION`, and a reader refuses a
  version that it does not know.
- `SummaryRecord` types the record, and `read_summary(dest)` reads it. A record comes from a
  store, so `read_summary` checks it at run time with `validate_summary`. That function raises
  `SummaryRecordError` for a missing, unknown, or wrongly typed field. `write` checks its record
  in the same way, before it writes anything.
- `gene_reference_version` is the label of the gene reference in the store's `var`, from
  `genes.REFERENCES`. No ontology names an annotation release. The label follows the GENCODE
  release that the pinned `single-cell-curation` tag declares.
- The store holds the organism and the curation schema version in `uns`, so the record does not
  copy them.

`write` adds the record after anndata writes the store and its consolidated metadata, so the
record is the last write.

`build` adds `dataset_id` to `obs` as the first column, so that each cell names its dataset when
stores are concatenated. The column is categorical with one category. A source `obs` that already
has a `dataset_id` column is an error: `build` raises `ValueError`, because that name belongs to the
store. `write` refuses an AnnData whose `dataset_id` column does not hold the summary record's
`dataset_id`.

## Receipts

A receipt is the JSON a pipeline writes beside a store, about one run of the pipeline. `Receipt`
types its shape, `provenance()` supplies the part this package knows, and the pipeline adds where
the dataset came from and what it did. For an admitted dataset it holds a copy of the store's
summary record, so that a report on the run opens no store. Nothing outside the pipeline reads a
receipt.

```python
census2_dataset_writer.provenance()
# {'recipe': 'v1', 'slot_policy': 'census', 'zarr_format': 3,
#  'census2_dataset_writer_version': ..., 'anndata_version': ...,
#  'zarr_version': ...}
```

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

See the repository root [`AGENTS.md`](../AGENTS.md) for the lint, format, and type-check policy,
and [`AGENTS.md`](AGENTS.md) here for the anndata and zarr behavior behind this code.
