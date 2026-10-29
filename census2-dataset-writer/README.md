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
| Admission — does the dataset belong in Census 2 | Minimal stub: an organism with a gene reference | `admit` |
| Cell filters — which cells of an admitted dataset a store keeps | Decided: the Census 1 filters on tissue type, assay, and primary data on a spatial assay | `cell_mask` |
| Slot policy — which AnnData slots a store carries | Decided: the Census set, or `all_slots`; the carried `uns` keys come from the `uns_keys` config | `build` |
| Harmonization — what makes stores queryable together | Partly made: the `obs` schema, the counts matrix, and the gene projection | `harmonize` |

A store holds the dataset's cells, with the counts on the organism's pinned gene reference:

- `X` holds the counts: `raw.X` when the source has `raw`, else `X`. The store has no `raw`.
- `var` is the reference (GENCODE v48 for human, 78,894 rows; GENCODE M37 for mouse), in reference
  order. Every store of one organism has the same genes in `var`, so a reader concatenates stores
  by position. A gene that the dataset did not measure is a zero column. The `measured` column of
  `var` marks the genes that the dataset measured, because a zero column cannot show that.
- The writer drops the features that are not genes of the organism, such as ERCC spike-ins. A
  gene that is not in the reference fails the dataset. The writer does not map genes across
  reference versions.
- The result reports the reference version and the dropped features. The store's summary record
  holds the reference version.

The organism criterion rejects a dataset whose organism has no gene reference.

The cell filters then drop the cells that Census 1 drops. A store keeps a cell when its
`tissue_type` is tissue or organoid and its assay is on the allowlist. On a spatial assay, the cell
must also be primary data (`is_primary_data`). The configuration holds the values (`tissue_types`,
`assay_allowlist`, `primary_data_assays`); the allowlist is the Census 1 assay list, from
`census2_dataset_writer.assays`. The result counts the cells that each filter drops. A dataset with
no cell left is rejected. The reason names the filters and, for the assay filter, the assays
present. The assay filter thus also keeps the bulk ingestion workflow away from modalities a Census
is unlikely to carry. On Census 1 LTS `2025-11-08`, the filters drop 32,069 cells from 6 of 1,258
datasets. They drop cell culture cells from 3 human datasets and Patch-seq cells from 3 mouse
datasets.

A rejection is a normal outcome, not an error: `build` returns no AnnData and a reason, and the
caller records that and succeeds.

Ingestion reads one dataset at a time and not the corpus. A harmonization step that needs a
corpus-wide aggregate must be redesigned or eliminated. This constraint lets many Census releases
reuse one store.

## obs schema

A store's `obs` holds a fixed set of columns in a fixed order. These are `dataset_id`, the 19 CxG
columns that Census 1 copies, the `tissue_general` pair, and five statistics of the counts.
`store_format.OBS_COLUMNS` lists them, each with its `ColType`.

- Each column has a fixed type. `is_primary_data` is `bool`, and `observation_joinid` is plain
  strings. `nnz` and `n_measured_vars` are `int64`, and `raw_sum`, `raw_mean_nnz`, and
  `raw_variance_nnz` are `float64`. Every other column is a categorical whose categories are the
  observed values, sorted, as strings.
- No value is null. A source that lacks a column, holds a null, has an `is_primary_data` that is
  not boolean, or has an `observation_joinid` that is not unique fails the build.
- A store does not yet hold any other source column. The result names each one that it drops,
  in `dropped_obs_columns`.
- `tissue_general` and `tissue_general_ontology_term_id` are the Census 1 roll-up of
  `tissue_ontology_term_id` to a high-level tissue, through `cellxgene-ontology-guide` at the
  version that Census 1 pins. The package holds its ontologies and reads no network. The receipt
  records its version.
- `raw_sum`, `nnz`, `raw_mean_nnz`, `raw_variance_nnz`, and `n_measured_vars` are the Census 1
  statistics of the counts. The writer computes them as Census 1 does, before the gene projection,
  so they are bit-identical to Census 1. `n_measured_vars` is the number of genes in `measured`.
  The writer reads `X` once for the statistics and once more for the write.
- Census 1 also holds `soma_joinid`, a corpus-wide number that a store for one dataset cannot
  assign.

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
# {'recipe': 'v2', 'slot_policy': 'census', 'zarr_format': 3,
#  'census2_dataset_writer_version': ..., 'anndata_version': ...,
#  'zarr_version': ..., 'cellxgene_ontology_guide_version': ...}
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
