# census2-sdk

Query Census 2 stores through a catalog and return one AnnData.

The SDK runs one query end to end over the stores that `census2-dataset-writer` writes. The
catalog excludes the datasets that cannot match. The SDK opens the remaining stores, reads the
selected cells, and concatenates them by gene position. The result is one AnnData. The filter
grammar, the per-chunk index, and private sources are later work.

## Install

```bash
pip install ./census2-sdk          # local catalog over local stores
pip install './census2-sdk[s3]'    # also read stores and catalogs on S3
```

zarr ships no S3 client. The `s3` extra adds `obstore`, and the SDK uses it only for a remote URI.

## Build the catalog

The catalog is a Parquet table with one row for each store. The Census 2 team publishes it, and the
SDK only reads it. `tools/build_catalog.py` is the prototype stand-in for the catalog builder.
It reads a prefix that holds `<dataset_version_id>.zarr` stores. It reads no receipt.

```bash
eval "$(aws configure export-credentials --profile <profile> --format env)"
python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet
```

Each row holds the identity of the dataset, the organism, the spatial flag, `n_obs`, and the
`gene_reference_version`. It also holds the distinct values of the 18 indexed columns of design 2.2.
Each store's summary record, in the attributes of its root group, gives `dataset_id`,
`dataset_version_id`, and `gene_reference_version`. The `census2-dataset-writer` README describes
the record. The builder reads the other fields from the store's `obs` and `uns`. It refuses a store
that has no summary record.

Each row also holds the cell counts of its store, in the `cell_counts` column. The column has one
entry for each combination of the counted values, with the number of cells in `n_cells`. The counted
values are assay, cell type, disease, sex, tissue, and `tissue_general` (each as the label and the
ontology term), and `is_primary_data`. The builder counts in the same pass that reads the values for the catalog row.
For the 1,258 stores written from Census 1 LTS `2025-11-08`, the column holds
117,712 entries, and the catalog is 1,014 KiB.

## Build a release

A release is a directory that pins the live catalog at one point in time (design 3.3).
`tools/build_release.py` is the prototype stand-in for the release builder. It reads the live
catalog one time, the `var` of each store, and the per-gene nnz file next to each store. It reads
no expression data.

```bash
python census2-sdk/tools/build_release.py catalog.parquet releases/ 2026-09-01
```

```text
releases/2026-09-01/
├── catalog.parquet               the live catalog content, unchanged
├── presence/<organism>.parquet   dataset_id, dataset_version_id, measured
├── aggregates/var.parquet        organism, feature_id, nnz, n_measured_obs
└── manifest.json                 the last write
```

- `catalog.parquet` is the content that the builder read. The snapshot keeps the `store_uri`
  column, so a release query opens the stores at the same locations.
- Each presence file has one row for each store of the organism. `measured` lists the positions in
  `var` of the genes that the dataset measured.
- `aggregates/var.parquet` has one row for each gene of each organism. `nnz` is the number of
  cells with a nonzero count of the gene. `n_measured_obs` is the number of cells of the datasets
  that measured the gene. For `nnz`, the builder adds the per-gene nnz files that the ingest task
  writes next to the stores (design 3.4). A store without its file stops the build.
- `manifest.json` pins the snapshot, the presence files, the aggregates, and the
  `dataset_version_id` of each store. The section "Query a release" lists its fields.

The caller gives the name of the release. The builder writes the manifest last, so a directory with
a manifest is a complete release. Each release is a new directory, and the builder refuses a
directory that exists. A build that stopped leaves a directory without a manifest. Remove that
directory before the next build. The builder writes the directory to the local file system. The
live catalog and the stores can be remote.

A release holds one gene reference version for each organism (design 1.5). The builder refuses a
catalog in which the rows of one organism record two versions.

The builder adds the nnz arrays by position, and the presence files hold positions. A position
names one gene in each store only when the stores hold the same genes in the same order. Thus the
builder compares the ordered gene identifiers of each store with those of the first store of the
organism. It stops at a store that differs. For this check, the builder reads the index of `var`
of each store, which is one string for each gene (78,894 for a human store).

The builder checks each nnz file before it adds the array. The array must be `int64` with one value
for each gene of `var`. A gene that the dataset did not measure must have zero, and no value can be
above the `n_obs` of the store.

This version of the release builder makes no dictionary file, compares no store with an inventory,
and adds no entry to a release index. The live catalog has no dictionary, and a store has no
inventory.

## Query

```python
import census2_sdk

census = census2_sdk.open_census(sources=["catalog.parquet"])
adata = census.get_anndata("homo_sapiens", obs_filter={"tissue": "lung"})
```

`obs_filter` maps a column to one value or to a list of values. The SDK selects a cell when every
column holds one of its values. The columns are the 18 indexed columns and `dataset_id`. The
`disease` columns are multi-valued, and a value matches when it is one of the members.

The result holds:

- `X`, the counts
- `layers["normalized"]`, with `layers=["normalized"]`: the Census 1 normalized values, which the
  SDK calculates from the counts on read. The SDK does not round them to 15 mantissa bits as
  Census 1 does, so they differ from Census 1 by less than 2^-15 relative.
- `var`, the pinned gene reference for the organism (GENCODE v48 for human, 78,894 rows)
- `obs`, the columns that all selected stores hold, which include `dataset_id` and
  `observation_joinid`
- `obs_names`, `<dataset_id>_<observation_joinid>`
- `uns["census2"]`, the sources, the releases among them, the count of catalog rows for the
  organism, and the count of stores opened

When no cell matches, `get_anndata` returns None.

## Query a release

A source is a catalog or the manifest of a release. A release is a directory that holds a snapshot
of the catalog, `catalog.parquet`, and a `manifest.json` (design 3.3). The manifest pins the
snapshot by its SHA-256 and lists the `dataset_version_id` of each store.

```python
census = census2_sdk.open_census(
    sources=["releases/2026-09-01/manifest.json", "my_lab_catalog.parquet"]
)
```

The SDK reads `catalog.parquet` from the directory of the manifest, and then queries that snapshot
as it queries a catalog. Thus `get_anndata`, `count_cells`, and the checks before concatenation are
the same for a release, a catalog, and a mix of the two. A later write to the live catalog does not
change the snapshot. `uns["census2"]["releases"]` maps the name of each release among the sources
to the SHA-256 of its manifest. The SDK refuses two sources that name the same release.

The manifest is UTF-8 JSON with five fields, and `census2_sdk.release` defines them:

| Field | Content |
| --- | --- |
| `version` | `1`, the format of the manifest |
| `release` | the name of the release, which is the name of its directory |
| `gene_reference_version` | the gene reference version of each organism |
| `stores` | one `{"dataset_version_id": ...}` for each store |
| `files` | the SHA-256 of each other file of the release, by its path in the directory |

The manifest pins the identity of each store and not its content. It holds no `writer_version` and
no content digest of a store, because no store and no catalog row records either value. The design
record reserves the format of the manifest, so these fields are the format of version 1.

`open_census` raises `release.ReleaseError`, a `ValueError`, in three cases: a malformed manifest,
an absent snapshot, and a snapshot that differs from the manifest. A source that does not exist
raises `FileNotFoundError`. `open_census` opens no store. A store that does not exist fails the
query that reads it, and the error names the store, the dataset, and the source.

## Count

```python
n = census.count_cells("homo_sapiens", {"cell_type": "neutrophil"})
n = census.count_cells("homo_sapiens", {"tissue": "blood"}, primary_only=True)
```

`count_cells` sums the cell counts in the catalog rows and opens no store. A cell matches by the
rule of `get_anndata`, and the catalog excludes the same datasets. Thus, without `primary_only`, the
count equals the `n_obs` of the result of `get_anndata` when there is one. A count joins no
matrices, so it does not apply the checks before concatenation. The filter takes the counted columns
and `dataset_id`, and a filter on another indexed column raises `ValueError`. With `primary_only`,
the count includes only the cells whose `is_primary_data` is true.

On the 1,258 stores, the call counts 374,463 human neutrophils, and Census 1 gives the same count.
On a laptop, `open_census` and the first count take about 40 ms together, and a later count takes
15 ms.

The writer projects each store onto the gene reference of its organism. Thus every store of one
organism has the same genes in `var`, in the same order. The SDK takes `var` from the stores and
downloads no reference. A store's `var` also marks the dataset's assayed genes in `measured`, and
the result leaves that column out.

## Checks before concatenation

The SDK applies the checks of design 4.3 to the selected catalog rows and reports the conflicting
values when a check fails:

- one organism, given by the query
- every row records a `gene_reference_version`, and it equals the version the SDK pins for the
  organism
- every opened store has the same genes in `var`; the SDK checks this when it reads the stores
- no spatial store together with a non-spatial store
- a `dataset_id` in one source only, checked when the sources open

## Compare with Census 1

`tools/compare_census1.py` runs one query through the SDK and through Census 1 release 2025-11-08.
It compares the cell sets, the gene sets, and the values on the intersection. A store written from a
re-curated dataset version can hold other cells, so a cell-set difference is a finding to record.

```bash
uv run --python 3.12 --with cellxgene-census==1.18.0 --with './census2-sdk[s3]' \
    census2-sdk/tools/compare_census1.py --catalog catalog.parquet --dataset-id <dataset_id> \
    --out census2-sdk/tools/results/compare_<dataset_id>.json
```

`--out` writes the figures as JSON. `tools/results/` holds the runs that a pull request cites.

## Validation harness

`tools/harness/` holds the pieces of the validation harness that the SDK's
differential run needs. They are the fixed query corpus, the canonical result form, the
comparator, and the oracle runner that records Census 1's answers.
`tools/harness/README.md` has the pins and the steps. `tools/harness/sdk_adapter.py` is the
SDK's adapter to the harness.

## Develop

```bash
cd census2-sdk
pip install -e ../census2-dataset-writer -e '.[dev,s3]'
pytest
mypy
```

The catalog builder in `tools/` imports `census2-dataset-writer` for the summary record: its type,
its location in the store, and the check of its fields. The release builder imports it for the name
of the `measured` column and the name of the per-gene nnz file. Thus the tests and mypy need the
writer installed. The SDK package itself does not depend on it.

Run both from this directory. `pytest` and `mypy` read their settings from `pyproject.toml` here. The tests write
small stores to a temporary directory and need no network and no credentials. Lint and format with `pipx run ruff check --fix .` and `pipx run ruff format .` from the
repository root. CI runs the unpinned `ruff check` and `ruff format --check` on this directory. See `AGENTS.md` here for the anndata and zarr
behavior that shaped this code.
