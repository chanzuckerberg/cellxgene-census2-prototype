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
- `uns["census2"]`, the sources, the count of catalog rows for the organism, and the count of stores opened

When no cell matches, `get_anndata` returns None.

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

## The obs table

This part is a spike. It is not intended to merge in its current form; see
[`spikes/obs-iceberg.md`](../spikes/obs-iceberg.md).

The writer puts an `obs` fragment next to each store: a Parquet file with the store's `obs`
columns, the store row ordinal, and the range of each cell's values in `X`. The `obs` table of an
organism is one [Apache Iceberg](https://iceberg.apache.org/) table whose data files are the
fragments of the catalog's stores, in place. `tools/build_obs_table.py` writes the table metadata
under the table location and needs no Iceberg catalog service:

```bash
pip install './census2-sdk[s3,iceberg]'
python census2-sdk/tools/build_obs_table.py catalog.parquet s3://my-bucket/census/obs/homo_sapiens
```

With the table, one scan selects the cells, and the SDK opens only the stores that hold a selected
cell. It reads the result `obs` from the scan and the `X` values of the selected cells from their
ranges, so it reads no store `obs` and no `X/indptr`:

```python
census = census2_sdk.open_census(
    ["catalog.parquet"], obs_tables=["s3://my-bucket/census/obs/homo_sapiens"]
)
adata = census.get_anndata(
    "homo_sapiens",
    "tissue IN ('lung', 'lung parenchyma') AND cell_type = 'macrophage' AND is_primary_data = true",
    max_workers=8,
)
obs = census.get_obs("homo_sapiens", {"tissue": "lung"}, columns=["donor_id"])
counts = census.count_obs("homo_sapiens", "disease = 'COVID-19'", by=["donor_id"])
```

The filter is the mapping filter or an Iceberg expression: a pyiceberg expression object, or a
string that pyiceberg parses. The catalog first excludes the datasets that cannot match the
expression, and the scan then reads the fragments of the rest. On the `disease` columns, `=`,
`!=`, `IN`, and `NOT IN` mean membership. `census2_sdk.expressions` states the rules, and
`expressions.from_value_filter` translates a Census 1 `value_filter` string. The result equals the
result through the stores, value for value, also with `layers=["normalized"]`. A categorical
column holds the labels that the result uses, where the result through the stores holds the labels
of every store that it opened. The snapshot of the table names the digest of its catalog, and
`open_census` refuses a table of another catalog.

A numeric column of the fragment, such as `nnz` or `raw_sum`, takes a range term:
`"nnz >= 500 AND tissue_general = 'lung'"`. The catalog records no numeric range, so only the table
tests such a term. `get_obs(..., dtype_backend="pyarrow")` returns `pd.ArrowDtype` columns over
the scan's Arrow data, which pandas does not copy; the default gives categoricals.

The table path does not read the `var` index of each store that it opens. Every catalog row must
record the pinned `gene_reference_version`, and the writer projects every store onto that
reference, so the SDK checks only that each store has the number of genes of the first store, and
takes `var` from the first store.

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
its location in the store, and the check of its fields. The tests and mypy therefore need the
writer installed. The SDK package itself does not depend on it.

Run both from this directory. `pytest` and `mypy` read their settings from `pyproject.toml` here. The tests write
small stores to a temporary directory and need no network and no credentials. Lint and format with `pipx run ruff check --fix .` and `pipx run ruff format .` from the
repository root. CI runs the unpinned `ruff check` and `ruff format --check` on this directory. See `AGENTS.md` here for the anndata and zarr
behavior that shaped this code.
