# census2-sdk

Query Census 2 stores through a catalog and return one AnnData.

The SDK is the skeleton of plan Phase 6, tracker 6.1. It runs one query end to end over the
stores that `census2-dataset-writer` writes. The catalog excludes the datasets that cannot match.
The SDK opens the remaining stores, reads the selected cells, and concatenates them by gene
position. The result is one AnnData. The filter grammar, the per-chunk index, and private sources
are later work.

## Install

```bash
pip install ./census2-sdk          # local catalog over local stores
pip install './census2-sdk[s3]'    # also read stores and catalogs on S3
```

zarr ships no S3 client. The `s3` extra adds `obstore`, and the SDK uses it only for a remote URI.

## Build the catalog

The catalog is a Parquet table with one row for each dataset. The Census 2 team publishes it, and
the SDK only reads it. `tools/build_catalog.py` is the Iteration 1 stand-in for the catalog builder.
It reads one or more prefixes that hold `<dataset_version_id>.zarr` stores, and reads up to
`--workers` stores at the same time. It reads no receipt.

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
values are assay, cell type, disease, sex, and tissue (each as the label and the ontology term), and
`is_primary_data`. The builder counts in the same pass that reads the values for the catalog row.
For the 1,258 stores of the Census 1 LTS `2025-11-08` run, the column holds 117,740 entries, and the
catalog is 925 KiB.

### Compacted stores

A compacted store is a second copy of several datasets in one store; see
[`compact-nf`](../compact-nf/README.md). `--compacted` names a prefix of compacted stores:

```bash
python census2-sdk/tools/build_catalog.py s3://bucket/stores/ catalog.parquet \
    --compacted s3://bucket/compacted/
```

The builder then fills three columns of each member's row: `compacted_store_uri`, and
`compacted_row_start` and `compacted_row_stop`, the rows of the dataset in that store. The row
keeps its per-dataset `store_uri`, and all its other fields come from the per-dataset store as
before. Exclusion therefore stays per dataset. The three columns are null for a dataset that no
compacted store holds, and a catalog without compacted stores has them too.

The builder refuses the catalog when a compacted store does not agree with the per-dataset stores:

- each member has a per-dataset store, with the same `dataset_id`, number of cells, gene reference,
  `observation_joinid` digest, and number of stored values of `X`
- the compacted store holds the member's `dataset_id` and `observation_joinid` values at the
  member's rows, and that number of stored values of `X`
- the members have one organism, one gene reference, and one spatial value
- a dataset is in one compacted store at most

These checks read two `obs` columns of each compacted store and some values of `X/indptr`, not the
matrix. The `compact-nf` pipeline can compare the whole matrix when it writes the store.

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
- `var`, the pinned gene reference for the organism (GENCODE v48 for human, 78,894 rows)
- `obs`, the columns that all selected stores hold, which include `dataset_id` and
  `observation_joinid`
- `obs_names`, `<dataset_id>_<observation_joinid>`
- `uns["census2"]`, the sources, the count of catalog rows for the organism
  (`datasets_in_catalog`), the count of datasets read (`datasets_read`), and the count of stores
  opened (`stores_opened`)

When no cell matches, `get_anndata` returns None.

The SDK reads a dataset from the compacted store that its catalog row names. It opens each store
once for all of its selected datasets, and reads only the rows of each dataset from `obs` and `X`.
The result is the same as from the per-dataset stores: the same cells in the same order, the same
`obs` columns with the same types and categories, and the same `X`.
`open_census(sources, compacted=False)` reads the per-dataset stores, which stay the primary copy.

`max_workers` applies twice: the SDK opens up to `max_workers` stores at the same time, and then
reads up to `max_workers` datasets at the same time. The datasets of one compacted store thus read
at the same time as well.

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
organism has the same `var`, in the same order. The SDK takes `var` from the stores and downloads
no reference.

## Checks before concatenation

The SDK applies the checks of design 4.3 to the selected catalog rows and reports the conflicting
values when a check fails:

- one organism, given by the query
- every row records a `gene_reference_version`, and it equals the version the SDK pins for the
  organism
- every opened store has the same `var`; the SDK checks this when it reads the stores
- no spatial store together with a non-spatial store
- a `dataset_id` in one source only, checked when the sources open

## Compare with Census 1

`tools/compare_census1.py` runs one query through the SDK and through Census 1 release 2025-11-08.
It compares the cell sets, the gene sets, and the values on the intersection. The Iteration 1
stores come from re-curated dataset versions, so a cell-set difference is a finding to record.

```bash
uv run --python 3.12 --with cellxgene-census==1.18.0 --with './census2-sdk[s3]' \
    census2-sdk/tools/compare_census1.py --catalog catalog.parquet --dataset-id <dataset_id> \
    --out census2-sdk/tools/results/compare_<dataset_id>.json
```

`--out` writes the figures as JSON. `tools/results/` holds the runs that a pull request cites.

## Validation harness

`tools/harness/` holds the pieces of the Phase 2 validation harness that the SDK's
differential run needs. They are the fixed query corpus, the canonical result form, the
comparator, and the oracle runner that records Census 1's answers.
`tools/harness/README.md` has the pins and the steps. The SDK's adapter to the harness
is tracker task 6.2.

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
