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
hs = census.homo_sapiens
adata = census.get_anndata(hs.obs.cell_type.isin(["neutrophil"]) & (hs.obs.nnz >= 500))
adata = census.get_anndata(hs.where({"tissue": "lung"}))
adata = census.get_anndata(hs)  # every human cell
```

A query selects cells of one organism. `census.homo_sapiens` and `census.mus_musculus` give the
organisms. Each column of `hs.obs` has a class for its type, with the pandas operators that suit
the type:

| Column | Operators |
|---|---|
| labels, such as `tissue` and `donor_id` | `==`, `!=`, `isin`, `.str.startswith`, `.endswith`, `.contains`, `.match`, `.fullmatch` |
| `disease`, `disease_ontology_term_id` (multi-valued) | `has_any([...])`, `has_all([...])` |
| `is_primary_data` | the column as a query, `~`, `==` a bool |
| `nnz` and the other count statistics | `==`, `!=`, `<`, `<=`, `>`, `>=`, `between` |
| every column | `isna`, `notna` |

`&`, `|`, and `~` combine queries, as for a pandas mask. `and`, `or`, and `not` raise `TypeError`.
`hs.where` takes a mapping of a catalog column to one value or to a list of values: the 18 indexed
columns and `dataset_id`. `hs.obs["name"]` takes a column that the static schema does not name.
`hs.expression(expr, columns=[...])` takes a `pyarrow.compute.Expression` as it is, with the names of
its columns, for a test that the builder does not have, such as arithmetic between columns.

Each operator gives a pyarrow expression, and pyarrow evaluates it on the columns of each store.
Thus the null semantics are pyarrow's, which are SQL's three-valued logic: a comparison with a null
gives null, and the SDK keeps only the cells where the query is true. `~(obs.x == "a")` drops the
cells where `x` is null, unlike pandas. `isin` follows the same logic. `isna` is true for a null and
for a float NaN, as in pandas. The writer stores no null in the columns of the store format, so the
rule matters only for other columns. A regular expression is RE2, not Python's `re`.

The query also records the columns that it reads and its hints: the allowed values of each catalog
column, from its top-level `&` of `==`, `isin`, and `has_any` terms. Two tests then run on the
catalog. The hints exclude a dataset whose value sets prove that no cell matches them. Then pyarrow
evaluates the query on every combination of each remaining dataset's value sets, with a null for a
column that the catalog does not know, such as `nnz`. A dataset where every combination gives false
is excluded, so `|`, `~`, and the string methods exclude datasets too. A dataset where every
combination gives true, with no unknown column, is whole: the SDK keeps every cell of the store and
reads no filter column. That proof relies on the store format's rule that an indexed column holds no
null. A query with `isna`, `notna`, or `hs.expression` keeps only the hint test.
`census.candidates(query)` shows the catalog rows that a query opens, with a `whole` column.

For each other store, the SDK reads only the query's columns, as codes for a categorical column,
and evaluates the expression with pyarrow. A store that keeps no cell reads nothing more.

The result holds:

- `X`, the counts
- `layers["normalized"]`, with `layers=["normalized"]`: the Census 1 normalized values, which the
  SDK calculates from the counts on read. The SDK does not round them to 15 mantissa bits as
  Census 1 does, so they differ from Census 1 by less than 2^-15 relative.
- `var`, the pinned gene reference for the organism (GENCODE v48 for human, 78,894 rows)
- `obs`, the columns that all stores with a selected cell hold, which include `dataset_id` and
  `observation_joinid`. A categorical column holds only the categories that the selected cells use.
- `obs_names`, `<dataset_id>_<observation_joinid>`
- `uns["census2"]`, the sources, the count of catalog rows for the organism, the count of stores
  opened, and the count of them that the catalog proved whole

When no cell matches, `get_anndata` returns None. `census.get_obs(query, columns=[...])` gives the
`obs` of the same cells and reads no `X`.

## Count

```python
n = census.count_cells(hs.obs.cell_type == "neutrophil")
n = census.count_cells((hs.obs.tissue == "blood") & hs.obs.is_primary_data)
```

`count_cells` sums the cell counts in the catalog rows and opens no store. A cell matches by the
rule of `get_anndata`, and the catalog excludes the same datasets. Thus the count equals the
`n_obs` of the result of `get_anndata` when there is one. A count joins no matrices, so it does
not apply the checks before concatenation. The query must be exactly its hints: a `&` of `==`,
`isin`, and `has_any` terms on the counted columns and `dataset_id`, and of `is_primary_data`.
Another query raises `ValueError`.

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
its location in the store, and the check of its fields. The tests and mypy therefore need the
writer installed. The SDK package itself does not depend on it.

Run both from this directory. `pytest` and `mypy` read their settings from `pyproject.toml` here. The tests write
small stores to a temporary directory and need no network and no credentials. Lint and format with `pipx run ruff check --fix .` and `pipx run ruff format .` from the
repository root. CI runs the unpinned `ruff check` and `ruff format --check` on this directory. See `AGENTS.md` here for the anndata and zarr
behavior that shaped this code.
