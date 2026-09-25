# census2-sdk

Query Census 2 stores through a catalog and return one AnnData.

The SDK is the skeleton of plan Phase 6, tracker 6.1. It runs one query end to end over the
Iteration 1 store set. The catalog excludes the datasets that cannot match. The SDK opens the
remaining stores, reads the selected cells, and concatenates them by gene position. The result is
one AnnData. The filter grammar, the per-chunk index, and private sources are later work.

## Install

```bash
pip install ./census2-sdk          # local catalog over local stores
pip install './census2-sdk[s3]'    # also read stores and catalogs on S3
```

zarr ships no S3 client. The `s3` extra adds `obstore`, and the SDK uses it only for a remote URI.

## Build the catalog

The catalog is a Parquet table with one row for each store. The Census 2 team publishes it, and the
SDK only reads it. `tools/build_catalog.py` is the Iteration 1 stand-in for the catalog builder. It
reads a prefix that holds `<dataset_version_id>.zarr` stores and the receipts next to them.

```bash
eval "$(aws configure export-credentials --profile <profile> --format env)"
python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet
```

Each row holds the identity of the dataset, the organism, the spatial flag, and `n_obs`. It also
holds the distinct values of the 18 indexed columns of design 2.2. The Iteration 1 stores hold no summary
record, so the builder reads `obs` and `uns` from each store and `dataset_id` from the receipt.

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

- `X`, the counts, from `raw/X` when the store has a `raw` group and from `X` when it does not
- `var`, the pinned gene reference for the organism (GENCODE v48 for human, 78,894 rows)
- `obs`, the columns that all selected stores hold, plus `dataset_id` and `observation_joinid`
- `obs_names`, `<dataset_id>_<observation_joinid>`
- `uns["census2"]`, the sources, the count of catalog rows for the organism, and the count of stores opened

The SDK downloads the gene reference once from the `single-cell-curation` tag `v7.0.0` into
`~/.cache/census2-sdk/v7.0.0/` and checks its SHA-256 digest against the pin on every load. Pass `refdata=<directory>` to `open_census` to read the file from a directory
that you staged.

## Checks before concatenation

The SDK applies the checks of design 4.3 to the selected catalog rows and reports the conflicting
values when a check fails:

- one organism, given by the query
- every recorded `gene_reference_version` equals the version the SDK pins for the organism; a row
  that records none passes
- no spatial store together with a non-spatial store
- a `dataset_id` in one source only, checked when the sources open

## Compare with Census 1

`tools/compare_census1.py` runs one query through the SDK and through Census 1 release 2025-11-08.
It compares the cell sets, the gene sets, and the values on the intersection. The Iteration 1
stores come from re-curated dataset versions, so a cell-set difference is a finding to record.

```bash
uv run --python 3.12 --with cellxgene-census==1.18.0 --with './census2-sdk[s3]' \
    census2-sdk/tools/compare_census1.py --catalog catalog.parquet --dataset-id <dataset_id>
```

## Develop

```bash
cd census2-sdk
pip install -e '.[dev,s3]'
pytest
mypy
```

Run both from this directory. `pytest` and `mypy` read their settings from `pyproject.toml` here. The tests write
small stores to a temporary directory and need no network and no credentials. Lint and format with `pipx run ruff check --fix .` and `pipx run ruff format .` from the
repository root. CI runs the unpinned `ruff check` and `ruff format --check` on this directory. See `CLAUDE.md` here for the anndata and zarr
behavior that shaped this code.
