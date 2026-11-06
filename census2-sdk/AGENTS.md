# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Prototype scaffolding

The `ponytail:` comments in the code mark the list columns.

- The catalog holds the distinct values of each indexed column as a list of strings. The bitmaps of
  design 2.2 replace the lists when the catalog size matters.
- The catalog holds no numeric ranges and the filter accepts no numeric column. The stores hold no
  `raw_sum`, `nnz`, or the other three columns of design 2.1.

## The catalog builder

`tools/build_catalog.py` imports `census2-dataset-writer`, which defines the summary record, reads
it, and checks its fields with `read_summary`. The builder thus shares one definition of the record
with the writer. The SDK package does not depend on the writer, so a user of the SDK does not
install it. The tests and mypy do need it, and CI installs it.

## The cell counts

`tools/build_catalog.py` writes the cell counts of each store into its catalog row, in the
`cell_counts` column, and `Census.count_cells` sums them. The catalog stays one file,
and the counts and `n_obs` of a store come from the same build. `catalog.COUNTED` names the counted
columns, and the builder and the SDK both read that tuple. A null is a group of its own, so the
entries of a store add up to its `n_obs`. A null matches no filter, as in `get_anndata`. A store
without a counted column holds a null in that column. For a store without `is_primary_data`, every
entry has `is_primary_data` false. The builder groups on categorical keys, because on object columns
a store of millions of cells takes much longer.

The builder writes the column with pyarrow and the declared type `build_catalog.COUNT_ENTRY`.
Type inference gives a null type to a counted column that no store holds. The counts of two
catalogs from different builds then differ in type. `catalog.load` reads the catalog without the
column, because pandas turns the column into Python dicts. For the 1,258 stores, that read takes
135 ms instead of 5 ms. `catalog.load_counts` flattens the column with the Arrow list and struct
accessors of pandas in about 20 ms. pyarrow ships no type information, so its imports carry
`type: ignore[import-untyped]`.

`open_census` reads the content of each catalog once and keeps it. The rows and the counts come
from that content, so they agree also when the catalog changes later. The live catalog of design
3.3 changes with each update of the catalog builder.

## The gene reference

The writer projects each store onto the pinned gene reference of its organism: `X` holds the
counts, and `var` is the reference in reference order. The store's summary record holds
`gene_reference_version`, and the catalog builder copies it into the catalog row. The SDK holds
only the version names, in `gene_reference.py`. These names must agree with the writer's pin in
`census2_dataset_writer.genes`. The SDK refuses a row that records another version or no version.
It also refuses a store whose `var` differs from the `var` of the first store that it opened. It
does not remap genes.

## Sharp edges: anndata and zarr

The SDK reads a row range of a store, because a compacted store holds several datasets and a query
reads only the rows of the selected ones. A per-dataset store is the range of all its rows, so one
read path serves both kinds of store. `anndata.io.sparse_dataset` over a CSR group reads the full
`indptr`, which in a compacted store holds every member's rows, so `_x_rows` reads `X` itself:
`indptr` from the first to the last selected row, then the `data` and `indices` of the selected
rows. When the selected rows follow one another, that is one slice; otherwise it is one integer
selection, as anndata does. `_obs_rows` reads a row range of each `obs` column, with the type that
`read_elem` gives the column, because anndata 0.13 has no public partial read of a dataframe.
`anndata.experimental.read_elem_lazy` reads a dataframe lazily, but needs xarray. A categorical
column keeps only the categories that its rows use: in a compacted store, that gives a member the
categories of its per-dataset store, and the result does not depend on where a dataset was read.
`read_lazy` is not used because it does not round-trip in anndata 0.13 (see the writer's
`AGENTS.md`).

`get_anndata` opens each store once, then reads each selected dataset. Both steps use the thread
pool of `max_workers`. A read does not wait for another task of the pool, so the two steps cannot
deadlock it. Each read of a dataset also reads its `obs` columns at the same time, with up to
`OBS_THREADS` threads of its own pool. Each column is a read of its own, and in a compacted store,
whose columns are sharded, a column costs a round trip for the shard index and one for the data.
Read one column after another, a member of 152,509 cells read in 4.6 s from a compacted store and
in 3.6 s from its per-dataset store, in region. Read at the same time, the two took 3.7 to 4.1 s
and 3.1 s in two runs of the matrix. On the laptop the member read is still about 1.4 times slower
from a compacted store, because each round trip costs more there.

Opening a store parses its consolidated metadata, about 100 KB of JSON for a store of either kind.
In region, opening 103 stores took 2.9 s with 8 workers and also with 32, and so did the laptop with
32 workers. That is consistent with a cost in the client's CPU of about 28 ms for each store, which
more threads do not divide. Fewer stores is the lever for that cost.

zarr has no S3 client. A remote URI goes through `zarr.storage.ObjectStore` over `obstore`, as in
the writer. `obstore` reads credentials from the environment and not from an AWS profile, so export
them first with `aws configure export-credentials --format env`. Set `AWS_REGION` when the bucket
region is not the default.

## Sharp edges: types

`AnnData.X` is typed as a wide union that includes `None`. The tests wrap it with a helper that
casts to `scipy.sparse.csr_matrix`. `AnnData.obs` is `DataFrame | Dataset2D`; narrow with `cast`
where a frame method is called.

`DataFrame.to_dict("records")` is typed `dict[Hashable, Any]`. Cast to `dict[str, Any]` at the call
site.
