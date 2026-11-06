# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Prototype scaffolding

The `ponytail:` comments in the code mark the list columns.

- The catalog holds the distinct values of each indexed column as a list of strings. The bitmaps of
  design 2.2 replace the lists when the catalog size matters.
- The catalog holds no numeric ranges and the mapping filter accepts no numeric column. The stores
  hold `raw_sum`, `nnz`, and the other three statistics of design 2.1 in `obs`. On the spike
  branch, only the obs table path takes a range on them.

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

`anndata.io.sparse_dataset` over a CSR group reads the full `indptr`. It then reads only the `data`
and `indices` of the selected rows, so a query does not load the whole `X`. `obs` and `var` are read
whole. `read_lazy` is not used because it does not round-trip in
anndata 0.13 (see the writer's `AGENTS.md`).

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

## The obs table (spike)

`census2_sdk.obs_table`, `census2_sdk.expressions`, and `tools/build_obs_table.py` are the obs
Iceberg spike. `spikes/obs-iceberg.md` at the repository root records the design, the measures,
and the open questions. pyiceberg is the optional `iceberg` extra, and the modules import it
lazily, so a census without a table does not need it.

## Sharp edges: pyiceberg and pyarrow

pyiceberg's terms (`In`, `EqualTo`, ...) take a column name and literals at run time, but their
pydantic signatures name other fields, so mypy refuses `In("tissue", [...])`. `expressions.is_in`,
`not_in`, and `not_null` state the run-time form once. `expr.term.name` needs an
`attr-defined` ignore for the same reason.

`DataScan.to_arrow()` plans the scan again. `obs_table.scan` plans once, times the plan, and gives
the tasks to `pyiceberg.io.pyarrow.ArrowScan`, which is not a documented interface. The columns
come back in the order of the table schema, not in the order of `selected_fields`.

`ArrowScan(dictionary_columns=...)` returns string columns dictionary-encoded, which pandas turns
into categoricals without one Python string for each cell. Two edges: pyarrow has no
`starts_with` kernel for a dictionary column, so a column of a `LIKE 'prefix%'` term must stay
plain (`expressions.prefix_columns`); and pyiceberg logs a warning for each file that it reads so,
which `obs_table` filters. Arrow's `sort_by` refuses a dictionary column, so `_scan` sorts with
numpy. `pd.concat` of categoricals with different categories gives plain strings, so `get_obs`
concatenates the Arrow tables first and unifies their dictionaries.

`StaticTable.from_metadata` needs the region in `s3.region`; it does not read `AWS_REGION`.
`obs_table._properties` passes it. A table location without a metadata file name opens through
`metadata/version-hint.text`, which pyiceberg reads as the Hadoop catalog convention does.

`add_files` derives the identity partition value of a file from its column bounds, so it needs
`dataset_id` to have one value in each fragment. It refuses a file that holds Parquet field IDs;
pyarrow writes none. It sets the table's name mapping on the first call. Its default metrics mode
is `truncate(16)`, so the bounds of a long label are its first 16 characters.

A delete of `dataset_id IN (...)` on a table partitioned by `dataset_id` drops whole files and
rewrites nothing. A predicate that does not follow the partition would make pyiceberg rewrite a
fragment into the table location, a copy that the Census does not want. `build_obs_table.py`
deletes by partition only.

## DuckDB over the fragments (spike)

`census2_sdk.obs_sql`, `census2_sdk.fragments`, `tools/bench/sql.py`, `tools/bench/variants.py`,
`tools/build_descendants.py`, and `tools/harness/filter_language_sql.py` are the second `obs`
spike, stacked on the Iceberg spike. `spikes/obs-duckdb.md` at the repository root records the
design, the measures, and the open questions. DuckDB is the optional `duckdb` extra, and
`census2_sdk.fragments` holds the fragment column names, so the SQL path does not import pyiceberg.

## Sharp edges: DuckDB

DuckDB 2.0 ships its own type information (`py.typed` and stubs), so its imports need no ignore.
A relation keeps no result: `fetchone()` is on a connection's result, and a relation has
`to_arrow_reader()` and `to_arrow_table()`. `json_serialize_sql` gives literals as
`{"kind": ..., "text": ...}` in 2.0; `obs_sql.parse` reads that form only.

An Arrow array that DuckDB exports keeps its whole result chunk alive. A batch whose label columns
the SDK dictionary-encodes would otherwise keep the plain strings of those columns in memory
through its other columns. `obs_sql._encoded` copies each kept column into pyarrow's memory
(`pa.concat_arrays`), which took the peak of a read of all `obs` from 10.7 to 5.9 GiB.

With `enable_external_access = false`, DuckDB's `httpfs` cannot read the CA bundle, so the TLS
connection to S3 fails. The bundle must be in `allowed_paths`; `obs_sql._CA_BUNDLES` lists the
usual places. `allowed_directories` allows writes as well as reads in those directories, so only
the filter language keeps a `COPY` out of a query. `lock_configuration` also locks
`enable_logging`, so the HTTP log is on or off for the life of a connection.

A table function cannot take a subquery, so `read_parquet((SELECT list(...) FROM catalog))` fails
to bind. The SDK writes the file list into the SQL; a SQL script can use `SET VARIABLE` and
`getvariable` instead. `read_parquet` without `union_by_name` takes the schema of the first file
and refuses a file without one of its columns.
