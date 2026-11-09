# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Prototype scaffolding

The `ponytail:` comments in the code mark the list columns.

- The catalog holds the distinct values of each indexed column as a list of strings. The bitmaps of
  design 2.2 replace the lists when the catalog size matters.
- The catalog holds no numeric ranges. The mapping filter accepts no numeric column. The typed
  filter builder accepts the five count statistics (`raw_sum`, `nnz`, and the others), but the
  catalog test proves nothing for a numeric term, so every dataset that the other terms keep is
  read.

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

The SDK reads a row range of a store, because a consolidated store holds several datasets and a
query reads only the rows of the selected ones. A per-dataset store is the range of all its rows, so
one read path serves both kinds of store. For all the rows of a store, `_obs_rows` reads each `obs`
column whole with `read_elem`. For part of the rows, it reads the row range of each column, with the
type that `read_elem` gives the column, because anndata 0.13 has no public partial read of a
dataframe. `anndata.experimental.read_elem_lazy` reads a dataframe lazily, but needs xarray. A
categorical column keeps only the categories that its rows use: in a consolidated store, that gives
a member the categories of its per-dataset store, and the result does not depend on where a dataset
was read. `read_lazy` is not used because it does not round-trip in anndata 0.13 (see the writer's
`AGENTS.md`).

`anndata.io.sparse_dataset` over a CSR group reads the full `indptr`, which in a consolidated store
holds every member's rows, so `_rows` reads `X` itself: `indptr` from the first to the last
selected row, then the `data` and `indices` of the selected rows. When the selected rows follow one
another, that is one slice; otherwise it is one integer selection of every selected stored value, as
anndata's row read makes. A read of runs of kept rows as ranges, one range for each run, would
generalize the slice; a member read would then pass row numbers from `row_start`.

`get_anndata` opens each store once, then reads each selected dataset. Both steps use the thread
pool of `max_workers`. A read does not wait for another task of the pool, so the two steps cannot
deadlock it. Each read of a dataset also reads its `obs` columns at the same time, with up to
`OBS_THREADS` threads of its own pool. Each column is a read of its own, and in a consolidated
store, whose columns are sharded, a partial read of a column costs a round trip for the shard index
and one for the data.

Opening a store parses its Zarr consolidated metadata, about 100 KB of JSON for a store of either
kind. `catalog.open_store` makes a new obstore client for each store, and zarr reads the metadata of
both Zarr formats when the format is not given, so an open makes 4 requests. See
`spikes/store-consolidation.md` for what that costs, and what a leaner open changes.

zarr has no S3 client. A remote URI goes through `zarr.storage.ObjectStore` over `obstore`, as in
the writer. `obstore` reads credentials from the environment and not from an AWS profile, so export
them first with `aws configure export-credentials --format env`. Set `AWS_REGION` when the bucket
region is not the default.

## The typed filter builder

`census.homo_sapiens` and `census.mus_musculus` are `Organism` namespaces. The builder is in four
modules: `obs_filter` (the columns, the expression tree, the normal form, and the label tests),
`obs_columns` (generated), `obs_plan` (the resolution and the catalog test), and `obs_arrow` (the
read from Zarr, the bind to each store, and the evaluation). `spikes/consolidation-filters.md`
records the measures and the reasons.

- `obs_columns.py` is generated by `tools/generate_obs_columns.py` from the writer's `OBS_COLUMNS`
  and `catalog.MULTIVALUED`. Do not edit it; run the tool. A test compares the columns of the copy
  with the store format, not the text, so a change of ruff's format does not fail it.
- mypy cannot make `==` a type error: for a refused `==`, it falls back to `object.__eq__` and
  gives `bool`. Each refused operator therefore carries `@deprecated` (PEP 702) and returns the
  type `Refused`. The error code `deprecated` is on in `pyproject.toml`, and the tests prove each
  refusal with a `# type: ignore[deprecated]` that `warn_unused_ignores` would report as unused.
  `typing_extensions.deprecated` is imported for the type checker only; at run time a no-op
  stands in, so the SDK needs no new dependency. An overload pair for `==` needs
  `# type: ignore[override]` on the first `@overload` line and `[overload-overlap]` on its `def`.
- Every term of the tree is two-valued (a null gives False), so `normalize` applies De Morgan's
  laws, and `~` turns the catalog's `NONE` into `ALL`. A term is `ALL` only when every value of the
  row's value set passes it and the column holds no null. The catalog test takes the store format's
  rule that no indexed column holds a null (the writer refuses one); a store format that admits
  nulls needs a null record in the catalog.
- One function, `obs_filter.label_test`, tests a label for the catalog and for the store, with
  Python's `str` and `re`, so the two cannot disagree. A regular expression is Python's, not RE2.
- pyarrow 26 decodes a dictionary array to one string for each cell for `equal` and `is_in`, and
  has no dictionary kernel for `starts_with` and the other string functions. `obs_arrow.bind`
  therefore tests each category once and maps the result to the cells by their codes (a code of
  -1 maps to False), and puts the result in the Arrow table as a boolean column.
- `Table.filter(expression)` takes every column of the table; `obs_arrow.project` evaluates the
  expression with an Acero `project` node in the calling thread, which costs about the time of the
  compute functions themselves.
- The builder path opens a store with `catalog.open_store(uri, lean=True)`: `zarr_format=3`, and
  one obstore client for each bucket. `bench.run.open_counted` takes the same keyword.
- A query reads in two rounds: the filter's columns over the rows of the read, then, for the kept
  rows only, the other columns of the result and the cell key. A read that keeps no cell reads
  nothing more. A candidate whose residual filter is True reads no column in the first round, and
  `count` takes its `n_obs` from the catalog with no store open.
- Adjacent members of one consolidated store read their filter columns in one read of their rows
  together; each member then evaluates its own residual filter on its part.
- A label that no catalog row holds raises `FilterError` when every label of the term is unknown,
  and warns when a near label exists beside known labels. An ontology list with absent terms
  therefore passes in silence.

## Sharp edges: types

`AnnData.X` is typed as a wide union that includes `None`. The tests wrap it with a helper that
casts to `scipy.sparse.csr_matrix`. `AnnData.obs` is `DataFrame | Dataset2D`; narrow with `cast`
where a frame method is called.

`DataFrame.to_dict("records")` is typed `dict[Hashable, Any]`. Cast to `dict[str, Any]` at the call
site.
