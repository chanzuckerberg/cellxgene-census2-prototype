# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Prototype scaffolding

- The catalog holds the distinct values of each indexed column as a list of strings. The bitmaps of
  design 2.2 replace the lists when the catalog size matters.
- The catalog holds no numeric ranges, so a numeric term of a query excludes no dataset.

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

## The query builder

This branch is the spike `spike/obs-arrow-native`; `spikes/obs-arrow-native.md` records it.
`query.py` types the operators of each column and gives a `pyarrow.compute.Expression` for each.
pyarrow evaluates the expression (`obs_read.evaluate`), so the SDK owns no filter language.

- A pyarrow expression has no public API to walk it. Thus the builder records, as it builds, the
  columns of a query and its hints, and `Organism.expression` takes the list of columns from the
  caller. The SDK reads only the recorded columns; a column that the expression names and the list
  does not makes pyarrow raise.
- The hints come only from a top-level `&` of `==`, `isin`, and `has_any` terms. `|`, `~`, a
  numeric term, and a string method give none, so they exclude no dataset. `catalog.covers` proves
  a store whole only when the query is exactly its hints (`Query.exact`), and that proof relies on
  the store format's rule that an indexed column holds no null. A test store with a null shows
  the sharp edge (`test_null_is_not_a_value`).
- pyarrow's `is_in` gives false for a null, not null. The builder wraps it in
  `if_else(is_valid(x), ..., null)`, so that `isin` follows the three-valued logic of `==`.
- pyarrow 26 decodes a dictionary array to one string for each cell for `==` and `is_in`, and it
  has no string kernel for a dictionary array. The string methods and `has_any` therefore cast the
  column to strings first. `has_any` is an RE2 expression over the joined value, about 65 ns for
  each cell; see the spike notes for its share of a query.
- A categorical column of the catalog takes its categories from the catalog row when the store has
  as many categories as the row lists, by the rule of `ColType.CATEGORY`.
  `dataset_id` reads the categories of the store, because the check of the dataset compares the
  store with the row.
- `obs_columns.py` is generated from the writer's `store_format.OBS_COLUMNS` by
  `tools/generate_obs_columns.py`. `tests/test_query.py` fails when the two differ.
- The refused operators carry `@deprecated` and the return type `Refused`, and mypy runs with the
  error code `deprecated` on, so mypy reports `obs.disease == "x"` on its line.

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
and `indices` of the selected rows, so a query does not load the whole `X`. The SDK reads `obs` one
column at a time with zarr (`obs_read`), and `var` whole. `read_lazy` is not used because it does
not round-trip in anndata 0.13 (see the writer's `AGENTS.md`).

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
