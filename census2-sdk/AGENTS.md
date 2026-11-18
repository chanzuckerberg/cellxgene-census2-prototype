# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Prototype scaffolding

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

## The release builder

`tools/build_release.py` imports the SDK and `census2-dataset-writer`. The SDK defines the manifest
in `census2_sdk.release`. The writer defines the name of the `measured` column and the name of the
per-gene nnz file. The SDK and the writer do not import each other, and neither imports the release
builder.

The snapshot is the content that the builder read with `catalog.read`, and not a frame that pandas
writes again. Thus the snapshot equals the live catalog, and the presence files and the manifest
come from the same catalog version. The builder writes the manifest last, under another name, and
then renames it. Thus `manifest.json` is complete when it exists. A build that stops leaves a
directory without a manifest, which `open_census` does not open. The builder refuses a directory
that exists, so a release holds no file of a previous build.

`open_census` finds a manifest by the content of the source. A Parquet file starts with `PAR1`, and
other content must be a manifest. The name of the file has no effect.

The manifest pins the identity of each store. It records no `writer_version` and no content digest
of a store, because no store and no catalog row records either value. A store that changes at its
location after the build changes the result of a release query.

The builder makes one presence file for each organism, and one set of aggregate rows. The catalog
also holds the `spatial` column, and the builder does not use it. Thus the file and the rows of an
organism cover its spatial and its non-spatial stores.

The builder adds the nnz arrays of an organism by position and reads the `feature_id` values from
the first store. It compares the length of `var` between the stores and not the gene identifiers.
One gene reference version for each organism is the reason that the positions agree.

The builder writes to the local file system only. To publish a release to object storage, copy
`manifest.json` last.

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
