# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Iteration 1 scaffolding

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

## The cell count table

`tools/build_catalog.py` writes the cell count table next to the catalog (tracker 5.9), and
`Census.count_cells` sums it. `catalog.COUNTED` names its columns, and the builder and the SDK both
read that tuple. A null is a group of its own in the table, so the rows of a store add up to its
`n_obs`. A null matches no filter, as in `get_anndata`. A store without a counted column holds a
null in that column. A store without `is_primary_data` counts as not primary. The builder groups
on categorical keys, because on object columns a store of millions of cells takes much longer.
Before it counts, the SDK compares the table total of each dataset with the catalog `n_obs`. A table
from another build gives a wrong count and no error, and this check finds it.

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
