# census2-sdk

The SDK skeleton: catalog exclusion, a read of the remaining stores, concatenation by gene position,
one AnnData out. `README.md` covers installing and using it. The module docstring in
`src/census2_sdk/__init__.py` states the result contract.

This file records the sharp edges of the code, and the scaffolding that the writer replaces.

## Iteration 1 scaffolding

The `ponytail:` comments in the code mark the source `var`, the `dataset_id` stamp, and the list
columns.

- The stores carry the source `var`. The SDK projects each store's genes onto the pinned reference
  before it concatenates (disposition rows 12, 14, 16, 20, 24). When the writer projects
  `var`, every store's `feature_id` set and order equal the reference, and the remap is an identity.
- The stores hold no `dataset_id` column, so the catalog builder reads it from the receipt and the
  SDK stamps it into `obs` (disposition row 13).
- The stores record no `gene_reference_version`, so the check of design 4.3 passes on the rows that
  record none.
- The catalog holds the distinct values of each indexed column as a list of strings. The bitmaps of
  design 2.2 replace the lists when the catalog size matters.
- The catalog holds no numeric ranges and the filter accepts no numeric column. The stores hold no
  `raw_sum`, `nnz`, or the other three columns of design 2.1.

## Sharp edges: anndata and zarr

anndata writes an absent `raw` as an array with `encoding-type: null`, not as a missing key. Test
the encoding, not the key, to decide whether `raw/X` holds the counts.

`anndata.io.sparse_dataset` over a CSR group reads the full `indptr`. It then reads only the `data`
and `indices` of the selected rows, so a query does not load the whole `X`. `obs` and `var` are read
whole. `read_lazy` is not used because it does not round-trip in
anndata 0.13 (see the writer's `CLAUDE.md`).

zarr has no S3 client. A remote URI goes through `zarr.storage.ObjectStore` over `obstore`, as in
the writer. `obstore` reads credentials from the environment and not from an AWS profile, so export
them first with `aws configure export-credentials --format env`. Set `AWS_REGION` when the bucket
region is not the default.

The CSR column remap in `_align_genes` rewrites `indices` through the reference positions. It then
sorts the indices, because a store's gene order is not the reference order.

## Sharp edges: types

`AnnData.X` is typed as a wide union that includes `None`. The tests wrap it with a helper that
casts to `scipy.sparse.csr_matrix`. `AnnData.obs` is `DataFrame | Dataset2D`; narrow with `cast`
where a frame method is called.

`DataFrame.to_dict("records")` is typed `dict[Hashable, Any]`. Cast to `dict[str, Any]` at the call
site.
