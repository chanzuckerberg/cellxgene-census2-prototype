# Spike: store compaction

This file records the store compaction spike on the branch `spike/store-compaction`: what it
built, what it measured, and what it found. The branch is a running-code reference and is not
intended to merge in its current form. The draft pull request gives the summary; this file holds
the details and the experience for a person or an agent who returns to this direction.

The operator's brief calls the idea "Zarr store consolidation". Zarr and anndata already use
"consolidated metadata" for another thing, and anndata writes it into every store. This spike
therefore uses the terms of the design record: compaction writes a compacted store, which is a
materialization, a second copy of data next to the per-dataset stores.

## 1. What the spike built

- `census2_dataset_writer.compaction`: `compact(sources, dest)` writes one compacted store from
  several per-dataset stores, with bounded memory. `verify(dest, sources)` compares the two copies
  value for value. `read_compaction` reads the compaction record. `store_format` defines the
  record (`CompactionRecord`, `CompactedMember`, `validate_compaction`) and `obs_digest`.
- `compact-nf/`: a Nextflow pipeline, patterned on `bulk-ingest-nf`. One task for each group of a
  plan, a HealthOmics profile, receipts, and a run report. `tools/compaction-plan.py` writes the
  plan from a catalog, outside the pipeline. The pipeline runs in the runtime image of
  `bulk-ingest-nf`.
- `census2-sdk/tools/build_catalog.py --compacted PREFIX`: the builder checks each compacted store
  against the per-dataset stores and writes three columns into the row of each member.
- The SDK: `get_anndata` opens each store once, reads only the rows of the selected datasets, and
  reads the datasets at the same time. The result equals the result of the per-dataset stores.
  `open_census(..., compacted=False)` reads the per-dataset stores.
- `census2-sdk/tools/bench`: named queries over one catalog (`matrix.py --catalog --queries`), a
  `filter` shape, `--per-dataset`, and `fanout.py`, which counts the datasets and the stores that
  each query reads from a catalog and a plan, with no store open.

## 2. Size distribution and strategy

The 103 datasets (the 100-dataset uniform sample and the 3 harness datasets) hold 13.3 million
cells in 55.3 GiB of stores. The full Census 1 LTS run holds 1,081 human stores in 629 GiB, of an
older writer; its receipts give the sizes at corpus scale.

| | stores | p25 | median | p75 | p90 | max |
|---|---|---|---|---|---|---|
| 103 sample, store MiB | 103 | 30 | 128 | 425 | 1,231 | 5,669 |
| 103 sample, cells | | 9,270 | 33,041 | 78,111 | 267,123 | 1,712,244 |
| LTS human, store MiB | 1,081 | 55 | 180 | 488 | 1,422 | 12,965 |
| LTS human, cells | | 11,177 | 34,833 | 108,440 | 304,989 | 4,062,980 |

The distribution is heavy-tailed. The stores below 256 MiB are 63% of the sample and 8.8% of its
bytes, and 58% of the LTS human stores and 8.4% of their bytes.

The pilot strategy: a store of less than 256 MiB is small; pack the small stores into groups of
at most 1 GiB; leave the other stores alone; drop a group of one. A group must compose: one
organism, one gene reference version, all spatial or all non-spatial, one writer version, and one
access grant. A compacted store of 1 GiB is near the 88th percentile of the per-dataset stores, so
no compacted store is among the largest stores, and a query still has many stores to read at the
same time.

Two packing orders vary the strategy:

- `size`: first-fit decreasing by store bytes, which gives the fewest groups.
- `tissue`: the stores in the order of their main tissue (the tissue with the most cells, from the
  catalog's cell counts), one group filled at a time. Stores of one tissue then share a compacted
  store, so a tissue query opens fewer stores. It packs a little less tightly.

| | stores before | small stores | compacted stores | stores after |
|---|---|---|---|---|
| 103 sample, `size` | 103 | 65 | 5 (5 to 33 members) | 43 |
| 103 sample, `tissue` | 103 | 64 | 5 (9 to 20 members) | 44 |
| LTS human and mouse, `size` (planned only) | 1,255 | 738 | 61 | 578 |
| LTS human, `size` (planned only) | 1,084 | | | 509 |
| LTS human, `tissue` (planned only) | 1,084 | | | 515 |

The store count at corpus scale falls by about half, and then the large stores dominate it. The
thresholds move it as follows (LTS human, `size`, from the receipts):

| small below | group at most | stores after (of 1,081) | groups | largest group |
|---|---|---|---|---|
| 64 MiB | 1 GiB | 784 | 9 | 1,024 MiB |
| 128 MiB | 1 GiB | 644 | 23 | 1,024 MiB |
| 256 MiB | 1 GiB | 506 | 53 | 1,024 MiB |
| 256 MiB | 2 GiB | 480 | 27 | 2,048 MiB |
| 512 MiB | 1 GiB | 382 | 121 | 1,024 MiB |
| 512 MiB | 2 GiB | 322 | 61 | 2,048 MiB |

A quarterly rollup, which compacts the small stores of each Discover publication quarter apart,
gives 519 stores (66 groups over 22 quarters) at 256 MiB and 1 GiB. That is near the 506 of
size-only packing. Its groups are less even: the 10th percentile of group size is 228 MiB, against
760 MiB for the smallest group of size-only packing. A time-based rollup thus costs few stores, and
it suits continuous ingestion; see section 9.

## 3. The compacted store and the format adjustments

The per-dataset store does not change. No adjustment of it was essential: compaction reads each
member store as it is, and the per-dataset writer needs nothing new to make compaction cheap. Two
changes would help later, and neither is made: a summary record that carries `obs_digest` and
`nnz` would let a catalog builder check a compacted store from metadata alone; and the `obs`
layout of the design record (chunks of 8,192 rows in shards of 65,536) would let a per-dataset
store serve partial reads as a compacted store does.

The compacted store is an AnnData-Zarr store, and `anndata.read_zarr` opens it. These are the
adjustments, each with its reason:

1. **Members are contiguous and in source order.** A member is one row range, `row_start` to
   `row_stop`, in the row order of its per-dataset store. A reader of one member reads one range,
   and a row ordinal of a per-dataset store maps to `row_start + ordinal`.
2. **A compaction record, and no summary record.** The root attributes hold
   `cellxgene-census2/compaction`: `version`, `gene_reference_version`, and the members with
   `dataset_id`, `dataset_version_id`, the row range, `nnz`, `obs_digest`, and the member's `uns`.
   A reader that expects one dataset in a store finds no summary record and refuses the store,
   so the two kinds of store cannot be confused. `obs_digest` is the SHA-256 of the
   `observation_joinid` values in row order; with `nnz`, it identifies the per-dataset store that a
   member came from, at the cost of one column read.
3. **The `obs` index is the cell key**, `<dataset_id>_<observation_joinid>`. The index of a
   per-dataset store is the author's cell name. In the sample, 535,874 author cell names occur in
   more than one store, and 25 of the 103 stores hold one. The SDK sets the cell key as
   `obs_names` in every result already.
4. **Categoricals hold the sorted union of the members' categories.** No category is unused in the
   store. A reader of one member removes the categories that its rows do not use, and gets the
   categories of the per-dataset store, so a result does not depend on which copy served it.
5. **Fixed chunks.** Each `obs` array and `X/indptr` has chunks of 8,192 rows in shards of 65,536
   rows, the layout that the design record gives for `obs`. `X/data` and `X/indices` have chunks of
   131,072 values in shards of 8,388,608 values. anndata's default gives a small per-dataset column
   one chunk, so a partial read of a compacted store with that layout would read the column whole.
6. **`uns`.** The store's `uns` holds the keys that all members share with one value (the organism,
   and the schema version in the sample). The record holds each member's `uns`, because a
   compacted store has no single-dataset home for it.
7. **Whole-shard writes.** `compact` writes `X/data` and `X/indices` one whole shard at a time. A
   write to part of a shard reads, merges, and writes the shard again.
8. **The address comes from the members.** A compacted store is `c-` and the first 16 hexadecimal
   digits of the SHA-256 of its sorted member dataset versions. One set of members has one address,
   and a new member set is a new store; nothing appends to a compacted store.

The catalog keeps one row for each dataset and adds three columns: `compacted_store_uri`,
`compacted_row_start`, and `compacted_row_stop`. They are null for a dataset that no compacted
store holds. `store_uri` stays the per-dataset store, and every other field still comes from the
per-dataset store, so exclusion stays per dataset. A catalog without compacted stores has the
columns, with every value null. A private catalog without the columns reads as before.

## 4. The compaction run

Both plans ran on HealthOmics with `--check`, one task for each group, 2 vCPU and a 4 GiB request.

| plan | groups | members | cells | MiB in | MiB out | compact s (each) | check s (each) | peak RSS MiB | run |
|---|---|---|---|---|---|---|---|---|---|
| `size` | 5 | 65 | 1,513,588 | 4,960 | 4,842 | 123 to 137 | 88 to 115 | 799 to 1,354 | 11.0 min |
| `tissue` | 5 | 64 | 1,491,620 | 4,871 | 4,753 | 115 to 138 | 90 to 106 | 843 to 1,136 | 10.5 min |

The run time includes the start of the HealthOmics run and its tasks, about 5 minutes. Every check
passed. A compacted store is 2 to 3% smaller than its members together. The ingestion of the 103
datasets, with its check, took 42.9 minutes and 394 task-minutes, so compaction of the small stores
costs about 5% of the task time of ingestion. The peak memory followed the member count more than
the cells: each member read its own copy of the 78,894-row `var`.

The catalog builder took 46 seconds for the 103 per-dataset stores on the in-region instance, with
16 workers, and 48 seconds with the 5 compacted stores. The compacted checks read two `obs`
columns and some values of `X/indptr` of each compacted store.

## 5. The SDK read path

`get_anndata` now reads a row range of each store. A per-dataset store is the range of all its
rows, so one read path serves both kinds of store, and the result of the compacted path equals the
result of the per-dataset path by construction as well as by test. The path:

1. Exclude per dataset, as before. Check the composition, as before.
2. Map each selected dataset to its read: the compacted store and the member's rows, or the
   per-dataset store.
3. Open each store once, up to `max_workers` at the same time. The open reads the consolidated
   metadata, the `var` index, and the categories of each categorical `obs` column.
4. Read each dataset, up to `max_workers` at the same time: the member's rows of each `obs`
   column, then the `indptr` from the first to the last selected row, then the `data` and
   `indices` of the selected rows. The datasets of one compacted store read at the same time.
5. Concatenate in catalog order, as before.

`anndata.io.sparse_dataset` reads the whole `indptr`, which in a compacted store holds every
member's rows, so the SDK reads `X` itself. When the selected rows follow one another, it reads one
slice; otherwise one integer selection, as anndata does. anndata 0.13 has no public partial read of
a dataframe, so `_obs_rows` reads each column's row range with the type that `read_elem` gives it.
`anndata.experimental.read_elem_lazy` reads a dataframe lazily, but needs xarray.

The new reader also changes the per-dataset path: a single-dataset query reads one slice of `X`
where anndata built an integer selection of every stored value. The measures therefore have three
configurations, so that the reader and compaction show apart: `main` (the SDK of `main`),
`per_dataset` (this branch, compaction off), and `compacted`.

The tests check equality of `obs` (with types and categories), `var`, and `X` between the two paths,
with 1 and 4 workers. On the real stores, every member of both plans read whole through the SDK
equals its per-dataset read, and so do the single, cell type, and compound queries.

## 6. Measures

MEASURES

## 7. Findings

FINDINGS

## 8. Deviations from the design record

DEVIATIONS

## 9. Relation to the obs Iceberg spike

ICEBERG

## 10. Open questions

OPEN

## 11. Reproduce

REPRODUCE

## 12. Sharp edges

SHARP
