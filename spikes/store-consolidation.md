# Spike: store consolidation

This file records the store consolidation spike on the branch `spike/store-compaction`: what it
built, what it measured, and what it found. The branch is a running-code reference and is not
intended to merge in its current form. The draft pull request gives the summary; this file holds
the details and the experience for a person or an agent who returns to this direction.

Consolidation writes a consolidated store: a lossless second copy of several small per-dataset
stores in one store. The per-dataset stores stay, unchanged, as the primary copy. The design
record's section on compaction and other materializations calls such a copy a materialization.
This consolidation of stores is not Zarr's consolidated metadata, which anndata writes into every
store, of either kind.

## 1. What the spike built

- `census2_dataset_writer.consolidation`: `consolidate(sources, dest)` writes one consolidated
  store from several per-dataset stores, with bounded memory. `verify(dest, sources)` compares the
  two copies value for value. `read_consolidation` reads the consolidation record. `store_format`
  defines the record (`ConsolidationRecord`, `ConsolidatedMember`, `validate_consolidation`) and
  `obs_digest`.
- `consolidate-nf/`: a Nextflow pipeline, patterned on `bulk-ingest-nf`. One task for each group
  of a plan, a HealthOmics profile, receipts, and a run report. `tools/consolidation-plan.py`
  writes the plan from a catalog, outside the pipeline. The pipeline runs in the runtime image of
  `bulk-ingest-nf`.
- `census2-sdk/tools/build_catalog.py --consolidated PREFIX`: the builder checks each consolidated
  store against the per-dataset stores and writes three columns into the row of each member.
- The SDK: `get_anndata` opens each store once, reads only the rows of the selected datasets, and
  reads the datasets at the same time. The result equals the result of the per-dataset stores,
  the normalized layer included. `open_census(..., consolidated=False)` reads the per-dataset
  stores.
- `census2-sdk/tools/bench`: named queries over one catalog (`matrix.py --catalog --queries`), a
  `filter` shape, `--per-dataset`, two store opens (`--reader sdk` and `--reader lean`), a split
  of the corpus shapes over processes (`--processes`), the CPU time of each query, and
  `fanout.py`, which counts the datasets and the stores that each query reads from a catalog and a
  plan, with no store open.

The branch holds `main` as of the merge of #38. The writer of that `main` applies the Census 1
cell filters, writes the `tissue_general` pair and the five count statistics into `obs`, marks
the dataset's assayed genes in the `measured` column of `var`, and checks the cell key on every
write. The SDK of that `main` skips the `obs` index, decodes the cell key of the kept cells only,
and calculates the normalized layer on read. All of these hold on the consolidated path as well.

## 2. Size distribution and strategy

The 103 datasets (the 100-dataset uniform sample and the 3 harness datasets) hold 13.3 million
cells in 55.4 GiB of stores. The LTS rows come from the receipts of an earlier writer run over the
1,081 human stores of Census 1 LTS `2025-11-08` (629 GiB), and give the shape at corpus scale.

| | stores | min | p25 | median | p75 | p90 | max |
|---|---|---|---|---|---|---|---|
| 103 sample, store MiB | 103 | 6 | 31 | 120 | 426 | 1,238 | 5,701 |
| 103 sample, cells | | 1,465 | 8,410 | 33,041 | 78,111 | 267,123 | 1,712,244 |
| LTS human, store MiB | 1,081 | | 55 | 180 | 488 | 1,422 | 12,965 |
| LTS human, cells | | | 11,177 | 34,833 | 108,440 | 304,989 | 4,062,980 |

The distribution is heavy-tailed. The stores below 256 MiB are 62% of the sample and 8.0% of its
bytes, and 58% of the LTS human stores and 8.4% of their bytes.

The pilot strategy: a store of less than 256 MiB is small; pack the small stores into groups of
at most 1 GiB, first-fit decreasing by store bytes (`--strategy size`); leave the other stores
alone; drop a group of one. A group must compose: one organism, one gene reference version, all
spatial or all non-spatial, one writer version, and one access grant. A consolidated store of 1 GiB
is near the 88th percentile of the per-dataset stores, so no consolidated store is among the
largest stores, and a query still has many stores to read at the same time. The planner also has a
`tissue` order, which packs the stores of one main tissue together; it is not measured again here.

The plan puts 64 small stores into 5 consolidated stores of 5, 7, 10, 17, and 25 members: 103 → 44
stores (2.3 times fewer). Four groups are full, at 1,020 to 1,024 MiB of members; the fifth holds
the rest, 480 MiB. After consolidation:

| | stores | min | p25 | median | p75 | p90 | max |
|---|---|---|---|---|---|---|---|
| store MiB | 44 | 258 | 402 | 586 | 1,112 | 4,128 | 5,701 |
| cells | | 29,343 | 67,805 | 154,684 | 304,652 | 815,910 | 1,712,244 |

The smallest store goes from 6 MiB to 258 MiB, and the median store from 120 MiB to 586 MiB.
Among the 39 stores that stay alone, the four full consolidated stores (993 to 1,010 MiB) are at
the 69th percentile by bytes and at the 67th to 82nd by cells; the fifth (444 MiB) is at the 36th
by bytes and the 64th by cells. A consolidated store thus stays inside the size range of the
per-dataset stores, near the upper third.

At corpus scale (LTS human, planned from the receipts, `size` packing), the same thresholds give
1,084 → 509 stores. Beyond that, the large stores dominate the count:

| small below | group at most | stores after (of 1,081) | groups | largest group |
|---|---|---|---|---|
| 64 MiB | 1 GiB | 784 | 9 | 1,024 MiB |
| 128 MiB | 1 GiB | 644 | 23 | 1,024 MiB |
| 256 MiB | 1 GiB | 506 | 53 | 1,024 MiB |
| 256 MiB | 2 GiB | 480 | 27 | 2,048 MiB |
| 512 MiB | 1 GiB | 382 | 121 | 1,024 MiB |
| 512 MiB | 2 GiB | 322 | 61 | 2,048 MiB |

A quarterly rollup, which consolidates the small stores of each Discover publication quarter
apart, gives 519 stores (66 groups over 22 quarters) at 256 MiB and 1 GiB, near the 506 of
size-only packing. A time boundary thus costs few stores, and it suits continuous ingestion; see
section 11.

## 3. The consolidated store and the format adjustments

The per-dataset store does not change. Consolidation reads each member store as it is. Two
changes would help later, and neither is made: a summary record that carries `obs_digest` and
`nnz` would let a catalog builder check a consolidated store from metadata alone; and the `obs`
layout of the design record (chunks of 8,192 rows in shards of 65,536) would let a per-dataset
store serve partial reads as a consolidated store does.

The consolidated store is an AnnData-Zarr store, and `anndata.read_zarr` opens it. These are the
adjustments, each with its reason:

1. **Members are contiguous and in source order.** A member is one row range, `row_start` to
   `row_stop`, in the row order of its per-dataset store. A reader of one member reads one range,
   and a row ordinal of a per-dataset store maps to `row_start + ordinal`.
2. **A consolidation record, and no summary record.** The root attributes hold
   `cellxgene-census2/consolidation`: `version`, `gene_reference_version`, and the members with
   `dataset_id`, `dataset_version_id`, the row range, `nnz`, `obs_digest`, and the member's `uns`.
   A reader that expects one dataset in a store finds no summary record and refuses the store, so
   the two kinds of store cannot be confused. `obs_digest` is the SHA-256 of the
   `observation_joinid` values in row order; with `nnz`, it identifies the per-dataset store that
   a member came from, at the cost of one column read.
3. **The assayed genes go to `varm["measured"]`.** A per-dataset `var` is the gene reference plus
   `measured`, the dataset's assayed genes. A shared `var` cannot hold a per-dataset column, so the
   consolidated `var` is the reference alone, and `varm["measured"]` is a boolean matrix of genes ×
   members, with column `k` for member `k` of the record (78,894 × 25 at most here, one chunk of
   47 to 177 KiB). It is the members' part of the per-release presence matrix that the design record's
   section on the gene reference plans, in the form that anndata already reads, and it stays out of
   the root attributes, which every open parses. `verify` and the catalog builder compare each
   column with the member's `measured`. The SDK result leaves `measured` out on both paths.
4. **The `obs` index is the cell key**, `<dataset_id>_<observation_joinid>`. The index of a
   per-dataset store is the author's cell name. In the sample, 535,874 author cell names occur in
   more than one store, and 25 of the 103 stores hold one. The SDK sets the cell key as
   `obs_names` in every result already.
5. **Categoricals hold the sorted union of the members' categories.** No category is unused in the
   store. A reader of one member removes the categories that its rows do not use, and gets the
   categories of the per-dataset store, so a result does not depend on which copy served it.
6. **Fixed chunks.** Each `obs` array and `X/indptr` has chunks of 8,192 rows in shards of 65,536
   rows, the layout that the design record gives for `obs`. `X/data` and `X/indices` have chunks of
   131,072 values in shards of 8,388,608 values. anndata's default gives a small per-dataset column
   one chunk, so a partial read with that layout would read the column whole.
7. **`uns`.** The store's `uns` holds the keys that all members share with one value (the
   organism, and the schema version in the sample). The record holds each member's `uns`.
8. **Whole-shard writes.** `consolidate` writes `X/data` and `X/indices` one whole shard at a time.
   A write to part of a shard reads, merges, and writes the shard again.
9. **The address comes from the members.** A consolidated store is `c-` and the first 16
   hexadecimal digits of the SHA-256 of its sorted member dataset versions. One set of members has
   one address, and a new member set is a new store; nothing appends to a consolidated store.

The catalog keeps one row for each dataset and adds three columns: `consolidated_store_uri`,
`consolidated_row_start`, and `consolidated_row_stop`. They are null for a dataset that no
consolidated store holds. `store_uri` stays the per-dataset store, and every other field still
comes from the per-dataset store, so exclusion stays per dataset. A catalog without consolidated
stores has the columns, with every value null.

Alternatives for `measured`: a bitmap for each member in the record (about 13 KB of base64 for each
member, in the root attributes that every open parses), or a `var` column for each member (one more
array in the metadata for each member). The `varm` matrix costs one small array in the metadata,
whatever the member count, and a reader reads it only when it asks for the assayed genes.

Member boundaries at shard boundaries are not tried. Zarr's regular chunk grid gives one chunk size
for each array, so an aligned member would need padding rows, which are not cells, or a chunk grid
with chunks of different sizes, which zarr-python does not yet write for these arrays.

## 4. The consolidation run

The plan ran on HealthOmics with `--check`, one task for each group, 2 vCPU and a 4 GiB request.

| groups | members | cells | MiB in | MiB out | consolidate s (each) | check s (each) | peak RSS MiB | run |
|---|---|---|---|---|---|---|---|---|
| 5 | 64 | 1,350,867 | 4,571 | 4,457 | 103 to 148 | 97 to 133 | 779 to 978 | 11.6 min |

The run time includes the start of the HealthOmics run and its tasks, about 5 minutes. Every check
passed. A consolidated store is 1 to 7% smaller than its members together, and the copy adds 7.9%
to the bytes of the sample. The ingestion of the 103 datasets, with its check, took 46.9 minutes
and 398 task-minutes, and every check passed; consolidation of the small stores, with its check,
took 20 task-minutes, about 5% of that. The peak memory followed the member count more than the
cells, because each member reads its own `var` for the check of the members.

The catalog builder took 46 seconds for the 103 per-dataset stores on the in-region instance, with
16 workers, and 48 seconds with the 5 consolidated stores. The consolidated checks read two `obs`
columns, `varm["measured"]`, and some values of `X/indptr` of each consolidated store.

## 5. The SDK read path

`get_anndata` reads a row range of each store. A per-dataset store is the range of all its rows,
so one read path serves both kinds of store, and the result of the consolidated path equals the
result of the per-dataset path by construction as well as by test. The path:

1. Exclude per dataset, as before. Check the composition, as before.
2. Map each selected dataset to its read: the consolidated store and the member's rows, or the
   per-dataset store.
3. Open each store once, up to `max_workers` at the same time. The open reads the consolidated
   metadata and the `var` index; for a consolidated store, also the record and the categories of
   each categorical `obs` column.
4. Read each dataset, up to `max_workers` at the same time. `obs` reads up to 8 columns at the
   same time, without the index and without `observation_joinid`: each column whole for a
   per-dataset store, as `read_elem` reads it, and the member's rows for a consolidated store. The
   filter then gives the kept rows, and the read takes `observation_joinid` for the kept rows only.
   `X` reads `indptr` from the first to the last kept row, then the `data` and `indices` of the
   kept rows. `layers=["normalized"]` calculates the layer from these counts, with
   `feature_length` from the store's `var` for a full-gene assay. The datasets of one consolidated
   store read at the same time.
5. Concatenate in catalog order, as before, and take `var` from the first store, without
   `measured`.

`anndata.io.sparse_dataset` reads the whole `indptr`, which in a consolidated store holds every
member's rows, so the SDK reads `X` itself. When the kept rows follow one another, it reads one
slice; otherwise one integer selection of every kept stored value, as anndata's row read makes. The
open pull request #40 reads each run of adjacent kept rows as one range of `data` and `indices`,
and keeps anndata's read above 256 runs. The two overlap: a member read that keeps all its rows is
one run, and #40's read, given row numbers from `row_start` and an `indptr` read over the member's
rows only, would replace the slice case here and add the case of a few runs. The branch does not
depend on #40.

anndata 0.13 has no public partial read of a dataframe, so `_obs_rows` reads each column's row
range with the type that `read_elem` gives it. `anndata.experimental.read_elem_lazy` reads a
dataframe lazily, but needs xarray.

The branch reader also changes the per-dataset path: it reads the `obs` columns at the same time,
and a single-dataset query reads one slice of `X` where anndata built an integer selection of every
stored value. The measures therefore have three configurations, so that the reader and
consolidation show apart: `main` (the SDK of `main`), `per_dataset` (this branch, consolidation
off), and `consolidated`.

The tests check the equality of `obs` (with types and categories), `var`, `X`, and
`layers["normalized"]` between the two paths, with 1 and 4 workers, over two members that assayed
different genes, one of them from a full-gene assay. On the real stores, digests of each result
show the same equality between `main`, `per_dataset`, and `consolidated`; see section 6.4.

## 6. Measures

The rig is `census2-sdk/tools/bench`: each run is a new process, and the record holds the elapsed
time of the query, the CPU time of the process during the query (user and system), the requests and
the payload at the zarr store interface, and the peak RSS. The numbers are medians of 3 runs, with
the configurations of one repeat run one after another. The machines are an in-region EC2 instance
(us-west-2, m6in.4xlarge, 16 vCPU, 64 GiB) and a remote laptop (Apple silicon, about 6 MB/s from the
bucket, with a large spread between runs). Both run anndata 0.13.4, zarr 3.4.2, and obstore 0.12.1.

The configurations:

- `main`: the SDK of `main` at the merge point, over the per-dataset catalog. This is the baseline.
- `per_dataset`: this branch, with `open_census(..., consolidated=False)`.
- `consolidated`: this branch, over the catalog of the plan (44 stores).

Two store opens apply to both layouts (`bench.run --reader`):

- `sdk`: the open of `catalog.open_store`, a new obstore client for each store and no Zarr format
  (4 requests).
- `lean`: `zarr_format=3` and one obstore client for the bucket (1 request).

The corpus shapes (`open`, `obs` of every store) also run over 8 processes with the `lean` open.
`open` and `obs` of every store do not call `get_anndata`, so `main` and `per_dataset` are the same
for them; they run as `per_dataset`.

The queries are in `census2-sdk/tools/results/consolidation-queries.json`:

| query | filter | datasets kept | with a matching cell | stores, per-dataset → consolidated | cells |
|---|---|---|---|---|---|
| one dataset, a member | `dataset_id` of an 89,849-cell member (160 MiB store) | 1 | 1 | 1 → 1 | 89,849 |
| one dataset, alone | `dataset_id` of a 59,738-cell dataset (314 MiB store) | 1 | 1 | 1 → 1 | 59,738 |
| `cell_type` | neutrophil, the cell type nearest to 15% of datasets (14 of 103) | 14 | 14 | 14 → 9 | 35,795 |
| lung macrophages | six lung tissue terms × six macrophage terms | 5 | 5 | 5 → 5 | 311,063 |
| catalog false positive | `assay` 10x 3' v3 and `tissue` mesenteric lymph node | 6 | 0 | 6 → 6 | 0 |
| `obs` of every store | the bench's `corpus_obs` | 103 | | 103 → 44 | 13.3 M |
| open every store | the bench's `open` | 103 | | 103 → 44 | |

The lung macrophage query has no `is_primary_data` term, because the SDK filter does not yet take
that column; its 311,063 cells include the non-primary ones. Of its 5 datasets, one small one is a
member.

### 6.1 Open every store, in region

Seconds (CPU seconds), per-dataset → consolidated:

| open | 1 thread | 8 threads | 32 threads |
|---|---|---|---|
| `sdk` | 6.95 → 2.97 (3.2 → 1.4) | 2.88 → 1.24 (3.3 → 1.4) | 3.05 → 1.32 (3.4 → 1.5) |
| `lean` | 4.29 → 1.94 (1.6 → 0.7) | 1.62 → 0.73 (1.6 → 0.7) | 1.69 → 0.72 (1.7 → 0.7) |
| `lean`, 8 processes | | 0.77 → 0.46 (2.6 → 1.6) | 0.43 → 0.30 (2.9 → 1.6) |

With processes, the threads are in all: 8 processes of 1 or 4 threads. Requests: 412 → 176 (`sdk`),
103 → 44 (`lean`). The payload is 11.3 → 4.9 MB, the consolidated metadata.

One open, one store after another, in region (`opencost.py`, medians over the stores):

| part | per-dataset store | consolidated store |
|---|---|---|
| root `zarr.json` with the consolidated metadata | 112.5 KB | 112.6 KB |
| new obstore client, CPU | 12.9 ms | 12.4 ms |
| parse of the consolidated metadata, CPU | 13.9 ms | 13.4 ms |
| `sdk` open, wall (CPU) | 68 ms (31 ms) | 63 ms (29 ms) |
| `lean` open, wall (CPU) | 41 ms (15 ms) | 38 ms (14 ms) |

### 6.2 `obs` of every store, in region

Seconds (CPU seconds), per-dataset → consolidated:

| open | 1 thread | 8 threads | 32 threads |
|---|---|---|---|
| `sdk` | 160.1 → 83.7 (23.9 → 17.8) | 27.2 → 19.0 (23.5 → 18.3) | 19.0 → 15.0 (23.4 → 19.0) |
| `lean` | 145.5 → 80.2 (21.7 → 19.1) | 27.0 → 18.7 (21.7 → 17.8) | 18.2 → 14.6 (21.7 → 18.3) |
| `lean`, 8 processes | | 21.8 → 16.7 (23.9 → 20.7) | 8.37 → 8.15 (26.9 → 24.8) |

Requests: 6,021 → 3,477 (`sdk`), 5,712 → 3,345 (`lean`). The payload is 452 → 450 MB: the same
values, in fewer objects. The shape reads each `obs` whole, the index included, with `read_elem`.

### 6.3 Selective queries, in region

8 threads, the `sdk` open. Seconds for `main` / `per_dataset` / `consolidated`; requests, CPU seconds,
and peak RSS for `per_dataset` → `consolidated` (and `main` for RSS):

| query | stores | seconds | requests | CPU s | MB read | peak RSS GB |
|---|---|---|---|---|---|---|
| one dataset, a member | 1 → 1 | 5.77 / 2.23 / 2.40 | 111 → 218 | 3.4 → 3.5 | 160 | 4.3 / 1.6 / 1.6 |
| one dataset, alone | 1 → 1 | 8.37 / 3.18 / 3.20 | 164 → 164 | 6.7 → 6.7 | 314 | 6.9 / 3.0 / 3.0 |
| `cell_type` neutrophil | 14 → 9 | 10.19 / 8.58 / 8.20 | 2,969 → 3,270 | 17.1 → 16.8 | 913 | 1.6 / 1.5 / 1.6 |
| lung macrophages | 5 → 5 | 35.93 / 32.84 / 32.62 | 2,239 → 2,282 | 97.9 → 98.7 | 5,212 | 24.1 / 20.0 / 19.3 |
| catalog false positive | 6 → 6 | 2.32 / 1.11 / 1.39 | 352 → 484 | 1.5 → 1.7 | 47 | 0.6 / 0.6 / 0.6 |

By threads, seconds, `per_dataset` → `consolidated`, with the `sdk` open and then the `lean` open:

| query | 1 | 8 | 32 |
|---|---|---|---|
| one dataset, a member | 2.42 → 2.80; 2.51 → 2.94 | 2.23 → 2.40; 2.13 → 2.35 | 2.17 → 2.35; 2.17 → 2.62 |
| `cell_type` neutrophil | 23.8 → 26.7; 22.1 → 22.1 | 8.58 → 8.20; 8.25 → 8.21 | 7.31 → 7.36; 7.10 → 7.14 |
| catalog false positive | 2.86 → 3.25; 2.55 → 2.95 | 1.11 → 1.39; 1.06 → 1.27 | 1.10 → 1.33; 1.11 → 1.32 |
| lung macrophages | 43.6 → 44.2 | 32.8 → 32.6 | 32.4 → 32.0 |

Before the SDK read the categories of a consolidated store at the same time, the member read took
3.35 s and the false positive 2.28 s at 8 threads (the first in-region batch of the results file).

For neutrophils, 4 large stores that stay alone hold 94% of the 35,795 cells; 10 small datasets in 5
consolidated stores hold the other 6%. Consolidation opens 9 stores in place of 14, but the read of
`X` (913 MB) sets the time, and each member read starts and stops inside shards, which costs about
10% more requests.

### 6.4 The same result

On the real stores, a digest of each result (the `obs_names`, each `obs` and `var` column with its
dtype and categories, `X`, and `layers["normalized"]`) is the same from `main`, `per_dataset`, and
`consolidated`, for the selective queries and for the members of each consolidated store read whole
(`spikes/store-consolidation/digest.py`, `consolidation-digests.json`).
The cases are the five selective queries (the false positive returns no result on every
path) and the members of each of the 5 consolidated stores, read whole: 1,350,867 cells. Every
case is the same.

### 6.5 Laptop

A subset, with the `sdk` open unless stated. Seconds, per-dataset → consolidated:

| query | threads | seconds | CPU seconds | requests |
|---|---|---|---|---|
| open every store | 1 | 68.1 → 28.5 | 4.7 → 2.2 | 412 → 176 |
| open every store, `lean` | 1 | 19.9 → 7.6 | 3.0 → 1.0 | 103 → 44 |
| open every store | 8 | 9.4 → 4.4 | 2.7 → 1.2 | 412 → 176 |
| open every store, `lean` | 8 | 3.4 → 1.6 | 1.3 → 0.6 | 103 → 44 |
| open every store | 32 | 3.1 → 1.8 | 2.0 → 0.9 | 412 → 176 |
| `obs` of every store | 8 | 125.7 → 98.4 | 35.8 → 24.4 | 6,021 → 3,477 |
| one dataset, a member | 8 | 26.1 → 22.3 | 4.9 → 5.4 | 111 → 218 |

The runs of one cell spread by up to 1.9 times. On the laptop, the CPU time is a small part of the
wall time: the link sets the time.

### 6.6 Corpus scale, planned only

`bench.fanout` over a catalog of the 1,084 human LTS stores of an earlier writer and plans of the
same rule: datasets kept, then stores opened per-dataset → `size` → `tissue`. No consolidated store
of the LTS is written, and nothing here is timed.

| query | datasets kept | with a matching cell | stores |
|---|---|---|---|
| `cell_type` vascular associated smooth muscle cell (15%) | 160 | 160 | 160 → 135 → 128 |
| `cell_type` neutrophil | 133 | 133 | 133 → 94 → 89 |
| `tissue` lung | 73 | 73 | 73 → 59 → 48 |
| lung macrophages | 34 | 34 | 34 → 34 → 31 |
| 10x 3' v3 in mesenteric lymph node | 27 | 3 | 27 → 26 → 26 |
| female in thalamic complex | 30 | 0 | 30 → 29 → 28 |
| every store | 1,084 | | 1,084 → 509 → 515 |

An extrapolation, not a measure: at the CPU costs of section 6.1, opening all 1,084 human stores on
one core takes about 34 s with the `sdk` open and 16 s with the `lean` open; opening 509 stores takes
about 16 s and 8 s.

## 7. Intrinsic and reader-dependent gains

A gain is intrinsic when it comes from the layout: fewer stores to open, fewer metadata documents to
parse, fewer objects to fetch, fewer round trips in series. A gain is reader-dependent when a better
reader of the per-dataset stores gets it too, or when it hides a cost of the reader that a better
reader removes. The two opens and the process split of section 6 separate them.

**The open of a store.** One open costs about 31 ms of CPU with the `sdk` open: 13 ms for a new
obstore client, 14 ms to parse 112 KB of consolidated metadata, and some for the 4 requests. The
`lean` open removes the client and 3 requests, and costs about 15 ms; the parse stays. A store of
either kind has the same metadata, so the parse is per store, and consolidation divides it by the
store ratio. In region, at 8 and 32 threads, the wall time of `open` is about its CPU time: one core
under the GIL is the limit. Thus:

- With the `sdk` open, consolidation opens the corpus 2.3 times faster (2.88 → 1.24 s at 8
  threads). With the `lean` open, it is still 2.2 to 2.3 times faster (1.62 → 0.73 s). The ratio
  survives the leaner reader, because the part that the leaner reader leaves, the parse, is per
  store. The saved time halves, from 1.6 s to 0.9 s.
- The `lean` open alone gives the per-dataset stores 1.8 times (2.88 → 1.62 s at 8 threads), less
  than consolidation alone (2.3 times). The two multiply: 2.88 → 0.73 s, 3.9 times.
- 8 processes spread the parse over cores. Then the per-dataset stores open in 0.43 s (8 processes
  of 4 threads), faster than the consolidated stores with any threaded open, and consolidation adds 1.4 times
  (0.43 → 0.30 s). The CPU in all stays per store (2.9 → 1.6 s), so consolidation still saves about
  40% of the CPU.
- On the laptop, at 1 thread, the time follows the round trips in series, not the CPU or the request
  count: about 0.66 s for each `sdk` open (a new client makes new connections) and 0.19 s for each
  `lean` open (one round trip on a warm connection, and the transfer). Consolidation gives 2.4 times
  with the `sdk` open (68.1 → 28.5 s) and 2.6 times with the `lean` open (19.9 → 7.6 s): each open
  has at least one round trip, and fewer stores is fewer round trips. The `lean` open alone gives
  3.4 times. Both together give 9 times.

**`obs` of every store.** Requests fall from 6,021 to 3,477, and the payload stays (452 MB): the
same values in fewer objects. That part is intrinsic. Its effect depends on what limits the read:

- At 1 thread, the round trips in series limit it, and consolidation gives 1.9 times (160 → 84 s)
  with the `sdk` open and 1.8 times with the `lean` open.
- At 8 to 32 threads, one core limits it. Consolidation lowers the CPU by 1.2 to 1.3 times (fewer
  column reads, each with a fixed cost in zarr and anndata), and the time by 1.3 to 1.4 times
  (27.2 → 19.0 s at 8 threads). The `lean` open changes neither, because the open is a small part.
- Over 8 processes at 32 workers in all, the gain is gone: 8.37 → 8.15 s. The largest stores then
  set the time, and they stay alone.

**Selective queries.** They open few stores and read many bytes of `X`, so neither the open nor
consolidation counts for much. The neutrophil query opens 9 stores in place of 14 and takes about the
same time (8.58 → 8.20 s at 8 threads), with either open. A member read costs about twice the
requests of a whole-store read, because each partial read of a sharded array first reads the shard
index: 111 → 218 requests and 2.23 → 2.40 s for one member. That cost is intrinsic to partial reads
of a sharded layout; in region it is small. The catalog false positive, which reads only `obs`,
shows it most: 1.11 → 1.39 s. The `lean` open does not change these.

In short: the gain on opens is intrinsic, about the store ratio, with any thread reader, and smaller
but real with processes. The gain on `obs` of every store is intrinsic in requests, but it shows in
time only while round trips or one core limit the read. A stronger reader gives as much as
consolidation does, and more with processes, and the two compose.

## 8. Findings

1. **Fewer stores open faster, with any reader.** 2.3 times fewer stores (103 → 44) open 2.3 times
   faster, with 2.3 times less CPU, in region at 1 to 32 threads: 2.88 → 1.24 s with the SDK's open
   and 1.62 → 0.73 s with the `lean` open, at 8 threads. The parse of each store's 112 KB of
   consolidated metadata (about 14 ms of CPU) is per store, and one core under the GIL limits a
   threaded open. At corpus scale the plan gives 1,084 → 509 stores.
2. **A leaner open gives about as much as consolidation, and the two compose.** `zarr_format=3` and
   one shared obstore client make the per-dataset stores open 1.8 times faster in region (2.88 →
   1.62 s) and 3.4 times faster from the laptop at 1 thread (68.1 → 19.9 s). 8 processes make them
   open in 0.43 s, 6.7 times faster than the SDK's open on 8 threads, and faster than consolidated
   stores with any threaded open. Together with consolidation: 3.9 times in region with threads, 9 times from
   the laptop.
3. **`obs` of every store gains while round trips or one core limit it.** Requests fall 42% (6,021 →
   3,477) for the same 452 MB. The time falls 1.9 times at 1 thread, 1.4 times at 8 threads
   (27.2 → 19.0 s), and 1.3 times at 32. Over 8 processes it does not fall (8.37 → 8.15 s): the
   largest stores, which stay alone, set the time.
4. **Selective queries do not gain.** The neutrophil query opens 9 stores in place of 14 for the
   same time; 4 large stores that stay alone hold 94% of its cells. A member read is 8% slower
   (2.23 → 2.40 s) and a catalog false positive 25% slower (1.11 → 1.39 s), because partial reads of
   sharded arrays read shard indexes. Member boundaries at shard boundaries are not yet tried.
5. **The reader matters more than the layout for selective queries.** From `main` to this branch's
   per-dataset path, one dataset reads in 2.23 s in place of 5.77 s, with 1.6 GB in place of 4.3 GB
   at peak, because the branch reads one slice of `X` where anndata's row read builds an index of
   every stored value; and a false positive takes 1.11 s in place of 2.32 s, because the branch
   reads the `obs` columns at the same time. These are gains for per-dataset stores; #40 covers
   part of the `X` read.
6. **The laptop is latency-bound, the instance CPU-bound.** On the laptop at 1 thread, the CPU is 7%
   of the wall time of `open`, and each open costs its round trips (0.66 s with the SDK's open, 0.19
   s with a warm shared client). In region at 8 threads and more, the wall time of `open` and of
   `obs` of every store is about their CPU time. Consolidation helps both regimes for the open, and
   for `obs` of every store only while one of them holds.
7. **Consolidated stores stay inside the size range.** The four full ones (about 1 GiB) are at the
   69th percentile by bytes and at the 67th to 82nd by cells of the 39 stores that stay alone. The
   smallest store goes from 6 MiB to 258 MiB, and the median from 120 MiB to 586 MiB.
8. **Consolidation is cheap to run and to check.** 5 groups took 103 to 148 s each to write and 97
   to 133 s to check, on 2 vCPU, with less than 1 GiB of memory: 20 task-minutes, about 5% of the
   398 task-minutes of the ingestion. The copy adds 7.9% to the bytes. The catalog builder's checks
   add 2 s to a 46 s build.
9. **The result is the same.** On the real stores, the results of `main`, `per_dataset`, and
   `consolidated` have the same digests for the selective queries and for every member, the
   normalized layer included. The assayed genes move from `var` to `varm["measured"]`, and the
   result leaves them out on every path.
10. **The sample is small for this question.** 103 stores are well below the thousands of the
    concern. The corpus-scale store counts are planned, not written, and the corpus open times are
    extrapolated from the measured cost of one open.

## 9. Deviations from the design record

1. **The tentpoles.** The repository tentpoles say one store per dataset, independent writes, and
   one store per access boundary. A consolidated store holds several datasets, and a consolidation
   run reads several stores. On this branch the tentpoles are overridden as a secondary
   optimization. The per-dataset stores stay the primary copy, every per-dataset store is still
   written alone, and the SDK keeps the per-dataset path, as the design record's section on
   compaction and other materializations requires.
2. **Small groups, not one store for each organism.** The design record's example of such a copy
   is a per-organism store for corpus-wide queries. The pilot consolidates groups of small stores
   of at most 1 GiB. Reason: one store for each organism is a colossal store, outside the envelope
   in which AnnData-Zarr is proven. Groups of 1 GiB stay inside the size distribution of the
   per-dataset stores, and keep enough stores for a client to read at the same time.
3. **An explicit address in the catalog.** The design record has the SDK compute a public store's
   URI from `dataset_version_id` and `writer_version`. A consolidated store serves many datasets,
   so its address comes from its member set, and the catalog row holds it explicitly
   (`consolidated_store_uri`) with the row range. The prototype catalog already holds `store_uri`
   for every row.
4. **A consolidation record beside the summary record.** The design record defines one summary
   record for each store. A consolidated store holds a consolidation record instead, under the same
   namespace, and no summary record. Reason: a reader that expects one dataset must refuse a store
   of many, and a separate record does that with no change to the per-dataset record or its
   version.
5. **The assayed genes in `varm`.** The design record keeps each store's assayed gene set in the
   store and builds the presence matrix for each release. A consolidated store holds the members'
   sets as one genes × members matrix in `varm["measured"]`, because its `var` is shared.
6. **The `obs` layout of the design record, in consolidated stores only.** The design record gives
   `obs` chunks of 8,192 rows in shards of 65,536 rows. The consolidated store uses that layout; the
   per-dataset writer does not yet.
7. **The consolidation runs now.** The design record defers this work. The spike builds it to
   measure it; nothing here makes it a decision.

## 10. Relation to the `obs` index spikes

The other spikes (`spike/obs-iceberg`, #33, and the DuckDB spike) answer the `obs` filter from a
release-wide table of `obs` rows, each with its dataset and its store row ordinal, and then read
only the stores with a hit, by ordinal.

The two directions compose without a rebuild of either:

- A member of a consolidated store holds the rows of its per-dataset store in order, so a row of
  the consolidated store is `consolidated_row_start + ordinal`. The `obs` table needs no change
  when a dataset is consolidated, and a catalog row tells the SDK where the ordinal points.
- An exact `obs` index removes the catalog false positives and the `obs` reads of the stores,
  which is where consolidation gains most today (section 7). What remains for the stores is the
  open and the read of `X`. Consolidation then lowers the count of stores that the hits open; for
  selective queries, that count changes little (section 6).
- A materialized consolidation of `obs` would be one Parquet fragment for each consolidated
  store, with its members in the same row order, written by the same task.
- Both need the same access rule: a fragment and a consolidated store each have one grant.

## 11. Open questions

- **Continuous ingestion: append or rebuild.** A consolidated store is immutable and addressed by
  its member set, so nothing appends to it. New small stores wait in their per-dataset form until
  enough of them fill a group, and then a new consolidated store holds them. A later run writes new
  stores and changes no earlier one, so the catalog of an earlier release keeps its stores. A third
  level (consolidating consolidated stores) is possible, but the store count falls slowly once the
  large stores dominate it (section 2). A quarterly rollup is the same scheme with the quarter as
  the group boundary.
- **Re-curation of a member.** A re-curated dataset gets a new `dataset_version_id` and a new
  per-dataset store. Its catalog row then names no consolidated store, and the old consolidated
  store still serves its other members by row range; the old member's rows are not read. A later
  run can write a new group when the share of rows that no catalog row reads passes a threshold,
  and retention can remove the old store when no supported release names it. The builder of this
  spike refuses a member with no matching per-dataset store; a production builder would instead
  keep the other members.
- **Release manifests and retention.** A release manifest that pins a consolidated store must pin
  it as it pins a per-dataset store, with a digest. At the pilot thresholds the copy adds about 8%
  of the corpus bytes.
- **Access domains.** A group must share one grant. A dataset whose grant changes (an embargo that
  ends) joins a later group. Private data stays in per-dataset stores.
- **Medium stores.** Consolidating stores of up to 512 MiB would leave 382 of the 1,081 human
  stores instead of 506. Section 7 suggests that a leaner open buys most of what this would.
- **The layout of a consolidated store.** A member read pays a shard index for each partial read
  (section 6). Larger `obs` shards, unsharded chunks of 65,536 rows, or shard boundaries at the
  boundaries of large members are not yet measured. The chunk of `X/data` (131,072 values) is a
  guess.

## 12. Reproduce

The commands use placeholders for the bucket, the account, and the image.

```bash
# 1. Ingest the 103 datasets with bulk-ingest-nf on HealthOmics (see its README), to
#    s3://my-bucket/census2/stores/, with --check. The manifest is the two files under
#    bulk-ingest-nf/manifests/, one after the other.

# 2. Build the per-dataset catalog and plan.
python census2-sdk/tools/build_catalog.py s3://my-bucket/census2/stores/ catalog-pd.parquet
consolidate-nf/tools/consolidation-plan.py catalog-pd.parquet --small-mib 256 --target-mib 1024 \
    --strategy size --sizes sizes.json --out plan-size.jsonl
aws s3 cp plan-size.jsonl s3://my-bucket/census2/consolidated/plan-size.jsonl

# 3. Consolidate on HealthOmics. The image is the bulk-ingest-nf runtime image.
consolidate-nf/tools/build-healthomics-bundle.sh consolidate-bundle.zip
aws omics create-workflow --engine NEXTFLOW --main main.nf \
    --definition-zip fileb://consolidate-bundle.zip \
    --parameter-template file://consolidate-nf/tools/healthomics-parameters.json \
    --storage-type DYNAMIC --name census2-consolidate
aws omics start-run --workflow-id <id> --role-arn <run role> --storage-type DYNAMIC \
    --output-uri s3://my-bucket/census2/consolidated/omics-run-output \
    --parameters '{"plan": "s3://my-bucket/census2/consolidated/plan-size.jsonl",
                   "outdir": "s3_//my-bucket/census2/consolidated/size",
                   "container": "<account>.dkr.ecr.<region>.amazonaws.com/<repo>:<tag>",
                   "check": "true"}'

# 4. Build the consolidated catalog, with its checks.
python census2-sdk/tools/build_catalog.py s3://my-bucket/census2/stores/ catalog-size.parquet \
    --consolidated s3://my-bucket/census2/consolidated/size/

# 5. Measure, from census2-sdk/tools. main is a checkout of main at the merge point.
Q=results/consolidation-queries.json
python -m bench.fanout catalog-pd.parquet --queries $Q --plan size=plan-size.jsonl
PYTHONPATH=. <venv of main>/bin/python -m bench.matrix --catalog catalog-pd.parquet \
    --queries $Q --workers 1,8,32 --repeats 3 --config main --where ec2 --out runs.jsonl
for reader in sdk lean; do
  python -m bench.matrix --catalog catalog-size.parquet --queries $Q --per-dataset \
      --reader $reader --workers 1,8,32 --repeats 3 --config per_dataset --where ec2 --out runs.jsonl
  python -m bench.matrix --catalog catalog-size.parquet --queries $Q \
      --reader $reader --workers 1,8,32 --repeats 3 --config consolidated --where ec2 --out runs.jsonl
done
python -m bench.matrix --catalog catalog-size.parquet --queries $Q --only open,corpus_obs \
    --reader lean --processes 8 --workers 8,32 --repeats 3 --config consolidated --where ec2 \
    --out runs.jsonl   # and the same with --per-dataset
```

## 13. Sharp edges

- `catalog.open_store` makes a new obstore client for each store, and opens without a Zarr format.
  The client costs about 13 ms of CPU, and the open makes 4 requests, 3 of them for Zarr v2 keys that
  do not exist. A new client also opens new connections, which on a remote link costs round trips.
  `bench.run --reader lean` shows the open without both.
- At 8 and more threads, the wall time of `open` and of `obs` of every store is about equal to the
  CPU time of the process: one core under the GIL is the limit, and more threads do not help. More
  processes do.
- pandas 3 infers the `str` type for strings, so a column of `observation_joinid` written from a
  pandas Series is a `nullable-string-array` in the store (values and mask), not a `string-array`.
  A partial reader must handle both, and build the same pandas type as `read_elem`.
- A frame built from a dict of Series aligns them on their index. A Series with a `RangeIndex` in a
  frame with a cell-name index becomes all nulls, without an error. `consolidation.union_obs` builds
  from arrays for that reason, and `consolidate` refuses a null `observation_joinid`.
- anndata's row selection on a backed CSR matrix builds an integer array of every selected stored
  value. Reading one slice when the rows follow one another avoids that array.
- A 1-D zarr array without explicit chunks gets one chunk for a small array. A partial read of it
  reads the whole column. The consolidated store sets its chunks.
- The open of a consolidated store reads the categories of each categorical `obs` column. Read one
  after another, they cost about 20 round trips in series; the SDK reads them at the same time.
- Nextflow's launcher computes a cache key with `md5`, which on macOS is in `/sbin`. Without `/sbin`
  on `PATH`, the launcher fails with `exec: Usage:: not found`.
- On GitHub, a rename of the head branch of an open pull request closes the pull request. This
  branch keeps its first name for that reason.
