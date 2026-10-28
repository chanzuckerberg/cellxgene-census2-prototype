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

1. **The tentpoles.** The repository tentpoles say one store per dataset, independent writes, and
   one store per access boundary. A compacted store holds several datasets, and a compaction run
   reads several stores. The operator's brief overrides the tentpoles on this branch, as a
   secondary optimization. The per-dataset stores stay the primary copy, every per-dataset store
   is still written alone, and the SDK keeps the uncompacted path, as the design record's section
   on compaction and other materializations requires.
2. **Small groups, not one store for each organism.** The design record's example of compaction is
   a per-organism store for corpus-wide queries. The pilot compacts groups of small stores of at
   most 1 GiB. Reason: one store for each organism is the colossal store that the brief rules out,
   and it is outside the envelope in which AnnData-Zarr is proven. Groups of 1 GiB stay inside the
   size distribution of the per-dataset stores, and keep enough stores for a client to read at the
   same time.
3. **An explicit address in the catalog.** The design record has the SDK compute a public store's
   URI from `dataset_version_id` and `writer_version`, and keeps an explicit `store_uri` for private
   stores only. A compacted store serves many datasets, so its address comes from its member set,
   and the catalog row holds it explicitly (`compacted_store_uri`) with the row range. The
   prototype catalog already holds `store_uri` for every row.
4. **A compaction record beside the summary record.** The design record defines one summary record
   for each store. A compacted store holds a compaction record instead, under the same namespace,
   and no summary record. Reason: a reader that expects one dataset must refuse a store of many,
   and a separate record does that with no change to the per-dataset record or its version.
5. **The `obs` layout of the design record, in compacted stores only.** The design record gives
   `obs` chunks of 8,192 rows in shards of 65,536 rows. The compacted store uses that layout; the
   per-dataset writer does not yet. Reason: a partial read needs it, and a per-dataset store is read
   whole today.
6. **The compaction runs now.** The design record defers compaction. The spike builds it to measure
   it; nothing here makes it a decision.

## 9. Relation to the obs Iceberg spike

The other half of the super-spike (`spike/obs-iceberg`) writes a Parquet `obs` fragment beside each
per-dataset store, with a store row ordinal, and queries a release-wide Iceberg table over the
fragments for `(dataset_id, observation_joinid, ordinal)`. The SDK then opens only the stores with
a hit and reads the rows of `X` by ordinal.

The two compose without a rebuild of either:

- A member of a compacted store holds the rows of its per-dataset store in order, so a row of the
  compacted store is `compacted_row_start + ordinal`. The Iceberg table needs no change when a
  dataset is compacted, and a catalog row tells the SDK where the ordinal points.
- The Iceberg table answers the `obs` filter exactly, so it removes the catalog false positives and
  the `obs` reads of the stores. What remains for the stores is the read of `X`. Compaction then
  lowers the count of stores that a hit opens. The two optimizations cut different costs: Iceberg
  cuts the `obs` bytes and the stores without a hit; compaction cuts the stores that each hold a
  few hits.
- The materialized consolidation that the brief defers would be one Parquet fragment for each
  compacted store, with its members in the same row order. The fragment ordinal is then the
  compacted row, and the table holds fewer, larger files. A compaction plan could write both in
  one task.
- Both need the same access rule: a fragment and a compacted store each have one grant.

## 10. Open questions

- **Continuous ingestion: append or rebuild.** A compacted store is immutable and addressed by its
  member set, so nothing appends to it. New small stores wait in their per-dataset form until
  enough of them fill a group, and then a new compacted store holds them. That is a size-tiered,
  LSM-like scheme of two levels: the per-dataset stores, then the compacted stores. A third level
  (merging compacted stores into larger ones) is possible, but the store count falls slowly once
  the large stores dominate it (section 2). A quarterly rollup is the same scheme with the
  quarter as the group boundary.
- **Re-curation of a member.** A re-curated dataset gets a new `dataset_version_id` and a new
  per-dataset store. Its catalog row then names no compacted store, and the old compacted store
  still serves its other members by row range; the old member's rows are never read. A later
  compaction can rewrite a group when the share of rows no catalog row reads passes a threshold,
  as an LSM tree drops dead entries. The builder of this spike refuses a member with no matching
  per-dataset store; a production builder would instead keep the other members.
- **Release manifests and retention.** A release manifest that pins a compacted store must pin it
  as it pins a per-dataset store, with a digest, and retention must keep it while a supported
  release names it. A compacted store is a second copy: at the pilot thresholds it adds about 8% of
  the corpus bytes.
- **Access domains.** A group must share one grant. A dataset whose grant changes (an embargo that
  ends) joins a later group. Private data stays in per-dataset stores.
- **The pull requests in flight.** PR #31 adds a per-dataset assayed-gene marker to `var`. A shared
  `var` cannot hold a per-dataset annotation, and `compact` refuses members whose `var` differs, so
  that marker must move: into the compaction record for each member, or into the per-release
  presence matrix that the design record plans and that is not built yet. Its `obs` columns
  (`tissue_general` and the five count statistics) compact as any other column. PR #30 drops cells
  in the writer, so a member holds fewer cells; compaction does not see the difference.
- **Medium stores.** Compacting stores of up to 512 MiB would leave 382 of the 1,081 human stores
  instead of 506. Whether the opening cost of a store justifies that depends on the measures of
  section 6 at corpus scale.
- **The `X` layout of a compacted store.** The chunk of `X/data` (131,072 values) is a guess. A
  smaller chunk lowers the over-read at member boundaries and for scattered rows; a larger one
  lowers the requests.

## 11. Reproduce

The commands use placeholders for the bucket, the account, and the image.

```bash
# 1. Ingest the 103 datasets with bulk-ingest-nf on HealthOmics (see its README), to
#    s3://my-bucket/census2/stores/, with --check.

# 2. Build the per-dataset catalog and plan.
python census2-sdk/tools/build_catalog.py s3://my-bucket/census2/stores/ catalog-pd.parquet
compact-nf/tools/compaction-plan.py catalog-pd.parquet --small-mib 256 --target-mib 1024 \
    --strategy size --sizes sizes.json --out plan-size.jsonl
aws s3 cp plan-size.jsonl s3://my-bucket/census2/compacted/plan-size.jsonl

# 3. Compact on HealthOmics. The image is the bulk-ingest-nf runtime image.
compact-nf/tools/build-healthomics-bundle.sh compact-bundle.zip
aws omics create-workflow --engine NEXTFLOW --main main.nf \
    --definition-zip fileb://compact-bundle.zip \
    --parameter-template file://compact-nf/tools/healthomics-parameters.json \
    --storage-type DYNAMIC --name census2-compact
aws omics start-run --workflow-id <id> --role-arn <run role> --storage-type DYNAMIC \
    --output-uri s3://my-bucket/census2/compacted/omics-run-output \
    --parameters '{"plan": "s3://my-bucket/census2/compacted/plan-size.jsonl",
                   "outdir": "s3_//my-bucket/census2/compacted/size",
                   "container": "<account>.dkr.ecr.<region>.amazonaws.com/<repo>:<tag>",
                   "check": "true"}'

# 4. Build the compacted catalog, with its checks.
python census2-sdk/tools/build_catalog.py s3://my-bucket/census2/stores/ catalog-size.parquet \
    --compacted s3://my-bucket/census2/compacted/size/

# 5. Measure. queries.json is in the pull request description; main is a checkout of main.
cd census2-sdk/tools
python -m bench.fanout catalog-pd.parquet --queries queries.json --plan size=plan-size.jsonl
PYTHONPATH=. <venv of main>/bin/python -m bench.matrix --catalog catalog-pd.parquet \
    --queries queries.json --workers 1,8,32 --repeats 3 --config main --where ec2 --out runs.jsonl
python -m bench.matrix --catalog catalog-size.parquet --queries queries.json --per-dataset \
    --workers 1,8,32 --repeats 3 --config per_dataset --where ec2 --out runs.jsonl
python -m bench.matrix --catalog catalog-size.parquet --queries queries.json \
    --workers 1,8,32 --repeats 3 --config compacted --where ec2 --out runs.jsonl
```

## 12. Sharp edges

- pandas 3 infers the `str` type for strings, so a column of `observation_joinid` written from
  a pandas Series is a `nullable-string-array` in the store (values and mask), not a
  `string-array`. A partial reader must handle both, and build the same pandas type as
  `read_elem`; `_column_rows` follows anndata's reader for it.
- A frame built from a dict of Series aligns them on their index. A Series with a `RangeIndex`
  in a frame with a cell-name index becomes all nulls, without an error.
  `compaction.union_obs` builds from arrays for that reason, and `compact` refuses a null
  `observation_joinid`.
- anndata's row selection on a backed CSR matrix builds an integer array of every selected stored
  value. For a single dataset of 152,509 cells, the `main` SDK peaked at about 5.9 GB on the
  instance. Reading one slice when the rows follow one another avoids that array.
- The equality check of all 65 compacted members at once needs two results of 1.5 million cells
  in memory, which passed 64 GB on the instance; one compacted store at a time peaked at 15 GB.
  A result of millions of cells costs memory in any layout.
- A 1-D zarr array without explicit chunks gets one chunk for a small array. A partial read of it
  reads the whole column. The compacted store sets its chunks.
- Nextflow's launcher computes a cache key with `md5`, which on macOS is in `/sbin`. Without
  `/sbin` on `PATH`, the launcher fails with `exec: Usage:: not found`.
