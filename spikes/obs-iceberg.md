# Spike: an Iceberg `obs` table over per-dataset fragments

**Status: a spike.** The branch `spike/obs-iceberg` is not intended to merge in its current form.
It is a running-code reference for one idea: publish a release-wide `obs` table beside the
catalog, so that one query selects the cells of the whole corpus, and the SDK opens only the
stores that hold a selected cell. The companion spike `spike/store-compaction` explores the other
idea of the pair. This file records what the branch does, what we measured, and what we learned,
for the people and agents who return to the idea.

## 1. The idea in one paragraph

The writer puts one Parquet fragment next to each store, with one row for each cell: the store's
`obs` columns, the store row ordinal, and the range of the cell's values in `X`. A builder
registers the fragments of a catalog as the data files of one Iceberg table for each organism,
in place, and writes only the Iceberg metadata under the Census prefix. A query first lets the
catalog exclude datasets, then scans the table once for the remaining datasets. The scan gives
the cell key `(dataset_id, observation_joinid)`, the store row ordinals, the result `obs`, and the
`X` ranges. The SDK then opens only the stores with a selected cell, reads the `X` ranges of the
selected cells, and reads no store `obs` and no `X/indptr`.

## 2. What the branch holds

| Part | Where | What it does |
|---|---|---|
| Fragment format | `census2_dataset_writer.store_format` | `fragment_location`, `FRAGMENT_ROW_COLUMNS`, `FragmentRecord`, `validate_fragment` |
| Fragment writer | `census2_dataset_writer.write`, `fragment_table`, `write_fragment` | writes `<name>.obs.parquet` after the store, then the fragment record beside the summary record |
| Pipeline | `bulk-ingest-nf` | publishes the fragment in local and remote modes; `--check` compares it with the store's `obs` and `X/indptr` |
| Catalog builder | `census2-sdk/tools/build_catalog.py` | copies the fragment's size and digest into the catalog row |
| Table builder | `census2-sdk/tools/build_obs_table.py` | one Iceberg commit per catalog: `add_files` for new rows, a partition delete for removed rows |
| SDK | `census2_sdk.obs_table`, `census2_sdk.expressions`, `Census.get_anndata`, `get_obs`, `count_obs` | the query path through the table |
| Filter language | `census2-sdk/tools/harness/filter_language.py`, `tools/results/filter_language.json` | the 74 bound corpus queries in the Iceberg expression language |
| Benchmark | `census2-sdk/tools/bench/` (`run`, `matrix`, `fragments`, `false_positives`, `readers`) | the measures below |

### 2.1 The fragment

```text
<prefix>/<dataset_version_id>.zarr/          the store, unchanged: plain AnnData-Zarr
<prefix>/<dataset_version_id>.obs.parquet    the obs fragment
```

- The columns are the store's `obs` columns, in the store's order, then `store_row`, `x_start`,
  and `x_end` (int64). A categorical column is plain strings with Parquet dictionary encoding,
  because Iceberg has no categorical type. The fragment takes whatever columns the store's `obs`
  holds, so `tissue_general` and the count statistics of the writer-obs-derived work arrive in
  the fragment with no change here.
- zstd, row groups of 65,536 rows (the `obs` shard of the store layout), delta-encoded row
  columns. 103 datasets, 13,326,922 cells: 187.5 MiB of fragments, 14.8 bytes per cell, 0.33
  percent of the 55.3 GiB of stores.
- The fragment record (`sha256`, `bytes`, `rows`) is in the store's root attributes, beside the
  summary record. The store is immutable, so the record proves which fragment belongs to it.
  The catalog row copies `bytes` and `sha256`; the table builder checks the size of each new
  fragment (one HEAD) and, with `--verify-digest`, the digest.

### 2.2 The table

```text
<census prefix>/obs/<organism>/metadata/
    00000-<uuid>.metadata.json, 00001-<uuid>.metadata.json, ...   one file for each commit
    snap-<id>-<uuid>.avro                                         manifest lists
    <uuid>-m0.avro                                                manifests
    version-hint.text                                             the current metadata file
```

- The data files are the fragments at their own locations. Nothing is copied.
- Partition spec: `identity(dataset_id)`. `add_files` derives each file's partition value from
  its column bounds, because each fragment holds one dataset.
- Metrics: the Iceberg default `truncate(16)` for strings, `counts` for `observation_joinid`.
- The snapshot summary holds `census2.catalog-sha256`, the digest of the catalog file that the
  table was built from. `open_census` refuses a table of another catalog, because a catalog row
  without a fragment in the table would silently lose its cells.

**No Iceberg catalog service is necessary.** The builder opens pyiceberg's `SqlCatalog` over a
SQLite file in a temporary directory, registers the current metadata file when the table exists
(`register_table`), commits, and then writes `metadata/version-hint.text`. Readers open a
`StaticTable` from a metadata file, or from the table location through the version hint, which
is the Hadoop-catalog convention that pyiceberg, DuckDB, and Polars read. A release pins the
metadata file itself, which is immutable. The cost: the version hint is a plain object write,
so two builders that commit at the same time can lose a commit. The catalog builder's one-writer
rule covers this, as it covers `catalog/current.parquet`; an S3 conditional write (`If-Match`)
on the hint would make the commit point a compare-and-swap.

### 2.3 A query

```python
census = census2_sdk.open_census(["catalog.parquet"], obs_tables=["s3://my-bucket/census/obs/homo_sapiens"])
adata = census.get_anndata(
    "homo_sapiens",
    "tissue IN ('lung') AND cell_type IN ('macrophage', 'alveolar macrophage') AND is_primary_data = true",
    max_workers=8,
)
```

1. `expressions.from_filter` gives an unbound Iceberg expression: from a pyiceberg expression,
   a string that pyiceberg parses, or the SDK's mapping filter.
2. `expressions.can_match` pushes the expression down to each catalog row. A term can match a
   dataset when some value of the row's value set satisfies it; `AND` and `OR` combine the
   results, and `NOT` moves to the terms first. A term on an unrecorded column proves nothing.
   For a mapping filter this gives exactly the rows of `catalog.exclude` (tested).
3. `expressions.with_members` rewrites the membership terms of the multi-valued `disease`
   columns into `IN` / `NOT IN` over the joined values that the remaining catalog rows record.
4. One `table.scan` with `expression AND dataset_id IN (remaining)`. Iceberg prunes the
   manifest entries by partition, then by column bounds; pyarrow prunes row groups by their
   statistics and filters rows. The SDK plans once and gives the tasks to pyiceberg's
   `ArrowScan`; pyiceberg reads the files in its thread pool.
5. The SDK groups the rows by dataset. A thread pool (`max_workers`) opens each store with a
   selected cell, checks its `X` shape against the catalog's `n_obs` and its `var`, and reads the
   selected `X` values by their ranges: one zarr selection of the value positions, which zarr
   fetches by chunk with its own concurrency.
6. The result `obs` comes from the scan. The result equals the result through the stores, value
   for value (tested on synthetic stores; checked on the 103 stores for the benchmark queries).

`get_obs` returns the selected `obs` from the table alone, and `count_obs` counts on any column,
for example `donor_id`, which the catalog's `cell_counts` does not count.

## 3. Measures

MEASURES

## 4. Findings

FINDINGS

## 5. The filter language

FILTERS

## 6. Readers

READERS

## 7. Deviations from the design record

DEVIATIONS

## 8. Relation to the compaction spike

COMPACTION

## 9. Open questions and future work

OPEN

## 10. Reproduce

REPRODUCE
