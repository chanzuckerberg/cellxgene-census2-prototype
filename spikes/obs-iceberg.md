# Spike: an Iceberg `obs` table over per-dataset fragments

**Status: a spike.** The branch `spike/obs-iceberg` is not intended to merge in its current form.
It is a running-code reference for one idea: publish a release-wide `obs` table beside the
catalog, so that one query selects the cells of the whole corpus, and the SDK opens only the
stores that hold a selected cell. The companion spike `spike/store-compaction` (#34) explores the
other idea of the pair. This file records what the branch does, what we measured, and what we
learned, for the people and agents who return to the idea.

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

The branch is `main` (at e54290c) plus the spike. The store path of `get_anndata` is the code of
`main`: the cell filters, the `tissue_general` pair and the count statistics in `obs`, the
`measured` column in `var`, the normalized layer on read, and the `obs` read without the barcode
index and with the cell key of the kept cells only.

| Part | Where | What it does |
|---|---|---|
| Fragment format | `census2_dataset_writer.store_format` | `fragment_location`, `FRAGMENT_ROW_COLUMNS`, `FragmentRecord`, `validate_fragment` |
| Fragment writer | `census2_dataset_writer.write`, `fragment_table`, `write_fragment` | writes `<name>.obs.parquet` after the store, then the fragment record beside the summary record |
| Pipeline | `bulk-ingest-nf` | publishes the fragment; `--check` compares it with the store's `obs` (the statistics bit for bit) and `X/indptr` |
| Catalog builder | `census2-sdk/tools/build_catalog.py` | copies the fragment's size and digest into the catalog row |
| Table builder | `census2-sdk/tools/build_obs_table.py` | one Iceberg commit per catalog: `add_files` for new rows, a partition delete for removed rows |
| SDK | `census2_sdk.obs_table`, `census2_sdk.expressions`, `Census.get_anndata`, `get_obs`, `count_obs` | the query path through the table |
| Filter language | `census2-sdk/tools/harness/filter_language.py`, `tools/results/filter_language.json` | the 74 bound corpus queries in the Iceberg expression language |
| Benchmark | `census2-sdk/tools/bench/` (`run`, `matrix`, `counting`, `pruning`, `fragments`, `false_positives`, `readers`) | the measures below |

### 2.1 The fragment

```text
<prefix>/<dataset_version_id>.zarr/          the store, unchanged: plain AnnData-Zarr
<prefix>/<dataset_version_id>.obs.parquet    the obs fragment
```

- The columns are the store's `obs` columns, in the store's order, then `store_row`, `x_start`,
  and `x_end` (int64). The fragment takes the columns that the store's `obs` holds, so the
  `tissue_general` pair and the five count statistics of `main` arrived in the fragment with no
  change to the fragment code.
- zstd, row groups of 65,536 rows (the `obs` shard of the store layout), delta-encoded row
  columns, and a dictionary for the other columns.
- The fragment record (`sha256`, `bytes`, `rows`) is in the store's root attributes, beside the
  summary record. The store is immutable, so the record proves which fragment belongs to it.
  The catalog row copies `bytes` and `sha256`; the table builder checks the size of each new
  fragment (one HEAD) and, with `--verify-digest`, the digest.

**Types and bounds.** Each column keeps the type of its store column:

| Store column type | Columns | Parquet | Iceberg | Iceberg bounds |
|---|---|---|---|---|
| categorical, string | the labels, `tissue_general` and its term ID, `observation_joinid` | `BYTE_ARRAY` (UTF8) | `string` | the first 16 characters (`truncate(16)`); `observation_joinid` counts only |
| bool | `is_primary_data` | `BOOLEAN` | `boolean` | exact |
| float64 | `raw_sum`, `raw_mean_nnz`, `raw_variance_nnz` | `DOUBLE` | `double` | exact |
| int64 | `nnz`, `n_measured_vars`, and the row columns | `INT64` | `long` | exact |

No statistic is NaN, so no float bound is lost to NaN. Over the 103 files, the bounds of `nnz` run
from 1 to 38,085; the median file has `nnz` in [305, 8,795]. `n_measured_vars` is one value in each
dataset (from 458 to 60,606), so its bounds are exact for each file.

**Size.** 103 datasets, 13,316,710 cells: 476.9 MiB of fragments, 37.5 bytes per cell, 0.84
percent of the 55.4 GiB of stores. That is 2.5 times the 14.8 bytes per cell of round 1, before
the statistics. Two derived statistics are about half of the bytes. Per cell, in 92 of the
fragments (4.4 million cells): `observation_joinid` 10.6 B, `raw_variance_nnz` 9.4 B,
`raw_mean_nnz` 9.1 B, `raw_sum` 2.6 B, `nnz` 1.8 B, `x_start` and `x_end` 1.5 B each, and all 21
label columns together less than 1 B. A scan reads only the columns that it selects and filters,
so the statistics cost bytes only in a scan of all columns and in storage. `BYTE_STREAM_SPLIT` on
the three float columns makes those fragments about 14 percent smaller (32.4 B per cell); the
writer does not use it yet.

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
- Metrics: the Iceberg default, `truncate(16)`: the bounds of a string are its first 16
  characters, and the bounds of a number are exact.
- The snapshot summary holds `census2.catalog-sha256`, the digest of the catalog file that the
  table was built from. `open_census` refuses a table of another catalog, because a catalog row
  without a fragment in the table would silently lose its cells.

**No Iceberg catalog service is necessary.** The builder opens pyiceberg's `SqlCatalog` over a
SQLite file in a temporary directory, registers the current metadata file when the table exists
(`register_table`), commits, and then writes `metadata/version-hint.text`. Readers open a
`StaticTable` from a metadata file, or from the table location through the version hint, which
is the Hadoop-catalog convention that pyiceberg reads. A release pins the metadata file itself,
which is immutable. The cost: the version hint is a plain object write, so two builders that
commit at the same time can lose a commit. The catalog builder's one-writer rule covers this, as
it covers `catalog/current.parquet`; an S3 conditional write (`If-Match`) on the hint would make
the commit point a compare-and-swap.

### 2.3 A query

```python
census = census2_sdk.open_census(
    ["catalog.parquet"], obs_tables=["s3://my-bucket/census/obs/homo_sapiens"]
)
adata = census.get_anndata(
    "homo_sapiens",
    "tissue IN ('lung') AND cell_type IN ('macrophage', 'alveolar macrophage') AND nnz >= 500",
    layers=["normalized"],
    max_workers=8,
)
```

1. `expressions.from_filter` gives an unbound Iceberg expression: from a pyiceberg expression,
   a string that pyiceberg parses, or the SDK's mapping filter.
2. `expressions.can_match` pushes the expression down to each catalog row. A term can match a
   dataset when some value of the row's value set satisfies it; `AND` and `OR` combine the
   results, and `NOT` moves to the terms first. A term on a column that the catalog does not
   record, such as `nnz`, proves nothing. For a mapping filter this gives exactly the rows of
   `catalog.exclude` (tested).
3. `expressions.with_members` rewrites the membership terms of the multi-valued `disease`
   columns into `IN` / `NOT IN` over the joined values that the remaining catalog rows record.
4. One `table.scan` of the expression. Iceberg prunes the manifest entries by the column
   bounds, and the SDK then keeps the files of the remaining datasets by their partition value;
   pyarrow prunes row groups by their statistics and filters rows. The SDK plans once and gives
   the tasks to pyiceberg's `ArrowScan`, which reads the files in its thread pool. (A
   `dataset_id IN (remaining)` term in the scan gives the same files, but pyiceberg plans it in
   time that grows with files times values: 23 s for 2,060 files and 2,060 values.)
5. The SDK groups the rows by dataset. A thread pool (`max_workers`) opens each store with a
   selected cell, checks its `X` shape against the catalog's `n_obs` and the number of genes of
   the first store, and reads the selected `X` values by their ranges: one zarr selection of the
   value positions, which zarr fetches by chunk with its own concurrency. It does not read the
   store's `var` index (section 7, item 9). With `layers=["normalized"]` it calculates the layer
   from those counts, as `main` does.
6. The result `obs` comes from the scan. The result equals the result of `main`'s store path,
   value for value, the normalized layer included (tested on synthetic stores with a full-gene
   assay; checked on the 103 stores for the benchmark queries).

`get_obs` returns the selected `obs` from the table alone, and `count_obs` counts on any column,
for example `donor_id`, which the catalog's `cell_counts` does not count.
`get_obs(dtype_backend="pyarrow")` returns `pd.ArrowDtype` columns over the scan's Arrow data.

## 3. Measures

### 3.1 Method

**Data.** The 103 human datasets of the prototyping manifest and the harness manifest, ingested
with the merged writer on AWS HealthOmics (103 admitted, 103 round-trip checks passed, the fragment
check included, in 47 minutes). The cell filters of `main` keep 13,316,710 cells (they drop
10,212); 55.4 GiB of stores; 476.9 MiB of fragments.

**Machines.** In-region: one EC2 m6in.4xlarge in us-west-2 (16 vCPU, 64 GiB), with the stores
and the table in an S3 bucket of the same region. Laptop: an Apple-silicon laptop on a home link
to us-west-2, shared with the other spikes under a lock that gives one spike the laptop at a time.

**Runs.** `bench.matrix` runs each query in a new process; the values are medians of 3 runs.
"Stores" is the store path of `main` at the merge point (`get_anndata` through the stores, which
this branch does not change), and "table" is this branch's path. The catalog is a local file on
both paths; the table is in S3.

**CPU time.** Each record holds the user and system CPU time of the process (`getrusage`, all
threads) during the query, beside the elapsed time. CPU divided by elapsed time is the mean
number of busy cores.

**Requests.** Each record counts the requests for every object that the query reads, by kind:

- *Store objects*, at the zarr store interface (`bench.counting.CountingStore`). One group of byte
  ranges that zarr fetches together counts as one request.
- *The version hint*, at `census2_sdk.catalog.read`: one GET.
- *The Iceberg metadata and the fragments*, through a counting `pyarrow.fs.FileSystemHandler`
  under pyiceberg's `PyArrowFileIO` (`bench.counting.CountingHandler`). Each file open is one
  HEAD, because pyarrow's S3 file system asks for the size of each file that it opens, and each
  read of the open file is one GET, because pyarrow reads a range of an S3 object with one GET.

The handler puts Python in each read of a fragment, and the reads of one file take turns: at 8
workers it made (d1) 0.66 s to 1.13 s, (b) 5.57 s to 6.66 s, and (e) 5.08 s to 5.38 s
(`ec2-103-count-overhead.jsonl`). The times therefore come from runs without it, and the request
counts of the table path from one more run of each query with it (`--count-arrow`). The counts do
not depend on the timing: the Parquet reader plans its coalesced reads from the file metadata. The
store path needs no second run. "Requests" below are the requests of the query; opening the table
adds 3 more (the version hint, and a HEAD and a GET of the metadata file).

**Queries** (`tools/results/obs_iceberg/queries.json`). The compaction spike uses its own values
for (a), (c), and (d); (b) and (e) are the same in both.

| Name | Filter | Shape |
|---|---|---|
| (a) one dataset | `dataset_id` of the first catalog row | `get_anndata` |
| (b) neutrophil | `cell_type = neutrophil`, the cell type whose share of the datasets is nearest 15 percent (14 of 103) | `get_anndata` |
| (c) lung macrophages, primary | `tissue` in 6 lung terms, `cell_type` in 9 macrophage terms, `is_primary_data` | `get_anndata` |
| (d1) catalog false positive, none | `tissue = mesenteric lymph node AND cell_type = megakaryocyte`: 4 datasets pass the catalog, none holds a cell | `get_anndata` |
| (d2) catalog false positive, some | `tissue = liver AND cell_type = mesothelial cell`: 4 pass, 2 hold cells | `get_anndata` |
| (e) all `obs` | none | stores: read each store's `obs`; table: `get_obs` |
| (f) `tissue_general` | `tissue_general = pancreas` (3 datasets) | `get_anndata` |
| count by donor | the filter of (b) | cells by `donor_id`; stores: read each candidate store's `obs` |
| (g) WMG count | the harness `wmg_production` filter: `is_primary_data`, `nnz >= 500`, `cell_type_ontology_term_id != 'unknown'`, six `tissue_general_ontology_term_id !=` terms | cells by tissue × cell type; stores: read each candidate store's `obs` and apply the filter in pandas (`bench.run.evaluate`) |

The mapping filter of `main` cannot express (g), so the store path of (g) is what a user of
`main` would write: the catalog test of the filter, then each candidate store's `obs` in pandas.

### 3.2 In-region, 103 datasets

Seconds at 1 / 8 / 32 workers. Stores opened, CPU time, requests, and MB: store path / table
path. CPU time is at 8 workers; the requests and the MB are of the query, from the counted runs
(`ec2-103.jsonl`, `ec2-103-counted.jsonl`; the store path of (g) from `ec2-103-wmg.jsonl`).

| Query | Cells | Stores opened | Store path (`main`) | Table path | CPU s | Requests | MB read |
|---|---|---|---|---|---|---|---|
| (a) one dataset | 3,799 | 1 / 1 | 2.37 / 2.20 / 2.28 | 1.10 / 1.15 / 1.08 | 0.96 / 0.65 | 77 / 32 | 17.0 / 16.9 |
| (b) neutrophil | 35,795 | 14 / 14 | 33.8 / 9.78 / 8.01 | 9.53 / 5.57 / 5.54 | 17.7 / 16.4 | 3,038 / 1,383 | 913 / 929 |
| (c) lung macrophages | 207,185 | 5 / 5 | 40.3 / 25.9 / 26.0 | 27.1 / 21.9 / 22.1 | 82.9 / 80.5 | 2,084 / 1,254 | 4,131 / 4,134 |
| (d1) no matching cell | 0 | 4 / 0 | 7.01 / 2.13 / 2.09 | 0.94 / 0.66 / 0.61 | 1.34 / 1.50 | 283 / 19 | 59 / 107 |
| (d2) 2 of 4 hold cells | 111 | 4 / 2 | 9.64 / 4.44 / 4.67 | 2.23 / 2.00 / 1.95 | 2.45 / 2.20 | 908 / 231 | 90 / 121 |
| (f) pancreas | 12,815 | 3 / 3 | 7.88 / 4.45 / 4.39 | 3.48 / 2.80 / 2.75 | 5.03 / 4.28 | 511 / 150 | 122 / 124 |
| (e) all `obs` | 13,316,710 | 103 / 0 | 153 / 23.6 / 15.9 | 5.46 / 5.08 / 5.17 | 18.7 / 11.3 | 5,801 / 322 | 380 / 483 |
| count by donor, (b) | 35,795 | 14 / 0 | 19.0 / 4.36 / 3.48 | 0.56 / 0.59 / 0.56 | 3.46 / 0.63 | 837 / 131 | 84 / 2.8 |
| (g) WMG count | 8,365,927 | 57 / 0 | 83.2 / 15.4 / 11.3 | 2.96 / 2.76 / 2.66 | 14.0 / 2.87 | 3,311 / 788 | 304 / 24 |

- **Requests by kind.** The table path's requests are its store objects, the fragments (one HEAD
  and 2 to 13 GETs for each fragment that it scans: the reader coalesces only near column chunks,
  so a scan of a few columns over many row groups takes more GETs), and four Iceberg metadata
  reads (the manifest list and the one manifest, a HEAD and a GET each). For (b): 1,334 store,
  45 fragment, and 4 metadata requests. For (e): 103 HEADs and 215 GETs of fragments. Opening the
  table is 3 more requests on the table path.
- **The table phase** (plan and scan) took 0.2 to 0.7 s for (a) to (d) and (f); the plan read a
  manifest of 88 KB in 0.1 s.
- **CPU.** The store path of (e) uses 19 s of CPU at any number of workers. At 32 workers it takes
  15.9 s, so about 1.2 cores are busy: the read is bound by the GIL. The table path uses 11.3 s
  of CPU in 5.1 s, 2.2 cores, most of it in Arrow outside the GIL. For (b) and (c), the CPU of the
  two paths is within 10 percent: decoding `X` is most of it.
- **Peak memory** follows the result: 13.5 GiB (table) and 16.0 GiB (stores) for (c), about 2 GiB
  for (b). For (e), the table path holds a 13.3-million-row frame (7.4 GiB peak), and the store
  path reads each `obs` and drops it (1.1 GiB at 8 workers).

### 3.3 All `obs`: the scan and the frame

`get_obs` for (e) at 8 workers, median of 3 (`ec2-103.jsonl`, `ec2-103-arrow-obs.jsonl`):

| Result | Seconds | Table (plan, scan, group) | Frame | CPU s | Peak GiB |
|---|---|---|---|---|---|
| categoricals (`dtype_backend="numpy"`, the default) | 5.08 | 1.79 | 3.22 | 11.3 | 7.4 |
| Arrow-backed (`dtype_backend="pyarrow"`) | 2.32 | 1.67 | 0.41 | 8.6 | 5.2 |
| store path, for comparison | 23.6 | | | 18.7 | 1.1 |

The scan of all 483 MB of fragments is 1.4 s. The frame of categoricals is most of the rest: one
remap of the integer codes of each of the 21 label columns, with the labels that the result uses,
sorted, as `get_anndata` gives them. The Arrow-backed frame keeps the scan's dictionary columns,
and pandas copies nothing; only the cell key index is built (0.4 s, one Arrow join for each
dataset in a thread pool). The Arrow-backed option is the cheap remedy for a reader that does not
need numpy categoricals.

### 3.4 The store path with a stronger column read

The store path reads every `obs` column of each candidate store, as `get_anndata` does. A store
reader that reads only the columns of the filter and of the count (`--narrow`, still through
anndata's `read_elem`, one column at a time) changes the obs-only comparison
(`ec2-103-wmg.jsonl`, `ec2-103-narrow.jsonl`, `ec2-103-narrow-1.jsonl`). Seconds at
1 / 8 / 32 workers, and the requests:

| Query | Store path, every column | Store path, needed columns | Table path |
|---|---|---|---|
| count by donor, (b) | 19.0 / 4.36 / 3.48 (837) | 2.53 / 0.75 / 0.62 (121) | 0.56 / 0.59 / 0.56 (131) |
| (g) WMG count | 83.2 / 15.4 / 11.3 (3,311) | 14.5 / 3.47 / 3.23 (631) | 2.96 / 2.76 / 2.66 (788) |

With only the needed columns, the store path reads 3 MB instead of 84 MB for the count by donor
and 31 MB instead of 304 MB for (g), and at 8 and 32 workers it comes within 1.1 to 1.3 times of
the table. So most of the table's gain on these obs-only questions is the column selection,
which a Zarr store also allows. The table keeps a lead at low concurrency (4 to 5 times at 1
worker), where its few large requests count. A stronger Zarr reader would narrow the gaps
further: read only the needed columns, read `zarr_format=3` metadata, and share one obstore client
across stores. What a store read cannot avoid is the open and the read of each store that the
catalog keeps, and for the catalog's false positives that is the whole cost.

### 3.5 In-region, growth with N

Copies of the 103 fragments with new `dataset_id`s (`bench.fragments`): 500, 1,000, and 2,060
fragments, one table each. A copy's catalog row points at the source store, so a query through
the table reads the source store again for each copy; a query through the stores does the same,
with the `dataset_id` check relaxed for copies. Both paths therefore read the same objects many
times, which can favour both. These are measured, not extrapolated; the store count is still
below the corpus. Seconds at 32 workers, median of 3 (`ec2-scale.jsonl`); the store path ran at
103 and 2,060 only. Requests at 2,060 from `ec2-scale-counted.jsonl`.

| N | Manifest | Plan | (b) stores / table | (d1) stores / table | (d2) stores / table | count by donor, stores / table |
|---|---|---|---|---|---|---|
| 103 | 88 KB | 0.12 s | 8.0 / 5.5 | 2.1 / 0.61 | 4.7 / 2.0 | 3.5 / 0.56 |
| 500 | 408 KB | 0.15 s | – / 22.1 | – / 1.6 | – / 3.3 | – / 1.9 |
| 1,000 | 811 KB | 0.21 s | – / 45.5 | – / 2.7 | – / 5.7 | – / 3.5 |
| 2,060 | 1.66 MB | 0.29 s | 140 / 99 | 16.3 / 6.7 | 27.6 / 11.8 | 52.5 / 7.4 |
| 2,060, requests | | | 60,475 / 27,299 | 5,660 / 304 | 17,875 / 4,259 | 16,740 / 2,544 |

- The plan stays below 0.3 s. The time of the table path follows the stores that a query opens
  and the bytes that it scans. (b) is `X`-bound: 16.4 GB of `X` at 2,060 on both paths.
- (d1) at 2,060: the table path opens none of the 80 stores that the catalog keeps, but its scan
  reads 2.1 GB of 80 atlas fragments, against 1.2 GB of `obs` on the store path; both use about
  27 s of CPU. The key-first scan reads only the filter columns and uses a quarter of the CPU
  (7 s), but at 32 workers it took 7.9 s against 6.7 s (`ec2-scale-keyfirst.jsonl`).
- At 2,060 the files that the column bounds alone keep for (b) are 1,420; the catalog keeps 280.

### 3.6 Pruning by the column statistics

`bench.pruning` plans each filter, then asks pyarrow which row groups of the kept files can match
by their Parquet statistics (`pruning-103.jsonl`, `pruning-2060.jsonl`). The catalog keeps every
dataset for a numeric term, because it records no numeric range.

| Filter, 103 datasets | Cells that match | Catalog keeps | Files after the bounds | Row groups after the statistics |
|---|---|---|---|---|
| `nnz >= 500` | 12,511,453 | 103 | 101 | 268 of 268 |
| `nnz >= 5000` | 1,164,209 | 103 | 87 | 197 of 199 |
| `nnz < 200` | 113,645 | 103 | 36 | 82 of 151 |
| `raw_sum >= 100000` | 86,569 | 103 | 49 | 108 of 156 |
| `raw_sum < 300` | 18,044 | 103 | 21 | 56 of 61 |
| `n_measured_vars < 30000` | 4,658,136 | 103 | 31 | 94 of 94 |
| (g) WMG filter | 8,365,927 | 57 | 58 | 184 of 201 |
| `cell_type = 'neutrophil'` | 35,795 | 14 | 71 | 221 of 224 |
| `tissue_general = 'lung'` | 1,993,595 | 17 | 28 | 105 of 115 |

- **Numeric bounds prune files that the catalog cannot.** A selective range keeps about a fifth to
  a half of the files: `raw_sum >= 100000` (0.65 percent of the cells) keeps 49 of 103 files.
  The plan still takes about 0.1 s. A column with one value in each dataset, such as
  `n_measured_vars`, prunes files exactly.
- **Row groups add little.** The rows of a fragment are in source order, so the range of a
  statistic in a row group of 65,536 rows is most of the range of its file: `raw_sum >= 100000`
  keeps 108 of the 156 row groups of those 49 files. Fragments sorted by a common filter column,
  with the store row ordinal as the link to the store, would make row groups prune; the spike does
  not sort them.
- **Label bounds stay weak.** `cell_type = 'neutrophil'` keeps 71 files by its 16-character
  lexicographic bounds, and the catalog keeps 14; at 2,060 datasets, 1,420 against 280. The
  selection on labels is the catalog's work, applied through the `dataset_id` partition.
- **In a compound filter the two levels add up.** For (g) the catalog keeps 57 datasets; the
  bounds keep 58 files, and the intersection is 55 files: two of the catalog's datasets have no
  cell with `nnz >= 500`.
- At 2,060 datasets the shares are the same (`raw_sum >= 100000`: 980 of 2,060 files), and the
  plan takes 0.2 to 0.5 s.

### 3.7 Laptop, 103 datasets

A small subset at 8 workers, median of 3 (`laptop-103.jsonl`, `laptop-103-keyfirst.jsonl`). The
link read at about 5 to 10 MB/s, and repeats of one query varied by up to a factor of 2, so these
are indications. Seconds:

| Query | Stores | Table | Table, key first |
|---|---|---|---|
| (a) one dataset | 11.7 | 7.3 | |
| (d1) no matching cell | 18.2 | 32.7 | 3.7 |
| (d2) 2 of 4 hold cells | 44.6 | 21.7 | 18.9 |
| count by donor, (b) | 20.6 | 3.4 | |

On this link, bytes decide. A count reads two columns of 14 fragments (about 2 MB) where the
store path reads 84 MB of `obs`. But `get_anndata` scans every column of each candidate fragment,
and the catalog's false positives are large atlases: for (d1) the scan reads about 107 MB of four
fragments to find no cell, more than the 59 MB of `obs` that the store path reads. The
fragments of this round are 2.5 times larger than in round 1 (the statistics), so this case got
worse. The key-first scan reads the filter columns first and the other columns only of the
fragments with a match: (d1) took 3.7 s. For (d2), whose two matching datasets are atlases, it
gains less.

### 3.8 Other measures

- **False positives of the catalog** (`false_positives.json`). Over the 10,852 (tissue, cell type)
  pairs that some dataset passes in the catalog, the catalog passes 12,569 datasets, and 6,944
  (55 percent) hold no cell of the pair; 6,500 pairs hold no cell anywhere. Over the 4,352 pairs
  that some cell holds, 261 of 5,886 passing datasets (4.4 percent) hold none. The false
  positives come from the multi-tissue atlases, which are among the largest stores.
- **Correctness on the real stores** (`ec2-equal.jsonl`). For (a) to (d) and (f), with
  `layers=["normalized"]`, the table path and the store path give equal `obs_names`, `obs` values
  and dtype kinds, `var`, `X`, and normalized layer.
- **Ingest, catalog, and table costs.** The HealthOmics run of the 103 datasets with `--check`
  took 47 minutes from start to end (the statistics read `X` once more). The catalog of 103
  stores took 127 s on EC2. The table took 33 s with `--verify-digest`, which reads and hashes all
  477 MB of fragments (29 s), and the commit 2.2 s. The table of 2,060 copies took 171 s: a HEAD
  of each new fragment, and a commit of 51 s, in which `add_files` reads each footer.
- **The key-first scan in region** (`ec2-103-keyfirst.jsonl`, `ec2-scale-keyfirst.jsonl`). At 103
  datasets and 8 workers: (d1) 0.66 to 0.51 s, (d2) 2.00 to 2.08 s, (b) 5.57 to 5.90 s. At 2,060
  and 32 workers: (d1) 6.7 to 7.9 s, (d2) 11.8 to 13.4 s. It uses less CPU and reads fewer bytes
  when the candidates hold no match, and costs a second plan and scan when they do.

## 4. Findings

1. **The table makes the catalog's false positives cheap.** The catalog tests each term alone, so
   it keeps a dataset that has some cell of the tissue and some cell of the type but no cell of
   both. The table path opens none of those stores: (d1) 2.09 s to 0.61 s at 32 workers, and 283
   to 19 requests; at 2,060 datasets, 16.3 s to 6.7 s, with 80 stores not opened and 5,660 to
   304 requests. Over all (tissue, cell type) pairs that the catalog passes, 55 percent of the
   passes hold no matching cell (4.4 percent over the pairs that some cell holds). The false
   positives are mostly multi-tissue atlases, which are large stores.
2. **But the scan reads every column of each candidate fragment, so on a slow link bytes
   decide.** With the count statistics a fragment is 2.5 times larger than in round 1. For (d1)
   the scan reads 107 MB of four atlas fragments to find no cell, against 59 MB of `obs` on the
   store path. In region the table path is still faster (fewer requests, decoding in Arrow); on
   the laptop it was slower (33 s against 18 s). The key-first scan (the filter columns first,
   then every column of the fragments with a match) took (d1) to 3.7 s on the laptop. In region it
   gains a little at 103 datasets and loses a little at 2,060 (section 3.8), so it is an option and
   not the default.
3. **An `obs`-only question needs no store, but most of that gain is column selection.** Counts by
   donor: 3.48 s to 0.56 s and 837 to 131 requests at 32 workers. A store reader that reads only
   the needed columns gets 0.62 s and 121 requests. For the WMG count (g), with nine terms and a
   range on `nnz`, the table takes 2.66 s, the narrow store read 3.23 s, and the store read of
   every column 11.3 s. The table's lead over the narrow read is 1.1 to 1.3 times at 8 and 32
   workers, and 4 to 5 times at 1 worker.
4. **All `obs` is faster on the table path because of fewer objects and decoding outside the
   GIL, not fewer bytes.** 15.9 s to 5.2 s at 32 workers, with 483 MB of fragments against 380 MB
   of Zarr `obs`, but 322 requests against 5,801. The store path keeps about 1.2 of 16 cores busy
   (19 s of CPU in 15.9 s); the table path keeps 2.2 busy. Of the table path's 5.1 s, the scan
   is 1.8 s and the frame of pandas categoricals is 3.2 s. An Arrow-backed frame
   (`dtype_backend="pyarrow"`) takes 0.4 s, for 2.3 s in all and 5.2 GiB peak instead of 7.4 GiB.
5. **`X` dominates `get_anndata`, and the table does not change `X`.** (b) reads about 820 MB of
   `X` chunks for 35,795 scattered cells on both paths, and both use about 17 s of CPU. The table
   path saves the 14 stores' `obs` and `indptr` reads: it is 1.45 times faster at 32 workers (8.0 s
   to 5.5 s), 1.75 times at 8, and 3.5 times at 1 worker, and requests fall from 3,038 to 1,383.
   The gain shrinks as the result grows: (c), 207,185 cells and 4.1 GB of `X`, is 1.18 times
   faster at 32 workers. Every store that the catalog keeps for (b) and (c) holds a hit, so the
   table skips no store there.
6. **Iceberg's own pruning is weak on labels and real but small on numbers.** File selection comes
   from the catalog, through the `dataset_id` partition: Iceberg's 16-character string bounds keep
   71 of 103 files for (b) where the catalog keeps 14 (1,420 against 280 at 2,060). Numeric bounds
   prune files that the catalog cannot: `raw_sum >= 100000` keeps 49 of 103 files, and for (g)
   they drop two of the catalog's 57 datasets, whose largest `nnz` is 276 and 311. Row-group
   statistics prune little, because the rows are in source order. What made the table path fast
   is Arrow's Parquet reader and the catalog's exclusion, not Iceberg's planning. Planning stays
   below 0.3 s up to 2,060 files, with one pyiceberg edge: a `dataset_id IN` term of N values plans
   in time that grows with files times values, so the SDK plans the filter alone and keeps files by
   partition value.
7. **The filter language covers the corpus.** All 50 corpus queries with an `obs` filter bind to
   the columns that the stores hold now. The 9 `tissue_general` queries go through the catalog on
   both paths. The 5 `wmg_production` queries need `!=` and `nnz >= 500`, which only the
   expression filter expresses, and `nnz` only the table tests.
8. **No Iceberg catalog service is necessary** for one writer: an ephemeral SQLite catalog at
   build time, static tables, and a version hint (section 2.2).
9. **The access boundary does not hold for a shared table.** The manifests carry the column bounds
   of every fragment, now with numeric ranges (section 7, item 6).

**What changed from round 1.** `main`'s store path became faster (the `obs` read skips the
barcode index and decodes the cell keys of kept cells only): at 32 workers, (b) 11.0 s to 8.0 s,
(d1) 4.4 s to 2.1 s, (e) 20.5 s to 15.9 s, counts by donor 4.2 s to 3.5 s, and at 2,060 datasets
(d1) 64.6 s to 16.3 s. The table path stayed about the same or became a little slower with the
larger fragments: (b) 6.1 s to 5.5 s, (d1) 0.51 s to 0.61 s, (e) 4.7 s to 5.2 s, and (d1) at 2,060
2.1 s to 6.7 s. So the ratios shrank: (b) 1.8 to 1.45 times, (d1) 8.6 to 3.4 times, (e) 4.4 to
3.1 times. Round 1 reported 0 requests for the obs-only queries of the table path, because it
counted only Zarr reads; with the fragment and metadata reads counted, a count by donor is 131
requests and all `obs` is 322. The new narrow store read (section 3.4) shows that the obs-only
gain is mostly column selection.

## 5. The filter language

pyiceberg's expression language has `AND`, `OR`, `NOT`, the six comparisons, `IN`, `NOT IN`,
`IS [NOT] NULL`, `IS [NOT] NAN`, and `[NOT] LIKE 'prefix%'` (Iceberg's starts-with). Literals are
strings, numbers, and booleans. `pyiceberg.expressions.parser.parse` reads the string form, and
`expressions.to_text` writes it back (tested round trip). `expressions.from_value_filter`
translates a Census 1 `value_filter` string through Python's `ast`.

**The corpus.** `tools/harness/filter_language.py` translates each of the 74 bound queries of the
validation harness and binds it to the fragment schema, which is now the store's `obs` schema
(`tools/results/filter_language.json`):

| Status | Queries | Entries |
|---|---|---|
| Expressible | 49 | `cell_type_slice`, `dataset_slice`, `tissue_slice`, `documented_narrow`, `documented_wide`, `wmg_production`, `cell_type_gene_subset`, `both_layers`, `obs_metadata_only`, `spatial_cell_type` |
| Expressible, with the SDK's membership rewrite | 1 | `multivalued_disease` |
| A gene filter, outside the `obs` table | 5 | `single_gene_all_cells` |
| No `obs` filter (`get_var`, `get_presence_matrix`) | 19 | `var_all_genes`, `var_symbol_collision`, `presence_matrix` |

All 50 queries with an `obs` filter bind to the columns that the stores hold now. The 14 that
waited for a column in round 1 split two ways:

- The 9 `tissue_general` queries (`documented_narrow`, `documented_wide`, `tissue_slice`): the
  catalog indexes `tissue_general`, so the mapping filter expresses them and the catalog prunes
  them, on both paths.
- The 5 `wmg_production` queries: `is_primary_data = true AND nnz >= 500 AND
  cell_type_ontology_term_id != 'unknown'` and six `tissue_general_ontology_term_id != ...`
  terms. The mapping filter cannot express them (`!=` and the range). The catalog tests the
  `is_primary_data` and `!=` terms; it records no numeric range, so only the table tests
  `nnz >= 500`. `nnz` is the one column of the corpus that only the table tests.

The SDK's mapping filter expresses 45 of the 50.

**Gaps.**

1. *Membership on a multi-valued column.* Iceberg has no membership term for a string, and
   pyiceberg has no predicate on the elements of a list column. The SDK therefore gives `=`, `!=`,
   `IN`, and `NOT IN` on the `disease` columns a membership meaning, and rewrites each such term
   into `IN` or `NOT IN` over the joined values that the remaining catalog rows record. The
   rewrite is exact, because the catalog records each dataset's joined values. A joined literal,
   as Census 1 filters write, matches exactly its set of values. Two costs: the design record asks
   for explicit membership syntax, and here the syntax is implicit, on those columns only; and a
   user who reads the table with DuckDB or Polars directly sees string semantics on those columns.
   DuckDB SQL can state membership (`list_contains(string_split(disease, ' || '), 'COVID-19')`), but
   that does not push down into Iceberg's file pruning. A `list<string>` copy of the column in the
   fragment would serve SQL readers; pyiceberg could not filter on it.
2. *No substring, case-insensitive, or pattern match* beyond a prefix. Census 1 has none either.
3. *No ontology closure.* "All T cells" needs the descendants of `CL:0000084`, expanded on the
   client into an `IN` list (for example with `cellxgene_ontology_guide`). This is the same in
   Census 1, and an SDK helper that expands a term gives the same expression in both languages.
4. *No column-to-column comparison and no arithmetic* (`raw_sum / nnz > 2`). Rare in Census use.
5. *The literal types follow the parser*: `is_primary_data = true`, not `True`; numbers are
   longs or doubles. A filter on a column that the table does not hold fails at bind time with a
   clear error. That is better than the store path, where a missing column selects nothing.

**Push-down.** One expression reaches four levels:

1. *Catalog*: `can_match` evaluates each term against the row's value set, so compound boolean
   logic, `NOT IN`, and the `is_primary_data` tri-state all exclude datasets. This is the only
   level that knows the value sets; it is exact for one term and conservative for a
   conjunction (section 4). It records no numeric range, so it keeps every dataset for
   `nnz >= 500`.
2. *Iceberg manifests*: the column bounds of each file prune the files, and the SDK then keeps
   the files of the remaining datasets by their partition value, exactly. The bounds of a label
   column are a lexicographic range, which excludes a file only when the wanted label sorts
   outside it. The bounds of a statistic are exact, and a selective range keeps a fifth to a half
   of the files (section 3.6).
3. *Parquet row groups*: pyarrow prunes row groups of 65,536 rows by their statistics. In source
   order, few row groups are excluded, on labels and on statistics alike.
4. *Rows*: pyarrow applies the expression. The stores get no filter at all: they get row
   ranges. So the per-chunk `obs_index` of the store layout is unnecessary on this path.

## 6. Readers

All four readers read the same table. These measures are from the round 1 stores
(`tools/results/obs_iceberg/r1/ec2-readers.jsonl`, in-region, 16 threads, median of 3, each time
with the open of the table in a new process); this round does not repeat them, and the DuckDB
spike stacks on this branch for SQL readers.

| Reader | How it reads the table | Notes |
|---|---|---|
| pyiceberg 0.12 | `StaticTable.from_metadata`, `table.scan(row_filter, selected_fields)` | The SDK's choice: the expression objects are the filter language, so the catalog push-down, the membership rewrite, and the scan share one tree. Adds pyiceberg (3.4 MB) and its dependencies (`pydantic` and `pydantic-core`, about 9 MB, `mmh3`, `pyroaring`, `strictyaml`, `tenacity`, `cachetools`, `rich`, `click`, `requests`, `fsspec`); pyarrow is already an SDK dependency. Import adds 0.15 to 0.2 s. |
| DuckDB 2.0 `iceberg` | `iceberg_scan('<metadata.json>')` | SQL, with membership through `list_contains`; the extension downloads on first use. An S3 secret through the `aws` extension's credential chain. |
| Polars 2.1 | `pl.scan_iceberg(table)` | Lazy frames; plans through pyiceberg. Brings a 150 MB runtime. |
| pyarrow dataset | `ds.dataset([fragment URIs])` from the catalog | No Iceberg at all: the catalog already lists the fragments and prunes datasets. Reads each footer at scan time. |

Select (c), 207,185 rows: pyarrow over the fragment list 0.57 s, pyiceberg 0.62 s, Polars 0.92 s,
DuckDB 1.20 s. Count by donor over all cells: Polars 0.90 s, DuckDB 1.06 s, pyiceberg 1.06 s,
pyarrow 2.36 s. All `obs` into Arrow: Polars 1.13 s, pyiceberg 1.74 s, DuckDB 2.23 s, pyarrow
3.14 s. Polars planned through pyiceberg, because its own planner is a separate package that was
not installed.

**R.** The `duckdb` R package 1.5.2 reads the table with `iceberg_scan` from local files and from
S3, through the `iceberg`, `httpfs`, and `aws` extensions. No other R Iceberg reader was
evaluated. `arrow` for R reads the fragments as Parquet without Iceberg, as pyarrow does. anndataR
and Rarr remain the path for the stores themselves.

## 7. Deviations from the design record

1. **The fragment holds `x_start` and `x_end`, not only the store row ordinal.** anndata 0.13
   writes `X/indptr` in chunks of the dask row blocks, one shard per chunk, so a full `indptr` read
   is many requests. With the ranges in the fragment, the SDK reads the `X` values of the
   selected cells and never reads `indptr`.
2. **The fragment record is in the store's root attributes**, beside the summary record, and
   holds the fragment's digest, size, and row count. The design record does not say how a builder
   proves that a fragment belongs to its store; the immutable store is the natural place.
3. **A separate table builder** (`tools/build_obs_table.py`) runs after `tools/build_catalog.py`.
   The design record makes one catalog update one Iceberg commit by the catalog builder. The
   spike keeps the two programs apart for clarity, and binds them by the catalog digest in the
   snapshot summary, which `open_census` checks.
4. **One table for each organism**, partitioned by `dataset_id`, as Census 1 has one `obs` for
   each experiment. A query names one organism, and a schema can then differ by organism.
5. **A release pins the metadata file**, which is immutable, rather than a tag of a snapshot in
   the live table. The builder can still tag (`--tag`).
6. **The table does not keep the access model unchanged.** The design record says that a private
   fragment uses the grant of its store and that the access model does not change. The fragment
   does, but the manifests hold each file's column bounds (lower and upper values, null counts,
   row counts) in metadata that every reader of the table reads. A private fragment in a shared
   table discloses its bounds, and with the count statistics the bounds now include numeric
   ranges such as the largest `raw_sum` of a dataset. Options: one table for each access domain
   (a private source brings its own catalog and its own table, which `open_census(obs_tables=[...])`
   already accepts); `write.metadata.metrics.column.<c> = none` for private data, which loses the
   pruning; or a private table that holds only private fragments. The pilot uses public data only.
7. **The result's categories differ.** The design record lists "keep all category labels in a
   result, or keep only labels that the result uses" as a detail that needs approval. The store
   path keeps the labels of every opened store; the table path keeps the labels that the result
   uses. The spike does not choose; the equality tests compare values.
8. **`is_primary_data` joins the mapping filter** on both paths, because the compound query
   needs it and the catalog already records the tri-state. The read of the store path is that of
   `main`.
9. **The table path does not read each store's `var` index.** It checks the number of genes, and
   relies on the catalog's `gene_reference_version` check for the identity of the genes. `main`
   reads and compares the index of every store that it opens.

## 8. Relation to the compaction spike (#34) and to the run read (#40)

- **The store row ordinal stays per dataset.** The compaction spike keeps each member
  contiguous and in its source row order, and its catalog builder writes `compacted_store_uri`,
  `compacted_row_start`, and `compacted_row_stop` into each member's row. A member's row in the
  compacted store is then `compacted_row_start + store_row`, so the fragments and the obs table
  need no rebuild when compaction happens. The member's values in the compacted `X` are at
  `value offset + x_start`, where the value offset is `indptr[compacted_row_start]` of the
  compacted store. One more catalog column for that offset would let the ranges translate with
  no `indptr` read.
- **The table removes the opens of stores without a hit; compaction removes the opens of many
  small stores.** The compaction spike measured a CPU cost of about 28 ms for each store open
  that more threads did not hide. A store that the table proves to hold no hit is thus worth
  more than another thread. With both, a query plans over the fragments, opens only the
  compacted stores that hold a hit, and reads only the hit ranges.
- **They can coalesce.** A compacted store could carry one fragment for all its members (the
  member fragments concatenated, with `store_row`, `x_start`, and `x_end` rebased). One commit
  then deletes the member fragments from the table and adds the compacted fragment: the same
  rows, fewer files to plan and to open. The access boundary limits both in the same way: a
  compacted store and its fragment can hold only datasets of one grant.
- **Iceberg's own compaction** (`rewrite_data_files`) is not the same thing: it writes new data
  files under the table location, a Census-owned copy of the `obs` data, which breaks "one store
  is one access boundary" for the `obs` columns.
- **#40 (open) reads each run of adjacent kept rows as one range** of `X/data` and `X/indices`,
  after one read of `X/indptr`, when a store has at most 256 runs. The fragment gives the same
  ranges with no `indptr` read: a run of adjacent store rows spans from the `x_start` of its
  first row to the `x_end` of its last row. The table path does not yet read by runs. It reads
  all the selected values of a store with one zarr selection of their positions, which builds
  one int64 position for each selected stored value, as anndata's row selection does; it reads
  one slice only when the whole selection is one run. The two compose: #40's run read over the ranges of
  the scan removes both the `indptr` read and, for each run, the position array. #40 measured
  that the cells of one cell type are mostly not adjacent in today's stores, so for those
  queries the gain of either is small; a dataset, donor, or tissue slice is long runs.

## 9. Open questions and future work

- **Materialized consolidation** of fragments into fewer, larger files, or of member fragments
  into a compacted store's fragment (section 8). The plan stays short up to 2,060 files, but the
  scan opens one file for each dataset that it reads, and a scan of the whole table opens every
  fragment.
- **Manifest growth under continuous ingestion.** Each catalog update is one commit and one
  manifest. pyiceberg's fast append keeps every manifest, so a live table after a thousand
  updates has a thousand manifests in its manifest list. Manifest merging, or a periodic rewrite
  of the manifests, is necessary; pyiceberg 0.12 has `merge_append` for the first.
- **Concurrent commits** without a catalog service: the version hint is a plain write. S3
  conditional writes, or the catalog builder's one-writer rule.
- **Private data**: section 7, item 6.
- **The catalog's value sets and the table's statistics overlap.** The catalog excludes datasets
  exactly for one label term, which Iceberg's string bounds cannot. The table makes compound
  queries exact, removes the store `obs` read, and is the only level that prunes a numeric
  range. A design could keep only one of the two exclusion structures; the spike keeps both.
- **Explicit membership syntax** for the multi-valued `disease` columns (section 5, gap 1).
- **The size of the statistics in the fragment.** `raw_mean_nnz` and `raw_variance_nnz` are about
  half of the fragment bytes (section 2.1). A scan that does not select them does not read them,
  but a scan of all `obs` does. `raw_mean_nnz` is `raw_sum / nnz`.
- **The per-chunk `obs_index`** of the store layout is unnecessary on the table path.
- **The `X` read by runs** (section 8, #40) and its position array for scattered cells.
- **`get_anndata_lazy()`** of the design record maps onto the scan: the lazy result can be the
  scan's `obs` plus the per-store ranges, with `X` read on demand.

## 10. Reproduce

Placeholders throughout: `s3://my-bucket/...`, `<account>`, `<region>`, `<repo>`, `<tag>`.

```bash
# 1. The image and the HealthOmics workflow (bulk-ingest-nf/README.md).
docker buildx build --platform linux/amd64 -t <account>.dkr.ecr.<region>.amazonaws.com/<repo>:<tag> \
    -f bulk-ingest-nf/docker/Dockerfile --push .
bulk-ingest-nf/tools/build-healthomics-bundle.sh bundle.zip
# Ingest both prototyping manifests to one outdir, with --check: s3_//my-bucket/spike/ingest/out

# 2. The catalog and the table. The table builder needs pyiceberg[sql-sqlite].
python census2-sdk/tools/build_catalog.py s3://my-bucket/spike/ingest/out/ catalogs/catalog-103.parquet
python census2-sdk/tools/build_obs_table.py catalogs/catalog-103.parquet \
    s3://my-bucket/spike/tables/obs-103 --verify-digest

# 3. Copies of the fragments for the growth with N, and their tables.
cd census2-sdk/tools
python -m bench.fragments catalogs/catalog-103.parquet --out s3://my-bucket/spike/scale/fragments \
    --n-copies 19 --sizes 500,1000,2060 --out-dir catalogs/
for n in 500 1000 2060; do
  python build_obs_table.py catalogs/catalog-$n.parquet s3://my-bucket/spike/tables/obs-$n
done

# 4. The matrix: timed runs, then one counted run of each query (--count-arrow).
Q=results/obs_iceberg/queries.json
SHAPES=single,cell_type,query:compound,query:fp_none,query:fp_some,query:tissue_general,corpus_obs
SHAPES=$SHAPES,count_donor:cell_type,count:wmg:tissue_general_ontology_term_id+cell_type_ontology_term_id
python -m bench.matrix --catalogs catalogs --sizes 103 --shapes $SHAPES --vias stores,table \
    --obs-table 's3://my-bucket/spike/tables/obs-{size}' --workers 1,8,32 --repeats 3 \
    --cell-type neutrophil --queries $Q --where ec2 --out ec2-103.jsonl
python -m bench.matrix --catalogs catalogs --sizes 103 --shapes $SHAPES --vias stores,table \
    --obs-table 's3://my-bucket/spike/tables/obs-{size}' --workers 8 --repeats 1 --count-arrow \
    --cell-type neutrophil --queries $Q --where ec2 --out ec2-103-counted.jsonl
# Variants: --arrow-obs (corpus_obs), --narrow (the store path of the counts), --key-first.

# 5. Growth with N, the pruning, and the false positives.
python -m bench.matrix --catalogs catalogs --sizes 500,1000,2060 \
    --shapes cell_type,query:fp_none,query:fp_some,count_donor:cell_type --vias table \
    --obs-table 's3://my-bucket/spike/tables/obs-{size}' --workers 32 --repeats 3 \
    --cell-type neutrophil --queries $Q --where ec2 --out ec2-scale.jsonl
python -m bench.pruning catalogs/catalog-103.parquet s3://my-bucket/spike/tables/obs-103 \
    --filter "nnz >= 500" --filter "raw_sum >= 100000" --filter "cell_type = 'neutrophil'"
python -m bench.false_positives catalogs/catalog-103.parquet s3://my-bucket/spike/tables/obs-103 \
    --out results/obs_iceberg/false_positives.json
python -m harness.filter_language --out results/filter_language.json
```

`tools/results/obs_iceberg/` holds the records of this round; `r1/` holds those of round 1 (the
stores before the cell filters and the statistics), which section 6 cites.
