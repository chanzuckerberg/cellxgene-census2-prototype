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
census = census2_sdk.open_census(
    ["catalog.parquet"], obs_tables=["s3://my-bucket/census/obs/homo_sapiens"]
)
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

pyiceberg's expression language has `AND`, `OR`, `NOT`, the six comparisons, `IN`, `NOT IN`,
`IS [NOT] NULL`, `IS [NOT] NAN`, and `[NOT] LIKE 'prefix%'` (Iceberg's starts-with). Literals are
strings, numbers, and booleans. `pyiceberg.expressions.parser.parse` reads the string form, and
`expressions.to_text` writes it back (tested round trip). `expressions.from_value_filter`
translates a Census 1 `value_filter` string through Python's `ast`; that is the translator that
the design record's migration error and migration guide need.

**The corpus.** `tools/harness/filter_language.py` translates each of the 74 bound queries of the
validation harness and binds it to the fragment schema (`tools/results/filter_language.json`):

| Status | Queries | Entries |
|---|---|---|
| Expressible now | 35 | `cell_type_slice`, `dataset_slice`, `cell_type_gene_subset`, `both_layers`, `obs_metadata_only`, `spatial_cell_type` |
| Expressible, with the SDK's membership rewrite | 1 | `multivalued_disease` |
| Expressible when the store holds the column | 14 | `documented_narrow`, `documented_wide`, `tissue_slice` (`tissue_general*`); `wmg_production` (`tissue_general*`, `nnz`) |
| A gene filter, outside the `obs` table | 5 | `single_gene_all_cells` |
| No `obs` filter (`get_var`, `get_presence_matrix`) | 19 | `var_all_genes`, `var_symbol_collision`, `presence_matrix` |

Every one of the 50 queries with an `obs` filter is expressible. The SDK's mapping filter
expresses 45 of them: it cannot express `wmg_production`, which needs `!=` and `nnz >= 500`.
`wmg_production` is the one compound query in the corpus, and it is a conjunction; the corpus has
no `OR`, `NOT IN`, or range other than that `nnz` term, so the corpus alone does not exercise most
of what the language adds. The tests do (`tests/test_obs_table.py`).

**Gaps.**

1. *Membership on a multi-valued column.* Iceberg has no membership term for a string, and
   pyiceberg has no predicate on the elements of a list column. The SDK therefore gives `=`, `!=`,
   `IN`, and `NOT IN` on the `disease` columns a membership meaning, and rewrites each such term
   into `IN` or `NOT IN` over the joined values that the remaining catalog rows record. The
   rewrite is exact, because the catalog records each dataset's joined values. A joined literal,
   as Census 1 filters write, matches exactly its set of values. Two costs: the design record asks
   for explicit membership syntax, and here the syntax is implicit; and a user who reads the table
   with DuckDB or Polars directly sees string semantics on those columns. DuckDB SQL can state
   membership (`list_contains(string_split(disease, ' || '), 'COVID-19')`), but that does not
   push down into Iceberg's file pruning. A `list<string>` copy of the column in the fragment
   would serve SQL readers; pyiceberg could not filter on it.
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
   conjunction (section 4).
2. *Iceberg manifests*: the scan adds `dataset_id IN (remaining)`, which the identity partition
   prunes exactly. The column bounds of each file then prune further, but the bounds of a label
   column are a lexicographic range, which excludes a file only when the wanted label sorts
   outside it. They earn more on numeric columns, such as the count statistics, and on sorted
   columns.
3. *Parquet row groups*: pyarrow prunes row groups of 65,536 rows by their statistics; in arrival
   order, few row groups are excluded on label columns.
4. *Rows*: pyarrow applies the expression. The stores get no filter at all: they get row
   ranges. So the per-chunk `obs_index` of the store layout is unnecessary on this path.

## 6. Readers

All four readers read the same table (laptop test on local files; EC2 measures in section 3):

| Reader | How it reads the table | Notes |
|---|---|---|
| pyiceberg 0.12 | `StaticTable.from_metadata`, `table.scan(row_filter, selected_fields)` | The SDK's choice: the expression objects are the filter language, so the catalog push-down, the membership rewrite, and the scan share one tree. Adds pyiceberg (3.4 MB) and its pure-Python dependencies (pydantic, which anndata's dependencies already bring, `mmh3`, `pyroaring`, `strictyaml`, `tenacity`, `cachetools`, `rich`, `click`, `requests`); pyarrow is already an SDK dependency. Import adds 0.15 to 0.2 s. |
| DuckDB 2.0 `iceberg` | `iceberg_scan('<metadata.json>')` | SQL, with membership through `list_contains`; the extension downloads on first use (21 s once). An S3 secret through the `aws` extension's credential chain. |
| Polars 2.1 | `pl.scan_iceberg(table)` | Lazy frames; plans through pyiceberg. Brings a 150 MB runtime. |
| pyarrow dataset | `ds.dataset([fragment URIs])` from the catalog | No Iceberg at all: the catalog already lists the fragments and prunes datasets. Reads each footer at scan time. |

**R.** The `duckdb` R package 1.5.2 (a binary from the rocker image's repository, installed in
7 to 12 s) reads the table with `iceberg_scan` from local files and from S3, through the `iceberg`,
`httpfs`, and `aws` extensions (`r-check` result in section 3). No other R Iceberg reader was
evaluated. `arrow` for R reads the fragments as Parquet without Iceberg, as pyarrow does. anndataR
and Rarr remain the path for the stores themselves (design record, existing user support).

## 7. Deviations from the design record

1. **The fragment holds `x_start` and `x_end`, not only the store row ordinal.** anndata 0.13
   writes `X/indptr` in chunks of the dask row blocks, one shard per chunk: 180 objects for a store
   of 566,546 cells, 3 s to read from a laptop. With the ranges in the fragment, the SDK reads the
   `X` values of the selected cells and never reads `indptr` (measured in section 3).
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
   the live table. The builder can still tag (`--tag`). A static reader reads one metadata file,
   so the pin needs no live pointer.
6. **The table does not keep the access model unchanged.** The design record says that a private
   fragment uses the grant of its store and that the access model does not change. The fragment
   does, but the manifests hold each file's column bounds (lower and upper labels, null counts,
   row counts) in metadata that every reader of the table reads. A private fragment in a shared
   table discloses its bounds. Options: one table for each access domain (a private source brings
   its own catalog and its own table, which `open_census(obs_tables=[...])` already accepts); or
   `write.metadata.metrics.column.<c> = none` for private data, which loses the pruning; or a
   private table that holds only private fragments. The pilot uses public data only.
7. **The result's categories differ.** The design record lists "keep all category labels in a
   result, or keep only labels that the result uses" as a detail that needs approval. The store
   path keeps the labels of every opened store, an accident of which stores the catalog keeps; the
   table path keeps the labels that the result uses. The spike does not choose; the equality
   tests compare values.
8. **`is_primary_data` joins the mapping filter** on both paths, because the compound query
   needs it and the catalog already records the tri-state. The read of the store path is that of
   `main`.

## 8. Relation to the compaction spike

- **The store row ordinal stays per-dataset.** A member of a compacted store keeps its rows
  contiguous and in order, so its row in the compacted store is `member row offset + store_row`.
  Its values in the compacted `X` are at `member value offset + x_start`. The compaction spike
  should record both offsets for each member, the row offset and the value offset
  (`compacted indptr[member row offset]`), so that the ranges in the fragments translate without
  a read of the compacted `indptr`. Then the obs table needs no change when compaction happens,
  and the SDK chooses the compacted store or the member store per query.
- **The table removes the open cost of non-hit stores; compaction removes the open cost of many
  small stores.** They compose: with both, a query plans over the fragments, opens only the
  compacted stores that hold a hit, and reads only the hit ranges.
- **They can coalesce.** A compacted store could carry one materialized fragment for all its
  members (the member fragments concatenated, with `store_row`, `x_start`, and `x_end` rebased).
  One commit then deletes the member fragments from the table and adds the compacted fragment:
  the same rows, fewer files to plan and to open. That is the materialized consolidation of
  fragments that this spike holds short of. The access boundary limits both in the same way: a
  compacted store and its fragment can hold only datasets of one grant.
- **Iceberg's own compaction** (`rewrite_data_files`) is not the same thing: it writes new data
  files under the table location, a Census-owned copy of the `obs` data, which breaks "one store
  is one access boundary" for the `obs` columns.

## 9. Open questions and future work

- **Materialized consolidation** of fragments into fewer, larger files, or of member fragments
  into a compacted store's fragment (section 8). The plan time grows with N (section 3), and the
  scan opens one file for each dataset that it reads.
- **Manifest growth under continuous ingestion.** Each catalog update is one commit and one
  manifest. pyiceberg's fast append keeps every manifest, so a live table after a thousand
  updates has a thousand manifests in its manifest list. Manifest merging, or a periodic rewrite
  of the manifests, is needed; pyiceberg 0.12 has `merge_append` for the first.
- **Concurrent commits** without a catalog service: the version hint is a plain write. S3
  conditional writes, or the catalog builder's one-writer rule.
- **Private data**: section 7, item 6.
- **Statistics columns.** When the stores hold `raw_sum`, `nnz`, and the other statistics, the
  column bounds in the manifests give numeric range pruning without more work, and the catalog's
  five numeric ranges become a second copy of the same facts.
- **The catalog's per-column test and the table overlap.** The catalog excludes datasets exactly
  for one term, which Iceberg's bounds cannot. The table makes compound queries exact and removes
  the store `obs` read. A design could keep only one of the two exclusion structures; the spike
  keeps both.
- **The per-chunk `obs_index`** of the store layout is unnecessary on the table path.
- **`get_anndata_lazy()`** of the design record maps onto the scan: the lazy result can be the
  scan's `obs` plus the per-store ranges, with `X` read on demand.

## 10. Reproduce

Placeholders throughout: `s3://my-bucket/...`, `<account>`, `<region>`, `<repo>`, `<tag>`.

```bash
# 1. The image and the HealthOmics workflow (bulk-ingest-nf/README.md).
docker buildx build --platform linux/amd64 -t <account>.dkr.ecr.<region>.amazonaws.com/<repo>:<tag> \
    -f bulk-ingest-nf/docker/Dockerfile --push .
bulk-ingest-nf/tools/build-healthomics-bundle.sh bundle.zip
# Ingest both prototyping manifests to one outdir: s3_//my-bucket/spike/ingest/out

# 2. The catalog and the table.
python census2-sdk/tools/build_catalog.py s3://my-bucket/spike/ingest/out/ catalog-103.parquet
python census2-sdk/tools/build_obs_table.py catalog-103.parquet s3://my-bucket/spike/tables/obs-103

# 3. Copies of the fragments for the growth with N, and their tables.
cd census2-sdk/tools
python -m bench.fragments ../../catalog-103.parquet --out s3://my-bucket/spike/scale/fragments \
    --n-copies 19 --sizes 500,1000,2060 --out-dir catalogs/
for n in 500 1000 2060; do
  python build_obs_table.py catalogs/catalog-$n.parquet s3://my-bucket/spike/tables/obs-$n
done

# 4. The matrix (queries.json holds the named mapping filters of section 3).
python -m bench.matrix --catalogs catalogs --sizes 103 \
    --shapes single,cell_type,query:compound,query:fp_none,query:fp_some,corpus_obs,count_donor:cell_type \
    --vias stores,table --obs-table 's3://my-bucket/spike/tables/obs-{size}' --workers 1,8,32 \
    --repeats 3 --cell-type neutrophil --queries queries.json --where ec2 --out results/ec2-103.jsonl
python -m bench.false_positives catalogs/catalog-103.parquet s3://my-bucket/spike/tables/obs-103
python -m bench.readers catalogs/catalog-103.parquet s3://my-bucket/spike/tables/obs-103 \
    --reader duckdb --query-shape select --query "$(jq -c .compound queries.json)"
```
