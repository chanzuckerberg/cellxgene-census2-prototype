# Spike: DuckDB over the `obs` fragments, without Iceberg

**Status: a spike.** The branch `spike/obs-duckdb` is not intended to merge in its current form.
It is a running-code reference, stacked on `spike/obs-iceberg` (#33, `spikes/obs-iceberg.md`).
The two spikes meet one requirement: one logical query over a release-wide `obs` gives the
selected cells of the whole corpus, so that the SDK opens only the stores with a selected cell
and reads only their rows. The Iceberg spike registers the per-dataset fragments as one Iceberg
table. This spike lets DuckDB read the same fragments directly, through the list of fragments that
the catalog already holds, and it uses SQL as the filter language. The two branches share the
fragment writer, the round-2 ingest, and the benchmark, so `main`'s store path, the Iceberg path,
and this path ran side by side on one instance in one session.

## 1. The idea in one paragraph

The writer puts one Parquet fragment next to each store (#33, section 2.1). The catalog has one
row for each store of a release, so it also lists the fragments of the release and pins their
versions. A query first lets DuckDB evaluate a catalog test over the catalog's value sets, then
reads the fragments of the remaining datasets with one `read_parquet` and the filter as its
`WHERE` clause. The scan gives the cell key, the result `obs`, and the `X` range of each selected
cell. The SDK then opens only the stores with a selected cell and reads their `X` ranges, with the
code of the Iceberg spike. There is no table metadata, no commit, and no table builder.

**Take-home figures** (in-region, 103 datasets and 13.3 million cells, 8 workers, median of 3):

- WMG's production filter (8.4 million cells), counted by tissue and cell type: **0.64 s**
  through DuckDB, 2.66 s through the Iceberg table, and 3.47 s through `main`'s stores with a
  reader that reads only the needed `obs` columns (#33's session); about 5 times the store
  baseline.
- Neutrophils counted by donor: **0.35 s**, 0.54 s, and 0.75 s with the same narrow store read.
- Count all 13.3 million cells by donor: **0.71 s** through DuckDB and 1.78 s through the table,
  at 0.5 and 1.7 GiB peak. The narrow store read of this count is unmeasured.
- A catalog false positive (four datasets pass the catalog, none holds a cell): **0.56 s**,
  0.71 s, 2.29 s in-region; **1.95 s**, 23 s, 27 s on a laptop.
- `get_anndata` with `X`: DuckDB and Iceberg tie (neutrophils 5.65 s and 5.58 s, against 9.24 s
  through the stores). `X` decides.
- The catalog test in DuckDB takes at most 7 ms at 2,060 catalog rows, and it is exact for every
  single-column term, patterns, explicit membership, and ontology closure included.
- The SQL filter language expresses all 50 corpus queries with an `obs` filter, and keeps the
  Census 1 meaning of each.

## 2. What the branch holds

| Part | Where | What it does |
|---|---|---|
| SQL path | `census2_sdk.obs_sql`; `open_census(..., obs_engine="duckdb")` | the filter check, the catalog test, the DuckDB engine, and `get_anndata`, `get_obs`, `count_obs` |
| Fragment names | `census2_sdk.fragments` | the fragment location and row columns, with no engine import; both `obs` paths use it |
| Ontology closure | `tools/build_descendants.py` | the `descendants` table from `cellxgene-ontology-guide` |
| Filter language | `tools/harness/filter_language_sql.py`, `tools/results/obs_duckdb/filter_language.json` | the 74 bound corpus queries and 13 probes, in the mapping, Iceberg, and SQL forms |
| Benchmark | `tools/bench/sql.py`, `variants.py`, `catalog_test.py`, `obs_sql_check.R`; the `sql` via of `bench.matrix` | the measures below; the request counts and the CPU time use #33's `bench.counting` |
| Tests | `tests/test_obs_sql.py`, `tests/test_bench.py` | equality with the store path, normalized layer included; the language, its refusals, and the file-access barrier |

The branch changes no writer code and needs no ingest. It reads the stores, fragments, catalogs,
fragment copies, and tables of #33's round-2 ingest in place.

### 2.1 A query

```python
census = census2_sdk.open_census(["catalog.parquet"], obs_engine="duckdb")
adata = census.get_anndata(
    "homo_sapiens",
    "tissue_general = 'lung' AND cell_type IN ('macrophage', 'alveolar macrophage') "
    "AND is_primary_data AND nnz >= 500",
    layers=["normalized"],
    max_workers=8,
)
```

1. `obs_sql.parse` reads the filter with DuckDB's parser and checks it (section 2.2). A mapping
   filter becomes SQL first (`obs_sql.from_mapping`).
2. A subquery of the `descendants` table becomes the `IN` list of its terms (section 2.4).
3. DuckDB runs the catalog test (section 2.3) over the catalog, which the engine holds as a table.
4. One query reads the fragments of the remaining rows:
   `SELECT * FROM read_parquet([<fragments>], union_by_name = true) WHERE <filter>`. A row's
   fragment is its `obs_fragment_uri`, or the file next to its store.
5. The SDK groups the rows by dataset and gives each group to #33's `_read_cells`, which opens
   the store, reads the `X` ranges of the selected cells, and calculates the normalized layer.
   The result `obs` comes from the scan.

`get_obs` returns the scan as a frame; DuckDB joins the cell key of the index in its own threads.
`count_obs` lets DuckDB group and count, so no per-cell value reaches Python. Neither opens a
store.

The result equals the result of `main`'s store path, value for value, with the normalized layer:
in the tests on synthetic stores, and on the 103 real stores for (d1), (d2), (b), the membership
query, and (c) (`ec2-equal.jsonl`). The three paths gave the same cells, groups, and catalog
candidates for every benchmark query. A categorical column holds the labels that the result uses,
as on the Iceberg path.

### 2.2 The filter language

**What a user writes.** One SQL boolean expression over the `obs` columns: the part of a query
after `WHERE`. A mapping filter still works and becomes SQL.

**How the SDK checks it.** `parse` puts the text in `SELECT 1 FROM obs WHERE <text>` and calls
DuckDB's `json_serialize_sql`, which parses a statement and does not bind or run it. The SDK then
walks the parse tree and refuses every node outside a short list:

- The statement must be one plain `SELECT` of the constant `1` from the base table `obs`, with a
  `WHERE` clause and no other clause. A filter that closes the expression and adds `UNION`,
  `ORDER BY`, `LIMIT`, `GROUP BY`, or a second statement fails this test.
- Expressions: `AND`, `OR`, `NOT`, the comparisons, `BETWEEN`, `IN`, `IS [NOT] NULL`,
  `IS [NOT] DISTINCT FROM`, `LIKE`, `ILIKE`, `GLOB`, `SIMILAR TO`, arithmetic, `CAST` to a scalar
  type, `coalesce`, string, number, boolean, and null literals, and list literals.
- Column names without a table name. The binder of the scan checks that a column exists; an
  unknown column gives `FilterError` with DuckDB's message.
- Functions from `obs_sql.FUNCTIONS` only: about 60 string, list, number, and null functions that
  read only their arguments. Not allowed: table functions (`read_csv`, `read_text`, `query`),
  functions that read state (`getenv`, `current_setting`), functions with side effects
  (`nextval`, `setseed`), aggregates, window functions, lambdas, parameters, and subqueries, with
  one exception for ontology closure (section 2.4).

`render` writes the SQL again from the accepted tree, with a parenthesis around each operation,
and the SDK runs that text, not the text of the user. What the SDK checked is what DuckDB runs.

**The second barrier.** Each engine connection loads `httpfs`, makes its S3 secret, and loads the
catalog. It then sets `allowed_directories` to the fragment directories,
`enable_external_access = false`, and `lock_configuration = true`. A statement can then read no
other file and change no setting (tested). Two limits: `allowed_directories` allows writes too, so
only the language keeps a `COPY` out; and `httpfs` must read the CA bundle, so the bundle is an
allowed path.

**Why the check.** The SDK runs on the user's machine with the user's credentials, so a user who
types a filter gains nothing that DuckDB does not already give them. The check is for two other
cases. First, a filter that does not come from the person who runs the SDK: a web form, a shared
notebook, a service that runs queries for its users, or a language model. Without the check, such
a filter could read local files or other S3 objects with the credentials of the process, write
files, or read an environment variable one bit at a time. Second, the contract: a predicate over
`obs` columns can go to the catalog test, can be translated for another engine, and keeps its
meaning across DuckDB releases only in the small part of SQL that the SDK names. The whole SQL
dialect of one DuckDB version is a much larger contract.

**Membership.** A filter means what SQL means. `disease = 'COVID-19'` compares the whole joined
value, as Census 1 does, so a translated Census 1 filter keeps its meaning. Membership is
explicit:

```sql
list_contains(string_split(disease, ' || '), 'COVID-19')
list_has_any(string_split(disease_ontology_term_id, ' || '), ['MONDO:0100096', 'MONDO:0005109'])
list_has_all(string_split(disease, ' || '), ['HIV infectious disease', 'visceral leishmaniasis'])
```

This is the explicit syntax that the design record asks for, in functions that any DuckDB client
knows. Its cost is small: over the 13.3 million cells of local fragments, `disease = 'normal'`
took 11 ms and the membership form 16 ms, because DuckDB evaluates the function once for each
dictionary entry of a column chunk. A `list<string>` column in the fragment would give shorter
syntax but no speed, so the spike did not make that fragment variant. The split uses the exact
separator of the CxG schema; `catalog.members` also strips spaces around each member, so the two
differ only on values that do not follow the schema.

### 2.3 Push-down

1. **The catalog.** `catalog_test` puts the filter in negation normal form, so that each `NOT`
   applies to one term. A term on one column with a value set in the catalog becomes:

   ```sql
   CASE WHEN tissue IS NULL THEN true ELSE coalesce(list_bool_or(list_transform(
       list_append(tissue, NULL), lambda v: v LIKE 'lung%')), false) END
   ```

   That is: the dataset can match when some value of its value set satisfies the term. DuckDB
   evaluates the term with the same functions as the scan, so `LIKE`, regular expressions,
   membership, `lower`, and the expanded closure push down to the catalog with no second
   implementation in Python. The value set holds no null, so the test adds one, and `IS NULL`
   proves nothing. A row that records no value set keeps its dataset. A term on two columns
   (`raw_sum / nnz > 3`), or on a column that the catalog does not record (`donor_id`, `nnz`),
   proves nothing. `AND` and `OR` combine the results, so the test is exact for one term and
   conservative for a conjunction, as #33's `expressions.can_match` is. For a mapping filter it
   keeps exactly the rows of `catalog.exclude` (tested).
2. **Parquet.** DuckDB pushes comparisons, `IN` lists, and `IS NULL` into the reader of each
   fragment. It skips a row group when the statistics exclude the term, and it uses a bloom filter
   when the file holds one. It reads the filter columns first and fetches the other columns only
   for the row groups with a match, so a scan that selects no cell reads few bytes.
3. **The rows.** DuckDB evaluates the whole filter on each remaining row.

### 2.4 Ontology closure

The language accepts one subquery form, a semi-join on the `descendants` table:

```sql
cell_type_ontology_term_id IN (SELECT descendant FROM descendants WHERE ancestor = 'CL:0000084')
```

The inner `SELECT` must name `descendant` from `descendants` and nothing more, and its `WHERE`
follows the same rules over `ontology`, `ancestor`, and `descendant`. `Options.descendants` names
the Parquet file of the table. Before the catalog test, the engine runs each such subquery and puts
the list of terms in its place, so the catalog test and the Parquet statistics see a plain `IN`
list. `tools/build_descendants.py` builds the table from `cellxgene-ontology-guide` at the ontology
versions of CxG schema 7.0.0, which the writer pins for `tissue_general`: CL, UBERON, MONDO, EFO,
HsapDv, MmusDv, PATO, and HANCESTRO, 918,866 pairs in 2.4 MiB, in 1.4 s. A release would publish it
next to its catalog. The same SQL runs in any DuckDB client that has the table.

### 2.5 How DuckDB finds the fragments, and fragment variants

`Options.files` selects one of four ways, for the comparison:

| `files` | The source of the scan |
|---|---|
| `catalog` (the default) | `read_parquet` of the fragments of the rows that the catalog test keeps |
| `all` | `read_parquet` of every fragment of the organism; only the Parquet statistics prune |
| `glob` | `read_parquet('<prefix>/*.obs.parquet')` and `filename IN (<the catalog's fragments>)` |
| `iceberg` | `iceberg_scan('<metadata.json>')` of #33's table and `dataset_id IN (<the kept datasets>)` |
| `iceberg_all` | `iceberg_scan('<metadata.json>')` with no catalog test: the manifests alone prune |

`tools/bench/variants.py` writes a variant of each fragment with DuckDB, under another prefix, and
a catalog whose `obs_fragment_uri` points at the variants. `bloom` keeps the row order and adds
DuckDB's bloom filters; `sorted` also sorts the rows by `cell_type`, `tissue`, and `store_row`.
The SDK sorts each dataset's selected rows by `store_row` after the scan, so a sorted fragment
gives the same result (tested). Each variant of the 103 fragments took about 5 s in-region.

## 3. Measures

**Setup.** #33's round-2 data: 103 human datasets, ingested with the round-2 writer (the Census 1
cell filters, `tissue_general`, and the five count statistics); 103 of 103 round-trip checks
passed; 13,316,710 cells; 476.9 MiB of fragments (37.5 bytes per cell). Its catalogs at 103, 500,
1,000, and 2,060 rows, its fragment copies, and its Iceberg tables. In-region: one EC2
m6in.4xlarge in us-west-2 (16 vCPU, 64 GiB), with the data in an S3 bucket of the same region.
Laptop: an Apple-silicon laptop on a home connection, shared with two other spikes under a lock.
Each run is a new process; the values are medians of 3 runs (2 on the laptop). "Stores" is
`main`'s store path, "table" the Iceberg path of #33, and "DuckDB" this branch, with DuckDB's
default of one thread for each core. Requests and MB count every object read: the zarr reads, the
catalog, the fragment reads from DuckDB's HTTP log, and, for the table path, the fragment and
metadata reads from separate counting runs, because that counting slows the table path.

The fragment columns differ in size. `observation_joinid` is 28 percent of the bytes, the three
float statistics 56 percent, `nnz` and the `X` ranges 13 percent, and the 21 label columns together
3 percent. A label filter thus reads a small part of a fragment, and the bytes of a selection come
mostly from the columns of the result.

### 3.1 In-region, 103 datasets

Seconds at 1 / 8 / 32 workers (`tools/results/obs_duckdb/ec2-main.jsonl`):

| Query | Cells | Stores opened: stores / table / DuckDB | Stores | Table | DuckDB |
|---|---|---|---|---|---|
| (a) one dataset | 3,799 | 1 / 1 / 1 | 2.22 / 2.22 / 2.17 | 1.03 / 1.01 / 1.10 | 1.01 / 0.94 / 0.99 |
| (b) neutrophils | 35,795 | 14 / 14 / 14 | 28.8 / 9.24 / 7.65 | 9.37 / 5.58 / 5.36 | 9.02 / 5.65 / 5.11 |
| (c) lung macrophages, primary | 207,185 | 5 / 5 / 5 | 35.2 / 24.6 / 24.7 | 25.6 / 21.2 / 21.1 | 25.8 / 21.2 / 21.3 |
| (d1) false positive, no cell | 0 | 4 / 0 / 0 | 6.76 / 2.29 / 2.35 | 0.87 / 0.71 / 0.75 | 0.61 / 0.56 / 0.45 |
| (d2) false positive, 2 of 4 | 111 | 4 / 2 / 2 | 9.78 / 4.42 / 4.40 | 2.50 / 1.89 / 2.13 | 1.93 / 1.63 / 1.60 |
| (f) `tissue_general = 'pancreas'` | 12,815 | 3 / 3 / 3 | 8.27 / 4.34 / 4.41 | 3.36 / 2.80 / 2.99 | 3.29 / 2.66 / 2.71 |
| numeric: (b) with `nnz >= 1000` | 4,113 | – / 11 / 11 | – | 8.62 / 3.51 / 3.22 | 6.32 / 3.50 / 3.08 |
| membership: COVID-19 neutrophils | 1,258 | 6 / 6 / 6 | 11.3 / 3.00 / 3.12 | 3.01 / 1.88 / 1.97 | 3.00 / 1.90 / 1.79 |
| (e) all `obs` | 13,316,710 | 103 / 0 / 0 | 136 / 22.7 / 15.3 ⁴ | 4.96 / 5.04 / 4.98 | 3.76 / 3.73 / 3.84 |
| count (b) by donor | 35,795 | 14 / 0 / 0 | 18.4 / 4.66 / 3.52 | 0.58 / 0.54 / 0.53 | 0.36 / 0.35 / 0.46 |
| count (b) by donor, stores read only the needed columns ¹ | 35,795 | 14 / – / – | 2.53 / 0.75 / 0.62 | | |
| count all by donor | 13,316,710 | 103 / 0 / 0 | 130 / 23.2 / 15.6 ² | 1.75 / 1.78 / 1.68 | 0.68 / 0.71 / 0.65 |
| (g) WMG count by tissue and cell type | 8,365,927 | 57 / 0 / 0 | 83.2 / 15.4 / 11.3 ³ | 2.53 / 2.66 / 2.51 | 0.60 / 0.64 / 0.60 |
| (g), stores read only the needed columns ¹ | 8,365,927 | 57 / – / – | 14.5 / 3.47 / 3.23 | | |

¹ The store path of #33's `bench.run --narrow`, which reads only the `obs` columns of the filter
and of the count, in #33's session on its own instance (`obs_iceberg/ec2-103-narrow*.jsonl`,
`ec2-103-wmg.jsonl`). This is the fair store baseline for an `obs`-only question.
² All `obs` columns of every store; a narrow read of this count was not measured.
⁴ The store path reads each store's `obs` and drops it, and builds no combined frame, so it is
not a like-for-like baseline for (e).
³ #33's session, after its fix of the bench evaluator (`03b2ec7`), which `main`'s path does not
use. This session ran before that fix and measured 82.8 / 20.4 / 16.4 s. No other store-path
shape here uses that evaluator: the other counts take a mapping filter.

At 8 workers: every object read, the MB read, the process CPU time (user and system), and the
peak memory (`ec2-main.jsonl`, `ec2-table-counts.jsonl`):

| Query | Requests: stores / table / DuckDB | MB | CPU s | Peak GiB |
|---|---|---|---|---|
| (a) | 77 / 32 / 27 | 17 / 17 / 17 | 0.9 / 0.6 / 0.6 | 0.7 / 0.5 / 0.5 |
| (b) | 3,038 / 1,383 / 1,489 | 913 / 929 / 919 | 17.2 / 16.2 / 15.9 | 1.5 / 1.9 / 2.1 |
| (c) | 2,084 / 1,254 / 1,390 | 4,131 / 4,134 / 4,134 | 82.5 / 79.4 / 79.4 | 17.5 / 14.3 / 14.5 |
| (d1) | 283 / 19 / 96 | 59 / 108 / 2.1 | 1.3 / 1.4 / 0.9 | 0.6 / 0.8 / 0.4 |
| (d2) | 908 / 231 / 299 | 90 / 121 / 57 | 2.4 / 2.1 / 1.7 | 0.6 / 0.9 / 0.6 |
| (f) | 511 / 150 / 159 | 122 / 124 / 124 | 4.9 / 4.2 / 4.2 | 2.6 / 2.4 / 2.4 |
| numeric | – / 764 / 832 | – / 426 / 425 | – / 7.9 / 7.8 | – / 1.1 / 1.3 |
| membership | 578 / 220 / 218 | 155 / 155 / 154 | 3.7 / 3.2 / 3.2 | 0.6 / 0.7 / 0.7 |
| (e) | 5,801 / 322 / 594 | 380 / 483 / 477 | 18.2 / 11.1 / 21.9 | 1.1 / 7.4 / 7.7 |
| count (b) | 837 / 131 / 113 | 84 / 2.8 / 9.7 | 3.4 / 0.6 / 0.9 | 0.7 / 0.5 / 0.4 |
| count (b), narrow stores ¹ | 121 / – / – | 3.1 / – / – | 0.7 / – / – | 0.3 / – / – |
| count all | 5,801 / 649 / 419 | 380 / 8.9 / 50 | 18.4 / 1.9 / 1.2 | 1.1 / 1.7 / 0.5 |
| (g) | 3,311 / 788 / 593 | 304 / 24 / 50 | 14.0 ³ / 2.7 / 1.8 | 1.4 / 1.4 / 0.5 |
| (g), narrow stores ¹ | 631 / – / – | 31 / – / – | 4.3 / – / – | 0.5 / – / – |

- **The `obs` engine is a small part of a selection with `X`.** For (b), DuckDB's query took about
  0.6 s of the 5.6 s; both engine paths read the same 821 MB of `X` and used 16 CPU-seconds. For
  (c), all three paths read 4.1 GB of `X` and used about 80 CPU-seconds in 21 to 25 s, about 3.7
  cores.
- **Fragment requests differ in kind.** pyarrow coalesces the column chunks of a file into a few
  large ranges, and DuckDB reads each column chunk on its own. A selective filter favours DuckDB,
  which reads the filter columns first: (d1) 2.1 MB in 96 requests against 107 MB in 15 through
  the table; (d2) 29 MB against 93 MB. A narrow scan of all fragments favours pyarrow: a count of
  all cells read 50 MB through DuckDB and 9 MB through the table, because DuckDB reads a remote
  file of up to about 2 MiB whole, in one GET (measured on single fragments), and 47 of the 103
  fragments are below 1 MiB. Where a result needs most columns, both read about the same.
- **(e), all `obs` into pandas,** took 3.7 s through DuckDB: 2.7 s for DuckDB's scan, its Arrow
  export, and the dictionary encoding of the labels, and 0.8 s for the pandas frame. Through the
  table it took 5.0 s: 1.8 s scan and 3.2 s frame. DuckDB used 22 CPU-seconds against 11, and both
  peaked at 7.4 to 7.7 GiB. With Arrow-backed columns (`dtype_backend="pyarrow"`), the table path
  is the faster one, 2.08 s at 5.2 GiB against 3.10 s at 6.3 GiB (`ec2-arrow-obs.jsonl`): pyarrow
  keeps the Parquet dictionaries, and DuckDB decodes them, so the SDK must encode them again.

### 3.2 How DuckDB finds the fragments

DuckDB path, 8 workers. Seconds, and the fragment requests and MB (`ec2-files-*.jsonl`):

| Query | Catalog list (default) | All fragments | Glob over the prefix | `iceberg_scan` |
|---|---|---|---|---|
| (d1) | 0.56 s, 96, 2.1 MB | 0.57 s, 362, 51 MB | 6.55 s, 263, 50 MB | 0.46 s, 96, 2.2 MB |
| (d2) | 1.63 s, 87, 29 MB | 2.02 s, 394, 81 MB | 6.88 s, 256, 77 MB | 1.64 s, 87, 29 MB |
| (b) | 5.65 s, 155, 99 MB | 5.66 s, 473, 141 MB | 11.9 s, 300, 139 MB | 5.61 s, 142, 99 MB |
| numeric | 3.50 s, 117, 107 MB | 3.34 s, 613, 217 MB | 9.48 s, 264, 148 MB | 3.02 s, 107, 107 MB |
| (e) | 3.73 s, 594, 477 MB | 3.60 s, 594, 477 MB | 9.88 s, 555, 477 MB | 3.38 s, 495, 477 MB |
| count all | 0.71 s, 419, 50 MB | 0.61 s, 419, 50 MB | 6.84 s, 380, 50 MB | 0.55 s, 452, 51 MB |

- **Without the catalog test**, DuckDB reads every footer and prunes by row-group statistics alone.
  The time in-region barely changes, but the fragment requests rise 3 to 5 times and the bytes 1.4
  to 24 times for the selective queries.
- **A glob lists the store prefix first.** The fragments share their prefix with the stores, so
  the listing pages through all 31,739 objects of the 103 stores: 32 sequential requests, about
  6 s in-region and 8.6 s from the laptop. With `union_by_name`, DuckDB then reads the footer of
  every listed fragment, so the file-name test saves no request. A glob also does not pin a
  release: it reads what the prefix holds now.
- **`iceberg_scan` of #33's table** is as fast as the catalog list, with a few fewer requests,
  because the manifests give the file sizes and DuckDB sends no HEAD. It needs the `iceberg`
  and `avro` extensions (section 6).

### 3.3 DuckDB threads and memory

DuckDB path, 8 workers, DuckDB threads 4 / 16 / 64 (`ec2-threads.jsonl`). Seconds, and peak GiB
where it changes:

| Query | 4 threads | 16 (the cores) | 64 |
|---|---|---|---|
| (d1) | 0.53 | 0.42 | 0.34 |
| (b) | 5.74 | 5.32 | 5.29 |
| (c) | 21.8 | 21.7 | 21.8 |
| (e) | 5.51 (7.3 GiB) | 3.79 (7.7 GiB) | 3.83 (8.3 GiB) |
| count all | 1.58 | 0.59 | 0.57 |
| (g) | 1.27 | 0.55 | 0.63 |

- The default of one thread for each core is close to the best in-region. More threads add
  parallel range requests but no speed here, and a little memory.
- `memory_limit` did not change the peak (2 GB and 4 GB on round-1 data): the result lives in
  Arrow and pandas, outside DuckDB's buffer manager.
- **One memory edge cost 5 GiB.** An Arrow array that DuckDB exports keeps its whole result chunk
  alive, so the plain strings of the label columns stayed in memory after the SDK had encoded
  them. A copy of each kept column into pyarrow's memory took the peak of (e) from 10.7 to
  5.9 GiB on round-1 data (`ec2-r1-memory*.jsonl`). Round 2 adds five numeric columns, so (e)
  peaks at 7.7 GiB.

### 3.4 Fragment variants

DuckDB path, 8 workers, in-region (`ec2-variants-*.jsonl`; the variant catalogs come from a catalog
build that equals #33's catalog row for row). Fragment MB and requests, and seconds:

| Query | Original (pyarrow) | `bloom` (DuckDB, same order) | `sorted` by cell type and tissue |
|---|---|---|---|
| (b) | 98.5 MB, 152, 6.0 s | 89.1 MB, 189, 6.0 s | 16.3 MB, 48, 5.4 s |
| (c) | 65.4 MB, 156, 22.9 s | 59.5 MB, 184, 22.0 s | 26.2 MB, 84, 22.2 s |
| (d1) | 2.1 MB, 96, 0.38 s | 0.6 MB, 72, 0.40 s | 0.4 MB, 63, 0.33 s |
| (d2) | 28.8 MB, 87, 1.7 s | 25.7 MB, 119, 1.9 s | 3.3 MB, 55, 2.2 s |
| numeric | 107 MB, 117, 3.3 s | 89.1 MB, 151, 3.1 s | 16.3 MB, 45, 2.8 s |
| membership | 7.3 MB, 20, 2.0 s | 6.6 MB, 21, 1.9 s | 4.5 MB, 15, 1.8 s |
| count (b) by donor | 9.7 MB, 113, 0.40 s | 9.1 MB, 149, 0.40 s | 7.4 MB, 40, 0.22 s |

Sorting cuts the fragment bytes of a label query by up to 9 times ((b) 6 times, (d2) 9 times), and
the fragment requests by a third to two thirds, with the same results. In-region the time changes little, because `X` dominates. Bloom
filters save bytes only where a filter value is absent from most row groups, and each filter is
one more request. The variants are 435 MiB against 477 MiB, because DuckDB's writer encodes the
float statistics more compactly.

### 3.5 Laptop

Seconds at 8 workers, median of 2 (`laptop-r2.jsonl`, `laptop-r2-sorted.jsonl`). The link is
shared and noisy; repeats varied by up to a factor of 2, so these are indications:

| Query | Stores | Table | DuckDB | DuckDB, sorted fragments |
|---|---|---|---|---|
| (a) one dataset | | 11.9 | 15.2 | |
| (d1) no cell | 26.7 | 23.1 | 1.95 | 2.0 |
| (d2) 2 of 4 | | 28.2 | 15.1 | 14.9 |
| count (b) by donor | 44.7 | 3.8 | 2.8 | 2.2 |

(d1) is where the engines differ most on this link: the table path reads 107 MB of four atlas
fragments to find no cell, and DuckDB reads 2.1 MB, which gives without a second scan what #33's
key-first scan gives (3.7 s there). (a) and (d2) read `X`, which decides here.

### 3.6 Growth with N

Copies of the 103 fragments (#33, section 3.5): seconds at 32 workers, DuckDB / table
(`ec2-scale.jsonl`). The store path at 2,060 is in #33.

| N | (b) neutrophils | (d1) no cell | numeric | count (b) by donor | count all by donor | peak GiB, count all |
|---|---|---|---|---|---|---|
| 103 | 5.11 / 5.36 | 0.45 / 0.75 | 3.08 / 3.22 | 0.46 / 0.53 | 0.65 / 1.68 | 0.5 / 1.7 |
| 500 | 19.4 / 19.8 | 0.45 / 1.17 | 10.7 / 10.9 | 0.61 / 1.78 | 2.74 / 7.75 | 0.8 / 6.2 |
| 1,000 | 41.0 / 41.8 | 0.85 / 2.22 | 21.2 / 22.1 | 1.25 / 3.46 | 4.72 / 15.5 | 1.1 / 11.8 |
| 2,060 | 89.0 / 90.5 | 1.53 / 3.89 | 44.3 / 46.1 | 1.99 / 6.99 | 10.9 / 31.2 | 1.8 / 23.7 |

- A selection with `X` grows with the stores that it opens on both paths: (b) at 2,060 reads
  16.4 GB of `X` from 280 stores.
- An `obs`-only question grows with the fragments that it reads, and DuckDB stays 2.5 to 3.5
  times faster: a count of all 266 million cells by donor took 10.9 s at 1.8 GiB peak against
  31.2 s at 23.7 GiB, because DuckDB aggregates as it reads.
- (d1) at 2,060: the catalog keeps 80 atlas fragments, and DuckDB reads 42 MB of them in 1,920
  requests and 2.3 CPU-seconds, against 27 CPU-seconds through the table.
- The catalog test stays below 7 ms at 2,060 rows (`catalog_test.jsonl`, in process on the
  laptop); Python's `can_match` takes up to 180 ms there, and `catalog.exclude` about 30 ms.

### 3.7 R

The `duckdb` R package 2.0.0 installed in 4 s as a binary from Posit Package Manager on the
instance. `tools/bench/obs_sql_check.R` runs the SDK's catalog-test SQL and filter SQL over the
round-2 catalog and fragments (`ec2-r-check.txt`): (c) by donor, 5 candidates and 207,185 cells,
in 0.27 s; COVID-19 membership in lung, 1 candidate and 125,530 cells, in 0.06 s; a count of all
13,316,710 cells, with the catalog's file list in a SQL variable, in 0.20 s; the connection with
its extensions in 0.30 s. The counts equal the Python counts.

### 3.8 File-level planning at scale

Can the engine plan the files, with no catalog test? Four ways for DuckDB to find the fragments of
an `obs`-only query (`get_obs`, or `count_obs` for the count), at N = 103, 500, 1,000, and 2,060
catalog rows (#33's fragment copies, catalogs, and tables):

- **catalog**: the catalog test, then `read_parquet` of the kept fragments (the default).
- **footers**: `read_parquet` of all N fragments, no catalog test. DuckDB reads each footer and
  prunes row groups by their statistics. This is with `union_by_name = false`, the faster
  setting: with `union_by_name`, DuckDB reads every footer before the scan, and the same queries
  took 2 to 4 times as long (`ec2-plan-footers-union.jsonl`).
- **iceberg**: `iceberg_scan` of #33's table of N fragments, no catalog test; DuckDB prunes files by
  the manifest bounds. No pyiceberg in the client.
- **sorted footers**: as footers, over copies of the `sorted` variant (section 3.4) at each N.

The queries: (d1), (d2), a rare cell type (`motor neuron`, which one source dataset holds, 2,955
cells), the numeric query, the membership query, and (b) counted by donor. In-region, 8 workers,
median of 3 (`ec2-plan-*.jsonl`). Seconds and requests at N = 2,060:

| Query | catalog | footers | iceberg | sorted footers |
|---|---|---|---|---|
| (d1) | **2.10 s**, 1,920 | 4.33 s, 7,240 | 2.38 s, 3,385 | 8.22 s, 6,340 |
| (d2) | **2.46 s**, 1,720 | 5.09 s, 7,897 | 3.03 s, 4,120 | 8.97 s, 6,430 |
| rare cell type | **1.50 s**, 1,017 | 5.48 s, 8,347 | 5.26 s, 5,675 | 5.94 s, 4,800 |
| numeric | **4.13 s**, 3,275 | 7.84 s, 12,598 | 6.42 s, 9,824 | 4.89 s, 4,959 |
| membership | **1.15 s**, 370 | 3.98 s, 6,589 | 3.27 s, 3,857 | 5.02 s, 4,822 |
| count (b) by donor | **2.27 s**, 2,270 | 5.57 s, 8,744 | 4.31 s, 6,010 | 5.30 s, 4,820 |

The laptop, 2 runs, at N = 103 / 500 / 1,000 / 2,060 (`laptop-plan-*.jsonl`); seconds:

| Query | catalog | footers | iceberg |
|---|---|---|---|
| (d1) | **1.9 / 4.9 / 6.7 / 13.7** | 5.6 / 17.2 / 49.7 / 93.4 | 3.8 / 7.8 / 14.9 / 30.7 |
| rare cell type | **4.1 / 14.8 / 31.4 / 54.6** | 12.3 / 41.2 / 80.1 / 149.8 | 11.9 / 37.1 / 55.5 / 151.5 |
| count (b) by donor | **2.4 / 8.9 / 14.5 / 27.6** | 7.2 / 27.8 / 53.4 / 110.5 | 12.1 / 17.5 / 41.0 / 85.4 |

**A cost model.** The least-squares slope of the seconds and of the requests against N, per 1,000
catalog rows, with the median over the queries (`tools/bench/plan_cost.py`):

| | catalog | footers | iceberg | sorted footers |
|---|---|---|---|---|
| in-region | 0.78 s, 880 requests | 2.29 s, 3,940 | 1.59 s, 2,370 | 2.60 s, 2,380 |
| laptop | 12.6 s, 930 requests | 52.8 s, 4,050 | 39.1 s, 2,750 | |

- **The catalog test still earns its place.** At 2,060 rows, which is more than the 1,258 stores
  of Census 1 LTS, the footer plan took 2 to 4 times as long in-region and 3 to 7 times as long on
  the laptop, with 4 to 18 times the requests and up to 24 times the bytes. Each fragment that the
  catalog can exclude still costs a HEAD, a footer GET or two, and, for a fragment below about
  2 MiB, a whole-file read.
- **Iceberg's manifests are a middle ground.** They save the HEAD and the footer of each file that
  the bounds exclude, so `iceberg_scan` came within 1.1 to 1.2 times of the catalog for the false
  positives (d1) and (d2) in-region. But label bounds are lexicographic ranges, so for a rare cell
  type, membership, and the count the plan kept most files: 1.9 to 3.5 times the catalog's time
  in-region, and 2.2 to 3.1 times on the laptop at 2,060.
- **Sorting the fragments made the footer plan cheaper only where the filter is on a sort column
  and few datasets match**: the numeric query (on `cell_type`) took 4.9 s against 7.8 s. For (d1)
  and (d2) the sorted copies were slower in this session (8.2 s against 4.3 s), with fewer
  requests; the copies are new objects under a new prefix, and the spike did not find the cause.
- **CPU and memory follow the files read**: at 2,060, (d1) used 2.6 CPU-seconds and 0.6 GiB through
  the catalog, 8.4 and 1.7 GiB through footers, and 3.9 and 0.9 GiB through `iceberg_scan`.
- The catalog's only gap is numeric terms; there DuckDB's footers or the manifests must prune.

## 4. Findings

1. **For a selection with `X`, the two `obs` engines tie.** In-region, the DuckDB path and the
   Iceberg path took the same time for every `get_anndata` query, within the spread of repeats:
   (b) 5.65 s and 5.58 s, (c) 21.2 s both, at 8 workers. Both open no store for a catalog false
   positive and read no store `obs`, which helps most at low concurrency: (b) 9.0 s against
   28.8 s through the stores at 1 worker. A large selection is bound by `X` on all three paths.
2. **For an `obs`-only question, DuckDB over the files is the fastest of the three and the
   smallest in memory.** It groups and counts in the engine. Against a store reader that reads
   only the needed `obs` columns (#33's session, section 3.4), WMG's filter, counted, took 0.64 s
   against 3.47 s, about 5 times faster, and the table 2.66 s; neutrophils by donor 0.35 s
   against 0.75 s, and the table 0.54 s. The narrow store read still opens every store that the
   catalog keeps. A count of all cells by donor took 0.71 s against 1.78 s through the table,
   and at 2,060 datasets 10.9 s at 1.8 GiB against 31.2 s at 23.7 GiB; a narrow store read of
   that count is unmeasured. The full store read (23.2 s at 103) is not a fair baseline for an
   `obs`-only question.
3. **The catalog is a sufficient file index.** It lists the fragments, pins their versions, and
   records each dataset's value sets. A catalog test that DuckDB evaluates over the value sets
   takes 1 to 7 ms at 2,060 rows, and it is exact for every single-column term (section 5). The
   SDK's query needs no table metadata, no builder, and no commit. What the table adds for the
   SDK is file-level bounds of the numeric columns, which the catalog does not record yet (#33,
   section 3.6).
4. **On a laptop, bytes decide, and DuckDB reads fewer of them for a selective filter.** DuckDB
   reads the filter columns first and the other columns only where a row group has a match: (d1)
   in 1.95 s from 2.1 MB, against 23 s through the table and 27 s through the stores. For a
   narrow scan of every fragment it reads more bytes than pyarrow, because it reads a small file
   whole.
5. **Sorted fragments cut the bytes of a label query by up to 9 times** ((b) 6 times), with no
   change to the stores, and the same results. The time in-region changes little. The sort would help the table
   path too, which reads the same files.
6. **A glob is not a release pin and is slow.** It lists every object under the store prefix (32
   requests, about 6 s) and, with `union_by_name`, reads every footer anyway.
7. **The defaults suffice for concurrency.** One DuckDB thread for each core was close to the
   best; DuckDB used 1.7 to 6 cores on `obs`-only questions, and `X` decoding set the CPU of the
   rest. `memory_limit` does not bound the result. The SDK needed one fix for DuckDB's memory: a
   copy of the kept Arrow columns, which saved 5 GiB on all `obs`.
8. **SQL is a practical filter language, with a cost.** It expresses every corpus query, gives
   explicit membership, patterns, arithmetic, and ontology closure, and its catalog test needs no
   rule for each kind of term. The cost: the language is a DuckDB dialect, so the SDK must pin
   DuckDB and keep an allow list (section 5).
9. **The access boundary holds.** A fragment read uses the grant of its store, and the only shared
   structure is the catalog, as on `main`. The Iceberg table's manifests disclose every fragment's
   column bounds to every reader (#33, section 7, item 6).
10. **File-level planning cannot yet be left to the engine.** Without the catalog test, DuckDB's
    footers took 2 to 4 times the catalog plan's time in-region and 3 to 7 times on the laptop at
    2,060 rows, and Iceberg's manifests through `iceberg_scan` 1.1 to 3.5 times and 2.2 to 3.1
    times (section 3.8). Per 1,000 catalog rows, the laptop pays 12.6 s through the catalog,
    39 s through the manifests, and 53 s through the footers. The catalog's per-dataset value
    sets exclude a dataset with no request at all; an engine must still touch each file or each
    manifest entry, and its label bounds are weak.

## 5. The filter language

**The corpus.** `tools/harness/filter_language_sql.py` translates each of the 74 bound queries of
the validation harness with `obs_sql.from_value_filter` and binds it to the fragment schema of the
writer (`tools/results/obs_duckdb/filter_language.json`):

| | Queries with an `obs` filter | Expressible |
|---|---|---|
| SQL (this branch) | 50 | 50 |
| Iceberg expressions (#33) | 50 | 50 |
| Mapping filter (`main`) | 50 | 45: not `wmg_production`, which needs `!=` and `nnz >= 500` |

The other 24 bound queries have no `obs` filter (19) or filter genes only (5). Each SQL translation
keeps the Census 1 meaning, also of `multivalued_disease`, which names a joined value. The corpus
is almost all conjunctions of equality, so it does not separate the languages. The probes do.
Each ran over the 103 round-2 fragments. "Kept" is the number of catalog rows that each
language's catalog test keeps, and "with cells" the number of datasets with a selected cell:

| Kind | Mapping | Iceberg | SQL | Kept: mapping / Iceberg / SQL | With cells | Cells |
|---|---|---|---|---|---|---|
| compound: `(tissue = 'lung' OR tissue_general = 'lung') AND NOT cell_type = 'macrophage'` | no | yes | yes | – / 17 / 17 | 17 | 1,885,447 |
| `IN` | yes | yes | yes | 22 / 22 / 22 | 22 | 122,097 |
| `NOT IN`, with `tissue_general = 'lung'` | no | yes | yes | – / 17 / 17 | 17 | 1,822,235 |
| `raw_variance_nnz IS NULL` | no | yes | yes | – / 103 / 103 | 0 | 0 |
| one member of `disease` | implicit | implicit, rewritten | explicit | 9 / 9 / 9 | 9 | 637,638 |
| a joined value, as Census 1 | no | yes | yes | – / 0 / 0 | 0 | 0 |
| numeric range with a label: (b) with `nnz >= 2000 AND raw_sum < 20000` | no | yes | yes | – / 14 / 14 | 6 | 1,077 |
| numeric range only: `nnz >= 8000` | no | yes | yes | – / 103 / 103 | 58 | 279,441 |
| prefix: `LIKE 'CD8%'` | no | yes | yes | – / 24 / 24 | 24 | 427,977 |
| substring, case-insensitive: `ILIKE '%macrophage%'` | no | no | yes | – / – / 36 | 36 | 581,205 |
| regular expression | no | no | yes | – / – / 27 | 27 | 1,051,329 |
| ontology closure: descendants of T cell | no | no; a client-side list | yes, a semi-join | – / – / 46 | 46 | 3,658,513 |
| arithmetic: `raw_sum / nnz > 3.0` with `tissue_general = 'blood'` | no | no | yes | – / – / 15 | 15 | 2,047,708 |

- **Every single-column term is exact at the catalog** in SQL, `ILIKE`, regular expressions,
  membership, and the closure included, because DuckDB evaluates the term itself over the value
  set. #33's `can_match` is exact for its terms too, but each new kind of term needs a Python rule.
  Numeric terms prove nothing at the catalog, which records no numeric range: it keeps 14 datasets
  for the numeric range, and 6 hold a cell; it keeps all 103 for `nnz >= 8000`, and 58 hold a
  cell. There the fragment statistics, or #33's manifest bounds, must do the work.
- **Membership is explicit and exact**, and a Census 1 filter that names a joined value keeps its
  meaning. The cost is length: `list_contains(string_split(disease, ' || '), 'COVID-19')` against
  `disease = 'COVID-19'`. A user who reads the fragments with DuckDB, Polars, or R sees the
  semantics of the SDK, which is not so on the Iceberg path, where only the SDK rewrites `=`.
- **Nulls.** The fragment holds a NaN of the store as a null (pyarrow's `from_pandas`), so
  `IS NULL` selects what `isna()` selects in the store. No round-2 cell has a null.
- **Gaps of the SQL subset**: no lambdas (so no pattern test on each member), no `CASE`, no
  aggregates, and no subquery other than the closure. Each could join the allow list; a `CASE` is
  safe, and a lambda needs scope rules in the check.

**Who owns the language.** With SQL, DuckDB owns the grammar, the parser, and the meaning of each
function, and the SDK owns the allow list. A DuckDB release can change the parse-tree format (the
literal form of `json_serialize_sql` changed in 2.0, and `parse` reads only that form) or a
function, so the SDK must pin a DuckDB range and test the language at each upgrade, as it does for
anndata. With Iceberg expressions, pyiceberg owns a small grammar, and the SDK owns the membership
meaning on top. A third option is a Census grammar that the SDK parses itself and renders to SQL,
to Iceberg expressions, and to R: more work, and no engine in the contract.

## 6. Packaging and deployment

- **The dependency.** `duckdb` 2.0 is one wheel with no Python dependency: 24 MB to download for
  Linux x86-64 (65 MB installed), and a 50 MB library on macOS arm64. It ships type information,
  so `mypy --strict` needs no ignore for it. pyiceberg 0.12 is a 0.8 MB wheel, but it brings
  pydantic, pydantic-core, mmh3, pyroaring, strictyaml, tenacity, cachetools, rich, click, requests,
  and fsspec. pyarrow is already an SDK dependency on both paths. The SQL path imports no
  pyiceberg.
- **Extensions.** Parquet, JSON, and ICU are in the library. A remote fragment needs `httpfs`
  (25.7 MiB for linux_amd64), DuckDB's credential chain needs `aws` (68.1 MiB), and `iceberg_scan`
  needs `iceberg` (54.4 MiB) and `avro` (27.0 MiB). DuckDB downloads each one from its extension
  repository on first use, into `~/.duckdb/extensions/<version>/<platform>`: in-region in 0.24 to
  0.48 s each; on the laptop, round 1 took 21 s once for the Iceberg set. A load then takes 20 to
  250 ms. The engine installs and loads them before it disables external access, because a locked
  connection cannot load an extension.
- **Offline and air-gapped.** An image can run `INSTALL httpfs` at build time, copy the files into
  the extension directory (`SET extension_directory`), `LOAD` a file by path, or point
  `custom_extension_repository` at a mirror. The extensions are signed, so a copied file loads
  without `allow_unsigned_extensions`. A community project publishes the extension files as wheels
  (`duckdb-extension-httpfs` and others); the spike did not test it. Without `httpfs`, the SDK
  could give DuckDB the fragments through an fsspec file system or pyarrow, at a speed that the
  spike did not measure.
- **Credentials.** `Options.credentials`: `env` makes a DuckDB secret from the `AWS_*` variables,
  which obstore also reads for the stores; `chain` is DuckDB's credential chain through the `aws`
  extension (it found the instance profile on EC2); `auto` takes `env` when the variables are set.
  A chain that finds nothing fails at `CREATE SECRET`, not at the first read. A secret holds the
  credentials of its creation, so expired temporary credentials need a new connection, or the
  chain's refresh option, which the spike did not test.
- **R.** The `duckdb` R package runs the same SQL (section 3.7). A table function cannot take a
  subquery, so `read_parquet((SELECT list(...)))` fails to bind; a SQL script puts the catalog's
  file list in a variable instead:

  ```sql
  SET VARIABLE files = (SELECT list(replace(store_uri, '.zarr', '.obs.parquet'))
                        FROM read_parquet('catalog.parquet')
                        WHERE organism = 'homo_sapiens' AND <catalog test>);
  SELECT donor_id, count(*) FROM read_parquet(getvariable('files')) WHERE <filter> GROUP BY 1;
  ```

  So the SDK's `obs` query, catalog pruning included, is two SQL statements that any DuckDB client
  can run. The R package shares its extension and secret directory, `~/.duckdb`, with the CLI and
  the other clients.

## 7. Iceberg against DuckDB over files

| | Iceberg table over the fragments (#33) | DuckDB over the catalog's fragment list (this branch) |
|---|---|---|
| What a release publishes | the catalog, and one table (metadata file, manifest list, manifests) for each organism | the catalog |
| Builder and commits | a table builder after the catalog builder; one commit for each catalog update; manifests grow with updates | none |
| Catalog service | none needed with a static table and a version hint, for one writer | none |
| Release pinning | the immutable metadata file, checked against the catalog digest | the catalog: each row names its store, so its fragment |
| Schema evolution | Iceberg's schema and name mapping; an added column is one commit | `union_by_name` at read time: a fragment without a column gives nulls. Nothing records the schema, and the scan reads every footer before the first row |
| File-level statistics | manifest bounds for every fragment, without a footer read; strings truncated to 16 characters; numeric bounds prune files | none before the footer; the catalog's per-dataset value sets do that work, exactly for one term, but hold no numeric range |
| Row groups | pyarrow's statistics | DuckDB's statistics and bloom filters, and late fetch of the other columns |
| Access boundary | the manifests hold the column bounds of every fragment, which every reader of the table sees | nothing beyond the catalog; a footer read needs the fragment's grant |
| Readers | every Iceberg engine (Spark, Trino, Snowflake, DuckDB, Polars, pyiceberg); R through DuckDB | every Parquet reader, given the catalog's list; the catalog test is SQL that any DuckDB client runs |
| Filter language | Iceberg expressions; the SDK gives `=` a membership meaning on `disease` | a checked subset of DuckDB SQL; membership is explicit |
| Python dependencies | pyiceberg and its dependencies; pyarrow | duckdb (one wheel) and its `httpfs` extension; pyarrow |

The two do not exclude each other. The fragments are the same files, and the results through both
are equal. A Census could keep the catalog as the contract that the SDK reads, and also publish an
Iceberg table over the same public fragments for engines that want one, without the SDK reading
it.

## 8. Deviations from the design record

1. **No Iceberg.** The design record names Iceberg as the candidate package for the per-cell `obs`
   frame, for manifest bounds, snapshot tags, and schema evolution. This branch tests the simpler
   form that the same section starts from: engines that read a list of fragments as one table. The
   catalog already lists and pins the fragments, and its value sets select datasets more exactly
   than label bounds.
2. **The filter interface is SQL**, a checked subset of the DuckDB dialect. The design record
   leaves the interface open and asks for explicit membership syntax and a translator from Census 1
   filters; both are here (`list_contains`, `obs_sql.from_value_filter`).
3. **A new release artifact, the `descendants` table**, and one subquery form in the language, for
   ontology closure. The design record does not name closure.
4. **The catalog test runs in DuckDB**, over the catalog's list columns. The design record plans
   bitmaps for the value sets; the test would then need bitmap functions or a list projection.
5. **The result schema.** With `union_by_name`, a column that only some fragments hold has nulls for
   the others, where the store path keeps only the columns that every opened store holds. The
   round-2 fragments are of one writer version, so the results agree; a mixed-vintage release needs
   a rule.
6. **The result's categories** are the labels that the result uses, as on the Iceberg path. The
   design record lists this as a detail that needs approval.

## 9. Open questions

- **One filter language or two.** SQL suits `obs`-only questions and readers outside the SDK.
  `get_anndata` still needs a selection interface (the design record notes this). A Census-owned
  grammar that renders to SQL and to Iceberg expressions would keep the contract independent of one
  engine.
- **Numeric ranges in the catalog, or not.** Numeric terms prove nothing at the catalog today;
  DuckDB's footers or Iceberg's manifests prune those files. Section 3.8 suggests that the
  catalog's value sets should stay for the label terms; whether five numeric ranges join them is
  a separate choice.
- **Fragment layout.** Which columns to sort by, and whether the writer or a backfill sorts. The
  sort helps both engine paths.
- **Mixed vintages** (section 8, item 5).
- **Credential refresh** for long sessions with temporary credentials (section 6).
- **Private sources.** The access boundary holds here, but a shared catalog still lists the value
  sets of a private dataset: the design record's open item on private discovery.

## 10. Reproduce

Placeholders throughout: `s3://my-bucket/...`, `<region>`. The stores, fragments, catalogs, copies,
and tables are those of `spikes/obs-iceberg.md`, section 10.

```bash
pip install -e ./census2-dataset-writer -e './census2-sdk[dev,s3,iceberg,duckdb]'
export AWS_REGION=<region>   # credentials in the AWS_* variables, or an instance profile
cd census2-sdk/tools

# The head-to-head: the store path, the Iceberg table, and DuckDB over the fragments.
python -m bench.matrix --catalogs catalogs --sizes 103 \
    --shapes single,cell_type,query:compound,query:fp_none,query:fp_some,query:tissue_general,query:numeric,query:membership,corpus_obs,count_donor:cell_type,count_donor:all,count:wmg:tissue_general_ontology_term_id+cell_type_ontology_term_id \
    --vias stores,table,sql --obs-table 's3://my-bucket/spike/tables/obs-{size}' \
    --workers 1,8,32 --repeats 3 --cell-type neutrophil \
    --queries results/obs_duckdb/queries.json --where ec2 --out ec2-main.jsonl
# The same shapes again: the table path's requests, in runs of their own.
python -m bench.matrix <as above> --vias table --workers 8 --repeats 1 --count-arrow --out ec2-table-counts.jsonl
# DuckDB threads, and the ways to find the files.
python -m bench.matrix <as above> --vias sql --workers 8 --duckdb-threads 4,16,64 --out ec2-threads.jsonl
python -m bench.matrix <as above> --vias sql --workers 8 --files glob --glob 's3://my-bucket/spike/ingest/out/*.obs.parquet'
python -m bench.matrix <as above> --vias sql --workers 8 --files iceberg --iceberg 's3://my-bucket/spike/tables/obs-{size}'
# Growth with N.
python -m bench.matrix --catalogs catalogs --sizes 500,1000,2060 --vias sql,table --workers 32 <shapes, queries, table>
# Fragment variants under another prefix, and their catalogs.
python -m bench.variants catalogs/catalog-103.parquet --variant sorted \
    --out s3://my-bucket/spike/variants/sorted --catalog-out variants/sorted/catalog-103.parquet
python -m bench.matrix --catalogs variants/sorted --sizes 103 --vias sql --workers 8 <shapes, queries>
# The catalog tests at each size, the descendants table, and the filter languages.
python -m bench.catalog_test catalogs/catalog-{103,500,1000,2060}.parquet --queries results/obs_duckdb/queries.json
python ../tools/build_descendants.py descendants.parquet
python -m harness.filter_language_sql --catalog catalogs/catalog-103.parquet \
    --descendants descendants.parquet --out results/obs_duckdb/filter_language.json
# File-level planning at each N (section 3.8), and its cost model.
python -m bench.matrix --catalogs catalogs --sizes 103,500,1000,2060 \
    --shapes obs:fp_none,obs:fp_some,obs:rare,obs:numeric,obs:membership,count_donor:cell_type \
    --vias sql --workers 8 --queries results/obs_duckdb/queries.json --files all --no-union-by-name
python -m bench.matrix <as above> --files iceberg_all --iceberg 's3://my-bucket/spike/tables/obs-{size}'
python -m bench.plan_cost catalog=ec2-plan-catalog.jsonl footers=ec2-plan-footers.jsonl iceberg=ec2-plan-iceberg.jsonl
# R.
Rscript bench/obs_sql_check.R s3://my-bucket/spike/catalogs/catalog-103.parquet queries-sql.json
```
