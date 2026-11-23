# Spike: store consolidation at full scale, and the whole `obs`

This file records the full-scale round of the store consolidation spike. `store-consolidation.md`
holds the pilot on 103 stores, the design of a consolidated store, and the sharp edges.

The question: Census 2 reads the whole `obs` of an organism much more slowly than Census 1. Does
consolidation make that read faster, at the pilot's rule (about half the stores) and at two more
aggressive rules (about 200 and about 100 stores)?

The answer is no. On the full corpus, the whole human `obs` takes about 420 s from the
per-dataset stores and 16 s from Census 1. The pilot's rule, with a whole read of each consolidated
store, takes 397 to 414 s; the two aggressive rules take 455 to 617 s. The consolidated layout does
not address the costs that set the time: the concatenation of one part for each dataset (about
147 s, in every layout), and per-dataset work that uses about one core of 64.

A scan-heavy `get_anndata` (section 5) is different: there `X` sets the time, Census 2 from the
per-dataset stores takes 1.2 to 1.6 times the time of Census 1, and consolidation makes the scan of
part of the corpus slower.

## 1. Setup

- **Corpus.** The 1,258 datasets of Census 1 LTS `2025-11-08`: 1,086 human and 172 mouse,
  202,636,280 cells. `bulk-ingest-nf` from `main` at 2ab7c32 wrote one store for each, with
  `--check`, and every check passed. The stores hold 875 GiB.
- **The read.** Census 1: `cellxgene_census.get_obs(census, organism)`, every column. Census 2:
  `Census.get_obs(organism)`, every column. `get_obs` is in #51, and `main` does not yet hold it.
  The branch holds #51 and reads a dataset from the consolidated store that its catalog row names.
  `bench.run --shape obs` runs the Census 2 read; `census1_obs.py` runs the Census 1 read.
- **The machine.** EC2 m6in.16xlarge in us-west-2 (64 vCPU, 256 GiB). Each run is a new process.
  The numbers are medians of 3 runs at 32 workers; the runs at 8 workers are single runs. Each
  configuration ran on an instance of its own, and each instance also read the human corpus from the
  per-dataset stores once: 420.7 to 425.0 s, so the instances agree to about 1%.

## 2. Census 1

| stack | organism | rows | seconds | read to Arrow | to pandas | CPU s | peak GiB |
|---|---|---|---|---|---|---|---|
| tiledbsoma 2.3.0, pandas 2.3.3 | human | 158,982,719 | 43.4 | 15.4 | 27.9 | 190 | 47.7 |
| tiledbsoma 2.3.0, pandas 2.3.3 | mouse | 43,653,561 | 11.5 | 4.7 | 6.8 | 52 | 13.7 |
| tiledbsoma 2.2.0, pandas 3.1.1 | human | 158,982,719 | 16.3 | 15.0 | 1.2 | 163 | 30.5 |
| tiledbsoma 2.2.0, pandas 3.1.1 | mouse | 43,653,561 | 4.6 | 4.3 | 0.3 | 45 | 9.6 |

The two organisms hold the 202,636,280 cells of the corpus. The pandas version sets most of the
Census 1 time. With pandas 2, the conversion to pandas makes a Python object for each string; with
pandas 3, the strings stay in Arrow buffers, and the conversion takes 1.2 s in place of 27.9 s. The
Census 2 SDK runs pandas 3.1.1, so the pandas 3 rows are the like-for-like baseline. A plain install
of `cellxgene-census` 1.18.0 resolves tiledbsoma 2.3.0 and pandas 2.3.3; with `pandas>=3`, the
resolver takes tiledbsoma 2.2.0. The read to Arrow is the same with either: 15 s for human, on
about 10 cores of CPU.

## 3. The consolidation policies

The planner packs, first-fit decreasing by store bytes, the stores of each organism, with at most
64 members in a group. `--max-members` is new; each member adds about 345 bytes to the
consolidation record, which every open parses. `halve` is the pilot's rule.

| policy | small below | group at most | stores | consolidated stores | members | copy, share of corpus | task-hours | longest task | peak memory |
|---|---|---|---|---|---|---|---|---|---|
| per-dataset | | | 1,258 | | | | | | |
| `halve` | 256 MiB | 1 GiB | 578 | 61 | 741 | 59 GiB, 6.8% | 4.0 | 5.6 min | 1.5 GiB |
| `t3` | 3 GiB | 3 GiB | 203 | 142 | 1,197 | 425 GiB, 48.5% | 20.8 | 11.4 min | 1.7 GiB |
| `t8` | 8 GiB | 8 GiB | 101 | 84 | 1,241 | 635 GiB, 72.5% | 28.8 | 24.2 min | 3.1 GiB |

Every consolidated store passed `verify` in its task, and the catalog builder's checks. The memory
of a task is about 0.7 GiB for each million cells, well under the request of the pipeline.

The first run of each policy started all of its groups at once, the three runs together, 287 tasks.
S3 returned `SlowDown` (503) on the prefix of the stores, and obstore's 10 retries ended within
seconds, so 15 groups failed twice and were ignored. A second run of those groups alone wrote them.
A full-scale consolidation therefore needs fewer tasks at the same time, or a longer retry budget in
the obstore client.

## 4. Census 2

Seconds of `get_obs` for the whole `obs`, with the time in `anndata.concat`, the CPU seconds, and
the requests at the zarr store interface. `members` reads each member of a consolidated store
apart; `whole` reads the store's `obs` columns once (below).

| configuration | stores | human, 32 workers | 8 workers | concat s | CPU s | requests | mouse, 32 | 8 |
|---|---|---|---|---|---|---|---|---|
| Census 1, pandas 3 | | 16.3 | | | 163 | | 4.6 | |
| `main` with #51 | 1,086 | 430.6 | 481.9 | 146 | 506 | 42,408 | 106.9 | 113.4 |
| per-dataset (branch) | 1,086 | 422.6 | 405.8 | 147 | 499 | 42,408 | 104.9 | 100.0 |
| `halve`, members | 510 | 433.1 | 413.6 | 145 | 523 | 65,703 | 103.8 | 101.0 |
| `halve`, whole | 510 | 414.4 | 397.0 | 147 | 495 | 30,612 | 101.0 | 97.5 |
| `t3`, members | 176 | 545.6 | 484.7 | 145 | 674 | 109,228 | 110.6 | 105.1 |
| `t3`, whole | 176 | 506.7 | 454.7 | 145 | 619 | 50,620 | 108.7 | 105.2 |
| `t8`, members | 85 | 617.2 | 546.1 | 149 | 766 | 130,889 | 114.8 | 109.7 |
| `t8`, whole | 85 | 581.0 | 495.9 | 149 | 715 | 65,595 | 109.6 | 107.4 |

The stores and the requests are those of the human read. Every configuration gives the same frame:
the digests of the index and of each column with its dtype are equal for both organisms. The peak
memory of a human read is 66 to 81 GiB.

The branch's per-dataset path differs from `main` with #51 in that it reads the `obs` columns of a
dataset at the same time; at 32 workers the two are within 2%.

### 4.1 The whole read of a consolidated store

A member read cuts the member's rows from each column, and each partial read of a sharded array
also reads its shard index. For the whole corpus, that made 1.5 to 3 times the requests of the
per-dataset stores. The branch therefore reads a consolidated store's `obs` columns whole, once,
when the selected datasets hold at least half of its rows (`WHOLE_STORE_SHARE`), and each member
cuts its rows from memory. The result is the same, and the tests compare both paths with the
per-dataset stores.

The whole read halves the requests of a member read, and makes the read 4 to 9% faster for human
and 2 to 5% for mouse. It is not enough: only `halve` is faster than the per-dataset stores, by 2
to 4%.

The obs arrays of a consolidated store have chunks of 8,192 rows, the layout of the design record,
so that a member read reads few rows that it does not need. A whole read then decodes many small
chunks. For the largest `t8` store (3 datasets, 4.59 M cells), a whole read of the `obs` columns
decodes 15,708 chunks with 5.0 to 5.2 CPU seconds; the same columns from the three per-dataset
stores, which have anndata's chunks, are 458 chunks and 1.74 CPU seconds. That is about 0.22 ms of
CPU for each chunk, and about 120 CPU seconds over the `t8` stores of the corpus: most of what the
`t8` whole read adds to the per-dataset read.

### 4.2 Where the time goes

- **The concatenation: about 147 s of 420 s, in every layout.** `get_obs` concatenates one part for
  each dataset, so consolidation cannot change it. A sampling profile of the mouse read (py-spy,
  only the thread that holds the GIL) puts 66% of its samples in `anndata.concat`. Of all samples,
  53% are the uniqueness check of `obs_names` that `AnnData.__init__` makes on the result, and 23%
  are the check of `get_obs` itself on the same index (`result.index.is_unique`). Each turns the
  Arrow strings of the cell key into Python objects for a hash table (about 30% of the samples).
  The cell key is unique by construction: each store's `observation_joinid` is checked unique, and
  the catalog refuses a `dataset_id` in two rows. The profiler fell behind, so these shares are
  approximate.
- **The reads: about 270 s, on about one core.** The load average of the instance was 1.1 to 1.4 of
  64 cores during the reads. 8 workers were as fast as 32 or faster in every configuration except
  `main` with #51, which reads the columns of a dataset one after another. The payload is 4.6 GB,
  about 17 MB/s. In the profile, the read functions hold the GIL in few samples, so decoding does
  not appear to be the limit. What limits the reads is not yet known.
- **The pandas version of Census 1.** Census 1 converts to pandas in 1.2 s with pandas 3 and in
  27.9 s with pandas 2. Census 2 already runs pandas 3.

### 4.3 Opportunities

Not tried here; each one applies to every layout.

- Skip the uniqueness checks that the cell key makes redundant, in `get_obs` and in the
  `AnnData` that `anndata.concat` builds, or concatenate the frames without an `AnnData`. The two
  checks hold most of the GIL samples of the mouse read; the concatenation alone is a third of the
  time.
- Find what limits the reads. 32 workers keep about one core busy and move about 17 MB/s.
- If whole reads are the common case, give the `obs` arrays of a consolidated store larger chunks,
  for example one chunk for each shard of 65,536 rows. The per-dataset read sets the bound: the
  consolidated layout saves at most the opens, 15 to 30 ms of CPU each in the pilot, or about 20 s
  for the 1,000 stores fewer of `t8`.

## 5. A scan-heavy `get_anndata`

`get_obs` of the whole `obs` is one shape. The other is a `get_anndata` query that the catalog keeps
for many datasets, but that matches few cells, so the read scans much `obs` and keeps little.

The query (`consolidation-anndata-queries.json`, `rare_19`): 19 rare cell types, chosen from the
catalog's cell counts for many datasets kept per matching cell. The catalog keeps 530 of the 1,086
human datasets, which hold 108.7 M cells (68% of the human cells); 992,451 cells match (0.6%). The
same filter also ran through `get_obs`, which makes the same exclusion and the same `obs` scan, and
reads neither `X` nor `var`. Medians of 3 runs at 32 workers, single runs at 8 workers:

| layout | stores | `get_anndata`, 32 workers | 8 workers | `get_obs`, 32 workers | 8 workers | `get_anndata` requests | read |
|---|---|---|---|---|---|---|---|
| Census 1 (either pandas) | | 302.3 | | | | | |
| `main` | 530 | 452.4 | 374.3 | | | 215,821 | 84.8 GB |
| per-dataset (branch) | 530 | 483.7 | 372.0 | 60.3 | 55.2 | 205,076 | 84.8 GB |
| `halve` | 317 | 510.2 | 371.6 | 71.0 | 61.9 | 209,002 | 83.8 GB |
| `t3` | 152 | 609.9 | 440.3 | 137.4 | 114.1 | 212,691 | 79.8 GB |
| `t8` | 79 | 623.3 | 478.4 | 205.4 | 153.5 | 215,950 | 77.9 GB |

`main` at 2ab7c32 has no `get_obs`. Every layout gives the same result (the digests of `obs` and
of `X` are equal), and the `get_obs` frames are equal too. Census 1 (`get_anndata` with
`obs_value_filter`) takes 302 s with pandas 2.3.3 and 306 s with pandas 3.1.1, for the same 992,451
cells, on about 15 cores (4,460 CPU seconds); Census 2 uses 2 to 3 cores.

- **On this path Census 2 is close to Census 1.** From the per-dataset stores, Census 2 takes 1.2
  to 1.6 times the time of Census 1: 372 s at 8 workers and 484 s at 32, against 302 s; `main`
  takes 1.24 and 1.50 times. The whole `obs` read is 26 times. Census 1 also spends most of this
  query in `X`.

- **`X`, not `obs`, sets the time.** The `obs` scan of the 530 stores takes 55 to 60 s; the read of
  `X` for 0.6% of the cells takes the other 320 to 420 s and reads about 85 GB. The rare cells are
  scattered, and each run of kept rows reads whole chunks of `X/data` and `X/indices`. In a
  per-dataset store, anndata's chunks hold about 150,000 values, about 30 cells at about 5,000
  stored values per cell, so a run of one rare cell reads and decodes two chunks of about 0.6 MB
  for about 40 KB of values. A consolidated store's `X` chunks hold 131,072 values, so the layout
  does not change this.
- **The `obs` scan of few kept cells is cheap.** 60 s for 108.7 M cells, where the whole `obs` of
  the corpus took 422 s through `get_obs`. The difference is the per-cell work on kept cells: the
  cell key strings, the uniqueness checks, and the concatenation.
- **Consolidation makes the partial scan slower.** `get_obs` for the filter: 60 s per-dataset,
  71 s `halve`, 137 s `t3`, 205 s `t8`. A consolidated store reads whole when the selected datasets
  hold half of its rows, so it decodes the rows of the other members too, in small chunks; a store
  below that share reads each member apart, with a shard index for each partial read. The requests
  grow from 21,430 to 68,887.
- **The shared `var` index of a consolidated store does not show.** Each per-dataset store reads
  its own `var` index (289 KB, 78,894 strings), and a consolidated store reads one for all of its
  members. By a local timing, the saving is about 40 ms of CPU for each store fewer, which the
  slower scan hides.
- **8 workers are faster than 32.** 372 s against 484 s for the per-dataset stores. Each dataset
  read has up to `RUN_THREADS` (16) span reads in flight, so 32 workers keep up to 512 in flight.
  `main`, which reads the `obs` columns of a dataset one after another, is faster at 32 workers
  than the branch (452 s against 484 s), and the two are equal at 8.

## 6. Reproduce

```bash
# Census 1, with the default resolve and with pandas 3:
uv run --with cellxgene-census python spikes/store-consolidation/census1_obs.py homo_sapiens
uv run --with cellxgene-census --with 'pandas>=3' \
    python spikes/store-consolidation/census1_obs.py homo_sapiens

# The plans (consolidate-nf/tools), from the per-dataset catalog:
consolidation-plan.py catalog-pd.parquet --small-mib 256 --target-mib 1024 --max-members 64 \
    --sizes sizes.json --out plan-halve.jsonl          # and 3072/3072, 8192/8192
# Consolidate each plan with consolidate-nf and --check, then build a catalog for each:
python census2-sdk/tools/build_catalog.py s3://my-bucket/stores/ catalog-halve.parquet \
    --consolidated s3://my-bucket/consolidated/halve/

# Census 2, from census2-sdk/tools, one configuration per instance:
python -m bench.run catalog-pd.parquet --shape obs --organism homo_sapiens --workers 32 \
    --per-dataset --digest
python -m bench.run catalog-halve.parquet --shape obs --organism homo_sapiens --workers 32

# The table:
python spikes/store-consolidation/obs_summary.py census2-sdk/tools/results/consolidation-obs-runs.jsonl

# Section 5: the query is rare_19 in census2-sdk/tools/results/consolidation-anndata-queries.json.
python -m bench.run catalog-pd.parquet --shape filter --organism homo_sapiens \
    --filter "$(jq -c .rare_19.filter consolidation-anndata-queries.json)" --workers 32 --per-dataset --digest
python -m bench.run catalog-pd.parquet --shape obs --organism homo_sapiens \
    --filter "$(jq -c .rare_19.filter consolidation-anndata-queries.json)" --workers 32 --per-dataset
uv run --with cellxgene-census python spikes/store-consolidation/census1_anndata.py \
    homo_sapiens rare_19.json      # a file that holds the filter
```

The `members` rows come from the branch at 0aade26, before the whole read of a consolidated store;
on the current branch, `census2_sdk.WHOLE_STORE_SHARE` above 1 gives that read. The `main` rows come
from #51 at 7d578cd, and the `main` rows of section 5 from 2ab7c32. The records are in
`census2-sdk/tools/results/`: `consolidation-obs-runs.jsonl`, `consolidation-anndata-runs.jsonl`,
and `consolidation-full-plans.json`.
