# Spike: store consolidation at full scale, and the whole `obs`

This file records the full-scale round of the store consolidation spike. `store-consolidation.md`
holds the pilot on 103 stores, the design of a consolidated store, and the sharp edges.

The question: Census 2 reads the whole `obs` of an organism much more slowly than Census 1. Does
consolidation make that read faster, at the pilot's rule (about half the stores) and at two more
aggressive rules (about 200 and about 100 stores)?

## 1. Setup

- **Corpus.** The 1,258 datasets of Census 1 LTS `2025-11-08`: 1,086 human and 172 mouse,
  202,636,280 cells. `bulk-ingest-nf` from `main` at 2ab7c32 wrote one store for each, with
  `--check`, and every check passed. The stores hold 875 GiB.
- **The read.** Census 1: `cellxgene_census.get_obs(census, organism)`, every column. Census 2:
  `Census.get_obs(organism)`, every column. `get_obs` is in #51, and `main` does not yet hold it.
  The branch holds #51 and reads a dataset from the consolidated store that its catalog row names.
  `bench.run --shape obs` runs the Census 2 read; `census1_obs.py` runs the Census 1 read.
- **The machine.** EC2 m6in.16xlarge in us-west-2 (64 vCPU, 256 GiB). Each run is a new process.
  The numbers are medians of 3 runs.

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

Not yet measured. The measurement instance had no read access to the private bucket of the stores.
The environments, the catalogs of the four layouts, and the driver are ready.

Two facts from small reads predict part of the result. A `get_obs` of a per-dataset store makes
about 34 requests. The same 3 datasets cost 102 requests from their per-dataset stores and 197 from
their consolidated store, because each partial read of a sharded array also reads its shard index
(as in the pilot, section 7 of `store-consolidation.md`). For the whole `obs`, the SDK reads each
member of a consolidated store apart. Consolidation thus removes opens but adds requests.

## 5. Reproduce

```bash
# Census 1:
python spikes/store-consolidation/census1_obs.py homo_sapiens --census-version 2025-11-08

# Census 2: one catalog for each layout (the catalog builder, with --consolidated), then
PYTHONPATH=census2-sdk/tools python -m bench.run catalog.parquet --shape obs \
    --organism homo_sapiens --workers 32 [--per-dataset] [--digest]

# The table:
python spikes/store-consolidation/obs_summary.py census2-sdk/tools/results/consolidation-obs-runs.jsonl
```

The records are in `census2-sdk/tools/results/consolidation-obs-runs.jsonl` and
`consolidation-full-plans.json`.
