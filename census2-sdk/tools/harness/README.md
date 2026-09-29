# The validation harness

Census 1 is the oracle. The harness records its answer to a fixed set of queries, asks a
stack under test the same queries, and reports every difference. The Census 2 SDK is
the stack under test (tracker task 6.2). This directory holds the harness pieces that
the SDK's run needs. The rest of the Phase 2 code stays in the `cellxgene-census`
repository, branch `spike/validation-harness`, `tools/validation_harness/`.

| File | Task | Role |
| --- | --- | --- |
| `corpus.py`, `query_values.json` | 2.1 | 16 entries, 74 bound queries over the seven experiments, with pinned filter values |
| `oracle_runner.py`, `requirements-oracle.lock`, `snapshot_manifest.json` | 2.2 | Records the Census 1 answers in the canonical form |
| `canonical.py` | 2.3 | The canonical result form: rows by identity, decoded categoricals, sorted CSR, `.npy` parts, a digest |
| `comparator.py` | 2.3 | The diff between two canonical results. Its corruption suite is `tests/test_harness.py` |

Not ported: the fixture reader over converted TileDB data (`zarr_adapter.py`,
`make_fixtures.py`, `compare_stacks.py`), the benchmark rig (2.5), and the catalog
prototype. The writer's stores and the SDK replace the first group. The rig returns in
Iteration 3.

## The pins

| What | Pin | Where recorded |
| --- | --- | --- |
| Census release | `2025-11-08`, the last LTS build | `query_values.json`, `_source` |
| `cellxgene-census` | 1.18.0 | `requirements-oracle.lock` |
| `tiledbsoma` | 1.15.3, the builder's pin | `requirements-oracle.lock` |
| `anndata` | 0.13.2 | `requirements-oracle.lock` |
| Python | 3.12 | this page |
| Query corpus | v1.0, 16 entries, 74 bound queries | `corpus.py` |
| Snapshot bound | 3 smallest datasets per experiment | `snapshot_manifest.json` |

Every serialized result embeds this record in its `manifest.json`, so a snapshot names
the oracle that produced it (evidence 11.1).

## Produce the oracle snapshot

Reads are anonymous. From `census2-sdk/tools/`:

1. `python3.12 -m venv harness/.venv-oracle`
2. `harness/.venv-oracle/bin/pip install -r harness/requirements-oracle.lock`
3. `harness/.venv-oracle/bin/python -m harness.oracle_runner run`. This writes
   `harness/oracle_snapshot/`, one directory per query with `manifest.json` and `.npy`
   parts, about 13 GB in total.
4. `harness/.venv-oracle/bin/python -m harness.oracle_runner run --verify` executes
   every query again and compares the result digests. The last line must read
   `PASS: 0 digest mismatches`.

The snapshot is not committed. Steps 3 and 4 gave 74 of 74 identical digests on
2026-08-24 from a clean environment on one host. That run proves the steps and the
pins on that host.

`snapshot_manifest.json` is committed. A change to the dataset sample changes every
digest, so regenerate it (`oracle_runner build-manifest`) only with a corpus version
bump.

## Scope limits

- The top-level SDK accessors do not reach the spatial experiments, so a spatial query
  runs through `experiment.axis_query` on the same handle. Each `manifest.json`
  records which path ran.
- An obs-filtered query is clamped to the manifest's datasets. A full-corpus X read is
  the benchmark rig's subject (task 2.5).
- The deferred corpus entries (identity coordinates, spatial scenes) have no oracle
  serialization. They wait on decision 7 and Phase 3.
- The filters are TileDB-SOMA `value_filter` strings. The Census 2 SDK takes a mapping
  of indexed column to values. Its adapter answers the entries it can express and
  reports the rest as unsupported.
