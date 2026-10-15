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
| Query corpus | v1.1, 16 entries, 74 bound queries | `corpus.py` |
| Snapshot bound | 3 datasets per experiment: the 3 smallest, or the human pick below | `snapshot_manifest.json` |

Every serialized result records the release, its URI, the package versions, and the
corpus version in its `manifest.json`. A result thus names the oracle that produced it
(evidence 11.1).

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

The snapshot is not committed. For corpus v1.1, steps 3 and 4 gave 74 of 74 identical
digests on 2026-10-02 from a clean environment on one host. That run proves the steps and
the pins on that host. 64 of the 74 digests equal those of v1.0. The 10 others are the
human queries that the new human sample changes.

`snapshot_manifest.json` is committed. A change to the dataset sample can change the digest
of every clamped query, so regenerate it (`oracle_runner build-manifest`) only with a
corpus version bump.

## The human sample of corpus v1.1

The SDK run (task 6.2) compares Census 1 with the Iteration 1 stores, so the human sample
must hold datasets that both stacks have. 73 of the 96 stores are in release 2025-11-08.
The human sample is three of them, chosen so that the silent-failure entries select
cells on both stacks:

| Dataset | Cells | Reason |
| --- | --- | --- |
| `5b25a46f-f1cf-4caa-9e10-6ab81dfe9174` | 12,822 | The only one of the 73 with a multi-valued disease. All its cells are primary data. |
| `6e00ccf7-0749-46ef-a999-dba785630d52` | 5,499 | The smallest dataset with primary data and brain or lung cells |
| `62315937-e268-4fa5-a032-8f7d776f3a3f` | 4,778 | The control: the SDK and Census 1 agreed on every value of this dataset (task 6.1) |

To rebuild the manifest with this pick:

```sh
harness/.venv-oracle/bin/python -m harness.oracle_runner build-manifest \
  --pick census_data/homo_sapiens=5b25a46f-f1cf-4caa-9e10-6ab81dfe9174,6e00ccf7-0749-46ef-a999-dba785630d52,62315937-e268-4fa5-a032-8f7d776f3a3f
```

The other six experiments keep the 3 smallest datasets. Four human pins in
`query_values.json` come from the Census 1 cells of this sample, and their `_cells`
counts are counts within the sample. The other human pins (the assay, the genes, and the
colliding feature name) keep their v1.0 values.

- `cell_type_ontology_term_id` is the cell type that the most datasets hold, then the most cells.
- `tissue_general_ontology_term_id` follows the same rule.
- `dataset_id` is the control dataset.
- `disease_ontology_term_id` is one term of the multi-valued disease. Census 1 compares the
  whole string. `disease_ontology_term_id_matches` thus lists every joined value in the
  sample that holds the term, and the corpus filter selects those values. A stack with
  membership semantics takes the bare term.

No dataset of the 73 holds vasculature cells, so `documented_narrow` selects no cells on
either stack. The steps of the previous section produce the whole v1.1 snapshot, 74
results. To regenerate the human results only, add `--only census_data/homo_sapiens`.

## Scope limits

- The top-level SDK accessors do not reach the spatial experiments. A spatial
  `get_anndata` query runs through `experiment.axis_query` on the same handle, and a
  spatial `obs`, `var`, or presence read goes to the SOMA object directly. Each
  `manifest.json` records which path ran.
- An obs-filtered query is clamped to the manifest's datasets. A full-corpus X read is
  the benchmark rig's subject (task 2.5).
- The deferred corpus entries (identity coordinates, spatial scenes) have no oracle
  serialization. They wait on decision 7 and Phase 3.
- The filters are TileDB-SOMA `value_filter` strings. The Census 2 SDK takes a mapping
  of indexed column to values. Its adapter answers the entries it can express and
  reports the rest as unsupported.
