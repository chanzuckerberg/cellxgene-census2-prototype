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
| `sdk_adapter.py` | 6.2 | Answers a bound query through the Census 2 SDK, or names why the SDK cannot |
| `compare_stacks.py` | 6.2 | Runs the corpus through the SDK, compares with the recorded Census 1 answers, and names the task that owns each finding |

Not ported: the fixture reader over converted TileDB data (`zarr_adapter.py`,
`make_fixtures.py`, and the fixture driver of the same name as `compare_stacks.py`), the
benchmark rig (2.5), and the catalog prototype. The writer's stores and the SDK replace the
first group. The rig returns in Iteration 3.

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

The SDK run (task 6.2) compares Census 1 with 103 stores. `bulk-ingest-nf/manifests/`
lists their datasets: 100 human datasets sampled from release 2025-11-08, and the three
datasets of this sample. The three datasets let the silent-failure entries select cells
on both stacks:

| Dataset | Cells | Reason |
| --- | --- | --- |
| `5b25a46f-f1cf-4caa-9e10-6ab81dfe9174` | 12,822 | The only one of the 103 with a multi-valued disease. All its cells are primary data. |
| `6e00ccf7-0749-46ef-a999-dba785630d52` | 5,499 | Primary data with brain or lung cells |
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

No dataset of the sample holds vasculature cells, so `documented_narrow` selects no cells
on either stack. The steps of the previous section produce the whole v1.1 snapshot, 74
results. To regenerate the human results only, add `--only census_data/homo_sapiens`.

## Run the SDK against Census 1

`compare_stacks.py` reads the catalog of the stores under test and the oracle snapshot.
`obstore` reads S3 credentials from the environment. From `census2-sdk/tools/`:

```sh
eval "$(aws configure export-credentials --profile single-cell-dev --format env)"
export AWS_REGION=us-west-2
../.venv/bin/python -m harness.compare_stacks CATALOG --out results/sdk_vs_census1.json
```

The driver prints one line for each query:

| Status | Meaning |
| --- | --- |
| `pass` | No finding, other than the genes that only the SDK holds |
| `fail` | At least one finding. The first five follow the line. |
| `vacuous` | Neither side holds a cell. The findings stay on record. |
| `unsupported` | The SDK cannot express the query. The line names the reason. |
| `skip` | The catalog holds no store of the experiment's clamp sample |

The adapter translates equality and membership clauses into the SDK's filter mapping and
adds the clamp as a `dataset_id` filter. The SDK has no var filter, so a `feature_id`
filter selects genes from the result after the read. A numeric or `!=` clause, a column
that the catalog does not index, and the `get_obs`, `get_var`, and presence accessors are
unsupported.

The SDK's `var` is the whole pinned reference, and the comparator compares no X value
while the gene sets differ. The driver thus restricts both sides to the shared genes. The
genes that only the SDK holds are an expected finding, and a Census 1 gene that the SDK
does not hold fails the query. Each finding names the task that owns it (the `OWNERS`
table in `compare_stacks.py`).

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
