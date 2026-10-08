# Getting the Census onto AnnData-Zarr: a plan

**Status:** draft, for discussion. **Last edited:** 2026-10-08.

A delivery order in iterations is under discussion:
[census2-iterative-plan-1pager](https://app.notion.com/p/68916a75378982a3908e018e928d7eb7). It keeps the design
record and the deliverables and changes the order of the phases below.

This plan says in what order the Census 2 design gets built and what must be true before each piece starts.
The design is in the [design record](./census_zarr_design_v3_ste100.md) and its measurements are in the
[evidence](./census_zarr_design_evidence.md). "Design N.N" below names a section of the design record. The
assumptions below are what the plan rests on. A and C are decided. The rest stay revisable, and each says what
happens if it is wrong. The `tools/` and `dev-docs/` paths in this plan name files in the
[cellxgene-census](https://github.com/chanzuckerberg/cellxgene-census) repository or on its spike branches. They
are not in this repository.

The plan has no dates and no effort estimates. The phases are ordered by dependency, and several run alongside
each other where nothing forces them apart.

---

## Assumptions this plan rests on

Each assumption states what happens if it is wrong, so the cost of a revision is visible.

**A. Datasets are published individually and composed through a catalog.** Decided (design summary, [5.1](./census_zarr_design_v3_ste100.md#51-one-store-for-each-organism-and-modality)). The
end of the central build and private access both require it. A single store cannot provide either.

**B. Compaction is additive and deferred.** Everything is built as though per-dataset stores are the permanent
arrangement. Private data stays in per-dataset stores, so the uncompacted read path must exist ([design 6.8](./census_zarr_design_v3_ste100.md#68-compaction-and-other-materializations)).
*If wrong:* nothing built earlier is wasted. Compaction becomes Phase 14 sooner.

**C. `var` publishes the full pinned reference**, GENCODE v48 for human and M37 for mouse. Decided ([design
1.5](./census_zarr_design_v3_ste100.md#15-gene-reference)). It makes a gene's position independent of corpus
membership, so the writer writes each dataset without reference to the others. Assumption A depends on it.

**D. A release is a manifest, not a rebuild.** A release pins dataset versions and writer versions with
content digests. It sits next to a catalog snapshot, its dictionary, the presence matrices, and the
corpus-wide aggregates ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)). Data objects are
immutable. A fix lands at a new `writer_version`. The five-year commitment holds once retention on every
referenced version is agreed with the bucket owners, which gates Phase 10.
*If wrong:* Phase 10 needs a different answer for how citation and the support commitment survive.

**E. Harmonization does not gate publication.** A dataset is queryable once the writer harmonizes it at
ingestion. Only the per-gene `nnz`, `n_measured_obs`, and the per-category counts need every dataset. The
release builder derives them from files next to each store ([design
3.4](./census_zarr_design_v3_ste100.md#34-derived-per-dataset-files)).
*If wrong:* Phase 9 becomes a gate again and ingestion latency returns to build latency, which is the
problem this direction exists to solve.

**F. Compaction, if it happens, produces the layout already measured**: chunked and sharded `obs` and the
per-chunk index ([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout)). The `X` encoding is open
([design 6.10](./census_zarr_design_v3_ste100.md#610-x-encoding)). The measurements in the evidence hold either
way.

---

## The shape of the work

Phase 0 is done. The layout is measured and the store organization is decided. Three phases then gate
everything else and are independent of each other, so they start together. Approvals come first because
nothing lands upstream without them. The validation harness comes first because every later phase is checked
with it, and it stops being possible once the old stack retires. Spatial representation comes first because
the writer cannot write spatial experiments without an answer.

After those the chain is: writer, catalog, SDK, the two ingestion phases, publish, migrate, retire. Backfill
and ongoing ingestion share the writer and run concurrently. Private datasets and R parity hang off the side
and gate nothing on the corpus path. Compaction hangs off the end.

| #   | Phase                        | Gated on                              |
| --- | ---------------------------- | ------------------------------------- |
| 0   | The design                   | Done                                  |
| 1   | Approvals and amendments     | Decisions in [design section 7](./census_zarr_design_v3_ste100.md#7-decisions)         |
| 2   | The validation harness       | Nothing                               |
| 3   | Spatial representation       | Nothing                               |
| 4   | The per-dataset writer       | 1, and 3 for spatial                  |
| 5   | The catalog                  | Nothing to design, 4 to populate      |
| 6   | The v2 Python SDK            | 4, 5                                  |
| 7   | The identity proof           | 6, and the mapping (decision 7)       |
| 8   | Backfill ingestion           | 4 to write, 6 to verify               |
| 9   | Ongoing ingestion            | 4, 5                                  |
| 10  | Publish                      | 1, 5, 7, 8, 9                         |
| 11  | Private datasets             | 5, 6                                  |
| 12  | Migrate the ecosystem        | 10                                    |
| 13  | R parity                     | Bioconductor tooling, 10              |
| 14  | Compaction (optional)        | Evidence it is needed (design 6.1)    |
| 15  | Cutover                      | Everything, and the support window    |

---

## Phase 0: Design

**Delivers** the design record and the measurements behind it.

**Done.** The layout is decided and measured, compaction reproduces it unchanged (assumption F), and the store
organization is per dataset ([design
5.1](./census_zarr_design_v3_ste100.md#51-one-store-for-each-organism-and-modality)).

| Settled                                                           | Where                                                                                   |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| The read path's cost problem, measured against the current format | [Evidence section 9](./census_zarr_design_evidence.md#9-read-path-feasibility-evidence) |
| The layout within a store, chosen against measured alternatives   | [Evidence section 10](./census_zarr_design_evidence.md#10-layout-comparison)            |
| How cells and genes are identified                                | [Design 1.6](./census_zarr_design_v3_ste100.md#16-cell-and-gene-identity)                                                                              |
| The SDK surface                                                   | [Design 4.1](./census_zarr_design_v3_ste100.md#41-api), [4.2](./census_zarr_design_v3_ste100.md#42-selection-and-filter-interface)                                                                         |
| The gene space is pinned upstream, and its cost measured          | [Design 1.5](./census_zarr_design_v3_ste100.md#15-gene-reference), [evidence section 10.6](./census_zarr_design_evidence.md)                   |
| One store for each dataset                                        | [Design 5.1](./census_zarr_design_v3_ste100.md#51-one-store-for-each-organism-and-modality)                                                                              |

**Watch for** what the design does not cover. Spatial representation is Phase 3. Every corpus-scale figure is
an extrapolation until Phase 8 replaces it with a measurement.

---

## Phase 1: Approvals and amendments

**Delivers** Census schema 3.0.0, an amended storage and release policy, and a documentation site that no
longer defines a release as a TileDB-SOMA collection.

**Needs first** the decisions in [design section 7](./census_zarr_design_v3_ste100.md#7-decisions). Decisions
1 and 2 are the blocking pair. Decision 6 is this phase's own scope. Decision 5 is complete: review confirmed
on 2026-09-09 that `soma_joinid` is per-build, and this phase records the formal approval.

**Done when** all three documents are merged upstream, and the policy amendment carries a check that fails a
release which does not satisfy it.

**Changed by this plan.** The schema defines seven things instead of one (decision 6):

- the catalog and its dictionary
- the manifest (assumption D)
- the writer summary record
- the inventory
- the per-chunk index
- the per-dataset aggregate files
- the artifact address from `dataset_version_id` and `writer_version` ([design 1.2](./census_zarr_design_v3_ste100.md#12-store-location-and-address))

The release policy describes a release that is a pinned list and not a folder of built output.

**Watch for** two things. The approver of the release policy is unidentified. The five-year commitment lives
on the documentation site and not in the policy, so an amendment to the policy alone does not touch it.

---

## Phase 2: The validation harness

**Delivers** three things. A differential checker runs a query through both stacks and compares the results.
The query corpus is fixed. A benchmark rig reports data volume and request count as its reproducible metrics.
It adds controlled elapsed-time and concurrency runs where a phase calls for them (Phase 6 does, for the wide
filters).

**Needs first** nothing. It starts first, because everything after it is checked with it.

**Done when** the harness runs both stacks over a set of datasets and agrees on cells, genes, values, and
metadata columns. A deliberate corruption of a result must make it fail. A check that cannot fail is not a
check.

**Watch for** the oracle expiring. The Census 1 SDK over TileDB is the oracle, so every differential proof
must run before the TileDB corpus retires. The method is in [evidence section
11](./census_zarr_design_evidence.md#11-validation-methodology).

---

## Phase 3: Spatial representation

**Delivers** a decision on how scenes, point clouds, images, and the cell-to-scene mapping are stored, and the
accessor that reads them ([design 6.3](./census_zarr_design_v3_ste100.md#63-spatial-representation)).

**Needs first** nothing. It blocks the writer for spatial datasets.

**Done when** the mapping to SpatialData is in the schema, or spatial is deferred and the corpus ships
non-spatial first.

Per-dataset publication makes a deferral cheap. A spatial dataset is not ingested yet, and no store is
missing.

---

## Phase 4: The per-dataset writer

**Delivers** the census dataset writer: a library that takes one curated dataset and writes one immutable store. The
store has chunked and sharded `obs` and the per-chunk index
([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout)). The `X` encoding is open ([design
6.10](./census_zarr_design_v3_ste100.md#610-x-encoding)), and writer version 1 writes anndata's default CSR. It also has the assayed gene set and the
per-cell statistics. The writer encodes genes against the pinned reference (assumption C) and writes to the
address derived from `dataset_version_id` and `writer_version` ([design
1.2](./census_zarr_design_v3_ste100.md#12-store-location-and-address)). It leaves two records at the store
root ([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout),
[2.3](./census_zarr_design_v3_ste100.md#23-discover-metadata)). The inventory lists the path, size, and
checksum of every object. The summary record is what the catalog builder reads. The writer also writes the
per-dataset aggregate files next to the store in the asset bucket ([design
3.4](./census_zarr_design_v3_ste100.md#34-derived-per-dataset-files)). These are the per-gene `nnz` and the
per-category counts. The catalog builder and the release builder are the only writers under the census prefix.

The library has two callers ([design 3.6](./census_zarr_design_v3_ste100.md#36-writer-operation)): the
curation pipeline step (Phase 9) and the reconciliation job (Phases 8 and 9). Both must write byte-identical
stores. The library receives three staged inputs and fetches nothing at run time ([design
3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). The inputs are the dataset, the reference data
(gene reference, tissue roll-up), and the policy configuration. The configuration holds the experiments, the
`tissue_type` rule, the full-gene assays, the declared schema version, the `uns` allowlist, and the blocklist.
The summary record carries the version or digest of each input.

**Needs first** Phase 1, because the schema defines where the index lives before the writer can emit one. Also
Phase 3 for spatial datasets. Non-spatial datasets proceed without it.

**Done when** three checks hold. A dataset written this way is value-identical to the same dataset read
through the Census 1 stack. The writer records the gene reference version it used. It asserts `(dataset_id,
observation_joinid)` uniqueness on every dataset it writes. Also done only when the [disposition
table](./census_zarr_writer_dispositions.md) is exhaustive against the builder source. Every transformation in
`tools/cellxgene_census_builder` is then implemented by the writer, moved to the catalog builder or the release builder, or dropped
with its reason.

**Simpler than it was.** The previous plan required a writer that streamed an entire organism because whole
matrices did not fit in memory. One dataset at a time removes that constraint. The largest dataset is a known,
bounded quantity.

**State.** Task 4.1, the disposition table, is done. Task 4.2, the gene-reference projection, is in progress.
Task 4.10, the pipeline entry point, is in review as PR #2 on `cellxgene-census2-prototype` (Mike). The PR
holds the `census2-writer` package with `create()`, `build()`, and `write()`, and the Nextflow reconciliation
job. Task 4.11 defines the library contract for both callers. Its five items are:

- address derivation
- atomic finalize, with the inventory as the last object
- the summary record typed in the package
- aggregate files next to the store
- input resolution

The review of PR #2 asks for that contract.

**Watch for** the gene reference version recorded per dataset and not assumed. A dataset curated against a
different GENCODE pin than the corpus is a real condition. The writer refuses it and names both pins. It does
not remap or drop the dataset.

Two findings from the 2026-08-24 local staging runs (`tools/backfill_inventory/`):

- **Contributor `uns` carries dtypes with no stable Zarr v3 encoding.** One dataset in a 25-dataset sample holds scanpy `rank_genes_groups` output in `uns`. The output has numpy record arrays and fixed-length UTF-32 strings, which zarr-python flags as unstable and Rarr cannot read. The writer keeps only the `uns` keys that the
  curation schema defines ([design 1.3](./census_zarr_design_v3_ste100.md#13-store-contents), disposition table).
- **The largest dataset fits a normal worker when the writer streams.** The bound on "a known, bounded
  quantity" is a 53.2 GB source H5AD. The writer in PR #2 reads the matrix as dask row blocks. A 100 dataset
  run on HealthOmics measured the sparse peak at about 0.90 GiB plus 0.16 times the file size. The largest
  dataset then peaks near 10 GiB. Dense matrices were the only cause of memory failure, and admission now
  excludes them on the assay allowlist. Conversion compute is not the constraint: the staging runs measured
  about 11 seconds per source GB on a laptop.

---

## Phase 5: The catalog

**Delivers** the catalog builder and the SDK's exclusion path ([design
2.1](./census_zarr_design_v3_ste100.md#21-catalog-fields),
[2.2](./census_zarr_design_v3_ste100.md#22-dataset-selection)). The catalog has one row per dataset. A row
carries identity, partition (organism and spatial), conformance, the per-column value bitmaps, and the numeric
ranges. The global dictionary sits next to it. The builder reads the Discover API for the metadata that
Discover owns ([design 2.3](./census_zarr_design_v3_ste100.md#23-discover-metadata)). It reads the summary
record and the inventory of each store for the rest. It reads no array data.

**Needs first** nothing to design. The dictionary and the bitmap schema come from the existing TileDB corpus,
where the categories and their cardinalities are known. The conformance and summary-record fields need real
stores, so those wait on Phase 4.

**Done when** a predicate resolved against the catalog excludes no dataset that a full scan selects, proven by
running both over the harness queries. Over-selection is allowed. An opened store that yields no rows costs
data transfer and not correctness, and the test measures it. The catalog for the full corpus stays under about
3 MB.

**State.** A first catalog construction from the TileDB corpus passed 55 of 55 harness queries with no false
negative (evidence section 17).

**This is the critical path piece.** Without the catalog, a query opens every store and pays about 50 kB each.
That is 54 MB across the human corpus against 3.4 MB for a whole documented query today. With it, the
documented filters cost 0.5 MB to 22 MB depending on how many datasets they touch. The catalog is what makes
assumption A viable.

**Watch for** the column set being treated as fixed. It follows the curation schema ([design
2.3](./census_zarr_design_v3_ste100.md#23-discover-metadata)). Schema 7.1.0's perturbation fields, which
postdate the builder's 7.0.0 pin, are the first addition the catalog and the writer's projection pick up. The
catalog and the per-chunk `obs` index are the same structure at two scales, so one implementation serves both.
`donor_id` is not indexed, at 11,633 distinct values, so a predicate on it reads the column. The publication
protocol is in [design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents). One builder makes one
serialized atomic write, and a release pins the dictionary and catalog it was encoded against.

---

## Phase 6: The v2 Python SDK

**Delivers** the surface in [design 4.1](./census_zarr_design_v3_ste100.md#41-api) and 4.2.
`open_census()` opens a set of sources, catalogs and release manifests. The SDK keeps the accessors and adds
the new filter interface, selection by the new identity, and lazy access. The SDK refuses a TileDB-SOMA filter
string and puts the translated equivalent in the error message and in the migration guide ([design
4.2](./census_zarr_design_v3_ste100.md#42-selection-and-filter-interface), decided 2026-09-08). The grammar
itself is open ([design section 8](./census_zarr_design_v3_ste100.md#8-required-design-details)).

**Needs first** Phases 4 and 5.

**Done when** three things hold. Every query in the harness corpus returns identical results through the old
stack and the new one. A query over two catalogs returns the union of both. The Phase 2 rig has measured the
harness corpus end-to-end over real per-dataset stores. It reports data volume and request count against the
same queries on the old stack. [Design
5.1](./census_zarr_design_v3_ste100.md#51-one-store-for-each-organism-and-modality) measures opening cost
only. For the wide filters this includes controlled elapsed-time and maximum-concurrency runs. To open
hundreds of stores at once is a latency, throttling, and partial-failure question that data volume and request
count do not show. The same rig compares the three exclusion arms of the per-cell `obs` frame item ([design
6.2](./census_zarr_design_v3_ste100.md#62-per-cell-obs-frame)). The arms are the bitmaps plus the in-store
per-chunk index, per-chunk summaries published with the catalog, and per-dataset Parquet `obs` fragments read
as one frame. Apache Iceberg is the candidate packaging for the third. One input is not taken: how many
datasets pass a two-predicate filter at dataset grain and return no cells. Two decisions follow from this
measurement. Phase 4 decides whether the writer emits the fragment. Phase 1, decision 6, decides whether the
per-chunk index stays in the schema.

**Build on anndata's lazy stack first.** anndata reads Zarr lazily and concatenates lazily, some of it over
dask with lighter paths available. Implement the read path on that, measure it in the Phase 2 rig, and write
custom reader code only where the measurement falls short. The pinned reference removes the cost that this
evaluation hinges on. Within an organism the gene axes match by construction, so concatenation is by position
and does not join gene lists across stores ([design 1.5](./census_zarr_design_v3_ste100.md#15-gene-reference),
[4.3](./census_zarr_design_v3_ste100.md#43-additive-source-composition)).

**Changed by this plan.** `open_census()` becomes a call over a source set. That is the visible face of
additive queries, and it is what lets private data participate in Phase 11.

**Watch for** the new grammar giving `disease` and `disease_ontology_term_id` explicit membership syntax and
refusing bare equality on them ([design
4.2](./census_zarr_design_v3_ste100.md#42-selection-and-filter-interface)). The old translator behavior fails
silently if missed, and it now lives only in the error-path translator and the harness. The `experimental`
namespace is unscoped, in particular `get_embedding`, which is keyed on the identifier being retired. The
supply side of embeddings is an open item ([design 6.7](./census_zarr_design_v3_ste100.md#67-embeddings)). Two
SDK details come from the use-case review. A result must carry the release or source identity it came from,
with a warning on an unpinned live-catalog source. Without that, reproducibility depends on user discipline
alone. The mapping from `(dataset_id, observation_joinid)` to `obs_names` must be unique and lossless, because
AnnData and Seurat tooling key on cell names. One result-schema case was measured during Phase 2. Census 1's
dtype for nullable integer columns (the spatial `array_row`, `array_col`, `in_tissue`) depends on whether the
selection includes nulls (evidence section 18.1). The SDK picks a stable representation and states it.

---

## Phase 7: The identity proof

**Delivers** a real contributed embedding, mapped through the published mapping and re-joined against the new
corpus, on the same biological cells.

**Needs first** the identity mapping and the SDK from Phase 6. The mapping does not come from the writer. With
no full-scale build to derive it from, its source is the last TileDB release (decision 7), and Phase 10
derives it.

**Done when** the four assertions in [evidence section
12](./census_zarr_design_evidence.md#12-identity-transition-proof) hold, each paired with a corruption that
makes it fail.

**More important under this plan.** `(dataset_id, observation_joinid)` stops being only a replacement for a
row number and becomes the key that makes results from different sources combinable ([design
1.6](./census_zarr_design_v3_ste100.md#16-cell-and-gene-identity)). The uniqueness assertion, verified on 3.7
percent of datasets so far, moves into the writer in Phase 4 (task 4.6). Every ingestion then checks it. What
remains here is the end-to-end re-join against real external work.

**Watch for** this being a release gate, demonstrated against a real embedding and not a synthetic one.

---

## Phase 8: Backfill ingestion

**Delivers** every existing dataset that passes the builder's inclusion gate, written as an individual store,
with a catalog over them. The queue is live and not fixed. The published Census holds 1,845 datasets, which
counts the last build (2025-11-17). On 2026-08-24 the same gate over the Discover corpus admitted 2,105 of
2,216 public datasets, 2.93 TB of H5AD, 660 of them spatial (`tools/backfill_inventory/`). Size the backfill
against a fresh inventory and not either snapshot. The assay gate stays for Census 2, decided 2026-08-24, at
dataset granularity for admission. Inside an admitted dataset the writer applies the builder's cell-level
filters unchanged (disposition table, row 11), so the Census 1 corpus semantics hold. A wider gate over the
full Discover corpus stays a later option. The delta is 111 datasets and 0.44 TB, gated on gene-reference pins
and non-RNA conformance and not on cost. Queue semantics and the failure ledger are in
`dev-docs/census_zarr_backfill_queue.md`.

**Needs first** Phase 4 to write. Verification needs Phase 6, because the differential harness reads through
the SDK. Backfill starts before the SDK is finished and is signed off after it.

**Done when**, in two steps tracked separately:

1. *Written.* Every dataset has a store and a catalog row, and every dataset that failed is listed with a
   reason.
2. *Verified.* The differential harness passes over a sample that spans every organism, and the
   corpus-scale figures replace the extrapolations the design carries.

**Replaces the full-scale build.** There is no single build to run. Backfill is a queue of independent
per-dataset jobs, run by the reconciliation job (Phase 9, tracker 9.1). The queue is parallel, resumable, and
pausable without a broken state. A failure affects one dataset and not the corpus, which is the property the
stalled build does not have.

**Watch for** the store layout. Writer version 1 writes plain CSR and default chunks. The `obs` layout in
[design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout) lands with task 4.4, and the `X` encoding waits
on [design 6.10](./census_zarr_design_v3_ste100.md#610-x-encoding) and task 4.3. Every store written before
them is rewritten at a new `writer_version`. Keep the backfill small until the layout settles. Also watch object count, but do not carry the prototype's figure across. Its stores were unsharded
Zarr v2 and averaged about 1,850 objects each. This layout shards, and no spike store exceeded about 700
objects at any scale tested. Record the real number here.

**Watch for** the gate count. Admission in PR #2 passes 2,175 of 2,228 datasets on the assay allowlist. The
gate above passed 2,105 of 2,216. The catalog builder holds one definition, and the count is re-measured when
the backfill starts.

**Watch for** the failure modes of the 2026-09-21 run of 100 datasets. A task that raises leaves a partial
store at its final address. Ten such tasks stopped the report step, so the run published no report. Task 4.8
removes the partial store with the version in the path and the inventory as the last object. The report step
must collect the receipts that exist (tracker 9.1).

---

## Phase 9: Ongoing ingestion

**Delivers** the standing path: curation publishes a dataset, the writer writes it, its catalog row appears,
and it becomes queryable without a build.

**Needs first** Phases 4 and 5. It shares the writer with Phase 8 and is built alongside it. Backfill is the
same machinery pointed at a finite queue.

**Done when** three things hold. A dataset that lands upstream is queryable without manual intervention. A
failed ingestion is visible and retryable without blocking others. The catalog update is atomic, so that no
reader sees a dataset whose store is incomplete.

**This is what unstalls releases**, the first of the two reasons this direction exists. Without it,
per-dataset publication is a different shape of the same manual build.

**Operation ([design 3.6](./census_zarr_design_v3_ste100.md#36-writer-operation), decision 8).** One writer
library, two callers.

- The curation pipeline gains a step that writes the store at ingestion and registers the Zarr artifact.
  The Census 2 team owns this step (decided 2026-09-14). It lands on Step Functions, the portal's
  orchestrator (tracker 9.2).
- A Nextflow reconciliation job, on AWS HealthOmics for the prototype with the production runtime open, scheduled and on demand, calls the same library. It covers backfill (Phase 8), a new `writer_version` after a writer or codec fix, and repair after an ingestion failure. Each run compares the published dataset versions with the catalog rows (tracker 9.1, PR #2).

Both callers end at the catalog builder's one serialized atomic write. A HealthOmics task has no internet
access by default, so the writer reads staged inputs in both runtimes ([design
3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). Three items are open under [design
3.6](./census_zarr_design_v3_ste100.md#36-writer-operation). How does the reconciliation job register an
artifact? What allows only one catalog write at a time? Who owns the continuing
operation (decision 8)?

**Watch for** the catalog as a write-contention point, because the builder rewrites it after every ingestion.
Also watch for conformance drift. The writer records a dataset that arrives under a newer curation schema or
gene reference than the corpus declares. It holds the dataset with a reason (disposition table, row 9). The
reference-transition rule lands here too ([design
3.2](./census_zarr_design_v3_ste100.md#32-gene-reference-update)). The catalog holds the old reference until
every applicable dataset converts. It then changes the organism in one atomic operation, and a dataset that
arrives during the migration enters at the change.

---

## Phase 10: Publish

**Delivers** the first release ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents),
[3.5](./census_zarr_design_v3_ste100.md#35-release-discovery)). The release builder assembles it. The release has:

- a manifest that pins dataset versions and writer versions with content digests
- the catalog and dictionary snapshots
- the presence matrices and aggregates next to it
- the entry in the release index
- the identity mapping

**Needs first** Phases 1, 5, 7, 8, and 9, and the retention agreement with the bucket owners ([design
6.4](./census_zarr_design_v3_ste100.md#64-storage-locations-and-retention)).

**Done when** the SDK reads a published manifest from S3 and resolves the right dataset versions. Each store's
inventory verifies against the manifest digest, and the mapping resolves for a real embedding.

**Changed by this plan.** A release is a manifest and not a folder of built output (assumption D). That is
what lets releases be cut cheaply and often while the five-year commitment attaches to immutable objects. The
release contents are in [design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents). The retention
agreement is part of the storage conversation in [design
6.4](./census_zarr_design_v3_ste100.md#64-storage-locations-and-retention).

**Watch for** this being the point where the move becomes mandatory, so the timeline is announced before the
release and not alongside it. The identity mapping comes from the last TileDB release.

---

## Phase 11: Private datasets

**Delivers** the second reason this direction exists: datasets outside the public Census, queryable in the
same call as those inside it ([design 4.5](./census_zarr_design_v3_ste100.md#45-private-data)).

**Needs first** Phases 5 and 6. It gates nothing on the corpus path and sits here for that reason, not because
it depends on publication.

**Done when** three things hold. A caller queries a public catalog and a private one together and receives the
union. A caller without rights to the private store gets a clear failure and not silently missing rows. A
private dataset with a category absent from the public dictionary resolves through a dictionary extension.

**Delivers specifically:** the access-control model, how a private catalog is produced and where it lives,
credential handling in the SDK, and the dictionary-extension mechanism. The starting point is the defaults in
[design 4.5](./census_zarr_design_v3_ste100.md#45-private-data). One grant covers a catalog and its stores. An
access-controlled file lists private catalog URIs. The SDK resolves credentials for each source. One
inaccessible source fails the whole query and names the source. Displace a default only by decision ([design
6.5](./census_zarr_design_v3_ste100.md#65-private-discovery-and-errors)).

**Watch for** access control staying in storage. The catalog describes and storage enforces. A catalog that
carries permissions discloses them.

Also watch for results that mix conformance levels. A private dataset can lack columns that a Census dataset
carries. So either results carry per-row provenance, or the query declares what it requires and excludes the
rest ([design 6.6](./census_zarr_design_v3_ste100.md#66-mixed-conformance)).

---

## Phase 12: Migrate the ecosystem

**Delivers** notebooks, the documentation site, and tutorials that run against the new corpus, and a published
migration path for external users.

**Needs first** Phase 10.

**Done when** the notebooks and vignettes execute green, and canary users outside the team have migrated.

**Watch for** the conceptual change being larger than the API change. There is no corpus object to open, and
the documentation teaches source sets instead of a handle.

---

## Phase 13: R parity

**Delivers** an R client that reads the published corpus.

**Needs first** a working sharded read in the R stack, and Phase 10 for something to read. The early check ran
on 2026-08-18 and failed. anndataR over the development Rarr reads an unsharded store of a real dataset
correctly. It returns incorrect data for the same store sharded ([design
4.6](./census_zarr_design_v3_ste100.md#46-existing-user-support)). So the gate is more than the Bioconductor
release around October 2026. An upstream fix must land before that release cuts from the development branch,
or this phase slips past it.

**Done when** R and Python agree on the harness query corpus.

**Watch for** the interim options if the upstream fix misses the release: a reticulate bridge over the Python
SDK, or a documented query-in-Python-then-`write_h5ad` workflow. Both avoid sharded reads. Until R can read
Zarr, R users see no new data. The catalog is Parquet, which R reads without new dependencies, so the Zarr
side gates.

---

## Phase 14: Compaction (optional)

**Delivers** a background process that folds settled datasets into per-organism stores, and the catalog fields
that point into them ([design 6.8](./census_zarr_design_v3_ste100.md#68-compaction-and-other-materializations)).
A gene-major (CSC) copy of `X` is a per-dataset materialization. It follows the `X` encoding decision ([design
6.10](./census_zarr_design_v3_ste100.md#610-x-encoding), task 4.3) and does not wait for this phase.

**Needs first** evidence that it is needed. The two measurements of [design
6.1](./census_zarr_design_v3_ste100.md#61-end-to-end-measurement) supply it: the fan-out distribution of real
queries, and the open cost as a function of store count. The team names an elapsed-time limit for an
interactive query. This phase starts when common queries exceed that limit. Neither measurement needs a full
ingest, so this decision can come before Phase 8.

**Done when** a compacted store returns results identical to the datasets it replaced. A query over hundreds
of datasets then reads one store instead of hundreds.

**Optional, and gates nothing on the corpus path.** Assumption B holds it here. The catalog brings the
narrower documented filters below what they cost today. Compaction is worth about 400x on metadata for the
widest filter and nothing for the narrow ones. Private data stays uncompacted, so the SDK supports the
per-dataset path in every case.

**Watch for** the trigger being wide-query cost, measured and not assumed. Object count looked like a second
trigger and does not survive inspection: that figure came from unsharded prototype stores, and this layout
shards. Also watch for continuous ingestion. Datasets arrive after a materialization exists, so this phase must
define an append to a large sharded store or a rebuild.

---

## Phase 15: Cutover

**Delivers** retirement of the TileDB corpus and the Census 1 SDK.

**Needs first** everything above, and the support window.

**Done when** the support matrix is retired per the availability commitment.

**Watch for** this being policy-gated and not readiness-gated, and constrained by R: a TileDB retirement
before R can read Zarr removes access for those users. The support matrix this retires is also weaker than it
looks. The nightly LTS compatibility workflow covers Python only, and as of 2026-08-19 it had failed for over
a month on its old `tiledbsoma` pin cells. Re-establish what the commitment holds before a cutover date is
set.

---

## Running alongside

These are not phases and gate nothing, but they need an owner.

- **Keeping the design record current.** The design record is
  [census_zarr_design_v3_ste100.md](./census_zarr_design_v3_ste100.md). The earlier v3, v2, and v1 are
  superseded; the old v3 page hosts the open review threads. The catalog write-up is archived at
  `dev-docs/archive/census_zarr_catalog.md`.
- **Scope for the `experimental` namespace**, which is SOMA-coupled and outside the design.
- **Embeddings as an ongoing product** ([design 6.7](./census_zarr_design_v3_ste100.md#67-embeddings), evidence section 16). If a posture beyond downstream-consumers'-work is chosen, three pieces follow. The one-time conversion of the hosted tables (about 1.36 TB, newest per model) hangs off Phase 7's mapping. Ingestion-time inference under a frozen model hangs off Phase 9. The per-release rerun that produced them was manual and has no owner.
- **Upstream documentation defects** found while measuring, listed in
  [evidence section 14](./census_zarr_design_evidence.md#14-findings-that-need-a-decision-elsewhere). They
  are cheap to fix and build goodwill with the team whose approval Phase 1 needs.
- **The dead links in the published copies.** The Notion spec pointed at a deleted branch until
  2026-08-05. Anything else that references that branch is still broken.

---

## Open questions for this plan

- **Are the remaining assumptions accepted?** A and C are decided. B, D, E, and F are open.
- Can Phases 4, 5, and 6 overlap, with the catalog and SDK built against a set of stores? Or does the writer need to finish first?
- Does Phase 11 move earlier? It sits after publication because it gates nothing there. But it is one of the two reasons for this direction, and late placement is how such things get dropped.
- Is spatial in the first published corpus, or deferred so that non-spatial ships sooner? Per-dataset
  publication makes a deferral cheap.
- Does backfill need infrastructure that must be requested before Phase 8 starts?
- Does a full-scale ingest run before tasks 4.3 and 4.4 land in the writer? The run measures the open cost
  and the read cost on real stores, after the measurements of design 6.1. The writer rewrites these stores at
  a new `writer_version` after those tasks.
- Who owns each phase? The Census 2 team owns the curation pipeline step (9.2) and Mike owns the
  reconciliation job (9.1). The continuing ingestion operation and every other phase are unassigned.
