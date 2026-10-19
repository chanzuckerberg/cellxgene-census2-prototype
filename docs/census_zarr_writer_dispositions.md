# The harmonization disposition table (task 4.1)

**Status:** draft, for review. **Last edited:** 2026-09-23.

This table assigns one disposition to each transformation of the Census 1 builder.
[Evidence section 15](https://app.notion.com/p/3b216a75378981a8b76fe2955c600c4d) lists the
39 steps and the validation suite, read from `tools/cellxgene_census_builder` in the
[cellxgene-census](https://github.com/chanzuckerberg/cellxgene-census) repository. Step numbers and section numbers
follow evidence section 15. Design section numbers refer to the
[design record](./census_zarr_design_v3_ste100.md).

The dispositions are:

- writer: the census dataset writer performs the step for one dataset at ingestion.
- caller: the curation pipeline step or the reconciliation job performs the step before it calls the writer
  ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)).
- catalog builder: the catalog builder performs the step when it updates the live catalog
  ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)).
- release builder: the release builder performs the step when it assembles a release
  ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)).
- replaced (task): another mechanism covers the step. The task that builds the mechanism is named.
- dropped: no consumer remains. The note names what breaks if that is wrong.

The writer rows define the scope of tasks 4.2 to 4.8.

The Iteration column says when a step lands under the proposed delivery order
([census2-iterative-plan-1pager](https://app.notion.com/p/68916a75378982a3908e018e928d7eb7)). Iteration 1 is the
end-to-end skeleton over the first store set, and it takes the steps that make stores concatenate into one
AnnData. Iteration 2 is full harmonization, and it takes the steps that make values match Census 1. Later is
Iteration 3 onward. A dropped step has none.

Two rules that this table sets:

- The writer keeps the Census 1 cell filters (step 11). These are `tissue_type` in tissue or organoid, the assays
  that the experiment permits, and `is_primary_data == True` for spatial data. Admission is per dataset. The writer admits
  or refuses a dataset whole, and it filters cells inside an admitted dataset exactly as the Census 1 builder
  does. Task 4.9 then compares equal cell sets.
- The writer keeps only the `uns` keys that the curation schema defines, such as `citation`, `title`, and the
  organism fields. It drops contributor analysis artifacts.

## 15.1 Admission

| # | Builder behavior | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 1 | Discover manifest + schema-version gate | writer + reconciliation job | 1 | New datasets arrive through the curation pipeline step. The reconciliation job compares published dataset versions with catalog rows and writes the missing stores ([design 3.6](./census_zarr_design_v3_ste100.md#36-writer-operation)). Schema conformance is writer admission, recorded as `curation_schema_version` and `census_eligible` ([design 2.1](./census_zarr_design_v3_ste100.md#21-catalog-fields)). The gate follows the curation schema version that the corpus declares ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). The writer records drift and does not ingest it (plan Phase 9). |
| 2 | Blocklist; duplicate `dataset_id` collapse | writer + catalog builder | 2 | The writer applies the reconciliation blocklist from the policy configuration ([design 1.3](./census_zarr_design_v3_ste100.md#13-store-contents), [3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)) and records the refusal as `census_eligible: false` with a reason. Duplicates cannot arise under `dataset_version_id` identity. The catalog builder asserts `dataset_id` uniqueness, and the SDK reports an error when a `dataset_id` occurs in more than one source ([design 4.4](./census_zarr_design_v3_ste100.md#44-query-execution)). |
| 3 | Filesize ordering of datasets | dropped | none | This step fixed the total order of `soma_joinid`. No corpus order exists and identity is intrinsic ([design 1.6](./census_zarr_design_v3_ste100.md#16-cell-and-gene-identity)), so no consumer remains. |
| 4 | Zero-cell dataset exclusion | writer | 2 | The writer refuses a dataset with zero cells after the filters and records the reason. It writes no store. |
| 5 | Dataset `soma_joinid` assignment | dropped | none | The old mapping is published once (decision 7). A wrong mapping breaks re-joins of old work. Phase 7 proves the mapping. |
| 6 | `census_info/datasets` table | catalog builder | 1 | The catalog's description fields carry it ([design 2.1](./census_zarr_design_v3_ste100.md#21-catalog-fields)). |

## 15.2 Reading and filtering

| # | Builder behavior | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 7 | Staging with size check | caller | 1 | The caller stages the dataset and checks its size against the manifest ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). The writer fetches nothing. |
| 8 | `raw/X` selection when a raw group exists | writer | 1 | The writer keeps this step. The corpus carries raw counts. |
| 9 | Per-file `uns/schema_version` and encoding assertion | writer | 2 | This step has three checks. The writer refuses the dataset when the file's `uns/schema_version` differs from the version that Discover declares for that dataset version. It also refuses an AnnData encoding version that it cannot parse. The third check records `curation_schema_version` ([design 2.1](./census_zarr_design_v3_ste100.md#21-catalog-fields)) and compares it with the schema version that the corpus declares in the policy configuration ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). The writer holds a dataset on a newer version and records the reason until the corpus declares that version (plan Phase 9). The live datasets are on schema 7.1.0. The Census 1 builder pins 7.0.0. |
| 10 | Duplicate-coordinate `sum_duplicates()` | writer | 2 | The writer keeps this step. |
| 11 | Cell filter (tissue_type, assays, spatial primary) | writer | 2, or 1 if the Iteration 1 gate compares against Census 1 | The writer applies the same criteria as the Census 1 builder, from the same allowlists ([design 1.3](./census_zarr_design_v3_ste100.md#13-store-contents), [3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). |
| 12 | Gene filter (`feature_biotype == "gene"`, organism match) | writer (4.2) | 1 | Projection onto the pinned reference covers this step. A gene outside the organism's reference does not map. |
| 13 | `obs` projection, dtype normalization, `dataset_id` stamp | writer (4.4) | 1 | The writer projects the column set of the curation schema version that the corpus declares. Schema 7.1.0 adds perturbation columns, which a frozen list omits. |
| 14 | `var` projection to the 4-column form | writer (4.2) | 1 | The pinned reference supplies the feature metadata. The dataset's `var` maps onto it. |

## 15.3 The corpus axes

| # | Builder behavior | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 15 | `obs` concatenation in dataset order | dropped | none | The SDK concatenates stores by position at read time ([design 4.3](./census_zarr_design_v3_ste100.md#43-additive-source-composition), [4.4](./census_zarr_design_v3_ste100.md#44-query-execution)). No corpus frame exists. |
| 16 | `var` first-seen union | replaced (4.2) | 1 | The pinned full reference replaces the union ([design 1.5](./census_zarr_design_v3_ste100.md#15-gene-reference)). Positions do not change when a dataset enters the corpus. |
| 17 | Cell and gene `soma_joinid` assignment | dropped | none | Identity is `(dataset_id, observation_joinid)` and `feature_id` ([design 1.6](./census_zarr_design_v3_ste100.md#16-cell-and-gene-identity)). Decision 7 publishes the old mapping. |
| 18 | `tissue_general` derivation | writer | 2 | Per cell through `cellxgene_ontology_guide`, and a catalog column ([design 2.1](./census_zarr_design_v3_ste100.md#21-catalog-fields)). The tissue roll-up is staged reference data ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)), with its version in the summary record. |
| 19 | Sentinel placeholders for missing columns | dropped | none | This step served the concatenation of a heterogeneous corpus. The writer computes statistics before the write and projects nullable dtypes, so a store has no sentinel values. |
| 20 | Joinid bookkeeping and the feature_id→joinid map | replaced (4.2) | 1 | Reference positions are the gene coordinates. No per-build map exists. |
| 21 | Contiguity assertion | dropped | none | No corpus axis exists. |
| 22 | Categorical casting and dtype assertion on write | writer (4.4) | 1 | This step becomes the store's dtype and encoding policy. The differential harness set two constraints: a stable representation for nullable integers (evidence section 18.1) and Zarr v3 stable dtypes only (evidence section 18.3). |
| 23 | Dataframe write with explicit TileDB domain | dropped | none | This step is TileDB-specific. The Zarr counterpart is the layout in [design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout): chunked and sharded `obs` and `obs_index` (tasks 4.3, 4.4). |

## 15.4 Expression and derived values

| # | Builder behavior | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 24 | Gene remap into corpus space; unmappable gene fails the build | writer (4.2) | 1 | The remap targets the pinned reference. The writer refuses a dataset on a different pin and names both pins. The refusal covers one dataset and not a build. |
| 25 | Single-modality assertion | writer | 2 | This check is admission. |
| 26 | X coercion (CSR, float32, non-negative-integer assertion) | writer (4.3) | 1 | The writer keeps this step. |
| 27 | Per-cell statistics (`raw_sum`, `nnz`, `raw_mean_nnz`, `raw_variance_nnz`) | writer | 2 | The statistics are per dataset. Their ranges go in the summary record ([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout), [2.3](./census_zarr_design_v3_ste100.md#23-discover-metadata)). |
| 28 | `n_measured_vars` | writer | 2 | `n_measured_vars` is the count of the dataset's assayed genes, from its filtered source `var`, and not the reference length. The meaning survives the longer `var`. |
| 29 | Per-gene `nnz` and `n_measured_obs` | writer + release builder | 2 | The writer writes the per-gene `nnz` for each dataset as a file next to the store, in the asset bucket ([design 3.4](./census_zarr_design_v3_ste100.md#34-derived-per-dataset-files)). A new aggregate definition rewrites these files and not the store. The release builder sums them at release. `n_measured_obs` is not stored. The release builder derives it from the assayed gene set and `n_obs`. |
| 30 | Stats merge into the corpus frame | dropped | none | The writer writes statistics into its own `obs`. No corpus frame and no sentinel exist. |
| 31 | `X["raw"]` write at corpus coordinates | writer (4.3) | 1 (plain CSR); the encoding is later (design 6.10) | The writer writes `X` at the reference gene coordinates. The store's rows are the dataset's rows, so no row shift occurs. |
| 32 | `X["normalized"]` derivation | writer | 2 | The writer applies the Census 1 algorithm: gene-length correction for the 21 full-gene assay terms ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)), division by the row sum, 15-mantissa-bit round-half-to-even, and a subnormal value for zero. The target is a bit-identical result. Task 4.9 compares the result under the harness tolerance. |
| 33 | Presence matrix assembly | writer + release builder | 2 | The writer records the assayed gene set for the dataset. `X` cannot recover the set, because a stored zero does not say whether the gene was measured ([design 1.3](./census_zarr_design_v3_ste100.md#13-store-contents)). The release builder builds the presence matrix for each organism from the recorded sets ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)). |
| 34 | Spatial extraction (scene, point cloud, image) | writer | later | The writer extracts spatial data after Phase 3 defines the representation ([design 6.3](./census_zarr_design_v3_ste100.md#63-spatial-representation)). The 660 spatial datasets wait until then. |

## 15.5 Corpus tables and metadata

| # | Builder behavior | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 35 | Per-category summary counts | release builder | 2 | The release builder sums the per-category count files next to each store at release ([design 3.4](./census_zarr_design_v3_ste100.md#34-derived-per-dataset-files)). A new definition rewrites the files and not the store. The release snapshot is the citable copy. |
| 36 | `census_info/summary` | release builder | 2 | This table becomes release metadata in the manifest ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)). |
| 37 | `census_info/organisms` | release builder | 2 | This table becomes release metadata in the manifest ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)). |
| 38 | Provenance (`created_on`, `git_commit_sha`) | writer + catalog builder | 1 | The summary record and the inventory carry the store's provenance ([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout)). The release manifest carries its own provenance and digests ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents)). |
| 39 | Bounding-box metadata scrub | dropped | none | This step was a TileDB-SOMA workaround. Nothing remains to scrub. |

## 15.6 The validation suite

| Check | Disposition | Notes |
| --- | --- | --- |
| Per-dataset passes: obs equals filtered source, `X["raw"]` equals source on assayed genes, stats recomputed, assayed set correct | writer | Self-checks at ingestion. A failed check fails the write. Task 4.9 runs the corpus-scale form through the differential harness. |
| Joinids unique and exactly `[0, n)` | replaced (4.6) | The writer asserts `(dataset_id, observation_joinid)` uniqueness on every write. |
| `var` union check | replaced (4.2) | Projection against the pinned reference is the check. |
| No duplicate coordinates in either X layer | writer | The writer checks each store at ingestion. |
| Normalization recomputed, row sums near 1, no stored zero | writer | Per-dataset self-check. The catalog builder checks a store only against its inventory ([design 3.3](./census_zarr_design_v3_ste100.md#33-release-contents), step 4). |
| `dataset_id`/`dataset_version_id` unique in the manifest | release builder | The release builder asserts uniqueness. The SDK reports an error when a `dataset_id` occurs in more than one source ([design 4.4](./census_zarr_design_v3_ste100.md#44-query-execution)). |

## New: contributor `uns`

| Input | Disposition | Iteration | Notes |
| --- | --- | --- | --- |
| Contributor `uns` (not among the 39 steps; the Census 1 corpus does not carry it) | writer: allowlist | 1 | The writer keeps the `uns` keys that the curation schema defines and drops contributor analysis artifacts ([design 1.3](./census_zarr_design_v3_ste100.md#13-store-contents)). The allowlist is policy configuration ([design 3.7](./census_zarr_design_v3_ste100.md#37-writer-inputs)). This removes the unstable-dtype exposure (evidence section 18.3: one dataset in a 25-dataset sample carried scanpy record arrays). The Census 1 corpus carries no contributor `uns`, so nothing regresses. |

## New: `obs` column names with a slash

| # | Input | Disposition | Iteration | Notes |
| --- | --- | --- | --- | --- |
| 40 | An `obs` column whose name holds a `/` (clinical units, compound cell type labels, cell2location factors; 6 of 100 sampled human datasets, 1 or 2 columns each) | writer: drop, prototype only | 1 (drop); the rename is later | Zarr gives `/` its path meaning, so the write fails on the column. The prototype writer (PR #2) drops the column and records the name in the receipt. The option when this moves past the prototype is a rename to a Zarr-safe name with the mapping recorded ([design 1.4](./census_zarr_design_v3_ste100.md#14-store-layout), task 4.4). |

## The writer rows, by task

| Task | Rows |
| --- | --- |
| Admission and the entry path | 1, 4, 8, 9, 10, 11, 25; row 7 is the caller's |
| 4.2 Build the gene-reference projection module | 12, 14, 16, 20, 24 |
| 4.3 Port the spike's X encoding | 23, 26, 31 |
| 4.4 Encode obs and var | 13, 22, 40 |
| Derived values | 18, 27, 28, 29 (writer half), 32, 33 (writer half), and the contributor `uns` allowlist |
| 4.6 Assert identity uniqueness on every write | the joinid check in 15.6 |
| 4.7 Emit the store-root records: inventory and summary | 27, 38 (writer half) |
| Spatial, after Phase 3 | 34 |
| Ingestion self-checks | the writer rows in 15.6 |

Every other row is the catalog builder's, the release builder's, replaced, or dropped with its reason above.

## Iteration, by row

| Iteration | Rows |
| --- | --- |
| 1 | 1, 6, 7, 8, 12, 13, 14, 16, 20, 22, 24, 26, 31 (plain CSR), 38, 40 (the drop), and the contributor `uns` allowlist |
| 2 | 2, 4, 9, 10, 11, 18, 25, 27, 28, 29, 32, 33, 35, 36, 37. Row 11 moves to 1 if the Iteration 1 gate compares against Census 1 |
| later | 34; the `X` encoding of 31 (design 6.10); the rename of 40 |
| none | 3, 5, 15, 17, 19, 21, 23, 30, 39 (dropped) |
