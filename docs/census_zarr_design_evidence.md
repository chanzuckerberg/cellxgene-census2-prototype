# Census AnnData-Zarr: measurements and methods

The evidence behind the [design record](./census_zarr_design_v3_ste100.md). This document is the source of truth for every figure quoted there: what was measured, how, and what each measurement does not establish.

The raw results are in `tools/zarr_layout_spike/results/*.csv`. This path and the other `tools/`, `docs/`, and `dev-docs/` paths in this document name files in the [cellxgene-census](https://github.com/chanzuckerberg/cellxgene-census) repository or on its spike branches. They are not in this repository. Sections are numbered 9 onward because they keep the numbers that version 1 of the design cited. References written "spec section N" point into version 1.

Two terms recur. A **rung** is one scale point of the measurement series, named for how many source datasets it holds: rung 01 is one dataset, rung 25 is twenty-five. The **ladder** is the series of rungs, run end to end so a figure can be read as a curve against dataset count rather than as a single number.

**Last edited:** 2026-08-04.

---

## 9. Read-path feasibility evidence

**Delivered.** 28 measured queries across rungs 01, 05 and 25 in two experiments, plus an `obs`-predicate probe. Every non-empty query agreed on cells, genes and `X`, including `normalized` within tolerance, with the only differences being the axis-metadata encoding class of 9.5. Per-query bytes and requests are in `tools/zarr_layout_spike/results/bench.csv`.

The read-path question resolves to a mechanism rather than a ratio:

> A single-column predicate locates 565 rows in 158,982,719 for **0.4 MB and 9 requests**. Adding one clause that
> cannot change the selection defeats TileDB's fragment pruning, after which the `obs` axis is scanned. Every other
> cost result in this section follows from that.

9.3 is the subsection to read if you read only one.

**Outstanding**, none of it blocking the layout decision: rung 100 and the extrapolation to 1,845 datasets, which would confirm rather than discover, since 9.2 gives the scaling law; in-region wall clock, deliberately not attempted (10.1); and cold versus warm, not separated: every figure here is a cold read.

### 9.1 What the ladder compares, and what it cannot

**The two stacks do not read the same corpus, so their byte columns are not a head-to-head.** The SOMA side answers every query from the whole published corpus; the Zarr side answers from a store holding only the rung's datasets, because a full-corpus Zarr copy is 1.4 TiB and does not exist. The ratio between the two columns is dominated by corpus scope, not by layout. What the run legitimately supports is the Zarr scaling curve across rungs, the SOMA cost of retrieving a small subset from the full corpus, and parity at every rung. Layout-versus-layout at equal scope is Section 10, where every candidate holds the same data.

**Size-ascending rungs select secondary data, which silently disables the WMG query.** Of the 100 smallest datasets in the release, 73 carry no primary cells, and rungs 01 and 05 carry none between them, because secondary datasets are typically subsets or re-annotations of larger ones. Since the WMG filter requires `is_primary_data == True` it selects zero cells below rung 25, and the rungs also skew shallow: median `nnz` per cell at rung 01 is 292 against the filter's floor of 500. A defect in the manifest selection criteria (11.5), not in the queries; fixing it invalidates measurements already taken.

**Two metric floors.** Amplification divides compressed bytes read by an uncompressed payload (`bytes_wanted` is `nnz * 12`), so below about 3x it measures compression rather than over-fetch: two rung-25 queries report 0.5x and 0.6x. Never quote it without the payload beside it. And `bench.py` re-probes filter values from each subset, which is what stops the queries selecting nothing, so each rung binds a *different* cell type. The rungs are the same query *shape* at growing scale, not the same query: request counts and byte totals compare across rungs, anything payload-normalised does not.

### 9.2 The current stack's cost is fixed per query and proportional to the experiment

Seven queries at three rungs against the published corpus, `homo_sapiens`. Every rung reads the real corpus and only the manifest scope in the predicate changes, so this carries no scope confound.

| Query | Cells, 01 → 05 → 25 | SOMA MB, 01 → 05 → 25 | SOMA requests |
|---|---|---|---|
| `cell_type_slice` | 391 → 3,799 → 11,470 | 150.1 → 154.7 → 163.6 | 294 → 294 → 294 |
| `dataset_slice` | 565 → 3,799 → 15,502 | 99.1 → 103.7 → 107.8 | 287 → 287 → 287 |
| `tissue_slice` | 565 → 3,799 → 30,714 | 101.6 → 106.2 → 171.5 | 287 → 287 → 291 |
| `both_layers` | 391 → 3,799 → 11,470 | 248.8 → 259.5 → 277.1 | 335 → 335 → 335 |
| `multivalued_disease` | 565 → 3,158 → 35,808 | 103.2 → 104.7 → 173.2 | 291 → 291 → 295 |
| `cell_type_gene_subset` | 391 → 3,799 → 11,470 | 147.3 → 153.6 → 163.0 | 294 → 294 → 294 |
| `wmg_production` | 0 → 0 → 2,728 | 625.9 → 625.9 → 451.9 | 353 → 353 → 331 |

**Request count is a constant of the query shape**, identical to the unit at rungs 01 and 05 and almost invariant at rung 25 across a 100x change in result size: `cell_type_slice` issued exactly 294 requests at 391, 3,799 and 11,470 cells. For the metric 10.1 identifies as dominant, TileDB's does not respond to how much data is returned.

**Cost is not monotonic in result size.** `wmg_production` returned the same 625,879,749 bytes at rungs 01 and 05, the same number to the byte, with five times as many datasets in scope, then **451.9 MB for 2,728 cells** at rung 25, 174 MB less to do strictly more. A zero-cell query is the most expensive on its rung: with nothing selected there is no `X` to fetch, so 625.9 MB in 353 requests is `obs`-axis predicate evaluation and nothing else, at 4x the cost of the queries that return data. Bytes do eventually move, at rung 25 `tissue_slice` rose about 65%, where rungs 01 → 05 gave 3--6%, but at 11x the cells for 1.6x the bytes the crossover is a knee, not a corner.

**The fixed cost is proportional to the experiment.** Rung 25 contains one *Macaca mulatta* dataset alongside its 24 human ones, so the same queries ran against two experiments differing **22.7x** in row count: 158,982,719 against 7,010,229. Cut within the rung, human costs 19.3x to 30.9x the macaque figure; cut instead against human at rung 05, where result sizes are close at 3,799 against 3,190 cells, 18.6x to 28.3x. Both bracket 22.7x, and the cleanest case is `wmg_production` matched on cells, where **both sides select exactly zero** and the cost still differs by 28.3x. Normalised, human reads 0.6 bytes per `obs` row for a single-column predicate and 3.9 for the WMG predicate; macaque 0.8 and 3.2. So:

> Bytes fetched are set by the number of rows in the `obs` axis and by which columns the predicate touches. They
> are very nearly independent of how many cells the query selects.

**This is the one property that gets worse as the Census grows.** A per-query cost proportional to corpus size means every selective read pays for every cell ever added. The corpus is 218M observations across 1,845 datasets today (spec section 2) and grows with each release; no chunk or shard tuning addresses that, because the `obs` axis is being scanned. 9.3 pins down why.

**Where the bytes go on the `X` side.** Instrumented with `tiledbsoma_stats`, a 167-cell read of `pan_troglodytes`, 0.11% of cells, 0.034% of values, wants 1.61 MB and fetches **510.67 MB**, an amplification of **317x** in 312 requests. The tile geometry accounts for it exactly: `X` is stored in **2048 x 2048** tiles under Zstd (`CENSUS_DEFAULT_X_LAYERS_PLATFORM_CONFIG`) and a compressed tile is atomic. The experiment is 78 row-bands by 12 column-bands; 167 scattered cells land in roughly 69 row-bands, and because the query reads every gene each touched band pulls all 12 column-bands: roughly 828 of 936 tiles, about 88% of the array, consistent with 88% of a ~580 MB compressed array. Request count is not the problem: 312 requests for 510 MB is 1.8 MB each, which is good coalescing. Read 317x with its payload beside it; since the numerator is a fixed cost the ratio falls by whatever factor the payload grows, so 317x is what a 167-cell read costs, not what the corpus costs.

Bounding these claims: two experiments, seven queries, one SDK, one laptop. The metric is deterministic counters, which reproduce exactly (10.1), so the ratios carry no network variance. The macaque subset is one dataset against human's 24, so their predicates differ in the length of the `dataset_id` scope clause: measured to contribute nothing, which is what makes either cut usable. This describes TileDB on the published corpus; that a full-axis scan is avoidable in principle is not evidence that the mandated stack avoids it.

### 9.3 The mechanism: the second predicate column defeats fragment pruning

`obs`-only reads against `homo_sapiens`, projecting `soma_joinid` alone so the returned payload is identical throughout, varying nothing but the predicate.

| Predicate | Rows | Bytes | Requests |
|---|---:|---:|---:|
| `dataset_id in ['d7476…']` | 565 | **0.4 MB** | 9 |
| the same **+ `nnz >= 0`** | 565 | 258.2 MB | 26 |
| the same **+ two further always-true clauses** | 565 | **259.9 MB** | 50 |
| the same **+ `cell_type_ontology_term_id == 'CL:0000128'`** | 0 | 51.4 MB | 11 |
| the same **+ the WMG clauses** | 0 | **312.2 MB** | 58 |
| unscoped `nnz >= 500` | 150,300,999 | 259.8 MB | 33 |
| unscoped `cell_type_ontology_term_id == 'CL:0000128'` | 8,105,477 | 53.1 MB | 89 |

**A single-column predicate on a dictionary-encoded column is nearly free**: 0.4 MB and 9 requests to locate 565 rows in 158,982,719, because TileDB prunes to the fragments that can contain the value. The fixed cost of 9.2 is not inherent to selective reads against this corpus: it has a trigger.

**The trigger is the second column, and it is not a per-column cost.** Adding `nnz >= 0`, a clause that cannot change the selection, and the same 565 rows come back, takes the read from 0.4 MB to 258.2 MB. Adding two *further* always-true clauses brings it to 259.9 MB, 1.007x rather than 3x. So the second column does not add its own scan: it **defeats the pruning**, after which the axis is read anyway and further columns cost only their own compressed size.

Per-column full-axis scan cost varies an order of magnitude with compressibility: `nnz`, a wide-range integer, is 259.8 MB (1.6 bytes/row) against dictionary-encoded `cell_type_ontology_term_id` with 898 levels at 53.1 MB (0.33 bytes/row). The arithmetic closes on the headline figure: WMG references `nnz` and three dictionary columns, and 258 + 51 ≈ 309 MB against 312.2 MB measured. The ladder's larger 626 MB is that scan plus the eight `obs` columns the query projects, where the probe projects one.

**How large the cliff is depends entirely on which column is named, and the SDK's documented filters sit at the bottom of the range.**

| Second clause added to the `dataset_id` baseline | Bytes | vs 0.4 MB |
|---|---:|---:|
| `is_primary_data == True`: SDK | **0.9 MB** | 2.5x |
| `tissue_general_ontology_term_id == X`: SDK | **2.8 MB** | 7.7x |
| both SDK columns together | **3.4 MB** | **9.2x** |
| `cell_type_ontology_term_id == X` (dictionary, 898 levels) | 51.4 MB | 139x |
| `nnz >= 0` (numeric) | 258.2 MB | 697x |
| WMG's four clauses | 312.2 MB | 843x |

**The cliff spans 2.5x to 843x, and documented SDK usage costs 9.2x: 0.4 MB to 3.4 MB.** The mechanism is the same at both ends; the consequence is not. The expensive columns are the ones that do not correlate with row order, at 1.6 and 0.33 bytes per row against roughly 0.01 for `is_primary_data` and `tissue_general` (11.2). **The `nnz` and WMG figures are the worst case in the corpus, not the workload, and are not the headline.**

Two smaller results:

- **Selecting zero rows is expensive once pruning is defeated**, at 51.4 MB and 312.2 MB for the two zero-row predicates, which is why 9.2's zero-cell WMG query is the most expensive on its rung. The cheap experiment misses during the ladder run were single-column `dataset_id` predicates whose value was absent from the dictionary, so pruning eliminated everything.
- **`is_primary_data in [True, False]` raises** `PyQueryCondition.create_uint8() not found`: the `in` operator is unsupported against boolean columns. 6.3 has to accept that as a gap in the dialect it promises compatibility with, since a v2 filter that *did* support it would accept filters the oracle rejects.

**What this asks of the layout:**

> Multi-column predicate evaluation must not degrade into a full scan of the `obs` axis. A predicate over k columns
> must cost on the order of the k-column single-value case, not the k-column full-axis scan.

Per-chunk statistics would do it, since a Zarr reader can skip chunks whose statistics exclude the predicate; 10.2 measures how far that gets. It also gives 11.4's threshold a defensible form, and the figure to hold v2 to is the **SDK** one: 0.4 MB and 9 requests for a single-column scoped read, 3.4 MB for a documented two-column read. A threshold set at the 258 MB `nnz` case would clear a bar users never approach, though the layout should still fix that case, because a numeric-threshold filter is legal, WMG issues one, and 10.2 measures it as the residual per-chunk min/max cannot prune.

Bounding the claim: one experiment, one SDK, one session. The always-true clauses are synthetic by design, since a real clause would change the selection and confound the measurement. The per-column costs are specific to these columns' encodings.

### 9.4 Sparse `X` in AnnData has no gene axis

`cell_type_gene_subset` asks for 100 of 61,497 genes on the same cells as the all-genes query. Against the Zarr store it read **29,054,752 bytes in 206 requests**: byte-identical and request-identical to the all-genes query, not merely similar. On the SOMA side it was within 0.35%, and at rung 01 restricting genes bought 1.9%: 147.3 MB against 150.1 MB for a payload 149x smaller.

The store explains it. AnnData writes sparse `X` in **CSR as flat one-dimensional arrays**: `data` and `indices` of 51,960,590 elements in 256 chunks of 202,972, plus `indptr`. The chunk grid is over *flattened nonzeros*, so **there is no gene axis to chunk on**, and no chunk or shard size can make a gene-selective read cheap: the reader fetches every chunk spanning the selected rows and discards unwanted genes from `indices` in memory.

Three consequences:

1. **TileDB's 2048 x 2048 tiling has a gene dimension and AnnData's sparse encoding does not.** On gene-selective reads the current stack is structurally more capable than the mandated target's default encoding. Recovering the axis takes a layout decision with real cost, a gene-major (CSC) copy, a dense chunked-2D `X`, or per-gene-band arrays, and doing nothing means accepting the regression. 10.3 measures all four.
2. **The column-locality matched pair of 11.2 cannot discriminate** under this encoding, on either stack: it varies an axis that does not exist. The pair measured nothing, which is the vacuous-comparison failure class this harness was built to surface.
3. **Few objects does not mean few requests.** 256 chunks over two arrays is 512 range requests before `indptr` and `obs`, which is why reading 55% of the store's rows took **617 requests** where TileDB took **295** for the same cells. Auto-sharding delivered the object-count win, 88 objects for 213.6 MB, but a partial read still issues a range request per chunk.

That last figure is the sharpest cross-stack statement the run supports, and requests are not confounded in the direction that flatters Zarr: for this configuration, holding more data cannot require fewer requests for the same rows. The configuration is anndata's automatic default, chosen by `auto_shard_zarr_v3` rather than by anyone, so it measures the default and not Zarr, and 295 requests is the target Section 10 was given.

Store footprint across the rungs, 565 cells at 2.1 MB in 83 objects, 5,716 at 18.3 MB in 84, and 64,798 at 213.6 MB in **88**, is 115x the cells and 100x the bytes for five more objects, because auto-sharding packs each array into essentially one shard. Against the SOMA baseline of 1.4 TiB across 325,695 objects (spec section 2), object count is not the constraint for a sharded Zarr corpus. Bytes per cell is 3.3 kB here against SOMA's 6.9 kB corpus-wide, but the ladder's cells are shallower than the corpus mean, so treat 3.3 kB as a floor rather than an extrapolation.

### 9.5 Parity holds, and two encoding differences need a decision

Every non-empty query agreed on cells by source identity, genes by `feature_id`, and `X` on every layer the experiment carries: `raw` exactly, `normalized` within tolerance. The first cross-stack comparison, a 167-cell slice of `pan_troglodytes`, matched all 140,587 `raw` values exactly across 167 cells and 23,534 genes.

Every difference the checker reports falls in one class, axis-metadata encoding, and neither half raises. SOMA carries the experiment's full category dictionary for `obs.cell_type` and its ontology column: 25 cell types, of which 24 are absent from the Zarr round trip, which retains only the levels present in the slice. And `var.feature_name` is `str` via SOMA and `category` via Zarr, because AnnData's writer converts string columns to categorical on write. Both are the failure class this harness exists for: values identical, nothing raised, and the user still receives something different. A comparison that decategorised before comparing, as the builder's own validation does, would report a clean pass.

Two decisions this leaves open (spec section 10). **Does the published corpus preserve the full category dictionary or only observed levels?** A caller who filters and then inspects `.cat.categories` gets a different answer under v1 and v2; both behaviours are defensible, differing silently is not. **Is `var.feature_name` a string or a categorical in v2?** That is AnnData's writer behaviour rather than a design choice, which is why it is worth settling before it becomes the contract.

Bounding the claim: the Zarr side was converted **from** the SOMA read, so this tests format fidelity and reader equivalence, not independent derivation: that arrives when the builder emits Zarr. The missing categories are partly an artefact of converting a slice rather than a whole experiment, so what the finding establishes is that category handling differs between the paths and needs specifying, not that the corpus would lose 24 levels.

## 10. Layout comparison

**Complete for `obs` and for the `X` encoding; one item open.** Every decision, with the number that decided it:

| Axis | Chosen | Rejected, with the number that rejected it |
|---|---|---|
| `X` encoding | **Banded CSR, 16 gene bands** (10.3) | Plain CSR: reads the whole 98.7 MB matrix for one gene. CSC: reads it all for any row query. Dense 2-D: works and costs only 1.38x, but banding matches it more cheaply on today's queries |
| `obs` chunking | **Uniform 8192-row chunks** (10.2, 10.5) | AnnData's default one-chunk-per-column leaves nothing to prune. Dataset-aligned chunks prune 10x better and read 3x more bytes |
| Sharding | **Shard** (10.2) | Unsharded: 28x the objects. Sharding needs adjacent candidate chunks to help requests |
| `obs` index | **A per-chunk statistics index**, applied adaptively below a ~1/20 survivor threshold (10.2) | No index: reads the whole column where a real `dataset_id` scope needs 2.6% of chunks. Always-on: costs ~3 requests per surviving chunk, which loses on broad predicates |
| Row order | **Assignment order, unchanged** (10.4) | Clustering by cell type: gains 6% on its own query's `X` read, costs 210% on another's |
| Store organization | **One store** (10.5) | Per-dataset: ~50 kB of metadata per store opened, 4.6 MB for a documented lung filter and 54 MB corpus-wide, against a 3.4 MB baseline |
| `X` chunk edges | **Padded to dataset boundaries** (10.5) | Unaligned: 17--35% more bytes and up to 2.2x the floor. Finer chunks instead: more bytes, 1.8x the requests, a larger store |

**Nothing is open on the layout.** The last item, aligning `X`'s nonzero chunk edges to dataset boundaries, where the byte saving attributed to per-dataset stores actually lived, is settled in 10.5: padding recovers 17--35% of the bytes and changes the store by −0.1%.

The two corpus-level metrics this section was asked to record, against the SOMA baseline of 1.4 TiB / 325,695 objects (spec section 2): **object count is not a constraint for a sharded layout**, 88 objects held 213.6 MB and 115x more cells added five (9.4), and **bytes per cell measured 3.3 kB against SOMA's 6.9 kB**, a floor rather than an estimate. On this floor the Zarr corpus would not exceed the SOMA one it replaces.

### 10.1 How candidates are ranked: counters, not wall clock

**Request count is ranked first, on prior art.** `chanzuckerberg/multimodal-slicing` converted CELLxGENE datasets to Zarr on S3: reading 16,527 cells originally issued **33,054 S3 requests**, two per cell, and restructuring the read to fetch contiguous ranges reduced that to **2**: a reported ~150x improvement, taking the read from over ten minutes to 4--7 seconds. With requests minimised, adding 8 Dask workers produced no speedup at all (3% slower) and the project raised its parallelism threshold to 50,000 rows, concluding that "I/O dominates over CPU processing time." Two consequences: bytes and CPU throughput are secondary, which is why object count is a first-class metric here; and request count is a property of layout *and* access pattern together, since those 2 requests came from coalescing rather than from storing fewer objects, so every candidate must be read through a path that coalesces or the comparison measures the reader. That work aligned datasets to a 79,948-gene universe of its own, covers human only, has no spatial data and is not census-harmonized, so its numbers are directional.

**Wall clock from outside `us-west-2` is unusable, and the counters are exact.** One 167-cell query, run repeatedly from a laptop in one session, took ~475 s, then ~27 s as the mean of three consecutive runs, then 39 s, 207 s and 349 s: an 18x spread with no cold-then-warm pattern. The control is stronger than the spread: the rung 01 run was executed twice an hour apart on the same laptop, and all seven queries returned **byte-identical counters**: the same `VFS.read_byte_num` to the byte, the same request count, the same LIST/HEAD count, while wall clock for those same reads ranged from 0.17x to 5.34x of the first run, a **31x spread**. `wmg_production` moved 625,879,749 bytes in 353 requests both times, taking 49.5 s and then 264.1 s; `multivalued_disease` moved 103,180,185 bytes in 291 requests, taking 128.2 s and then 21.4 s. The identical byte counts *prove* the two runs did the same work, so the time difference cannot be a difference in what was read.

Three consequences, all adopted by 11.4: deterministic counters are the primary metric and wall clock secondary; a no-regression threshold is stated in requests and bytes, which is possible precisely because those are exact rather than merely reproducible; and any wall-clock comparison must run **in-region** on an EC2 instance in `us-west-2`, which is why no timing claim appears in this document.

**A naive reader makes every layout look bad.** A scattered-row `X` read initially reported **36 GB from a 98.7 MB store**, because it issued one range per contiguous row run and re-fetched the same chunk for every row in it. Every figure in this section comes from a reader that coalesces by chunk, with a regression test asserting that no scattered read exceeds a full read.

### 10.2 `obs`: chunk it, shard it, and index the multi-column predicates

The sweep that tests whether a Zarr layout can hold a multi-column predicate near the cost of a single-column one (9.3). Four axes, chunk size, shard size, row order, chunk statistics on or off, over an `obs` axis scaled to **5,000,000 rows**, all local, counters exact. 64 measurements in `results/sweep.csv`.

Three preliminaries that shaped it. **AnnData writes each `obs` column as one chunk**, measured `chunks=(64798, )` on every column of the rung 25 store, so there is no granularity to prune at and chunking `obs` is a precondition, not a tuning. **Zarr carries no chunk statistics**, unlike Parquet, so they live in a file, which makes the index an *extension* to the AnnData-Zarr layout needing a defined home in the schema (spec section 3.1, section 5). And **the rung 25 store is too small to measure this on**: at 64,798 rows a compressed `obs` column is a few kB and the index is comparable in size to the data it prunes. The 5M-row axis is built by repeating whole dataset blocks so `dataset_id` stays contiguous and cell types stay distributed within a dataset as they really are; it is unsuitable for any claim about corpus content.

**The standing configuration and what it costs.** Assignment order, 8192-row chunks, 65536-row shards: the layout 3.1 specifies. Predicates are a single dictionary column, that column plus `nnz >= 500`, WMG's four clauses, and a value absent from the axis:

| Predicate | Index on | Index off | Bytes saved | At 20M rows |
|---|---:|---:|---:|---:|
| 1 column | 0.53 MB, 626 reqs | 1.23 MB, 100 reqs | 2.3x | 2.10 / 4.90 MB, 2.3x |
| 2 columns (`+ nnz >= 500`) | **4.19 MB**, 1,005 reqs | 10.38 MB, 181 reqs | 2.5x | 16.75 / 41.49 MB, 2.5x |
| 4 columns (WMG) | 8.88 MB, 4,664 reqs | 10.44 MB, 346 reqs | 1.2x | 35.50 / 41.70 MB, 1.2x |
| absent value | **0.02 MB, 16 reqs** | 1.23 MB, 100 reqs | 61x | 0.03 / 4.90 MB, 157x |

**The index buys bytes and costs requests.** Pruning reads 34% of the chunks for the one- and two-column predicates, saving 2.3--2.5x on bytes, but under assignment order those surviving chunks are scattered across shards and each costs its own range request, where an unpruned read fetches whole shards. That is 3.0's second qualification measured: sharding's request benefit needs adjacent candidates. So the index is worth applying where pruning is *decisive* rather than partial, and the absent-value row is the clean case, better on both metrics by 61x and 6x. It is not a toy either: 9.3 shows a predicate naming a value the corpus does not hold is a real and expensive query on the current stack. Which predicates take the index is therefore a reader policy the layout should state rather than leave to the implementation (Section 14).

**The trade does not improve with scale, and the reason bounds what this sweep can answer.** Repeated at **20,000,000 rows** (`results/sweep_20m.csv`, same 64 configurations), every absolute figure grows about 4x and every ratio holds: the byte saving is 2.33x on one column against 2.31x at 5M, 2.48x on two columns against 2.48x, and the request penalty widens slightly, 6.3x to 7.5x on one column. Only the decisive case improves, from 80x to 157x on bytes.

That flatness is structural rather than a result. The ratio is set by `fraction_read`, the share of chunks the filtered value occupies, and it came out **0.340 at 20M against 0.339 at 5M**: identical to three decimals, because the axis is scaled by repeating whole dataset blocks, so a value's chunk occupancy is a property of the 64,798-row block cycle and cannot vary with row count. **This sweep therefore measures the index at one fixed value distribution**, and that distribution is the thing the answer depends on. So it was measured on the real corpus.

#### What a real predicate actually occupies, measured on the published corpus

`probe_chunk_occupancy.py`, over `census_data/homo_sapiens` in assignment order, 158,982,719 rows, 19,408 chunks of 8192, asking for each predicate how many chunks contain a match. `soma_joinid` was verified contiguous from 0 first, since positional chunking is only assignment order if it is. Reads are one column at a time and total about 310 MB. Results in `results/chunk_occupancy.csv`; the same run against `macaca_mulatta` is there for contrast.

| Predicate | Rows matched | Of the axis | Chunks touched | Spread |
|---|---:|---:|---:|---:|
| `cell_type == CL:0000540`: the most common | 10,035,322 | 6.3% | **0.124** | 2.0x |
| `cell_type == CL:0002365`: median rank | 16,371 | 0.010% | **0.0030** | 29x |
| `cell_type == CL:0008024`, rarest | 8 | ~0 | **0.0001** (1 chunk) | n/a |
| `tissue_general == UBERON:0000955`: brain, SDK-documented | 48,581,856 | 30.6% | 0.328 | 1.1x |
| `dataset_id == 37a17b78…`: the largest | 4,140,453 | 2.6% | **0.0261** | 1.0x |
| `nnz >= 500`: WMG's numeric clause | 150,300,999 | 94.5% | 0.999 | 1.1x |

"Spread" is chunks touched divided by rows matched, so 1.0x is a perfectly contiguous value and higher numbers mean the rows are scattered. Three things follow.

**The synthetic axis was pessimistic, and by more than a little.** The hardest real dictionary case, the single most common cell type in the corpus, touches **0.124** of chunks where the synthetic axis reports 0.339, and a median-frequency cell type touches **0.0030**, a hundredfold better. Sweep figures should be read as an upper bound on cost, not as the expected case.

**Assignment order is dataset-contiguous, now measured at corpus scale.** `dataset_id` touches 0.0261 of chunks to return 2.6% of rows, a spread of 1.0x: a dataset's rows occupy exactly the chunks they have to and no others. That is the property 10.4 relies on when it keeps assignment order, and it had not previously been checked above 64,798 rows.

**The `nnz` residual is confirmed and is partly not the layout's fault.** `nnz >= 500` keeps 99.9% of chunks, but it also *matches* 94.5% of rows, so most of that is a genuinely unselective predicate rather than a pruning failure. A filter that wants almost every row cannot be helped by skipping chunks.

#### The policy this settles: prune adaptively, on a threshold the reader can evaluate for free

The index is read first and is cheap, 0.012 MB and 16 requests at 5M rows, 0.15 MB at 20M, and reading it tells the reader **how many chunks survive before it fetches any column data**. That makes the choice measurable at query time rather than fixed at write time.

The break-even follows from the sweep's own counters. A pruned read costs about 3 requests per surviving chunk (2,447 requests for 830 chunks at 20M; 626 for 207 at 5M), while an unpruned read costs about one request per shard, and a shard holds 8 chunks. Pruning therefore wins on requests when survivors fall below roughly **(1/8) / 3 ≈ 4-5% of chunks**, which the measured configurations put at 33 of 611 chunks at 5M and 111 of 2,442 at 20M. Below that line pruning wins on both metrics at once; above it, it trades bytes for requests and 10.1's ranking says to decline the trade.

Applied to the real predicates above, at human scale the line sits at about 889 of 19,408 chunks:

| Predicate | Chunks touched | Prune? |
|---|---:|---|
| `dataset_id ==` one dataset | 506 | **yes**, and this is the scoping clause most queries carry |
| `cell_type ==` median rank | 58 | **yes**, by two orders of magnitude |
| `cell_type ==` rarest | 1 | **yes** |
| `cell_type ==` most common | 2,411 | no: read whole shards |
| `tissue_general ==` brain | 6,361 | no |
| `nnz >= 500` | 19,396 | no |

So the index earns its keep on the selective end, which is where 9.3 located the read-path problem, and the reader declines it on the broad end where a whole-column read is cheaper. **3.1 adopts this as the index scope rule.** What remains unmeasured is the threshold's exact position on S3 rather than a local filesystem, since it balances a request count against a byte count and only their relative latency decides the crossover (10.1).

**Sharding buys the object-count property outright.** At 2048-row chunks the unsharded store is **20,642 objects** against **726** sharded, a 28x reduction, with bytes unchanged (0.11 vs 0.12 MB). Against the SOMA baseline of 325,695 objects that removes object count as a concern for a sharded layout, while showing how an unsharded one would blow past the baseline.

**8192-row chunks are chosen on requests, not bytes.** Against 2048 under the same order and pruning, 8192 issues 1,005 requests to 2048's 1,850 on the two-column predicate and 4,664 to 7,704 on WMG, roughly 1.7x fewer, while 2048 reads 25--35% fewer bytes (3.15 against 4.19 MB; 6.54 against 8.88). 10.1's ranking decides it. Both directions hold at 20M rows, 3,951 requests against 7,303 on two columns, and 12.57 MB against 16.75, so the choice is not an artefact of the smaller axis. Do not go finer than 8192 without re-measuring.

**Does it meet 9.3's target?** Partly. A single-column predicate costs 0.53 MB and two columns 4.19 MB, a **7.9x cliff**, where TileDB's on the same `nnz`-threshold shape is 697x. The two axes differ in scale, so treat that as an order-of-magnitude statement rather than a like-for-like ratio: the cliff is roughly 90-fold smaller, and it is not gone. The residual is `nnz >= 500` itself: min/max statistics cannot exclude a chunk when the threshold falls inside almost every chunk's range, so the `nnz` column is read for every surviving candidate. That is a property of the *predicate shape*, not the layout: equality on a dictionary column prunes to almost nothing, a threshold on a high-cardinality numeric column prunes to almost everything, which is why 3.2's diagram cannot show the `nnz` clause being pruned. Materialising the threshold as a stored boolean would reduce that clause to a dictionary predicate, and is rejected because the corpus should not accrete a column per consumer.

Bounding these claims: the axis is synthetic above 64,798 rows and 5M is 1/32nd of the human experiment. Reads are local, so request counts are real and their S3 latency consequence is not measured and cannot be from here (10.1). Only `obs` is swept: `X` is 10.3 and 10.4, and 9.4 established `X` has no gene axis to sweep. The pruning correctness check, pruned and unpruned reads returning identical rows, passed in all 32 configurations, which is the precondition for any of these numbers meaning anything.

### 10.3 `X` encoding: banded CSR recovers the gene axis

Four candidate encodings of the same matrix, 64,798 cells x 61,497 genes, 51,960,590 nonzeros, 1.30% dense, each written from the rung 25 store, with every encoding required to return identical values for every read.

| Encoding | Store | Objects | 1 gene | 100 contiguous genes | all genes, contiguous rows |
|---|---:|---:|---:|---:|---:|
| `bands` (16 gene bands, CSR each) | **92.7 MB** | 217 | **0.51 MB** | 0.31 MB | **16.62 MB** |
| `csr` (AnnData's default) | 98.7 MB | 13 | 98.71 MB | 16.68 MB | 16.68 MB |
| `csc` | 115.1 MB | 13 | 4.65 MB | 4.65 MB | 115.13 MB |
| `dense` (2-D, 2048x2048 chunks) | 136.5 MB | 26 | 1.04 MB | **0.20 MB** | 22.02 MB |

**CSR reads the entire matrix to return one gene**: 98.71 MB of a 98.7 MB store for 1 of 61,497 genes. This is 9.4's finding at the store boundary, and it is total rather than partial. Banding costs 0.51 MB for the same query, a **194x** reduction.

**Banding is the best all-rounder and has the smallest store.** Its row-selective read is within 0.4% of CSR's (16.62 against 16.68 MB) and its contiguous gene-block read is 54x better. Splitting genes into 16 bands narrows each band's column indices from 0--61,496 to 0--3,843 and `indices` compresses better as a result, which is why 217 objects still cost fewer bytes than CSR's 13. Its price is request count: 148--264 against CSR's 18--20.

**A dense 2-D `X` costs 1.38x the sparse store, not the 15.9 GB the arithmetic suggests.** At 1.30% density the zeros compress to almost nothing, so the objection that has kept dense off the table is wrong by two orders of magnitude. Dense buys true two-dimensional chunking, at 2048x2048, exactly the TileDB tile geometry of 9.2, for 38% more storage, and it wins outright on contiguous gene blocks. It is rejected only because banding matches it more cheaply on the queries that exist today (spec section 3.3).

**CSC is the mirror-image failure and is dropped.** It reads the whole array for any row-selective query (115.13 MB), the exact pathology CSR has on the gene axis, moved to the other axis, and is never the best choice except for a query no consumer issues.

**Gene-axis locality only pays for *contiguous* gene selections, which no current consumer asks for.** A contiguous 100-gene block costs 0.31 MB banded and 0.20 MB dense; 100 *scattered* genes cost 16.62 and 21.59, some 50x and 100x worse, because 100 random genes hit all 16 bands and most column chunks. A real gene query is a marker panel scattered through `feature_id` order, so recovering the benefit needs the genes *reordered* so co-queried ones are adjacent: a second clustering decision, on the other axis, with no obvious key (spec section 3.3). Without it, banding's advantage applies to single-gene lookups and little else.

Bounding these claims: one matrix, one layer, 64,798 cells: a rung 25 store, not the corpus. Local filesystem. `raw` only; the `normalized` layer would roughly double every store figure. Band count fixed at 16, not swept.

### 10.4 Row order: assignment order is kept

Row order is a layout choice the current corpus does not exercise, and clustering by a commonly-filtered attribute is the obvious candidate. Measured on `X` in the banded encoding, same matrix, two real predicates, two row orders, against the theoretical floor of store size times the fraction of nonzeros the query needs:

| Predicate | Rows | Floor | Assignment order | Clustered by cell type |
|---|---:|---:|---:|---:|
| `cell_type_ontology_term_id == 'CL:0002627'` (17.5% of nonzeros) | 11,470 | 16.26 MB | 19.55 MB (1.20x), 250 reqs | **18.42 MB** (1.13x), 206 reqs |
| `dataset_id == 'd3566d6a…'` (3.7% of nonzeros) | 15,502 | 3.42 MB | **5.92 MB** (1.73x), 152 reqs | 18.34 MB (5.4x), 314 reqs |

**Clustering gains 6% on the query it is built for and costs 210% on the other.** That is not a trade worth making, and the reason is the useful part: **real predicates already have row locality under assignment order, because datasets are biologically homogeneous.** A cell type concentrates in the few datasets that assayed it, so its rows are not scattered through the corpus even though nothing sorted them. Assignment order is within **1.2x--1.7x of optimal** for both predicates with no clustering at all, which leaves clustering very little to win and plenty to lose.

**The `obs`-side gain from clustering is real, and larger than it first looked.** Measured on the 20M-row axis (10.2), clustering cuts a single-column predicate from 2.10 MB to 0.11 MB and WMG from 35.50 MB to 1.87 MB, saving 2.0 MB and 33.6 MB: which scales to roughly 16 MB and 270 MB per query at the human experiment's 162M rows. It cuts requests hard too, 2,447 to 84 on one column, because clustering is what makes the surviving chunks adjacent.

So the rejection does not rest on the `obs` term being small. It rests on the `X` term being larger, and on which way each scales. **A clustered corpus makes a `dataset_id`-shaped selection scattered through the whole axis**, so its `X` read touches a corpus-proportional share rather than one contiguous block: at 64,798 cells that took 5.92 MB to 18.34 MB, essentially the whole store. The `obs` saving and the `X` penalty therefore both grow with the corpus, and `X` is three orders of magnitude the size of one `obs` column (10.5). A few hundred megabytes saved on `obs` does not pay for a corpus-proportional share of a terabyte-scale `X`.

**That argument is not measured at corpus scale**, and it is the weakest link in this decision. The `X` penalty is one matrix at 64,798 cells; the extrapolation is from the mechanism, not from a number.

**Assignment order is kept.** It is dataset-contiguous, already near-optimal for both predicate shapes measured, and free. The effort belongs on the `obs` index instead, which attacks the predicate-evaluation cost 9.3 identifies as the dominant one and still prunes under assignment order, at 2.3--2.5x on bytes (10.2).

Bounding these claims: one matrix, 64,798 cells, one clustering key, two predicates, `raw` only. The floor calculation assumes compression is uniform across the store, which it is not exactly. And these are the two predicates the manifests make available at this rung; a query whose rows genuinely are scattered, a rare cell type present in hundreds of datasets, would behave differently, and none was tested. That is the one case that could reopen this decision (spec section 3.3).

### 10.5 One store, uniform `obs` chunks, dataset-aligned `X` chunks

Two related questions: whether to split the corpus into per-dataset stores, and whether to place chunk boundaries at dataset starts within one store. Both are measured. Per-dataset stores lose; alignment loses on `obs` and **wins on `X`**, which is where the byte saving that motivated it was coming from all along. The per-dataset half is then rechecked against a corpus already deployed in that layout, which is the only measurement in this section taken on S3 rather than local disk.

**Per-dataset stores.** One store holding all 24 rung-25 datasets against 24 per-dataset stores plus a shared index of per-dataset summaries, so a predicate can rule a store out without opening it. Both were required to select the same cells for every query.

| Query | Stores touched | Monolithic | Per-dataset |
|---|---:|---:|---:|
| `dataset_id ==` one dataset | 1 of 24 | 5.93 MB, 167 reqs | **3.22 MB, 132 reqs** |
| `cell_type ==` | 2 of 24 | 19.56 MB, 264 reqs | **15.64 MB**, 276 reqs |
| WMG | 3 of 24 | 17.39 MB, 335 reqs | **11.58 MB**, 400 reqs |
| all datasets (`nnz >= 0`) | 24 of 24 | 92.79 MB, **275 reqs** | 93.92 MB, 3,021 reqs |
| Store total | | **93.6 MB, 268 objects** | 108.2 MB, 3,449 objects |

**Monolithic wins the moment a query fans out.** At full fan-out per-dataset issues **3,021 requests against 275**, an 11x penalty for the same bytes, which is **114 extra requests per store opened**. That per-store figure describes *this* build and does not extrapolate: these stores carry no consolidated metadata, so opening one costs a separate request for every group and array document. The deployed measurement below deals with a per-dataset corpus that does consolidate, where the same cost is **10.1 requests per store**, and shows that consolidation moves the cost rather than removing it. Object count scales with dataset count either way: 144 objects per dataset store implies about 265,000 at 1,845 datasets, close to the SOMA baseline, where a monolithic store grows with data volume and sits at 268 objects for 93.6 MB.

**And the fan-out is measured, not assumed.** 11.2 records the fan-out of the SDK's own documented filters over the 1,086 human datasets: the baseline is `get_anndata`, not any one downstream pipeline. `tissue_general == 'vasculature'` touches 10 datasets, `is_primary_data == True and tissue_general == 'lung'` touches 45, `tissue_general == 'lung'` touches 92, and `tissue_general in ['brain', 'lung']` touches 447. The measured crossover is low: per-dataset is cheaper at one store, level at two, behind at three, so every documented filter sits past it. What each of those fan-outs costs is priced below, against a deployed per-dataset corpus rather than this rung.

**Dataset-aligned chunk boundaries within one store.** The per-dataset byte advantage is boundary alignment: in a monolithic store a dataset's rows share their first and last chunks with its neighbours', so those rows are fetched and discarded. Alignment promises that saving without the per-store request cost. Measured on the same 64,798-row `obs`:

| Metric | `aligned_chunks` (exact, unsharded) | `aligned_shards` (rounded, sharded) | `uniform` (8192 / 65536) |
|---|---:|---:|---:|
| Rows scanned, `dataset_id ==` | **15,502**: exactly the dataset | 21,504 | 24,576 |
| Rows scanned, WMG | **3,752** | 34,816 | 40,960 |
| kB read, `dataset_id ==` | 16.7 | 27.1 | **5.5** |
| kB read, SDK 2-column | 17.4 | 33.2 | **6.4** |
| kB read, WMG | 46.1 | 146.6 | 118.0 |
| Objects | 229 | 232 | **50** |
| Index size | 23.5 kB | 39.3 kB | **9.6 kB** |

**Alignment prunes an order of magnitude better and still reads 3x more bytes.** It does exactly what it was designed to on rows scanned, but dataset-sized chunks are small, median 1,324 rows against uniform's 8,192, and small chunks compress badly. This is 10.2's chunk-size result restated: on the `obs` axis compression efficiency dominates, and a layout that shrinks chunks to gain precision loses more to the compressor than it saves in rows. `aligned_shards` is the useful negative result: rounding each dataset up to a multiple of 512 rows keeps sharding but gives up exactness, and it scans nearly as much as uniform while paying alignment's object and index costs.

A second reason not to ship it: variable-width chunks exist as Zarr v3's *rectilinear chunk grid*, but in zarr 3.2.1 they are **experimental and disabled by default** (`zarr.config.set({'array.rectilinear_chunks': True})`) and **cannot be combined with sharding**. Zarr directs you instead to rectilinear *shards* with uniform chunks inside, each shard a whole multiple of its inner chunk; since real dataset widths at this rung run 146 to 15,502 rows, exact alignment and sharding are mutually exclusive. A recommendation resting on an experimental feature and forfeiting either sharding or exactness is not one this design makes (spec section 5).

**Where the 20--46% byte saving actually lived: `X`.** The store-organization figures are megabytes and the `obs` alignment figures kilobytes, because those queries read `obs` **and** `X` and `X` dominates by three orders of magnitude. The saving was on the `X` side, a dataset's nonzero ranges in a banded CSR array stop sharing chunks with its neighbours', and almost none of it was `obs`. So the decision splits by axis, and the `X` half is measured below.

Bounding the `obs` half: 24 datasets, one rung, one organism. The fan-out extrapolation assumes the per-store cost is constant, which it should be since it is metadata. The per-dataset store duplicates every category dictionary, which is most of its 16% size overhead and would grow with dataset count. Dataset widths here (146 to 15,502 rows) are far smaller than the corpus median, and alignment's compression penalty shrinks as datasets grow, so a corpus-scale test could land differently, though it would have to overcome a 3x gap. Local filesystem throughout.

#### `X`: aligning nonzero chunk edges to datasets, by padding

Exact alignment needs Zarr v3's rectilinear chunk grid, which is experimental in zarr 3.2.1 and **cannot be combined with sharding**, so it is not a candidate whatever it reads. The buildable form keeps uniform chunks and **pads** each dataset's nonzero range up to a chunk multiple, so no chunk holds two datasets' nonzeros while sharding survives. Padding entries carry `indices = -1`, which is outside every band's gene range, so the *same* reader reads both layouts and the comparison is layout-only; `test_x_aligned.py` asserts value-identical results and mutates the sentinel to prove the assertion bites. Same matrix as 10.3: 64,798 cells, 51,960,590 nonzeros, banded CSR, 16 bands, 125,000-element chunks per band. "Floor" is the store size times the fraction of nonzeros the query needs, as in 10.4.

| Query | Uniform | Padded | Change | vs floor |
|---|---:|---:|---:|---:|
| `dataset_id ==` one dataset | 5.92 MB, 152 reqs | **3.86 MB, 142 reqs** | **−34.8%** | 1.73x → **1.13x** |
| `cell_type ==` cross-dataset | 19.55 MB, 250 reqs | **16.19 MB**, 282 reqs | **−17.2%** | 1.20x → **1.00x** |
| WMG | 17.27 MB, 266 reqs | **12.17 MB**, 268 reqs | **−29.6%** | 2.23x → 1.57x |
| Store | 92.7 MB, 217 objects | **92.6 MB**, 239 objects | −0.1% | |

**Padding recovers the saving and costs nothing to store.** 17--35% fewer bytes on every query, and the cell-type read lands exactly on its floor. The storage objection, the reason this was expected to be a trade, does not survive contact: 19.16 million padding elements, 37% as many as there are real nonzeros, changed the store by **−0.1%**, because a run of identical zeros and −1s compresses to almost nothing and aligned chunks hold one dataset's values rather than two datasets' mixed together. Object count rises 10%, and requests are a wash: −10 on the dataset query, +32 on the cell-type one, +2 on WMG.

**It also beats the obvious alternative, which is simply using a smaller chunk.** Shrinking the chunk reduces edge spillover too, and the same run at 25,000 elements per band shows what that costs:

| Layout | `dataset_id` read | Requests | Store |
|---|---:|---:|---:|
| uniform, 125,000/band | 5.92 MB | 152 | 92.7 MB |
| uniform, 25,000/band | 4.31 MB | 272 | 98.4 MB |
| **padded, 125,000/band** | **3.86 MB** | **142** | **92.6 MB** |

Padding at the coarse chunk dominates the fine uniform chunk on all three metrics at once. That is the same compression effect 10.2 found on `obs`, in the other direction: fine chunks buy precision and lose it again to the compressor, where alignment buys the precision without shrinking anything. **3.1 therefore pads `X` nonzero ranges to dataset boundaries at the standing 125,000-element chunk.**

Bounding these claims: one rung, one organism, `raw` only, 24 datasets, local filesystem. The magnitude is scale-sensitive in a way the direction is not: the saving is edge spillover, which is roughly the chunk size over a dataset's nonzeros per band, and at this rung that ratio is about 0.9 because a dataset barely fills one chunk. Corpus datasets are much larger, so the percentage should fall; the 25,000-element run shows exactly that dependency, with the saving dropping to 4--8% when the ratio drops fivefold. **Do not quote 17--35% as a corpus-scale figure.** What does carry is the dominance argument above, which is structural, and the storage cost, which is nil at both chunk sizes.

#### The same layout as deployed: consolidation moves the cost into bytes

The build above carries no consolidated metadata, so most of its 114 requests per store are one-per-document metadata reads. A deployed per-dataset Census corpus does consolidate, and `probe_datalake.py` measures it: the multimodal-slicing POC's own query code against its own 662 per-dataset stores of Census human data on S3, over 11.2's documented filters, counted through the same `CountingStore` as everything above.

| Query | Stores | Matched | Requests | Per store | Bytes |
|---|---:|---:|---:|---:|---:|
| `tissue_general == 'vasculature'` | 200 | 2 | 2,019 | 10.1 | 21.0 MB |
| `tissue_general == 'lung'` | 200 | 20 | 2,019 | 10.1 | 21.0 MB |
| `is_primary_data == True and tissue_general == 'lung'` | 200 | 11 | 3,038 | 15.2 | 31.0 MB |

**Consolidation cuts requests by 11x and converts the cost into bytes.** Opening one store costs 4 requests and **49.7 kB**, because the consolidated document is fetched whole. It is fixed in size, 49,579 bytes for a 5.5 MB store against 49,772 for a 4.9 GB one, does not shrink with a more selective predicate, and is paid once per store opened. The per-store cost that decides this question is therefore **~50 kB, not 114 requests**:

| Documented filter (11.2) | Stores touched | Metadata alone |
|---|---:|---:|
| `tissue_general == 'vasculature'` | 10 | 0.5 MB |
| `is_primary_data == True and tissue_general == 'lung'` | 45 | 2.2 MB |
| `tissue_general == 'lung'` | 92 | 4.6 MB |
| `tissue_general in ['brain', 'lung']` | 447 | 22.2 MB |
| every human dataset, which is what a store-pruning index avoids | 1,086 | 54.0 MB |

Against 9.3, where the current SDK answers the documented two-column filter in **3.4 MB total**, the lung filter spends 1.3x that baseline on metadata before reading a cell, and a reader with no store-pruning index spends 16x it on every query regardless of what the query asks.

**Fan-out is priced, and does not reach these figures.** Driving the same work with threads rather than the POC's Dask removes scatter, serialization and scheduling, so this measures the ceiling rather than its implementation: 8 threads return 6.9--7.9x, 32 return 16.7--19.3x, and 64 return 19.3--25.0x. Concurrency is a real answer to latency and no answer to bytes.

Bounding it: 200 of the 662 stores, predicate only at that scale, laptop and off-region, so the fan-out ratios are indicative where the byte and request counts are not. Two costs in the POC are its own rather than the layout's and are excluded from every figure here: `filter_obs` opens each store twice, once in `_get_n_obs` and once in `_load_obs_field`, doubling the metadata it fetches; and it holds no index over its stores, which is what puts the 1,086-store row in reach. The layout is charged one open per store throughout.

### 10.6 The gene space is pinned upstream, and how much of it a build sees

`var` is the column space of `X`, so how it is defined decides whether a dataset can be written without touching its neighbours. Four measurements settle most of the question.

**The vocabulary is not ours.** The Census schema pins gene references to the CELLxGENE dataset schema, which publishes them per organism as `genes_<organism>.csv.gz` in `single-cell-curation`, currently GENCODE v48 for human and M37 for mouse, with an automated bump workflow driven by `gene_info.yml`. No curated dataset can carry a gene outside that reference, so the identifier space is already fixed and versioned upstream.

**What varies is how much of it a build has seen.** Today `var` holds the genes assayed by at least one included dataset. Measured against the `2025-11-08` LTS corpus and the current reference files:

| Experiment | Reference | Assayed | Ratio | Never assayed |
|---|---:|---:|---:|---:|
| `macaca_mulatta` | 35,432 | 34,676 | 1.02 | 756 |
| `homo_sapiens` | 78,894 | 61,497 | 1.28 | 17,397 |
| `callithrix_jacchus` | 40,188 | 28,348 | 1.42 | 11,840 |
| `pan_troglodytes` | 33,729 | 23,534 | 1.43 | 10,195 |
| `mus_musculus` | 78,317 | 53,384 | 1.47 | 24,933 |
| `homo_sapiens`, spatial | 78,894 | 43,386 | 1.82 | 35,508 |
| `mus_musculus`, spatial | 78,317 | 35,162 | 2.23 | 43,155 |
| non-spatial total | 266,560 | 201,439 | 1.32 | 65,121 |

**The spread is 1.02x to 2.23x, and the widest cases are the smallest stores.** The two spatial experiments assay targeted panels, so they see least of the reference, and they also hold the least data, so a wider `var` costs least there in absolute terms. Human, which carries most of the corpus, sits at 1.28.

**Assayed is not detected.** Every row in `var` has `n_measured_obs > 0`, without exception, so membership tracks assay coverage rather than expression. But some genes are assayed and never nonzero anywhere in the corpus: **119** in human, **2,868** in mouse, **1,014** in human spatial. Those rows carry no values at all, and without the presence matrix they are indistinguishable from genes nothing ever measured. That is spec section 2's measured-versus-zero requirement, visible at corpus scale rather than only per dataset.

**Most of `var` is the reference copied verbatim.** Of its seven columns, four (`feature_id`, `feature_name`, `feature_type`, `feature_length`) come straight from the pinned reference and depend on nothing else. Two, `nnz` and `n_measured_obs`, are aggregates over every cell in the corpus. The seventh is `soma_joinid`, retired by spec section 4. So the only parts of `var` that require knowing the whole corpus are those two aggregates, and they are exactly what a writer holding one dataset cannot produce.

**The cost of publishing the full reference is shape, not storage.** A gene assayed nowhere contributes no nonzeros, so it adds nothing to `data` or `indices`. What grows is the `var` table, by tens of thousands of rows, and the presence matrix, at one bit per dataset and gene. What it does change is band geometry, since bands are contiguous ranges of `var`.

**What it bears on.** Publishing the full reference makes a gene's index a function of `(reference version, feature_id)` alone, so adding a dataset never shifts another dataset's indices. Publishing the assayed subset makes indices a function of corpus membership. No decision is recorded here; the choice is open and is stated in spec section 3.3.

## 11. Validation methodology

The metrics, query corpus and pass conditions were fixed before the measurements of Sections 9 and 10 were taken, so those results are reported against criteria that were not chosen after seeing them. One criterion did change: 10.1's finding on wall clock, and 11.4 states the revised form.

### 11.1 v1-over-TileDB is the oracle

A correct implementation of every claim in this document already exists and keeps running throughout the transition: the v1 SDK reading the TileDB-SOMA corpus. Equivalence is therefore established **differentially**, by running the same query through both stacks and comparing results rather than against expectations written by hand. That is what makes validation meaningful at one to five datasets instead of requiring full corpus scale, and it is why the published TileDB releases have testing value beyond continuity: for as long as they remain available, every claim remains checkable.

The oracle expires at cutover. Until then, two conditions keep it trustworthy:

- **Every recorded result names the oracle version.** `tiledbsoma` is on an active release cadence and is pinned at 1.15.3 in the builder (11.7). A parity result measured against an unrecorded oracle version cannot be reproduced. The same applies to `anndata`, whose `~=0.12` pin floats across minor releases (11.7).
- **A disagreement is not automatically a v2 defect.** Where the stacks differ, the source datasets are the tiebreaker, since both derive from them. The oracle can be wrong: a discrepancy may be a pre-existing TileDB-path bug, which is a finding to report upstream rather than a v2 failure to fix.

### 11.2 Query corpus

A fixed set of queries, versioned alongside the harness, each earning its place by exercising something that could plausibly break independently.

**The baseline is the v1 Python SDK, not any downstream service.** What v2 must not regress against is `cellxgene_census.get_anndata()`: the call the API contract and the 5-year LTS commitment are made about (Sections 6, 7). Every SOMA-side figure in Sections 9 and 10 is that path: `readers.read_soma` builds an `axis_query` and calls `to_anndata()`, which is what `get_anndata` does internally, so the oracle of 11.1 and the performance baseline are the same code.

**The primary anchor is what the SDK's own documented usage does.** Filters taken from this repo's notebooks, with their measured fan-out over the 1,086 human datasets:

| Documented SDK filter | Cells | Datasets touched |
|---|---:|---:|
| `tissue_general == 'vasculature'` | 431,901 | 10 (1%) |
| `is_primary_data == True and tissue_general == 'lung'` | 6,167,731 | 45 (4%) |
| `tissue_general == 'lung'` | 10,181,125 | 92 (8%) |
| `tissue_general in ['brain', 'lung']` | 58,762,981 | 447 (41%) |

These are one- and two-column predicates over corpus-wide attributes, which matters twice over. The two-column cliff of 9.3 is not an artificial construction: `is_primary_data == True and tissue_general == 'lung'` is exactly that shape, it is in the documentation, and it costs 9.2x, and the fan-out range that decides store organization (10.5) is **10 to 447 datasets**, measured rather than inferred. The same measurement gives the per-column costs 9.3 rests on: `dataset_id`, `tissue_general` and `is_primary_data` together cost **4 MB to scan across 158,982,719 rows**, 0.025 bytes per row, because all three correlate with row order and compress to almost nothing. Against `nnz` at 1.6 bytes/row and `cell_type` at 0.33, predicate cost varies about 65x with *which* columns are named, and the expensive clauses are the numeric ones.

WMG sits in the corpus as what it is: a real filter shape from a real consumer, useful because it is the widest fan-out anyone has asked for and the only 4-column predicate available. It is not the baseline, and its own cost today is unknown: no WMG wall-clock baseline exists in recorded form, and 10.1 shows one taken from outside the region would have been unusable anyway.

| Query | What it exercises |
|---|---|
| The documented SDK filters above | The baseline the no-regression threshold is set against |
| WMG's production filter, verbatim below | A real 4-column predicate at maximum fan-out |
| Single-gene lookup across all cells | Worst-case sparse column access |
| Broad tissue slice | Large contiguous read |
| Cell-type slice | Selective read with a dictionary-encoded predicate |
| Identity-coordinate access | `soma_joinid` on v1, its replacement on v2 (spec section 4) |
| Presence-matrix read | The measured-vs-zero distinction (spec section 2) |
| Spatial read: scene, point cloud, image | The four object kinds with no non-spatial counterpart (spec section 2) |
| Multi-valued `disease` filter | Membership rather than string equality (spec section 2) |
| Symbol-collision `feature_name` lookup | Non-unique gene symbols (spec section 2) |
| `raw` and `normalized` on non-spatial; `raw` only on spatial | Layer presence varying by experiment (spec section 2) |

**The WMG filter, verbatim**, from `backend/wmg/pipeline/constants.py` in `single-cell-data-portal`, fully interpolated. Human and mouse use an identical string:

```
is_primary_data == True and nnz >= 500 and cell_type_ontology_term_id != 'unknown'
and tissue_general_ontology_term_id != 'UBERON:0001017'
and tissue_general_ontology_term_id != 'UBERON:0001007'
and tissue_general_ontology_term_id != 'UBERON:0002405'
and tissue_general_ontology_term_id != 'UBERON:0000990'
and tissue_general_ontology_term_id != 'UBERON:0001004'
and tissue_general_ontology_term_id != 'UBERON:0001434'
```

The field is `cell_type_ontology_term_id`, not `cell_type`, and the threshold is `nnz >= 500`. The read is `census_version="latest"`, two organisms sequentially, measurement `RNA`, **no `var` filter at all**, and `X["raw"]` consumed blockwise on axis 0 at `row_stride = 50_000` as CSR. Two inefficiencies in it affect fairness: the pipeline materialises the full filtered `obs` to pandas **twice**, in `expression_summary.py` and again in `cell_counts.py`, and passes no `column_names` projection, so every `obs` column is read. A comparison must either reproduce both or state plainly that it did not.

Two rules on coverage. **The corpus iterates experiment specs rather than hardcoding human and mouse**, because there are seven experiments (spec section 2) and a result covering two of them is not a result about the Census. And **the last four rows exist because they fail silently**: a filter that treats a multi-valued field as a string, or resolves a gene by a colliding symbol, returns a plausible subset and raises nothing.

**Matched pairs, for the benchmark rather than for parity.** Amplification is query-dependent, so ranking layouts needs pairs that vary one axis and hold the rest fixed, each sized to make any difference attributable to locality rather than volume:

| Pair | Varies | Status |
|---|---|---|
| `dataset_id == X` vs `cell_type_ontology_term_id == Y` | row locality at similar result size | Needs a third variant: its halves differ in predicate column count as well as locality (9.2) |
| all genes vs a ~100-gene subset, same cells | column locality | Cannot discriminate as written: sparse `X` has no gene axis (9.4). Becomes a measurement only when the *encoding* varies, which is 10.3 |
| 0.1% vs ~50% selective, same experiment | selectivity | Live |

### 11.3 Parity checker

Runs a query through both stacks and asserts:

1. **Identical cell selection**, compared by source `dataset_id` plus the source H5AD's own row identity, not by either stack's internal join key: comparing by `soma_joinid` would assume the mapping under test (spec section 4).
2. **Identical gene selection**, compared by `feature_id` (spec section 2).
3. **`X` values equal within tolerance**, for every layer the experiment carries. `raw` holds integer counts and must match exactly; `normalized` is floating point and is compared within a stated tolerance.
4. **`obs` and `var` columns equal**, including dtype and categorical level sets, since a dictionary-encoded column that loses its levels still compares equal value-by-value (9.5).

Parity is a precondition, not a tradeoff: a faster stack that fails parity fails the gate.

### 11.4 Benchmark rig

Reported per query, per dataset-count, and per layout candidate. **The metrics are ranked, not a flat list**, because wall clock from outside `us-west-2` varies by more than an order of magnitude on an identical query and cannot carry a comparison (10.1).

**Primary: deterministic counters.** Reproducible across machines and networks, so a difference between two layouts is a real difference.

| Metric | Why | How measured |
|---|---|---|
| Request count | Dominates cold reads (10.1), and the metric object-per-chunk layouts regress | Zarr: wrap the store and count `get`/`get_partial_values`. TileDB: `Context.VFS.read_ops_num` |
| Bytes read | Distinguishes a slow stack from an over-fetching one | Zarr: sum returned buffer lengths. TileDB: `Context.VFS.read_byte_num` |
| Corpus bytes and object count | Stored size against the baseline of 1.4 TiB / 325,695 objects (spec section 2) | Store enumeration |

On the TileDB side, count **LIST and HEAD operations too**: `Context.VFS.ls_num`, `is_object_num` and `file_size_num` are real S3 round-trips that `read_ops_num` excludes. Counting reads alone flatters TileDB.

**Secondary: wall clock, and only in-region.** Filter time and `X` slice throughput, split cold from warm, on an EC2 instance in `us-west-2` rather than a laptop. A laptop number is reported only as a sanity check on ordering.

Runs are reported at **1, 5, 25, and 100 datasets** so the output is a scaling curve rather than a single number, and any extrapolation to 1,845 datasets is stated with the curve it came from. The ladder takes the smallest datasets available, so it scales dataset *count* rather than bytes; the corpus median dataset is far larger than the ladder mean, and a large-dataset case is needed before the extrapolation is trusted. Selectivity is a dimension rather than an accident: each query records what fraction of the experiment it selects, and the corpus deliberately spans both ends.

**The threshold's baseline is `cellxgene_census.get_anndata()` on the published TileDB corpus**, measured in requests and bytes rather than seconds, against the documented SDK filters of 11.2. 9.3 gives it a concrete form: the current SDK answers a one-column scoped read in **0.4 MB and 9 requests** and a documented two-column read in **3.4 MB**. Those are the figures v2 is held to. The 258 MB an `nnz` clause costs is the worst case rather than the target: worth fixing, since WMG issues one, but a threshold set there would clear a bar users never approach.

**The tolerance band: no more than 2x the v1 figure on either metric, per query, and no regression at all on the corpus-level totals.** Two-fold rather than parity because v2 is a different engine with different metadata overheads, and because the counters are exact (10.1), so a 2x band is a real bound rather than noise absorption. Corpus bytes and object count get no allowance because 10.2 measured them as strictly better.

### 11.5 Dataset subsets

Canonical manifests at 1, 5, 25, and 100 datasets, each pinned by explicit `dataset_id` list rather than by selection criteria, so results are comparable across runs and across people. The subsets must include spatial datasets and at least one dataset from each of the five organisms; a random sample would draw almost entirely from human and mouse (spec section 2). The current manifests have a known defect: size-ascending selection concentrates secondary data, which disables the WMG query below rung 25 (9.1).

### 11.6 What this harness does not establish

- **Full-scale behaviour.** Every result below 1,845 datasets is an extrapolation. Section 13 addresses this and only a full build settles it.
- **Anything the sweep's synthetic `obs` axis implies about value distribution.** It is grown past 64,798 rows by repeating whole dataset blocks (10.2), which keeps `dataset_id` contiguous and cell types realistically distributed *within* a dataset, but fixes how widely a value spreads *across* chunks: the fraction of chunks a predicate must read is 0.339 at 5M rows and 0.340 at 20M. So sweep ratios hold at one distribution and no row count changes them, and that distribution is pessimistic: the real corpus puts the same predicate shape at 0.124 and a median-frequency one at 0.0030 (10.2). **Read sweep costs as an upper bound.** Occupancy itself is now measured on the real corpus, so the gap this bullet used to describe is closed; what is not closed is where the 4--5% pruning threshold sits on S3 rather than a local filesystem, since it trades requests against bytes and only their relative latency fixes the crossover (10.1).
- **Errors common to both stacks.** Differential testing cannot detect a mistake made identically on both sides, such as a shared misreading of the source schema.
- **Real cold-S3 latency at scale.** Small-subset cold reads approximate it; a full corpus in its published location does not exist until the first full build.

### 11.7 Toolchain versions, and why each result must record them

| Dependency    | Builder                | Reader (`cellxgene_census`)                   |
| ------------- | ---------------------- | --------------------------------------------- |
| `tiledbsoma`  | `==1.15.3`             | `>=1.15.3` (must not fall behind the builder) |
| `anndata`     | `~=0.12`               | unpinned                                      |
| `spatialdata` | `>=0.5.0` (test extra) | `>=0.5.0` (`spatial` extra)                   |

The builder already requires `anndata~=0.12`. The Zarr v3 / lazy-read generation this design depends on is already a census build dependency.

Two qualifications on that, both relevant to how much reassurance it offers:

- **The lazy APIs are experimental.** `read_lazy` and `read_elem_lazy` live in `anndata.experimental` and are documented as subject to change. Every result must therefore record the resolved `anndata` version, the same discipline Section 11.1 applies to the oracle.
- **But `read_elem_lazy` is already in production use in this repo**, at `build_soma/build_soma.py:475`, where the builder lazy-reads spatial images before handing them to a dask store. So the dependency on experimental lazy reads is not a new risk this design introduces; it is one the census build already carries.

Note what `~=0.12` actually permits. Under PEP 440 the two-component compatible-release form means `>=0.12, <1.0`, **not** `>=0.12, <0.13`; only `~=0.12.0` would restrict to the 0.12 series. So the builder already accepts every anndata minor release below 1.0, and a resolution performed on 2026-07-29 selects **0.13.2**.

The consequence runs opposite to the obvious reading: the builder's anndata version is not pinned in any practical sense, it floats across minor releases. Two builds months apart can therefore be produced by different anndata minor versions with different lazy-read behaviour, and nothing records which was used. The version-recording discipline Section 11.1 imposes on the oracle applies to the builder for the same reason.

## 12. Identity transition proof

`TBD` until the mapping artifact exists. A real contributed embedding from `cellxgene-contrib-public`, mapped through the identity model (spec section 4) and re-joined against the new corpus, landing on the same biological cells it did under v1.

### 12.1 What the proof has to assert

The identity model (spec section 4) is a design, and its claim is falsifiable: an embedding row computed against v1 must, after mapping, attach to the same cell in v2. The gate is a differential test with the same shape as the parity checker (11.3), and the same discipline applies: it must be possible for it to fail.

Four assertions, in the order they can go wrong:

1. **The mapping is total over the release.** Every `soma_joinid` in the v1 experiment appears exactly once in the 4.4 mapping. A mapping that silently omits cells produces an embedding that loads and is short.
2. **The mapping is injective.** No two `soma_joinid` values map to the same `(dataset_id, observation_joinid)`. 4.1 verifies the pair is unique in the corpus, which makes this a check on the artifact rather than the data.
3. **A mapped embedding lands on the same cells.** Take an embedding's `soma_joinid` coordinates, map them, resolve them in v2, and the resulting cells must equal the v1 cells by source identity. This reuses the parity checker's cell key directly, since it is already `(dataset_id, observation_joinid)`.
4. **An unmappable identity raises.** A `soma_joinid` absent from the mapping, a cell dropped between releases, must produce an error, not a silently shorter result (spec section 4.3).

The paired mutations that prove the gate is live: drop one row from the mapping (1 must fail), duplicate one target pair (2 must fail), shift the mapping by one position (3 must fail), and pass a `soma_joinid` beyond the experiment's range (4 must raise).

### 12.2 What is available to test with

`cellxgene-contrib-public` publishes contributed embeddings with a `census_version` in their metadata, and `get_all_available_embeddings` enumerates them, so a real embedding pinned to the 2025-11-08 release can be selected without guessing. The mapping artifact itself does not exist yet, it is produced by the builder (spec section 4.3), so the proof needs either a builder change or a mapping derived from the published corpus, which is the same two-column projection and costs 1--19 MB per small experiment (spec section 4.3). Choosing the smallest experiment that has a contributed embedding keeps the whole proof laptop-scale.

## 13. Build resource envelope

`TBD` until the builder emits Zarr. Memory and runtime as dataset count grows, and the projection to full scale that sizes the build host.

### 13.1 What is already known, from the spike's own conversions

The spike converted census slices to Zarr repeatedly, and those conversions are the same operation the builder performs, at small scale: 565 cells and a 2.1 MB store at negligible peak memory, 5,716 cells and 18.3 MB likewise, and 64,798 cells and 213.6 MB at roughly **1.2 GB resident**.

The conversion path holds the whole matrix in memory as CSR before writing (`convert.py`), so peak memory scales with nonzeros rather than with the store. Extrapolating the rung-25 figure on nonzeros alone puts a 340,000-cell rung at roughly 4 GB and the full human experiment far beyond a laptop.

Two things this establishes and one it does not. It establishes that **the builder must stream rather than materialise**, since a whole-experiment CSR is not resident-memory-shaped, and that **the banded encoding of 10.3 makes streaming natural**, because each gene band is written independently and a band is a bounded slice of the nonzeros. The `X` padding of 10.5 helps here rather than complicating it: since every dataset's nonzero range begins on a chunk boundary, the builder can append one dataset at a time without ever reading back or rewriting a chunk another dataset shares. It does not establish the envelope: these are conversions of a corpus that already exists, not builds from source H5ADs, and the builder additionally holds harmonization state this path never touches.

`TBD`: the actual envelope, which needs the builder emitting Zarr rather than the spike converting it.

---

## 14. Findings that need a decision elsewhere

The decisions themselves are listed in spec section 10, which is the authoritative set. This section holds what was found while measuring: defects in existing documentation, and two questions only the census team can answer.

### Defects in existing documentation

In `docs/cellxgene_census_schema.md`:

1. The obs-to-scene mapping is named `obs_spatial_presence` in the implementation and at line 952, but `obs_presence_matrix` at lines 998 and 1033. The latter object does not exist.
2. Line 935 states spatial experiments live in `census_obj["census_data"][organism]`; should be `census_spatial_sequencing`.

In `dev-docs/cellxgene_census_storage_and_release_policy.md`, five points of drift, tabulated in spec section 5. The filename and bucket are wrong, the `V#` alias requirement is unmet in production, the `stable` alias is undocumented despite being the SDK default, and `do_not_delete` and `relative_uri` are undocumented. This document is the one spec section 5 amends, so reconciling it is part of that work rather than separate from it.

On the `V#` alias specifically, which of the five is least obvious how to settle: the policy states one "MUST exist that points to the release" for every LTS release, six releases carry `flags.lts`, and no such alias exists. Nothing consumes it either. The SDK resolves aliases generically, `stable` is its default, and the docsite directs users to pin by build date. So implementing it, deleting the requirement, or downgrading it to `MAY` all appear to cost nothing, which is what makes it cheap to settle while the document is already open.

### Questions for the census team

1. **Who approves an amendment to the storage and release policy?** `dev-docs/cellxgene_census_storage_and_release_policy.md` must be amended before a Zarr corpus can be published (spec section 5). The file records "Last edited: April, 2023" and has no substantive commit in the repository's visible history; the only change touching it is a link fix. Sustained recent ownership of the governing documents in `dev-docs/` and of the Census schema sits with **Brian Raymor** (7 of the last 10 `dev-docs/` commits, plus the 2.4.0 schema work through Nov 2025), which is the best available signal but not a confirmation.
2. **Is** `soma_joinid` **intended to be stable across releases?** It effectively is today, which is why users join on it, but the guarantee does not appear to be documented. This materially changes the identity model (spec section 4), and is decision 5 there.

---

## 15. The current builder's harmonization catalogue

This backs the table in v2 section 6. Every transformation the builder applies between curated H5ADs and the published corpus, read from `tools/cellxgene_census_builder` on 2026-08-10; paths below are relative to its `src/cellxgene_census_builder/`. Each step is tagged **per-dataset** (needs only one dataset) or **corpus-wide** (needs more than one). Build mechanics with no effect on data, logging, parallelism, TileDB plumbing and consolidation, are excluded. Phase 4 of the plan is done only when every step here is implemented by the per-dataset writer, moved to release assembly or the catalog build, or explicitly dropped.

Two findings stand above the list. **There is no external gene reference**: the corpus `var` is `pd.concat(per-dataset vars).drop_duplicates()`, the union of eligible genes in first-seen order across filesize-sorted datasets (`build_soma/experiment_builder.py:403-405`), so the total order depends on which datasets are in the build and in what order, which is the renumbering v2 section 4 removes. And **no step deduplicates cells across datasets**: `is_primary_data` is a pass-through column used only to count (`build_soma/experiment_builder.py:458`), and no cross-dataset uniqueness check exists outside the validation suite's coordinate check.

One update since the catalogue was read: the builder gates on dataset schema 7.0.0 (step 1), and curation schema 7.1.0 has since shipped, with every live dataset re-versioned onto it (all 2,216 as of 2026-08-19). 7.1.0 adds perturbation fields, `perturbation_types` and `genetic_perturbation_strategy` at the dataset level plus per-cell experimental-condition and genetic-perturbation columns, with CRISPR screen strategies among the allowed values. Verified live 2026-08-19: one dataset populates `perturbation_types` (chemical), none populate `genetic_perturbation_strategy` yet. These columns are the first concrete case for v2's rule (section 3) that the writer's obs projection and the catalog's column set track the curation schema rather than freeze at this catalogue's list, which the builder's fixed projection (step 13) would silently drop.

### 15.1 Admission, all corpus-wide

1. **Manifest and schema gate.** Datasets from the Discover API filtered to `schema_version == CXG_SCHEMA_VERSION` (7.0.0), one H5AD asset required; a CSV manifest can substitute. `build_soma/manifest.py:42-104`, `build_soma/globals.py:15`.
2. **Blocklist and duplicate collapse.** Blocked `dataset_id`s dropped (a missing blocklist URI is a hard error); duplicate `dataset_id`s collapse to one. `build_soma/manifest.py:107-175`.
3. **Filesize ordering.** Datasets sorted ascending by `asset_h5ad_filesize`; this sort is what ultimately fixes the total order of obs and var joinids. `build_soma/build_soma.py:207-217`.
4. **Utilized filter.** Post-filter cell counts summed across experiments; zero-cell datasets are excluded and their staged H5ADs deleted. `build_soma/build_soma.py:261-282,151-161`.
5. **Dataset `soma_joinid` assignment**, 0..n-1 in the sorted order. `build_soma/datasets.py:63-65`.
6. **`census_info/datasets`.** Utilized datasets projected onto the 11-column table spec with the post-filter cell count. `build_soma/datasets.py:68-85`, `build_soma/globals.py:18-33`.

### 15.2 Reading and filtering, all per-dataset

7. **Staging.** Each H5AD copied locally as `{dataset_id}.h5ad` with a size check. `build_soma/source_assets.py:36,61`.
8. **Raw matrix selection.** If a `raw` group exists the corpus uses `raw/X` and `raw/var`, so the Census always carries raw counts. `build_soma/anndata.py:216-267`.
9. **Per-file enforcement.** `uns/schema_version` must equal the pinned version or the build raises; the AnnData encoding version is asserted too. `build_soma/anndata.py:240-252`.
10. **Duplicate-coordinate coalescing.** Non-canonical sparse X gets `sum_duplicates()` on slice. `build_soma/anndata.py:127-140`.
11. **Cell filter.** Keeps `tissue_type in {tissue, organoid}` and the experiment's allowed assays; requires `uns/organism_ontology_term_id` to match the experiment or the whole file is dropped; spatial experiments additionally require `is_primary_data == True`. `build_soma/anndata.py:311-366`, `build_soma/experiment_specs.py:20-108`, `build_soma/globals.py:211-333`.
12. **Gene filter.** Keeps `feature_biotype == "gene"` and `feature_reference ==` the experiment organism, so cross-organism genes are dropped. `build_soma/anndata.py:335,350`.
13. **`obs` projection.** Only the spec's term columns read (19, plus 3 spatial); dtypes normalized to nullable pandas types so datasets missing a column still concatenate; `dataset_id` stamped. `build_soma/experiment_builder.py:356-383`.
14. **`var` projection.** Reduced to `feature_id, feature_name, feature_length, feature_type`, index renamed to `feature_id`. `build_soma/experiment_builder.py:377-381`.

### 15.3 The corpus axes, all corpus-wide per experiment

15. **`obs` concatenation.** Per-dataset frames concatenated in dataset order, so each dataset's cells are contiguous. `build_soma/experiment_builder.py:397-411`.
16. **`var` union.** The first-seen-union finding above; dedup is on the full 4-tuple, so a `feature_id` whose name, length or type disagrees between datasets makes two rows, caught only by validation. `build_soma/experiment_builder.py:403-405`, `build_soma/validate_soma.py:288-290`.
17. **Cell and gene `soma_joinid` assignment**, `range(len)` over the concatenated frames. `build_soma/experiment_builder.py:445-446`.
18. **`tissue_general` derivation.** Each `tissue_ontology_term_id` mapped to a high-level UBERON term via `cellxgene_ontology_guide`; unmappable IDs fall back to themselves. Computed on the concatenated frame but each value depends only on its own cell. `build_soma/experiment_builder.py:824-842`, `build_soma/tissue_mapper.py:15-37`.
19. **Sentinel placeholders.** Missing obs columns filled with int64 min or NaN, later overwritten by the stats merge; missing var numerics filled with 0. `build_soma/experiment_builder.py:423-428,450-454`.
20. **Joinid bookkeeping.** `dataset_obs_joinid_start`/`dataset_n_obs` per dataset, and the `feature_id` to var-joinid map used to remap every dataset's genes. `build_soma/experiment_builder.py:459-463`.
21. **Contiguity assertion.** Each dataset's obs joinids must form one contiguous block. `build_soma/experiment_builder.py:432`.
22. **Categorical handling on write.** Dictionary-spec fields cast to `category` with the narrowest index type the observed codes allow; observed Arrow types asserted against the spec. `build_soma/schema_util.py:122-175`, `build_soma/globals.py:84-150`.
23. **Dataframe write with explicit domain**, index `soma_joinid`, domain the observed min and max. `build_soma/experiment_builder.py:227-286`.

### 15.4 Expression and derived values

24. **Gene remap into corpus space**, per-dataset, but against the corpus-wide var map of step 20; an unmappable gene fails the build. `build_soma/experiment_builder.py:593-594,679`.
25. **Single-modality assertion** per dataset, no mixing spatial and non-spatial assays. `build_soma/experiment_builder.py:596-599`.
26. **X coercion and validation**, per-dataset. Dense to CSR, zeros eliminated, values to float32, and `is_nonnegative_integral` asserted, so non-integer or negative raw values abort. `build_soma/experiment_builder.py:612-640`, `build_soma/util.py:25-36`.
27. **Per-cell statistics**, per-dataset. `raw_sum`, `nnz`, `raw_mean_nnz` (non-finite to 0.0), `raw_variance_nnz` (over stored values, ddof=1). `build_soma/stats.py:12-38,63-131`.
28. **`n_measured_vars`**, per-dataset: the size of that dataset's filtered var, identical for all its cells. `build_soma/experiment_builder.py:529-550`.
29. **Per-gene `nnz` and `n_measured_obs`**, corpus-wide sums of per-dataset reductions; each dataset contributes its full post-filter cell count to every gene it assayed. `build_soma/stats.py:41-60`, `build_soma/experiment_builder.py:513-521,542,801-805`.
30. **Stats merge** into the corpus obs frame by joinid, replacing the step 19 sentinels. `build_soma/experiment_builder.py:798-800`.
31. **`X["raw"]` write**, per-dataset values at corpus coordinates, rows shifted by the dataset's joinid start. `build_soma/experiment_builder.py:661-683`.
32. **`X["normalized"]` derivation**, per-dataset, non-spatial only. Values for cells whose assay is in the 21-term `FULL_GENE_ASSAY` list are divided by `feature_length` first; every row is divided by its sum; precision is reduced to 15 mantissa bits by round-half-to-even; resulting exact zeros are replaced with the smallest float32 subnormal so no stored value is zero. `build_soma/experiment_builder.py:601-607,666-693,862-899`, `build_soma/globals.py:336-358`.
33. **Presence matrix**, corpus-wide, assembled from per-dataset rows. True iff the gene appears in the dataset's filtered var, assayed rather than nonzero; False is implicit. `build_soma/experiment_builder.py:310-341,543-549`.
34. **Spatial extraction**, per-dataset. A `Scene` per dataset with a pixel coordinate space, `obsm/spatial` as a point cloud joined to corpus joinids with radius from `spot_diameter_fullres`, and the `hires` image as a multiscale CHW array with its scale transform; plus the corpus-wide `obs_spatial_presence` concatenation. `build_soma/build_soma.py:313-441,467-509,417-441`.

### 15.5 Corpus tables and metadata, all corpus-wide

35. **Per-category counts.** For 8 categories plus an `all` row, counts pivoted on `is_primary_data` into `unique_cell_count` and `total_cell_count`, accumulated dataset by dataset, then concatenated across experiments and summed by `(organism, category, ontology_term_id)`; spatial and non-spatial rows for one organism merge here. `build_soma/summary_cell_counts.py:14-131`.
36. **`census_info/summary`.** Schema version, build date, dataset schema version, total and unique cell counts summed over every experiment. `build_soma/census_summary.py:14-38`, `build_soma/experiment_builder.py:817-821`.
37. **`census_info/organisms`**, one row per organism with spatial and non-spatial collapsed. `build_soma/census_summary.py:41-66`.
38. **Provenance metadata**, `created_on` and `git_commit_sha` on the root collection. `build_soma/build_soma.py:164-179`, `build_soma/util.py:39-52`.
39. **Bounding-box metadata scrub**, a post-write deletion of `soma_dim_*_domain_*` keys working around TileDB-SOMA#1969. TileDB-specific; nothing to carry forward. `build_soma/build_soma.py:533-565`.

### 15.6 The validation suite, `build_soma/validate_soma.py`

Per-dataset passes: the corpus obs slice equals the filtered source obs column for column (`:304-407`); `X["raw"]` equals the source matrix on assayed genes and is empty elsewhere, with all five obs statistics recomputed and `nnz > 0` asserted per cell (`:685-866`); the presence row equals the intersection of corpus `feature_id`s with the source var (`:580-682`). Corpus-wide passes: joinids unique and exactly `[0, n)`, the var `feature_id` set matched against an independently computed union (`:232-301`); no duplicate coordinates in either X layer, the closest thing to a cross-dataset duplicate-cell check (`:528-577`); the normalization recomputed with row sums asserted near 1 and no stored zero in either layer (`:410-525`); `dataset_id` and `dataset_version_id` unique in the written manifest (`:897-909`).

What this establishes for v2: the corpus-wide steps are exactly the joinid system (dropped, v2 section 5), the `var` union (replaced by the pinned reference, v2 section 4), sums and concatenations over per-dataset reductions (release assembly, v2 section 6), and admission bookkeeping (the catalog build, v2 section 3). No corpus-wide step touches expression data except through a per-dataset reduction, which is what makes a build-free corpus possible.

---

## 16. The hosted embeddings: inventory and generation mechanism

Probed live 2026-08-17: `aws s3 ls --no-sign-request` against `s3://cellxgene-contrib-public/contrib/cell-census/soma/`, joined with the SDK's manifest at `contrib.cellxgene.cziscience.com/contrib/cell-census/contributions.json` (the source `get_embedding_metadata_by_name` reads).

### 16.1 Inventory

23 embedding arrays across four release prefixes, 2,995 GB in total: 2,290 GB of embedding tables and 705 GB of `indexes/` prefixes, the embedding-search indexes, present for the first three releases only.

| Release | Id | Model | Organism | Cells | Dims | Size |
| --- | --- | --- | --- | --- | --- | --- |
| 2023-12-15 | CxG-contrib-1 | scgpt | human | 62,998,417 | 512 | 46.5 GB |
| 2023-12-15 | CxG-contrib-2 | uce | human | 62,996,451 | 1280 | 116.9 GB |
| 2023-12-15 | CxG-contrib-3 | uce | mouse | 5,684,805 | 1280 | 10.5 GB |
| 2023-12-15 | CxG-contrib-4 | nmf | human | 28,534,015 | 200 | 3.8 GB |
| 2023-12-15 | CxG-contrib-6 | nmf | human | 60,664 | 200 | <0.1 GB |
| 2023-12-15 | CxG-czi-1 | geneformer | human | 62,998,417 | 512 | 46.7 GB |
| 2023-12-15 | CxG-czi-4 | scvi | mouse | 5,684,805 | 50 | 0.4 GB |
| 2023-12-15 | CxG-czi-5 | scvi | human | 62,998,417 | 50 | 4.8 GB |
| 2024-07-01 | CxG-contrib-7 | scgpt | human | 74,322,510 | 512 | 53.9 GB |
| 2024-07-01 | CxG-czi-6 | geneformer | human | 74,322,510 | 512 | 54.4 GB |
| 2024-07-01 | CxG-czi-7 | scvi | mouse | 41,233,630 | 50 | 3.2 GB |
| 2024-07-01 | CxG-czi-8 | scvi | human | 74,322,510 | 50 | 5.6 GB |
| 2025-01-30 | CxG-czi-9 | scvi | mouse | 42,776,863 | 50 | 3.3 GB |
| 2025-01-30 | CxG-czi-10 | scvi | human | 106,118,167 | 50 | 8.0 GB |
| 2025-01-30 | CxG-czi-11 | geneformer | human | 106,118,167 | 512 | 77.6 GB |
| 2025-01-30 | CxG-czi-12 | tf-sapiens | human | 106,118,167 | 2048 | 316.7 GB |
| 2025-01-30 | CxG-czi-13 | tf-exemplar-human | human | 106,118,167 | 2048 | 316.2 GB |
| 2025-01-30 | CxG-czi-14 | tf-exemplar-mouse | mouse | 42,776,863 | 2048 | 127.5 GB |
| 2025-11-08 | CxG-czi-15 | scvi | human | 158,982,719 | 50 | 12.0 GB |
| 2025-11-08 | CxG-czi-16 | scvi | mouse | 43,653,561 | 50 | 3.4 GB |
| 2025-11-08 | CxG-czi-17 | tf-sapiens | human | 158,982,719 | 2048 | 474.5 GB |
| 2025-11-08 | CxG-czi-18 | tf-exemplar-human | human | 158,982,719 | 2048 | 473.9 GB |
| 2025-11-08 | CxG-czi-19 | tf-exemplar-mouse | mouse | 43,653,561 | 2048 | 130.1 GB |

What the inventory establishes:

- **The maintained embeddings were regenerated for every LTS release while builds ran.** The newest set is pinned to 2025-11-08 and covers 158,982,719 human and 43,653,561 mouse cells, the full corpus for those organisms at the final build. Geneformer's last run is 2025-01-30 (106,118,167 cells).
- **The community contributions are frozen earlier**: scgpt at 2024-07-01, uce and nmf at 2023-12-15.
- **A one-time conversion of the newest table per model is about 1.36 TB** (all four releases: about 2.29 TB), dominated by the two 2048-dimension tables at roughly 950 GB between them. Each array is at most 9 objects, so object counts are negligible against the corpus's 325,695.

### 16.2 The generation mechanism, in this repository

- `tools/models/{scvi,geneformer,transcriptformer}` are the pipelines; the `tf-` tables are the transcriptformer pipeline's output. The scVI pipeline is three scripts: `scvi-prepare.py` (primary cells with `nnz >= 300`, top 8,000 highly variable genes, a batch key concatenating `dataset_id`, `assay`, `suspension_type`, `donor_id`), `scvi-train.py`, `scvi-create-latent-update.py`.
- **The per-release rerun is inference, not retraining.** `scvi-create-latent-update.py` pushes all cells of a new Census version through `load_query_data()` on the existing trained model; the README states this is "required so that the model doesn't need to be re-trained from scratch on each Census version".
- `tools/census_contrib` validates and packages both maintained and community embeddings into the bucket's SOMA arrays; `tools/census_embeddings_indexer` builds the `indexes/` prefixes.

Two boundaries on the dataset-local claim: the scVI batch key includes `donor_id` and the highly-variable-gene selection was computed corpus-wide, so embedding a new dataset under a frozen model is dataset-local, while refreshing the model itself is not. And the embeddings describe one release's dataset versions, so a dataset re-curated since an embedding's release may be covered only partially; a per-dataset layout should record the covered `dataset_version_id`.

## 17. The catalog prototyped over the validation-harness fixtures

Until 2026-08-25 the catalog (v2 section 3) had never been constructed at any scale: the data-lake probe of
10.5 opened stores without one, and the plan's Phase 5 records that no surveyed prototype carries one. This
section is the first construction, at fixture scale, with a correctness proof. Code:
`tools/validation_harness/catalog_proto.py`, run against the 19 per-dataset Zarr fixture stores the harness
derives from the pinned 2025-11-08 release (three smallest datasets per experiment, one for chimp).

**Built as specified.** One row per store carrying all 18 indexed columns (every one present in the census
`obs`), packed bitmaps over a 361-value dictionary, the five numeric ranges over the per-cell statistics, and
the `is_primary_data` tri-state; written as Parquet (24,964 bytes) beside its dictionary (7,115 bytes). The
bitmap width came to 54 bytes per dataset, which is the dictionary size divided by eight plus per-column
padding; v2's corpus figure of 479 bytes at 3,834 values is the same arithmetic, so the sizing formula is
confirmed even though the corpus-scale total is not measured here.

**The correctness rule holds on every harness query.** The pruner excludes a dataset only when the catalog
proves no row can match: bitmap membership for equality and in-lists, with `' || '` splitting on the
multi-valued disease pair so membership semantics cannot wrongly exclude; `!=` excludable only when the named
value is the only value the dataset holds; numeric comparisons against the ranges; the tri-state for
`is_primary_data`; exact identity for `dataset_id`; and anything unproven keeps the dataset. Answering all 55
obs-bearing corpus queries (11.2) from only the candidate stores, every result compares identical to the v1
oracle snapshot under the 11.3 comparator, WMG's four-column filter and the silent-failure traps included.

**The check can fail, which is what makes it a check.** Deliberately removing one true candidate from a pruned
set is reported by the comparator as a `strict-subset` finding, the exact shape a pruning bug produces. Phase
5's done-when, no false negatives proven by running both pruned and full reads over the harness queries, is
therefore mechanically checkable today.

**What this does not establish.** Pruning power: with 19 stores and dataset-clamped queries, the 45-of-149
stores-opened figure is illustrative only. Nor the 2.8 MB corpus-scale catalog size, the dictionary built from
the full corpus, or build cost. And the builder here reads each store's `obs` directly because fixture stores
carry no writer summary record; v2's builder reads summary records, which arrive with the Phase 4 writer.

## 18. Findings from building the validation harness

Three facts surfaced while building Phase 2 (harness code in `tools/validation_harness/`, branch
`spike/validation-harness`), each recorded here because it outlives the harness that found it.

### 18.1 The v1 stack's dtype for nullable integer columns depends on the selection

Probed 2026-08-25 during the first two-stack parity run. The spatial experiments' `array_row`, `array_col` and
`in_tissue` columns are integers that are null for non-Visium datasets (a Slide-seqV2 dataset carries them as
all-null). Arrow-to-pandas conversion in the v1 read path yields `int64` when the selected rows happen to
contain no nulls and `float64` with NaN when they do, so the same column's dtype differs between two v1 queries
over the same experiment. Dtype equality on these columns is therefore not a property the baseline has, and the
harness comparator tolerates int/float disagreement when values compare equal (NaN equal to NaN), with the
policy stated in `comparator.py`. Two consequences: the parity checker's "obs columns equal, including dtype"
(11.3) is amended by that stated policy for nullable integer columns; and a v2 SDK that returns a stable
nullable dtype (for example Int64 always) will legitimately differ from v1 here, a filter-translator and
result-schema case Phase 6 has to decide rather than inherit.

### 18.2 TileDB's VFS counters undercount wire traffic

Validated 2026-08-25 by placing a logging proxy between `tiledbsoma` and S3 (`vfs.s3.endpoint_override`;
`tools/validation_harness/trace_check.py`). On the probe query (the chimp experiment's `var`, one column), the
VFS counters reported 46 read ops and 3,417 bytes while the proxy logged 52 GETs and 35,618 body bytes: the
counters are contained in, not equal to, the wire traffic, excluding some schema and listing requests. The
same-boundary Zarr counter matched its access log exactly (3,386 counted, 3,386 logged, against a
Range-supporting server; the stock `python -m http.server` ignores Range headers and corrupts sharded reads).
Consequence for this document: the SOMA request and byte figures in sections 9 and 10 were produced by the same
VFS counters and carry the same downward bias; cross-stack request comparisons flatter the SOMA side. The bias
is recorded with each benchmark run in `bench_soma_config.json`.

### 18.3 anndata 0.13 exposes X as `layers[None]`

An announced 0.13.0 behavior change (release notes, 2026-07-07: "None appears in `.layers.items()`/`.keys()`"),
met in practice 2026-08-25: code iterating `layers.items()` receives X itself under the key `None`, and the
harness serializer written against the 0.12 mental model wrote every X twice, which the corruption suite caught
and byte-identity verification could not (both runs of the same code duplicated identically). The guard is two
lines in `canonical.py`. Anything in the migration that iterates layers, the Phase 4 writer above all, needs the
same guard, and the v1 SDK's floating `anndata~=0.12` pin (11.7) resolves to 0.13.x today, so the behavior is
already live in freshly built environments.
