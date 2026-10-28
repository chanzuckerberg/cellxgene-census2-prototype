# census2-compact

Compaction of small Census 2 stores, as a Nextflow pipeline. `README.md` covers what it does and
how to run it. The repository root `AGENTS.md` gives the Census 2 tentpoles; compaction is a
secondary materialization beside them, and the per-dataset stores stay the primary copy.
`bulk-ingest-nf/AGENTS.md` records the Nextflow and HealthOmics sharp edges, and every one of them
applies here; this pipeline follows the same patterns. This file records what is particular to
compaction.

## Seams

`census2_dataset_writer.compaction` decides what a compacted store holds; its module docstring
gives the layout and the memory contract. `compact-group.py` handles the I/O: where the store goes,
the check, and the receipt. `tools/compaction-plan.py` decides the groups, outside the pipeline,
from a catalog. The catalog builder (`census2-sdk/tools/build_catalog.py --compacted`) checks the
compacted stores against the per-dataset stores and names them in the catalog. No receipt reaches
the catalog.

The pipeline runs in the runtime image of `bulk-ingest-nf`. A change to the compaction module needs
a rebuild of that image, as a change to the writer does.

## Sharp edges

`compact` writes `X/data` and `X/indices` one whole shard at a time. A write to part of a shard
reads the shard, merges, and writes it again; anndata's append of each dask block does that. A
shard written once, whole, is written once.

The source `obs` frames are held until the union frame is written, then dropped. Memory therefore
follows the cell count of the group. On the 103-dataset pilot, groups of about 1 GiB and 190,000
to 340,000 cells peaked at 0.8 to 1.4 GiB, and the group of 33 members peaked highest. Each member
reads its own copy of the 78,894-row `var`. `compact` now keeps one copy after the check of the
members, which is not yet measured. The
`memory` directive requests 4 GiB or 2 GiB plus 2 GiB for each million cells, whichever is more.
The planner bounds a group by bytes, not cells, so a group of many small datasets with few stored
values for each cell can have more cells than its bytes suggest.

A plan names stores by URI. A local plan names local paths, which a task container sees only
through a mount; `tests/nextflow.config` adds one. On HealthOmics, the stores are on S3.
