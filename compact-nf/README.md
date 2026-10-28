# census2-compact

A Nextflow pipeline that compacts small Census 2 stores: it writes one compacted store for each
group of per-dataset stores. It runs locally under Docker, and on AWS HealthOmics.

A compacted store is a materialization, a second copy. The per-dataset stores stay the primary
copy, and stay available. The catalog builder names a compacted store in the catalog rows of its
datasets, and the SDK then reads those datasets from it: one store opens for several datasets.
Exclusion stays per dataset, so a query opens a compacted store only when one of its datasets can
match. The `census2-dataset-writer` package decides what a compacted store holds, in
`census2_dataset_writer.compaction`; this pipeline is the I/O around it.

A run starts from a plan that `tools/compaction-plan.py` prepares from a catalog, one group per
line. The pipeline runs one task per group: read the member stores, write the compacted store and
a receipt to `--outdir`, and with `--check`, compare the compacted store with each member. A final
task reads the receipts and reports on the run.

## Requirements

The requirements of [`bulk-ingest-nf`](../bulk-ingest-nf/README.md): Nextflow 25.10.0, Docker, and
AWS credentials in the environment. The pipeline runs in the runtime image of `bulk-ingest-nf`,
which holds the `census2-dataset-writer` package. No other image is necessary.

## Plan

Build a catalog of the per-dataset stores, then plan from it:

```bash
python census2-sdk/tools/build_catalog.py s3://my-bucket/stores/ catalog.parquet

compact-nf/tools/compaction-plan.py catalog.parquet \
    --small-mib 256 --target-mib 1024 --strategy size --out plan.jsonl
```

The planner packs the stores of less than `--small-mib` MiB into groups of at most
`--target-mib` MiB. A larger store stays alone. The members of a group have one organism, one gene
reference version, and one spatial value. `--strategy size` packs first-fit decreasing by size,
for the fewest groups. `--strategy tissue` packs in the order of each store's main tissue, so
that the stores of one tissue share a compacted store. The planner lists each store for its size;
`--sizes sizes.json` keeps the sizes for a second plan.

All members of a group must share one access grant, because a compacted store has one. Plan only
from a public catalog. A private dataset stays in its per-dataset store.

## Compact

```bash
docker build -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .

nextflow run compact-nf -profile docker \
    --plan plan.jsonl --outdir s3://my-bucket/compacted --check
```

Then build the catalog again, with the compacted stores:

```bash
python census2-sdk/tools/build_catalog.py s3://my-bucket/stores/ catalog.parquet \
    --compacted s3://my-bucket/compacted/
```

The builder checks each compacted store against the per-dataset stores, and refuses the catalog
when a check fails: see `census2-sdk/README.md`.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--plan` | none, required | JSONL plan from `tools/compaction-plan.py`. Local path or `s3://` URI. |
| `--outdir` | `./out` | Local path or `s3://` prefix for the compacted stores and the receipts. |
| `--check` | `false` | Compare each compacted store with each member store, value for value. A failure goes into the receipt, and `REPORT` fails the run. |
| `--report_dir` | none | Where the run report publishes. The `healthomics` profile sets it. |

## Output layout

```
<outdir>/<group>.zarr/    the compacted store
<outdir>/<group>.json     the receipt
```

The name of a group comes from its members: `c-` and 16 hexadecimal digits of the SHA-256 of the
sorted member dataset versions. One set of members has one name. A rerun overwrites the store, as
in `bulk-ingest-nf`.

A task that fails twice is ignored. Its datasets stay in their per-dataset stores, which is a
correct catalog, so the run continues.

## AWS HealthOmics

As for `bulk-ingest-nf`: build the bundle, push a `linux/amd64` runtime image to private ECR, then
create a workflow and start a run. Pass `--outdir` as `s3_//bucket/prefix`.

```bash
compact-nf/tools/build-healthomics-bundle.sh compact-bundle.zip

aws omics create-workflow --engine NEXTFLOW --main main.nf \
    --definition-zip fileb://compact-bundle.zip \
    --parameter-template file://compact-nf/tools/healthomics-parameters.json \
    --storage-type DYNAMIC --name census2-compact
```

The run role reads the per-dataset stores and writes `--outdir`.

## Tests

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:$PWD" -w "$PWD" census2-bulk-ingest:dev \
    python tools/seed-stores.py tests/stores   # from compact-nf/

nf-test test --tag validation,compact
```

`tools/seed-stores.py` writes three small stores in the store format and a plan of one group. The
mount keeps the paths of the stores the same inside and outside the container, because the plan
names them.
