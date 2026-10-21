# census2-bulk-ingest

A Nextflow pipeline that converts CELLxGENE Discover datasets into Census 2 stores in bulk: one
AnnData-Zarr store per dataset. It runs locally under Docker, and on AWS HealthOmics.

A run starts from a manifest that `tools/cxg-manifest.py` prepares from the Discover API, one
record per dataset. The pipeline runs one task per dataset: fetch the source h5ad from Discover's
asset bucket on S3, build the store, and write the store and a receipt to `--outdir`. A final task
reads the receipts and reports on the run.

The [`census2-dataset-writer`](../census2-dataset-writer/README.md) package decides whether a
dataset is admitted and what a store contains; the runtime image installs it as a dependency. This
pipeline is the I/O around it. [`AGENTS.md`](AGENTS.md) records how the pipeline is put together.

## Requirements

- Nextflow 25.10.0, the version that the pipeline pins (`manifest.nextflowVersion`). HealthOmics
  runs that version, and CI tests it. Another version gives a warning.
- Docker.
- AWS credentials that can read the source bucket.
- Internet access for `tools/cxg-manifest.py` and `tools/stage-refdata.py` only. The pipeline
  itself needs none.
- `aria2c` for `tools/stage-refdata.py`, when you run it outside the runtime image
  (`brew install aria2`).

**Credentials in the container.** Tasks read credentials from the environment. A profile that
resolves credentials through a helper, such as SSO or an external `credential_process`, does not
work inside the container. Resolve on the host and export the result:

```bash
eval "$(aws configure export-credentials --profile <your-profile> --format env)"
```

On HealthOmics the run role supplies credentials to each task.

## Quick start

Build the runtime image from the repository root:

```bash
docker build -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile .
```

Prepare a manifest. This step queries the Discover API, so run it on a host with internet access.
It needs only the Python standard library.

```bash
tools/cxg-manifest.py --limit 5 --sort size_asc --out manifest.jsonl
```

Stage the pinned gene references. The writer projects every dataset onto them, and a task cannot
fetch them. The tool needs the `census2-dataset-writer` package, which holds the pin, and
`aria2c`. The runtime image has both. `aria2c` downloads the files in parallel chunks, retries,
resumes a partial file, and checks each file against its pinned SHA-256 digest.

```bash
docker run --rm -v "$PWD:/w" -w /w census2-bulk-ingest:dev \
    python tools/stage-refdata.py refdata/
```

Ingest it:

```bash
eval "$(aws configure export-credentials --profile <your-profile> --format env)"

nextflow run . -profile docker \
    --manifest manifest.jsonl --refdata refdata/ --check \
    --outdir ./out
```

The stores go to `./out/<dataset_version_id>.zarr`, with one receipt JSON beside each store.

To write directly to S3, give `--outdir` an `s3://` prefix and read the manifest from S3 as well.
The task writes each store to its final location; there is no staging step.

```bash
aws s3 cp manifest.jsonl s3://my-bucket/census2-ingest/manifest.jsonl

nextflow run . -profile docker \
    --manifest s3://my-bucket/census2-ingest/manifest.jsonl --refdata refdata/ \
    --outdir s3://my-bucket/census2-ingest
```

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--manifest` | none, required | JSONL manifest from `tools/cxg-manifest.py`. Local path or `s3://` URI. |
| `--outdir` | `./out` | Local path or `s3://` prefix. |
| `--all_slots` | `false` | Carry `obsm`, `varm`, `obsp`, `varp`, `layers`, and the full `uns`, instead of the default Census slot set. |
| `--check` | `false` | Verify each output against its source: `obs`, and the counts at their gene reference positions. A failure goes into the receipt, and `REPORT` fails the run. |
| `--config` | none | Ingestion configuration file, staged into every task. The file format is not defined yet. |
| `--refdata` | none, required | Directory or `s3://` prefix that holds the gene references from `tools/stage-refdata.py`, staged into every task. |
| `--report_dir` | none | Where the run report publishes. The `healthomics` profile sets it. |

`tools/cxg-manifest.py` selects which datasets to ingest, by organism, size, and count. Run it
with `--help` for the options.

## Prototyping manifest

`manifests/census1-2025-11-08-homo_sapiens-100.jsonl` is the shared corpus for prototype work.
It holds 100 datasets from the `homo_sapiens` experiment of Census 1 LTS release `2025-11-08`,
sampled uniformly at random from the 1086 datasets there. It is 147 GiB of h5ad and 13.3M cells.

Each record names the dataset version that the release holds, at CxG schema 7.0.0, and not the
current Discover version. A store made from it therefore holds the same source data as Census 1.
You can then compare a query against the stores with the same query against the release,
filtered to these `dataset_id`s.

Census 1 does not hold every cell of every dataset, because it drops cells that fail its own
filters. For example, it drops cells from cell culture. In this sample, one dataset is partial:
Census 1 holds the 5594 organoid cells of `ef1a4b05-540c-4f73-8ceb-49dfa800645d` and drops its
10212 primary cell culture cells. A comparison must apply the same filters to the stores.

`tools/census1-sample.py` writes the manifest. The script reads the release with
`cellxgene-census`, so it runs under uv:

```bash
tools/census1-sample.py --census-version 2025-11-08 --seed 0 --limit 100 \
    --out manifests/census1-2025-11-08-homo_sapiens-100.jsonl
```

The script shuffles all 1086 datasets with the seed and takes the first `--limit`. Thus no dataset
occurs twice, and the first N records of the manifest are also a uniform sample. A larger sample
with the same seed starts with this one. `--limit 0` gives all 1086 datasets (1.7 TiB) in that
order.

`manifests/census1-2025-11-08-homo_sapiens-harness.jsonl` is a companion manifest. It holds the
three datasets of the human sample of the validation harness (`census2-sdk/tools/harness/`), at
the versions that the release holds. The harness chose them for their properties, so they are not
part of the uniform sample. Ingest both manifests to the same `--outdir`.

| Dataset | Property |
| --- | --- |
| `5b25a46f-f1cf-4caa-9e10-6ab81dfe9174` | A multi-valued disease |
| `6e00ccf7-0749-46ef-a999-dba785630d52` | Primary data, lung cells |
| `62315937-e268-4fa5-a032-8f7d776f3a3f` | The control: the SDK and Census 1 agree on every value |

`--dataset-id` names a dataset instead of a sample:

```bash
tools/census1-sample.py --census-version 2025-11-08 \
    --dataset-id 5b25a46f-f1cf-4caa-9e10-6ab81dfe9174 \
    --dataset-id 6e00ccf7-0749-46ef-a999-dba785630d52 \
    --dataset-id 62315937-e268-4fa5-a032-8f7d776f3a3f \
    --out manifests/census1-2025-11-08-homo_sapiens-harness.jsonl
```

## Profiles

| Profile | Use |
|---|---|
| `docker` | Run tasks locally in the runtime image. Credentials come from the environment. |
| `healthomics` | Run on AWS HealthOmics. `--container` is required and must be a private ECR URI. |

## Output layout

```
<outdir>/<dataset_version_id>.zarr/    the store, if the dataset was admitted
<outdir>/<dataset_version_id>.json     the receipt, always
```

The store holds its summary record, with the dataset's identifiers and the gene reference version,
in the attributes of its root `zarr.json`. The `census2-dataset-writer` README describes the
record. The receipt is the run log of the pipeline, and the run report reads it.

A rejection is not an error: the task succeeds, no store is written, and the receipt records
`admitted: false` with a reason.

This flat layout is a prototype convenience: a rerun overwrites the store it wrote before. The
run report publishes to `--report_dir`, not `--outdir`; a local run leaves it in the task's
working directory unless you set `--report_dir`. On HealthOmics, the run report goes to the
`pubdir/` prefix of the run output location. Nextflow's `report.html`, `timeline.html`, and
`trace.txt` go to the `output/` prefix. `trace.txt` records the peak memory and the runtime of each
task.

## Running a component alone

Each step is a plain script and runs without Nextflow. `tools/cxg-manifest.py` and
`tools/seed-mock-s3.py` need only the Python standard library, so they run anywhere.
`cxg-ingest.py` needs the image's Python dependencies, so run it in the image, or install the
package: `pip install '../census2-dataset-writer[s3]'`.

```bash
modules/local/ingest/resources/usr/bin/cxg-ingest.py \
    --record "$(head -1 manifest.jsonl)" --outdir ./out --refdata refdata/ --check

modules/local/report/resources/usr/bin/ingestion-report.py ./out/*.json \
    --out ingestion_report.json
```

## AWS HealthOmics

Build the bundle, push an `amd64` image to private ECR, then create a workflow and start a run:

```bash
tools/build-healthomics-bundle.sh healthomics-bundle.zip

docker buildx build --platform linux/amd64 \
    -t <account>.dkr.ecr.<region>.amazonaws.com/<repo>:census2-bulk-ingest-0.1.0 \
    -f docker/Dockerfile --push ..

aws omics create-workflow --engine NEXTFLOW --main main.nf \
    --definition-zip fileb://healthomics-bundle.zip \
    --parameter-template file://tools/healthomics-parameters.json \
    --storage-type DYNAMIC --name census2-bulk-ingest
```

Pass the destination as `s3_//bucket/prefix`, not `s3://`: HealthOmics validates an `s3://`
parameter value as an existing input object, and an output prefix is not one. `--outdir` must be
remote on HealthOmics, and the pipeline refuses a local path. Stage `--refdata` to an S3 prefix
with `tools/stage-refdata.py s3://bucket/prefix/refdata/`, and pass that prefix as an `s3://`
input.

## Tests

Tests use [nf-test](https://www.nf-test.com) and run against the pinned manifest in
`tests/fixtures/manifest.jsonl`. Every test needs the gene references in `tests/refdata`, which
is not tracked. Stage them once, as CI does:

```bash
docker run --rm -v "$PWD:/w" -w /w census2-bulk-ingest:dev \
    python tools/stage-refdata.py tests/refdata
```

```bash
nf-test test --tag validation         # fast; no network and no credentials after staging
nf-test test --tag validation,ingest  # plus the pinned ingests; what CI runs
nf-test test                          # everything, including the test that writes to S3
```

Install nf-test with `curl -fsSL https://get.nf-test.com | bash`. It needs a JRE on `PATH`, which
Nextflow does not provide.

The ingest tests need either AWS credentials in the environment or a mock endpoint. To run them
with no credentials, serve the fixture datasets from a mock S3, under the bucket and key the
manifest already names:

```bash
pipx install 'moto[server]'
moto_server -H 0.0.0.0 -p 9000 &

export AWS_ACCESS_KEY_ID=test AWS_SECRET_ACCESS_KEY=test
AWS_ENDPOINT_URL=http://127.0.0.1:9000 tools/seed-mock-s3.py

# A task runs in a container, so it needs an address that reaches the host. On
# Linux use the docker bridge gateway instead of this name.
export AWS_ENDPOINT_URL=http://host.docker.internal:9000
nf-test test --tag validation,ingest
```

Regenerate the pinned manifest to change the fixtures:

```bash
tools/cxg-manifest.py --limit 2 --sort size_asc --out tests/fixtures/manifest.jsonl
```
