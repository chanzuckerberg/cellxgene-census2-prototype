# census2-bulk-ingest

A Nextflow pipeline that converts CELLxGENE Discover datasets into Census 2 stores in bulk: one
AnnData-Zarr store per dataset, written independently and composed at query time through a
catalog. It runs locally under Docker, and on AWS HealthOmics.

A run starts from a manifest, one JSON record per dataset, which `tools/cxg-manifest.py` prepares
from the Discover API. The pipeline then runs one task per dataset: it fetches the source h5ad from
Discover's asset bucket, builds the store, and writes the store and a receipt to `--outdir`. A
final task reads every receipt and reports on the run as a whole.

The source is always on S3, and it refuses anonymous reads, so every run needs AWS credentials
that can read it. The destination is your choice: `--outdir` takes a local path or an `s3://`
prefix, and only a remote one also needs write access. HealthOmics requires a remote destination,
because its workers do not share the filesystem you launch from. Nothing in the pipeline reaches
the public internet, which is what lets it run there at all.

What a store contains is not decided here. That belongs to `census2-store`, a package at the
repository root: it decides whether a dataset is admitted, what a store carries, and writes it out,
and the runtime image installs it as a dependency. This pipeline is the I/O around it — fetching,
publishing, receipts, and per-dataset resources. Read that package's docs for the data model, and
`CLAUDE.md` here for how this pipeline is put together.

Two of those decisions are still stubs, which shapes what a run means today. Every dataset is
admitted, because which datasets Census 2 admits is undecided; a rejection, when there is one, is
not an error, and the receipt records `admitted: false` with a reason while the task succeeds.
Harmonization returns its input, so a store today still matches its source, which is what makes
`--check` able to compare the two. Both will change.

## Requirements

- Nextflow 24.04 or later.
- Docker.
- Internet access for `tools/cxg-manifest.py` only. The pipeline itself needs none.
- AWS credentials that can read the source bucket.

**Getting credentials into the container.** Tasks read them from the environment, and nothing is
mounted. If your profile resolves credentials through a helper such as SSO or an external
`credential_process`, that helper is not present inside the container, so resolve on the host
first and export the result:

```bash
eval "$(aws configure export-credentials --profile <your-profile> --format env)"
```

On HealthOmics this is handled for you: the run role supplies credentials to each task.

## Quick start

Build the runtime image:

```bash
docker build -t census2-bulk-ingest:dev -f bulk-ingest-nf/docker/Dockerfile ..
# ^ from the repository root: the image installs the census2-store package
```

Running the pipeline takes two steps. First prepare a manifest, then ingest it.

**Step 1 — prepare a manifest.** This queries the Discover API, so run it somewhere with
internet access. It needs only the Python standard library.

```bash
tools/cxg-manifest.py --limit 5 --sort size_asc --out manifest.jsonl
```

**Step 2 — ingest.** The pipeline reads only the manifest and S3, so it needs no internet access.

```bash
eval "$(aws configure export-credentials --profile <your-profile> --format env)"

nextflow run . -profile docker \
    --manifest manifest.jsonl --check \
    --outdir ./out
```

Results land in `./out/v1/<dataset_version_id>.zarr`, with one receipt JSON beside each store.

A local destination is published by Nextflow from the task's declared output, so any local path
works. A remote destination is written straight to object storage by the task, and the directory
output is absent.

To write straight to S3, give `--outdir` a remote prefix and read the manifest from S3 as well.
Stores go to their final location; there is no staging step.

```bash
aws s3 cp manifest.jsonl s3://my-bucket/census2-ingest/manifest.jsonl

nextflow run . -profile docker \
    --manifest s3://my-bucket/census2-ingest/manifest.jsonl \
    --outdir s3://my-bucket/census2-ingest
```

Policy lives in the configuration rather than in separate flags. `all_slots` is a config key
today, with `--all_slots` as a convenience override. Admission and harmonization will need more
configuration and some reference data; both are stubbed. `--config FILE` and `--refdata DIR` are
staged into every task and recorded in the receipt, but no config file format is defined yet.
They are staged rather than fetched because a HealthOmics task has no internet access.

Because a store is written in place, a reader can see it partially written. Consolidated metadata
is written last, so opening with `zarr.open_consolidated` fails on an incomplete store rather than
returning part of one.

## Pipeline parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--manifest` | none, required | JSONL manifest from `tools/cxg-manifest.py`. Local path or `s3://` URI. |
| `--outdir` | `./out` | Local path or `s3://` prefix. |
| `--all_slots` | `false` | Override the config's `all_slots`: carry `obsm`, `varm`, `obsp`, `varp`, `layers`, and the full `uns`. |
| `--check` | `false` | Verify each output against its source. A failure lands in the receipt and `REPORT` fails the run. Meaningful only while harmonization is a stub. |
| `--config` | none | Ingestion configuration file, staged into every task. Holds policy such as `all_slots`; the file format is not defined yet. |
| `--refdata` | none | Directory of reference data, staged into every task. Nothing reads it yet. |

## Manifest options

`tools/cxg-manifest.py` selects which datasets to ingest. Run `--help` for the full list.

| Option | Default | Meaning |
|---|---|---|
| `--organism` | `NCBITaxon:9606` | Ontology term ID. A dataset is included if any of its organisms matches. |
| `--limit` | `0` | Take the first N after sorting. `0` means all. |
| `--sort` | `size_asc` | `size_asc`, `size_desc`, or `random`. |
| `--seed` | `0` | Seed for `--sort random`. |
| `--schema-version` | none | Omit for the current corpus. See "Rakes: upstream data" in `CLAUDE.md`. |
| `--source-bucket` | `dataset-assets-public-prod` | Bucket holding the source `.h5ad` files. |

By default the pipeline writes the Census slot set: `obs`, `var`, `X`, `raw`, and two `uns` keys.
Use `--all_slots` when the consumer needs a complete CxG object.

## Profiles

| Profile | Use |
|---|---|
| `docker` | Run tasks locally in the runtime image. Credentials come from the environment. |
| `healthomics` | Run on AWS HealthOmics. `--container` is required and must be a private ECR URI. |

## Running a component alone

Each step is a plain script and runs without Nextflow. Nextflow puts the module binaries on
`PATH` inside a task; from a shell, call them by path. `cxg-ingest.py` needs the image's Python
dependencies, so run it in the image, or install the package: `pip install '../census2-store[s3]'`.

```bash
tools/cxg-manifest.py --limit 5 --sort size_asc --out manifest.jsonl

modules/local/ingest/resources/usr/bin/cxg-ingest.py \
    --record "$(head -1 manifest.jsonl)" --outdir ./out --check

modules/local/report/resources/usr/bin/ingestion-report.py ./out/v1/*.json \
    --out ingestion_report.json
```

## Using the store builder directly

The data model is an installable package, so a caller that wants the transformation without the
pipeline copies nothing out of `modules/`:

```bash
pip install './census2-store[s3]'      # the extra only if the destination is remote
```

`census2-store/README.md` has the API and the standalone example, and
`census2-store/CLAUDE.md` has the anndata and zarr behaviour behind it. One habit carries over
from there: `obstore` reads credentials from the environment and does not resolve AWS profiles, so
export them first, exactly as the pipeline requires.

## What runs where

`tools/cxg-manifest.py` and `tools/seed-mock-s3.py` are not part of the pipeline and are not
module binaries. They use only the Python standard library, so they run anywhere. `cxg-ingest.py`
and `ingestion-report.py` ship as Nextflow module binaries under
`modules/local/*/resources/usr/bin/`, and the data model arrives in the image as an installed
package rather than as a staged file.

## Output layout

```
<outdir>/<recipe>/<dataset_version_id>.zarr/    the store, if the dataset was admitted
<outdir>/<recipe>/<dataset_version_id>.json     the receipt, always
```

`<recipe>` is `census2_store.RECIPE`, `v1` today. It covers the anndata version, the Zarr format,
the X encoding, the slot policy, and harmonization. Raise it there when any of those change. Outputs are keyed on
`dataset_version_id`, which names an immutable snapshot, so a rerun can skip work by listing
the prefix and subtracting it from the manifest.

## AWS HealthOmics

Build the bundle, push an `amd64` image to private ECR, then create a workflow and start a run:

```bash
tools/build-healthomics-bundle.sh healthomics-bundle.zip

docker buildx build --platform linux/amd64 \
    -t <account>.dkr.ecr.<region>.amazonaws.com/<repo>:census2-bulk-ingest-0.1.0 \
    -f bulk-ingest-nf/docker/Dockerfile --push ..

aws omics create-workflow --engine NEXTFLOW --main main.nf \
    --definition-zip fileb://healthomics-bundle.zip \
    --parameter-template file://tools/healthomics-parameters.json \
    --storage-type DYNAMIC --name census2-bulk-ingest
```

Pass the destination as `s3_//bucket/prefix`, not `s3://`. HealthOmics validates an `s3://`
parameter value as an existing input object, and an output prefix is not one; the pipeline
normalizes the escaped form. `--outdir` must be remote on HealthOmics, and the pipeline refuses a
local path rather than discarding the output.

## Tests

Tests use [nf-test](https://www.nf-test.com). They run against a pinned manifest in
`tests/fixtures/manifest.jsonl`, so a run is reproducible rather than depending on whichever
datasets are smallest today.

```bash
nf-test test --tag validation       # parameter validation: fast, no network, no credentials
nf-test test --tag ingest           # ingests the pinned datasets from object storage
nf-test test --tag integration      # those ingests, plus the test that writes to S3
nf-test test --tag validation,ingest  # the union, which is what CI runs
nf-test test                        # everything
```

nf-test can only select tags, not exclude them, so a bare `nf-test test` includes every test that
reaches object storage.

Install nf-test with `curl -fsSL https://get.nf-test.com | bash`. It needs a JRE on `PATH`, which
Nextflow does not provide; set `JAVA_HOME` if `java` is missing.

The ingest tests need either AWS credentials in the environment or a mock endpoint, and they take
about a minute.

**Running the ingest tests with no credentials.** Serve the fixture datasets from a mock S3, under
the bucket and key the manifest already names. The bytes come from the public CDN.

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

The variable reaches the task because `conf/docker.config` passes it through, and the AWS CLI
honours it. `.github/workflows/bulk-ingest-nf.yml` is this same recipe with the Linux address.

Regenerate the pinned manifest when you want different fixtures:

```bash
tools/cxg-manifest.py --limit 2 --sort size_asc --out tests/fixtures/manifest.jsonl
```
