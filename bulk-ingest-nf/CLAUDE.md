# census2-bulk-ingest

Bulk ingestion of CELLxGENE Discover datasets into Census 2 stores, as a Nextflow pipeline.

`README.md` covers what it does and how to run it. The repository root `CLAUDE.md` gives the
Census 2 tentpoles and the documentation policy. This file records the rakes: the things that cost
us effort to find, and that are easy to step on again.

## Python style

See the repository root `CLAUDE.md`. Every script here is glue, so ruff applies and `mypy --strict`
does not.


## Seams, and how provisional they are

The `census2-store` package decides what a store contains. It holds three decisions, two of them
placeholders waiting on the harmonization design:

- `admit` decides whether a dataset belongs in Census 2. It currently admits everything. A
  rejection is a normal outcome rather than an error, so it returns a reason and the task
  succeeds.
- Slot policy decides which AnnData slots a store carries. It defaults to the set the Census 1
  builder reads, which drops embeddings and most of `uns`. `all_slots` carries everything.
- `harmonize` returns its input.

`cxg-ingest.py` handles what running this as a pipeline needs: fetching the source, publishing, and
the receipt. Policy that a caller might choose lives in the configuration rather than in function
arguments, so that there is one place to look.

These arrangements suit a stub. Expect to revisit them once the harmonization steps are known,
particularly if a step needs data the current read order does not fetch.

## Rakes: upstream data

The source bucket is `s3://dataset-assets-public-prod/<dataset_version_id>.h5ad` in `us-west-2`.
The public CloudFront URL uses the same key. The bucket refuses anonymous reads, and it belongs to
another AWS account. Its policy grants at the account level rather than to named roles, so a
HealthOmics run role reads it. Evidence for the name is `backend/scripts/mirror_s3_data.sh` in the
`single-cell-data-portal` repository; elsewhere there it is a deployment secret.

Nothing in the chain publishes a content hash. The Discover API gives `filesize`, and the S3 ETag
is multipart for most objects. File size is the only integrity check, and the ingest step asserts
it.

Each CxG schema version is a complete, separately versioned copy of the corpus. A migration gives
every dataset a new `dataset_version_id`, so it invalidates every output keyed on one. Budget for a
full reconversion at each schema release.

## Rakes: anndata and zarr

These moved to `census2-store/CLAUDE.md`, beside the code they describe. Read them before changing
how a store is written.


## Rakes: Nextflow

`publishDir` with `mode: 'move'` removes the file from the work directory, so a downstream process
that consumes the same output then fails. The store moves and the receipt is copied, because
`REPORT` reads the receipt.

A `publishDir` `enabled:` value must be a static expression. Only `path:` accepts a closure.

`file()` is unavailable to the configuration parser, and a top-level `def` in a module `.nf` parses
as a process directive. Normalizing a parameter in `nextflow.config` also does not work, because
command-line and params-file overrides are applied after the configuration is evaluated.

A task that writes outside its work directory is a side effect Nextflow does not track. An earlier
version did that for every destination, which needed a bind mount whose path had to match the
task's path exactly, symlinks included. A mismatch discarded every store while each task reported
success. The store is now a declared optional directory output, published by Nextflow for a local
destination and written by the task for a remote one. `verify_written` covers the remote case,
where nothing tracks the write.

A constraint that only one execution profile has is better as a declared parameter than as
something the workflow infers. `params.remote_outdir_required` travels with the profile and can be
tested without the platform. Inferring it from the profile name was untestable.

Tests read a pinned manifest from `tests/fixtures`, so a failure means the code changed rather than
the corpus, and they assert on the receipt, which has to stay correct anyway.

Nextflow publishes nothing from a failed task. A task that has to report a failure downstream must
therefore exit 0 and record the failure in its output. `cxg-ingest.py` exits 0 on a failed round-trip
check for that reason, and `REPORT` fails the run. Exiting non-zero there hid the receipt that
named the failure and stopped `REPORT` from running at all.

An empty manifest gives a successful run that writes nothing, because no task runs and
`collect()` on an empty channel emits nothing, so `REPORT` is skipped as well. `main.nf` fails
early instead.

`nf-test --tag` selects tests and cannot exclude them, so a bare `nf-test test` runs the
integration tests too. Every test carries a tag, and the offline set is `--tag validation`. A test
that checks a validation rule but still reaches a task belongs with the integration set, because
starting a task needs Docker. Assert that the tasks were created: asserting only that an error
message is absent also passes when the run fails for some other reason.

A comma-separated `--tag validation,ingest` selects the union of the two, which is what CI runs.
`ingest` is a second tag on the ingest suite, so CI can take those four tests without the `s3_//`
escape test: that one writes to object storage, and against a mock it would pass without proving
anything. Note that `nf-test list` takes `--tags` rather than `--tag`, and does not filter the way
`test` does.

## Rakes: AWS HealthOmics and credentials

HealthOmics validates any parameter whose value looks like an `s3://` URI as an existing input
object. An output prefix is not one. Pass the destination as `s3_//bucket/prefix`, which the
pipeline normalizes. Creating a placeholder object at the output key also works, but leaves a stray
object behind.

There is no way to pass `-profile`. Neither `create-workflow` nor `start-run` accepts Nextflow
arguments, so the platform settings have to be the defaults inside the bundle.
`tools/build-healthomics-bundle.sh` appends that include when it builds the zip, which keeps the
repository free of a HealthOmics-only default that would surprise a local run.

HealthOmics has no boolean parameter type, so every parameter arrives as a string. A non-empty
string is true in Groovy, so `check=false` would switch the check on. Read a boolean parameter
through `.toString().toBoolean()`.

Tasks have no public internet access; `StartRun` reports `networkingMode: RESTRICTED`. Images must
be private ECR URIs in the same region, and ARM images are unsupported, so build `linux/amd64`.

A profile that resolves credentials with a helper, such as SSO or an external `credential_process`,
cannot work inside a container, because the helper is absent. Mounting the AWS configuration
directory makes it worse: the CLI then finds the configuration and tries to run a binary that is
not there. Resolve on the host and export with `aws configure export-credentials --format env`. On
HealthOmics and Batch the task role supplies credentials through the environment, so this matches
the deployment target.

Naming a credential variable in the Nextflow `env` scope exports an empty value when the launcher
does not have it, and Nextflow warns once for each empty variable. `-e NAME` in
`docker.runOptions` passes the value through and omits the variable instead.
