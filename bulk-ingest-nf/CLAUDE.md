# census2-bulk-ingest

Bulk ingestion of CELLxGENE Discover datasets into Census 2 stores, as a Nextflow pipeline.

`README.md` covers what it does and how to run it. The repository root `CLAUDE.md` gives the
Census 2 tentpoles and the documentation policy. This file records the sharp edges: what cost us
effort to find, and what is easy to catch again.

## Python style

See the repository root `CLAUDE.md`. Every script here is glue, so ruff applies and `mypy --strict`
does not.

## Seams

The `census2-writer` package decides whether a dataset is admitted and what a store contains; its
README and module docstring cover the three decisions and their stubs, and its `CLAUDE.md` has the
anndata and zarr sharp edges. Read those before changing how a store is written.

`cxg-ingest.py` handles what running that package as a pipeline needs: fetching the source,
publishing, and the receipt. The receipt's shape is `census2_writer.Receipt`, and its store-level
facts come from `census2_writer.provenance`, so the script that writes a receipt and the script
that reads one agree on it. Policy that a caller might choose lives in the configuration rather
than in function arguments, so there is one place to look.

These arrangements suit a stub. Expect to revisit them once the harmonization steps are known,
particularly if a step needs data the current read order does not fetch.

## Sharp edges: upstream data

The source bucket is `s3://dataset-assets-public-prod/<dataset_version_id>.h5ad` in `us-west-2`.
The public CloudFront URL uses the same key. The bucket refuses anonymous reads and belongs to
another AWS account. Its policy grants at the account level rather than to named roles, so a
HealthOmics run role reads it. Evidence for the name is `backend/scripts/mirror_s3_data.sh` in the
`single-cell-data-portal` repository; elsewhere in that repository it is a deployment secret.

Nothing in the chain publishes a content hash. The Discover API gives `filesize`, and the S3 ETag
is multipart for most objects. File size is the only integrity check, and the ingest step asserts
it.

Each CxG schema version is a complete, separately versioned copy of the corpus. A migration gives
every dataset a new `dataset_version_id`, so it invalidates every output keyed on one. Budget for
a full reconversion at each schema release.

## Sharp edges: Nextflow

`publishDir` with `mode: 'move'` removes the file from the work directory, so a downstream process
that consumes the same output then fails. The store moves and the receipt is copied, because
`REPORT` reads the receipt.

A `publishDir` `enabled:` value must be a static expression; only `path:` accepts a closure.
`file()` is unavailable to the configuration parser. A top-level `def` in a module `.nf` parses as
a process directive. Normalizing a parameter in `nextflow.config` does not work, because
command-line and params-file overrides apply after the configuration is evaluated.

A task that writes outside its work directory is a side effect Nextflow does not track, and a
mistake there discards stores while every task reports success. The store is therefore a declared
optional directory output: Nextflow publishes it for a local destination, and the task writes it
for a remote one, where `verify_written` covers the untracked write.

A constraint that only one execution profile has is better as a declared parameter than as
something the workflow infers. `params.remote_outdir_required` travels with the profile and can be
tested without the platform; inferring it from the profile name was untestable.

Nextflow publishes nothing from a failed task. A task that must report a failure downstream
therefore exits 0 and records the failure in its output: `cxg-ingest.py` exits 0 on a failed
round-trip check, and `REPORT` fails the run. A non-zero exit would hide the receipt that names
the failure and stop `REPORT` from running.

An empty manifest gives a successful run that writes nothing: no task runs, `collect()` on an
empty channel emits nothing, and `REPORT` is skipped. `main.nf` fails early instead.

Nextflow escapes a `path` value when it renders it into a task script, and renders every other
value verbatim. A `path` input therefore goes in bare, and quoting one turns the backslashes
Nextflow added into literal characters. A string is the opposite: quote it, and escape any quote
inside it, or the shell reads it as syntax; a dataset title with an apostrophe fails the task.
`--record` and `--outdir` are strings and are escaped; `--config` and `--refdata` are paths and
are bare. A `script:` block cannot call a helper closure, so the escape is written out at each
use.

nf-test's `--tag` selects tests and cannot exclude them, so a bare `nf-test test` runs the tests
that reach object storage too. Every test carries a tag, and the offline set is `--tag validation`.
A test that checks a validation rule but still starts a task belongs with the integration set,
because starting a task needs Docker. Assert that the tasks were created: asserting only that an
error message is absent also passes when the run fails for another reason. Tests read the pinned
manifest in `tests/fixtures`, so a failure means the code changed rather than the corpus, and they
assert on the receipt, which has to stay correct anyway.

CI runs `--tag validation,ingest`, the union of the two tags. `ingest` is a second tag on the
ingest suite so that CI takes those tests without the `s3_//` escape test, which writes to object
storage and proves nothing against a mock. `nf-test list` takes `--tags` rather than `--tag`, and
does not filter the way `test` does.

## Sharp edges: dependency tooling

The dependency set is declared once, in `census2-writer/pyproject.toml`. pip cannot install only a
project's dependencies, so the image stages that file beside an empty module, installs it, and a
later layer replaces the module with the real source. Docker then reuses the large dependency
layer when only the package source changes. Nothing pins an exact version; the bound that matters
is
`anndata<0.14`, and every receipt records what wrote a store.

uv and pixi were evaluated, and the image stays on pip. uv makes the Dockerfile simpler but adds
another build tool. pixi works: an image built from a `pixi.lock` ran the pipeline end to end,
locks the Linux and macOS platforms in one file, and removes the apt layer and the AWS CLI
install, because conda-forge carries every dependency here. But its image is larger, it manages
the local environment and CI as well as the image, and every platform must be declared. Reconsider
pixi when a lock file is appropriate, or when dependencies outside Python arrive. Two points are
untested: TLS to real S3 with no apt `ca-certificates`, and whether the environment relocates.

## Sharp edges: AWS HealthOmics and credentials

HealthOmics validates any parameter whose value looks like an `s3://` URI as an existing input
object, and an output prefix is not one. Pass the destination as `s3_//bucket/prefix`, which the
pipeline normalizes.

There is no way to pass `-profile`: neither `create-workflow` nor `start-run` accepts Nextflow
arguments, so the platform settings must be the defaults inside the bundle.
`tools/build-healthomics-bundle.sh` appends that include when it builds the zip, which keeps a
HealthOmics-only default out of the configuration that a local run reads.

HealthOmics has no boolean parameter type, so every parameter arrives as a string, and a non-empty
string is true in Groovy: `check=false` would switch the check on. Read a boolean parameter
through `.toString().toBoolean()`.

Tasks have no public internet access (`networkingMode: RESTRICTED`). Images must be private ECR
URIs in the same region, and ARM images are unsupported, so build `linux/amd64`.

A profile that resolves credentials with a helper, such as SSO or an external
`credential_process`, cannot work inside a container, because the helper is absent. Mounting the
AWS configuration directory also fails: the CLI then finds the configuration and tries to run a
binary that is not there. Resolve on the host and export with
`aws configure export-credentials --format env`. On HealthOmics and Batch the task role supplies
credentials through the environment, so this matches the deployment target.

Naming a credential variable in the Nextflow `env` scope exports an empty value when the launcher
does not have it, and Nextflow warns once for each empty variable. `-e NAME` in
`docker.runOptions` passes the value through and omits the variable instead.
