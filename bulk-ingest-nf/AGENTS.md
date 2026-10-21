# census2-bulk-ingest

Bulk ingestion of CELLxGENE Discover datasets into Census 2 stores, as a Nextflow pipeline.

`README.md` covers what it does and how to run it. The repository root `AGENTS.md` gives the
Census 2 tentpoles and the documentation policy. This file records the sharp edges: what cost us
effort to find, and what is easy to catch again.

## Python style

See the repository root `AGENTS.md`. Every script here is glue, so ruff applies and `mypy --strict`
does not.

## Seams

The `census2-dataset-writer` package decides whether a dataset is admitted and what a store
contains; its README and module docstring cover the three decisions and their stubs, and its
`AGENTS.md` has the anndata and zarr sharp edges. Read those before changing how a store is
written.

`cxg-ingest.py` handles what running that package as a pipeline needs: fetching the source,
publishing, and the receipt. It gives the writer the dataset's identifiers from the manifest
record, and the writer puts the facts about a store into the store, in its summary record. The
receipt is the pipeline's run log. Only `ingestion-report.py` reads it, and the catalog reads the
store. For an admitted dataset, the receipt holds a copy of the summary record, so that the report
opens no store. The receipt's shape is `census2_dataset_writer.Receipt`, and its store-level facts
come from `census2_dataset_writer.provenance`, so the script that writes a receipt and the script
that reads one agree on it. Policy that a caller might choose lives in the configuration rather than
in function arguments, so there is one place to look.

These arrangements suit the harmonization steps that exist: the `obs` projection and the gene
projection. Expect to revisit them when more steps are known, particularly if a step needs data
that the current read order does not fetch.

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

The superseded versions stay in the source bucket under their own `dataset_version_id`. The
curation API's `/curation/v1/datasets/<dataset_id>/versions` gives their metadata. A manifest can
therefore pin an older generation. The prototyping manifest pins the 7.0.0 versions that Census 1
LTS `2025-11-08` holds.

Each Census 1 release also keeps its own copy of each h5ad in a public bucket, keyed by
`dataset_id`. For each of the 1086 human datasets in `2025-11-08`, that copy and the Discover
version have the same size. `tools/census1-sample.py` asserts this for each record it writes.

`cellxgene-census` needs tiledbsoma, and tiledbsoma 2.3.0 publishes no wheel for Python 3.14. The
script's uv header therefore bounds Python below 3.14. Without the bound, uv selects 3.14 and tries
a source build, which fails.

Some datasets name an `obs` column with a `/` in it. Zarr gives `/` its path meaning, so the
write fails on that column. The names seen are clinical units (`SCr_(mg/dl)`,
`eGFR_(ml/min/1.73m2)`), compound cell types (`T/NK_cycling`), and deconvolution factors
(`q05_spot_factorsStem-like/TA`); the last two are cell2location output on Visium datasets, one
column per cell type. 6 of 100 randomly sampled human datasets hold one, 1 or 2 columns out of 60
to 102.

`census2_dataset_writer.project_obs` removes those columns, and `harmonize` calls it. The
receipt records the names in `dropped_obs_columns`, and the run report totals them, because a store
that holds less than its dataset should not be invisible. The round-trip check compares `obs`
against the source minus exactly those names, so an unreported drop still fails.

The slashes measured were all in `obs`. The store's `var` is the pinned gene reference, not the
source `var`, so a `/` in a source `var` column does not reach the store.

`census2_dataset_writer.project_genes` puts the counts (`raw/X`, else `X`) in the store's `X`, on
the pinned gene reference. The store has no `raw`. The store's summary record holds
`gene_reference_version`. The receipt records the dropped non-gene features in
`dropped_var_features`. Thus the round-trip check does not compare the matrix
byte for byte. It reads the store's `X` at the reference position of each kept source feature, and
compares the values with the source. It also requires equal nonzero counts. It does this one block
of rows at a time. A gene that is not in the reference fails the task, because the writer does not
map genes across reference versions.

anndata 0.12 added a `disallow_forward_slash_in_h5ad` setting, on by default, so anndata also
refuses to *write* such a name to an h5ad as well. The corpus datasets that carry one predate it.
A test fixture that needs one must turn the setting off; zarr has no equivalent opt-out, because
there the restriction belongs to the format.

Peak memory depends on the encoding of the matrix that the store holds. The store's `X` holds the
counts: `raw.X` when the source has `raw`, else `X`. Census 1 reads the same matrix. Some datasets
also hold a dense matrix beside sparse counts, for example a dense normalized `X` beside a sparse
`raw.X`. The writer does not read that dense matrix. Of the 1086 human datasets in Census 1 LTS
`2025-11-08`, 9 hold a dense `X` beside a sparse `raw.X`, and none holds dense counts.

The DNA methylation datasets (snmC-seq) hold dense counts: `X` is dense, and there is no `raw` and
no `layers`. Admission rejects them on their assay: `census2_dataset_writer.admit` requires one cell
on an assay allowlist, and no methylation assay is on it. Of 2228 datasets, 53 have no allowed
assay, and 18 dense datasets are among the 53. The criterion is the assay and not the encoding, so
admission accepts a dataset that holds dense counts on an allowed assay.

The memory request is still `size * 2.5` for every dataset. On the prototyping manifest, the
largest peak was 4.0 GiB, on a 16.4 GiB h5ad, and no task used more than 40% of its reservation.
If you tune the request, start from the measured fit for sparse counts: peak is about
`0.90 GiB + 0.20 * filesize`.

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

That covers a failed check. It does not cover a task that raises: `errorStrategy = 'finish'` lets
the running tasks end but starts no new one, so `REPORT` never runs and the run publishes no
after-action report. A 100 dataset run lost its report to 10 such tasks. The report is most wanted
when tasks fail, so this needs a design that collects the receipts that do exist.

A task that raises part way through a write leaves a partial store. anndata writes consolidated
metadata last, so the absence of `consolidated_metadata` in the root `zarr.json` identifies one.
The same run left 8 partial stores, and one complete store whose task then failed its check.

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

An ingest test fails intermittently on macOS with `OSError: Cannot call rmtree on a symbolic
link`, raised from `zarr.storage.LocalStore.clear()` when anndata starts a write. Nothing is a
symbolic link. The message comes from `shutil._rmtree_safe_fd`, which raises it when
`os.path.samestat(orig_st, os.fstat(topfd))` is false: a guard against a directory being swapped
between the stat and the open. A Docker Desktop bind mount trips that guard on its own, so the
failure follows no particular test and does not reproduce. It has never appeared in CI, which
runs on Linux. Re-run before investigating.

## Sharp edges: dependency tooling

The dependency set is declared once, in `census2-dataset-writer/pyproject.toml`. pip cannot install
only a project's dependencies, so the image stages that file beside an empty module, installs it,
and a later layer replaces the module with the real source. Docker then reuses the large dependency
layer when only the package source changes. Nothing pins an exact version; the bound that matters
is `anndata<0.14`, and every receipt records what wrote a store.

uv and pixi were evaluated, and the image stays on pip. uv makes the Dockerfile simpler but adds
another build tool. pixi works: an image built from a `pixi.lock` ran the pipeline end to end,
locks the Linux and macOS platforms in one file, and removes the apt layer and the AWS CLI
install, because conda-forge carries every dependency here. But its image is larger, it manages
the local environment and CI as well as the image, and every platform must be declared. Reconsider
pixi when a lock file is appropriate, or when dependencies outside Python arrive. Two points are
untested: TLS to real S3 with no apt `ca-certificates`, and whether the environment relocates.

## Sharp edges: AWS HealthOmics and credentials

HealthOmics exports two directories of a Nextflow run. A task publishes files to
`/mnt/workflow/pubdir` with `publishDir`, and HealthOmics exports them to the `pubdir/` prefix of
the run output location. The run report goes there. The engine writes its own files to
`/mnt/workflow/output/`, and HealthOmics exports them to the `output/` prefix, on Nextflow 25.10 and
later. Nextflow's `report.html`, `timeline.html`, and `trace.txt` are engine files, so
`conf/healthomics.config` writes them to `/mnt/workflow/output/`. HealthOmics does not export an
engine file from any other directory, `pubdir` included, and it logs no warning. The trace records
the `peak_rss` and the `realtime` of each task. Nextflow runs `ps` in each task container to collect
them, so the image installs `procps`.

HealthOmics also records the task metrics itself, in the run manifest log: CloudWatch log group
`/aws/omics/WorkflowLog`, stream `manifest/run/<run id>/<uuid>`. It holds one JSON record for each
task, and the `metrics` of a record include `memoryMaximumGiB`, `memoryReservedGiB`, `cpusMaximum`,
and `runningSeconds`. The reservation is larger than the request, because HealthOmics rounds a
request up to the next of 4, 8, 16, 32, and 64 GiB: 5 GiB to 8, and 41 to 64.

HealthOmics validates any parameter whose value looks like an `s3://` URI as an existing input
object, and an output prefix is not one. Pass the destination as `s3_//bucket/prefix`, which the
pipeline normalizes.

HealthOmics selects the Nextflow version from `manifest.nextflowVersion`. It runs an exact
version that it supports, such as `'25.10.0'`. For a range, its documentation gives an order of
preference (23.10.0, 22.04.0, 24.10.8, 25.10.0, 26.04.0), but `'>=24.04'` gave 25.10.0 and not
24.10.8. The pipeline therefore pins an exact version, and CI sets `NXF_VER` to the same version,
so the tests run the engine that HealthOmics runs. `GetRun` gives the version of a run as
`engineVersion`. On 25.10, a `healthomics` profile test failed with `Unknown config attribute
params.container` until the test set `container`, and 26.x ran it without that, so both profile
tests set it.

`start-run` selects a Nextflow profile through `--engine-settings '{"profile": "<name>"}'`. The
HealthOmics documentation says that this works for every Nextflow version it supports. This
pipeline does not use it. Instead, `tools/build-healthomics-bundle.sh` appends
`includeConfig 'conf/healthomics.config'` when it builds the zip. The platform settings are thus
the defaults inside the bundle, and the configuration that a local run reads holds no
HealthOmics-only default. A `start-run` that selects the `healthomics` profile would make that
step unnecessary.

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
