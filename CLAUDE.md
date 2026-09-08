# chanzuckerberg/cellxgene-census2-prototype

This repository prototypes a successor to [CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census), a scRNA data corpus and SDK.

Census 1 reprocesses Discover's scRNA AnnData datasets into one monolithic TileDB structure on Amazon S3. Its client libraries are thin wrappers around TileDB's embedded query engine.

Census 2 replaces the monolithic data structure with an S3 data lake. It also moves from TileDB to AnnData's maturing Zarr materialization. **The tentpoles are: one AnnData-Zarr store per dataset, written independently, composed at query time through a catalog.**

- **One store per dataset.** Each dataset that the Census admits becomes one Zarr store. The store goes into CELLxGENE Discover's asset storage, beside the dataset's other artifacts. Its address is the `dataset_version_id` and the `artifact_revision`.
- **Independent writes.** We write the stores one dataset at a time. There is no corpus-wide data build. A store depends only on its dataset and the schema version.
- **Harmonization without aggregation.** Cross-dataset queries need harmonized data. We apply some of Census 1's harmonization steps during Zarr ingestion (so a store is *not* merely the Discover h5ad saved in a different format). Ingestion reads only one dataset at a time, so steps that needed corpus-wide aggregates are redesigned or eliminated.
- **Catalog and releases.** The Census owns both. The catalog holds one summary row per dataset, projected from Discover's curation API. A query uses these summaries to skip irrelevant stores without opening them. A release is a manifest and a catalog snapshot, which together pin the exact dataset versions. Because ingestion does not use corpus-wide aggregates, many releases can pin the same store.
- **Query-time composition.** A query reads the catalog and opens only the stores that survive pruning. It combines those stores through stable cell and gene identity. One store is one access boundary. A private dataset therefore joins the same query as the public corpus as soon as its reader holds a grant, and the data does not move.

## Rough consensus and running code

Beneath those tentpoles, our design and roadmap are subject to iterative refinement with experience and feedback.

When we change an interface, do not add backward-compatibility shims (unless the operator so directs). Instead, update callers and other dependents to the revised interface.

## Python style

Lint and format the Python code with `pipx run ruff` before you commit a change to it. Run it in
this order, and let ruff fix as much as it can before you touch anything yourself:

```bash
pipx run ruff check --fix .     # ruff fixes what it can
pipx run ruff format .
pipx run ruff check .           # then work through what is left, by hand
```

`pipx run` takes the current ruff release, and we take its default rules and its default width.
There is no ruff configuration file, and nothing pins a version. This suits a small codebase at the
start of a project: the defaults cost us no decisions, and a release that moves them costs little
to absorb. CI runs these same unpinned commands, so a runner and a workstation cannot disagree. The
cost is that a ruff release can fail a commit that changed nothing, and the fix is to run the
commands again and commit the result.

Every Python file carries a `.py` extension, the Nextflow module binaries included. Ruff finds
files by extension, so a script without one is skipped without a word. That is also what keeps a
configuration file unnecessary today.

`ruff format` also reformats the Python blocks inside the markdown files, so expect `README.md` to
appear in the diff.

**A standalone package also type-checks under `mypy --strict`.** `census2-store` is one; the
Nextflow scripts are not. The line is what a caller can depend on. A package publishes an API, so
its types are part of the contract, and it ships a `py.typed` marker to say so. The scripts are
glue around `argparse` and `subprocess`, read top to bottom, and imported by nothing, so
annotating them would cost more than it catches.

```bash
pip install './census2-store[dev]' && mypy census2-store/src
```

The settings live in the package's own `pyproject.toml` rather than in a repository-wide file, and
mypy runs from an environment where the package's dependencies are installed. It cannot run
isolated the way `pipx run ruff` does, because it has to resolve the types of what the package
imports.

## Documentation policy

This policy applies to all project documentation: markdown files, docstrings, and code comments.

- MANDATORY DIRECTIVE: write in the style of ASD-STE100 Simplified Technical English.
- Assume the reader has a general background in scRNA bioinformatics. Use the domain terms and acronyms without verbose introduction.
- Disclose complex ideas and designs progressively. The first pass can be simplified, or imperfectly accurate, if it gives a *readable* big picture. Add the necessary details and nuance after it.
- Keep the whole sentence/section/document readable and coherent when you edit it. The smallest correct edit is not always the best edit to make.
- State the current truth as if you knew it from the start. Avoid verbosely documenting the journey there; git retains that history.
- AVOID sentence structures and rhetorical flourishes that intend to elicit human emotional reactions. This isn't a dramatic novel—it's technical documentation.

Humility note: the scverse core team will review this project under a funded collaboration, and the project will eventually open for public consumption. What we write about the capabilities and the performance of scverse dependencies, TileDB-SOMA, and other related projects must be honest. It must also be sensitive to human social dynamics.
