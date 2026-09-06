# chanzuckerberg/cellxgene-census2-prototype

In this repository we're prototyping a successor to [CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census), a scRNA data corpus and SDK.

The original Census reprocessed Discover scRNA AnnData datasets into a monolithic TileDB structure published on Amazon S3, with client libraries lightly wrapping TileDB's embedded query engine.

In Census 2 we're migrating from TileDB to AnnData's maturing Zarr materialization, and replacing the monolithic data structure with a Data Lake approach. The tentpoles are: **one AnnData-Zarr store per dataset, written independently, composed at query time through a catalog.**

- Each dataset the Census admits is one Zarr store in CELLxGENE Discover's asset storage, beside the dataset's existing artifacts, addressed by `dataset_version_id` and `artifact_revision`.
- Stores are written one dataset at a time; there is no corpus-wide data build.
- Each store is preprocessed with some (but not all) Census1-insipired harmonization steps, facilitating cross-dataset queries.
- The Census owns a catalog and releases. The catalog has one summary row per dataset, projected from Discover's curation API, so a query can skip a store without opening it. Each release is a manifest and a catalog snapshot that pin exact dataset versions.
- A query reads the catalog, opens only the stores that survive pruning, and combines them through stable cell and gene identity. One store is one access boundary, so a private dataset joins the same query as the public corpus once its reader holds a grant, without the data moving.

**Rough consensus and running code:** beneath the tentpoles above, our design and roadmap are subject to iterative refinement with experience and feedback.

## Documentation policy

For all project documentation (markdown files, docstrings, and code comments):

- MANDATORY DIRECTIVE: Write in the style of ASD-STE100 Simplified Technical English.
- Assume the reader has general background in scRNA bioinformatics, making domain terminology and acronyms available without verbose introduction.
- Communicate complex idea/designs using progressive disclosure. The first take may be simplified or otherwise imperfectly accurate, in order to paint a *readable* big picture. Then follow up with necessary details & nuance.
- Always mind overall readability and coherence when editing documentation; the smallest correct edit is not always the best one to make. Prefer to state the current truth as if you knew it from the start, instead of exhaustively documenting the winding road there; git naturally retains edit history.
- AVOID sentence structures and rhetorical flourishes intending to elicit human emotional reactions. This isn't a dramatic novel--it's technical documentation.

Note on humility: our project will be reviewed by the scverse core team under a funded collaboration, and eventually opened for public consumption. Our documentation must be honest, but also sensitive to human social dynamics, when we discuss the capabilities and/or performance of scverse dependencies, TileDB-SOMA, and other related projects.
