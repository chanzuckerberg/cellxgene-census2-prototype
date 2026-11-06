# cellxgene-census2-prototype

This repository holds the prototype of Census 2, the successor to
[CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census). The prototype is
in development and is not stable.

Census 2 publishes one AnnData-Zarr store for each dataset that it admits. The census dataset writer writes each
store independently, next to the other artifacts of the dataset. A catalog holds one summary row
for each dataset. At query time, the SDK uses the catalog to exclude the stores that cannot match.
It opens the stores that remain and returns one AnnData. A release pins the exact dataset versions
and a catalog snapshot. The [design record](https://github.com/chanzuckerberg/cellxgene-census2-prototype/wiki/Design-record) on the project wiki states the
design in full.

## Components

| Directory | What it is |
|---|---|
| [`census2-dataset-writer`](census2-dataset-writer/README.md) | A Python package that converts one Discover dataset into one store. It decides admission, the slot set, and harmonization. |
| [`bulk-ingest-nf`](bulk-ingest-nf/README.md) | A Nextflow pipeline that runs the writer over many datasets, locally under Docker or on AWS HealthOmics. |
| [`census2-sdk`](census2-sdk/README.md) | A Python package that queries stores through a catalog and returns one AnnData. |
| [`compact-nf`](compact-nf/README.md) | A Nextflow pipeline that writes compacted stores: a second copy of several small per-dataset stores in one store. A spike; see [`spikes/store-compaction.md`](spikes/store-compaction.md). |

Each component README holds the setup and use instructions for that component.

## Design documents

The design documents are on the [project wiki](https://github.com/chanzuckerberg/cellxgene-census2-prototype/wiki).

| Document | Wiki page |
|---|---|
| The design record, version 3 | [Design record](https://github.com/chanzuckerberg/cellxgene-census2-prototype/wiki/Design-record) |
| The harmonization disposition table | [Disposition table](https://github.com/chanzuckerberg/cellxgene-census2-prototype/wiki/Disposition-table) |

A README or a docstring cites a section of the design record in a fixed form: design 2.2 names
section 2.2.

## Develop

[`AGENTS.md`](AGENTS.md) holds the lint, format, type-check, and documentation policy. CI runs
`ruff` on each Python component and the tests of each package. The documentation policy applies to
every markdown file, docstring, and code comment in this repository.

## Code of conduct

This project follows the Contributor Covenant
[code of conduct](https://github.com/chanzuckerberg/.github/blob/main/CODE_OF_CONDUCT.md). All
participants must obey it. Send reports of unacceptable behavior to
[opensource@chanzuckerberg.com](mailto:opensource@chanzuckerberg.com).

## Security

Do not open a public issue for a security issue. Send an email to
[security@chanzuckerberg.com](mailto:security@chanzuckerberg.com).
[SECURITY.md](SECURITY.md) states this policy.

## License

The contents of this repository and of its wiki are available under the [MIT license](LICENSE).
