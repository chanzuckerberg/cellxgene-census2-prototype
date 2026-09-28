"""The gene reference that this SDK pins for each organism.

The writer projects every store onto the reference of its organism: the store's
`var` is the reference, in reference order, and `uns["gene_reference_version"]`
names it. Stores of one organism therefore concatenate by position.
The SDK holds only the version names. It refuses a store whose catalog row records
another version, and it does not remap genes during a query.
"""

from typing import TypedDict


class Organism(TypedDict):
    taxon: str
    version: str


# The versions of the single-cell-curation tag v7.0.0, which the writer pins.
ORGANISMS: dict[str, Organism] = {
    "homo_sapiens": {"taxon": "NCBITaxon:9606", "version": "GENCODE v48"},
    "mus_musculus": {"taxon": "NCBITaxon:10090", "version": "GENCODE M37"},
}
TAXON_TO_ORGANISM = {spec["taxon"]: name for name, spec in ORGANISMS.items()}
