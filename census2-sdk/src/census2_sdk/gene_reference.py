"""The gene reference that this SDK pins for each organism.

The writer projects every store onto the reference of its organism: the store's
`var` is the reference, in reference order, and the store's summary record names
its version. Stores of one organism therefore concatenate by position (design 1.5).
The SDK holds only the version names. It refuses a store whose catalog row records
another version, and it does not remap genes during a query (design 4.3).
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

# The store's `var` column of the dataset's assayed genes. The writer names it in
# `census2_dataset_writer.store_format.MEASURED_COLUMN`.
MEASURED_COLUMN = "measured"
