"""The format of a Census 2 store: the names and records that a reader relies on.

A store is an AnnData-Zarr group. This module names what Census 2 adds to that
format, so that the writer and a reader of a store use one definition:

- `DATASET_ID_COLUMN`, the `obs` column that names each cell's dataset.
- `OBS_COLUMNS`, the columns of a store's `obs`, in order, each with its
  `ColType`.
- The summary record in the root group's attributes, under
  `attrs[SUMMARY_NAMESPACE][SUMMARY_KEY]`: `SummaryRecord`, `SUMMARY_VERSION`,
  and `validate_summary`, which raises `SummaryRecordError`.

The module imports only the standard library, so a reader could depend on it
without the writer's dependencies. The package re-exports these names. The SDK
repeats some of them today, because it does not depend on the writer.
"""

from enum import Enum
from typing import TypedDict, cast

# Where the summary record is: a key of the root group's attributes, and a key
# under it. The namespace keeps the record clear of anndata's own attributes.
SUMMARY_NAMESPACE = "cellxgene-census2"
SUMMARY_KEY = "summary"

# The format of the summary record. Raise it when a field changes meaning or
# goes away. A reader refuses a version that it does not know.
SUMMARY_VERSION = 1

# The `obs` column that names each cell's dataset. `build` adds it after
# harmonization, so a cell keeps its dataset when stores are concatenated.
DATASET_ID_COLUMN = "dataset_id"


class SummaryRecord(TypedDict):
    """The facts about a store that a catalog needs, held in the store.

    The store holds everything else that a catalog reads, such as the organism
    and the curation schema version in `uns`, so the record does not copy it.
    """

    version: int
    dataset_id: str
    dataset_version_id: str
    # The label of the gene reference in the store's `var`, as in
    # `genes.REFERENCES`, for example "GENCODE v48".
    gene_reference_version: str


class SummaryRecordError(ValueError):
    """A summary record that is malformed: a missing, unknown, or wrongly typed field."""


def validate_summary(obj: object) -> SummaryRecord:
    """Return `obj` as a summary record, or raise `SummaryRecordError` saying why not.

    A summary record is read from a store, which is external input, so its type is
    checked here at run time and not only asserted. The record must have exactly
    the fields of `SummaryRecord`, `version` must be `SUMMARY_VERSION`, and the
    other fields must be strings that are not empty.
    """
    if not isinstance(obj, dict):
        raise SummaryRecordError(
            f"a summary record is a JSON object, not {type(obj).__name__}"
        )
    fields = set(SummaryRecord.__annotations__)
    if set(obj) != fields:
        missing, extra = sorted(fields - set(obj)), sorted(set(obj) - fields)
        raise SummaryRecordError(
            f"summary record fields: missing {missing}, unknown {extra}"
        )
    version = obj["version"]
    if type(version) is not int or version != SUMMARY_VERSION:
        raise SummaryRecordError(
            f"summary record version {version!r}; this package reads {SUMMARY_VERSION}"
        )
    for field in sorted(fields - {"version"}):
        if not isinstance(obj[field], str) or not obj[field]:
            raise SummaryRecordError(
                f"summary record {field} is {obj[field]!r}, not a string"
            )
    return cast("SummaryRecord", obj)


class ColType(Enum):
    """The type of a store column, which a reader can rely on in every store.

    It is the store format's contract, not a pandas dtype. Each comment says how
    pandas holds the type.
    """

    # A pandas categorical whose categories are strings: the observed values,
    # sorted, with no unused category.
    CATEGORY = "category"
    # NumPy bool.
    BOOL = "bool"
    # Plain strings, not categorical.
    STRING = "string"


# The columns of a store's `obs`, in order. They are the columns that Census 1
# copies from every CxG h5ad, in Census 1's order, after `dataset_id`. A store
# holds no other `obs` column yet. Census 1 also writes the `tissue_general`
# pair and five statistics of the counts, which a store does not hold yet, and
# `soma_joinid`, a corpus-wide number that a store for one dataset cannot
# assign. No value is null: the CxG schema requires each of these.
OBS_COLUMNS: tuple[tuple[str, ColType], ...] = (
    (DATASET_ID_COLUMN, ColType.CATEGORY),
    ("assay", ColType.CATEGORY),
    ("assay_ontology_term_id", ColType.CATEGORY),
    ("cell_type", ColType.CATEGORY),
    ("cell_type_ontology_term_id", ColType.CATEGORY),
    ("development_stage", ColType.CATEGORY),
    ("development_stage_ontology_term_id", ColType.CATEGORY),
    ("disease", ColType.CATEGORY),
    ("disease_ontology_term_id", ColType.CATEGORY),
    ("donor_id", ColType.CATEGORY),
    ("is_primary_data", ColType.BOOL),
    ("observation_joinid", ColType.STRING),
    ("self_reported_ethnicity", ColType.CATEGORY),
    ("self_reported_ethnicity_ontology_term_id", ColType.CATEGORY),
    ("sex", ColType.CATEGORY),
    ("sex_ontology_term_id", ColType.CATEGORY),
    ("suspension_type", ColType.CATEGORY),
    ("tissue", ColType.CATEGORY),
    ("tissue_ontology_term_id", ColType.CATEGORY),
    ("tissue_type", ColType.CATEGORY),
)
