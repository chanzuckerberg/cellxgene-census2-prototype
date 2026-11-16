"""The format of a Census 2 store: the records that a reader relies on.

A store is an AnnData-Zarr group. The schema of its `obs`, and the `var` column of
the assayed genes, are in the SDK (`census2_sdk.schema`), which the writer and every
reader import. This module names the rest of what Census 2 adds to the format:

- `DERIVED_OBS_COLUMNS`, the `obs` columns that the writer computes.
- The summary record in the root group's attributes, under
  `attrs[SUMMARY_NAMESPACE][SUMMARY_KEY]`: `SummaryRecord`, `SUMMARY_VERSION`,
  and `validate_summary`, which raises `SummaryRecordError`.

The package re-exports these names.
"""

from typing import TypedDict, cast

from census2_sdk.schema import DATASET_ID_COLUMN

# Where the summary record is: a key of the root group's attributes, and a key
# under it. The namespace keeps the record clear of anndata's own attributes.
SUMMARY_NAMESPACE = "cellxgene-census2"
SUMMARY_KEY = "summary"

# The format of the summary record. Raise it when a field changes meaning or
# goes away. A reader refuses a version that it does not know.
SUMMARY_VERSION = 1


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


# The writer computes these columns. Each other column copies the source.
DERIVED_OBS_COLUMNS = frozenset(
    {DATASET_ID_COLUMN, "tissue_general", "tissue_general_ontology_term_id"}
    | {"raw_sum", "nnz", "raw_mean_nnz", "raw_variance_nnz", "n_measured_vars"}
)
