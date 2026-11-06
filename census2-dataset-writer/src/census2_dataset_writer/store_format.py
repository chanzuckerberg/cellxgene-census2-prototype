"""The format of a Census 2 store: the names and records that a reader relies on.

A store is an AnnData-Zarr group. This module names what Census 2 adds to that
format, so that the writer and a reader of a store use one definition:

- `DATASET_ID_COLUMN`, the `obs` column that names each cell's dataset.
- `OBS_COLUMNS`, the columns of a store's `obs`, in order, each with its
  `ColType`, and `DERIVED_OBS_COLUMNS`, the ones that the writer computes.
- `MEASURED_COLUMN`, the `var` column of the dataset's assayed genes.
- The summary record in the root group's attributes, under
  `attrs[SUMMARY_NAMESPACE][SUMMARY_KEY]`: `SummaryRecord`, `SUMMARY_VERSION`,
  and `validate_summary`, which raises `SummaryRecordError`.
- The compaction record of a compacted store, under
  `attrs[SUMMARY_NAMESPACE][COMPACTION_KEY]`: `CompactionRecord`,
  `COMPACTION_VERSION`, and `validate_compaction`. A compacted store is a second
  copy of several per-dataset stores in one store; see
  `census2_dataset_writer.compaction`.
- `obs_digest`, the identity of a store's cells: a digest of its
  `observation_joinid` values in row order.

The module imports only the standard library, so a reader could depend on it
without the writer's dependencies. The package re-exports these names. The SDK
repeats some of them today, because it does not depend on the writer.
"""

import hashlib
from collections.abc import Iterable
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


# Where the compaction record is: a key beside SUMMARY_KEY. A compacted store
# holds a compaction record and no summary record, so a reader that expects one
# dataset in a store refuses a compacted store.
COMPACTION_KEY = "compaction"

# The format of the compaction record. Raise it when a field changes meaning or
# goes away. A reader refuses a version that it does not know.
COMPACTION_VERSION = 1


class CompactedMember(TypedDict):
    """One dataset in a compacted store: its identity, and where its rows are.

    The rows of a member are `row_start` to `row_stop` (exclusive) of the
    compacted store, in the row order of its per-dataset store. `nnz` and
    `obs_digest` identify the per-dataset store that the rows came from, so a
    catalog builder can check the two copies against each other cheaply.
    """

    dataset_id: str
    dataset_version_id: str
    row_start: int
    row_stop: int
    # The stored values of the member's rows of `X`.
    nnz: int
    # `obs_digest` of the member's `observation_joinid` values.
    obs_digest: str
    # The `uns` of the per-dataset store. A compacted store has one `uns` for
    # all of its members, so the values of each member are kept here.
    uns: dict[str, str]


class CompactionRecord(TypedDict):
    """The facts about a compacted store that a catalog needs, held in the store."""

    version: int
    # The gene reference of every member, as in `SummaryRecord`.
    gene_reference_version: str
    # In row order: the first member starts at row 0, and each member starts
    # where the one before it stops.
    members: list[CompactedMember]


def validate_compaction(obj: object) -> CompactionRecord:
    """Return `obj` as a compaction record, or raise `SummaryRecordError` saying why not.

    A compaction record is read from a store, which is external input, so its type
    is checked here at run time. The record must have exactly the fields of
    `CompactionRecord`, and each member exactly the fields of `CompactedMember`.
    The members must cover the rows from 0 without a gap or an overlap, each with at
    least one row, and no dataset may occur twice.
    """
    if not isinstance(obj, dict) or set(obj) != set(CompactionRecord.__annotations__):
        raise SummaryRecordError(
            f"a compaction record is an object with the fields "
            f"{sorted(CompactionRecord.__annotations__)}"
        )
    version = obj["version"]
    if type(version) is not int or version != COMPACTION_VERSION:
        raise SummaryRecordError(
            f"compaction record version {version!r}; "
            f"this package reads {COMPACTION_VERSION}"
        )
    if not isinstance(obj["gene_reference_version"], str):
        raise SummaryRecordError(
            "compaction record gene_reference_version is not a string"
        )
    members = obj["members"]
    if not isinstance(members, list) or not members:
        raise SummaryRecordError("a compaction record lists one or more members")
    fields = set(CompactedMember.__annotations__)
    stop = 0
    for member in members:
        if not isinstance(member, dict) or set(member) != fields:
            raise SummaryRecordError(
                f"a compacted member is an object with the fields {sorted(fields)}"
            )
        for field in ("dataset_id", "dataset_version_id", "obs_digest"):
            if not isinstance(member[field], str) or not member[field]:
                raise SummaryRecordError(f"compacted member {field} is not a string")
        for field in ("row_start", "row_stop", "nnz"):
            if type(member[field]) is not int:
                raise SummaryRecordError(f"compacted member {field} is not an integer")
        if member["row_start"] != stop or member["row_stop"] <= member["row_start"]:
            raise SummaryRecordError(
                f"compacted member {member['dataset_version_id']} has rows "
                f"{member['row_start']}:{member['row_stop']}; the next row is {stop}"
            )
        stop = member["row_stop"]
        uns = member["uns"]
        if not isinstance(uns, dict) or not all(
            isinstance(k, str) and isinstance(v, str) for k, v in uns.items()
        ):
            raise SummaryRecordError("compacted member uns is not an object of strings")
    for field in ("dataset_id", "dataset_version_id"):
        values = [m[field] for m in members]
        if len(set(values)) != len(values):
            raise SummaryRecordError(f"a {field} occurs twice in a compaction record")
    return cast("CompactionRecord", obj)


def obs_digest(joinids: Iterable[str]) -> str:
    """The identity of a sequence of cells: SHA-256 of the `observation_joinid` values.

    Each value is UTF-8, followed by a newline. The digest depends on the order, so
    two stores with the same digest hold the same cells in the same order.
    """
    digest = hashlib.sha256()
    for value in joinids:
        digest.update(value.encode())
        digest.update(b"\n")
    return digest.hexdigest()


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
    # NumPy int64.
    INT64 = "int64"
    # NumPy float64.
    FLOAT64 = "float64"


# The columns of a store's `obs`, in Census 1's order: `dataset_id`, the columns
# that Census 1 copies from every CxG h5ad, the `tissue_general` pair, and five
# statistics of the counts. A store holds no other `obs` column yet. Census 1 also
# writes `soma_joinid`, a corpus-wide number that a store for one dataset cannot
# assign. No value is null.
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
    ("tissue_general", ColType.CATEGORY),
    ("tissue_general_ontology_term_id", ColType.CATEGORY),
    ("raw_sum", ColType.FLOAT64),
    ("nnz", ColType.INT64),
    ("raw_mean_nnz", ColType.FLOAT64),
    ("raw_variance_nnz", ColType.FLOAT64),
    ("n_measured_vars", ColType.INT64),
)

# The writer computes these columns. Each other column copies the source.
DERIVED_OBS_COLUMNS = frozenset(
    {DATASET_ID_COLUMN, "tissue_general", "tissue_general_ontology_term_id"}
    | {"raw_sum", "nnz", "raw_mean_nnz", "raw_variance_nnz", "n_measured_vars"}
)

# The dataset's assayed genes. `X` cannot show them, because an unmeasured gene is 0.
# It is the only part of `var` that differs between stores of one organism.
MEASURED_COLUMN = "measured"
