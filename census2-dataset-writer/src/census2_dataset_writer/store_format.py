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
- The `obs` fragment next to the store: `fragment_location`, the row columns
  that the fragment adds to `obs` (`FRAGMENT_ROW_COLUMNS`), and the fragment
  record under `attrs[SUMMARY_NAMESPACE][FRAGMENT_KEY]`: `FragmentRecord`,
  `FRAGMENT_VERSION`, and `validate_fragment`.

The module imports only the standard library, so a reader could depend on it
without the writer's dependencies. The package re-exports these names. The SDK
repeats some of them today, because it does not depend on the writer.

The `obs` fragment
------------------
The fragment is one Parquet file next to the store, with one row for each cell.
It holds the store's `obs` columns, in the store's order, and three more int64
columns: the store row ordinal, and the start and end of the cell's values in
`X` (`X/indptr[i]` and `X/indptr[i + 1]`). A query engine can then select
cells from the fragments of many stores, and the SDK reads only the `X` ranges
of the selected cells. The fragment is not part of the Zarr hierarchy, so the
store stays a plain AnnData-Zarr store.

The fragment record in the store names the digest, the size, and the row count
of the fragment. The record is in the immutable store, so a builder can prove
that a fragment is the one that the writer wrote with the store.
"""

import re
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

# Where the fragment record is: a key under SUMMARY_NAMESPACE, beside the summary
# record.
FRAGMENT_KEY = "obs_fragment"

# The format of the fragment and its record. Raise it when a column or a field
# changes meaning or goes away. A reader refuses a version that it does not know.
FRAGMENT_VERSION = 1

# The fragment of `<name>.zarr` is `<name>.obs.parquet`, in the same directory.
FRAGMENT_SUFFIX = ".obs.parquet"

# The columns that the fragment adds after the `obs` columns, each int64. The
# store row ordinal is the cell's row in the store, from 0. `x_start` and `x_end`
# give the half-open range of the cell's stored values in `X/data` and
# `X/indices`.
STORE_ROW_COLUMN = "store_row"
X_START_COLUMN = "x_start"
X_END_COLUMN = "x_end"
FRAGMENT_ROW_COLUMNS = (STORE_ROW_COLUMN, X_START_COLUMN, X_END_COLUMN)


class FragmentRecord(TypedDict):
    """The fragment that the writer wrote with a store, held in the store."""

    version: int
    # The SHA-256 digest of the fragment file, as 64 lowercase hex digits.
    sha256: str
    # The size of the fragment file.
    bytes: int
    # The number of rows, which is the store's number of cells.
    rows: int


def fragment_location(store: str) -> str:
    """The location of the fragment of the store at `store`: a path or a URL.

    The store's name must end with `.zarr`. The fragment is in the same directory,
    with `FRAGMENT_SUFFIX` in place of `.zarr`.
    """
    base = store.rstrip("/")
    if not base.endswith(".zarr"):
        raise ValueError(f"a store name ends with .zarr: {store!r}")
    return base.removesuffix(".zarr") + FRAGMENT_SUFFIX


def validate_fragment(obj: object) -> FragmentRecord:
    """Return `obj` as a fragment record, or raise `SummaryRecordError` saying why not.

    The record must have exactly the fields of `FragmentRecord`. `version` must be
    `FRAGMENT_VERSION`, `sha256` a hex digest, and `bytes` and `rows` integers that
    are not negative.
    """
    if not isinstance(obj, dict):
        raise SummaryRecordError(
            f"a fragment record is a JSON object, not {type(obj).__name__}"
        )
    fields = set(FragmentRecord.__annotations__)
    if set(obj) != fields:
        missing, extra = sorted(fields - set(obj)), sorted(set(obj) - fields)
        raise SummaryRecordError(
            f"fragment record fields: missing {missing}, unknown {extra}"
        )
    if type(obj["version"]) is not int or obj["version"] != FRAGMENT_VERSION:
        raise SummaryRecordError(
            f"fragment record version {obj['version']!r}; "
            f"this package reads {FRAGMENT_VERSION}"
        )
    if not isinstance(obj["sha256"], str) or not re.fullmatch(
        "[0-9a-f]{64}", obj["sha256"]
    ):
        raise SummaryRecordError(f"fragment record sha256 is {obj['sha256']!r}")
    for field in ("bytes", "rows"):
        if type(obj[field]) is not int or obj[field] < 0:
            raise SummaryRecordError(f"fragment record {field} is {obj[field]!r}")
    return cast("FragmentRecord", obj)
