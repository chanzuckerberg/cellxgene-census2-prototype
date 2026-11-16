"""The schema of a store's `obs`: one typed attribute for each column, in store order.

The writer writes these columns, and the SDK reads and queries them, so both take
them from this module, and nothing else lists them. `Obs` gives each column its
class of the query builder (`census2_sdk.query`). Everything else here comes from its
annotations:

- `OBS_COLUMNS`, the names in store order, each with its `ColType`.
- `MULTIVALUED`, the columns whose cells join several values (`MultiValueColumn`).
  The store holds them as `ColType.CATEGORY`.
- `DATASET_ID_COLUMN`, the first column, which names each cell's dataset.

`MEASURED_COLUMN` is the `var` column of the dataset's assayed genes.
"""

from enum import Enum

from census2_sdk.query import (
    BoolColumn,
    CategoryColumn,
    Column,
    FloatColumn,
    IntColumn,
    MultiValueColumn,
    ObsNamespace,
    StringColumn,
)

__all__ = [
    "COLUMNS",
    "DATASET_ID_COLUMN",
    "MEASURED_COLUMN",
    "MULTIVALUED",
    "OBS_COLUMNS",
    "ColType",
    "Obs",
]


class ColType(Enum):
    """The type of a store column, which a reader can rely on in every store.

    It is the store format's contract, not a pandas dtype. Each comment says how
    pandas holds the type. No value is null.
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


class Obs(ObsNamespace):
    """The `obs` columns of one organism: `census.homo_sapiens.obs`.

    The columns are in Census 1's order: `dataset_id`, the columns that Census 1
    copies from every CxG h5ad, the `tissue_general` pair, and five statistics of
    the counts. A store holds no other `obs` column yet. Census 1 also writes
    `soma_joinid`, a corpus-wide number that a store for one dataset cannot assign.
    """

    dataset_id: CategoryColumn
    assay: CategoryColumn
    assay_ontology_term_id: CategoryColumn
    cell_type: CategoryColumn
    cell_type_ontology_term_id: CategoryColumn
    development_stage: CategoryColumn
    development_stage_ontology_term_id: CategoryColumn
    disease: MultiValueColumn
    disease_ontology_term_id: MultiValueColumn
    donor_id: CategoryColumn
    is_primary_data: BoolColumn
    observation_joinid: StringColumn
    self_reported_ethnicity: CategoryColumn
    self_reported_ethnicity_ontology_term_id: CategoryColumn
    sex: CategoryColumn
    sex_ontology_term_id: CategoryColumn
    suspension_type: CategoryColumn
    tissue: CategoryColumn
    tissue_ontology_term_id: CategoryColumn
    tissue_type: CategoryColumn
    tissue_general: CategoryColumn
    tissue_general_ontology_term_id: CategoryColumn
    raw_sum: FloatColumn
    nnz: IntColumn
    raw_mean_nnz: FloatColumn
    raw_variance_nnz: FloatColumn
    n_measured_vars: IntColumn

    def __init__(self, organism: str) -> None:
        super().__init__(organism)
        for name, cls in COLUMNS.items():
            setattr(self, name, cls(name, organism))


# The class of each column, in store order, from the annotations of `Obs`.
COLUMNS: dict[str, type[Column]] = dict(Obs.__annotations__)
_COL_TYPES: dict[type[Column], ColType] = {
    CategoryColumn: ColType.CATEGORY,
    MultiValueColumn: ColType.CATEGORY,
    BoolColumn: ColType.BOOL,
    StringColumn: ColType.STRING,
    IntColumn: ColType.INT64,
    FloatColumn: ColType.FLOAT64,
}
OBS_COLUMNS = tuple((name, _COL_TYPES[cls]) for name, cls in COLUMNS.items())
# Schema 2.4.0 made the `disease` columns multi-valued: a cell holds members joined
# by `query.MULTIVALUE_SEP`, and equality is membership.
MULTIVALUED = frozenset(n for n, cls in COLUMNS.items() if cls is MultiValueColumn)
# The writer adds it first, after harmonization, so a cell keeps its dataset when
# stores are concatenated.
DATASET_ID_COLUMN = next(iter(COLUMNS))

# The dataset's assayed genes. `X` cannot show them, because an unmeasured gene is 0.
# It is the only part of `var` that differs between stores of one organism.
MEASURED_COLUMN = "measured"
