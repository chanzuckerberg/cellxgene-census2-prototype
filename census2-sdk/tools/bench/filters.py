"""The named `obs` queries of the filter benchmarks, in builder form.

Each query builds its filter from the typed columns of one organism (`obs`), so the
same function runs on `census.homo_sapiens.obs` of this SDK and on the builder of
the earlier filter spike (`spike/consolidation-filters`), which has the same
attributes and operators. The module imports no SDK, so a checkout of another
branch can import it. `mapping` is the same filter as a mapping, for
`Organism.where` and for `main`'s `get_anndata`, when a mapping can express it.

`mode` tells how a cross-check selects the cells: `anndata` reads `X`, and `keys`
reads only the cell keys, because the query selects hundreds of thousands of cells
or more. The tags name the query sets: `pr` for the six queries that earlier pull
requests measured, `shared` for the shared set of the filter spikes, `e2e` for the
end-to-end set of the pyarrow spike, and `catalog` for a query whose terms the
catalog cannot use.
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

LUNG = [
    "epithelial cell of lung",
    "lower lobe of left lung",
    "lung",
    "lung epithelium",
    "middle lobe of right lung",
    "upper lobe of left lung",
]
MACROPHAGES = [
    "TREM2-positive macrophage",
    "alveolar macrophage",
    "central nervous system macrophage",
    "gastrointestinal tract (lamina propria) macrophage",
    "inflammatory macrophage",
    "kidney resident macrophage",
    "lung interstitial macrophage",
    "lung macrophage",
    "macrophage",
]
# The lung macrophages of the consolidation spike: six of the nine.
LUNG_MACROPHAGES = [
    "TREM2-positive macrophage",
    "alveolar macrophage",
    "inflammatory macrophage",
    "lung interstitial macrophage",
    "lung macrophage",
    "macrophage",
]
WMG_EXCLUDED = [
    "UBERON:0001017",
    "UBERON:0001007",
    "UBERON:0002405",
    "UBERON:0000990",
    "UBERON:0001004",
    "UBERON:0001434",
]
ONE = "f75f2ff4-2884-4c2d-b375-70de37a34507"
MEMBER = "6ceeaa86-9ceb-4582-b390-6d4dd6ff0572"
ALONE = "c6413315-7158-46c5-a5fd-7228385694f1"
CD4 = "CD4-positive, alpha-beta T cell"


@dataclass(frozen=True)
class Query:
    name: str
    build: Callable[[Any], Any]
    mapping: dict[str, Any] | None = None
    mode: str = "anndata"
    tags: tuple[str, ...] = field(default_factory=tuple)


def _wmg(o: Any) -> Any:
    return (
        o.is_primary_data
        & (o.nnz >= 500)
        & (o.cell_type_ontology_term_id != "unknown")
        & ~o.tissue_general_ontology_term_id.isin(WMG_EXCLUDED)
    )


def _one(name: str, column: str, value: str, *tags: str) -> Query:
    """A query of one value of one column, which a mapping also expresses."""
    multi = column in ("disease", "disease_ontology_term_id")

    def build(o: Any) -> Any:
        col = getattr(o, column)
        return col.has_any([value]) if multi else col == value

    return Query(name, build, {column: value}, tags=tags)


QUERIES: dict[str, Query] = {
    q.name: q
    for q in [
        _one("liver", "tissue", "liver", "pr", "e2e"),
        _one("crohn", "disease", "Crohn disease", "pr", "e2e"),
        _one("assay_v1", "assay", "10x 3' v1", "pr", "e2e"),
        _one("neutrophil", "cell_type", "neutrophil", "pr", "shared", "e2e"),
        _one("macrophage", "cell_type", "macrophage", "pr", "e2e"),
        _one("cd4", "cell_type", CD4, "pr", "e2e"),
        _one("one_dataset", "dataset_id", ONE, "shared"),
        _one("single_member", "dataset_id", MEMBER),
        _one("single_alone", "dataset_id", ALONE),
        Query(
            "lung_macrophage_primary",
            lambda o: (
                o.tissue.isin(LUNG) & o.cell_type.isin(MACROPHAGES) & o.is_primary_data
            ),
            tags=("shared",),
        ),
        Query(
            "lung_macrophage",
            lambda o: o.tissue.isin(LUNG) & o.cell_type.isin(LUNG_MACROPHAGES),
            {"tissue": LUNG, "cell_type": LUNG_MACROPHAGES},
        ),
        Query(
            "fp_none",
            lambda o: (
                (o.tissue == "mesenteric lymph node") & (o.cell_type == "megakaryocyte")
            ),
            {"tissue": "mesenteric lymph node", "cell_type": "megakaryocyte"},
            tags=("shared", "e2e"),
        ),
        Query(
            "fp_some",
            lambda o: (o.tissue == "liver") & (o.cell_type == "mesothelial cell"),
            {"tissue": "liver", "cell_type": "mesothelial cell"},
            tags=("shared", "e2e"),
        ),
        Query(
            "false_positive",
            lambda o: (o.assay == "10x 3' v3") & (o.tissue == "mesenteric lymph node"),
            {"assay": "10x 3' v3", "tissue": "mesenteric lymph node"},
        ),
        _one("tissue_general", "tissue_general", "pancreas"),
        Query(
            "numeric",
            lambda o: (o.cell_type == "neutrophil") & (o.nnz >= 1000),
            tags=("e2e",),
        ),
        Query(
            # A narrow range: nnz >= 8000 selects 279,441 dense cells, and their X
            # does not fit in the memory of the bench instance.
            "numeric_only",
            lambda o: o.nnz.between(8000, 8100),
            tags=("e2e", "catalog"),
        ),
        Query(
            "membership",
            lambda o: o.disease.has_any(["COVID-19"]) & (o.cell_type == "neutrophil"),
            {"disease": "COVID-19", "cell_type": "neutrophil"},
            tags=("e2e",),
        ),
        Query(
            # Cells of unknown sex that are not primary data: two negated terms,
            # which give no hint.
            "negation",
            lambda o: ~o.sex.isin(["female", "male"]) & ~o.is_primary_data,
            tags=("e2e", "catalog"),
        ),
        Query(
            "prefix_small",
            lambda o: o.cell_type.str.startswith("mesothelial"),
            tags=("e2e", "catalog"),
        ),
        Query(
            "or_columns",
            lambda o: (o.cell_type == "neutrophil") | (o.tissue_general == "pancreas"),
            tags=("e2e", "catalog"),
        ),
        Query(
            "compound",
            lambda o: (
                ((o.tissue == "lung") | (o.tissue_general == "lung"))
                & ~(o.cell_type == "macrophage")
            ),
            mode="keys",
            tags=("catalog",),
        ),
        Query(
            "has_any",
            lambda o: o.disease_ontology_term_id.has_any(
                ["MONDO:0100096", "MONDO:0005109"]
            ),
            mode="keys",
        ),
        Query(
            "has_all",
            lambda o: o.disease.has_all(
                ["HIV infectious disease", "visceral leishmaniasis"]
            ),
        ),
        Query(
            "numeric_range",
            lambda o: (
                (o.cell_type == "neutrophil") & (o.nnz >= 2000) & (o.raw_sum < 20000)
            ),
        ),
        Query(
            "between",
            lambda o: o.nnz.between(2000, 3000) & (o.tissue_general == "blood"),
            mode="keys",
        ),
        Query(
            "not_in",
            lambda o: (
                (o.tissue_general == "lung")
                & ~o.cell_type.isin(["alveolar macrophage", "type II pneumocyte"])
            ),
            mode="keys",
        ),
        Query("prefix", lambda o: o.cell_type.str.startswith("CD8"), mode="keys"),
        Query(
            "substring",
            lambda o: o.cell_type.str.contains("macrophage", case=False, regex=False),
            mode="keys",
        ),
        Query(
            "regex", lambda o: o.cell_type.str.match(r"CD(4|8)-positive"), mode="keys"
        ),
        Query("isna", lambda o: o.raw_variance_nnz.isna()),
        Query("wmg", _wmg, mode="keys", tags=("shared",)),
        Query("all_obs", lambda o: None, mode="keys", tags=("shared",)),
    ]
}


def tagged(tag: str) -> list[Query]:
    return [q for q in QUERIES.values() if tag in q.tags]
