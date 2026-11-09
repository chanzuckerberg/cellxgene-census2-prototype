"""The named queries of the filter builder benchmark, in builder form.

Each query builds its filter from the typed columns of one organism, and names its
mode: `anndata` (`get_anndata`), `obs` (`get_obs`), `count`, or `count_by` with its
columns. A filter that selects hundreds of thousands of cells or more counts them by
donor, as the `count_donor` shape of the DuckDB spike does, so that the measure is the
`obs` work and not a read of `X` of many gigabytes. `mapping` is the same filter as a mapping filter of `Census.get_anndata`,
when one can express it, so that `main`'s path runs it too. `sql` is the filter of
the DuckDB path of the `spike/obs-duckdb` branch, for the cross-check of the cell
keys; `spikes/consolidation-filters/duckdb_keys.py` reads it.

The shared queries take their values from the DuckDB and consolidation spikes:
(a) one dataset, (b) neutrophils, (c) lung macrophages of primary data, and the
catalog false positives (d1) and (d2).
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from census2_sdk.obs_columns import Obs
from census2_sdk.obs_filter import FilterLike

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


def _sql_list(values: list[str]) -> str:
    return "(" + ", ".join("'" + v.replace("'", "''") + "'" for v in values) + ")"


@dataclass(frozen=True)
class Query:
    name: str
    build: Callable[[Obs], FilterLike | None]
    mode: str = "anndata"
    by: tuple[str, ...] = ()
    mapping: dict[str, Any] | None = None
    sql: str | None = None
    tags: tuple[str, ...] = field(default_factory=tuple)


def _wmg(obs: Obs) -> FilterLike:
    return (
        obs.is_primary_data
        & (obs.nnz >= 500)
        & (obs.cell_type_ontology_term_id != "unknown")
        & ~obs.tissue_general_ontology_term_id.isin(WMG_EXCLUDED)
    )


WMG_SQL = (
    "is_primary_data AND nnz >= 500 AND cell_type_ontology_term_id <> 'unknown' AND "
    f"tissue_general_ontology_term_id NOT IN {_sql_list(WMG_EXCLUDED)}"
)
ONE = "f75f2ff4-2884-4c2d-b375-70de37a34507"
MEMBER = "6ceeaa86-9ceb-4582-b390-6d4dd6ff0572"
ALONE = "c6413315-7158-46c5-a5fd-7228385694f1"

QUERIES: dict[str, Query] = {
    q.name: q
    for q in [
        Query(
            "one_dataset",
            lambda o: o.dataset_id == ONE,
            mapping={"dataset_id": ONE},
            sql=f"dataset_id = '{ONE}'",
            tags=("shared",),
        ),
        Query(
            "single_member",
            lambda o: o.dataset_id == MEMBER,
            mapping={"dataset_id": MEMBER},
            sql=f"dataset_id = '{MEMBER}'",
        ),
        Query(
            "single_alone",
            lambda o: o.dataset_id == ALONE,
            mapping={"dataset_id": ALONE},
            sql=f"dataset_id = '{ALONE}'",
        ),
        Query(
            "neutrophil",
            lambda o: o.cell_type == "neutrophil",
            mapping={"cell_type": "neutrophil"},
            sql="cell_type = 'neutrophil'",
            tags=("shared",),
        ),
        Query(
            "lung_macrophage_primary",
            lambda o: (
                o.tissue.isin(LUNG) & o.cell_type.isin(MACROPHAGES) & o.is_primary_data
            ),
            sql=(
                f"tissue IN {_sql_list(LUNG)} AND cell_type IN {_sql_list(MACROPHAGES)} "
                "AND is_primary_data"
            ),
            tags=("shared",),
        ),
        Query(
            "lung_macrophage",
            lambda o: o.tissue.isin(LUNG) & o.cell_type.isin(LUNG_MACROPHAGES),
            mapping={"tissue": LUNG, "cell_type": LUNG_MACROPHAGES},
            sql=f"tissue IN {_sql_list(LUNG)} AND cell_type IN {_sql_list(LUNG_MACROPHAGES)}",
        ),
        Query(
            "fp_none",
            lambda o: (
                (o.tissue == "mesenteric lymph node") & (o.cell_type == "megakaryocyte")
            ),
            mapping={"tissue": "mesenteric lymph node", "cell_type": "megakaryocyte"},
            sql="tissue = 'mesenteric lymph node' AND cell_type = 'megakaryocyte'",
            tags=("shared",),
        ),
        Query(
            "fp_some",
            lambda o: (o.tissue == "liver") & (o.cell_type == "mesothelial cell"),
            mapping={"tissue": "liver", "cell_type": "mesothelial cell"},
            sql="tissue = 'liver' AND cell_type = 'mesothelial cell'",
            tags=("shared",),
        ),
        Query(
            "false_positive",
            lambda o: (o.assay == "10x 3' v3") & (o.tissue == "mesenteric lymph node"),
            mapping={"assay": "10x 3' v3", "tissue": "mesenteric lymph node"},
            sql="assay = '10x 3'' v3' AND tissue = 'mesenteric lymph node'",
        ),
        Query(
            "tissue_general",
            lambda o: o.tissue_general == "pancreas",
            mapping={"tissue_general": "pancreas"},
            sql="tissue_general = 'pancreas'",
        ),
        Query(
            "numeric",
            lambda o: (o.cell_type == "neutrophil") & (o.nnz >= 1000),
            sql="cell_type = 'neutrophil' AND nnz >= 1000",
        ),
        Query(
            "membership",
            lambda o: o.disease.has_any(["COVID-19"]) & (o.cell_type == "neutrophil"),
            mapping={"disease": "COVID-19", "cell_type": "neutrophil"},
            sql=(
                "list_contains(string_split(disease, ' || '), 'COVID-19') AND "
                "cell_type = 'neutrophil'"
            ),
        ),
        Query(
            "or_columns",
            lambda o: (o.cell_type == "neutrophil") | (o.tissue_general == "pancreas"),
            sql="cell_type = 'neutrophil' OR tissue_general = 'pancreas'",
            tags=("new",),
        ),
        Query(
            "compound",
            lambda o: (
                ((o.tissue == "lung") | (o.tissue_general == "lung"))
                & ~(o.cell_type == "macrophage")
            ),
            sql="(tissue = 'lung' OR tissue_general = 'lung') AND NOT cell_type = 'macrophage'",
            tags=("new",),
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "has_any",
            lambda o: o.disease_ontology_term_id.has_any(
                ["MONDO:0100096", "MONDO:0005109"]
            ),
            sql=(
                "list_has_any(string_split(disease_ontology_term_id, ' || '), "
                "['MONDO:0100096', 'MONDO:0005109'])"
            ),
            tags=("new",),
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "has_all",
            lambda o: o.disease.has_all(
                ["HIV infectious disease", "visceral leishmaniasis"]
            ),
            sql=(
                "list_has_all(string_split(disease, ' || '), "
                "['HIV infectious disease', 'visceral leishmaniasis'])"
            ),
            tags=("new",),
        ),
        Query(
            "numeric_range",
            lambda o: (
                (o.cell_type == "neutrophil") & (o.nnz >= 2000) & (o.raw_sum < 20000)
            ),
            sql="cell_type = 'neutrophil' AND nnz >= 2000 AND raw_sum < 20000",
            tags=("new",),
        ),
        Query(
            "between",
            lambda o: o.nnz.between(2000, 3000) & (o.tissue_general == "blood"),
            sql="nnz BETWEEN 2000 AND 3000 AND tissue_general = 'blood'",
            tags=("new",),
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "not_in",
            lambda o: (
                (o.tissue_general == "lung")
                & ~o.cell_type.isin(["alveolar macrophage", "type II pneumocyte"])
            ),
            sql=(
                "tissue_general = 'lung' AND cell_type NOT IN "
                "('alveolar macrophage', 'type II pneumocyte')"
            ),
            tags=("new",),
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "prefix",
            lambda o: o.cell_type.str.startswith("CD8"),
            sql="cell_type LIKE 'CD8%'",
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "substring",
            lambda o: o.cell_type.str.contains("macrophage", case=False, regex=False),
            sql="cell_type ILIKE '%macrophage%'",
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "regex",
            lambda o: o.cell_type.str.match(r"CD(4|8)-positive"),
            sql="regexp_matches(cell_type, '^CD(4|8)-positive')",
            mode="count_by",
            by=("donor_id",),
        ),
        Query(
            "isna",
            lambda o: o.raw_variance_nnz.isna(),
            sql="raw_variance_nnz IS NULL OR isnan(raw_variance_nnz)",
        ),
        Query(
            "count_b_donor",
            lambda o: o.cell_type == "neutrophil",
            mode="count_by",
            by=("donor_id",),
            sql="cell_type = 'neutrophil'",
            tags=("shared",),
        ),
        Query(
            "count_all_donor",
            lambda o: None,
            mode="count_by",
            by=("donor_id",),
            sql="true",
            tags=("shared",),
        ),
        Query(
            "count_b",
            lambda o: o.cell_type == "neutrophil",
            mode="count",
            mapping={"cell_type": "neutrophil"},
            sql="cell_type = 'neutrophil'",
        ),
        Query(
            "wmg_count",
            _wmg,
            mode="count_by",
            by=("tissue_general_ontology_term_id", "cell_type_ontology_term_id"),
            sql=WMG_SQL,
            tags=("shared",),
        ),
        Query("wmg", _wmg, mode="obs", sql=WMG_SQL),
        Query("all_obs", lambda o: None, mode="obs", sql="true"),
    ]
}


def run_query(organism: Any, query: Query, max_workers: int) -> int:
    """Run `query` by its mode over `organism`; return the number of cells."""
    obs_filter = query.build(organism.obs)
    if query.mode == "anndata":
        adata = organism.get_anndata(obs_filter, max_workers=max_workers)
        return 0 if adata is None else int(adata.n_obs)
    if query.mode == "obs":
        return len(organism.get_obs(obs_filter, max_workers=max_workers))
    if query.mode == "count":
        return int(organism.count(obs_filter, max_workers=max_workers))
    if query.mode == "count_by":
        counts = organism.count_by(list(query.by), obs_filter, max_workers=max_workers)
        return int(counts.sum())
    raise ValueError(f"no mode {query.mode!r}")


def matrix_queries(kind: str) -> Any:
    """The queries in the form of a `--queries` file.

    `builder` and `mapping` give the list of `bench.matrix`; `duckdb` gives the
    mapping of names to filters of the `spike/obs-duckdb` branch's `bench.matrix`, in
    which a shape `query:<name>`, `obs:<name>`, or `count_donor:<name>` names a filter.
    """
    if kind == "builder":
        return [
            {"query": q.name, "shape": "builder", "builder": q.name}
            for q in QUERIES.values()
        ]
    if kind == "mapping":
        return [
            {"query": q.name, "shape": "filter", "filter": q.mapping}
            for q in QUERIES.values()
            if q.mapping is not None and q.mode == "anndata"
        ]
    if kind == "duckdb":
        return {q.name: {"sql": q.sql} for q in QUERIES.values() if q.sql is not None}
    raise ValueError(f"no kind {kind!r}")


if __name__ == "__main__":
    import json
    import sys

    print(json.dumps(matrix_queries(sys.argv[1]), indent=1))
