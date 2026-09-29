"""The fixed query corpus (task 2.1, evidence section 11.2).

Every later phase is checked against this set of reads. Phase 4 signs off value
identity over it. Phase 5 proves that exclusion drops no dataset over it. Phase 6
proves result identity over it. Phase 13 requires R and Python to agree on it.

Each entry names a Census 1 SDK accessor (`get_anndata`, `get_obs`, `get_var`,
`get_presence_matrix`) and the phases that rely on it. The filter values are pinned in
`query_values.json` from release 2025-11-08. A filter that matches nothing passes
vacuously, so no value is written inline.

The filters are TileDB-SOMA `value_filter` strings, the language the oracle speaks. A
stack under test translates them. The Census 2 SDK takes a mapping of indexed column to
values, and its grammar is open (design section 8), so its adapter answers the entries it
can express and reports the rest as unsupported.
"""

from __future__ import annotations

import json
import pathlib
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

CORPUS_VERSION = "1.0"
VALUES_PATH = pathlib.Path(__file__).parent / "query_values.json"

NON_SPATIAL = (
    "homo_sapiens",
    "mus_musculus",
    "macaca_mulatta",
    "callithrix_jacchus",
    "pan_troglodytes",
)
SPATIAL = ("homo_sapiens", "mus_musculus")

# WMG's production filter, verbatim from backend/wmg/pipeline/constants.py in
# single-cell-data-portal: a real consumer's 4-column predicate at maximum fan-out.
WMG_SYSTEM_TISSUES = (
    "UBERON:0001017",
    "UBERON:0001007",
    "UBERON:0002405",
    "UBERON:0000990",
    "UBERON:0001004",
    "UBERON:0001434",
)
WMG_FILTER = (
    "is_primary_data == True and nnz >= 500 and cell_type_ontology_term_id != 'unknown'"
    + "".join(
        f" and tissue_general_ontology_term_id != '{t}'" for t in WMG_SYSTEM_TISSUES
    )
)


@dataclass(frozen=True)
class Entry:
    """One corpus entry, before it is bound to an experiment.

    `obs_filter` and `var_filter` are SOMA `value_filter` strings. A `{placeholder}` is
    filled for each experiment from `query_values.json`.
    """

    key: str
    accessor: str  # get_anndata | get_obs | get_var | get_presence_matrix
    exercises: str
    relied_on_by: tuple[str, ...]
    obs_filter: str | None = None
    var_filter: str | None = None
    layers: tuple[str, ...] = ("raw",)
    genes: str = "all"  # "all" | "subset"
    scope: str = "both"  # "non_spatial" | "spatial" | "both" | "human_non_spatial"
    obs_columns: tuple[str, ...] = ()
    pair: str | None = None


CORPUS: tuple[Entry, ...] = (
    # ---- documented-usage anchors (evidence 11.2: the baseline the no-regression
    # threshold is set against). Filters verbatim from the Census 1 notebooks, so they
    # bind to label columns with inline values and are scoped to human, where the
    # documented usage lives.
    Entry(
        key="documented_narrow",
        accessor="get_anndata",
        exercises="the narrowest documented filter: 431,901 cells across 10 datasets (1%), "
        "the low end of the fan-out range that decides store organization (evidence 10.5, 11.2)",
        relied_on_by=("4.9", "5", "6", "13"),
        obs_filter="tissue_general == 'vasculature'",
        scope="human_non_spatial",
    ),
    Entry(
        key="documented_wide",
        accessor="get_anndata",
        exercises="the widest documented filter: 58,762,981 cells across 447 datasets (41%), "
        "the high end of the same fan-out range and the wide-filter case Phase 6's "
        "concurrency runs target",
        relied_on_by=("5", "6"),
        obs_filter="tissue_general in ['brain', 'lung']",
        scope="human_non_spatial",
    ),
    Entry(
        key="client_documented_two_column",
        accessor="get_anndata",
        exercises="the two-column shape the SDK's own documentation uses (is_primary_data "
        "== True and tissue_general == X, 9.2x cost cliff measured in evidence 9.3)",
        relied_on_by=("4.9", "5", "6", "13"),
        obs_filter="is_primary_data == True and "
        "tissue_general_ontology_term_id == '{tissue_general_primary}'",
        scope="non_spatial",
    ),
    # ---- predicate and locality shapes
    Entry(
        key="wmg_production",
        accessor="get_anndata",
        exercises="a real 4-column predicate at maximum fan-out; reads every gene, raw only",
        relied_on_by=("5", "6"),
        obs_filter=WMG_FILTER,
        scope="non_spatial",
    ),
    Entry(
        key="cell_type_slice",
        accessor="get_anndata",
        exercises="selective read with a dictionary-encoded predicate, maximally scattered "
        "because rows are ordered by soma_joinid, not cell type",
        relied_on_by=("4.9", "5", "6", "13"),
        obs_filter="cell_type_ontology_term_id == '{cell_type_ontology_term_id}'",
        pair="row_locality",
    ),
    Entry(
        key="dataset_slice",
        accessor="get_anndata",
        exercises="the same shape against a near-contiguous row range, isolating scatter "
        "from volume; also the per-dataset identity read that Census 2 stores answer natively",
        relied_on_by=("4.9", "6"),
        obs_filter="dataset_id == '{dataset_id}'",
        pair="row_locality",
    ),
    Entry(
        key="tissue_slice",
        accessor="get_anndata",
        exercises="large contiguous-ish read over a broad predicate",
        relied_on_by=("5", "6"),
        obs_filter="tissue_general_ontology_term_id == '{tissue_general_ontology_term_id}'",
    ),
    # ---- gene-axis shapes
    Entry(
        key="single_gene_all_cells",
        accessor="get_anndata",
        exercises="worst-case sparse column access: one gene, every cell, no obs filter "
        "(evidence 11.2 corpus table)",
        relied_on_by=("6",),
        var_filter="feature_id == '{feature_id_sample}'",
        scope="non_spatial",
    ),
    Entry(
        key="cell_type_gene_subset",
        accessor="get_anndata",
        exercises="the cell_type_slice cells over ~100 genes, so column-spanning "
        "amplification separates from row selection",
        relied_on_by=("6",),
        obs_filter="cell_type_ontology_term_id == '{cell_type_ontology_term_id}'",
        genes="subset",
        pair="column_locality",
    ),
    # ---- silent-failure traps (evidence 11.2: these fail by returning a plausible
    # subset and raising nothing)
    Entry(
        key="multivalued_disease",
        accessor="get_anndata",
        exercises="disease fields hold ' || '-joined values, so equality is membership; "
        "a translator doing string equality returns a subset silently (schema 2.4.0 "
        "made disease_ontology_term_id multi-valued too)",
        relied_on_by=("6",),
        obs_filter="disease_ontology_term_id == '{disease_ontology_term_id}'",
        scope="non_spatial",
    ),
    Entry(
        key="var_symbol_collision",
        accessor="get_var",
        exercises="gene lookup by a non-unique feature_name must surface every colliding "
        "feature_id, not silently pick one",
        relied_on_by=("6", "13"),
        var_filter="feature_name == '{colliding_feature_name}'",
        scope="non_spatial",
    ),
    # ---- layer and modality coverage
    Entry(
        key="both_layers",
        accessor="get_anndata",
        exercises="raw and normalized together; spatial carries raw only, so layer sets "
        "genuinely differ by experiment",
        relied_on_by=("4.9", "6"),
        obs_filter="cell_type_ontology_term_id == '{cell_type_ontology_term_id}'",
        layers=("raw", "normalized"),
        scope="non_spatial",
    ),
    Entry(
        key="spatial_cell_type",
        accessor="get_anndata",
        exercises="the spatial experiments' obs and X: three extra obs columns, raw only",
        relied_on_by=("3", "4.9", "6"),
        obs_filter="cell_type_ontology_term_id == '{cell_type_ontology_term_id}'",
        scope="spatial",
    ),
    # ---- metadata-only accessors
    Entry(
        key="obs_metadata_only",
        accessor="get_obs",
        exercises="a metadata-only read with a column projection: no X touched, which is "
        "the read the catalog excludes and the shape evidence 9.3 costs per column",
        relied_on_by=("5", "6"),
        obs_filter="cell_type_ontology_term_id == '{cell_type_ontology_term_id}'",
        obs_columns=(
            "cell_type_ontology_term_id",
            "tissue_general_ontology_term_id",
            "dataset_id",
        ),
    ),
    Entry(
        key="var_all_genes",
        accessor="get_var",
        exercises="the full var frame: under Census 2 this is the pinned reference itself "
        "(plan assumption C), so its content changes meaning between stacks",
        relied_on_by=("4.9", "6"),
    ),
    Entry(
        key="presence_matrix",
        accessor="get_presence_matrix",
        exercises="the datasets x genes presence matrix: without it, 'not measured' and "
        "'measured as zero' are indistinguishable",
        relied_on_by=("4", "5", "6"),
    ),
)

# Required by evidence 11.2 and not expressible yet. Listed here so the gap is visible.
DEFERRED: tuple[tuple[str, str], ...] = (
    (
        "identity_coords",
        (
            "positional access by published identity: soma_joinid on Census 1, its "
            "replacement depends on decision 7, so the Census 2 half cannot be written yet"
        ),
    ),
    (
        "spatial_scene",
        (
            "scenes, point clouds, images: no AnnData form, needs the Phase 3 "
            "spatialdata decision and its own comparison path"
        ),
    ),
)


@dataclass(frozen=True)
class BoundQuery:
    """One corpus entry bound to one experiment, with its placeholders filled."""

    key: str
    accessor: str
    modality: str
    experiment: str
    relied_on_by: tuple[str, ...]
    obs_filter: str | None
    var_filter: str | None
    layers: tuple[str, ...]
    obs_columns: tuple[str, ...] = ()
    label: str = ""


def load_values(path: pathlib.Path = VALUES_PATH) -> dict[str, Any]:
    values: dict[str, Any] = json.loads(path.read_text())
    return values


def clamp_filter(obs_filter: str | None, dataset_ids: Iterable[str]) -> str:
    """Bound a query to the snapshot's datasets.

    The oracle runner and the stack-under-test runner share this one definition. Two
    clamps that differ would report phantom row differences on every comparison.
    """
    clause = "dataset_id in [" + ", ".join(f"'{d}'" for d in dataset_ids) + "]"
    return f"({obs_filter}) and {clause}" if obs_filter else clause


def _applies(entry: Entry, modality: str, experiment: str) -> bool:
    if entry.scope == "both":
        return True
    if entry.scope == "human_non_spatial":
        return modality == "census_data" and experiment == "homo_sapiens"
    spatial = modality == "census_spatial_sequencing"
    return entry.scope == ("spatial" if spatial else "non_spatial")


def expand(values: dict[str, Any] | None = None) -> list[BoundQuery]:
    """Bind every entry to every experiment it applies to.

    All seven experiments are visited. An entry whose placeholders cannot be filled for
    an experiment is skipped, so no query goes out with a filter that selects nothing.
    """
    values = values if values is not None else load_values()
    out: list[BoundQuery] = []
    modalities = (("census_data", NON_SPATIAL), ("census_spatial_sequencing", SPATIAL))
    for modality, experiments in modalities:
        for experiment in experiments:
            pinned = values.get(f"{modality}/{experiment}")
            if pinned is None:
                continue
            for e in CORPUS:
                if not _applies(e, modality, experiment):
                    continue
                try:
                    obs_filter = e.obs_filter.format(**pinned) if e.obs_filter else None
                    var_filter = e.var_filter.format(**pinned) if e.var_filter else None
                except KeyError:
                    continue
                if e.genes == "subset":
                    subset = pinned.get("gene_subset")
                    if not subset:
                        continue
                    var_filter = (
                        "feature_id in [" + ", ".join(f"'{g}'" for g in subset) + "]"
                    )
                out.append(
                    BoundQuery(
                        key=e.key,
                        accessor=e.accessor,
                        modality=modality,
                        experiment=experiment,
                        relied_on_by=e.relied_on_by,
                        obs_filter=obs_filter,
                        var_filter=var_filter,
                        layers=e.layers,
                        obs_columns=e.obs_columns,
                        label=f"{e.key}:{modality}/{experiment}",
                    )
                )
    return out


if __name__ == "__main__":
    queries = expand()
    print(
        f"corpus v{CORPUS_VERSION}, source release {load_values()['_source']['release']}"
    )
    print(
        f"{len(CORPUS)} entries -> {len(queries)} bound queries; {len(DEFERRED)} deferred"
    )
    by_accessor: dict[str, int] = {}
    for q in queries:
        by_accessor[q.accessor] = by_accessor.get(q.accessor, 0) + 1
    print("by accessor:", by_accessor)
