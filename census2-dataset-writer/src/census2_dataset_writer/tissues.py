"""The tissue roll-up that gives a store's `tissue_general` columns.

Census 1 maps each `tissue_ontology_term_id` to a high-level tissue from the
curated `TISSUE_GENERAL` list of `cellxgene_ontology_guide`, and this module
repeats that mapping (its `TissueMapper`, in `build_soma/tissue_mapper.py` of
the `chanzuckerberg/cellxgene-census` repository). A term with no high-level
tissue maps to itself.

The package holds its ontologies in its own files and reads no network, so a
HealthOmics task can run it. `pyproject.toml` pins the version that Census 1
pins, and the receipt records it.
"""

from __future__ import annotations

import functools
from collections.abc import Iterable

from cellxgene_ontology_guide import curated_ontology_term_lists
from cellxgene_ontology_guide.entities import CuratedOntologyTermList
from cellxgene_ontology_guide.ontology_parser import OntologyParser

# The CxG schema whose ontology versions the mapping reads: the schema of
# Census 1 LTS `2025-11-08` (its `CXG_SCHEMA_VERSION`).
SCHEMA_VERSION = "7.0.0"


@functools.cache
def _ontology() -> tuple[OntologyParser, list[str]]:
    """The parser and the high-level tissues. Loading them takes a few seconds."""
    parser = OntologyParser(f"v{SCHEMA_VERSION}")
    high_level = curated_ontology_term_lists.get_curated_ontology_term_list(
        CuratedOntologyTermList.TISSUE_GENERAL
    )
    return parser, list(high_level)


def tissue_general(term_ids: Iterable[str]) -> dict[str, tuple[str, str]]:
    """Map each tissue term to its high-level tissue term and that term's label.

    Raises `ValueError` for a term that the ontology does not hold. Census 1
    fails its build on such a term.
    """
    parser, high_level = _ontology()
    mapping: dict[str, tuple[str, str]] = {}
    for term in term_ids:
        try:
            general = parser.get_high_level_terms(term, high_level)[0]
        except (IndexError, KeyError, ValueError):
            general = term
        try:
            mapping[term] = (general, parser.get_term_label(general))
        except KeyError as e:
            raise ValueError(f"tissue term {term} is not in the ontology") from e
    return mapping
