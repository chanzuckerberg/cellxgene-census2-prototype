"""The uniform category lists of a store's `obs`.

A column of `store_format.UNIFORM_OBS_ONTOLOGIES` holds one category list in
every store. The list holds each term of its ontologies that is not deprecated,
as an identifier or as a label, sorted. The lists come from the pinned
`cellxgene_ontology_guide`, so a store depends only on its dataset and on that
version. `write` records a digest of each list in the store (`list_digests`), so
a reader can tell which stores hold equal lists.

The list of a label column is sorted by label, so it is not parallel to the list
of the identifier column. Two ontologies of one pair can share a label, and the
label list holds it once.
"""

from __future__ import annotations

import functools
import hashlib
import json
from typing import Any, cast

import pandas as pd

from . import tissues
from .store_format import UNIFORM_OBS_ONTOLOGIES, UNIFORM_OBS_UNKNOWN

_ID_SUFFIX = "_ontology_term_id"


@functools.cache
def uniform_categories(column: str) -> tuple[str, ...] | None:
    """The category list of `obs` column `column`, or None for a column without one."""
    pair = column.removesuffix(_ID_SUFFIX)
    ontologies = UNIFORM_OBS_ONTOLOGIES.get(pair)
    if ontologies is None:
        return None
    parser, _ = tissues._ontology()
    values: set[str] = set()
    for ontology in ontologies:
        terms = cast("dict[str, dict[str, Any]]", parser.cxg_schema.ontology(ontology))
        for term, entry in terms.items():
            if not entry.get("deprecated"):
                values.add(term if column.endswith(_ID_SUFFIX) else str(entry["label"]))
    if pair in UNIFORM_OBS_UNKNOWN:
        values.add("unknown")
    return tuple(sorted(values))


def list_digests(obs: pd.DataFrame) -> dict[str, str]:
    """A digest of the category list of each uniform column of `obs`.

    `write` records the digests in the store. Two stores with one digest for a
    column hold equal lists for it, so a code has one meaning in both. The digest
    comes from the list that `obs` holds, and not from the ontology package, so
    it is true for each `obs` that a caller gives to `write`. A list that holds
    another type than strings gets no digest.
    """
    digests = {}
    for name in obs.columns:
        column = obs[name]
        if uniform_categories(str(name)) is None or not isinstance(
            column.dtype, pd.CategoricalDtype
        ):
            continue
        categories = list(column.cat.categories)
        # A list that holds another type than strings gets no digest, because
        # its JSON form can equal the form of a list of strings.
        if not all(isinstance(c, str) for c in categories):
            continue
        # A JSON array has one form for one list, whatever the strings hold.
        encoded = json.dumps(categories).encode()
        digests[str(name)] = hashlib.sha256(encoded).hexdigest()
    return digests
