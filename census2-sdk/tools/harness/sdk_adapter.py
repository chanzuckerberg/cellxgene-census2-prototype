"""The Census 2 SDK as the stack under test of the validation harness (task 6.2).

`answer` takes one bound corpus query and the clamp of its recorded Census 1 result. It
translates the query's TileDB-SOMA `value_filter` into the SDK's filter mapping, runs
`get_anndata` over the catalog, and returns the canonical arrays of the result. The
comparator then reads the two sides in the same form.

The SDK takes equality and membership on the indexed catalog columns and on
`dataset_id`. A query that it cannot express raises `Unsupported` with the reason, so
the driver lists the query instead of a false difference. The SDK has no `get_obs`,
`get_var`, or presence accessor and no var filter. A var filter on `feature_id` is
applied to the result after the read, and the notes of the answer say so.
"""

import ast
import re
from collections.abc import Sequence
from typing import Any

import pandas as pd

import census2_sdk
from census2_sdk.gene_reference import ORGANISMS
from harness.canonical import Arrays, canonicalize_anndata
from harness.corpus import BoundQuery

CLAUSE = re.compile(r"^(?P<column>\w+) (?P<op>==|in) (?P<value>.+)$")


class Unsupported(Exception):
    """The SDK cannot express this query. The message names the reason."""


class NoCells(Exception):
    """The SDK selected no cell, so `get_anndata` returned None."""


def _values(text: str) -> list[str]:
    try:
        value = ast.literal_eval(text)
    except (SyntaxError, ValueError) as error:
        raise Unsupported(f"the filter value is not a literal: {text!r}") from error
    return [str(v) for v in (value if isinstance(value, list) else [value])]


def clauses(value_filter: str) -> list[str]:
    """The clauses of a conjunction, split on ` and ` outside quotes and brackets."""
    parts, start, quote, depth, i = [], 0, "", 0, 0
    while i < len(value_filter):
        c = value_filter[i]
        if quote:
            if c == "\\":
                i += 1  # the escaped character belongs to the string
            elif c == quote:
                quote = ""
        elif c in "'\"":
            quote = c
        elif c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif depth == 0 and value_filter.startswith(" and ", i):
            parts.append(value_filter[start:i])
            start = i + len(" and ")
            i = start
            continue
        i += 1
    parts.append(value_filter[start:])
    return [p.strip() for p in parts]


def translate(value_filter: str, pins: dict[str, Any]) -> dict[str, list[str]]:
    """The SDK filter mapping of a conjunction of equality and membership clauses.

    A multi-valued column takes the bare term of the pin. The Census 1 filter lists the
    joined values that hold the term (`<column>_matches` in the pins), because Census 1
    compares the whole string.
    """
    mapping: dict[str, list[str]] = {}
    for clause in clauses(value_filter):
        match = CLAUSE.match(clause)
        if match is None:
            raise Unsupported(
                f"the SDK filter takes equality and membership only: {clause!r}"
            )
        column, values = match["column"], _values(match["value"])
        if column in census2_sdk.catalog.MULTIVALUED:
            if values != pins.get(f"{column}_matches"):
                raise Unsupported(f"no bare term is pinned for {column} = {values}")
            values = [pins[column]]
        # A conjunction keeps the values that every clause of the column accepts.
        if column in mapping:
            values = [v for v in mapping[column] if v in values]
        mapping[column] = values
    return mapping


def _genes(var_filter: str) -> list[str]:
    match = CLAUSE.match(var_filter.strip())
    if match is None or match["column"] != "feature_id":
        raise Unsupported(f"the SDK has no var filter: {var_filter!r}")
    return _values(match["value"])


def answer(
    query: BoundQuery,
    census: census2_sdk.Census,
    clamp: Sequence[str],
    pins: dict[str, Any],
) -> tuple[Arrays, list[str]]:
    """The canonical arrays of the SDK's answer to `query`, and notes on how it ran."""
    if query.accessor != "get_anndata":
        raise Unsupported(f"the SDK has no {query.accessor} accessor")
    if query.modality != "census_data":
        raise Unsupported("the SDK reads no spatial experiment")
    organism = query.experiment
    if organism not in ORGANISMS:
        raise Unsupported(f"the SDK pins no gene reference for {organism}")

    obs_filter = translate(query.obs_filter, pins) if query.obs_filter else {}
    clamped = set(clamp)
    if "dataset_id" in obs_filter:
        clamped &= set(obs_filter["dataset_id"])
    obs_filter["dataset_id"] = sorted(clamped)
    try:  # the SDK refuses a column that the catalog does not index
        census2_sdk.catalog.canonicalize(obs_filter)
    except ValueError as error:
        raise Unsupported(str(error)) from error
    # Any other error of the read is a defect of the stores or the SDK, so it fails the run.
    layers = [name for name in query.layers if name != "raw"]
    adata = census.get_anndata(organism, obs_filter, layers=layers)
    if adata is None:
        raise NoCells("the SDK selected no cell")

    notes = []
    if query.var_filter:
        genes = set(_genes(query.var_filter))
        keep = pd.Series(adata.var["feature_id"]).isin(genes).to_numpy()
        adata = adata[:, keep].copy()
        notes.append("genes selected after the read; the SDK has no var filter")
    adata.uns["x_layer_name"] = "raw"
    return canonicalize_anndata(adata), notes
