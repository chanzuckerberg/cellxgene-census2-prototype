"""The prototype catalog: one row for each store, as Parquet (design section 2).

`build` reads a prefix that holds `<dataset_version_id>.zarr` stores and the receipts
next to them. `exclude` applies the one correctness rule of design 2.2: a query
excludes a dataset only when the catalog proves that no cell of the dataset can match.
An absent column or an unknown value keeps the dataset.

Command line, for the Iteration 1 store set:

    census2-catalog s3://bucket/prefix/ catalog.parquet
"""

import io
import json
import os
import pathlib
import sys
from collections.abc import Iterable, Mapping
from typing import Any, cast

import pandas as pd
import zarr
from anndata.io import read_elem

from census2_sdk.gene_reference import TAXON_TO_ORGANISM

# obstore is the optional `s3` extra. It is imported here, at module level, when it
# is installed, so that a local catalog over local stores installs and runs
# without it. Each remote path calls `_require_obstore` first.
try:
    import obstore
    import obstore.store
except ImportError:  # pragma: no cover
    obstore = None  # type: ignore[assignment]


def _require_obstore() -> None:
    if obstore is None:
        raise ImportError("a remote URI needs obstore: pip install 'census2-sdk[s3]'")


Filter = Mapping[str, str | Iterable[str]]
# A filter after `normalize`: one set of accepted values for each column.
Normalized = dict[str, frozenset[str]]

# Design 2.2: the 18 indexed columns.
PAIRS = (
    "assay",
    "cell_type",
    "development_stage",
    "disease",
    "self_reported_ethnicity",
    "sex",
    "tissue",
    "tissue_general",
)
INDEXED = tuple(c for p in PAIRS for c in (p, f"{p}_ontology_term_id")) + (
    "suspension_type",
    "tissue_type",
)
# Schema 2.4.0 made these two columns multi-valued. Equality is membership.
MULTIVALUED = frozenset({"disease", "disease_ontology_term_id"})
MULTIVALUE_SEP = " || "
# Assays that the curation schema lists as spatial.
SPATIAL_ASSAYS = frozenset({"EFO:0010961", "EFO:0030062"})

REMOTE_SCHEMES = ("s3://", "gs://", "az://", "abfs://", "http://", "https://")


def is_remote(uri: str | os.PathLike[str]) -> bool:
    return str(uri).startswith(REMOTE_SCHEMES)


def open_store(uri: str | os.PathLike[str]) -> zarr.Group:
    """Open a store for reading. A remote URI goes through obstore."""
    if is_remote(uri):
        _require_obstore()
        store = zarr.storage.ObjectStore(
            obstore.store.from_url(str(uri)), read_only=True
        )
        return zarr.open_group(store, mode="r")
    return zarr.open_group(str(uri), mode="r")


def _list(prefix: str) -> list[str]:
    """Names of the objects directly under `prefix`."""
    if is_remote(prefix):
        _require_obstore()
        store = obstore.store.from_url(prefix)
        listing = obstore.list_with_delimiter(store, "")
        return [str(o["path"]) for o in listing["objects"]]
    return [p.name for p in pathlib.Path(prefix).iterdir()]


def _read_json(uri: str) -> dict[str, Any]:
    if is_remote(uri):
        _require_obstore()
        base, _, name = uri.rpartition("/")
        store = obstore.store.from_url(base)
        data = bytes(obstore.get(store, name).bytes())
    else:
        data = pathlib.Path(uri).read_bytes()
    return cast("dict[str, Any]", json.loads(data))


def members(value: str, column: str) -> list[str]:
    """The values that one cell holds in `column`, one for a single-valued column."""
    if column in MULTIVALUED:
        return [m.strip() for m in value.split(MULTIVALUE_SEP)]
    return [value]


def normalize(obs_filter: Filter | None) -> Normalized:
    """One set of accepted values for each column. Refuse a column that is not indexed."""
    normalized: Normalized = {}
    for column, value in (obs_filter or {}).items():
        if column != "dataset_id" and column not in INDEXED:
            raise ValueError(
                f"{column!r} is not a catalog column; use one of {INDEXED}"
            )
        wanted = {value} if isinstance(value, str) else {str(v) for v in value}
        normalized[column] = frozenset(wanted)
    return normalized


def build(prefix: str, out: str | os.PathLike[str] | None = None) -> pd.DataFrame:
    """One catalog row for each admitted receipt under `prefix`.

    The receipt supplies `dataset_id`; the store supplies everything else.
    """
    # ponytail: dataset_id comes from the receipt until the writer stamps it into obs
    # (disposition row 13).
    prefix = prefix.rstrip("/")
    rows: list[dict[str, Any]] = []
    for name in sorted(n for n in _list(prefix) if n.endswith(".json")):
        receipt = _read_json(f"{prefix}/{name}")
        if not receipt.get("admitted"):
            continue
        version_id = str(receipt["dataset_version_id"])
        store_uri = f"{prefix}/{version_id}.zarr"
        rows.append(
            {
                "dataset_id": str(receipt["dataset_id"]),
                "dataset_version_id": version_id,
                "store_uri": store_uri,
                "curation_schema_version": receipt.get("schema_version"),
                # The Iteration 1 stores carry the source var and record no version.
                "gene_reference_version": None,
                **summarize(open_store(store_uri)),
            }
        )
    frame = pd.DataFrame(rows)
    if out is not None:
        frame.to_parquet(out, index=False)
    return frame


def summarize(group: zarr.Group) -> dict[str, Any]:
    """The catalog fields that one store supplies."""
    obs = group["obs"]
    assert isinstance(obs, zarr.Group)
    row: dict[str, Any] = {}
    n_obs: int | None = None
    # ponytail: list<string> columns; bitmaps over a dictionary (design 2.2) when the
    # catalog size matters.
    for column in INDEXED:
        if column not in obs:
            row[column] = None
            continue
        values = pd.Series(read_elem(obs[column]))
        n_obs = len(values)
        row[column] = sorted(values.astype(str).unique())
    if "is_primary_data" in obs:
        primary = pd.Series(read_elem(obs["is_primary_data"])).astype(bool)
        n_obs = len(primary)
        row["is_primary_data"] = (
            "all_true"
            if primary.all()
            else "all_false"
            if not primary.any()
            else "mixed"
        )
    else:
        row["is_primary_data"] = "mixed"
    if n_obs is None:
        n_obs = len(pd.Series(read_elem(obs[str(obs.attrs["_index"])])))
    row["n_obs"] = n_obs

    uns = group["uns"]
    assert isinstance(uns, zarr.Group)
    taxon = str(read_elem(uns["organism_ontology_term_id"]))
    row["organism"] = TAXON_TO_ORGANISM.get(taxon, taxon)
    assays = row.get("assay_ontology_term_id") or []
    row["spatial"] = bool(SPATIAL_ASSAYS & set(assays))
    return row


def load(uri: str | os.PathLike[str]) -> pd.DataFrame:
    if is_remote(uri):
        _require_obstore()
        base, _, name = str(uri).rpartition("/")
        data = bytes(obstore.get(obstore.store.from_url(base), name).bytes())
        return pd.read_parquet(io.BytesIO(data))
    return pd.read_parquet(uri)


def excludes(row: Mapping[str, Any], obs_filter: Normalized) -> bool:
    """True when the row proves that no cell of the dataset matches `obs_filter`."""
    for column, wanted in obs_filter.items():
        if column == "dataset_id":
            if str(row["dataset_id"]) not in wanted:
                return True
            continue
        values = row[column]
        if values is None or (not isinstance(values, list) and pd.isna(values).all()):
            continue  # the store has no such column; unproven keeps the dataset
        present = {m for v in values for m in members(str(v), column)}
        if not (present & wanted):
            return True
    return False


def exclude(
    catalog: pd.DataFrame, organism: str, obs_filter: Filter | None
) -> pd.DataFrame:
    """The catalog rows that can hold a matching cell."""
    scoped = catalog[catalog["organism"] == organism]
    normalized = normalize(obs_filter)
    if not normalized:
        return scoped
    records = cast("list[dict[str, Any]]", scoped.to_dict("records"))
    keep = [not excludes(row, normalized) for row in records]
    return scoped[keep]


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    frame = build(argv[0], argv[1])
    print(f"{len(frame)} rows -> {argv[1]}")
    return 0


def cli() -> None:
    """The `census2-catalog` console script."""
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    cli()
