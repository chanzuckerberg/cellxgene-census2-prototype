#!/usr/bin/env python3
"""Build the prototype catalog from the Iteration 1 stores: one Parquet row for each store.

The catalog builder is not part of the SDK. The SDK reads a catalog; the Census 2
team publishes one (design 2.3, 3.3). This script is the Iteration 1 stand-in for
that builder. It reads a prefix that holds `<dataset_version_id>.zarr` stores and the
receipts next to them. It takes the indexed values from each store's `obs`, because the
stores have no summary record yet. The real builder reads the Discover API and the
summary records and opens no array data.

    python census2-sdk/tools/build_catalog.py s3://bucket/prefix/ catalog.parquet

A remote prefix needs the `s3` extra and credentials in the environment.
"""

import json
import pathlib
import sys
from typing import Any, cast

import pandas as pd
import zarr
from anndata.io import read_elem

from census2_sdk.catalog import INDEXED, is_remote, open_store
from census2_sdk.gene_reference import TAXON_TO_ORGANISM

# Assays that the curation schema lists as spatial.
SPATIAL_ASSAYS = frozenset({"EFO:0010961", "EFO:0030062"})


def _list(prefix: str) -> list[str]:
    """Names of the objects directly under `prefix`."""
    if is_remote(prefix):
        import obstore
        import obstore.store

        store = obstore.store.from_url(prefix)
        listing = obstore.list_with_delimiter(store, "")
        return [str(o["path"]) for o in listing["objects"]]
    return [p.name for p in pathlib.Path(prefix).iterdir()]


def _read_json(uri: str) -> dict[str, Any]:
    if is_remote(uri):
        import obstore
        import obstore.store

        base, _, name = uri.rpartition("/")
        store = obstore.store.from_url(base)
        data = bytes(obstore.get(store, name).bytes())
    else:
        data = pathlib.Path(uri).read_bytes()
    return cast("dict[str, Any]", json.loads(data))


def build(prefix: str, out: str | pathlib.Path | None = None) -> pd.DataFrame:
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
        # A null is not a value. astype(str) would turn it into the string "nan".
        row[column] = sorted(values.dropna().astype(str).unique())
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
    # The writer records the gene reference that it projected the store onto.
    row["gene_reference_version"] = (
        str(read_elem(uns["gene_reference_version"]))
        if "gene_reference_version" in uns
        else None
    )
    assays = row.get("assay_ontology_term_id") or []
    row["spatial"] = bool(SPATIAL_ASSAYS & set(assays))
    return row


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    frame = build(argv[0], argv[1])
    print(f"{len(frame)} rows -> {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
